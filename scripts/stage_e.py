#!/usr/bin/env python
"""Stage E: simulated null and injected skies through Stages C-D.

    python scripts/stage_e.py --precompute            # pixel templates + truths; one process, <= 8 GB, ~10 min
    python scripts/stage_e.py --pilot                 # one sim per truth, serial, timed
    python scripts/stage_e.py --n-sims 50 --workers 6 # resumable: one JSON per sim
    python scripts/stage_e.py --summarise

THE QUESTION (STATE_2026-09-22 §7 item 24): Totani's residual bubble
construction is a data-driven template, and the structure selection is chosen
on the data. Neither he nor we have a look-elsewhere calibration. Four truths,
each a real-data background-only fit (`simsky.build_truth`) of the best
physical row of one no-residual selection, with and without his halo:

  null_structured  SO_Z6_R30_T150_C5 1-gas + the 1-gas no-residual selection
                   (structured -> eROSITA unif N/S -> Haslam -> solar IC)
  halo_structured  the same + NFW-rho^2 on Totani's fitted spectrum
  null_flat        SS_Z10_R30_T150_C5 5-ring + the 5-ring no-residual selection
                   (flat -> eROSITA unif outer N/S -> Haslam -> solar IC)
  halo_flat        the same + the halo

Every simulated sky goes through the Stage C selection (fraction-of-rows
rule, `selection.fraction_rule_step`, per component split, with and without
the residual candidate -- rebuilt from the simulated 4.3 GeV pixels for every
row) and the scoring of every row with the selected set (per-bin halo, TS,
quadrant gate, and the joint b-bbar 500 GeV <sigma v>).

THE FIT ENSEMBLE is a stratified 16-model subset of the 2012 grid, in both
splits, plus A2020 (33 rows) -- not the 257 of the data run (15 min per sky).
The two truth models are NOT in it: no fit row is the true model, as none is
for the real sky.

"2 sigma" here is TS > 4 with the sign stated -- NOT sigmav_lo68 > 0, which
is a 1-sigma bound (rank_stage_d's "fraction_positive_2sigma" was that, found
2026-09-22).
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani import dmspectra, ranking, selection, simsky  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-22"
OUT = ROOT / "results" / "fits" / "stage_e"
FRACTION_CUTOFF = 0.8                   # as rank_stage_c
K20 = 5                                 # 20.8 GeV
INJECTION_JSON = ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json"

TRUTHS = {
    "null_structured": {"model": "SO_Z6_R30_T150_C5", "variant": "gas1_ics1", "halo": False,
                        "structures": ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic")},
    "halo_structured": {"model": "SO_Z6_R30_T150_C5", "variant": "gas1_ics1", "halo": True,
                        "structures": ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic")},
    "null_flat": {"model": "SS_Z10_R30_T150_C5", "variant": "gas5_ics3", "halo": False,
                  "structures": ("bubble_flat", "erosita_unif_outer_NS", "loopI_haslam_raw", "sun_moon_ic")},
    "halo_flat": {"model": "SS_Z10_R30_T150_C5", "variant": "gas5_ics3", "halo": True,
                  "structures": ("bubble_flat", "erosita_unif_outer_NS", "loopI_haslam_raw", "sun_moon_ic")},
}
TRUTH_SEED = {name: i + 1 for i, name in enumerate(TRUTHS)}


def fit_codes() -> list:
    """16 models: every (source distribution, z_h) pair once, with R_h, T_S and
    the E(B-V) cut assigned by XORs of their bits so each is 8/8 and balanced
    within every level of the other two factors."""
    srcs, zs = ("SS", "SL", "SY", "SO"), (4, 6, 8, 10)
    out = []
    for si, zi in itertools.product(range(4), range(4)):
        s1, s2, z1, z2 = si & 1, si >> 1, zi & 1, zi >> 1
        r, t, c = (20, 30)[s1 ^ z1], (150, 100000)[s2 ^ z2], (2, 5)[s1 ^ s2 ^ z2]
        out.append(f"{srcs[si]}_Z{zs[zi]}_R{r}_T{t}_C{c}")
    truth_codes = {t["model"] for t in TRUTHS.values()}
    assert not truth_codes & set(out), "a truth model is in the fit ensemble"
    return out


def fit_rows(variant: str) -> list:
    return [(c, variant) for c in fit_codes()] + [("a2020", "-")]


def truth_path(name: str) -> Path:
    return ranking.CACHE / f"stage_e_truth__{name}.npz"


def precompute() -> None:
    t0 = time.time()
    simsky.precompute_shared(EDGES)
    rows = [("a2020", "-")] + [(c, v) for v in ranking.VARIANTS for c in fit_codes()] + \
           sorted({(t["model"], t["variant"]) for t in TRUTHS.values()})
    for i, (c, v) in enumerate(rows, 1):
        simsky.precompute_row(c, v, EDGES)
        print(f"  [{i}/{len(rows)}] pixels {c} {v} ({time.time()-t0:.0f}s)", flush=True)
    amp = simsky.halo_amplitudes_totani(INJECTION_JSON)
    for name, t in TRUTHS.items():
        tr = simsky.build_truth(t["model"], t["variant"], t["structures"], EDGES, amp if t["halo"] else None)
        np.savez(truth_path(name), cells=tr["cells"], pix2=tr["pix2"], halo_cells=tr["halo_cells"],
                 halo_pix2=tr["halo_pix2"], amplitudes=tr["amplitudes"],
                 _meta=json.dumps({"names": tr["names"], **{k: (list(v) if isinstance(v, tuple) else v)
                                                            for k, v in t.items()}}))
        keep = ranking_keep()
        print(f"  truth {name}: {tr['cells'][:, keep].sum(axis=1).round(0).tolist()} counts per bin at |b|>=10; "
              f"halo {tr['halo_cells'][:, keep].sum(axis=1).round(0).tolist()}", flush=True)


def ranking_keep():
    from totani.cells import halo_cell_mask
    return halo_cell_mask()


# ---- one simulated sky --------------------------------------------------------------

_W: dict = {}


def _init():
    _W["sh"] = ranking.shared_cells(EDGES)
    _W["shared"] = simsky.load_shared()
    _W["truths"] = {n: dict(np.load(truth_path(n))) for n in TRUTHS}
    _W["rows"] = {r: simsky.load_row(*r) for v in ranking.VARIANTS for r in fit_rows(v)}


def consensus(rows, shr, cands):
    """The Stage C consensus path on one sky (rank_stage_c step 1, serial)."""
    stacks = {r: selection.base_stack(r[0], r[1], EDGES, shr[r]) for r in rows}
    cur = {r: selection.bg_neg2lnl(stacks[r]) for r in rows}
    selected, trace, remaining = [], [], list(cands)
    while remaining:
        per_row, added = {}, {}
        for r in rows:
            d = {}
            for c in remaining:
                st, v, _ = selection.add_candidate(stacks[r], c, r[0], r[1], EDGES, shr[r])
                d[c] = float((cur[r] - v).sum()); added[(r, c)] = (st, v)
            per_row[r] = d
        rule = selection.fraction_rule_step(list(per_row.values()), remaining, FRACTION_CUTOFF)
        trace.append({"step": len(selected) + 1, "best": rule["best"], "accepted": rule["accepted"],
                      "clears_fraction": rule["clears_fraction"],
                      "best_fraction_among_eligible": rule["best_fraction_among_eligible"],
                      "per_row": {f"{r[0]}__{r[1]}": d for r, d in per_row.items()}})
        if not rule["accepted"]:
            break
        best = rule["best"]
        selected.append(best)
        for r in rows:
            stacks[r], cur[r] = added[(r, best)]
        remaining = [c for c in remaining if c != best and not any(
            best in g and c in g for g in selection.GROUPS.values())]
    return selected, trace, stacks


def score_row(r, st, J) -> dict:
    """`ranking.score` without intervals or the injection gate, and with the
    quadrant fits at 20.8 GeV only (the gate reads nothing else; all nine bins
    was 4x the cost); plus the joint b-bbar 500 GeV fit, best fit and TS only."""
    from totani.cells import halo_cell_mask
    from totani.likelihood import fit_all_bins, fit_all_bins_profile, quadrant_masks, restrict_stack
    keep = halo_cell_mask()
    sr, Jr = restrict_stack(st, keep, J)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    fit = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False)["per_bin"]
    n_bins = st["n_bins"]
    ts = [bg[k]["neg2lnl"] - fit[k]["neg2lnl"] for k in range(n_bins)]
    halo = [fit[k]["signal_amplitude"] * float(Jr[k][keep].sum()) for k in range(n_bins)]
    one = {"counts": st["counts"][K20:K20 + 1], "n_bins": 1,
           "templates": {n: t[K20:K20 + 1] for n, t in st["templates"].items()}}
    quads = {}
    for qn, qm in quadrant_masks().items():
        stq, Jq = restrict_stack(one, qm, J[K20:K20 + 1])
        quads[qn] = fit_all_bins_profile(stq, Jq, halo_mask=qm, with_uncertainty=False)["per_bin"][0]["signal_amplitude"] \
            * float(Jq[0][qm].sum())
    d = {"ts": ts, "halo_counts": halo, "quadrant_counts_20": quads,
         "quadrant_any_negative_20": bool(any(q < 0 for q in quads.values())),
         "halo_north_20": quads["l+b+"] + quads["l-b+"], "halo_south_20": quads["l+b-"] + quads["l-b-"]}
    sv = dmspectra.joint_fit(st, J, "b", 500, EDGES, n_scan=31, with_intervals=False)
    return {"model": r[0], "variant": r[1], **d, "neg2lnl_total": float(sum(b["neg2lnl"] for b in bg)),
            "favor": bool(ts[K20] > 4 and halo[K20] > 0 and not d["quadrant_any_negative_20"]),
            "sigmav_best": sv["sigmav_best"], "sigmav_ts": sv["ts_vs_zero"]}


def run_one(task) -> str:
    name, i = task
    p = OUT / f"{name}__{i:03d}.json"
    if p.exists():
        return f"{p.name} exists"
    t0 = time.time()
    if not _W:
        _init()
    sh, shared = _W["sh"], _W["shared"]
    rng = np.random.default_rng([TRUTH_SEED[name], i, 20260922])
    cells, pix = simsky.draw(_W["truths"][name], sh["counts"], shared, rng)
    sh_sim = dict(sh); sh_sim["counts"] = cells
    out = {"date": DATE, "truth": name, "sim": i, "truth_spec": {k: (list(v) if isinstance(v, tuple) else v)
                                                                 for k, v in TRUTHS[name].items()},
           "counts_per_bin": cells[:, ranking_keep()].sum(axis=1).tolist(), "splits": {}}
    all_c = list(selection.CANDIDATES)
    for variant in ranking.VARIANTS:
        rows = fit_rows(variant)
        shr = {r: {**sh_sim, **simsky.residual_cells(pix, _W["rows"][r], shared, EDGES)} for r in rows}
        out["splits"][variant] = {}
        for branch, cands in (("residual", all_c), ("no_residual", [c for c in all_c if c != "bubble_totani"])):
            selected, trace, stacks = consensus(rows, shr, cands)
            scored = [score_row(r, stacks[r], sh["J"]) for r in rows]
            out["splits"][variant][branch] = {"selected": selected, "trace": trace, "rows": scored}
    out["seconds"] = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp"); tmp.write_text(json.dumps(out)); tmp.rename(p)
    return f"{p.name} {out['seconds']:.0f}s"


# ---- summary -------------------------------------------------------------------------

def summarise() -> None:
    sims = [json.loads(p.read_text()) for p in sorted(OUT.glob("*.json"))]
    summary = {"date": DATE, "fit_codes": fit_codes(), "truths": {k: {**v, "structures": list(v["structures"])}
                                                                    for k, v in TRUTHS.items()}, "cells": {}}
    lines = [f"# Stage E — simulated skies ({DATE})", "",
             f"{len(sims)} simulated skies. Fit ensemble per split: {len(fit_codes())} models + A2020 "
             f"(truth models excluded). Per sky: Stage C (fraction-of-rows, {FRACTION_CUTOFF:.0%}) with and without "
             "the residual construction, then every row scored. '2σ' = TS > 4 at 20.8 GeV with the sign stated; "
             "⟨σv⟩ is b b̄ 500 GeV, joint, signed; 'σv 2σ+' = TS > 4 and ⟨σv⟩ > 0.", "",
             "| truth | split | branch | n | selection (most common, share) | bubble chosen | "
             "rows favor (median, range) | rows 2σ+ @20.8 (median) | median halo @20.8 | median ⟨σv⟩ | rows σv 2σ+ (median) |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for tn in TRUTHS:
        for v in ranking.VARIANTS:
            for br in ("residual", "no_residual"):
                S = [s["splits"][v][br] for s in sims if s["truth"] == tn]
                if not S:
                    continue
                sel = [" → ".join(x["selected"]) for x in S]
                top = max(set(sel), key=sel.count)
                bub = [next((c for c in x["selected"] if c.startswith("bubble_")), "none") for x in S]
                bub_share = {b: bub.count(b) / len(bub) for b in sorted(set(bub))}
                fav = [np.mean([r["favor"] for r in x["rows"]]) for x in S]
                sig = [np.mean([r["ts"][K20] > 4 and r["halo_counts"][K20] > 0 for r in x["rows"]]) for x in S]
                halo = [np.median([r["halo_counts"][K20] for r in x["rows"]]) for x in S]
                sv = [np.median([r["sigmav_best"] for r in x["rows"]]) for x in S]
                svp = [np.mean([r["sigmav_ts"] > 4 and r["sigmav_best"] > 0 for r in x["rows"]]) for x in S]
                key = f"{tn}__{v}__{br}"
                summary["cells"][key] = {"n": len(S), "selections": {k: sel.count(k) for k in set(sel)},
                                         "bubble_share": bub_share, "favor_fraction": fav, "sig2_fraction": sig,
                                         "median_halo_20": halo, "median_sigmav": sv, "sigmav_2sig_pos_fraction": svp}
                lines.append(f"| {tn} | {v} | {br} | {len(S)} | {top} ({sel.count(top)}/{len(S)}) | "
                             + ", ".join(f"{b.replace('bubble_', '')} {100*s:.0f}%" for b, s in bub_share.items())
                             + f" | {np.median(fav):.2f} [{min(fav):.2f}, {max(fav):.2f}] | {np.median(sig):.2f} | "
                             f"{np.median(halo):,.0f} | {np.median(sv):+.2e} | {np.median(svp):.2f} |")
    (ROOT / "results" / "fits" / f"stage_e_summary_{DATE}.json").write_text(json.dumps(summary, indent=2))
    (ROOT / "results" / f"stage_e_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--precompute", action="store_true")
    ap.add_argument("--pilot", action="store_true", help="sim 0 of every truth, serial, timed")
    ap.add_argument("--n-sims", type=int, default=0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--summarise", action="store_true")
    args = ap.parse_args()
    if args.precompute:
        precompute(); return 0
    if args.summarise:
        summarise(); return 0
    tasks = [(n, 0) for n in TRUTHS] if args.pilot else \
        [(n, i) for i in range(args.n_sims) for n in TRUTHS]           # interleaved: partial runs stay balanced
    t0 = time.time()
    if args.pilot:
        for t in tasks:
            print(run_one(t), flush=True)
    else:
        with Pool(args.workers, initializer=_init) as pool:
            for k, msg in enumerate(pool.imap_unordered(run_one, tasks), 1):
                print(f"[{k}/{len(tasks)}] {msg} ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
