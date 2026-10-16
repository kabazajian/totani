#!/usr/bin/env python
"""Structure selection and scoring at any cell size, over all nine bins or
over the response-valid bins only (Oscar, 2026-09-29, item 1: 'a
restricted-energy check should restrict both inference and selection').
Generalizes stage_2deg.py (kept as the record of the 09-26/27 runs).

    FINE_TAG=psf python scripts/stage_select.py --cell 2 --bins resp --workers 5
    FINE_TAG=psf python scripts/stage_select.py --cell 1 --bins all  --workers 5
    FINE_TAG=psf python scripts/stage_select.py --cell 1 --bins resp --workers 5

'resp' = bins 5.6-130 GeV (PSF R68 <= 1/4 of a 1-deg cell), fixed before any
result. With a bin subset, every step uses only those bins: the per-candidate
Delta(-2lnL) sums, the BIC threshold p ln n with p = K x (templates) and
n = (fitted cells) x K, the own greedy paths, and the scoring (per-bin TS,
Q^+ over the K bins with the Macias et al. 2018 mixture, the 20.8 GeV halo
with its 68% interval, the quadrant gate, and the joint b-bbar 500 GeV
<sigma v> fitted over the same bins). The fraction-of-rows rule (80%) and the
candidate list are those of stage_2deg. Voters: the 16 stratified 2012 models
of the decomposition + A2020; Galp21 rows present at this cell size are
scored and given own paths but do not vote (their input checks are pending).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine, selection, stats  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-29"
CUTOFF = 0.8
BINSETS = {"all": list(range(9)), "resp": [3, 4, 5, 6, 7, 8]}
G21 = [f"G21_{m}{ic}" for m in ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12") for ic in ("", "+iso")]
_W: dict = {}


def _init(cell, bins):
    _W["cell"] = cell; _W["bins"] = bins; _W["lib"] = fine.load_library(cell)


def _row(r):
    if r not in _W:
        _W[r] = fine.load_row(_W["cell"], *r)
    return _W[r]


def _sub(st):
    """The stack restricted to the chosen bins."""
    b = _W["bins"]
    return {"counts": st["counts"][b], "n_bins": len(b), "templates": {n: t[b] for n, t in st["templates"].items()}}


def _bg(st):
    return fine.bg_neg2lnl(_sub(st), _W["lib"]["keep"])


def _add(st, cand, row):
    """fine.add_candidate with the sign choice made on the chosen bins only."""
    lib = _W["lib"]; keep = lib["keep"]
    new = {n: fine.template(lib, row, n) for n in selection.CANDIDATES[cand]}
    free = [n for n in new if n in selection.SIGN_FREE]
    if not free:
        s2 = {**st, "templates": {**st["templates"], **new}}
        return s2, _bg(s2)
    nf = free[0]; fixed = {n: v for n, v in new.items() if n != nf}
    plus = {**st, "templates": {**st["templates"], **fixed, nf: new[nf]}}
    minus = {**st, "templates": {**st["templates"], **fixed, nf: -new[nf]}}
    vp, vm = _bg(plus), _bg(minus)
    signs_sub = np.where(vp <= vm, 1, -1)
    signs = np.ones(9); signs[_W["bins"]] = signs_sub
    s2 = {**st, "templates": {**st["templates"], **fixed, nf: new[nf] * signs[:, None, None]}}
    return s2, np.minimum(vp, vm)


def _stack_with(r, selected):
    lib = _W["lib"]; row = _row(r); st = fine.base_stack(lib, row)
    for c in selected:
        st, _ = _add(st, c, row)
    return st


def _deltas(args):
    r, selected, cands = args
    row = _row(r); st = _stack_with(r, selected); cur = _bg(st)
    return r, {c: float((cur - _add(st, c, row)[1]).sum()) for c in cands}


def _thr(c):
    K = len(_W["bins"]); return fine.bic_threshold(K * len(selection.CANDIDATES[c]), int(_W["lib"]["keep"].sum()) * K)


def _own(args):
    r, cands = args
    row = _row(r); st = fine.base_stack(_W["lib"], row); cur = _bg(st)
    chosen, steps, rem = [], [], list(cands)
    while rem:
        d = {c: float((cur - _add(st, c, row)[1]).sum()) for c in rem}
        best = max(d, key=d.get); thr = _thr(best)
        steps.append({"best": best, "delta": d[best], "threshold": thr, "accepted": d[best] > thr})
        if d[best] <= thr:
            break
        st, cur = _add(st, best, row); chosen.append(best)
        rem = [c for c in rem if c != best and not any(best in g and c in g for g in selection.GROUPS.values())]
    return {"model": r[0], "variant": r[1], "chosen": chosen, "steps": steps}


def _score(args):
    r, selected = args
    from totani import dmspectra
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]; b = _W["bins"]
    st = _stack_with(r, selected); sub = _sub(st); Jb = J[b]
    sr, Jr = restrict_stack(sub, keep, Jb)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=True, n_refine=6)["per_bin"]
    ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(len(b))]; a = [pf[k]["signal_amplitude"] for k in range(len(b))]
    k20 = b.index(5); roi20 = float(Jr[k20][keep].sum())
    qa = {}
    for q, m in lib["quads"].items():
        stq, Jq = restrict_stack(st, m, J)
        one = {"counts": stq["counts"][5:6], "n_bins": 1, "templates": {n: t[5:6] for n, t in stq["templates"].items()}}
        qa[q] = fit_all_bins_profile(one, Jq[5:6], halo_mask=m, with_uncertainty=False)["per_bin"][0]["signal_amplitude"] * float(Jq[5][m].sum())
    edges = EDGES[b[0]:b[-1] + 2]
    sv = dmspectra.joint_fit(sub, Jb, "b", 500, edges, halo_mask=keep, n_scan=31)
    return {"model": r[0], "variant": r[1], "structures": list(selected), "ts": ts, "a": a, "Q": stats.summarize(ts, a),
            "halo20": a[k20] * roi20, "halo20_lo": pf[k20]["signal_amplitude_lo"] * roi20, "halo20_hi": pf[k20]["signal_amplitude_hi"] * roi20,
            "quadrants20": qa, "favor": bool(ts[k20] > 4 and a[k20] > 0 and all(v > 0 for v in qa.values())),
            "neg2lnl_H0": float(sum(x["neg2lnl"] for x in bg)),
            "sigmav": {k: sv[k] for k in ("sigmav_best", "sigmav_lo68", "sigmav_hi68", "ts_vs_zero")}}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, required=True); ap.add_argument("--bins", choices=BINSETS, default="all")
    ap.add_argument("--workers", type=int, default=5); ap.add_argument("--render", action="store_true"); a = ap.parse_args()
    bins = BINSETS[a.bins]; tag = f"cells{a.cell:g}{fine.tag()}_{a.bins}"
    js = ROOT / "results" / "fits" / f"stage_select_{tag}_{DATE}.json"
    if a.render:
        render(json.loads(js.read_text()), tag); return 0
    t0 = time.time(); out = {"date": DATE, "cell_deg": a.cell, "bins": bins, "psf": fine.tag(), "cutoff": CUTOFF, "splits": {}}
    g21 = [g for g in G21 if fine.row_path(a.cell, g, "gas1_ics1").exists()]
    with Pool(a.workers, initializer=_init, initargs=(a.cell, bins)) as pool:
        for v in ("gas1_ics1", "gas5_ics3"):
            voters = [(c, v) for c in stage_e.fit_codes()] + [("a2020", "-")]
            extra = [(g, v) for g in g21 if fine.row_path(a.cell, g, v).exists()]
            for branch in ("residual", "no_residual"):
                cands = [c for c in selection.CANDIDATES if branch == "residual" or c != "bubble_totani"]
                selected, trace, rem = [], [], list(cands)
                _init(a.cell, bins)
                while rem:
                    rows = dict(pool.map(_deltas, [(r, tuple(selected), tuple(rem)) for r in voters]))
                    phys = list(rows.values())
                    thr = {c: _thr(c) for c in rem}
                    clears = {c: float(np.mean([p[c] > thr[c] for p in phys])) for c in rem}
                    elig = [c for c in rem if clears[c] >= CUTOFF]
                    if elig:
                        tops = [max(elig, key=lambda c: p[c]) for p in phys]
                        bf = {c: float(np.mean([t == c for t in tops])) for c in elig}; best, acc = max(bf, key=bf.get), True
                    else:
                        bf = {}; best, acc = max(clears, key=clears.get), False
                    trace.append({"best": best, "accepted": acc, "clears": clears, "best_fraction": bf,
                                  "per_row": {f"{r[0]}__{r[1]}": d for r, d in rows.items()}})
                    print(f"{v} {branch} step {len(trace)}: {'ACCEPT' if acc else 'stop'} {best} (clears {100*clears[best]:.0f}%, "
                          f"best {100*bf.get(best, 0):.0f}%) ({time.time()-t0:.0f}s)", flush=True)
                    if not acc:
                        break
                    selected.append(best)
                    rem = [c for c in rem if c != best and not any(best in g and c in g for g in selection.GROUPS.values())]
                own = pool.map(_own, [(r, tuple(cands)) for r in voters + extra])
                scored = pool.map(_score, [(r, tuple(selected)) for r in voters + extra])
                out["splits"].setdefault(v, {})[branch] = {"selected": selected, "trace": trace, "own_paths": own, "scored": scored}
                js.write_text(json.dumps(out, indent=1, default=float))
                print(f"{v} {branch}: selected {selected} ({time.time()-t0:.0f}s)", flush=True)
    render(out, tag)
    return 0


def render(out, tag):
    K = len(out["bins"])
    lines = [f"# Structure selection and scoring — {tag} ({DATE})", "",
             f"Bins used for selection AND inference: {out['bins']} (K = {K}). Fraction-of-rows rule, cutoff {out['cutoff']}; BIC p ln n with n = cells × K. "
             "Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, "
             "all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).", ""]
    for v, br in out["splits"].items():
        for b, d in br.items():
            sc = d["scored"]; s12 = [x for x in sc if not x["model"].startswith("G21")]; sg = [x for x in sc if x["model"].startswith("G21")]
            def summ(S):
                if not S:
                    return "—"
                h = [x["halo20"] for x in S]; sv = [x["sigmav"]["sigmav_best"] * 1e25 for x in S]; z = [x["Q"]["Z_2018"] for x in S]
                return (f"favor {sum(x['favor'] for x in S)}/{len(S)}; halo20 median {np.median(h):,.0f} [{min(h):,.0f}–{max(h):,.0f}]; "
                        f"Z₂₀₁₈ median {np.median(z):.1f} [{min(z):.1f}–{max(z):.1f}]; ⟨σv⟩ median {np.median(sv):+.1f}, 2σ+ "
                        f"{sum(x['sigmav']['ts_vs_zero'] > 4 and x['sigmav']['sigmav_best'] > 0 for x in S)}/{len(S)}")
            own_b = {}
            for o in d["own_paths"]:
                bub = next((c for c in o["chosen"] if c.startswith("bubble_")), "none")
                key = "G21" if o["model"].startswith("G21") else "2012"
                own_b.setdefault(key, {}); own_b[key][bub] = own_b[key].get(bub, 0) + 1
            lines += [f"## {v}, {b.replace('_', ' ')}", "", f"Selected: **{' → '.join(d['selected']) or 'nothing'}**", "",
                      f"- 2012 subset: {summ(s12)}", f"- Galp21: {summ(sg)}", f"- own paths, bubble chosen: {own_b}", ""]
            t1 = d["trace"][0]
            bub = {c: t1["clears"].get(c) for c in ("bubble_flat", "bubble_structured", "bubble_totani") if c in t1["clears"]}
            per = t1["per_row"]
            fvs = [p["bubble_flat"] - p["bubble_structured"] for p in per.values() if "bubble_flat" in p and "bubble_structured" in p]
            lines += [f"- step 1: clears {bub}; flat − structured improvement per row median {np.median(fvs):+.0f} "
                      f"(flat larger in {sum(x > 0 for x in fvs)}/{len(fvs)} rows)", ""]
    (ROOT / "results" / f"stage_select_{tag}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
