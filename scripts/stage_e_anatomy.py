#!/usr/bin/env python
"""Anatomy of the null-sky halo: which ingredient of the residual path makes
Totani's halo appear on a sky that has none.

    python scripts/stage_e_anatomy.py --n-skies 20        # ~30 min, one process

On the `null_structured` skies of Stage E (same seeds, so the same skies),
three diffuse rows -- the TRUTH row itself (SO_Z6_R30_T150_C5, no diffuse
mismatch at all), Totani's baseline SL_Z6_R30_T150_C2, and A2020 -- are each
fitted with five fixed structure sets, and the 20.8 GeV halo, its spectrum,
TS, the quadrant gate and the joint b-bbar 500 GeV <sigma v> are recorded:

  truth      structured + eROSITA unif N/S + Haslam + solar IC   (the sky's own set)
  construct  residual pair in place of the structured map, else the truth set
  consensus  residual pair + shells+NPS + eROSITA brightness      (Stage E's residual selection)
  flat_truth flat mask in place of the structured map, else the truth set
  flat_fig1  flat mask + two shells + eROSITA unif N/S            (Figure 1's flat panel)

The truth row under the truth set measures Poisson noise alone; the same row
under `construct` isolates the residual construction; the other rows add the
diffuse mismatch; `consensus` adds the companion mismatch; the two flat sets
give the flat mask's bias on a structured-bubble sky.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from totani import ranking, selection, simsky  # noqa: E402

EDGES = stage_e.EDGES
DATE = "2026-09-23"
CONFIG = {
    "null_structured": {
        "rows": [("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas1_ics1"), ("a2020", "-")],
        "sets": {
            "truth": ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
            "construct": ("bubble_totani", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
            "consensus": ("bubble_totani", "loopI_shells+nps", "erosita_bright_NS"),
            "flat_truth": ("bubble_flat", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
            "flat_fig1": ("bubble_flat", "loopI_shells", "erosita_unif_NS")}},
    # the mirror image: a FLAT-bubble sky; the other 5-ring row is from the Stage E subset
    "null_flat": {
        "rows": [("SS_Z10_R30_T150_C5", "gas5_ics3"), ("SL_Z6_R20_T150_C5", "gas5_ics3"), ("a2020", "-")],
        "sets": {
            "truth": ("bubble_flat", "erosita_unif_outer_NS", "loopI_haslam_raw", "sun_moon_ic"),
            "construct": ("bubble_totani", "erosita_unif_outer_NS", "loopI_haslam_raw", "sun_moon_ic"),
            "consensus": ("bubble_totani", "loopI_shells+nps", "sun_moon_ic"),
            "structured_truth": ("bubble_structured", "erosita_unif_outer_NS", "loopI_haslam_raw", "sun_moon_ic"),
            "structured_fig1": ("bubble_structured", "loopI_shells", "erosita_unif_NS")}},
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="null_structured")
    ap.add_argument("--n-skies", type=int, default=20)
    args = ap.parse_args()
    ROWS, SETS = CONFIG[args.truth]["rows"], CONFIG[args.truth]["sets"]
    t0 = time.time()
    sh = ranking.shared_cells(EDGES)
    shared = simsky.load_shared()
    truth = dict(np.load(stage_e.truth_path(args.truth)))
    rows_pix = {r: simsky.load_row(*r) for r in ROWS}
    out = {"date": DATE, "truth": args.truth, "sets": {k: list(v) for k, v in SETS.items()}, "skies": []}
    for i in range(args.n_skies):
        rng = np.random.default_rng([stage_e.TRUTH_SEED[args.truth], i, 20260922])
        cells, pix = simsky.draw(truth, sh["counts"], shared, rng)
        sh_sim = dict(sh); sh_sim["counts"] = cells
        sky = {"sim": i, "rows": {}}
        for r in ROWS:
            shr = {**sh_sim, **simsky.residual_cells(pix, rows_pix[r], shared, EDGES)}
            st0 = selection.base_stack(r[0], r[1], EDGES, shr)
            res = {}
            for sname, cands in SETS.items():
                st = st0
                for c in cands:
                    st, _, _ = selection.add_candidate(st, c, r[0], r[1], EDGES, shr)
                res[sname] = stage_e.score_row(r, st, sh["J"])
            sky["rows"][f"{r[0]}__{r[1]}"] = res
        out["skies"].append(sky)
        print(f"[{i + 1}/{args.n_skies}] " + "  ".join(
            f"{s}:{np.median([sky['rows'][k][s]['halo_counts'][stage_e.K20] for k in sky['rows']]):6.0f}" for s in SETS)
            + f"  ({time.time() - t0:.0f}s)", flush=True)
    p = ROOT / "results" / "fits" / f"stage_e_anatomy_{DATE}_{args.truth}.json"
    p.write_text(json.dumps(out))
    # summary
    lines = [f"# Null-sky halo anatomy — {args.truth}, {args.n_skies} skies ({DATE})", "",
             "Median over skies of the 20.8 GeV halo (counts), TS, favor fraction and joint b b̄ 500 GeV ⟨σv⟩, per row and structure set.", "",
             "| row | set | halo @20.8 (16–84%) | TS | favor | ⟨σv⟩ | σv TS |", "|---|---|---|---|---|---|---|"]
    for k in out["skies"][0]["rows"]:
        for s in SETS:
            h = np.array([sk["rows"][k][s]["halo_counts"][stage_e.K20] for sk in out["skies"]])
            ts = np.median([sk["rows"][k][s]["ts"][stage_e.K20] for sk in out["skies"]])
            fav = np.mean([sk["rows"][k][s]["favor"] for sk in out["skies"]])
            sv = np.median([sk["rows"][k][s]["sigmav_best"] for sk in out["skies"]])
            svt = np.median([sk["rows"][k][s]["sigmav_ts"] for sk in out["skies"]])
            lines.append(f"| {k} | {s} | {np.median(h):,.0f} ({np.percentile(h, 16):,.0f}–{np.percentile(h, 84):,.0f}) | {ts:.1f} | {fav:.2f} | {sv:+.2e} | {svt:.1f} |")
    lines += ["", "Median halo spectrum (counts per bin) for the truth row:", "",
              "| set | " + " | ".join(f"{e:g}" for e in [1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100]) + " |", "|---|" + "---|" * 9]
    k0 = f"{ROWS[0][0]}__{ROWS[0][1]}"
    for s in SETS:
        med = np.median([sk["rows"][k0][s]["halo_counts"] for sk in out["skies"]], axis=0)
        lines.append(f"| {s} | " + " | ".join(f"{m:,.0f}" for m in med) + " |")
    (ROOT / "results" / f"stage_e_anatomy_{DATE}_{args.truth}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
