#!/usr/bin/env python
"""Oscar's items 4c and 9: (A) the halo fitted OUTSIDE the geometric Fermi-bubble
footprint versus inside, each with its own injection check; (B) the eROSITA
brightness-vs-uniform comparison and the halo restricted to the native
western half of the eROSITA map (Galactic West, l < 0) versus the
reconstructed eastern half.

    python scripts/spatial_split_tests.py

Cell space, 10 deg. 'Outside' cells are those in which the flat bubble
template contributes < 2% of its per-cell maximum in the 20.8 GeV bin; at
this resolution that is a conservative footprint (the PSF is < 0.5 deg here).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking, selection  # noqa: E402
from totani.cells import cell_centers, halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, injection_slope, restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-24"
K20 = 5
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SL_Z10_R30_T150_C5", "gas1_ics1"),
        ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-"), ("G21_SA0_R12+iso", "gas1_ics1")]
ROWS = [r for r in ROWS if ranking.available(*r)]          # the Galp21 row only when its cubes are ingested
PANELS = {"residual": ("bubble_totani", "loopI_shells", "erosita_unif_NS"),
          "flat": ("bubble_flat", "loopI_shells", "erosita_unif_NS"),
          "structured": ("bubble_structured", "loopI_shells", "erosita_unif_NS")}


def halo_on(st, J, mask):
    sr, Jr = restrict_stack(st, mask, J)
    bg = fit_all_bins(sr, halo_mask=mask)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=mask, with_uncertainty=True)["per_bin"]
    roi = float(Jr[K20][mask].sum())
    slope = injection_slope(sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=mask)["slopes"][K20]
    return {"halo20": pf[K20]["signal_amplitude"] * roi, "lo": pf[K20]["signal_amplitude_lo"] * roi, "hi": pf[K20]["signal_amplitude_hi"] * roi,
            "ts20": bg[K20]["neg2lnl"] - pf[K20]["neg2lnl"], "n_cells": int(mask.sum()), "injection_slope20": slope,
            "halo_bins": [pf[k]["signal_amplitude"] * float(Jr[k][mask].sum()) for k in range(9)]}


def main() -> int:
    sh = ranking.shared_cells(EDGES); J = sh["J"]; keep = halo_cell_mask(); Lc, _ = cell_centers()
    fb = sh["bubble_flat"][K20]; outside = keep & (fb < 0.02 * fb.max()); inside = keep & ~outside
    west = keep & (Lc < 0); east = keep & (Lc > 0)
    out = {"date": DATE, "n_cells": {"outside": int(outside.sum()), "inside": int(inside.sum())}, "A": [], "B": []}
    lines = [f"# Spatial split tests ({DATE})", "", f"## A. Halo inside vs outside the flat-bubble footprint ({int(inside.sum())} / {int(outside.sum())} of 120 cells)", "",
             "| row | panel | region | halo@20.8 [68%] | TS@20.8 | injection slope @20.8 |", "|---|---|---|---|---|---|"]
    for row in ROWS:
        for pname, cands in PANELS.items():
            try:
                st = selection.base_stack(row[0], row[1], EDGES, sh)
                for c in cands:
                    st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
            except (FileNotFoundError, KeyError):
                continue
            for reg, m in (("all", keep), ("inside", inside), ("outside", outside)):
                r = halo_on(st, J, m); out["A"].append({"model": row[0], "variant": row[1], "panel": pname, "region": reg, **r})
                lines.append(f"| {row[0]} {row[1]} | {pname} | {reg} | {r['halo20']:,.0f} [{r['lo']:,.0f}, {r['hi']:,.0f}] | {r['ts20']:.1f} | {r['injection_slope20']:.2f} |")
                print(lines[-1], flush=True)
    lines += ["", "## B. eROSITA: native West (l < 0) vs reconstructed East (l > 0)", "",
              "Δ(−2lnL) = uniform − brightness (H0, 9 bins; positive favors the brightness map), and the halo under each bubble template with the brightness lobes, per half.", "",
              "| row | half | Δ(−2lnL) unif − bright | halo@20.8, residual | flat | structured |", "|---|---|---|---|---|---|"]
    for row in ROWS:
        if row[1] == "gas5_ics3":
            continue
        for half, m in (("west", west), ("east", east)):
            rec = {"model": row[0], "variant": row[1], "half": half}
            vals = {}
            for e in ("erosita_unif_NS", "erosita_bright_NS"):
                st = selection.base_stack(row[0], row[1], EDGES, sh)
                for c in ("bubble_structured", "loopI_shells", e):
                    st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
                vals[e] = float(sum(b["neg2lnl"] for b in fit_all_bins(restrict_stack(st, m), halo_mask=m)["per_bin"]))
            rec["delta_unif_minus_bright"] = vals["erosita_unif_NS"] - vals["erosita_bright_NS"]
            for pname, cands in PANELS.items():
                st = selection.base_stack(row[0], row[1], EDGES, sh)
                for c in (cands[0], "loopI_shells", "erosita_bright_NS"):
                    st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
                rec[pname] = halo_on(st, J, m)
            out["B"].append(rec)
            lines.append(f"| {row[0]} {row[1]} | {half} | {rec['delta_unif_minus_bright']:+.0f} | " + " | ".join(
                f"{rec[p]['halo20']:,.0f} (TS {rec[p]['ts20']:.0f})" for p in PANELS) + " |")
            print(lines[-1], flush=True)
    (ROOT / "results" / "fits" / f"spatial_split_tests_{DATE}.json").write_text(json.dumps(out, indent=1))
    (ROOT / "results" / f"spatial_split_tests_{DATE}.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
