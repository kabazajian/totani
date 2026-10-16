#!/usr/bin/env python
"""Do the Loop I shells double-count Oscar's structured map? Fitted shell
counts over the halo cells, background-only, next to each bubble template.

If the structured map (a 1-10 GeV residual template, Ackermann+ 2017 Fig. 12)
already contains Loop I / North Polar Spur emission, the two geometric shells
should fit much LOWER beside it than beside the flat mask. Cell space; one
stack build. 2026-09-18, for the no-Loop-I baseline decision.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
from astropy.io import fits
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from systematics import replace_component  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, restrict_stack  # noqa: E402
from totani.stack import _strip_to_2d, build_background_stack, cell_template_from_shape  # noqa: E402
from totani.templates import bubble_on_target  # noqa: E402
RUN = ROOT / "runs" / "fff2de5bc7e0"; EDGES = TOTANI_EDGES_GEV[:10]

def main():
    t0 = time.time()
    bg = build_background_stack(RUN, EDGES, bubble="macias")
    with fits.open(RUN / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    flat = replace_component(bg, "bubble", cell_template_from_shape(
        RUN, EDGES, (bubble_on_target(target) > 0).astype(float)))
    keep = halo_cell_mask()
    out = {"E_centers_gev": list(TOTANI_CENTERS_GEV[:9]), "counts_over_halo_cells_bg_only": {}}
    for tag, st in (("structured", bg), ("flat", flat)):
        sr = restrict_stack(st, keep)
        fit = fit_all_bins(sr, halo_mask=keep)["per_bin"]
        res = {}
        for name in ("loopI_shell1", "loopI_shell2", "bubble", "isotropic"):
            res[name] = [fit[k]["amplitudes"][name] * float(sr["templates"][name][k][keep].sum())
                         for k in range(9)]
        res["loopI_total"] = [a + b for a, b in zip(res["loopI_shell1"], res["loopI_shell2"])]
        out["counts_over_halo_cells_bg_only"][tag] = res
        print(f"{tag:10s} Loop I total: " + " ".join(f"{v:8.0f}" for v in res["loopI_total"]))
        print(f"{'':10s} bubble      : " + " ".join(f"{v:8.0f}" for v in res["bubble"]))
    s, f = (out["counts_over_halo_cells_bg_only"][k]["loopI_total"] for k in ("structured", "flat"))
    out["loopI_ratio_structured_over_flat"] = [a / b if b else None for a, b in zip(s, f)]
    print("Loop I ratio structured/flat: " + " ".join(f"{r:8.2f}" for r in out["loopI_ratio_structured_over_flat"]))
    p = ROOT / "results" / "fits" / "loopi_amplitudes_2026-09-18.json"
    p.write_text(json.dumps(out, indent=2)); print(f"saved {p} ({time.time()-t0:.0f}s)")

if __name__ == "__main__":
    main()
