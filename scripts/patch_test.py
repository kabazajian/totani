#!/usr/bin/env python
"""Oscar's matched-angular-distance test on the Zhang+ (2408.06312) patches.

    python scripts/patch_test.py

🚨 THE NULL IS A SYMMETRY, NOT A MODEL. Oscar's proposal (2026-09-05) is that
"a spherical halo should give the same intensity at the same angular distance;
a lobe or shell need not." Reflecting a patch in l -> -l or in b -> -b leaves
theta_GC EXACTLY invariant, so each reflected pair is a matched-theta
comparison that needs no bubble boundary, no digitised curve and no
assumption about where the eROSITA edge runs. That is the whole reason this
test can run before Oscar's templates arrive.

THE PATCHES ARE THE ONLY ONES ZHANG+ SPECIFY NUMERICALLY (their sec 3.1):

    R1  20 < l < 40,   45 < b < 60    northern "outer region"
    R3  20 < l < 40,  -45 < b < -30   southern "outer region"

("outer region" = outside the Fermi Bubbles, inside the eROSITA Bubbles.)
Their comparison patches R2/R4, the background R0, the Fermi-cap R5 and both
bubble boundary curves are given only as figure elements -- no coordinates, no
released region files -- which is why this module reflects the two specified
boxes instead of trying to reproduce their R1-R2 / R3-R4 differences.

⚠ THE PATCHES DO NOT LAND ON THE 10-DEGREE CELL GRID. Their longitudes (20,
40) are cell edges; their latitudes (45, -45, 60, -30) cut the b in [40,50]
and b in [-50,-40] cells in half. So this test CANNOT be done on the cell
arrays the fit uses -- it has to sum pixels. What it reuses from the fit is
the fitted normalisations, which are one scalar per component per energy bin
and are therefore valid at any resolution. `mu_patch = sum_l f_l T_l(patch)`
is exact given those f_l; nothing is re-fitted here.

MEMORY: each pixel-space template is reduced to its patch sums and freed
before the next is built, because the machine is shared and the naive version
(hold all 19 in pixel space) peaks near 7 GB. That step now lives in
`totani.patches.region_sums`, shared with `scripts/boundary_test.py`; the
boxes here are the 0/1 special case of its weight maps. The refactor was
verified 2026-09-08 by re-running this script and diffing the JSON against the
committed `patch_test_2026-09-06.json`.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402
from totani.stack import build_background_stack, halo_cell_template  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]

#: The two Zhang+ boxes and all four reflections of each. Reflection in l or b
#: preserves theta_GC exactly, so within a family every patch samples the same
#: angular distances from the Galactic Centre.
FAMILIES = {
    "R1": {"R1":        (20.0, 40.0, 45.0, 60.0),     # Zhang+ north, in outer region
           "R1_lmir":  (-40.0, -20.0, 45.0, 60.0),    # mirrored in longitude
           "R1_bmir":   (20.0, 40.0, -60.0, -45.0),   # mirrored in latitude
           "R1_lbmir": (-40.0, -20.0, -60.0, -45.0)},
    "R3": {"R3":        (20.0, 40.0, -45.0, -30.0),   # Zhang+ south, in outer region
           "R3_lmir":  (-40.0, -20.0, -45.0, -30.0),
           "R3_bmir":   (20.0, 40.0, 30.0, 45.0),
           "R3_lbmir": (-40.0, -20.0, 30.0, 45.0)},
}

#: A box that IS an exact union of cells, used to prove the pixel-sum path here
#: agrees with the cell-space stack the fit actually uses.
VALIDATION_BOX = (20.0, 40.0, 30.0, 40.0)


def _box_mask(L, B, box):
    l0, l1, b0, b1 = box
    return (L >= l0) & (L < l1) & (B >= b0) & (B < b1)


def patch_sums(run_dir, edges_gev, boxes) -> dict:
    """`region_sums` on the boxes, which are the 0/1 case of its weight maps.

    Point-source membership stays the explicit box test rather than the
    nearest-pixel lookup `region_sums` defaults to: a source at l = 19.99 sits
    in the pixel centred on 20.0625, so the two rules disagree for sources
    within half a pixel of an edge, and the box test is the one every committed
    number here was produced with.
    """
    from totani.patches import region_sums
    from totani.stack import _strip_to_2d
    from totani.templates import grid_lb

    with fits.open(Path(run_dir) / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    L, B = grid_lb(target)
    weights = {n: _box_mask(L, B, b).astype(float) for n, b in boxes.items()}

    def source_weight(name, sl, sb):
        l0, l1, b0, b1 = boxes[name]
        return ((sl >= l0) & (sl < l1) & (sb >= b0) & (sb < b1)).astype(float)

    return region_sums(run_dir, edges_gev, weights, source_weight)


def main() -> int:
    import resource

    t0 = time.time()
    boxes = {n: b for fam in FAMILIES.values() for n, b in fam.items()}
    boxes["_validation"] = VALIDATION_BOX

    ps = patch_sums(RUN, EDGES, boxes)
    print(f"patch sums built ({time.time()-t0:.0f}s)", flush=True)

    bg = build_background_stack(RUN, EDGES)
    J = halo_cell_template(RUN, EDGES)
    names = sorted(bg["templates"])
    print(f"cell stack built ({time.time()-t0:.0f}s)", flush=True)

    # ---- validation: the pixel-sum path must reproduce the cell-space stack --
    # VALIDATION_BOX is exactly cells (row b=35) x (cols l=25, l=35).
    cell_sum = {n: [float(bg["templates"][n][k, 9, 8:10].sum()) for k in range(9)]
                for n in names}
    worst = 0.0
    for n in names:
        a = np.array(ps["templates"]["_validation"][n])
        b = np.array(cell_sum[n])
        good = b != 0
        if good.any():
            worst = max(worst, float(np.abs(a[good] / b[good] - 1).max()))
    cw = np.array(ps["counts"]["_validation"])
    cc = np.array([float(bg["counts"][k, 9, 8:10].sum()) for k in range(9)])
    print(f"VALIDATION worst template mismatch vs cell stack: {worst:.2e}; "
          f"counts {cw.sum():.0f} vs {cc.sum():.0f}")

    # ---- the fits whose normalisations we evaluate on the patches -----------
    halo_keep = halo_cell_mask()
    sr, Jr = restrict_stack(bg, halo_keep, J)
    bg_only = fit_all_bins(sr, halo_mask=halo_keep)
    with_halo = fit_all_bins_profile(sr, Jr, halo_mask=halo_keep,
                                     with_uncertainty=False)
    print(f"fits done ({time.time()-t0:.0f}s)", flush=True)

    results = {"date": "2026-09-06", "run": "fff2de5bc7e0",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]),
               "boxes": {n: list(b) for n, b in boxes.items()},
               "validation_worst_rel_mismatch": worst,
               "patches": {}}

    for pn in boxes:
        if pn == "_validation":
            continue
        n_obs = np.array(ps["counts"][pn])
        mu_bg = np.zeros(9)
        mu_halo = np.zeros(9)
        for k in range(9):
            for nm in names:
                t = ps["templates"][pn][nm][k]
                mu_bg[k] += bg_only["per_bin"][k]["amplitudes"][nm] * t
                mu_halo[k] += with_halo["per_bin"][k]["bg_amplitudes"][nm] * t
        results["patches"][pn] = {
            "observed": n_obs.tolist(),
            "mu_background_only": mu_bg.tolist(),
            "mu_background_with_halo_fitted": mu_halo.tolist(),
            "sigma_background_only": ((n_obs - mu_bg) / np.sqrt(mu_bg)).tolist(),
        }

    # ---- the reflection tests -------------------------------------------
    # Each pair is a theta_GC-preserving reflection, so a spherical component
    # can only produce a difference through exposure/masking, which the
    # `halo_diff` term removes. What is left is a limit on a lobe.
    roi_J = [float(Jr[k][halo_keep].sum()) for k in range(9)]
    halo_counts = [with_halo["per_bin"][k]["signal_amplitude"] * roi_J[k]
                   for k in range(9)]
    results["fitted_halo_counts_roi"] = halo_counts
    results["reflections"] = {}
    for tag, (a, b) in (("R1_vs_l_mirror", ("R1", "R1_lmir")),
                        ("R3_vs_l_mirror", ("R3", "R3_lmir")),
                        ("R1_vs_b_mirror", ("R1", "R1_bmir")),
                        ("R3_vs_b_mirror", ("R3", "R3_bmir"))):
        rows = []
        for k in range(9):
            na = ps["counts"][a][k]; ma = results["patches"][a]["mu_background_only"][k]
            nb = ps["counts"][b][k]; mb = results["patches"][b]["mu_background_only"][k]
            excess = (na - ma) - (nb - mb)
            sd = float(np.sqrt(ma + mb))
            hA = halo_counts[k] * ps["halo"][a][k] / roi_J[k]
            hB = halo_counts[k] * ps["halo"][b][k] / roi_J[k]
            lobe = excess - (hA - hB)
            rows.append({"bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
                         "excess": excess, "sigma_excess": sd,
                         "halo_predicted_diff": hA - hB,
                         "lobe": lobe, "lobe_sigma": lobe / sd,
                         "halo_counts_in_A": hA,
                         "lobe_over_halo_95pct_UL": (lobe + 2 * sd) / hA if hA > 0 else None})
        results["reflections"][tag] = rows
        b5 = rows[5]
        print(f"  {tag:16s} bin5: lobe {b5['lobe']:+6.0f} +- {b5['sigma_excess']:.0f} "
              f"({b5['lobe_sigma']:+.2f}s)  95% UL / halo = "
              f"{b5['lobe_over_halo_95pct_UL']:+.2f}")

    out = ROOT / "results" / "fits" / "patch_test_2026-09-06.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"saved {out}")
    print("peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
