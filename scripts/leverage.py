#!/usr/bin/env python
"""Where the halo's constraining power actually lives, in latitude.

    python scripts/leverage.py

🚨 THE QUESTION THIS ANSWERS IS YOUYOU'S, NOT OURS. She proposed (2026-09-05)
forecasting what the |b| <= 60 Galp21 extension buys by measuring the spread
across the six models inside 10 < |b| < 30 -- the band her current cubes
already cover -- before anyone spends compute on the extension. That is the
right instinct, and this script measures how much of the answer such a
forecast can see.

The finding is that it can see less than half of it, and structurally so: the
fit uses ONE set of background normalizations for the whole ROI, so a
background combination that imitates the halo at 10-30 deg must imitate it at
30-60 deg too, with the same numbers. That shared-coefficient constraint --
not the halo's brightness at any latitude -- supplies the majority of the
orthogonal power. A forecast confined to |b| <= 30 cannot measure it, because
the contrast it rests on is precisely the part the current maps are missing.

Reported for every energy bin and under two weightings, because a claim this
structural should not depend on either.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani.cells import cell_centers, halo_cell_mask  # noqa: E402
from totani.degeneracy import lever_arm_share, leverage  # noqa: E402
from totani.likelihood import (  # noqa: E402
    fit_all_bins_profile,
    injection_slope,
    restrict_stack,
)
from totani.stack import build_background_stack, halo_cell_template  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]

#: The overlap band is set by Youyou's current cubes (|b| <= 30) against
#: Totani's ROI floor (|b| >= 10), so it is a fact about the delivery, not a
#: tuning choice. The outer region is the rest of the halo mask.
BAND_MAX_DEG = 30.0


def main() -> int:
    import resource

    t0 = time.time()
    bg = build_background_stack(RUN, EDGES)
    J = halo_cell_template(RUN, EDGES)
    names = sorted(bg["templates"])
    T = np.stack([bg["templates"][n] for n in names])
    print(f"stack built ({time.time()-t0:.0f}s)", flush=True)

    _, Bc = cell_centers()
    full = halo_cell_mask()
    band = full & (np.abs(Bc) < BAND_MAX_DEG)
    outer = full & (np.abs(Bc) >= BAND_MAX_DEG)
    if not np.array_equal(band | outer, full) or (band & outer).any():
        raise ValueError("band and outer do not partition the halo cells")

    results = {"date": "2026-09-06", "run": "fff2de5bc7e0",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]),
               "band_max_deg": BAND_MAX_DEG,
               "n_cells": {"full": int(full.sum()), "band": int(band.sum()),
                           "outer": int(outer.sum())},
               "component_order": names, "per_bin": []}

    for k in range(9):
        # `full` orders cells the same way for h and A, which is all that is
        # needed -- the cell ordering itself is arbitrary and cancels.
        h = J[k][full]
        A = T[:, k][:, full].T
        sub_band = band[full]
        sub_outer = outer[full]
        counts_full = bg["counts"][k][full]
        w = 1.0 / np.sqrt(np.maximum(counts_full, 1.0))

        row = {"bin": k, "E_gev": TOTANI_CENTERS_GEV[k]}
        for tag, weights in (("unweighted", None), ("poisson", w)):
            la = lever_arm_share(h, A, [sub_band, sub_outer], weights)
            row[tag] = {
                "leverage_full": la["leverage_common"],
                "leverage_band": la["leverage_per_region"][0],
                "leverage_outer": la["leverage_per_region"][1],
                "leverage_if_fitted_separately": la["leverage_separate"],
                "lever_arm_share": la["share"],
            }
        row["halo_template_fraction_in_band"] = float(
            J[k][band].sum() / J[k][full].sum())
        row["counts_fraction_in_band"] = float(
            bg["counts"][k][band].sum() / bg["counts"][k][full].sum())
        results["per_bin"].append(row)

    print(f"{'bin':>3}{'E':>7}{'lev_full':>10}{'lev_band':>10}{'lev_outer':>11}"
          f"{'separate':>10}{'share':>8}{'J in band':>11}")
    for row in results["per_bin"]:
        u = row["unweighted"]
        print(f"{row['bin']:>3}{row['E_gev']:7.1f}{u['leverage_full']:10.4f}"
              f"{u['leverage_band']:10.4f}{u['leverage_outer']:11.4f}"
              f"{u['leverage_if_fitted_separately']:10.4f}"
              f"{100*u['lever_arm_share']:7.1f}%"
              f"{100*row['halo_template_fraction_in_band']:10.1f}%")

    # ---- can the band measure the halo at all? -----------------------------
    # The forecast is only worth running if a band-restricted fit is a fit and
    # not a shrug, so it is gated on the same injection test as everything else.
    print("\nHalo fit restricted to each region (bin 5, 20.8 GeV):")
    results["restricted_fits"] = {}
    for tag, keep in (("full", full), ("band", band), ("outer", outer)):
        st, Jr = restrict_stack(bg, keep, J)
        fit = fit_all_bins_profile(st, Jr, halo_mask=keep)["per_bin"]
        slopes = injection_slope(st, Jr, [0.0, 100e-35, 300e-35],
                                 halo_mask=keep)["slopes"]
        per_bin = []
        for k in range(9):
            tot = float(Jr[k][keep].sum())
            per_bin.append({
                "bin": k,
                "counts": fit[k]["signal_amplitude"] * tot,
                "lo": fit[k]["signal_amplitude_lo"] * tot,
                "hi": fit[k]["signal_amplitude_hi"] * tot})
        results["restricted_fits"][tag] = {
            "per_bin": per_bin,
            "injection_slope_min": float(min(slopes)),
            "injection_slope_max": float(max(slopes))}
        b5 = per_bin[5]
        half = 50 * (b5["hi"] - b5["lo"]) / abs(b5["counts"])
        print(f"  {tag:6s} {b5['counts']:8.0f} [{b5['lo']:7.0f}, {b5['hi']:7.0f}]"
              f"  +-{half:5.1f}%   injection slope "
              f"{min(slopes):.3f}-{max(slopes):.3f}")

    out = ROOT / "results" / "fits" / "leverage_2026-09-06.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out}")
    print("peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
