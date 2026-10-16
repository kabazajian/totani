#!/usr/bin/env python
"""Oscar's item 5: the residual-template construction, written out and checked
on a noise-free synthetic sky with known flat bubbles and a correctly
specified background.

    python scripts/reconstruction_check.py        # pixel level, ~8 GB, a few minutes

The sequence implemented in totani.bubble_totani (his Secs. 3.1-3.2):
  1. at 4.3 GeV, fit the 144 cells (disk included) with the diffuse model,
     the FLAT bubble mask, Loop I shells, isotropic, point sources and a
     GC-excess template (NFW rho^2.5);
  2. image = counts - sum over every fitted component EXCEPT the bubble
       = fitted flat-bubble contribution + fit residual      (the add-back);
  3. divide by exposure x solid angle (flux), mask extended sources;
  4. smooth with a 1 deg Gaussian (weighted by the valid-pixel mask);
  5. split into positive and negative parts; each becomes an
     energy-independent template with a free normalization per bin (the
     negative part's sign free per bin).

Check: build a noise-free 4.3 GeV sky = the fitted model of Totani's
configuration (flat bubbles at their fitted amplitude, no halo), run steps
1-5 on it, and compare the positive image with the smoothed flat mask at the
fitted amplitude: agreement inside the footprint, leakage outside, and the
size of the negative part.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani.bubble_totani import bubble_image, fit_all_cells, gaussian_smooth_sky, pixel_model_for_bin  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-24"


def main() -> int:
    pm = pixel_model_for_bin(RUN, EDGES, 2, diffuse="galprop_totani", with_gce=True)
    fit = fit_all_cells(pm)
    amps = fit["amplitudes"]
    model = sum(amps[n] * pm["templates"][n] for n in fit["names"])            # noise-free sky, flat bubbles at the fitted amplitude
    pm2 = dict(pm); pm2["counts"] = model
    fit2 = fit_all_cells(pm2)
    img = bubble_image(pm2, fit2)
    valid = img["valid"]
    # expected: the flat bubble's flux, smoothed the same way
    with np.errstate(divide="ignore", invalid="ignore"):
        flat_flux = np.where(valid, amps["bubble"] * pm["templates"]["bubble"] / pm["flat_sr"], 0.0)
    expected = gaussian_smooth_sky(flat_flux, valid.astype(float), pm["B"], 1.0); expected[~valid] = 0.0
    foot = pm["templates"]["bubble"] > 0
    pos, neg = img["positive"], img["negative"]
    inside = foot & valid; outside = ~foot & valid
    rel_in = np.abs(pos - expected)[inside].max() / expected[inside].max()
    out = {"date": DATE,
           "refit_amplitude_ratios": {n: fit2["amplitudes"][n] / amps[n] for n in fit["names"] if amps[n] > 0},
           "positive_vs_expected_max_rel_dev_inside_footprint": float(rel_in),
           "positive_flux_outside_over_inside": float(pos[outside].sum() / pos[inside].sum()),
           "negative_over_positive_total": float(neg[valid].sum() / pos[valid].sum()),
           "negative_max_over_positive_max": float(neg[valid].max() / pos[inside].max())}
    (ROOT / "results" / "fits" / f"reconstruction_check_{DATE}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
