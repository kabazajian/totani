#!/usr/bin/env python
"""Two questions Oscar's 2026-09-18 email raises, answered in cell space.

    python scripts/bubble_lowE_check.py     # ~6 min; five stack builds, one at a time

Oscar (2026-09-18): the structured map is Ackermann+ 2017 Fig. 12, derived from
1-10 GeV data, so it should work especially well at LOW energies -- "I wonder
if the poorer fit at 1.5-2.6 GeV might be related to some other templates."
And: would a SkyFACT-style controlled spatial flexibility of the bubble
template help, tested on simulated skies with and without an injected halo?

(a) IS THE LOW-ENERGY DEFICIT A POINT-SOURCE EFFECT? The structured map fits
    worse than the flat mask at 1.5-2.6 GeV (delta -2lnL = +243, +45, bg-only,
    `bubble_macias_2026-09-07.md`) and masking point sources instead of
    modelling them removes 97% of the bin-0 halo. So: repeat the
    structured-vs-flat comparison with point sources MASKED (95% containment
    per bin, and 1 deg fixed) instead of modelled. If the +243 goes away, the
    deficit was the point-source template, not the bubble map. Also with the
    dust residuals and Loop I dropped, since those are the other low-energy
    suspects.

(b) A CONTROLLED FLEXIBILITY LADDER, gated by injection. Give the structured
    template 1, 2 (N/S), 4 (quadrants) and 8 (quadrants x |b| <= 30 / > 30)
    independent normalisations per bin. Report the halo, TS, the quadrant
    diagnostic and the injection slope at each rung. The slope is the gate:
    once the bubble freedom can absorb an injected spherical halo the rung is
    telling us about degeneracy, not about the sky -- which is exactly the
    failure Oscar's simulated-sky test is designed to catch, measured here on
    the real stack in one minute.

MEMORY: five `build_background_stack` calls (~7 GB each), strictly
sequential, each freed before the next. Everything else is cell-space.
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
from systematics import drop, quadrant_diagnostic, replace_component  # noqa: E402
from totani.cells import cell_centers, halo_cell_mask  # noqa: E402
from totani.likelihood import (  # noqa: E402
    fit_all_bins,
    fit_all_bins_profile,
    injection_slope,
    quadrant_masks,
    restrict_stack,
)
from totani.pointsources import psf95_radius_deg  # noqa: E402
from totani.stack import (  # noqa: E402
    _strip_to_2d,
    build_background_stack,
    cell_template_from_shape,
    halo_cell_template,
)
from totani.templates import bubble_on_target  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9


def bg_only_neg2lnl(stack):
    keep = halo_cell_mask()
    sr = restrict_stack(stack, keep)
    return [row["neg2lnl"] for row in fit_all_bins(sr, halo_mask=keep)["per_bin"]]


def split_regions(stack, name, regions: dict):
    """Replace component `name` by one copy per (cell-mask) region."""
    t = dict(stack["templates"]); b = t.pop(name)
    for tag, m in regions.items():
        t[f"{name}_{tag}"] = b * m[None, :, :]
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def ladder_regions():
    Lc, Bc = cell_centers()
    q = quadrant_masks()
    # quadrant_masks already excludes the disk; the disk cells are not fitted
    # anyway (restrict_stack), so the split is exact on the fitted cells.
    rungs = {
        "1": {"all": np.ones_like(Bc, dtype=bool)},
        "2_NS": {"N": Bc > 0, "S": Bc < 0},
        "4_quadrants": q,
        "8_quadrants_x_lat": {f"{k}_{lat}": m & ((np.abs(Bc) <= 30) if lat == "lo" else (np.abs(Bc) > 30))
                              for k, m in q.items() for lat in ("lo", "hi")},
    }
    return rungs


def score(stack, J):
    d = quadrant_diagnostic(stack, J)
    keep = halo_cell_mask()
    sr, Jr = restrict_stack(stack, keep, J)
    d["injection_slopes"] = injection_slope(
        sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
    fit = fit_all_bins_profile(sr, Jr, halo_mask=keep)["per_bin"]
    roi_J = [float(Jr[k][keep].sum()) for k in range(N_BINS)]
    d["halo_intervals"] = [{"bin": k, "counts": fit[k]["signal_amplitude"] * roi_J[k],
                            "lo": fit[k]["signal_amplitude_lo"] * roi_J[k],
                            "hi": fit[k]["signal_amplitude_hi"] * roi_J[k]}
                           for k in range(N_BINS)]
    d["n_bubble_params"] = sum(1 for n in stack["templates"] if n.startswith("bubble"))
    return d


def main() -> int:
    import resource
    t0 = time.time()
    results = {"date": "2026-09-18", "run": "fff2de5bc7e0",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]),
               "a_structured_minus_flat_bg_only": {}, "b_flexibility_ladder": {}}

    # ---- (a) structured minus flat, under different treatments ------------
    from astropy.io import fits
    with fits.open(RUN / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    flat_shape = (bubble_on_target(target) > 0).astype(float)

    def compare(tag, st_structured):
        flat = replace_component(st_structured, "bubble", flat_cells[tag])
        a, b = bg_only_neg2lnl(st_structured), bg_only_neg2lnl(flat)
        results["a_structured_minus_flat_bg_only"][tag] = {
            "structured": a, "flat": b, "delta": [x - y for x, y in zip(a, b)]}
        print(f"  (a) {tag:26s} d(-2lnL) structured-flat: "
              + " ".join(f"{x - y:+6.0f}" for x, y in zip(a, b)) + f"  ({time.time()-t0:.0f}s)",
              flush=True)

    flat_cells = {}
    bg = build_background_stack(RUN, EDGES, bubble="macias")
    J = halo_cell_template(RUN, EDGES)
    flat_cells["ps_modelled"] = cell_template_from_shape(RUN, EDGES, flat_shape)
    for tag in ("no_dust", "no_loopI", "no_dust_no_loopI"):
        flat_cells[tag] = flat_cells["ps_modelled"]
    print(f"baseline stack ({time.time()-t0:.0f}s)", flush=True)
    compare("ps_modelled", bg)
    compare("no_dust", drop(bg, {"dust_pos", "dust_neg"}))
    compare("no_loopI", drop(bg, {"loopI_shell1", "loopI_shell2"}))
    compare("no_dust_no_loopI", drop(bg, {"dust_pos", "dust_neg", "loopI_shell1", "loopI_shell2"}))

    # Point-source-masked stacks: the mask is per bin, so the flat template
    # must be rebuilt through the same per-bin mask -- take it from the flat
    # stack itself rather than from `cell_template_from_shape`.
    for tag, radii in (("ps_masked_95pct", [psf95_radius_deg(e) for e in TOTANI_CENTERS_GEV[:9]]),
                       ("ps_masked_fixed_1deg", 1.0)):
        st_s = build_background_stack(RUN, EDGES, ps_mask_radii_deg=radii, bubble="macias")
        st_f = build_background_stack(RUN, EDGES, ps_mask_radii_deg=radii, bubble="ackermann")
        flat_cells[tag] = st_f["templates"]["bubble"]
        del st_f
        compare(tag, st_s)
        del st_s

    # ---- (b) flexibility ladder on the structured template ----------------
    for rung, regions in ladder_regions().items():
        st = split_regions(bg, "bubble", regions)
        d = score(st, J)
        results["b_flexibility_ladder"][rung] = d
        h5 = d["halo_intervals"][5]; q5 = d["summary"]["bin5"]
        print(f"  (b) rung {rung:18s} halo@20.8 {h5['counts']:6.0f} [{h5['lo']:.0f}, {h5['hi']:.0f}] "
              f"TS {d['ts'][5]:5.1f}  quad sum/base {q5['sum_over_baseline']:+.2f} neg={q5['any_negative']}  "
              f"slope min {min(d['injection_slopes']):.2f} bin5 {d['injection_slopes'][5]:.2f}  "
              f"({time.time()-t0:.0f}s)", flush=True)
        print("        halo/bin " + " ".join(f"{h['counts']:7.0f}" for h in d["halo_intervals"]))
        print("        slope    " + " ".join(f"{s:7.2f}" for s in d["injection_slopes"]))

    out = ROOT / "results" / "fits" / "bubble_lowE_check_2026-09-18.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out}\ntotal {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
