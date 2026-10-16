#!/usr/bin/env python
"""Does the eROSITA-bounded residual move the halo? The systematic variant.

    python scripts/erosita_variant.py            # ~3 min, one 7 GB stack build

🚨 WHAT THIS ADDS AND WHY. `boundary_test_2026-09-08` found +404 +- 200 counts
at 20.8 GeV (2.0 sigma; 14.7 sigma at 1.5 GeV) inside Oscar's eROSITA edge and
outside the Fermi bubbles that the baseline model does not account for. The
systematics envelope (`systematics_2026-09-07.json`, 3.5x at 20.8 GeV) is built
by SWAPPING background choices; none of its variants contains a template for
that region, so however the halo responds to that unmodelled emission is not
in the envelope. Kev's call, 2026-09-18: measure it.

The variant adds ONE template to the stack: the region inside the eROSITA
boundary and outside the Ackermann flat Fermi-bubble footprint -- Zhang et
al.'s "outer region" -- as a flat-intensity shape with a free normalisation
per bin, through the same pipeline as every other diffuse template
(`stack.cell_template_from_shape`). A north/south split of it is also run,
because the bubbles want separate N/S normalisations at ~3 sigma
(`bubble_NS`) and Kev adopted that split the same day.

⚠ WHAT THIS CAN AND CANNOT SAY. A flat-intensity lobe template of this size is
partly degenerate with a spherical halo over the |b| >= 10 cells (both are
smooth, both peak toward the GC in cell space). So the halo's shift is an
upper bound on the bias the unmodelled emission induces, not a measurement of
the emission itself; the fitted amplitude of the added template is reported
alongside so the two can be told apart. If the shift is small, the residual is
not a halo systematic. If it is large, that is a degeneracy statement as much
as a sky statement, and `boundary_breakdown.py` is what says which.

Every variant is scored the way `systematics.py` scores its variants: halo
counts with 68% interval, TS, the quadrant diagnostic, injection slopes. The
quadrant diagnostic is the interesting one here -- if the eROSITA template is
what the four quadrants have been disagreeing about (open item 7), adding it
should bring them together.

MEMORY: one `build_background_stack` (~7 GB peak); everything else is
cell-space. The machine is shared -- check the other queue's RSS first.
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
from systematics import quadrant_diagnostic, replace_component, split_ns  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.erosita import check_handedness, erosita_fraction_on_target  # noqa: E402
from totani.likelihood import (  # noqa: E402
    fit_all_bins,
    fit_all_bins_profile,
    injection_slope,
    restrict_stack,
)
from totani.stack import (  # noqa: E402
    _strip_to_2d,
    build_background_stack,
    cell_template_from_shape,
    halo_cell_template,
)
from totani.templates import bubble_on_target, grid_lb, pixel_solid_angle  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9


def erosita_outer_shapes(run_dir) -> dict:
    """The 'outer region' as fractional-coverage maps: inside the eROSITA edge,
    outside the flat Fermi-bubble footprint; whole, and split N/S."""
    with fits.open(Path(run_dir) / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    L, B = grid_lb(target)
    inside = erosita_fraction_on_target(target)
    check_handedness(L, B, inside, pixel_solid_angle(target))
    fb = bubble_on_target(target) > 0
    outer = inside * (~fb)
    return {"erosita_outer": outer,
            "erosita_outer_N": outer * (B > 0),
            "erosita_outer_S": outer * (B < 0),
            "_fb_footprint": fb.astype(float)}


def add_components(stack, new: dict):
    t = dict(stack["templates"]); t.update(new)
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def score(stack, J, base_bg_only=None) -> dict:
    """`systematics.py`'s scoring, plus what the ADDED components absorb."""
    d = quadrant_diagnostic(stack, J)
    keep = halo_cell_mask()
    sr, Jr = restrict_stack(stack, keep, J)
    d["injection_slopes"] = injection_slope(
        sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
    d["components"] = sorted(stack["templates"])
    fit = fit_all_bins_profile(sr, Jr, halo_mask=keep)["per_bin"]
    roi_J = [float(Jr[k][keep].sum()) for k in range(N_BINS)]
    d["halo_intervals"] = [
        {"bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
         "counts": fit[k]["signal_amplitude"] * roi_J[k],
         "lo": fit[k]["signal_amplitude_lo"] * roi_J[k],
         "hi": fit[k]["signal_amplitude_hi"] * roi_J[k],
         "neg2lnl_with_halo": fit[k]["neg2lnl"]}
        for k in range(N_BINS)]

    # Background-only fit: what each eROSITA / bubble component absorbs, in
    # counts over the halo cells, and the fit quality with no halo at all.
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    absorbed = {}
    for name in sr["templates"]:
        if name.startswith("erosita") or name.startswith("bubble"):
            absorbed[name] = [bg[k]["amplitudes"][name]
                              * float(sr["templates"][name][k][keep].sum())
                              for k in range(N_BINS)]
    d["component_counts_bg_only"] = absorbed
    d["neg2lnl_bg_only"] = [bg[k]["neg2lnl"] for k in range(N_BINS)]
    if base_bg_only is not None:
        d["delta_neg2lnl_bg_only_vs_baseline"] = [
            bg[k]["neg2lnl"] - base_bg_only[k] for k in range(N_BINS)]
    return d


def main() -> int:
    import resource

    t0 = time.time()
    bg = build_background_stack(RUN, EDGES, bubble="macias")
    J = halo_cell_template(RUN, EDGES)
    print(f"baseline stack built ({time.time()-t0:.0f}s)", flush=True)

    shapes = erosita_outer_shapes(RUN)
    cells = {n: cell_template_from_shape(RUN, EDGES, s)
             for n, s in shapes.items() if not n.startswith("_")}
    flat_bubble = cell_template_from_shape(RUN, EDGES, shapes["_fb_footprint"])
    print(f"eROSITA outer templates built ({time.time()-t0:.0f}s)", flush=True)

    keep = halo_cell_mask()
    geometry = {n: {"halo_cells_sum_bin5": float(c[5][keep].sum())}
                for n, c in cells.items()}

    ero = {"erosita_outer": cells["erosita_outer"]}
    ero_ns = {"erosita_outer_N": cells["erosita_outer_N"],
              "erosita_outer_S": cells["erosita_outer_S"]}
    flat = replace_component(bg, "bubble", flat_bubble)
    variants = {
        "baseline": bg,
        "erosita": add_components(bg, ero),
        "erosita_NS": add_components(bg, ero_ns),
        "bubble_NS": split_ns(bg),
        "bubble_NS_erosita_NS": add_components(split_ns(bg), ero_ns),
        "bubble_flat": flat,
        "bubble_flat_erosita": add_components(flat, ero),
    }

    results = {"date": "2026-09-18", "run": "fff2de5bc7e0", "bubble_baseline": "macias",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]),
               "added_template": "inside eROSITA edge (Macias reconstruction of "
                                 "Zhang+ ED Fig 5a) AND outside the Ackermann flat "
                                 "Fermi-bubble footprint; flat intensity; free per bin",
               "geometry": geometry, "variants": {}}

    base_bg_only = None
    for name, st in variants.items():
        d = score(st, J, base_bg_only)
        if name == "baseline":
            base_bg_only = d["neg2lnl_bg_only"]
        results["variants"][name] = d
        h5 = d["halo_intervals"][5]
        q5 = d["summary"]["bin5"]
        ero_counts = d["component_counts_bg_only"]
        e5 = sum(v[5] for k, v in ero_counts.items() if k.startswith("erosita"))
        print(f"  {name:22s} halo@20.8: {h5['counts']:6.0f} [{h5['lo']:.0f}, {h5['hi']:.0f}] "
              f"TS {d['ts'][5]:5.1f}  quad sum/base {q5['sum_over_baseline']:+.2f} "
              f"neg={q5['any_negative']}  erosita counts {e5:6.0f}  "
              f"slope {min(d['injection_slopes']):.2f}-{max(d['injection_slopes']):.2f} "
              f"({time.time()-t0:.0f}s)", flush=True)

    out = ROOT / "results" / "fits" / "erosita_variant_2026-09-18.json"
    out.write_text(json.dumps(results, indent=2))
    base = results["variants"]["baseline"]["halo_intervals"]
    print("\n  halo counts per bin, variant / baseline:")
    print("  " + " ".join(f"{c:>7.1f}" for c in TOTANI_CENTERS_GEV[:N_BINS]))
    for name, d in results["variants"].items():
        print(f"  {name:22s} " + " ".join(
            f"{d['halo_intervals'][k]['counts'] / base[k]['counts']:7.2f}"
            if base[k]["counts"] else "    n/a" for k in range(N_BINS)))
    print(f"\nsaved {out}")
    print(f"total {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
