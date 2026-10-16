#!/usr/bin/env python
"""Inside vs outside the eROSITA Bubble edge, matched in theta_GC.

    python scripts/boundary_test.py                  # both bubble templates
    python scripts/boundary_test.py --bubble macias  # just the baseline

🚨 WHAT THIS TEST IS. Oscar's boundary (`totani.erosita`) is a geometric
footprint reconstructed from Zhang et al. Extended Data Fig. 5(a). It is not
an X-ray brightness map and it never enters the fit as a template. It splits
the halo search region in two, and `matched_theta_weights` reweights the two
halves annulus by annulus so they carry the same theta_GC distribution. A
spherical halo then predicts the same response-corrected intensity on both
sides whatever its radial profile; a lobe bounded by the eROSITA edge does not.
The residual difference, minus the fitted halo's own predicted difference, is
a limit on such a lobe.

🚨 WHY IT IS RUN ON BOTH BUBBLE TEMPLATES, AND WHY THAT IS THE POINT. The
baseline changed on 2026-09-07 to Oscar's structured map on physical grounds,
and the flat Ackermann mask became a systematic variant -- a choice made
against the fit statistic, not against the sky. The two templates predict
DIFFERENT things in exactly this region: outside the Fermi bubbles, inside the
eROSITA edge. So this is the independent on-the-sky check of that choice, and
the two runs are meant to be read side by side.

⚠ THE FERMI BUBBLES ARE EXCLUDED FROM BOTH SIDES, USING THE SAME FIXED
FOOTPRINT IN BOTH RUNS. Zhang et al.'s "outer region" means outside the Fermi
Bubbles and inside the eROSITA Bubbles, and that is the region of interest. The
exclusion uses the Ackermann flat 0/1 footprint in BOTH runs deliberately: if
each run excluded its own template's support the two would be comparing
different patches of sky and the difference between them would be
uninterpretable. What varies between runs is which bubble template is FITTED,
nothing else.

⚠ THIS TEST HAS NO SUPPORT INSIDE theta_GC ~ 34.5 deg -- every pixel of the
halo region closer than that is inside the footprint, so there is nothing to
compare it with at the same angular distance. That is where most of the halo's
leverage lives (54-61% of it is in the 10-30 vs 30-60 deg contrast,
`band_forecast_2026-09-06.json`). This test bounds a lobe's morphology; it is
a weak constraint on the halo amplitude, and the JSON records the dropped area
so that cannot be forgotten.

MEMORY: ~7 GB peak per bubble variant, sequentially, via
`totani.patches.region_sums`. The machine is shared -- check the other queue's
RSS first.
"""
from __future__ import annotations

import argparse
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
from totani import GCE_FISHER_ROOT  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.erosita import (  # noqa: E402
    check_handedness,
    erosita_fraction_on_target,
    matched_theta_weights,
    spherical_contrast,
    theta_gc,
    unmatched_contrast,
)
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402
from totani.patches import region_sums  # noqa: E402
from totani.stack import _strip_to_2d, build_background_stack, halo_cell_template  # noqa: E402
from totani.templates import (  # noqa: E402
    bubble_on_target,
    grid_lb,
    pixel_solid_angle,
)

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9

#: Reused from `patch_test.py`: a box that IS an exact union of cells (row
#: b = 35, columns l = 25 and 35), carried through the pixel path to prove it
#: still reproduces the cell-space stack the fit actually uses.
VALIDATION_BOX = (20.0, 40.0, 30.0, 40.0)


def build_regions(run_dir, *, theta_step_deg=0.25, node_weight="min",
                  exclude_bubble="ackermann"):
    """Weight maps for the matched inside/outside comparison, per hemisphere.

    Returns `(weights, geometry)`. `weights` holds, for each hemisphere, an
    `_in` and `_out` stack and a `_var` stack equal to (w_in - w_out)^2, each
    (n_bins, ny, nx).

    🚨 THE WEIGHTS ARE PER ENERGY BIN BECAUSE THE MATCHING IS RESPONSE-WEIGHTED.
    A spherical component contributes `f(theta) x exposure x omega` counts, so
    cancelling it means normalising each annulus by exposure x solid angle --
    `flat_sr`, the same response every template in the stack is built through
    -- and the exposure changes with energy. Matching on solid angle alone left
    the fitted halo predicting 34% more counts on one side than the other; see
    `erosita.matched_theta_weights`.

    🚨 THE `_var` STACK IS NOT DECORATION. Boundary pixels are PARTIALLY inside,
    so they carry weight on BOTH sides and their Poisson fluctuation enters the
    difference with a cancelling sign. Var(D) = sum_i (a_i - b_i)^2 mu_i, which
    is NOT sum_i a_i^2 mu_i + sum_i b_i^2 mu_i. Summing the templates against
    (w_in - w_out)^2 gives the correct variance from the fitted normalisations
    directly; treating the sides as independent would overstate the error bar
    on every partially covered pixel.

    🚨 THE EXTENDED-SOURCE MASK IS FOLDED IN BEFORE THE MATCHING, NOT AFTER.
    `region_sums` also multiplies by it, which is then a no-op. If it were
    applied only afterwards the annulus normalisations would be computed over
    pixels that contribute nothing, and the theta_GC matching -- the entire
    null -- would be wrong by whatever the mask removes.
    """
    from totani.exposure import flat_spectrum_bin_exposure
    from totani.pointsources import default_catalog_path, extended_source_mask

    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        flat_sr = flat_spectrum_bin_exposure(
            h[0].data.astype(float),
            h["ENERGIES"].data["Energy"].astype(float), EDGES)

    L, B = grid_lb(target)
    omega = pixel_solid_angle(target)
    flat_sr = flat_sr * omega
    inside = erosita_fraction_on_target(target)
    handed = check_handedness(L, B, inside, omega)

    with fits.open(default_catalog_path()) as h:
        ext_mask = extended_source_mask(L, B, h["ExtendedSources"].data)

    if exclude_bubble == "ackermann":
        fb = bubble_on_target(target, GCE_FISHER_ROOT) > 0
    elif exclude_bubble == "none":
        fb = np.zeros_like(L, dtype=bool)
    else:
        raise ValueError(f"exclude_bubble must be 'ackermann' or 'none', "
                         f"got {exclude_bubble!r}")

    th = theta_gc(L, B)
    weights, geometry = {}, {"handedness_mean_lon_deg": handed,
                             "exclude_bubble": exclude_bubble,
                             "hemispheres": {}}

    for lobe, sign in (("north", +1), ("south", -1)):
        keep = (sign * B >= 10.0) & ~fb & ~ext_mask
        w_in = np.empty((N_BINS,) + L.shape)
        w_out = np.empty((N_BINS,) + L.shape)
        per_bin = []
        for k in range(N_BINS):
            a, b, info = matched_theta_weights(
                L, B, omega, inside, keep, flat_sr[k],
                theta_step_deg=theta_step_deg, node_weight=node_weight)
            w_in[k], w_out[k] = a, b
            info["spherical_contrast_matched"] = spherical_contrast(
                a, b, th, flat_sr[k])
            info["spherical_contrast_unmatched"] = unmatched_contrast(
                inside, keep, omega, th, flat_sr[k])
            per_bin.append(info)

        weights[f"{lobe}_in"] = w_in
        weights[f"{lobe}_out"] = w_out
        weights[f"{lobe}_var"] = (w_in - w_out) ** 2
        # The same two regions with no matching at all -- the naive number,
        # kept so the difference the matching makes is on the record.
        weights[f"{lobe}_raw_in"] = inside * keep
        weights[f"{lobe}_raw_out"] = (1.0 - inside) * keep

        # The geometry is energy independent; only the response is not. Report
        # it once, and the cancellation for every bin.
        g = {k: v for k, v in per_bin[0].items()
             if not k.startswith("spherical_contrast")}
        deg2 = (180.0 / np.pi) ** 2
        for key in ("inside_area", "outside_area",
                    "inside_area_dropped", "outside_area_dropped"):
            g[f"{key}_deg2"] = g[f"{key}_sr"] * deg2
        g["spherical_contrast_matched_per_bin"] = [
            i["spherical_contrast_matched"] for i in per_bin]
        g["spherical_contrast_unmatched_per_bin"] = [
            i["spherical_contrast_unmatched"] for i in per_bin]
        g["worst_matched_contrast"] = max(
            abs(v) for i in per_bin for v in i["spherical_contrast_matched"].values())
        g["worst_unmatched_contrast"] = max(
            abs(v) for i in per_bin for v in i["spherical_contrast_unmatched"].values())
        geometry["hemispheres"][lobe] = g

    l0, l1, b0, b1 = VALIDATION_BOX
    weights["_validation"] = ((L >= l0) & (L < l1) & (B >= b0) & (B < b1)).astype(float)
    return weights, geometry


def _validate_pixel_path(rs, bg, names):
    """The pixel path must still reproduce the cell-space stack exactly."""
    worst = 0.0
    for n in names:
        a = np.array(rs["templates"]["_validation"][n])
        b = np.array([float(bg["templates"][n][k, 9, 8:10].sum()) for k in range(N_BINS)])
        good = b != 0
        if good.any():
            worst = max(worst, float(np.abs(a[good] / b[good] - 1).max()))
    return worst


def run_one(bubble, weights, geometry, drop=()) -> dict:
    """The whole test for one bubble template. `drop` names components left out
    of the fit AND the region sums (2026-09-18 baseline: the Loop I shells)."""
    t0 = time.time()
    rs = region_sums(RUN, EDGES, weights, bubble=bubble, drop_components=drop)
    print(f"  [{bubble}] region sums built ({time.time()-t0:.0f}s)", flush=True)

    bg = build_background_stack(RUN, EDGES, bubble=bubble, drop_components=drop)
    J = halo_cell_template(RUN, EDGES)
    names = sorted(bg["templates"])
    worst = _validate_pixel_path(rs, bg, names)
    print(f"  [{bubble}] pixel path vs cell stack: {worst:.2e}", flush=True)

    halo_keep = halo_cell_mask()
    sr, Jr = restrict_stack(bg, halo_keep, J)
    bg_only = fit_all_bins(sr, halo_mask=halo_keep)
    with_halo = fit_all_bins_profile(sr, Jr, halo_mask=halo_keep,
                                     with_uncertainty=False)
    print(f"  [{bubble}] fits done ({time.time()-t0:.0f}s)", flush=True)

    def mu(region, k, amps):
        return sum(amps[n] * rs["templates"][region][n][k] for n in names)

    roi_J = [float(Jr[k][halo_keep].sum()) for k in range(N_BINS)]
    halo_counts = [with_halo["per_bin"][k]["signal_amplitude"] * roi_J[k]
                   for k in range(N_BINS)]

    out = {"bubble": bubble,
           "validation_worst_rel_mismatch": worst,
           "fitted_halo_counts_roi": halo_counts,
           "hemispheres": {}}

    for lobe in ("north", "south"):
        # How much a SPHERICAL component can leak into the statistic, in
        # counts: the annulus-width residual of the matching (1-2e-5, see
        # `tests/test_erosita.py`) times the halo's own counts inside. Reported
        # per bin so it can be compared against the error bar rather than
        # taken on trust as "negligible".
        leak = geometry["hemispheres"][lobe]["worst_matched_contrast"]
        rows = []
        for k in range(N_BINS):
            amps = bg_only["per_bin"][k]["amplitudes"]
            n_in = rs["counts"][f"{lobe}_in"][k]
            n_out = rs["counts"][f"{lobe}_out"][k]
            m_in = mu(f"{lobe}_in", k, amps)
            m_out = mu(f"{lobe}_out", k, amps)
            excess = (n_in - m_in) - (n_out - m_out)
            # Var(D) = sum_i (a_i - b_i)^2 mu_i, evaluated through the _var map.
            var = mu(f"{lobe}_var", k, amps)
            sd = float(np.sqrt(var))

            # The fitted halo's own predicted difference. Small by construction
            # -- the matching cancels any spherical INTENSITY -- but not zero,
            # because exposure and the extended-source mask are not spherical.
            h_in = halo_counts[k] * rs["halo"][f"{lobe}_in"][k] / roi_J[k]
            h_out = halo_counts[k] * rs["halo"][f"{lobe}_out"][k] / roi_J[k]
            lobe_counts = excess - (h_in - h_out)

            raw_in = rs["counts"][f"{lobe}_raw_in"][k] - mu(f"{lobe}_raw_in", k, amps)
            raw_out = rs["counts"][f"{lobe}_raw_out"][k] - mu(f"{lobe}_raw_out", k, amps)

            rows.append({
                "bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
                "n_inside": n_in, "mu_inside": m_in,
                "n_outside": n_out, "mu_outside": m_out,
                "excess": excess, "sigma_excess": sd,
                "halo_predicted_diff": h_in - h_out,
                "halo_counts_inside": h_in,
                "lobe": lobe_counts, "lobe_sigma": lobe_counts / sd if sd else None,
                "lobe_over_halo_95pct_UL": ((lobe_counts + 2 * sd) / h_in
                                            if h_in > 0 else None),
                "unmatched_residual_diff": raw_in - raw_out,
                "spherical_leak_bound_counts": leak * h_in,
            })
        out["hemispheres"][lobe] = rows

    # North and south weight maps are disjoint (|b| >= 10 either side), so the
    # combined statistic and its variance are exact sums.
    combined = []
    for k in range(N_BINS):
        n = out["hemispheres"]["north"][k]
        s = out["hemispheres"]["south"][k]
        lb = n["lobe"] + s["lobe"]
        sd = float(np.hypot(n["sigma_excess"], s["sigma_excess"]))
        hA = n["halo_counts_inside"] + s["halo_counts_inside"]
        combined.append({
            "bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
            "excess": n["excess"] + s["excess"], "sigma_excess": sd,
            "halo_predicted_diff": n["halo_predicted_diff"] + s["halo_predicted_diff"],
            "halo_counts_inside": hA,
            "lobe": lb, "lobe_sigma": lb / sd if sd else None,
            "lobe_over_halo_95pct_UL": (lb + 2 * sd) / hA if hA > 0 else None,
        })
    out["hemispheres"]["combined"] = combined
    return out


def main(argv=None) -> int:
    import resource

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bubble", choices=("ackermann", "macias", "both"),
                    default="both", help="bubble template to fit (default: both)")
    ap.add_argument("--theta-step", type=float, default=0.25,
                    help="theta_GC annulus width in degrees (default: 0.25)")
    ap.add_argument("--theta-weight", choices=("min", "harmonic"), default="min",
                    help="annulus weight; 'min' matches Oscar's own quadrature")
    ap.add_argument("--exclude-bubble", choices=("ackermann", "none"),
                    default="ackermann",
                    help="Fermi-bubble footprint excluded from BOTH sides")
    ap.add_argument("--no-loopI", action="store_true",
                    help="2026-09-18 baseline: fit without the two Loop I shells")
    ap.add_argument("--date", default="2026-09-08")
    args = ap.parse_args(argv)
    drop = ("loopI_shell1", "loopI_shell2") if args.no_loopI else ()

    t0 = time.time()
    weights, geometry = build_regions(
        RUN, theta_step_deg=args.theta_step, node_weight=args.theta_weight,
        exclude_bubble=args.exclude_bubble)
    for lobe, info in geometry["hemispheres"].items():
        print(f"{lobe}: theta {info['theta_min_deg']:.2f}-{info['theta_max_deg']:.2f} deg, "
              f"in {info['inside_area_deg2']:.0f} / out {info['outside_area_deg2']:.0f} deg2, "
              f"dropped {info['inside_area_dropped_deg2']:.0f} deg2 inside "
              f"(no outside at that theta)")
        print(f"  spherical profiles cancel to {info['worst_matched_contrast']:.1e} "
              f"matched, {info['worst_unmatched_contrast']:.1%} unmatched")

    results = {"date": args.date, "run": "fff2de5bc7e0",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]),
               "boundary": "Macias reconstruction of Zhang+ 2408.06312v2 ED Fig 5(a)",
               "theta_step_deg": args.theta_step,
               "theta_weight": args.theta_weight,
               "exclude_bubble": args.exclude_bubble,
               "dropped_components": list(drop),
               "geometry": geometry,
               "runs": {}}

    variants = ("ackermann", "macias") if args.bubble == "both" else (args.bubble,)
    for bubble in variants:
        results["runs"][bubble] = run_one(bubble, weights, geometry, drop)

    # Non-default options get their own file: a systematic variant must never
    # overwrite the baseline it is being compared against.
    tag = "_noLoopI" if args.no_loopI else ""
    if args.theta_weight != "min":
        tag += f"_{args.theta_weight}"
    if args.exclude_bubble != "ackermann":
        tag += f"_fb-{args.exclude_bubble}"
    if args.theta_step != 0.25:
        tag += f"_step{args.theta_step:g}"
    out = ROOT / "results" / "fits" / f"boundary_test_{args.date}{tag}.json"
    out.write_text(json.dumps(results, indent=2))

    print("\n  bin 5 (20.8 GeV), lobe = inside minus theta-matched outside:")
    for bubble in variants:
        for lobe in ("north", "south", "combined"):
            r = results["runs"][bubble]["hemispheres"][lobe][5]
            ul = r["lobe_over_halo_95pct_UL"]
            ul_s = f"{ul:+.2f}" if ul is not None else "n/a"
            print(f"    {bubble:9s} {lobe:8s}: {r['lobe']:+8.0f} +- {r['sigma_excess']:6.0f} "
                  f"({r['lobe_sigma']:+5.2f}s)   95% UL / halo = {ul_s}")
    print(f"\nsaved {out}")
    print(f"total {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
