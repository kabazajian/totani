#!/usr/bin/env python
"""WHERE inside the eROSITA edge does the boundary excess sit? Spatial breakdown.

    python scripts/boundary_breakdown.py                  # both bubble templates
    python scripts/boundary_breakdown.py --bubble macias  # baseline only

🚨 THE QUESTION. `boundary_test_2026-09-08` gives one number per energy: the
inside-minus-outside excess at matched theta_GC, +404 +- 200 counts at 20.8 GeV
(14.7 sigma at 1.5 GeV) on the structured baseline. It does not say where
inside the footprint that emission is. Two readings were left open:

  * it HUGS THE FERMI-BUBBLE EDGE -- spillover from a bubble template whose
    edge is too sharp or misplaced, i.e. a bubble-modelling defect;
  * it FILLS THE OUTER REGION EVENLY -- a component bounded by the X-ray edge
    (or a low-energy diffuse foreground wrong across the whole area).

And two open items look like the same thing: the quadrant failure (item 7) and
the antisymmetric |b| = 45-60 residual (item 10). So the inside region is cut
three ways, by distance from the Fermi-bubble edge, by sign of l, and by
latitude band, and each piece is matched in theta_GC against the SAME full
outside region (`erosita.matched_theta_weights_regions`). One common reference
makes the pieces comparable with each other; the statistic is reported per
matched square degree for the same reason.

🚨 THE 'all' SUB-REGION IS THE WHOLE INSIDE, AND MUST REPRODUCE THE 09-08 RESULT.
It goes through the generalised weights and is compared against
`boundary_test_2026-09-08.json` bin by bin; the JSON records the worst relative
mismatch and the script refuses to write a result if it is not at rounding.

⚠ THE SAME CAVEATS AS THE PARENT TEST. No support inside theta_GC ~ 34.5 deg;
the Fermi bubbles are excluded from both sides with the fixed Ackermann
footprint in both runs; a sub-region whose theta range has no outside partner
loses those annuli (recorded as dropped area). Sub-regions at the edge bands
are narrow in theta, so their matched outside is a narrow ring too -- the
error bars will be larger, and that is the honest size of the question.

MEMORY: ~2 GB of weight stacks per hemisphere (10 sub-regions x 3 stacks x
66 MB), doubled inside `region_sums`, on top of the template cubes, so
hemispheres are run as SEPARATE `region_sums` calls. MEASURED PEAK 12.2 GB
(2026-09-18, both templates) -- the largest job in this repo. The machine is
shared -- check the other queue's RSS first.
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
from totani.cells import halo_cell_mask  # noqa: E402
from totani.erosita import (  # noqa: E402
    check_handedness,
    erosita_fraction_on_target,
    matched_theta_weights_regions,
    theta_gc,
)
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402
from totani.patches import region_sums  # noqa: E402
from totani.stack import _strip_to_2d, build_background_stack, halo_cell_template  # noqa: E402
from totani.templates import bubble_on_target, grid_lb, pixel_solid_angle  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9
DEG2 = (180.0 / np.pi) ** 2
#: The parent boundary test this must reproduce: the 09-08 run for the 09-07
#: baseline, or the same-dated `_noLoopI` run for the 2026-09-18 baseline.
PARENT = {(): ROOT / "results" / "fits" / "boundary_test_2026-09-08.json",
          ("loopI_shell1", "loopI_shell2"):
              ROOT / "results" / "fits" / "boundary_test_2026-09-18_noLoopI.json"}

#: Distance-from-the-Fermi-bubble-edge bands, degrees.
EDGE_BANDS = ((0.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, 90.0))
#: Latitude bands, degrees of |b|. 45-60 is open item 10's band.
LAT_BANDS = ((10.0, 30.0), (30.0, 45.0), (45.0, 60.0))


def fb_edge_distance_deg(L, B, fb) -> np.ndarray:
    """Great-circle distance from every pixel to the nearest pixel on the
    Fermi-bubble footprint's edge. Exact on the sphere (unit vectors, KD-tree,
    chord -> angle), not a CAR pixel distance, which would be stretched by
    1/cos b at the latitudes the bubbles reach."""
    from scipy.ndimage import binary_erosion
    from scipy.spatial import cKDTree

    edge = fb & ~binary_erosion(fb)
    if edge.sum() == 0:
        raise ValueError("Fermi-bubble footprint has no edge pixels")

    def unit(l, b):
        l, b = np.radians(l), np.radians(b)
        return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], -1)

    tree = cKDTree(unit(L[edge], B[edge]))
    chord, _ = tree.query(unit(L.ravel(), B.ravel()), k=1)
    return np.degrees(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))).reshape(L.shape)


def sub_regions(L, B, dist_fb) -> dict:
    """Boolean pixel maps, applied on top of `inside * keep`."""
    out = {"all": np.ones_like(L, dtype=bool)}
    for lo, hi in EDGE_BANDS:
        tag = f"edge_{lo:g}-{hi:g}" if hi < 90 else f"edge_{lo:g}+"
        out[tag] = (dist_fb >= lo) & (dist_fb < hi)
    out["l_pos"] = L > 0
    out["l_neg"] = L < 0
    for lo, hi in LAT_BANDS:
        out[f"b_{lo:g}-{hi:g}"] = (np.abs(B) >= lo) & (np.abs(B) < hi)
    return out


def build_regions(run_dir, *, theta_step_deg=0.25, node_weight="min"):
    """Per hemisphere: weight stacks for every sub-region vs the full outside."""
    from totani.exposure import flat_spectrum_bin_exposure
    from totani.pointsources import default_catalog_path, extended_source_mask

    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        flat_sr = flat_spectrum_bin_exposure(
            h[0].data.astype(float), h["ENERGIES"].data["Energy"].astype(float), EDGES)

    L, B = grid_lb(target)
    omega = pixel_solid_angle(target)
    flat_sr = flat_sr * omega
    inside = erosita_fraction_on_target(target)
    handed = check_handedness(L, B, inside, omega)
    with fits.open(default_catalog_path()) as h:
        ext_mask = extended_source_mask(L, B, h["ExtendedSources"].data)
    fb = bubble_on_target(target) > 0
    dist_fb = fb_edge_distance_deg(L, B, fb)
    subs = sub_regions(L, B, dist_fb)

    per_hemi = {}
    geometry = {"handedness_mean_lon_deg": handed, "exclude_bubble": "ackermann",
                "edge_bands_deg": EDGE_BANDS, "lat_bands_deg": LAT_BANDS,
                "hemispheres": {}}
    for lobe, sign in (("north", +1), ("south", -1)):
        keep = (sign * B >= 10.0) & ~fb & ~ext_mask
        f_out = (1.0 - inside) * keep
        weights, geo = {}, {}
        for tag, S in subs.items():
            f_in = inside * keep * S
            if f_in.sum() == 0:
                continue
            w_in = np.empty((N_BINS,) + L.shape)
            w_out = np.empty((N_BINS,) + L.shape)
            infos = []
            try:
                for k in range(N_BINS):
                    a, b, info = matched_theta_weights_regions(
                        L, B, omega, f_in, f_out, flat_sr[k],
                        theta_step_deg=theta_step_deg, node_weight=node_weight)
                    w_in[k], w_out[k] = a, b
                    infos.append(info)
            except ValueError as e:      # no annulus with both sides
                geo[tag] = {"skipped": str(e)}
                continue
            weights[f"{tag}_in"] = w_in
            weights[f"{tag}_out"] = w_out
            weights[f"{tag}_var"] = (w_in - w_out) ** 2
            g = {k: v for k, v in infos[0].items()}
            for key in ("inside_area", "outside_area",
                        "inside_area_dropped", "outside_area_dropped"):
                g[f"{key}_deg2"] = g[f"{key}_sr"] * DEG2
            g["unmatched_inside_area_deg2"] = float((omega * f_in).sum() * DEG2)
            g["mean_theta_gc_deg"] = float((omega * f_in * theta_gc(L, B)).sum()
                                           / (omega * f_in).sum())
            geo[tag] = g
        per_hemi[lobe] = weights
        geometry["hemispheres"][lobe] = geo
    return per_hemi, geometry


def run_one(bubble, per_hemi, geometry, drop=()) -> dict:
    t0 = time.time()
    bg = build_background_stack(RUN, EDGES, bubble=bubble, drop_components=drop)
    J = halo_cell_template(RUN, EDGES)
    names = sorted(bg["templates"])
    halo_keep = halo_cell_mask()
    sr, Jr = restrict_stack(bg, halo_keep, J)
    bg_only = fit_all_bins(sr, halo_mask=halo_keep)
    with_halo = fit_all_bins_profile(sr, Jr, halo_mask=halo_keep, with_uncertainty=False)
    roi_J = [float(Jr[k][halo_keep].sum()) for k in range(N_BINS)]
    halo_counts = [with_halo["per_bin"][k]["signal_amplitude"] * roi_J[k]
                   for k in range(N_BINS)]
    del bg, sr
    print(f"  [{bubble}] stack + fits ({time.time()-t0:.0f}s)", flush=True)

    out = {"bubble": bubble, "fitted_halo_counts_roi": halo_counts, "hemispheres": {}}
    for lobe, weights in per_hemi.items():
        rs = region_sums(RUN, EDGES, weights, bubble=bubble, drop_components=drop)
        print(f"  [{bubble}] {lobe} region sums, {len(weights)//3} sub-regions "
              f"({time.time()-t0:.0f}s)", flush=True)

        def mu(region, k, amps):
            return sum(amps[n] * rs["templates"][region][n][k] for n in names)

        tags = sorted({n.rsplit("_", 1)[0] for n in weights})
        rows_by_tag = {}
        for tag in tags:
            g = geometry["hemispheres"][lobe][tag]
            area_in = g["inside_area_deg2"]
            rows = []
            for k in range(N_BINS):
                amps = bg_only["per_bin"][k]["amplitudes"]
                n_in, n_out = rs["counts"][f"{tag}_in"][k], rs["counts"][f"{tag}_out"][k]
                m_in, m_out = mu(f"{tag}_in", k, amps), mu(f"{tag}_out", k, amps)
                excess = (n_in - m_in) - (n_out - m_out)
                sd = float(np.sqrt(mu(f"{tag}_var", k, amps)))
                h_in = halo_counts[k] * rs["halo"][f"{tag}_in"][k] / roi_J[k]
                h_out = halo_counts[k] * rs["halo"][f"{tag}_out"][k] / roi_J[k]
                lobe_counts = excess - (h_in - h_out)
                rows.append({
                    "bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
                    "n_inside": n_in, "mu_inside": m_in,
                    "n_outside": n_out, "mu_outside": m_out,
                    "excess": excess, "sigma_excess": sd,
                    "halo_predicted_diff": h_in - h_out,
                    "halo_counts_inside": h_in,
                    "lobe": lobe_counts,
                    "lobe_sigma": lobe_counts / sd if sd else None,
                    "lobe_per_deg2": lobe_counts / area_in if area_in else None,
                    "lobe_sigma_per_deg2": sd / area_in if area_in else None,
                    "lobe_over_model_inside": lobe_counts / m_in if m_in else None,
                })
            rows_by_tag[tag] = rows
        out["hemispheres"][lobe] = rows_by_tag
        del rs

    # Combined N+S per sub-region: disjoint weight maps, exact sums.
    combined = {}
    for tag in out["hemispheres"]["north"]:
        if tag not in out["hemispheres"]["south"]:
            continue
        rows = []
        for k in range(N_BINS):
            n = out["hemispheres"]["north"][tag][k]
            s = out["hemispheres"]["south"][tag][k]
            lb = n["lobe"] + s["lobe"]
            sd = float(np.hypot(n["sigma_excess"], s["sigma_excess"]))
            area = (geometry["hemispheres"]["north"][tag]["inside_area_deg2"]
                    + geometry["hemispheres"]["south"][tag]["inside_area_deg2"])
            m_in = n["mu_inside"] + s["mu_inside"]
            rows.append({"bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
                         "lobe": lb, "sigma_excess": sd,
                         "lobe_sigma": lb / sd if sd else None,
                         "halo_predicted_diff": n["halo_predicted_diff"] + s["halo_predicted_diff"],
                         "inside_area_deg2": area,
                         "lobe_per_deg2": lb / area if area else None,
                         "lobe_sigma_per_deg2": sd / area if area else None,
                         "lobe_over_model_inside": lb / m_in if m_in else None})
        combined[tag] = rows
    out["hemispheres"]["combined"] = combined
    return out


def check_against_parent(results, parent_path) -> dict:
    """The 'all' sub-region must reproduce the parent boundary test."""
    parent = json.loads(parent_path.read_text())
    worst = {}
    for bubble, run in results["runs"].items():
        p = parent["runs"][bubble]["hemispheres"]["combined"]
        q = run["hemispheres"]["combined"]["all"]
        w = 0.0
        for k in range(N_BINS):
            for key in ("lobe", "sigma_excess"):
                a, b = p[k][key], q[k][key]
                w = max(w, abs(a - b) / max(abs(a), 1e-9))
        worst[bubble] = w
    return worst


def main(argv=None) -> int:
    import resource

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bubble", choices=("ackermann", "macias", "both"), default="both")
    ap.add_argument("--theta-step", type=float, default=0.25)
    ap.add_argument("--theta-weight", choices=("min", "harmonic"), default="min")
    ap.add_argument("--no-loopI", action="store_true",
                    help="2026-09-18 baseline: fit without the two Loop I shells")
    ap.add_argument("--date", default="2026-09-18")
    args = ap.parse_args(argv)
    drop = ("loopI_shell1", "loopI_shell2") if args.no_loopI else ()
    parent_path = PARENT[drop]

    t0 = time.time()
    per_hemi, geometry = build_regions(RUN, theta_step_deg=args.theta_step,
                                       node_weight=args.theta_weight)
    print(f"weights built ({time.time()-t0:.0f}s)", flush=True)
    for lobe, geo in geometry["hemispheres"].items():
        for tag, g in geo.items():
            if "skipped" in g:
                print(f"  {lobe} {tag}: skipped ({g['skipped']})")
            else:
                print(f"  {lobe} {tag:10s}: theta {g['theta_min_deg']:.1f}-{g['theta_max_deg']:.1f}, "
                      f"in {g['inside_area_deg2']:6.0f} deg2 matched "
                      f"(+{g['inside_area_dropped_deg2']:4.0f} dropped), "
                      f"<theta> {g['mean_theta_gc_deg']:.1f}")

    results = {"date": args.date, "run": "fff2de5bc7e0",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]),
               "parent": parent_path.name, "dropped_components": list(drop),
               "theta_step_deg": args.theta_step,
               "theta_weight": args.theta_weight, "geometry": geometry, "runs": {}}
    variants = ("ackermann", "macias") if args.bubble == "both" else (args.bubble,)
    for bubble in variants:
        results["runs"][bubble] = run_one(bubble, per_hemi, geometry, drop)

    if args.theta_step == 0.25 and args.theta_weight == "min":
        worst = check_against_parent(results, parent_path)
        results["reproduces_parent_worst_rel"] = worst
        print(f"\n  'all' vs {parent_path.name}: worst rel mismatch {worst}")
        if max(worst.values()) > 1e-6:
            raise SystemExit("the whole-inside sub-region does not reproduce the "
                             "09-08 result; refusing to write")

    tag = ("_noLoopI" if args.no_loopI else "") + ("" if args.theta_weight == "min" else f"_{args.theta_weight}")
    out = ROOT / "results" / "fits" / f"boundary_breakdown_{args.date}{tag}.json"
    out.write_text(json.dumps(results, indent=2))

    for k in (0, 5):
        print(f"\n  bin {k} ({TOTANI_CENTERS_GEV[k]:.1f} GeV), combined N+S, "
              f"lobe = inside minus theta-matched outside:")
        print(f"  {'sub-region':12s} " + "  ".join(f"{b:>28s}" for b in variants))
        for tag in results["runs"][variants[0]]["hemispheres"]["combined"]:
            cells = []
            for bubble in variants:
                r = results["runs"][bubble]["hemispheres"]["combined"][tag][k]
                cells.append(f"{r['lobe']:+8.0f} ± {r['sigma_excess']:5.0f} ({r['lobe_sigma']:+5.1f}σ) "
                             f"{r['lobe_per_deg2']:+6.2f}/deg²")
            print(f"  {tag:12s} " + "  ".join(cells))
    print(f"\nsaved {out}")
    print(f"total {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
