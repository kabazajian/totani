"""Per-region, per-bin sums of every template and of the counts, in pixel space.

Both on-the-sky tests need the same thing: the fitted normalisations evaluated
on a region that does not land on the 10-degree cell grid. `scripts/patch_test.py`
needs it for the Zhang+ boxes, whose latitudes cut cells in half;
`scripts/boundary_test.py` needs it for the eROSITA footprint, which cuts 22 of
the 120 cells. This module is that shared step, and it is shared rather than
copied so the two tests cannot drift in how they build a template.

What is reused from the fit is the fitted normalisations -- one scalar per
component per energy bin, valid at any resolution. `mu_region = sum_l f_l
T_l(region)` is exact given those f_l; nothing is re-fitted here.

🚨 REGIONS ARE WEIGHT MAPS, NOT MASKS. A boolean mask is the special case of a
0/1 map, and the box tests use exactly that. The eROSITA regions need real
weights: partial pixel coverage at the boundary, and the theta_GC matching
factors. Carrying both through one float multiply keeps one code path.

MEMORY: each pixel-space template is reduced to its region sums and freed
before the next is built, because the machine is shared and the naive version
(hold all 19 in pixel space) peaks near 7 GB.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits


def nearest_pixel_weight(weight, L, B):
    """Look a weight map up at arbitrary (l, b) -- used for point sources.

    Point sources are NOT a pixel map in this pipeline: the template deposits
    each source's expected counts at its own position (see
    `pointsources.point_source_cell_template`), so a source's region weight is
    the weight of the pixel it falls in, not an average over anything.
    """
    lv, bv = L[0, :], B[:, 0]

    def _idx(vec, val):
        # `vec` may run either way: grid_lb reads CDELT straight off the WCS
        # and CDELT1 is negative on every fermitools product.
        if vec[0] < vec[-1]:
            i = np.searchsorted(vec, val)
        else:
            i = len(vec) - np.searchsorted(vec[::-1], val)
        i = np.clip(i, 1, len(vec) - 1)
        return np.where(np.abs(val - vec[i - 1]) <= np.abs(val - vec[i]), i - 1, i)

    return lambda sl, sb: weight[_idx(bv, sb), _idx(lv, sl)]


def region_sums(run_dir, edges_gev, weights, source_weight_fn=None,
                bubble: str = "ackermann", drop_components=(),
                diffuse: str = "a2020", ics_split: bool = False,
                gas_rings: bool = False, galprop_model_set=None) -> dict:
    """Per-region, per-bin weighted sums of every template, the counts and the halo.

    `weights` maps a region name to a (ny, nx) float map on the target grid,
    or to a (n_bins, ny, nx) stack when the weights differ per energy bin.
    `source_weight_fn(name, source_l, source_b)` returns the per-source weight
    for one region; the default looks the map up at the source's own pixel.

    `bubble` selects the same two templates `build_background_stack` offers and
    MUST match the stack the normalisations came from -- evaluating a fit made
    with one bubble template on region sums built with the other is the
    silent-wrong-answer failure this argument exists to prevent.
    `drop_components` likewise mirrors the stack's: a component absent from
    the fit must be absent here, or `mu_region` would silently miss it (the
    callers sum over the FIT's component names, so an extra one here is
    ignored -- but keeping the two in step is the only safe habit).

    Same pixel-space construction as `stack.build_background_stack` -- same
    solid angle, same exposure bases, same extended-source mask -- reduced to
    weighted sums instead of cell sums, one template at a time.
    """
    from . import GCE_FISHER_ROOT
    from .exposure import flat_spectrum_bin_exposure, single_energy_bin_exposure
    from .halo import nfw_rho2_on_target
    from .pointsources import (default_catalog_path, extended_source_mask,
                               source_expected_counts)
    from .stack import _strip_to_2d
    from .templates import (bubble_macias_on_target, bubble_on_target,
                            gas_dust_on_target, grid_lb, ics_four_ring,
                            ics_on_target, isotropic_on_target, loopi_on_target,
                            pixel_solid_angle)

    if bubble not in ("ackermann", "macias"):
        raise ValueError(f"bubble must be 'ackermann' or 'macias', got {bubble!r}")

    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        counts_cube = h[0].data.astype(float)
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float)
        exposure_mev = h["ENERGIES"].data["Energy"].astype(float)
        exp_header = _strip_to_2d(h[0].header)

    L, B = grid_lb(target)
    with fits.open(default_catalog_path()) as h:
        cat = h["LAT_Point_Source_Catalog"].data
        extended = h["ExtendedSources"].data
    mask = extended_source_mask(L, B, extended)

    # 🚨 THE EXTENDED-SOURCE MASK IS APPLIED HERE, ONCE, TO EVERY REGION.
    # Its photons are cut from the data, so a region that kept those pixels in
    # its templates would compare a model including them against counts that
    # do not. The box test used to do this with `& ~mask`; a weight map does it
    # by multiplication, and the effect is identical for a 0/1 map.
    # 🚨 A REGION MAY BE PER-BIN. The eROSITA matching normalises annuli by
    # exposure x solid angle, and the exposure is energy dependent, so those
    # weights are a (n_bins, ny, nx) stack. A 2-D map means the same weights in
    # every bin, which is what every box region is -- and it is kept 2-D rather
    # than broadcast into a stack, because 9 x 960 x 960 x 8 bytes per region
    # is 66 MB and this runs beside a 16-28 GB queue.
    n_bins = len(edges_gev) - 1
    w = {}
    for n, m in weights.items():
        m = np.asarray(m, dtype=float)
        if m.shape not in (L.shape, (n_bins,) + L.shape):
            raise ValueError(f"region {n!r} has shape {m.shape}, expected "
                             f"{L.shape} or {(n_bins,) + L.shape}")
        w[n] = m * ~mask

    def wk(rw, k):
        return rw if rw.ndim == 2 else rw[k]
    omega = pixel_solid_angle(target)
    flat_sr = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev) * omega
    single_sr = single_energy_bin_exposure(exposure_cube, exposure_mev, edges_gev) * omega

    out = {n: {} for n in w}

    def reduce_and_free(name, cube):
        """cube is (n_bins, ny, nx); store its region sums, then drop it."""
        for rn, rw in w.items():
            out[rn][name] = [float((cube[k] * wk(rw, k)).sum()) for k in range(n_bins)]

    # `diffuse` / `ics_split` mirror `build_background_stack` exactly; see there.
    if diffuse not in ("a2020", "galprop_totani"):
        raise ValueError(f"diffuse must be 'a2020' or 'galprop_totani', got {diffuse!r}")
    if diffuse == "galprop_totani":
        from .galprop import TOTANI_SET, galprop_on_target
        for name, cube in galprop_on_target(target, edges_gev, ics_split=ics_split,
                                            gas_rings=gas_rings,
                                            model_set=galprop_model_set or TOTANI_SET).items():
            reduce_and_free(f"galprop_{name}", cube * single_sr)
    else:
        for name, shape_map in gas_dust_on_target(target, GCE_FISHER_ROOT).items():
            reduce_and_free(name, shape_map[None, :, :] * flat_sr)
    bub = (bubble_on_target(target, GCE_FISHER_ROOT) if bubble == "ackermann"
           else bubble_macias_on_target(target))
    reduce_and_free("bubble", bub[None, :, :] * flat_sr)
    del bub
    for name, shape_map in loopi_on_target(target).items():
        reduce_and_free(name, shape_map[None, :, :] * flat_sr)
    reduce_and_free("isotropic", isotropic_on_target(target)[None, :, :] * flat_sr)
    if diffuse == "a2020":
        six = ics_on_target(target, edges_gev, GCE_FISHER_ROOT)
        for name, cube in zip(("ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6"),
                              ics_four_ring(six)):
            reduce_and_free(name, cube * single_sr)
        del six

    drop_components = tuple(drop_components)
    unknown = [n for n in drop_components if n not in next(iter(out.values()))]
    if unknown:
        raise ValueError(f"drop_components names {unknown} that are not built here")
    for rn in out:
        for n in drop_components:
            del out[rn][n]

    counts = {n: [float((counts_cube[k] * wk(rw, k)).sum()) for k in range(n_bins)]
              for n, rw in w.items()}

    # 🚨 THE HALO TEMPLATE ON THE REGIONS, built by the SAME route as
    # `stack.halo_cell_template`. It is needed because a spherical halo does
    # NOT predict equal COUNTS in two theta-matched regions -- exposure and the
    # extended-source mask break the symmetry that theta_GC preserves.
    # Subtracting the fitted halo's own predicted difference is what turns the
    # residual into a limit on a component no spherical profile can produce.
    halo_cube = (nfw_rho2_on_target(target)[None, :, :] * flat_sr)
    halo = {n: [float((halo_cube[k] * wk(rw, k)).sum()) for k in range(n_bins)]
            for n, rw in w.items()}
    del halo_cube

    cat_l = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"])
    cat_b = cat["GLAT"]
    in_roi = (np.abs(cat_l) <= 60.0) & (np.abs(cat_b) <= 60.0)
    src_keep = in_roi & ~extended_source_mask(cat_l, cat_b, extended)
    sc = source_expected_counts(cat[src_keep], cat_l[src_keep], cat_b[src_keep],
                                edges_gev, exposure_cube, exposure_mev, exp_header)
    sl, sb = cat_l[src_keep], cat_b[src_keep]
    for rn, rw in w.items():
        if source_weight_fn is None:
            # Per bin: a source's weight follows its pixel's weight, which for
            # the eROSITA regions changes with the exposure.
            col = [nearest_pixel_weight(wk(rw, k), L, B)(sl, sb)
                   for k in range(n_bins)]
        else:
            col = [source_weight_fn(rn, sl, sb)] * n_bins
        out[rn]["point_sources"] = [float((sc[:, k] * col[k]).sum())
                                    for k in range(n_bins)]

    return {"templates": out, "counts": counts, "halo": halo, "n_bins": n_bins}
