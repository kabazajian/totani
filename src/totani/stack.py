"""Assemble every background template into per-bin cell-space counts basis.

    stack = build_background_stack(run_dir, edges_gev)
    stack["counts"]              # (n_bins, 12, 12) observed, masked
    stack["templates"]["bubble"] # (n_bins, 12, 12), one f_l per bin multiplies this

🚨 ONE ARRAY PER NAMED COMPONENT, EACH MEANT TO BE SCALED BY EXACTLY ONE FREE
NORMALIZATION PER BIN -- not per pixel, not per source, not per ring-sub-plane.
This is the assembled form of every decision recorded in CONVENTIONS.md and in
`totani.exposure`'s module docstring: gas x8 + dust x2 (A2020 convention, each
ring independently free), ICS x4 (4-ring collapse, CONFIRMED by Oscar), bubble
x1, Loop I x2 (Totani sec 2.3: "two independent fit parameters"), isotropic x1,
point sources x1 (Totani sec 2.3 + sec 3.4.3, CONFIRMED — see CONVENTIONS.md).
15 free normalizations per bin against 120 halo cells, comfortably determined.

⚠ THE HALO/DM TEMPLATE IS DELIBERATELY NOT HERE. This module builds the
BACKGROUND stack only. Adding a halo template and fitting it against real data
is gated by `validation/RULING_real_data.md` and `scripts/hooks/
rule6_halo_guard.py`; keeping this module halo-free means importing or running
it can never itself be the thing that crosses that line.

🚨 EXTENDED-SOURCE MASKING IS APPLIED TO EVERY TEMPLATE AND TO THE DATA, THE
SAME MASK, BEFORE CELL BINNING. Totani sec 2.3 masks circular regions (2x the
catalogued semi-major axis) and excludes them from the analysis. Masking the
data but not a template (or vice versa) would silently bias every cell touching
a masked region -- `bin_to_cells`'s `mask` argument excludes pixels from the
sum entirely, which is what "excluded from the analysis" means, not zeroed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

from .cells import bin_to_cells
from .exposure import flat_spectrum_bin_exposure, single_energy_bin_exposure
from .pointsources import (
    default_catalog_path,
    extended_source_mask,
    point_source_cell_template,
    point_source_mask,
    psf95_radius_deg,
)
from .templates import (
    bubble_macias_on_target,
    bubble_on_target,
    gas_dust_on_target,
    grid_lb,
    ics_four_ring,
    ics_on_target,
    isotropic_on_target,
    loopi_on_target,
    pixel_solid_angle,
)

__all__ = ["build_background_stack", "halo_cell_template"]


def _strip_to_2d(header):
    h = header.copy()
    for k in ("NAXIS3", "CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3"):
        h.remove(k, ignore_missing=True)
    h["NAXIS"] = 2
    return h


def build_background_stack(run_dir, edges_gev, root=None, catalog_path=None,
                           config: str = "baseline", ts: str = "150K",
                           half_deg: float = 60.0, cell_deg: float = 10.0,
                           ps_mask_radii_deg=None, bubble: str = "ackermann",
                           drop_components=(), diffuse: str = "a2020",
                           ics_split: bool = False, gas_rings: bool = False,
                           galprop_model_set=None) -> dict:
    """Every background template and the observed counts, on the cell grid.

    `diffuse` selects the cosmic-ray diffuse model (2026-09-18):
      * "a2020"          -- the Pohl et al. ring templates as delivered for
                            A2020: 4 HI + 4 H2 annuli, 2 dust residuals, 4 IC
                            ring groups (14 templates). Every result before
                            2026-09-18 used this.
      * "galprop_totani" -- TOTANI'S OWN MODEL, GALPROP SLZ6R30T150C2 from
                            Ackermann+ 2012 (`totani.galprop`): ONE `gas`
                            template (pi0 + bremsstrahlung, ratio fixed) and
                            ONE `ics` template, as his sec 2.3 -- or three
                            ICS components with `ics_split` (his sec 3.4.4).
                            2 (or 4) templates in place of 14; `gas_rings`
                            splits gas into his five Galactocentric groups.
    Bubble, Loop I, isotropic and point sources are the same in both.

    `drop_components` removes named templates from the returned stack. Added
    2026-09-18 to test structured bubbles WITHOUT the two Loop I shells as a
    baseline candidate: the structured map (Ackermann+ 2017 Fig. 12, a
    1-10 GeV residual template) and the shells disagree at low energy -- with
    the shells present the map fits worse than the flat mask at 1.5-2.6 GeV
    (delta -2lnL = +243), without them better by 180-260. ⚠ The candidate was
    run and NOT adopted: the shells are not double-counted by the map (they
    fit 8-15% lower beside it, and buy -2,986 at 1.5 GeV), see
    `results/rebaseline_candidate_noLoopI_2026-09-18.md`. Unknown names raise,
    so a typo cannot silently keep a component.

    `run_dir` is a `gcereduce.reduce.reduce()` output directory (must contain
    `ccube.fits` and `expcube.fits`). `edges_gev` should be restricted to the
    ICS-covered range (bins 0-8, 1.162-129.985 GeV for this analysis) -- both
    `flat_spectrum_bin_exposure` and `source_expected_counts` raise if asked for
    energies outside the tabulated exposure, and `ics_on_target` has no
    coverage above its own ceiling regardless.
    """
    from . import GCE_FISHER_ROOT

    run_dir = Path(run_dir)
    root = GCE_FISHER_ROOT if root is None else root
    catalog_path = default_catalog_path() if catalog_path is None else catalog_path

    with fits.open(run_dir / "ccube.fits") as h:
        counts_cube = h[0].data.astype(float)
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float)
        exposure_mev = h["ENERGIES"].data["Energy"].astype(float)
        exp_header = _strip_to_2d(h[0].header)
    if (target["CDELT1"], target["CDELT2"], target["NAXIS1"], target["NAXIS2"]) != (
       exp_header["CDELT1"], exp_header["CDELT2"], exp_header["NAXIS1"], exp_header["NAXIS2"]):
        raise ValueError("ccube and expcube grids disagree; refusing to guess "
                         "which one is right")

    L, B = grid_lb(target)

    with fits.open(catalog_path) as h:
        cat = h["LAT_Point_Source_Catalog"].data
        extended = h["ExtendedSources"].data
    mask = extended_source_mask(L, B, extended)
    cat_l_all = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"])
    cat_b_all = cat["GLAT"]

    # 🚨 EVERY DIFFUSE TEMPLATE IS AN INTENSITY -- a value per sky DIRECTION --
    # so its counts contribution per pixel is value x exposure x SOLID ANGLE.
    # The pixel solid angle varies by a factor of TWO across |b| <= 60 (cos b),
    # so omitting it is a latitude-dependent SHAPE error no free per-bin
    # normalization can absorb. It WAS omitted in the first version, and was
    # found by the validation step (2026-09-04): the fitted isotropic jumped
    # five orders of magnitude between adjacent bins, while point sources --
    # the one component that is already an integrated flux and needs NO solid
    # angle -- fit at f = 0.64-0.98. See `templates.pixel_solid_angle`.
    omega = pixel_solid_angle(target)
    flat_basis = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev)
    single_basis = single_energy_bin_exposure(exposure_cube, exposure_mev, edges_gev)
    flat_basis_sr = flat_basis * omega[None, :, :]
    single_basis_sr = single_basis * omega[None, :, :]

    if diffuse not in ("a2020", "galprop_totani"):
        raise ValueError(f"diffuse must be 'a2020' or 'galprop_totani', got {diffuse!r}")
    templates = {}
    if diffuse == "galprop_totani":
        # Already bin-integrated under GALPROP's own spectrum: the IC-ring
        # route (single-energy exposure), not the flat-spectrum one.
        from .galprop import TOTANI_SET, galprop_on_target
        for name, cube in galprop_on_target(target, edges_gev, ics_split=ics_split,
                                            gas_rings=gas_rings,
                                            model_set=galprop_model_set or TOTANI_SET).items():
            templates[f"galprop_{name}"] = cube * single_basis_sr
    else:
        for name, shape_map in gas_dust_on_target(target, root, config, ts).items():
            templates[name] = shape_map[None, :, :] * flat_basis_sr
    # `bubble`: "ackermann" is the flat 2014 mask every result before
    # 2026-09-07 used and stays the default so they all reproduce; "macias" is
    # Oscar's structured map (see `templates.bubble_macias_on_target`). Same
    # component name either way -- they are alternative templates for ONE
    # emission, never both in the same stack.
    if bubble == "ackermann":
        bub = bubble_on_target(target, root)
    elif bubble == "macias":
        bub = bubble_macias_on_target(target)
    else:
        raise ValueError(f"bubble must be 'ackermann' or 'macias', got {bubble!r}")
    templates["bubble"] = bub[None, :, :] * flat_basis_sr
    for name, shape_map in loopi_on_target(target).items():
        templates[name] = shape_map[None, :, :] * flat_basis_sr
    templates["isotropic"] = isotropic_on_target(target)[None, :, :] * flat_basis_sr

    if diffuse == "a2020":
        six = ics_on_target(target, edges_gev, root)
        for name, cube in zip(("ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6"),
                              ics_four_ring(six)):
            templates[name] = cube * single_basis_sr

    n_bins = len(edges_gev) - 1
    counts_used = counts_cube[:n_bins]
    if counts_used.shape[0] != n_bins:
        raise ValueError(f"ccube has {counts_cube.shape[0]} bins, need {n_bins} "
                         f"for the requested edges")

    # 🚨 EVERYTHING ABOVE IS PIXEL-SPACE (960x960); bin it into cells NOW,
    # BEFORE adding the point-source template, which is built ALREADY in
    # cell-space (`point_source_cell_template` sums each source directly into
    # its one cell -- running it back through `bin_to_cells` would be binning
    # an already-binned array against the wrong grid shape).
    # 🚨 POINT-SOURCE MASKING (Totani sec 3.4.3's ALTERNATIVE treatment) is
    # PER BIN, because the 95% containment radius shrinks by 7x from 1.5 to
    # 100 GeV -- 59% of pixels masked in bin 0, 2% in bin 8. When it is active
    # the point-source TEMPLATE is dropped: a masked source's photons are gone
    # from the data, so modelling it too would double-count its absence. The
    # two treatments are alternatives, never combined.
    ps_masked = ps_mask_radii_deg is not None
    if ps_masked:
        radii = np.atleast_1d(np.asarray(ps_mask_radii_deg, float))
        if radii.size == 1:
            radii = np.full(n_bins, float(radii[0]))
        if radii.size != n_bins:
            raise ValueError(f"need one PS mask radius per bin ({n_bins}), "
                             f"got {radii.size}")
        # sources slightly outside the ROI can still mask pixels inside it
        near = (np.abs(cat_l_all) <= half_deg + 3.0) & (np.abs(cat_b_all) <= half_deg + 3.0)
        per_bin_mask = [mask | point_source_mask(L, B, cat_l_all[near],
                                                 cat_b_all[near], r)
                       for r in radii]
    else:
        per_bin_mask = [mask] * n_bins

    def _bin_cube(cube):
        return np.stack([bin_to_cells(cube[k], L, B, half_deg, cell_deg,
                                      mask=per_bin_mask[k])
                        for k in range(n_bins)])

    drop_components = tuple(drop_components)
    unknown = [n for n in drop_components if n not in templates]
    if unknown:
        raise ValueError(f"drop_components names {unknown} that are not in the "
                         f"stack: {sorted(templates)}")
    for n in drop_components:
        del templates[n]

    cell_templates = {name: _bin_cube(arr) for name, arr in templates.items()}
    cell_counts = _bin_cube(counts_used)

    cat_l, cat_b = cat_l_all, cat_b_all
    in_roi = (np.abs(cat_l) <= half_deg) & (np.abs(cat_b) <= half_deg)

    # 🚨 SOURCES INSIDE THE EXTENDED-SOURCE MASK MUST LEAVE THE TEMPLATE.
    # `mask` above excludes those pixels from the DATA and from every diffuse
    # template; a point source sitting in them has its photons excluded too, so
    # depositing its model counts is a pure model/data inconsistency.
    #
    # This was WRONG until 2026-09-04 and it was not small: 240 of the 2,651
    # ROI sources fall inside the mask circles, but they carry 29.5% (bin 0) to
    # 45.1% (bin 8) of the total predicted point-source flux -- bright sources
    # are exactly what gets catalogued as extended, so the masked minority
    # dominates the flux. That is 5-30x larger than the PSF-spill effect this
    # module already documents, and it is the leading suspect for the fitted
    # f_ps = 0.65-0.96 running below unity.
    #
    # Totani sec 2.3 is explicit that masked regions are "excluded from the
    # analysis" -- both sides of the comparison, not just the data. The
    # criterion here is `extended_source_mask` itself, evaluated at the SOURCE
    # POSITIONS rather than on the pixel grid, so the model and the data can
    # never disagree about where the mask is.
    in_mask = extended_source_mask(cat_l, cat_b, extended)
    keep_src = in_roi & ~in_mask
    if not ps_masked:
        cell_templates["point_sources"] = point_source_cell_template(
            cat[keep_src], cat_l[keep_src], cat_b[keep_src], edges_gev,
            exposure_cube, exposure_mev, exp_header, half_deg, cell_deg)

    return {"counts": cell_counts, "templates": cell_templates,
           "n_bins": n_bins, "n_masked_pixels": int(mask.sum()),
           "n_sources_kept": 0 if ps_masked else int(keep_src.sum()),
           "n_sources_dropped_in_mask": int((in_roi & in_mask).sum()),
           "ps_masked": bool(ps_masked), "bubble": bubble, "diffuse": diffuse,
           "galprop_model_set": galprop_model_set,
           "dropped_components": list(drop_components),
           "masked_pixel_fraction_by_bin": [float(m.mean()) for m in per_bin_mask]}


def halo_cell_template(run_dir, edges_gev, halo_params=None,
                       half_deg: float = 60.0, cell_deg: float = 10.0) -> np.ndarray:
    """The NFW rho^2 halo template, exposure-weighted and cell-binned -- the
    SAME pipeline `build_background_stack` applies to bubble/Loop I/isotropic
    (energy-independent spatial shapes, flat-assumed-spectrum exposure
    integral, extended-source masked), kept as a SEPARATE function.

    🚨 SEPARATE ON PURPOSE. `build_background_stack`'s own docstring states why
    it carries no halo template: importing and running it must never be the
    thing that crosses the rule 6 boundary. This function is what does --
    calling it is the deliberate, auditable step where a halo/DM template
    enters real-data machinery, gated by `validation/RULING_real_data.md` and
    `scripts/hooks/rule6_halo_guard.py`. `halo_params` overrides
    `halo.TOTANI_NFW_PARAMS` for a systematic variant (e.g. NFW-rho^1); the
    default reproduces Totani's baseline exactly.
    """
    from .halo import nfw_rho2_on_target

    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float)
        exposure_mev = h["ENERGIES"].data["Energy"].astype(float)

    J = nfw_rho2_on_target(target, **(halo_params or {}))
    # J is per-steradian (GeV^2 cm^-5 sr^-1): the same intensity-to-counts
    # solid-angle factor as every diffuse template. See build_background_stack.
    omega = pixel_solid_angle(target)
    basis = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev)
    pixel_cube = J[None, :, :] * basis * omega[None, :, :]

    L, B = grid_lb(target)
    with fits.open(default_catalog_path()) as h:
        extended = h["ExtendedSources"].data
    mask = extended_source_mask(L, B, extended)

    return bin_to_cells(pixel_cube, L, B, half_deg, cell_deg, mask=mask)


def cell_template_from_shape(run_dir, edges_gev, shape_map,
                             half_deg: float = 60.0, cell_deg: float = 10.0) -> np.ndarray:
    """ONE energy-independent spatial shape through the diffuse pipeline --
    flat-spectrum exposure x solid angle, extended-source mask, cells -- WITHOUT
    rebuilding the stack. Exactly the route `build_background_stack` takes for
    bubble / Loop I / isotropic and `halo_cell_template` takes for the halo, so
    a template added to a stack this way is on the same footing as the ones
    already in it.

    Used for systematic variants that ADD a component (2026-09-18: the region
    inside the eROSITA edge and outside the Fermi bubbles, `scripts/
    erosita_variant.py`). The shape is an intensity-like map on the target
    grid: a 0/1 footprint, a fractional-coverage map, or a structured map.
    """
    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float)
        exposure_mev = h["ENERGIES"].data["Energy"].astype(float)

    L, B = grid_lb(target)
    shape_map = np.asarray(shape_map, float)
    if shape_map.shape != L.shape:
        raise ValueError(f"shape map is {shape_map.shape}, target grid is {L.shape}")
    if not np.isfinite(shape_map).all():
        raise ValueError("shape map has non-finite values")

    omega = pixel_solid_angle(target)
    basis = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev)
    pixel_cube = shape_map[None, :, :] * basis * omega[None, :, :]

    with fits.open(default_catalog_path()) as h:
        extended = h["ExtendedSources"].data
    mask = extended_source_mask(L, B, extended)
    return bin_to_cells(pixel_cube, L, B, half_deg, cell_deg, mask=mask)
