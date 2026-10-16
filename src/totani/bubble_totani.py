"""Totani's Fermi-bubble template, built the way his sec 3.1 says.

Kev, 2026-09-19: "Let's build that bubble template as he did." His recipe
(arXiv:2507.07209 sec 3.1, read from `docs/`):

  1. A flat template with sharp edges "within a rough area of the Fermi
     bubbles" is fitted to the data with every other component -- INCLUDING
     the GC GeV excess as a known component (NFW-rho^2.5, his sec 2.3) -- over
     the WHOLE ROI, disk included, at 1.5 and 4.3 GeV.
  2. The "best-fit image of the bubble" is the flat template plus the fit
     residuals, per pixel, converted to flux with the exposure and smoothed
     with a Gaussian of sigma = 1 deg. He iterates the flat boundary by eye
     against these images; here the boundary is the Ackermann+ 2014 flat mask
     and is not iterated (stated wherever the result is used).
  3. From the 4.3 GeV image: the POSITIVE-flux regions, wherever they are,
     become one energy-independent template; the NEGATIVE-flux regions
     another. Both are free per bin in the halo fits, in place of any bubble
     template ("the positive and negative residual templates").

His own words on what that does: "With this treatment, the model and data
maps should agree well without additional model components at around
4.3 GeV. Even if a halo component is present near this energy bin, it will be
absorbed by these templates." So a halo fit at 4.3 GeV must come out ~0 by
construction -- that is the check that the templates were built right -- and
the search is for the energy DEPENDENCE relative to 4.3 GeV.

🚨 THE RESIDUAL IS PER PIXEL BUT THE FIT IS PER CELL, in his analysis and
here. The normalisations f_l come from the 144-cell fit; the pixel residual
is counts(pixel) - sum_l f_l T_l(pixel). Point sources are deposited at
their own pixel with no PSF (this pipeline has no PSF-convolved source map);
at 4.3 GeV the 68% PSF is ~0.2 deg against the 1 deg smoothing, so a source's
deposit and its true PSF-spread counts smooth to nearly the same thing, but
the brightest sources will leave a positive ring / negative centre pattern
at the percent level -- which is also what his figure 1 shows ("residuals of
bright point sources are visible in some places").

Everything else -- diffuse model, Loop I shells, isotropic, catalogue,
extended-source mask, cells -- is the same machinery as the rest of the
project, selected by the same arguments.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

from .cells import bin_to_cells
from .exposure import flat_spectrum_bin_exposure, single_energy_bin_exposure
from .halo import TOTANI_NFW_PARAMS
from .likelihood import fit_bin
from .pointsources import (default_catalog_path, exposure_pixel, extended_source_mask,
                           source_expected_counts)
from .stack import _strip_to_2d
from .templates import (axis_coords, bubble_on_target, gas_dust_on_target, grid_lb,
                        ics_four_ring, ics_on_target, isotropic_on_target,
                        loopi_on_target, pixel_solid_angle)


def nfw_rho_power_on_target(target, power: float = 2.5, s_max: float = 200.0,
                            n_s: int = 800, row_chunk: int = 80) -> np.ndarray:
    """int rho^power ds with Totani's NFW parameters. power = 2.5 is his GC
    GeV excess template (sec 2.3: "epsilon ∝ rho^2.5 with the original NFW",
    emissivity ∝ r^-2.5 at the centre). Same row-chunked line-of-sight
    integral as `halo.nfw_rho2_on_target`, which is the power = 2 case."""
    from gcereduce.dm_templates import A2020_DEFAULTS, KPC_CM, _los_nodes, rho_gnfw

    p = dict(A2020_DEFAULTS)
    p.update({k: v for k, v in TOTANI_NFW_PARAMS.items() if k != "profile"})
    l_c, b_c = axis_coords(target)
    out = np.empty((len(b_c), len(l_c)))
    for r0 in range(0, len(b_c), row_chunk):
        r1 = min(r0 + row_chunk, len(b_c))
        r, jac, t = _los_nodes(l_c, b_c[r0:r1], p["R_sun"], s_max, n_s, p["q_b"], p["q_c"], xp=np)
        rho = rho_gnfw(r, p["gamma"], p["rs"], p["rho_sun"], p["R_sun"], xp=np)
        out[r0:r1] = np.trapezoid(rho ** power * jac, t, axis=-1) * KPC_CM
    return out


def gaussian_smooth_sky(image, weight, B, sigma_deg: float = 1.0, pixel_deg: float = 0.125):
    """Masked Gaussian smoothing on a plate-carree grid: sigma = `sigma_deg`
    on the sky, so along l the pixel sigma is sigma/(pixel cos b), per row.
    `weight` is 1 where the pixel is valid, 0 where masked; the result is
    G(w f) / G(w), so masked pixels neither contribute nor dilute."""
    from scipy.ndimage import gaussian_filter1d

    s_pix = sigma_deg / pixel_deg
    wf = image * weight
    w = weight.astype(float)
    # along b (columns): constant sigma
    wf = gaussian_filter1d(wf, s_pix, axis=0, mode="nearest")
    w = gaussian_filter1d(w, s_pix, axis=0, mode="nearest")
    # along l (rows): sigma / cos b, row by row
    out_f = np.empty_like(wf); out_w = np.empty_like(w)
    for i in range(image.shape[0]):
        s = s_pix / max(np.cos(np.radians(B[i, 0])), 0.2)
        out_f[i] = gaussian_filter1d(wf[i], s, mode="nearest")
        out_w[i] = gaussian_filter1d(w[i], s, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        sm = np.where(out_w > 1e-3, out_f / out_w, 0.0)
    return sm


#: In-process cache of everything in `pixel_model_for_bin` that does NOT
#: depend on the diffuse model: the counts, the exposure bases, the mask, the
#: flat bubble, Loop I, isotropic, the GC-excess profile and the point-source
#: deposit. Building the residual templates for the 128-model grid calls this
#: 256 times with identical values for all of it, and the rho^2.5
#: line-of-sight integral alone is seconds of that. Keyed by everything that
#: could change it. Arrays are shared, not copied -- no caller mutates them.
_COMMON: dict = {}


def _common_pixel_parts(run_dir, edges_gev, k: int, with_gce: bool):
    key = (str(run_dir), tuple(np.asarray(edges_gev, float).tolist()), int(k), bool(with_gce))
    if key in _COMMON:
        return _COMMON[key]
    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        counts = h[0].data[k].astype(float)
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float)
        exposure_mev = h["ENERGIES"].data["Energy"].astype(float)
        exp_header = _strip_to_2d(h[0].header)
    L, B = grid_lb(target)
    omega = pixel_solid_angle(target)
    flat_sr = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev)[k] * omega
    single_sr = single_energy_bin_exposure(exposure_cube, exposure_mev, edges_gev)[k] * omega

    with fits.open(default_catalog_path()) as h:
        cat = h["LAT_Point_Source_Catalog"].data
        extended = h["ExtendedSources"].data
    mask = extended_source_mask(L, B, extended)

    shared = {"bubble": (bubble_on_target(target) > 0).astype(float) * flat_sr,
              "isotropic": isotropic_on_target(target) * flat_sr}
    for name, m in loopi_on_target(target).items():
        shared[name] = m * flat_sr
    if with_gce:
        shared["gce_nfw25"] = nfw_rho_power_on_target(target, 2.5) * flat_sr

    cat_l = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"]); cat_b = cat["GLAT"]
    in_roi = (np.abs(cat_l) <= 60.0) & (np.abs(cat_b) <= 60.0)
    keep = in_roi & ~extended_source_mask(cat_l, cat_b, extended)
    sc = source_expected_counts(cat[keep], cat_l[keep], cat_b[keep], edges_gev,
                                exposure_cube, exposure_mev, exp_header)[:, k]
    row, col = exposure_pixel(cat_l[keep], cat_b[keep], target)
    ps = np.zeros(L.shape)
    np.add.at(ps, (row, col), sc)
    shared["point_sources"] = ps

    out = {"shared": shared, "counts": counts, "flat_sr": flat_sr, "single_sr": single_sr,
           "mask": mask, "L": L, "B": B, "target": target, "n_sources": int(keep.sum())}
    _COMMON[key] = out
    return out


def pixel_model_for_bin(run_dir, edges_gev, k: int, *, diffuse: str = "galprop_totani",
                        with_gce: bool = True, galprop_model_set=None, galprop_cubes=None) -> dict:
    """Per-pixel counts-basis maps of every component at analysis bin `k`,
    the counts map, the exposure basis and the extended-source mask -- the
    same construction as `stack.build_background_stack`, one bin, kept in
    pixel space. The bubble is the sharp-edged flat mask.

    The model-independent half is cached in `_COMMON` (see there)."""
    c = _common_pixel_parts(run_dir, edges_gev, k, with_gce)
    target, L, B = c["target"], c["L"], c["B"]
    flat_sr, single_sr, mask, counts = c["flat_sr"], c["single_sr"], c["mask"], c["counts"]

    T = {}
    if diffuse == "galprop_totani":
        from .galprop import TOTANI_SET, galprop_on_target
        # `galprop_cubes`: precomputed {name: (n_bins, ny, nx)} (the ranking
        # reads them from a tarball); otherwise built from the extracted files.
        cubes = galprop_cubes if galprop_cubes is not None else galprop_on_target(
            target, edges_gev, model_set=galprop_model_set or TOTANI_SET)
        for name, cube in cubes.items():
            T[f"galprop_{name}"] = cube[k] * single_sr
    elif diffuse == "a2020":
        from . import GCE_FISHER_ROOT
        for name, m in gas_dust_on_target(target, GCE_FISHER_ROOT).items():
            T[name] = m * flat_sr
        six = ics_on_target(target, edges_gev, GCE_FISHER_ROOT)
        for name, cube in zip(("ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6"), ics_four_ring(six)):
            T[name] = cube[k] * single_sr
    else:
        raise ValueError(diffuse)
    T.update(c["shared"])      # flat bubble, Loop I, isotropic, GCE, point sources

    return {"templates": T, "counts": counts, "flat_sr": flat_sr, "mask": mask,
            "L": L, "B": B, "target": target, "n_sources": c["n_sources"]}


def fit_all_cells(pm: dict) -> dict:
    """The 144-cell fit of every component at this bin, disk included."""
    names = sorted(pm["templates"])
    L, B, mask = pm["L"], pm["B"], pm["mask"]
    T = np.stack([bin_to_cells(pm["templates"][n], L, B, mask=mask) for n in names])
    n = bin_to_cells(pm["counts"], L, B, mask=mask)
    val, amp, r = fit_bin(T, n)
    return {"names": names, "amplitudes": dict(zip(names, amp)), "neg2lnl": float(val),
            "converged": bool(r.status == 0)}


def bubble_image(pm: dict, fit: dict, sigma_deg: float = 1.0) -> dict:
    """Flat template + residual, per pixel, in flux per bin per sr, smoothed.
    = counts - sum_{l != bubble} f_l T_l, divided by exposure x solid angle."""
    model_other = sum(fit["amplitudes"][n] * pm["templates"][n] for n in fit["names"] if n != "bubble")
    resid_counts = pm["counts"] - model_other
    valid = ~pm["mask"] & (pm["flat_sr"] > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        flux = np.where(valid, resid_counts / pm["flat_sr"], 0.0)
    smooth = gaussian_smooth_sky(flux, valid.astype(float), pm["B"], sigma_deg)
    smooth[~valid] = 0.0
    return {"flux_raw": flux, "flux_smooth": smooth, "valid": valid,
            "positive": np.clip(smooth, 0, None), "negative": np.clip(-smooth, 0, None),
            "bubble_amplitude": fit["amplitudes"]["bubble"]}
