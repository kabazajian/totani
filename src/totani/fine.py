"""The analysis on finer cells (default 2 deg): a cached template library and
row stacks, and the selection/scoring primitives, independent of the 10-deg
cache that ranking/selection use.

Why a separate module: every 10-deg cache in `ranking` is keyed by name
only, and the 10-deg selection code reads them implicitly; running the 2-deg
programme through the same paths would mix binnings. Here the library for
one cell size lives in `results/fits/ranking_cells/cells<deg>/`, built once
from the pixel maps (`finer_cells_test.pixel_library`, extended to every
Stage C candidate), and each row's diffuse templates and Totani residual pair
(built on his 144 cells of 10 deg, projected) are cached beside it.

BIC at this binning uses n = (number of |b|>=10 cells) x 9 bins -- 27,000 at
2 deg, not 1,080 -- so the per-template threshold is 9 ln 27000 = 92, not 63.
Expected counts are exposure-weighted templates without PSF convolution.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import ranking, selection
from .likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack

K20 = 5
STRUCTURE_NAMES = sorted({n for names in selection.CANDIDATES.values() for n in names
                          if not n.startswith("bubble_totani_")})


def cache_dir(cell_deg: float) -> Path:
    d = ranking.CACHE / f"cells{cell_deg:g}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def build_library(cell_deg: float, edges_gev) -> Path:
    """All model-independent templates at this cell size, from the pixel maps."""
    import sys
    sys.path.insert(0, str(ranking.ROOT / "scripts"))
    from finer_cells_test import pixel_library
    from . import structures
    from .cells import bin_to_cells
    from .stack import cell_template_from_shape

    p = cache_dir(cell_deg) / "library.npz"
    if p.exists():
        return p
    lib = pixel_library(cell_deg)
    target, L, B, single_sr, mask = lib["grid"]
    T = dict(lib["templates"])
    for n in STRUCTURE_NAMES:
        if n in T:
            continue
        arr = structures.STRUCTURES[n](target, edges_gev)
        if arr.ndim == 2:
            T[n] = cell_template_from_shape(ranking.RUN, edges_gev, arr, cell_deg=cell_deg)
        else:
            pix = arr * single_sr
            T[n] = np.stack([bin_to_cells(pix[k], L, B, cell_deg=cell_deg, mask=mask) for k in range(9)])
    missing = [n for n in STRUCTURE_NAMES if n not in T]
    if missing:
        raise KeyError(f"structures not built at {cell_deg} deg: {missing}")
    q = lib["quads"]
    np.savez(p, counts=lib["counts"], J=lib["J"], keep=lib["keep"],
             **{f"T__{n}": v for n, v in T.items()},
             **{f"A__{n}": v for n, v in lib["a2020"].items()},
             **{f"Q__{k}": v for k, v in q.items()})
    return p


def load_library(cell_deg: float) -> dict:
    z = np.load(cache_dir(cell_deg) / "library.npz")
    lib = {"counts": z["counts"], "J": z["J"], "keep": z["keep"].astype(bool), "templates": {}, "a2020": {}, "quads": {}}
    for k in z.files:
        if k.startswith("T__"):
            lib["templates"][k[3:]] = z[k]
        elif k.startswith("A__"):
            lib["a2020"][k[3:]] = z[k]
        elif k.startswith("Q__"):
            lib["quads"][k[3:]] = z[k].astype(bool)
    lib["n_data"] = int(lib["keep"].sum()) * 9
    return lib


def row_path(cell_deg, code, variant) -> Path:
    return cache_dir(cell_deg) / f"row__{code}__{variant}.npz"


def build_row(cell_deg, code, variant, edges_gev) -> Path:
    """One row's diffuse templates and residual pair at this cell size (pixel level, ~8 GB)."""
    import sys
    sys.path.insert(0, str(ranking.ROOT / "scripts"))
    from finer_cells_test import pixel_library, row_templates
    p = row_path(cell_deg, code, variant)
    if p.exists():
        return p
    lib = pixel_library(cell_deg)
    diffuse, res = row_templates((code, variant), cell_deg, lib)
    np.savez(p, **{f"D__{n}": v for n, v in diffuse.items()}, **res)
    return p


def load_row(cell_deg, code, variant) -> dict:
    z = np.load(row_path(cell_deg, code, variant))
    return {"diffuse": {k[3:]: z[k] for k in z.files if k.startswith("D__")},
            "residual": {k: z[k] for k in ("bubble_totani_pos", "bubble_totani_neg")}}


def base_stack(lib, row, counts=None) -> dict:
    t = {"isotropic": lib["templates"]["isotropic"], "point_sources": lib["templates"]["point_sources"], **row["diffuse"]}
    return {"counts": lib["counts"] if counts is None else counts, "n_bins": 9, "templates": t}


def template(lib, row, name):
    return row["residual"][name] if name.startswith("bubble_totani_") else lib["templates"][name]


def bg_neg2lnl(stack, keep) -> np.ndarray:
    return np.array([b["neg2lnl"] for b in fit_all_bins(restrict_stack(stack, keep), halo_mask=keep)["per_bin"]])


def add_candidate(stack, cand, lib, row, keep):
    """As selection.add_candidate: a sign-free template gets the better sign per bin."""
    new = {n: template(lib, row, n) for n in selection.CANDIDATES[cand]}
    free = [n for n in new if n in selection.SIGN_FREE]
    if not free:
        st = {**stack, "templates": {**stack["templates"], **new}}
        return st, bg_neg2lnl(st, keep)
    nf = free[0]; fixed = {n: v for n, v in new.items() if n != nf}
    plus = {**stack, "templates": {**stack["templates"], **fixed, nf: new[nf]}}
    minus = {**stack, "templates": {**stack["templates"], **fixed, nf: -new[nf]}}
    vp, vm = bg_neg2lnl(plus, keep), bg_neg2lnl(minus, keep)
    signs = np.where(vp <= vm, 1, -1)
    st = {**stack, "templates": {**stack["templates"], **fixed, nf: new[nf] * signs[:, None, None]}}
    return st, np.minimum(vp, vm)


def bic_threshold(n_params: int, n_data: int) -> float:
    return n_params * float(np.log(n_data))


def score(stack, lib, with_sigmav=True, edges_gev=None) -> dict:
    """TS per bin, the 20.8 GeV amplitude with interval, the spectrum, the
    quadrant gate, and optionally the joint b-bbar 500 GeV <sigma v>."""
    keep, J = lib["keep"], lib["J"]
    sr, Jr = restrict_stack(stack, keep, J)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=True)["per_bin"]
    roi = [float(Jr[k][keep].sum()) for k in range(9)]
    qa = {}
    for q, m in lib["quads"].items():
        stq, Jq = restrict_stack(stack, m, J)
        one = {"counts": stq["counts"][K20:K20 + 1], "n_bins": 1, "templates": {n: t[K20:K20 + 1] for n, t in stq["templates"].items()}}
        qa[q] = fit_all_bins_profile(one, Jq[K20:K20 + 1], halo_mask=m, with_uncertainty=False)["per_bin"][0]["signal_amplitude"] * float(Jq[K20][m].sum())
    ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(9)]
    out = {"neg2lnl_H0": float(sum(b["neg2lnl"] for b in bg)), "ts": ts,
           "halo_counts": [pf[k]["signal_amplitude"] * roi[k] for k in range(9)],
           "halo20_lo": pf[K20]["signal_amplitude_lo"] * roi[K20], "halo20_hi": pf[K20]["signal_amplitude_hi"] * roi[K20],
           "quadrants20": qa, "n_free_per_bin": len(stack["templates"]),
           "favor": bool(ts[K20] > 4 and pf[K20]["signal_amplitude"] > 0 and all(v > 0 for v in qa.values()))}
    if with_sigmav:
        from . import dmspectra
        sv = dmspectra.joint_fit(stack, J, "b", 500, edges_gev, halo_mask=keep, n_scan=31)
        out["sigmav"] = {k: sv[k] for k in ("sigmav_best", "sigmav_lo68", "sigmav_hi68", "ts_vs_zero")}
    return out


# ---- PSF-folded library (2026-09-27) --------------------------------------------------
# Every template's pixel map is convolved with the LAT PSF per energy bin
# (totani.psf) before binning to cells. The residual construction's pair is
# NOT convolved: it is built from the counts, which already carry the PSF.
# Selected with FINE_TAG=psf in the environment; caches in cells<deg>_psf/.

import os as _os


def tag() -> str:
    t = _os.environ.get("FINE_TAG", "")
    return f"_{t}" if t else ""


_cache_dir_plain = cache_dir


def cache_dir(cell_deg: float) -> Path:                     # noqa: F811 -- tag-aware
    d = ranking.CACHE / f"cells{cell_deg:g}{tag()}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _psf_grid(edges_gev):
    from .exposure import flat_spectrum_bin_exposure
    from .pointsources import default_catalog_path, extended_source_mask
    from .stack import _strip_to_2d
    from .templates import grid_lb, pixel_solid_angle
    from astropy.io import fits
    with fits.open(ranking.RUN / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(ranking.RUN / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float); exposure_mev = h["ENERGIES"].data["Energy"].astype(float)
        exp_header = _strip_to_2d(h[0].header)
    L, B = grid_lb(target); om = pixel_solid_angle(target)
    basis = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev) * om[None]
    _, _, _, single_sr = ranking._grid_and_basis(edges_gev)
    with fits.open(default_catalog_path()) as h:
        mask = extended_source_mask(L, B, h["ExtendedSources"].data)
    return {"target": target, "L": L, "B": B, "basis": basis, "single_sr": single_sr, "mask": mask,
            "exposure_cube": exposure_cube, "exposure_mev": exposure_mev, "exp_header": exp_header}


def _cells_psf(cube_pix, g, cell_deg, edges_gev):
    from . import psf
    from .cells import bin_to_cells
    conv = psf.convolve_cube(np.asarray(cube_pix, float), g["B"], edges_gev)
    return np.stack([bin_to_cells(conv[k], g["L"], g["B"], cell_deg=cell_deg, mask=g["mask"]) for k in range(conv.shape[0])])


def point_sources_cells(g, cell_deg, edges_gev, exact: bool = True):
    """The PSF-folded point-source template on cells. `exact` (2026-10-01, the
    production choice from that date): each source at its catalog position with
    the bin's kernel integrated over the pixels (psf.point_sources_exact), which
    the gtsrcmaps benchmark showed agrees with Fermitools to <= 2% of a source's
    counts per cell. `exact=False`: the pre-10-01 placement, each source
    deposited in the pixel containing it and the map convolved (1-14% of a
    source's counts in a different cell than Fermitools)."""
    from . import psf
    from .cells import bin_to_cells
    from .pointsources import default_catalog_path, exposure_pixel, extended_source_mask, source_expected_counts
    from astropy.io import fits
    with fits.open(default_catalog_path()) as h:
        cat = h["LAT_Point_Source_Catalog"].data; extended = h["ExtendedSources"].data
    cl = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"]); cb = cat["GLAT"]
    keep = (np.abs(cl) <= 60.0) & (np.abs(cb) <= 60.0) & ~extended_source_mask(cl, cb, extended)
    sc = source_expected_counts(cat[keep], cl[keep], cb[keep], edges_gev, g["exposure_cube"], g["exposure_mev"], g["exp_header"])
    if exact:
        conv = psf.point_sources_exact(cl[keep], cb[keep], sc, g["L"], g["B"], edges_gev)
        return np.stack([bin_to_cells(conv[k], g["L"], g["B"], cell_deg=cell_deg, mask=g["mask"]) for k in range(conv.shape[0])])
    row, col = exposure_pixel(cl[keep], cb[keep], g["target"])
    ps = np.zeros((len(edges_gev) - 1,) + g["L"].shape)
    for k in range(ps.shape[0]):
        np.add.at(ps[k], (row, col), sc[:, k])
    return _cells_psf(ps, g, cell_deg, edges_gev)


def build_library_psf(cell_deg: float, edges_gev) -> Path:
    """The model-independent templates, PSF-folded, at this cell size. Counts and
    masks are taken from the unconvolved library (the data are the data)."""
    from . import structures
    from .halo import nfw_rho2_on_target
    from .pointsources import (default_catalog_path, exposure_pixel, extended_source_mask,
                               source_expected_counts)
    from .templates import bubble_macias_on_target, bubble_on_target, isotropic_on_target, loopi_on_target
    from astropy.io import fits
    p = cache_dir(cell_deg) / "library.npz"
    if p.exists():
        return p
    plain = np.load(ranking.CACHE / f"cells{cell_deg:g}" / "library.npz")
    g = _psf_grid(edges_gev); target = g["target"]
    T = {}
    shapes = {"bubble_flat": (bubble_on_target(target) > 0).astype(float), "bubble_structured": bubble_macias_on_target(target),
              "isotropic": isotropic_on_target(target), **loopi_on_target(target)}
    for n, shp in shapes.items():
        T[n] = _cells_psf(shp[None] * g["basis"], g, cell_deg, edges_gev)
    for n in STRUCTURE_NAMES:
        if n in T:
            continue
        arr = structures.STRUCTURES[n](target, edges_gev)
        cube = arr[None] * g["basis"] if arr.ndim == 2 else arr * g["single_sr"]
        T[n] = _cells_psf(cube, g, cell_deg, edges_gev)
    T["point_sources"] = point_sources_cells(g, cell_deg, edges_gev)
    J = _cells_psf(nfw_rho2_on_target(target)[None] * g["basis"], g, cell_deg, edges_gev)
    # A2020 diffuse, PSF-folded
    from . import GCE_FISHER_ROOT
    from .templates import gas_dust_on_target, ics_four_ring, ics_on_target
    A = {n: _cells_psf(m[None] * g["basis"], g, cell_deg, edges_gev) for n, m in gas_dust_on_target(target, GCE_FISHER_ROOT).items()}
    six = ics_on_target(target, edges_gev, GCE_FISHER_ROOT)
    for n, cube in zip(("ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6"), ics_four_ring(six)):
        A[n] = _cells_psf(cube * g["single_sr"], g, cell_deg, edges_gev)
    missing = [n for n in STRUCTURE_NAMES if n not in T]
    if missing:
        raise KeyError(f"not built: {missing}")
    np.savez(p, counts=plain["counts"], J=J, keep=plain["keep"], **{f"T__{n}": v for n, v in T.items()},
             **{f"A__{n}": v for n, v in A.items()}, **{k: plain[k] for k in plain.files if k.startswith("Q__")})
    return p


def build_row_psf(cell_deg, code, variant, edges_gev) -> Path:
    """One row's diffuse templates PSF-folded; the residual pair copied from the
    unconvolved row cache (built from counts)."""
    p = row_path(cell_deg, code, variant)
    if p.exists():
        return p
    plain = np.load(ranking.CACHE / f"cells{cell_deg:g}" / f"row__{code}__{variant}.npz")
    if code == "a2020":
        lib = np.load(cache_dir(cell_deg) / "library.npz")
        diffuse = {k[3:]: lib[k] for k in lib.files if k.startswith("A__")}
    else:
        g = _psf_grid(edges_gev)
        cubes = ranking.galprop_cubes(code, variant, g["target"], edges_gev)
        diffuse = {f"galprop_{n}": _cells_psf(c * g["single_sr"], g, cell_deg, edges_gev) for n, c in cubes.items()}
    np.savez(p, **{f"D__{n}": v for n, v in diffuse.items()},
             bubble_totani_pos=plain["bubble_totani_pos"], bubble_totani_neg=plain["bubble_totani_neg"])
    return p
