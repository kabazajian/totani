"""Stage E: simulated skies through the structure selection (plan §5, Stage E).

A TRUTH is a real-data, background-only fit of one diffuse row plus a fixed
structure set, optionally with Totani's halo added on his fitted spectrum. A
SIMULATED SKY is a Poisson draw of that truth, made so that every step of the
analysis can run on it exactly as on the data -- including Totani's residual
bubble construction, which is built from PIXEL counts at 4.3 GeV:

  * the 120 analysis cells (|b| >= 10) of every bin are Poisson draws of the
    truth; bin 2 (4.3 GeV) is drawn per PIXEL and binned, so the residual
    construction sees a pixel map;
  * the 24 disk cells (|b| < 10), and the disk pixels at 4.3 GeV, carry the
    REAL counts (plus a Poisson draw of the injected halo there, when there is
    one, as `bubble_totani_injection` did). No analysis fit uses the disk; only
    the residual construction's 144-cell fit does, and it is dominated by the
    disk -- a truth fitted at |b| >= 10 would extrapolate a disk nobody fitted.

THE FAST RESIDUAL REBUILD. `selection.totani_residual_cells` costs ~47 s and
8 GB per row (GALPROP cubes from the tarball onto the 960x960 grid). Per
simulated sky that is impossible. But the only data-dependent parts are one
144-cell fit and the smoothed residual image, so the 4.3 GeV pixel templates
of each row are PRECOMPUTED once (`precompute_row`, float64 .npy, memory-mapped
so workers share the page cache), and a rebuild is then `fit_all_cells` +
`bubble_image` -- the same functions, on the same inputs -- plus a cached
shape-to-cells projection (`shape_to_cells`, equal to
`stack.cell_template_from_shape` without re-reading the FITS files per call).
`tests/test_simsky.py` checks the rebuild on the real counts reproduces the
cached residual templates.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import ranking
from .cells import bin_to_cells, halo_cell_mask

K_BUILD = 2                                         # 4.3 GeV, the bin the residual templates are built from
PIX = ranking.CACHE / "stage_e_pix2"                # gitignored with the rest of ranking_cells/
#: pixel templates every row shares (the residual construction's non-diffuse components)
SHARED_NAMES = ("bubble", "isotropic", "loopI_shell1", "loopI_shell2", "gce_nfw25", "point_sources")


# ---- pixel precompute -------------------------------------------------------------

def _row_dir(code: str, variant: str) -> Path:
    return PIX / f"{code}__{variant}"


def precompute_shared(edges_gev) -> None:
    """The model-independent 4.3 GeV pixel parts: shared templates, real counts,
    exposure basis, mask and grid. One process, once."""
    from .bubble_totani import _common_pixel_parts
    d = PIX / "_shared"
    if (d / "done").exists():
        return
    d.mkdir(parents=True, exist_ok=True)
    c = _common_pixel_parts(ranking.RUN, edges_gev, K_BUILD, True)
    for n in SHARED_NAMES:
        np.save(d / f"{n}.npy", np.asarray(c["shared"][n], float))
    for n in ("counts", "flat_sr", "L", "B"):
        np.save(d / f"{n}.npy", np.asarray(c[n], float))
    np.save(d / "mask.npy", np.asarray(c["mask"], bool))
    (d / "done").write_text("ok\n")


def precompute_row(code: str, variant: str, edges_gev) -> list:
    """One row's diffuse pixel templates at 4.3 GeV, exactly as
    `selection.totani_residual_cells` builds them. ~15 s, up to 8 GB: run one
    at a time. Returns the template names written."""
    from .bubble_totani import pixel_model_for_bin
    d = _row_dir(code, variant)
    if (d / "done").exists():
        return json.loads((d / "done").read_text())
    if code == "a2020":
        pm = pixel_model_for_bin(ranking.RUN, edges_gev, K_BUILD, diffuse="a2020", with_gce=True)
    else:
        target, _, _, _ = ranking._grid_and_basis(edges_gev)
        cubes = ranking.galprop_cubes(code, variant, target, edges_gev)      # Paper II tarball or Galp21
        pm = pixel_model_for_bin(ranking.RUN, edges_gev, K_BUILD, diffuse="galprop_totani", with_gce=True,
                                 galprop_cubes=cubes)
    d.mkdir(parents=True, exist_ok=True)
    names = sorted(n for n in pm["templates"] if n not in SHARED_NAMES)
    for n in names:
        np.save(d / f"{n}.npy", np.asarray(pm["templates"][n], float))
    (d / "done").write_text(json.dumps(names))
    return names


def load_shared() -> dict:
    d = PIX / "_shared"
    out = {n: np.load(d / f"{n}.npy", mmap_mode="r") for n in SHARED_NAMES + ("counts", "flat_sr", "L", "B")}
    out["mask"] = np.load(d / "mask.npy")
    return out


def load_row(code: str, variant: str) -> dict:
    d = _row_dir(code, variant)
    return {n: np.load(d / f"{n}.npy", mmap_mode="r") for n in json.loads((d / "done").read_text())}


# ---- the fast residual rebuild ------------------------------------------------------

_SHAPE: dict = {}


def _shape_basis(edges_gev) -> tuple:
    """(flat-spectrum exposure x solid angle per bin, L, B, extended-source
    mask), read once per process."""
    from .exposure import flat_spectrum_bin_exposure
    from .pointsources import default_catalog_path, extended_source_mask
    from .stack import _strip_to_2d
    from .templates import grid_lb, pixel_solid_angle

    key = tuple(np.asarray(edges_gev, float).tolist())
    if key not in _SHAPE:
        with fits.open(ranking.RUN / "ccube.fits") as h:
            target = _strip_to_2d(h[0].header)
        with fits.open(ranking.RUN / "expcube.fits") as h:
            exposure_cube = h[0].data.astype(float)
            exposure_mev = h["ENERGIES"].data["Energy"].astype(float)
        L, B = grid_lb(target)
        basis = flat_spectrum_bin_exposure(exposure_cube, exposure_mev, edges_gev) * pixel_solid_angle(target)[None]
        with fits.open(default_catalog_path()) as h:
            mask = extended_source_mask(L, B, h["ExtendedSources"].data)
        _SHAPE[key] = (basis, L, B, mask)
    return _SHAPE[key]


def shape_to_cells(shape_map, edges_gev) -> np.ndarray:
    """`stack.cell_template_from_shape(RUN, edges, shape)` without re-reading
    the FITS files per call."""
    basis, L, B, mask = _shape_basis(edges_gev)
    return bin_to_cells(np.asarray(shape_map, float)[None] * basis, L, B, mask=mask)


def residual_cells(counts_pix, row: dict, shared: dict, edges_gev) -> dict:
    """Totani's positive/negative residual templates for one row, from a 4.3 GeV
    pixel counts map, by `bubble_totani.fit_all_cells` + `bubble_image`."""
    from .bubble_totani import bubble_image, fit_all_cells
    T = {n: shared[n] for n in SHARED_NAMES}
    T.update(row)
    pm = {"templates": T, "counts": np.asarray(counts_pix, float), "flat_sr": shared["flat_sr"],
          "mask": shared["mask"], "L": shared["L"], "B": shared["B"]}
    img = bubble_image(pm, fit_all_cells(pm))
    return {"bubble_totani_pos": shape_to_cells(img["positive"], edges_gev),
            "bubble_totani_neg": shape_to_cells(img["negative"], edges_gev)}


# ---- truths -------------------------------------------------------------------------

def structure_pixels(name: str, edges_gev, shared: dict) -> np.ndarray:
    """A cell template's 4.3 GeV PIXEL map, by the route `ranking` bins it."""
    from . import structures
    from .templates import bubble_macias_on_target
    if name in ("isotropic", "point_sources"):
        return np.asarray(shared[name], float)
    if name == "bubble_flat":
        return np.asarray(shared["bubble"], float)
    target, _, _, single_sr = ranking._grid_and_basis(edges_gev)
    if name == "bubble_structured":
        return shape_to_basis(bubble_macias_on_target(target), edges_gev)
    arr = structures.STRUCTURES[name](target, edges_gev)
    if arr.ndim == 2:
        return shape_to_basis(arr, edges_gev)
    return arr[K_BUILD] * single_sr[K_BUILD]


def shape_to_basis(shape_map, edges_gev) -> np.ndarray:
    """An energy-independent shape's 4.3 GeV pixel counts per unit amplitude."""
    return np.asarray(shape_map, float) * _shape_basis(edges_gev)[0][K_BUILD]


def halo_amplitudes_totani(path) -> np.ndarray:
    """Per-bin NFW-rho^2 amplitudes of Totani's fitted spectrum (his
    configuration, `bubble_totani_injection_2026-09-19`), negative bins set
    to zero -- the 'totani_like' spectrum of that test."""
    z = json.loads(Path(path).read_text())
    counts = np.asarray(z["a0_halo_counts"], float)
    return np.clip(counts / np.asarray(z["roi_J"], float), 0.0, None)


def build_truth(model: str, variant: str, structures: tuple, edges_gev, halo_amp=None) -> dict:
    """Expected counts of a truth: the real-data background-only fit (|b| >= 10)
    of `model` + `structures`, as (n_bins, 12, 12) cells and a 4.3 GeV pixel
    map, plus the halo (cells and pixels) if `halo_amp` (per bin) is given."""
    from . import selection
    from .halo import nfw_rho2_on_target
    from .likelihood import fit_all_bins, restrict_stack

    sh = ranking.shared_cells(edges_gev)
    shared = load_shared()
    st = selection.base_stack(model, variant, edges_gev, sh)
    for c in structures:
        st, _, signs = selection.add_candidate(st, c, model, variant, edges_gev, sh)
        if signs is not None:
            raise ValueError(f"truth structure {c} is sign-free; not supported")
    keep = halo_cell_mask()
    fit = fit_all_bins(restrict_stack(st, keep), halo_mask=keep)
    names = fit["component_order"]
    amps = np.array([[fit["per_bin"][k]["amplitudes"][n] for n in names] for k in range(st["n_bins"])])
    cells = np.einsum("kn,nkij->kij", amps, np.stack([st["templates"][n] for n in names]))
    row = load_row("a2020" if model == "a2020" else model, variant)
    pix = np.zeros(shared["L"].shape)
    for i, n in enumerate(names):
        p = row[n] if n in row else structure_pixels(n, edges_gev, shared)
        pix += amps[K_BUILD, i] * p
    out = {"model": model, "variant": variant, "structures": list(structures), "names": names,
           "amplitudes": amps, "cells": cells, "pix2": pix, "halo_cells": np.zeros_like(cells),
           "halo_pix2": np.zeros_like(pix)}
    if halo_amp is not None:
        target, _, _, _ = ranking._grid_and_basis(edges_gev)
        halo_amp = np.asarray(halo_amp, float)
        out["halo_cells"] = halo_amp[:, None, None] * sh["J"]
        out["halo_pix2"] = halo_amp[K_BUILD] * nfw_rho2_on_target(target) * np.asarray(shared["flat_sr"])
        out["halo_amp"] = halo_amp
    return out


def disk_pixels(shared: dict) -> np.ndarray:
    """Pixels in the 24 disk cells (|b| < 10), by the cell assignment
    `bin_to_cells` uses (pixel centres; the boundary falls between pixels)."""
    return np.abs(np.asarray(shared["B"])) < 10.0


def draw(truth: dict, real_counts_cells, shared: dict, rng) -> tuple:
    """One simulated sky: (cells (n_bins, 12, 12), 4.3 GeV pixel counts)."""
    keep = halo_cell_mask()
    mu = truth["cells"] + truth["halo_cells"]
    cells = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), real_counts_cells).astype(float)
    disk = disk_pixels(shared)
    real_pix = np.asarray(shared["counts"])
    pix = np.where(disk, real_pix + rng.poisson(np.clip(truth["halo_pix2"], 0, None)),
                   rng.poisson(np.clip(truth["pix2"] + truth["halo_pix2"], 0, None))).astype(float)
    cells[K_BUILD] = bin_to_cells(pix, shared["L"], shared["B"], mask=shared["mask"])
    return cells, pix
