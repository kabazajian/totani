"""Pixels to 10-degree cells — the binning the whole analysis is cheap because of.

Totani sec 2.2 evaluates the likelihood not per pixel but on 12 x 12 cells of
10 deg, because "the halo-like component is expected to be much smoother than
the pixel scale of 0.125 deg". That turns a 921,600-pixel comparison into a
144-cell one, and it is why this ROI is tractable at all.

He also gives the second reason, which is not about speed: "smoothing at coarse
resolutions may reduce the impact of pixel-scale mismatches between a model map
and data on the likelihood estimate". The gas templates have real pixel-scale
structure that the model will never match exactly; cells stop that structure
from driving the fit.

⚠ **AND IT IS ALSO A COST, WHICH HE STATES.** Fine structure is exactly what
would otherwise SEPARATE the gas from the smooth components -- "some components
of the Galactic diffuse emission, such as those from interstellar gas, possess
pixel-scale fine structures, which may help separate them from others by
pixel-scale likelihood calculation". So cells buy robustness against
mismodelling at the price of discriminating power between templates. That
trade is his design and we reproduce it; it is not free and should not be
described as though it were.
"""
from __future__ import annotations

import numpy as np

from .pointsources import cell_edges

__all__ = ["bin_to_cells", "cell_centers", "halo_cell_mask"]


def bin_to_cells(data, l_grid, b_grid, half_deg: float = 60.0,
                 cell_deg: float = 10.0, mask=None) -> np.ndarray:
    """Sum `data` into (n_cell, n_cell) cells; leading axes are preserved.

    `data` is a 2-D map or an (n_bins, ny, nx) cube. `mask` is a boolean array
    of pixels to EXCLUDE (extended-source regions, per Totani sec 2.3).

    Summing is right for counts and for exposure-weighted model counts alike --
    both are extensive. It would be wrong for an intensity, which must be
    weighted by solid angle first; nothing here passes one.
    """
    data = np.asarray(data, float)
    l_grid = np.asarray(l_grid, float)
    b_grid = np.asarray(b_grid, float)
    if data.shape[-2:] != l_grid.shape:
        raise ValueError(f"data {data.shape} does not match grid {l_grid.shape}")

    edges = cell_edges(half_deg, cell_deg)
    n = len(edges) - 1
    # np.digitize on the pixel CENTRES; a pixel belongs to exactly one cell.
    col = np.digitize(l_grid, edges) - 1
    row = np.digitize(b_grid, edges) - 1
    good = (col >= 0) & (col < n) & (row >= 0) & (row < n)
    if mask is not None:
        good &= ~np.asarray(mask, bool)

    flat = (row * n + col)[good]
    cube = data.reshape(-1, *l_grid.shape)
    out = np.empty((cube.shape[0], n * n))
    for k in range(cube.shape[0]):
        out[k] = np.bincount(flat, weights=cube[k][good], minlength=n * n)
    out = out.reshape(*data.shape[:-2], n, n) if data.ndim > 2 else out.reshape(n, n)
    return out


def cell_centers(half_deg: float = 60.0, cell_deg: float = 10.0):
    """(l, b) of each cell centre, as (n_cell, n_cell) arrays."""
    e = cell_edges(half_deg, cell_deg)
    c = 0.5 * (e[:-1] + e[1:])
    return np.meshgrid(c, c)          # (L, B): rows are latitude


def halo_cell_mask(b_min: float = 10.0, half_deg: float = 60.0,
                   cell_deg: float = 10.0) -> np.ndarray:
    """Cells kept for the halo search: |b| >= `b_min`.

    With 10 deg cells and b_min = 10 the split is exact -- no cell straddles the
    boundary -- which is why Totani's disk cut and his cell grid agree. A
    different cell size would need the boundary handled explicitly.
    """
    if abs(b_min / cell_deg - round(b_min / cell_deg)) > 1e-9:
        raise ValueError(
            f"|b| >= {b_min} does not fall on a {cell_deg} deg cell boundary; "
            "the disk cut would split cells and this returns a mask that "
            "silently keeps or drops whole ones")
    _, B = cell_centers(half_deg, cell_deg)
    return np.abs(B) >= b_min
