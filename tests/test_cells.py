"""Pixel -> 10-degree cell binning."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from totani.cells import bin_to_cells, cell_centers, halo_cell_mask  # noqa: E402
from totani.templates import grid_lb, target_header  # noqa: E402

CCUBE = Path(__file__).resolve().parents[1] / "runs/fff2de5bc7e0/ccube.fits"


def _target_grid():
    t = target_header()
    return grid_lb(t)


def test_binning_conserves_the_total():
    L, B = _target_grid()
    rng = np.random.default_rng(0)
    d = rng.random((3,) + L.shape)
    c = bin_to_cells(d, L, B)
    assert c.shape == (3, 12, 12)
    assert c.sum() == pytest.approx(d.sum(), rel=1e-12)


def test_a_2d_map_stays_2d():
    L, B = _target_grid()
    c = bin_to_cells(np.ones(L.shape), L, B)
    assert c.shape == (12, 12)
    assert c.sum() == pytest.approx(L.size)


def test_mask_removes_flux_rather_than_moving_it():
    L, B = _target_grid()
    d = np.ones(L.shape)
    mask = np.abs(L) < 5.0
    c = bin_to_cells(d, L, B, mask=mask)
    assert c.sum() == pytest.approx(L.size - mask.sum())


def test_shape_mismatch_is_refused():
    L, B = _target_grid()
    with pytest.raises(ValueError, match="does not match grid"):
        bin_to_cells(np.ones((5, 5)), L, B)


def test_the_disk_cut_falls_on_a_cell_boundary():
    """With 10 deg cells and |b| >= 10 the split is exact -- 120 of 144 cells,
    latitude centres +/-15 ... +/-55. If it did not, whole cells would be kept
    or dropped without that being visible anywhere."""
    m = halo_cell_mask()
    assert m.sum() == 120
    _, B = cell_centers()
    assert sorted(set(np.abs(B[m]).tolist())) == [15.0, 25.0, 35.0, 45.0, 55.0]
    with pytest.raises(ValueError, match="cell boundary"):
        halo_cell_mask(b_min=12.0)


@pytest.mark.skipif(not CCUBE.exists(), reason="reduction not present")
def test_cells_reproduce_the_counts_cube():
    """Measured 2026-09-02 on the real reduction: 11,803,231 counts conserved
    exactly, and 3,736,576 in the halo cells -- identical to the pixel-space
    answer, which is the point."""
    from astropy.io import fits

    h = fits.open(CCUBE)
    d = h[0].data.astype(float)
    hd = h[0].header.copy()
    for k in ("NAXIS3", "CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3"):
        hd.remove(k, ignore_missing=True)
    hd["NAXIS"] = 2
    L, B = grid_lb(hd)
    c = bin_to_cells(d, L, B)
    assert int(c.sum()) == int(d.sum()) == 11803231
    halo_pixels = d[:, (np.abs(B) >= 10) & (np.abs(B) <= 60)].sum()
    assert int(c[:, halo_cell_mask()].sum()) == int(halo_pixels) == 3736576
