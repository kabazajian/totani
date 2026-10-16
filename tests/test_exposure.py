"""Exposure-weighting of spatial templates into a per-bin cell basis."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from totani.exposure import flat_spectrum_bin_exposure, single_energy_bin_exposure  # noqa: E402

RUN = Path(__file__).resolve().parents[1] / "runs/fff2de5bc7e0"
have_run = pytest.mark.skipif(not (RUN / "expcube.fits").exists(),
                              reason="reduction not present")


def _fake_exposure(n_planes=5, ny=4, nx=4):
    E = np.logspace(3, 5, n_planes)             # 1-100 GeV in MeV
    # a smoothly INCREASING exposure with energy, uniform across the tiny grid
    cube = np.stack([np.full((ny, nx), 1e10 * (e / E[0])) for e in E])
    return cube, E


def test_flat_spectrum_integral_is_monotonic_and_positive():
    cube, E = _fake_exposure()
    edges_gev = np.array([1.0, 3.0, 10.0, 30.0, 100.0])
    basis = flat_spectrum_bin_exposure(cube, E, edges_gev)
    assert basis.shape == (4, 4, 4)
    assert (basis > 0).all()
    # exposure rises with energy, and each successive bin is both higher-E
    # AND wider in log-space here, so the integral must rise bin to bin.
    totals = basis.mean(axis=(1, 2))
    assert (np.diff(totals) > 0).all()


def test_single_energy_matches_flat_at_a_constant_exposure():
    """If exposure is CONSTANT in energy, both conventions must agree exactly
    (up to the trivial factor of bin width), since there is no spectral shape
    left to disagree about."""
    ny, nx, n_planes = 3, 3, 6
    E = np.logspace(3, 5, n_planes)
    cube = np.ones((n_planes, ny, nx)) * 7.5e11
    edges_gev = np.array([1.0, 3.0, 10.0, 30.0])
    flat = flat_spectrum_bin_exposure(cube, E, edges_gev)
    single = single_energy_bin_exposure(cube, E, edges_gev)
    edges_mev = edges_gev * 1000.0
    widths = np.diff(edges_mev)
    assert np.allclose(flat, single * widths[:, None, None], rtol=1e-3)


def test_out_of_range_energy_is_refused():
    cube, E = _fake_exposure()
    with pytest.raises(ValueError, match="outside the tabulated"):
        flat_spectrum_bin_exposure(cube, E, np.array([0.5, 2.0]))
    # Center = sqrt(500*1000) = 707 GeV, genuinely outside [1,100] GeV --
    # unlike [1, 1000], whose center (31.6 GeV) is INSIDE range and must NOT
    # raise, since single_energy_bin_exposure only checks centers.
    with pytest.raises(ValueError, match="outside the tabulated"):
        single_energy_bin_exposure(cube, E, np.array([500.0, 1000.0]))
    single_energy_bin_exposure(cube, E, np.array([1.0, 1000.0]))  # must NOT raise


@have_run
def test_real_exposure_gives_a_sensible_basis():
    """Measured 2026-09-03 on the actual Totani reduction: exposure varies only
    ~6% from bin 0 to bin 8 over 15 years of near-uniform sky coverage, so the
    flat-spectrum basis should track the bin's own energy WIDTH far more than
    any exposure variation."""
    import sys as _s
    _s.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from astropy.io import fits
    from reduce_totani import TOTANI_EDGES_GEV

    with fits.open(RUN / "expcube.fits") as h:
        cube = h[0].data.astype(float)
        E = h["ENERGIES"].data["Energy"].astype(float)
    edges = TOTANI_EDGES_GEV[:10]
    basis = flat_spectrum_bin_exposure(cube, E, edges)
    assert basis.shape == (9,) + cube.shape[1:]
    assert (basis > 0).all()
