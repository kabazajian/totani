"""The GALPROP (Ackermann+ 2012 SLZ6R30T150C2) loader: the bin integration is
exact on a power law, and the delivered maps are what the docstring says."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.galprop import (  # noqa: E402
    _ROOT, GAS_PHASES, RING_BOUNDS_KPC, TOTANI_RING_GROUPS, TOTANI_SET, integrate_bins, load_skymap)

have_maps = pytest.mark.skipif(
    not (_ROOT / f"pi0_decay_healpix_{TOTANI_SET}.gz").exists(),
    reason="GALPROP archive not downloaded")


@pytest.mark.parametrize("index", [-2.0, -2.7, -1.5, -1.0])
def test_bin_integral_is_exact_on_a_power_law(index):
    """Nodes on the archive's own grid (factor 1.2997 from 50 MeV). The
    log-log interpolant of a power law is the power law, and the per-segment
    integral is analytic, so this must be exact to rounding."""
    E = 50.0 * 1.2996566 ** np.arange(38)
    A = 3e-7
    spec = (A * (E / 1000.0) ** index)[None, :] * np.array([[1.0], [2.0]])
    edges_gev = np.array([1.162, 1.962, 3.315, 5.598, 9.456, 15.971])
    got = integrate_bins(spec, E, edges_gev)
    e = edges_gev * 1e3
    if index == -1.0:
        exact = A * 1000.0 * np.log(e[1:] / e[:-1])
    else:
        exact = A * 1000.0 / (index + 1) * ((e[1:] / 1000.0) ** (index + 1) - (e[:-1] / 1000.0) ** (index + 1))
    assert np.allclose(got[0], exact, rtol=1e-9)
    assert np.allclose(got[1], 2 * exact, rtol=1e-9)


def test_bin_outside_tabulated_range_raises():
    E = 50.0 * 1.2996566 ** np.arange(38)
    with pytest.raises(ValueError, match="outside"):
        integrate_bins(np.ones((1, 38)), E, [0.01, 0.1])


@have_maps
def test_maps_are_dnde_with_the_pion_bump():
    from totani.healpix import ang2pix_ring

    spec, E, nside = load_skymap("pi0_decay")
    assert nside == 128 and E.size == 38 and np.isclose(E[0], 50.0)
    gc = ang2pix_ring(nside, 0.0, 0.0)
    e2 = spec[gc] * E ** 2
    peak = E[np.argmax(e2)]
    assert 300 < peak < 1500, f"E^2 dN/dE peaks at {peak} MeV; not a pi0 spectrum in dN/dE units"
    hi = ang2pix_ring(nside, 0.0, 45.0)
    assert spec[gc, 13] > 30 * spec[hi, 13]


@have_maps
def test_ics_components_sum_to_the_total():
    tot, E, _ = load_skymap("ics_isotropic")
    parts = sum(load_skymap(f"ics_isotropic_comp_{i}")[0] for i in (1, 2, 3))
    assert np.allclose(parts, tot, rtol=1e-6, atol=0)


@have_maps
def test_templates_on_target_are_positive_and_plane_bright():
    from totani.galprop import galprop_on_target
    from totani.templates import grid_lb, target_header

    t = target_header(npix=96, pixel_deg=1.25)          # coarse grid, fast
    edges = [1.162, 1.962, 3.315]
    d = galprop_on_target(t, edges)
    L, B = grid_lb(t)
    for name, cube in d.items():
        assert cube.shape == (2,) + L.shape and (cube > 0).all(), name
    gas, ics = d["gas"], d["ics"]
    assert gas[0][np.abs(B) < 2].mean() > 10 * gas[0][np.abs(B) > 40].mean()
    # ICS is the smoother component: a smaller plane/high-latitude contrast
    assert (ics[0][np.abs(B) < 2].mean() / ics[0][np.abs(B) > 40].mean()
            < gas[0][np.abs(B) < 2].mean() / gas[0][np.abs(B) > 40].mean())


def test_ring_groups_tile_the_galaxy_exactly():
    """Totani's five ranges are unions of whole GALPROP annuli, no gaps."""
    used = sorted(k for rings in TOTANI_RING_GROUPS.values() for k in rings)
    assert used == list(range(1, 18))
    for tag, rings in TOTANI_RING_GROUPS.items():
        lo = RING_BOUNDS_KPC[rings[0] - 1][0]; hi = RING_BOUNDS_KPC[rings[-1] - 1][1]
        assert tag == f"R{lo:g}-{hi:g}", (tag, lo, hi)
        for a, b in zip(rings[:-1], rings[1:]):
            assert RING_BOUNDS_KPC[a - 1][1] == RING_BOUNDS_KPC[b - 1][0]


@have_maps
def test_bremss_total_includes_ionized_and_rings_sum_to_it():
    """The finding that fixed the gas template: `bremss` is the whole thing."""
    br, _, _ = load_skymap("bremss")
    bri, _, _ = load_skymap("bremss_ionized")
    rings = sum(load_skymap(f"bremss_{p}_ring_{k}", allow_negative=True)[0]
                for p in GAS_PHASES for k in range(1, 18))
    hii = sum(load_skymap(f"bremss_HII_ring_{k}", allow_negative=True)[0] for k in range(1, 18))
    assert np.allclose(rings, br, rtol=1e-6, atol=1e-30)
    assert np.allclose(hii.sum(axis=0), bri.sum(axis=0), rtol=3e-3)
    assert 0.03 < bri.sum() / br.sum() < 0.08


@have_maps
def test_ring_groups_sum_to_the_gas_total_on_target():
    from totani.galprop import galprop_on_target
    from totani.templates import target_header

    t = target_header(npix=48, pixel_deg=2.5)
    edges = [1.162, 1.962]
    whole = galprop_on_target(t, edges)["gas"]
    parts = galprop_on_target(t, edges, gas_rings=True)
    assert sorted(parts) == sorted([f"gas_{k}" for k in TOTANI_RING_GROUPS] + ["ics"])
    rel = sum(v for k, v in parts.items() if k.startswith("gas_")) / whole - 1
    # clipping of residual negatives + nonlinear (log-log) bin integration
    # done per group instead of once: a ~1e-3 effect, not zero, by design
    assert np.abs(rel).max() < 1e-2, np.abs(rel).max()
    assert np.abs(np.median(rel)) < 2e-3
    for tag, frac in galprop_on_target.last_clipped.items():
        assert 0 <= frac < 5e-3, (tag, frac)


@have_maps
def test_ring_files_have_small_negatives_that_totals_do_not():
    """The reason `allow_negative` exists; pinned so a re-delivery that changes
    it is noticed."""
    with pytest.raises(ValueError, match="negative"):
        load_skymap("pi0_decay_H2R_ring_17")
    spec, _, _ = load_skymap("pi0_decay_H2R_ring_17", allow_negative=True)
    assert (spec < 0).any() and -spec[spec < 0].sum() / spec[spec > 0].sum() < 1e-2
    load_skymap("pi0_decay")           # strict, passes
