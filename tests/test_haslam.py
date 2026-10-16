"""The hand-written HEALPix RING pixelisation, and the Haslam map through it.

No healpy in this environment, so `haslam.ang2pix_ring` / `pix2ang_ring` are
transcriptions of the reference algorithm. The round trip over EVERY pixel of
an NSIDE-512 map is the test that they are transcribed correctly: a wrong
branch, a wrong shift or a wrong ring boundary breaks it for thousands of
pixels at once.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.haslam import (  # noqa: E402
    _HASLAM,
    ang2pix_ring,
    haslam_nps_on_target,
    nps_disc,
    haslam_on_target,
    load_haslam,
    pix2ang_ring,
)


@pytest.mark.parametrize("nside", [1, 2, 4, 16, 512])
def test_ring_round_trip_every_pixel(nside):
    npix = 12 * nside * nside
    pix = np.arange(npix)
    l, b = pix2ang_ring(nside, pix)
    back = ang2pix_ring(nside, l, b)
    bad = np.flatnonzero(back != pix)
    assert bad.size == 0, f"nside {nside}: {bad.size} pixels fail, first {bad[:5]}"


def test_ring_structure_nside1():
    """12 pixels: 4 in each cap ring and 4 on the equator, at the textbook
    positions (Gorski+ 2005 fig. 4)."""
    l, b = pix2ang_ring(1, np.arange(12))
    assert np.allclose(np.sort(np.unique(np.round(b, 6))),
                       np.round([-41.810315, 0.0, 41.810315], 6))
    assert np.allclose(np.sort(l[:4]), [-135, -45, 45, 135])
    assert np.allclose(np.sort(l[4:8]), [-180, -90, 0, 90])


def test_pixels_partition_the_sphere():
    """Every direction maps to a valid pixel, and a fine random sample hits
    every pixel of a small map about equally (equal-area property)."""
    rng = np.random.default_rng(0)
    n = 400_000
    z = rng.uniform(-1, 1, n); l = rng.uniform(-180, 180, n)
    b = np.degrees(np.arcsin(z))
    nside = 8
    pix = ang2pix_ring(nside, l, b)
    assert pix.min() >= 0 and pix.max() < 12 * nside * nside
    counts = np.bincount(pix, minlength=12 * nside * nside)
    assert counts.min() > 0
    assert counts.std() / counts.mean() < 0.06


def test_poles_and_seam():
    """The poles land in the first / last ring; l = 180 and l = -180 agree."""
    nside = 64
    assert ang2pix_ring(nside, 0.0, 89.9999) < 4
    assert ang2pix_ring(nside, 0.0, -89.9999) >= 12 * nside * nside - 4
    assert ang2pix_ring(nside, 180.0, 12.3) == ang2pix_ring(nside, -180.0, 12.3)


have_map = pytest.mark.skipif(not _HASLAM.exists(), reason="Haslam map not downloaded")


@have_map
def test_map_loads_as_expected():
    data, nside = load_haslam()
    assert nside == 512 and data.size == 12 * 512 * 512
    assert 10 < np.median(data) < 40                # K, high-latitude sky
    assert data.max() > 1000                        # the Galactic centre


@have_map
def test_sampling_sees_the_galaxy_and_the_north_polar_spur():
    from totani.templates import grid_lb, target_header

    t = target_header()
    m = haslam_on_target(t)
    L, B = grid_lb(t)
    assert m.shape == L.shape and np.isfinite(m).all() and (m > 0).all()
    plane = m[np.abs(B) < 2].mean()
    high = m[np.abs(B) > 50].mean()
    assert plane > 5 * high
    # The North Polar Spur: brighter than the mirror region in the south at
    # the same |b| and l. Classic NPS position ~ (l, b) = (+25, +30).
    nps = m[(np.abs(L - 25) < 5) & (np.abs(B - 30) < 5)].mean()
    mirror = m[(np.abs(L - 25) < 5) & (np.abs(B + 30) < 5)].mean()
    assert nps > 1.3 * mirror


def test_shell1_footprint_is_the_whole_sky_so_shell2_disc_is_used():
    """The reason `nps_disc` exists: shell 1 alone covers every pixel."""
    from totani.loopi import loopi_maps
    from totani.templates import grid_lb, target_header

    L, B = grid_lb(target_header())
    sh = loopi_maps(L, B)
    assert (sh["loopI_shell1"] > 0).all()
    disc = nps_disc(L, B)
    assert 0.3 < disc.mean() < 0.9
    # the disc is centred in the north at negative l: more north than south
    assert disc[B > 0].mean() > 1.5 * disc[B < 0].mean()
    assert disc[L < 0].mean() > disc[L > 0].mean()


@have_map
def test_restricted_form_is_inside_the_disc():
    from totani.templates import grid_lb, target_header

    t = target_header()
    r = haslam_nps_on_target(t)
    L, B = grid_lb(t)
    d = nps_disc(L, B)
    assert (r[~d] == 0).all() and (r[d] > 0).all()
