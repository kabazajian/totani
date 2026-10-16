"""Oscar Macias' structured bubble map: the loader and its Tier-1 properties.

Fast: builds the target grid from `target_header()` and regrids one 720x360
map. No run directory needed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.templates import (  # noqa: E402
    bubble_macias_on_target,
    bubble_on_target,
    grid_lb,
    pixel_solid_angle,
    target_header,
)

MAP = ROOT / "external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits"
pytestmark = pytest.mark.skipif(not MAP.exists(), reason="Macias map not present")


@pytest.fixture(scope="module")
def maps():
    t = target_header()
    return t, bubble_macias_on_target(t), bubble_on_target(t)


def test_shape_finite_nonnegative(maps):
    t, mac, _ = maps
    assert mac.shape == (int(t["NAXIS2"]), int(t["NAXIS1"]))
    assert np.isfinite(mac).all()
    assert (mac >= 0).all()
    assert mac.max() > 0


def test_it_is_structured_not_a_mask(maps):
    """The whole point: this must NOT be loadable as a 0/1 mask."""
    _, mac, _ = maps
    assert len(np.unique(mac)) > 100


def test_handedness_agrees_with_the_ackermann_mask(maps):
    """Pixel correlation must be higher as-is than with longitude flipped.
    This project has met four handedness traps; assert, do not assume."""
    _, mac, ack = maps

    def corr(a, b):
        a = a - a.mean(); b = b - b.mean()
        return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))

    assert corr(mac, ack) > corr(mac[:, ::-1], ack) + 0.02


def test_flat_footprint_is_inside_the_structured_map(maps):
    """94% of the Ackermann footprint is non-zero in Oscar's map (Tier 1)."""
    _, mac, ack = maps
    assert np.mean(mac[ack > 0] > 0) > 0.9


def test_flux_conservation_through_the_regrid(maps):
    """Flux over |b| >= 10 in the ROI must match the native map's over the
    same box to a few percent -- `regrid_car` conserves flux, and this is the
    check that the WCS was read the way the map was written."""
    from astropy.io import fits
    from astropy.wcs import WCS

    t, mac, _ = maps
    L, B = grid_lb(t)
    om = pixel_solid_angle(t)
    sel = (np.abs(B) >= 10) & (np.abs(B) <= 60) & (np.abs(L) <= 60)
    ours = float((mac * om)[sel].sum())

    with fits.open(MAP) as h:
        d = np.asarray(h[0].data, float); hdr = h[0].header
    w = WCS(hdr)
    yy, xx = np.mgrid[:d.shape[0], :d.shape[1]]
    ll, bb = w.wcs_pix2world(xx, yy, 0)
    ll = np.where(ll > 180, ll - 360, ll)
    nat_om = np.radians(0.5) ** 2 * np.cos(np.radians(bb))
    nsel = (np.abs(bb) >= 10) & (np.abs(bb) <= 60) & (np.abs(ll) <= 60)
    native = float((d * nat_om)[nsel].sum())
    assert ours == pytest.approx(native, rel=0.03)


def test_refuses_a_bad_path(tmp_path):
    with pytest.raises(FileNotFoundError):
        bubble_macias_on_target(target_header(), path=tmp_path / "nope.fits")
