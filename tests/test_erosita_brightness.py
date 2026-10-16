"""Oscar's full-sky eROSITA brightness package (delivered 2026-09-21).

The package ships its own `verify_package.py` for integrity and normalisation.
These tests cover what that script cannot: that OUR transfer onto the Totani
grid is exact and correctly handed, and that the one number the analysis will
lean on -- the East/West brightness contrast, whose eastern half is a
reconstruction -- is stable across the coherent alternatives Oscar supplies.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.erosita import (  # noqa: E402
    BRIGHTNESS_EXTENSIONS, LOBE_MEAN_LON_DEG, UNIFORM_MEAN_LON_DEG,
    brightness_on_target, _BRIGHTNESS,
)
from totani.templates import grid_lb, pixel_solid_angle, target_header  # noqa: E402

_SYS = ROOT / "external" / "eROSITA_fullsky_lat_handoff" / "systematics" \
    / "erosita_brightness_uncertainties.fits"
pytestmark = pytest.mark.skipif(not _BRIGHTNESS.exists(), reason="brightness package not present")


def _roi(L, B):
    return (np.abs(L) <= 60) & (np.abs(B) >= 10) & (np.abs(B) <= 60)


def test_uniform_controls_reproduce_the_geometric_boundary_handedness():
    """The uniform controls are pure Zhang et al. geometry, so their mean
    longitudes must match the independently delivered 2026-09-08 boundary --
    north leaning negative, south positive. A longitude flip in either
    delivery would swap both signs and nothing else in the pipeline would."""
    t = target_header()
    L, B = grid_lb(t)
    omega = pixel_solid_angle(t)
    k = _roi(L, B)
    for lobe, key in (("north", "unif_whole_N"), ("south", "unif_whole_S")):
        w = brightness_on_target(t, key) * omega * k
        ml = float((w * L).sum() / w.sum())
        assert np.sign(ml) == np.sign(LOBE_MEAN_LON_DEG[lobe]), (
            f"{key} mean l = {ml:+.3f} but the geometric boundary gives "
            f"{LOBE_MEAN_LON_DEG[lobe]:+.3f} -- a longitude handedness flip")
        assert abs(ml - UNIFORM_MEAN_LON_DEG[lobe]) < 0.05, ml
        # and it agrees with the boundary map to better than a degree
        assert abs(ml - LOBE_MEAN_LON_DEG[lobe]) < 1.0, (ml, LOBE_MEAN_LON_DEG[lobe])


def test_transfer_is_exact_and_refuses_a_shifted_grid():
    """Same pixel size and half-pixel phase -> nearest neighbour moves nothing.
    A grid on a different phase must raise, not interpolate silently."""
    t = target_header()
    a = brightness_on_target(t, "whole_N")
    assert a.shape == (t["NAXIS1"], t["NAXIS2"]) and np.isfinite(a).all()
    with fits.open(_BRIGHTNESS, memmap=True) as h:
        src = np.asarray(h["T_N"].data, float)
    # the target is the central 960x960 of the 2880x1440 delivery
    assert np.array_equal(a, src[240:1200, 960:1920])

    bad = target_header()
    bad["CRPIX1"] = bad["CRPIX1"] + 0.3
    with pytest.raises(ValueError, match="share a phase"):
        brightness_on_target(bad, "whole_N")


def test_each_lobe_keeps_its_own_unit_integral_and_is_not_renormalised():
    """Oscar's lobes integrate to one over their own support; on our smaller
    grid they must integrate to LESS than one. Renormalising here would
    redefine the flux convention the README warns about."""
    t = target_header()
    omega = pixel_solid_angle(t)
    for key in ("whole_N", "whole_S", "outer_N", "outer_S"):
        s = float((brightness_on_target(t, key) * omega).sum())
        assert 0.4 < s < 1.0, f"{key} integrates to {s} on the Totani grid"


@pytest.mark.skipif(not _SYS.exists(), reason="systematics file not present")
def test_east_west_contrast_survives_every_coherent_alternative():
    """The eastern half of the map is reconstructed from an RGB image, so the
    East/West contrast is the number most at risk. Measured 2026-09-21: 0.911
    north, 0.819 south, moving by only a few per cent across calibration
    percentiles, smoothing, registration shifts, residual-shape transport and
    -- the real test -- ED1, which is calibrated from a DIFFERENT published
    display. If a redelivery breaks this, the morphology has changed."""
    from totani.healpix import pix2ang_ring

    d = fits.open(_SYS)["UNCERTAINTY"].data
    nside = 256
    l, b = pix2ang_ring(nside, np.arange(12 * nside * nside))
    l = np.where(l > 180, l - 360, l)
    roi = _roi(l, b)

    def ratio(col, valid):
        v = np.asarray(d[col], float)
        m = np.asarray(d[valid])
        out = {}
        for lab, sel in (("north", b > 0), ("south", b < 0)):
            k = m & sel & (v > 0) & roi
            out[lab] = float(v[k & (l > 0)].mean() / v[k & (l < 0)].mean())
        return out

    nom = ratio("FULL_NOMINAL", "VALID_FULL")
    assert nom["north"] == pytest.approx(0.911, abs=0.02)
    assert nom["south"] == pytest.approx(0.819, abs=0.02)
    for col in ("FULL_CAL_LOW", "FULL_CAL_HIGH", "FULL_SMOOTH075", "FULL_SMOOTH150",
                "FULL_REG_X", "FULL_REG_Y", "FULL_RESIDUAL_SHAPE"):
        r = ratio(col, "VALID_FULL")
        assert abs(r["north"] - nom["north"]) < 0.02, (col, r)
        assert abs(r["south"] - nom["south"]) < 0.02, (col, r)
    ed1 = ratio("FULL_ED1", "VALID_ED1")     # independent display and calibration
    assert abs(ed1["north"] - nom["north"]) < 0.06, ed1
    assert abs(ed1["south"] - nom["south"]) < 0.06, ed1
    # East is dimmer than West in both lobes, and more so in the south
    assert ed1["south"] < ed1["north"] < 1.0
