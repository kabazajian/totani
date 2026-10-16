"""The halo signal template -- NFW rho^2, Totani's baseline (sec 2.3).

🚨 `test_j_map_wants_1d_vectors_not_a_meshgrid` is the load-bearing test. The
first version of `nfw_rho2_on_target` passed the FULL 2-D `grid_lb()` output
to `j_map`, which itself builds a 2-D grid internally from two 1-D vectors --
so a 2-D input silently picked up ANOTHER broadcast axis inside `_los_nodes`,
turning a (960,960) grid into a (960,960,960)-shaped intermediate. It ran for
77s and reached 13.85 GB RSS (still climbing) before being killed, on a machine
shared with another session's 30+ GB queue -- caught, not repeated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from totani.halo import TOTANI_NFW_PARAMS, nfw_rho2_on_target  # noqa: E402
from totani.templates import axis_coords, grid_lb, target_header  # noqa: E402


def test_totani_parameters_are_the_via_lactea_ii_values():
    """Sec 2.3, quoted in the module docstring. gamma=1 is not a separately
    stated number -- it IS what "NFW" (vs a generalized profile) means, per
    the quoted density formula rho = rho_s x^-1 (1+x)^-2."""
    assert TOTANI_NFW_PARAMS == dict(profile="gnfw", gamma=1.0, rs=21.0,
                                     rho_sun=0.42, R_sun=8.0)


def test_j_map_wants_1d_vectors_not_a_meshgrid():
    """Direct regression, against the SAME convention `gce-fisher`'s own tests
    use (`test_dm_profile.py`'s `GRID = np.linspace(...)`, a 1-D vector)."""
    from gcereduce.dm_templates import j_map

    l1d = np.linspace(-4.0, 4.0, 9)
    b1d = np.linspace(-4.0, 4.0, 7)
    out = j_map(l1d, b1d, n_s=50)
    assert out.shape == (7, 9), (
        "j_map(l, b) must return (len(b), len(l)) from 1-D inputs -- if this "
        "fails, the calling convention changed upstream and "
        "nfw_rho2_on_target's use of axis_coords() must be re-checked")

    # And the FAILURE MODE this test exists to catch: feeding it an
    # ALREADY-2-D grid gives a WRONG shape rather than raising -- which is
    # exactly how the bug went undetected until memory, not a shape error,
    # caught it.
    L2d, B2d = np.meshgrid(l1d, b1d)
    wrong = j_map(L2d, B2d, n_s=50)
    assert wrong.shape != (7, 9), (
        "if j_map now raises or returns the RIGHT shape for a 2-D input, "
        "nfw_rho2_on_target no longer needs axis_coords() -- update both")


def test_small_grid_shape_and_peak_at_the_galactic_centre():
    t = target_header(npix=40, pixel_deg=1.0)
    J = nfw_rho2_on_target(t, n_s=200, row_chunk=10)
    assert J.shape == (40, 40)
    assert np.isfinite(J).all()
    assert (J > 0).all()          # rho^2 integral, never zero or negative
    L, B = grid_lb(t)
    iy, ix = np.unravel_index(np.argmax(J), J.shape)
    assert abs(L[iy, ix]) <= 1.0 and abs(B[iy, ix]) <= 1.0, (
        "the brightest pixel must be within one pixel of the GC")


def test_matches_grid_lb_meshgrid_via_axis_coords():
    """`axis_coords`'s 1-D vectors, meshed, must reproduce `grid_lb`'s 2-D
    output exactly -- the guarantee that switching functions changes nothing
    else about the grid."""
    t = target_header(npix=40, pixel_deg=1.0)
    L, B = grid_lb(t)
    l_c, b_c = axis_coords(t)
    Lm, Bm = np.meshgrid(l_c, b_c)
    assert np.array_equal(L, Lm) and np.array_equal(B, Bm)


def test_row_chunking_does_not_change_the_result():
    t = target_header(npix=40, pixel_deg=1.0)
    whole = nfw_rho2_on_target(t, n_s=200, row_chunk=40)
    chunked = nfw_rho2_on_target(t, n_s=200, row_chunk=7)   # doesn't divide evenly
    assert np.array_equal(whole, chunked)


@pytest.mark.slow
def test_full_grid_peaks_at_the_galactic_centre():
    """The real 960x960 analysis grid. Measured 2026-09-04: 17.3 s, 3.73 GB
    peak RSS, peak pixel within 0.1 deg of (0,0)."""
    t = target_header()
    J = nfw_rho2_on_target(t)
    assert J.shape == (960, 960)
    assert np.isfinite(J).all() and (J > 0).all()
    L, B = grid_lb(t)
    iy, ix = np.unravel_index(np.argmax(J), J.shape)
    assert abs(L[iy, ix]) < 0.5 and abs(B[iy, ix]) < 0.5
