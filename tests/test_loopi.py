"""Loop I is the one template we GENERATE rather than receive, so it is the one
that can be silently wrong. These are the checks that would catch it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from totani.loopi import SHELLS, loopi_maps, shell_map, sphere_path_length  # noqa: E402

POHL_LOOPI = Path("/Users/aba/work/gce-fisher/data/external/Pohl_etal_2022/"
                  "LoopI_Normalized.fits")


def test_parameters_are_the_lat_team_values():
    """Ackermann et al. 2014 ApJ 793,64 sec 3.1, read from the PDF (rule 1a).
    Totani takes these rather than Wolleben's originals, and so do we."""
    s1, s2 = SHELLS
    assert (s1["l_deg"], s1["b_deg"], s1["d_pc"], s1["r_in_pc"], s1["r_out_pc"]) \
        == (341.0, 3.0, 78.0, 62.0, 81.0)
    assert (s2["l_deg"], s2["b_deg"], s2["d_pc"], s2["r_in_pc"], s2["r_out_pc"]) \
        == (332.0, 37.0, 95.0, 58.0, 82.0)


@pytest.mark.parametrize("cos_psi,d,R,expect", [
    ( 1.0, 95.0, 82.0, 2 * 82.0),   # origin OUTSIDE, straight at the centre
    (-1.0, 95.0, 82.0, 0.0),        # origin outside, pointing away
    ( 1.0, 78.0, 81.0, 78.0 + 81),  # origin INSIDE, straight at the centre
    (-1.0, 78.0, 81.0, 81.0 - 78),  # origin inside, pointing away
])
def test_path_length_analytic_limits(cos_psi, d, R, expect):
    """🚨 The two shells put the Sun on OPPOSITE sides of their outer spheres
    (d=78 < r_out=81 for shell 1; d=95 > r_out=82 for shell 2). A formula that
    forgets to clip the near root at s=0 gets the inside case wrong by half a
    chord and still returns a plausible map."""
    got = float(sphere_path_length(np.array(cos_psi), d, R))
    assert got == pytest.approx(expect, abs=1e-9)


def test_the_sun_is_inside_shell_one_and_outside_shell_two():
    """Guards the premise of the test above -- if a future parameter edit moves
    the Sun, the clipping test stops testing what it claims to."""
    s1, s2 = SHELLS
    assert s1["r_in_pc"] < s1["d_pc"] < s1["r_out_pc"]
    assert s2["d_pc"] > s2["r_out_pc"]


def test_shell_is_empty_where_the_sphere_is_missed():
    """Far from shell 2's centre the sight line misses it entirely."""
    s2 = SHELLS[1]
    anti = shell_map(np.array([s2["l_deg"] - 180.0]), np.array([-s2["b_deg"]]),
                     s2["l_deg"], s2["b_deg"], s2["d_pc"],
                     s2["r_in_pc"], s2["r_out_pc"])
    assert anti.item() == 0.0


def test_inner_radius_must_be_smaller():
    with pytest.raises(ValueError, match="inner radius"):
        shell_map(np.array([0.0]), np.array([0.0]), 0, 0, 78, 81, 62)


@pytest.mark.skipif(not POHL_LOOPI.exists(), reason="Pohl archive not present")
def test_matches_pohls_delivered_loopi_map():
    """🚨 THE FALSIFIABLE CHECK. Pohl et al. (2022) ship a Loop I map over
    40x40 deg built from this same geometric model. If our generated template is
    the same object, the two must agree over the shared footprint; if the
    parameters or the path-length model were wrong, they would not.

    Measured 2026-09-02: r = 0.9997 for the unweighted sum of the two shells,
    and a free two-shell fit recovers near-equal weights (0.01894, 0.01926),
    consistent with Pohl's map being their sum.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    h = fits.open(POHL_LOOPI)[0]
    ref = h.data.astype(float)
    ny, nx = ref.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    # ⚠ Always through the WCS -- this project has met four longitude-handedness
    # traps and indexing by hand is how each of them started.
    L, B = WCS(h.header).wcs_pix2world(xx, yy, 0)

    m = loopi_maps(L, B)
    total = m["loopI_shell1"] + m["loopI_shell2"]
    r = np.corrcoef(total.ravel(), ref.ravel())[0, 1]
    assert r > 0.999, f"generated Loop I disagrees with Pohl's map: r = {r:.4f}"

    A = np.column_stack([m["loopI_shell1"].ravel(),
                         m["loopI_shell2"].ravel(),
                         np.ones(ref.size)])
    c, *_ = np.linalg.lstsq(A, ref.ravel(), rcond=None)
    assert c[0] > 0 and c[1] > 0, f"a shell fits NEGATIVE: {c}"
    assert abs(c[0] - c[1]) / max(c[0], c[1]) < 0.05, f"weights not near-equal: {c}"
