"""The cutout grid reconstruction for the Galp21 band forecast.

`scripts/band_forecast.py` rebuilds the (l, b) of `galp21.cutout`'s output from
the documented WCS rather than reading a header, because the cutout is a
concatenation of two column slices and no single WCS describes it. This project
has met four longitude-handedness traps, so the reconstruction is asserted
against the real header whenever the map share is reachable, and against the
documented convention always.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, "/Users/aba/work/gce-fisher/src")

from band_forecast import FOUR_RING, HALF_B, HALF_L, native_grid  # noqa: E402

MAPS = Path("/private/nfs/Data2/gce-fisher/GCE_Galprop_maps_2026/maps")


def test_longitude_is_ascending_and_spans_the_cut():
    L, _ = native_grid()
    l = L[0]
    assert np.all(np.diff(l) > 0), "longitude must ascend across the concatenation"
    assert l.min() == pytest.approx(-HALF_L + 0.05)
    assert l.max() == pytest.approx(HALF_L - 0.05)
    assert l.size == int(2 * HALF_L * 10)


def test_latitude_is_ascending_and_spans_the_cut():
    _, B = native_grid()
    b = B[:, 0]
    assert np.all(np.diff(b) > 0)
    assert b.min() == pytest.approx(-HALF_B + 0.05)
    assert b.max() == pytest.approx(HALF_B - 0.05)


def test_the_galactic_centre_is_not_at_an_edge():
    """The GC sits on the native array's wrap; after the cut it must be central."""
    L, _ = native_grid()
    l = L[0]
    i = int(np.argmin(np.abs(l)))
    assert abs(i - l.size // 2) <= 1, "the seam was not rejoined about l = 0"


def test_four_ring_collapse_uses_every_group_once():
    flat = [g for grp in FOUR_RING for g in grp]
    assert sorted(flat) == [1, 2, 3, 4, 5, 6]
    assert FOUR_RING[0] == (1, 2, 3), "A2020 collapses rings 1-3, confirmed by Oscar"


@pytest.mark.skipif(not MAPS.exists(), reason="Galp21 map share not mounted")
def test_reconstruction_matches_the_delivered_wcs():
    """Rebuild against astropy's own inversion of the native header."""
    from astropy.io import fits
    from astropy.wcs import WCS

    import gcereduce.galp21 as g21
    g21.BASE = MAPS
    p = g21.component_path("SA0_F98", "ics_anisotropic", ring=1, isrf=0)
    with fits.open(p, memmap=True) as h:
        w = WCS(h[0].header).celestial

    lo, hi = g21.cut_columns(HALF_L)
    rows = g21.cut_rows(HALF_B)
    cols = np.concatenate([np.arange(lo.start, lo.stop),
                           np.arange(hi.start, hi.stop)])
    r = np.arange(rows.start, rows.stop)
    l_wcs, b_wcs = w.wcs_pix2world(cols, np.full(cols.size, r[0]), 0)
    l_wcs = np.where(l_wcs > 180, l_wcs - 360, l_wcs)

    L, B = native_grid()
    assert np.allclose(L[0], l_wcs, atol=1e-6)
    _, b_col = w.wcs_pix2world(np.zeros(r.size), r, 0)
    assert np.allclose(B[:, 0], b_col, atol=1e-6)


def test_file_enumeration_matches_the_delivery_structure():
    """918 native files: per model 51 IC (17 rings x 3 ISRF) and 17 per gas
    phase for each of pion decay and bremsstrahlung (3 phases each)."""
    import gcereduce.galp21 as g21
    from band_forecast import KINDS, _files, _phases

    per_model = 0
    for kind in KINDS:
        for phase in _phases(kind):
            files = _files("SA0_F98", kind, phase)
            assert sorted(files) == [1, 2, 3, 4, 5, 6]
            n = sum(len(v) for v in files.values())
            assert n == (51 if kind == "ics_anisotropic" else 17), (kind, phase, n)
            per_model += n
            # every native ring appears in exactly one group
            rings = [r for grp in g21.SIX_RING_GROUPS.values() for r in grp]
            assert sorted(rings) == list(range(1, 18))
    assert per_model == 153
    assert per_model * len(g21.MODELS) == 918


def test_gas_kinds_have_three_phases_and_ic_has_none():
    from band_forecast import _phases
    assert _phases("ics_anisotropic") == (None,)
    assert _phases("pion_decay") == ("H2R", "HIR", "HII")
    assert _phases("bremss") == ("H2R", "HIR", "HII")
