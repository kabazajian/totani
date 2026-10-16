"""Template construction on the 960x960 Totani grid.

🚨 `test_delivered_header_coordinates_are_finite` is the important one. It is the
regression for a bug that produced NO error: `astropy.wcs` returned NaN for every
latitude in A2020's delivered gas cube, so an ROI mask selected zero pixels and a
flux check reported "source flux = 0". Meanwhile `regrid_car` -- which reads the
same headers linearly and never calls astropy -- was perfectly correct. A
coordinate helper and a regridder that disagree about what a header means is the
quietest possible failure, and this file exists to keep them agreeing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import totani  # noqa: E402,F401
from totani.templates import (  # noqa: E402
    CONFIGS,
    bubble_on_target,
    delivered_paths,
    gas_dust_on_target,
    grid_lb,
    ics_four_ring,
    loopi_on_target,
    target_header,
    _load,
)

DELIVERY = totani.GCE_FISHER_ROOT / "data/external/Templates_arxiv2003.10416"
have_delivery = pytest.mark.skipif(
    not DELIVERY.is_dir(), reason="A2020 delivery not present")


def test_target_grid_geometry():
    t = target_header()
    assert (t["NAXIS1"], t["NAXIS2"]) == (960, 960)
    assert t["CDELT1"] < 0 < t["CDELT2"]        # longitude to the LEFT
    L, B = grid_lb(t)
    assert L[0, 0] > L[0, -1], "longitude must decrease with column index"
    assert np.isclose(L.min(), -59.9375) and np.isclose(L.max(), 59.9375)
    assert np.isclose(B.min(), -59.9375) and np.isclose(B.max(), 59.9375)


@have_delivery
def test_delivered_header_coordinates_are_finite():
    """THE REGRESSION. A2020's gas cube carries CRVAL2 = -89.75 with CRPIX2 = 1
    and spans 0-360 in longitude; with the -CAR code that `load_map` supplies,
    astropy reads it as a rotated projection with out-of-domain longitudes and
    returns NaN. `grid_lb` must reproduce the header's actual plate-carree grid.
    """
    _, hdr, _ = _load(delivered_paths(totani.GCE_FISHER_ROOT)["HI"], fix_h2=True)
    L, B = grid_lb(hdr)
    assert np.isfinite(L).all() and np.isfinite(B).all()
    assert np.isclose(B.min(), -89.75) and np.isclose(B.max(), 89.75)
    assert np.isclose(L.min(), -179.75) and np.isclose(L.max(), 179.75)


@have_delivery
def test_regrid_conserves_flux_into_the_roi():
    """Measured 2026-09-02: ratio 0.999997 for all four HI annuli. The 3e-6 is
    the edge of the check's own mask (source pixel CENTRES inside the ROI vs
    target pixels extending to the boundary), not a bias -- which is why it is
    identical across rings rather than varying with their morphology."""
    from gcereduce.regrid import regrid_car

    data, hdr, _ = _load(delivered_paths(totani.GCE_FISHER_ROOT)["HI"], fix_h2=True)
    t = target_header()
    L, B = grid_lb(hdr)
    Lt, Bt = grid_lb(t)
    mask = (np.abs(L) <= 59.9375) & (np.abs(B) <= 59.9375)
    for ring in range(4):
        src = (data[ring] * np.cos(np.deg2rad(B)))[mask].sum() * 0.5 * 0.5
        tgt = (regrid_car(data[ring], hdr, t)
               * np.cos(np.deg2rad(Bt))).sum() * 0.125 * 0.125
        assert tgt / src == pytest.approx(1.0, abs=1e-4), f"ring {ring + 1}"


@have_delivery
def test_gas_dust_stack_is_ten_templates():
    g = gas_dust_on_target(target_header())
    assert sorted(g) == sorted(
        [f"HI_ring{k}" for k in range(1, 5)]
        + [f"H2_ring{k}" for k in range(1, 5)] + ["dust_pos", "dust_neg"])
    assert all(v.shape == (960, 960) for v in g.values())
    # The negative-residual map is genuinely negative; a sign flip would be
    # absorbed by a free normalization and never noticed.
    assert g["dust_neg"].sum() < 0 < g["dust_pos"].sum()


def test_ts_has_no_default_and_config_is_checked():
    """Oscar said to pick one fixed T_S but not which; a silent default is
    exactly the plausible-looking wrong choice rule 1 exists to prevent."""
    with pytest.raises(ValueError, match="150K"):
        delivered_paths(totani.GCE_FISHER_ROOT, ts="160K")
    with pytest.raises(ValueError, match="unknown config"):
        delivered_paths(totani.GCE_FISHER_ROOT, config="nope")


def test_baseline_variant_is_pohl_mag5():
    """Oscar, 2026-08-11. Baseline = hydrodynamic gas + the mag5 dust cut."""
    assert CONFIGS["baseline"] == ("pohl", "mag5")
    p = delivered_paths(totani.GCE_FISHER_ROOT)
    assert p["HI"].name == "HI_pohl_4annuli_TS150K.fits"
    assert p["dust_pos"].name == "posres_pohl_mag5_150K.fits"


def test_ics_four_ring_grouping():
    """(1+2+3), 4, 5, 6 -- CONFIRMED by Oscar, was previously our inference."""
    six = [np.full((2, 2), float(k)) for k in range(1, 7)]
    four = ics_four_ring(six)
    assert len(four) == 4
    assert four[0][0, 0] == 1 + 2 + 3
    assert [f[0, 0] for f in four[1:]] == [4.0, 5.0, 6.0]
    with pytest.raises(ValueError, match="expected 6"):
        ics_four_ring(six[:5])


@have_delivery
def test_bubble_is_a_binary_mask_of_the_right_shape():
    """Sampled from the MOLLWEIDE copy: the HEALPix one has no ORDERING keyword,
    so RING vs NESTED is unstated and the wrong choice gives scattered noise a
    free normalization would absorb."""
    b = bubble_on_target(target_header())
    assert set(np.unique(b)) <= {0.0, 1.0}
    assert b.mean() == pytest.approx(0.192, abs=0.01)
    L, B = grid_lb(target_header())
    # Two lobes: present well off the plane in BOTH hemispheres, centred on l=0.
    assert (b[B > 15] > 0).any() and (b[B < -15] > 0).any()
    on = b > 0
    assert abs(np.average(L[on])) < 10.0, "bubble is not centred near l = 0"


def test_loopi_on_target_returns_two_shells():
    lo = loopi_on_target(target_header())
    assert sorted(lo) == ["loopI_shell1", "loopI_shell2"]
    assert all(v.shape == (960, 960) and (v >= 0).all() for v in lo.values())
