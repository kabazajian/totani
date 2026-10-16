"""The full background-template stack, on real reduction output."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

RUN = Path(__file__).resolve().parents[1] / "runs/fff2de5bc7e0"
have_run = pytest.mark.skipif(not (RUN / "expcube.fits").exists(),
                              reason="reduction not present")
have_catalog = pytest.mark.skipif(
    not Path("/Users/aba/data/fermi/catalog/gll_psc_v35.fit").exists(),
    reason="4FGL-DR4 not present")


def _edges():
    from reduce_totani import TOTANI_EDGES_GEV
    return TOTANI_EDGES_GEV[:10]


@have_run
@have_catalog
def test_stack_has_exactly_the_named_background_components():
    """15 free normalizations per bin -- no halo, no per-source freedom."""
    from totani.stack import build_background_stack

    st = build_background_stack(RUN, _edges())
    expected = {f"HI_ring{k}" for k in range(1, 5)} | {f"H2_ring{k}" for k in range(1, 5)}
    expected |= {"dust_pos", "dust_neg", "bubble", "loopI_shell1", "loopI_shell2",
                "isotropic", "ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6",
                "point_sources"}
    assert set(st["templates"]) == expected
    assert len(expected) == 19  # 8 gas + 2 dust + bubble + 2 loopI + isotropic + 4 ICS + PS
    assert "halo" not in st["templates"] and "dm" not in str(st["templates"]).lower()


@have_run
@have_catalog
def test_every_template_is_cell_shaped_and_finite():
    from totani.stack import build_background_stack

    st = build_background_stack(RUN, _edges())
    assert st["counts"].shape == (9, 12, 12)
    for name, arr in st["templates"].items():
        assert arr.shape == (9, 12, 12), name
        assert np.isfinite(arr).all(), name


@have_run
@have_catalog
def test_masking_reduces_but_does_not_wildly_distort_the_counts():
    """Measured 2026-09-03: masking removed 17,578 of the halo-box pixels and
    1.66% of bin-0 counts -- a small, physically reasonable correction, not a
    sign that the mask ate a large fraction of the ROI."""
    from totani.stack import build_background_stack

    st = build_background_stack(RUN, _edges())
    assert 0 < st["n_masked_pixels"] < 50000
    halo = np.abs(np.linspace(-55, 55, 12)) >= 10
    assert 1_500_000 < st["counts"][0][halo].sum() < 2_000_000


@have_run
@have_catalog
def test_dust_negative_stays_negative_after_the_full_pipeline():
    """A sign flip anywhere in the exposure/binning chain would be absorbed by
    a free f_l and never show up in a fit residual -- so it is asserted here
    directly, on the actual assembled template."""
    from totani.stack import build_background_stack

    st = build_background_stack(RUN, _edges())
    assert (st["templates"]["dust_neg"] <= 0).all()
    assert (st["templates"]["dust_pos"] >= 0).all()


@have_run
@have_catalog
def test_point_source_template_excludes_masked_sources():
    """🚨 REGRESSION FOR THE MASKED-SOURCE DEFECT (found by the 2026-09-04
    audit). The stack's point-source template must equal a direct call with
    sources inside the extended-source mask DROPPED, and must NOT equal the
    version that keeps them.

    The old version of this test asserted the opposite -- it compared against
    the unmasked direct call and passed, encoding the bug. 240 of the 2,651
    ROI sources sit in mask circles and carry 29-45% of the predicted flux;
    their photons are excluded from the data by the same mask, so depositing
    their model counts was a pure model/data inconsistency. Fixing it moved
    the fitted point-source normalization from 0.65-0.96 to 0.91-1.17.
    """
    from astropy.io import fits

    from totani.pointsources import extended_source_mask, point_source_cell_template
    from totani.stack import build_background_stack

    st = build_background_stack(RUN, _edges())
    with fits.open(RUN / "expcube.fits") as h:
        cube = h[0].data.astype(float)
        E = h["ENERGIES"].data["Energy"].astype(float)
        hdr = h[0].header.copy()
        for k in ("NAXIS3", "CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3"):
            hdr.remove(k, ignore_missing=True)
        hdr["NAXIS"] = 2
    with fits.open("/Users/aba/data/fermi/catalog/gll_psc_v35.fit") as h:
        cat = h["LAT_Point_Source_Catalog"].data
        ext = h["ExtendedSources"].data
    l = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"])
    b = cat["GLAT"]
    in_roi = (np.abs(l) <= 60) & (np.abs(b) <= 60)
    in_mask = extended_source_mask(l, b, ext)
    keep = in_roi & ~in_mask

    assert st["n_sources_dropped_in_mask"] == int((in_roi & in_mask).sum()) > 0
    assert st["n_sources_kept"] == int(keep.sum())

    masked_out = point_source_cell_template(cat[keep], l[keep], b[keep],
                                            _edges(), cube, E, hdr)
    assert np.array_equal(st["templates"]["point_sources"], masked_out)

    kept_all = point_source_cell_template(cat[in_roi], l[in_roi], b[in_roi],
                                          _edges(), cube, E, hdr)
    assert not np.array_equal(st["templates"]["point_sources"], kept_all), (
        "the stack is still including masked sources -- the defect is back")
    # and the difference is LARGE, not a rounding detail
    assert kept_all.sum() / st["templates"]["point_sources"].sum() > 1.2
