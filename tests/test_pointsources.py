"""4FGL-DR4 point sources on 10-degree cells.

🚨 `test_spectra_reproduce_the_catalogs_own_flux1000` is the load-bearing one.
The three 4FGL parameterizations -- especially the PLEC4 SuperExpCutoff -- are
easy to get subtly wrong in a way that yields entirely plausible fluxes, and
this project's rule 1 exists for exactly that. The check uses the catalog's OWN
model-derived integrals, so it is self-contained and needs no external number.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from totani.pointsources import (  # noqa: E402
    MEV_PER_ERG,
    bin_photon_flux,
    boundary_fraction,
    cell_edges,
    cell_index,
    dnde,
    extended_source_mask,
    misplaced_flux_fraction,
)

CATALOG = Path("/Users/aba/data/fermi/catalog/gll_psc_v35.fit")
have_catalog = pytest.mark.skipif(not CATALOG.exists(),
                                  reason="4FGL-DR4 not present")


def _catalog():
    from astropy.io import fits
    return fits.open(CATALOG)


def _halo_rows(c):
    l = np.where(c["GLON"] > 180, c["GLON"] - 360, c["GLON"])
    sel = (np.abs(l) <= 60) & (np.abs(c["GLAT"]) >= 10) & (np.abs(c["GLAT"]) <= 60)
    return c[sel], l[sel], c["GLAT"][sel]


@have_catalog
@pytest.mark.parametrize("kind", ["PowerLaw", "LogParabola", "PLSuperExpCutoff"])
def test_spectra_reproduce_the_catalogs_own_flux1000(kind):
    """Flux1000 is the 1-100 GeV photon flux OF THE FITTED MODEL, so integrating
    our dN/dE over the same range must return it. Measured 2026-09-02: median
    ratio within 1e-4 for all three forms over 400 sources each."""
    c = _catalog()["LAT_Point_Source_Catalog"].data
    rows = c[np.char.strip(c["SpectrumType"]) == kind][:400]
    assert len(rows) > 50, f"too few {kind} sources to test"
    E = np.logspace(3, 5, 256)                       # 1-100 GeV in MeV
    ratio = np.array([np.trapezoid(dnde(r, E), E) / r["Flux1000"] for r in rows])
    assert np.median(ratio) == pytest.approx(1.0, abs=1e-3), (
        f"{kind}: median {np.median(ratio):.5f}")


@have_catalog
def test_energy_flux100_also_reproduces():
    """A second, independent catalog integral -- different range and weighting
    (E dN/dE over 100 MeV-100 GeV), so it fails differently if a form is wrong."""
    c = _catalog()["LAT_Point_Source_Catalog"].data
    rows = c[:300]
    E = np.logspace(2, 5, 384)
    got = np.array([np.trapezoid(dnde(r, E) * E, E) / MEV_PER_ERG
                    / r["Energy_Flux100"] for r in rows])
    assert np.median(got) == pytest.approx(1.0, abs=1e-3)


def test_unknown_spectrum_type_is_refused():
    """A silently substituted power law is the plausible wrong number."""
    row = {"SpectrumType": "Exponential", "Pivot_Energy": 1000.0}
    with pytest.raises(ValueError, match="unknown SpectrumType"):
        dnde(row, np.array([1000.0]))


def test_cell_grid_divides_and_indexes():
    e = cell_edges()
    assert len(e) == 13 and e[0] == -60 and e[-1] == 60
    with pytest.raises(ValueError, match="does not divide"):
        cell_edges(half_deg=60.0, cell_deg=7.0)
    r, c = cell_index(np.array([-59.9, 0.0, 59.9, 61.0]),
                      np.array([-59.9, 0.0, 59.9, 0.0]))
    assert (r[:3] >= 0).all() and (c[:3] >= 0).all()
    assert r[3] == -1 and c[3] == -1, "a source outside the ROI must be flagged"


@have_catalog
def test_every_halo_source_lands_in_one_of_the_120_cells():
    c = _catalog()["LAT_Point_Source_Catalog"].data
    rows, l, b = _halo_rows(c)
    assert len(rows) == 1425
    r, cc = cell_index(l, b)
    assert (r >= 0).all() and (cc >= 0).all()
    # 12 columns of longitude x 10 rows of latitude at |b| >= 10.
    assert len(set(zip(r.tolist(), cc.tolist()))) == 120


@have_catalog
def test_bin_flux_shape_and_positivity():
    from reduce_totani import TOTANI_EDGES_GEV

    c = _catalog()["LAT_Point_Source_Catalog"].data
    rows, _, _ = _halo_rows(c)
    f = bin_photon_flux(rows[:50], TOTANI_EDGES_GEV)
    assert f.shape == (50, 13)
    assert (f > 0).all(), "a non-positive photon flux is unphysical"
    # Falling spectra: the lowest bin must dominate the highest for essentially
    # every source, and certainly in the median.
    assert np.median(f[:, 0] / f[:, -1]) > 1e3


@have_catalog
def test_misplaced_flux_is_far_smaller_than_the_boundary_count():
    """The distinction this module insists on: 27% of sources sit within a PSF
    width of a cell edge at 1 GeV, but only ~8% of their FLUX crosses, because
    a source 0.7 deg away contributes only its far tail. Quoting the source
    count as the systematic would overstate it by 3-4x."""
    c = _catalog()["LAT_Point_Source_Catalog"].data
    rows, l, b = _halo_rows(c)
    bf = boundary_fraction(l, b, psf_deg=0.727)
    mf = misplaced_flux_fraction(l, b, psf68_deg=0.727)
    assert 0.2 < bf < 0.4
    assert mf < bf / 2.0
    # and it must fall steeply with energy -- the PSF narrows by 7x
    mf_high = misplaced_flux_fraction(l, b, psf68_deg=0.10)
    assert mf_high < mf / 5.0


@have_catalog
def test_extended_source_mask_covers_known_extents():
    """Totani sec 2.3 masks 2x the semi-major axis of each catalogued extended
    source. The mask must be non-empty over the ROI and stay a small fraction."""
    h = _catalog()
    ext = h["ExtendedSources"].data
    l = np.linspace(-60, 60, 241)
    b = np.linspace(-60, 60, 241)
    L, B = np.meshgrid(l, b)
    m = extended_source_mask(L, B, ext)
    assert m.any(), "no extended source masked anywhere in the ROI"
    assert m.mean() < 0.25, f"masking {m.mean():.1%} of the ROI is too much"
