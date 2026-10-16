"""Stage A of the ranking: registries, the tarball reader, and that the
cached cell templates are the same numbers the stack builder makes."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from totani import ranking  # noqa: E402


def test_registry_is_the_full_factorial():
    assert len(ranking.MODEL_CODES) == 128 == len(set(ranking.MODEL_CODES))
    assert ranking.code_to_tag("SL_Z6_R30_T150_C2") == "54_Lorimer_z6kpc_R30kpc_Ts150K_EBV2mag"
    assert ranking.code_to_tag("SS_Z10_R20_T100000_C5") == "54_SNR_z10kpc_R20kpc_Ts100000K_EBV5mag"
    with pytest.raises(ValueError):
        ranking.code_to_tag("SL_Z7_R30_T150_C2")


have_baseline = pytest.mark.skipif(
    not (ranking.LOCAL_TARS / "SL_Z6_R30_T150_C2.tar").exists(), reason="baseline tarball not present")


@have_baseline
def test_tar_reader_matches_the_extracted_files():
    from totani.galprop import load_skymap
    src = ranking.TarSkymaps(ranking.LOCAL_TARS / "SL_Z6_R30_T150_C2.tar", ranking.code_to_tag("SL_Z6_R30_T150_C2"))
    try:
        a, Ea, na = src.load("ics_isotropic")
    finally:
        src.close()
    b, Eb, nb = load_skymap("ics_isotropic")
    assert na == nb and np.array_equal(Ea, Eb) and np.array_equal(a, b)


@have_baseline
def test_model_cells_equal_the_stack_builders_templates(tmp_path, monkeypatch):
    """The cache must be the same arithmetic as `build_background_stack(
    diffuse='galprop_totani')`, or Stage B ranks something else."""
    from reduce_totani import TOTANI_EDGES_GEV
    from totani.stack import build_background_stack
    monkeypatch.setattr(ranking, "CACHE", tmp_path)
    edges = TOTANI_EDGES_GEV[:10]
    cells = ranking.model_cells("SL_Z6_R30_T150_C2", "gas1_ics1", edges)
    st = build_background_stack(ranking.RUN, edges, bubble="macias", diffuse="galprop_totani")
    for n in ("galprop_gas", "galprop_ics"):
        assert np.allclose(cells[n], st["templates"][n], rtol=1e-10, atol=0)
    assert set(cells) == {"galprop_gas", "galprop_ics"}


@have_baseline
def test_assemble_reproduces_the_totani_diffuse_baseline(tmp_path, monkeypatch):
    """Assembled stack with the reference structures == the stack the
    totani_diffuse run used (structured bubbles there; flat here -> use flat)."""
    from reduce_totani import TOTANI_EDGES_GEV
    from totani.stack import build_background_stack
    monkeypatch.setattr(ranking, "CACHE", tmp_path)
    edges = TOTANI_EDGES_GEV[:10]
    st = ranking.assemble("SL_Z6_R30_T150_C2", "gas1_ics1", ("bubble_flat", "loopI_shell1", "loopI_shell2"), edges)
    ref = build_background_stack(ranking.RUN, edges, bubble="ackermann", diffuse="galprop_totani")
    assert np.array_equal(st["counts"], ref["counts"])
    for n, arr in st["templates"].items():
        rn = "bubble" if n == "bubble_flat" else n
        assert np.allclose(arr, ref["templates"][rn], rtol=1e-10, atol=0), n
    assert set(st["templates"]) == {"bubble_flat" if n == "bubble" else n for n in ref["templates"]}


def test_radio_loops_are_where_the_paper_puts_them():
    """Loop IV sits in the northern l < 0 ROI; Loop II only clips the eastern
    edge; Loop III is outside. Areas from the projected thin-shell template."""
    from totani.structures import RADIO_LOOPS, radio_loop
    from totani.templates import grid_lb, target_header

    assert RADIO_LOOPS["loopIV"] == (315.0, 48.5, 19.8)
    t = target_header(npix=192, pixel_deg=0.625)
    L, B = grid_lb(t)
    iv = radio_loop("loopIV")(t, None)
    assert (iv > 0).any() and (iv >= 0).all() and np.isfinite(iv).all()
    on = iv > 0
    assert L[on].mean() < -30 and B[on].mean() > 35            # north, negative l
    assert 0.02 < on.mean() < 0.15
    ii = radio_loop("loopII")(t, None)
    assert 0 < (ii > 0).mean() < 0.15 and L[ii > 0].mean() > 40   # eastern edge only
    assert (radio_loop("loopIII")(t, None) > 0).mean() < 0.01


_GIEM = ROOT / "external" / "fssc" / "gll_iem_v07.fits"
_SUN = Path("/Users/aba/work/gce-fisher/data/external/Pohl_etal_2022/Sun_and_Moon_Template_IC.fits")


@pytest.mark.skipif(not _GIEM.exists(), reason="gll_iem_v07 not downloaded")
def test_gll_iem_v07_samples_onto_the_target_grid():
    """Half-pixel latitude offset handled by averaging; plane bright; positive."""
    from totani.structures import gll_iem_v07
    from totani.templates import grid_lb, target_header

    t = target_header()
    cube = gll_iem_v07(t, [1.162, 1.962, 3.315])
    L, B = grid_lb(t)
    assert cube.shape == (2,) + L.shape and np.isfinite(cube).all() and (cube > 0).all()
    assert cube[0][np.abs(B) < 1].mean() > 10 * cube[0][np.abs(B) > 45].mean()
    # symmetric handling of the half-pixel offset: rows just above and below
    # b = 0 are nearly equal (the model is close to plane-symmetric there)
    i0 = np.argmin(np.abs(B[:, 0] + 0.0625)); i1 = np.argmin(np.abs(B[:, 0] - 0.0625))
    assert np.abs(cube[0][i1].mean() / cube[0][i0].mean() - 1) < 0.05


@pytest.mark.skipif(not _SUN.exists(), reason="Sun/Moon IC template not present")
def test_sun_moon_ic_is_an_ecliptic_band():
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from totani.structures import sun_moon_ic
    from totani.templates import grid_lb, target_header

    t = target_header(npix=192, pixel_deg=0.625)
    cube = sun_moon_ic(t, [1.162, 1.962])
    L, B = grid_lb(t)
    assert cube.shape == (1,) + L.shape and np.isfinite(cube).all() and (cube >= 0).all()
    ecl = SkyCoord(l=L.ravel() * u.deg, b=B.ravel() * u.deg, frame="galactic").barycentrictrueecliptic
    beta = np.abs(ecl.lat.deg).reshape(L.shape)
    near, far = cube[0][beta < 5].mean(), cube[0][beta > 40].mean()
    assert near > 3 * far, (near, far)


@have_baseline
def test_greedy_selection_takes_the_bubbles_first():
    """On Totani's baseline model with two candidates, the bubbles must be
    chosen before Loop I and both must clear the BIC threshold."""
    from reduce_totani import TOTANI_EDGES_GEV
    from totani import selection
    edges = TOTANI_EDGES_GEV[:10]
    sh = ranking.shared_cells(edges)
    path = selection.greedy_for_model("SL_Z6_R30_T150_C2", "gas1_ics1", edges, sh,
                                      candidates=["loopI_shells", "bubble_flat"])
    assert path["chosen"] == ["bubble_flat", "loopI_shells"]
    assert path["steps"][0]["delta"] > selection.bic_threshold(9) * 10
    assert selection.bic_threshold(9) == pytest.approx(9 * np.log(1080))


def test_sign_free_candidate_keeps_the_better_sign_per_bin():
    """A synthetic stack: a component whose sign flips at bin 3."""
    from totani import selection
    rng = np.random.default_rng(1)
    n_bins, n = 9, 12
    base = np.full((n_bins, n, n), 200.0)
    shape = np.zeros((n_bins, n, n)); shape[:, 2:10, 2:10] = 30.0
    sign = np.array([1, 1, 1, -1, -1, -1, 1, 1, 1])
    counts = base + sign[:, None, None] * shape
    counts = rng.poisson(np.clip(counts, 1, None)).astype(float)
    sh = {"isotropic": base, "point_sources": np.zeros_like(base), "counts": counts, "x": shape}
    stack = {"counts": counts, "n_bins": n_bins, "templates": {"isotropic": base}}
    selection.CANDIDATES["_test"] = ("x",)
    selection.SIGN_FREE.add("x")
    try:
        _, v, signs = selection.add_candidate(stack, "_test", "m", "v", None, sh)
    finally:
        selection.CANDIDATES.pop("_test"); selection.SIGN_FREE.discard("x")
    assert signs == sign.tolist()
