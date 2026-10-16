"""The pieces of Totani's bubble-template construction that can be tested
without a stack: the NFW-rho^2.5 GCE profile and the masked sky smoothing."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from totani.bubble_totani import gaussian_smooth_sky, nfw_rho_power_on_target  # noqa: E402
from totani.halo import nfw_rho2_on_target  # noqa: E402
from totani.templates import grid_lb, target_header  # noqa: E402


def test_rho25_is_the_rho2_profile_made_steeper():
    t = target_header(npix=48, pixel_deg=2.5)
    L, B = grid_lb(t)
    r2 = nfw_rho2_on_target(t, n_s=300)
    r25 = nfw_rho_power_on_target(t, 2.5, n_s=300)
    assert (r25 > 0).all() and np.isfinite(r25).all()
    theta = np.degrees(np.arccos(np.cos(np.radians(B)) * np.cos(np.radians(L))))
    ratio = r25 / r2
    # more centrally concentrated: the ratio falls monotonically with angle
    inner = ratio[theta < 10].mean(); mid = ratio[(theta > 20) & (theta < 30)].mean(); outer = ratio[theta > 50].mean()
    assert inner > mid > outer
    # and the power = 2 case reproduces halo.nfw_rho2_on_target
    assert np.allclose(nfw_rho_power_on_target(t, 2.0, n_s=300), r2, rtol=1e-12)


def test_masked_smoothing_keeps_a_constant_and_ignores_the_mask():
    t = target_header(npix=96, pixel_deg=1.25)
    L, B = grid_lb(t)
    img = np.full(L.shape, 3.0)
    w = np.ones(L.shape); w[40:44, 50:54] = 0                      # a masked hole
    img[40:44, 50:54] = 1e6                                          # garbage under the mask
    s = gaussian_smooth_sky(img, w, B, sigma_deg=2.0, pixel_deg=1.25)
    assert np.allclose(s[w > 0], 3.0, atol=1e-9)


def test_smoothing_width_is_one_degree_on_the_sky():
    """A delta at high latitude spreads wider in l (by 1/cos b) than in b."""
    t = target_header(npix=192, pixel_deg=0.625)
    L, B = grid_lb(t)
    img = np.zeros(L.shape); w = np.ones(L.shape)
    i, j = np.unravel_index(np.argmin(np.abs(B - 50) + np.abs(L - 0)), L.shape)
    img[i, j] = 1.0
    s = gaussian_smooth_sky(img, w, B, sigma_deg=1.0, pixel_deg=0.625)
    prof_b = s[:, j]; prof_l = s[i, :]
    sig_b = np.sqrt(np.sum(prof_b * (B[:, j] - 50) ** 2) / prof_b.sum())
    sig_l = np.sqrt(np.sum(prof_l * (L[i, :]) ** 2) / prof_l.sum())
    assert abs(sig_b - 1.0) < 0.15
    assert abs(sig_l - 1.0 / np.cos(np.radians(50))) < 0.25


_RUN = ROOT / "runs" / "fff2de5bc7e0"


@pytest.mark.skipif(not (_RUN / "ccube.fits").exists(), reason="reduction not present")
def test_cached_common_parts_change_nothing():
    """The in-process cache of the model-independent half must be a pure
    speedup: the same templates, the same counts, the same mask."""
    import time
    from totani import bubble_totani as bt
    from reduce_totani import TOTANI_EDGES_GEV

    edges = TOTANI_EDGES_GEV[:10]
    bt._COMMON.clear()
    t0 = time.time(); a = bt.pixel_model_for_bin(_RUN, edges, 2, diffuse="a2020"); cold = time.time() - t0
    t0 = time.time(); b = bt.pixel_model_for_bin(_RUN, edges, 2, diffuse="a2020"); warm = time.time() - t0
    assert set(a["templates"]) == set(b["templates"])
    for n in a["templates"]:
        assert np.array_equal(a["templates"][n], b["templates"][n]), n
    for k in ("counts", "flat_sr", "mask", "L", "B"):
        assert np.array_equal(a[k], b[k]), k
    assert a["n_sources"] == b["n_sources"]
    assert warm < 0.5 * cold, f"cache gave no speedup: {cold:.1f}s then {warm:.1f}s"
