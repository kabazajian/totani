"""PPPC4DMID reading and the joint DM fit."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from totani import dmspectra  # noqa: E402

have_pppc = pytest.mark.skipif(not dmspectra.PPPC_PATH.exists(), reason="PPPC table not present")


@have_pppc
def test_table_shape_and_channels():
    t = dmspectra.load_pppc()
    assert t["m"].size == 62 and t["log10x"].size == 179
    for c in ("b", "tau", "W", "mu"):
        assert c in t and t[c].shape == (62, 179) and (t[c] >= 0).all()


@have_pppc
def test_the_one_conversion_that_matters():
    """int dN/dE dE over E must equal int dN/dlog10x dlog10x over the same x
    range -- the E ln10 factor is the whole test."""
    t = dmspectra.load_pppc(); m = 100.0
    i = int(np.argmin(np.abs(t["m"] - m)))
    lx = t["log10x"]; row = t["b"][i]
    sel = lx >= -3                                        # x from 1e-3 to 1
    n_logx = np.trapezoid(row[sel], lx[sel])
    E = np.logspace(-3, 0, 4000) * t["m"][i]
    n_E = np.trapezoid(dmspectra.dnde("b", float(t["m"][i]), E), E)
    assert abs(n_E / n_logx - 1) < 1e-2
    assert dmspectra.dnde("b", 100.0, np.array([150.0]))[0] == 0.0        # nothing above m


@have_pppc
def test_bin_yields_are_positive_and_shift_with_mass():
    edges = [1.162, 1.962, 3.315, 5.598, 9.456, 15.971, 26.976, 45.564, 76.959, 129.985]
    n50 = dmspectra.bin_yields("b", 50.0, edges); n500 = dmspectra.bin_yields("b", 500.0, edges)
    assert (n50 >= 0).all() and (n500 >= 0).all()
    assert n50[-1] == 0.0                                  # no photons above 50 GeV
    ratio = n500[:-1] / n50[:-1]
    assert (np.diff(ratio) > 0).all()                     # the 500 GeV spectrum is harder: ratio rises with E
    dE = dmspectra.bin_widths_mev(edges)
    assert np.allclose(dmspectra.amplitudes(1e-26, 50.0, "b", edges), 1e-26 / (8 * np.pi * 2500.0) * n50 / dE)


_CELLS = ROOT / "results" / "fits" / "ranking_cells" / "shared__base.npz"


@pytest.mark.skipif(not (_CELLS.exists() and dmspectra.PPPC_PATH.exists()), reason="caches not built")
def test_joint_fit_recovers_an_injected_signal():
    """Inject a b-bbar 100 GeV halo at a known <sigma v> into the counts of
    the A2020 stack (cell space, cached) and recover it."""
    from reduce_totani import TOTANI_EDGES_GEV
    from totani import ranking
    edges = TOTANI_EDGES_GEV[:10]
    st = ranking.assemble("a2020", "-", ("bubble_flat", "loopI_shell1", "loopI_shell2"), edges)
    J = ranking.shared_cells(edges)["J"]
    sv_true = 3e-25
    a = dmspectra.amplitudes(sv_true, 100.0, "b", edges)
    st2 = dict(st); st2["counts"] = st["counts"] + a[:, None, None] * J
    r0 = dmspectra.joint_fit(st, J, "b", 100.0, edges, n_scan=21)
    r1 = dmspectra.joint_fit(st2, J, "b", 100.0, edges, n_scan=21)
    assert abs((r1["sigmav_best"] - r0["sigmav_best"]) / sv_true - 1) < 0.1
    assert r1["sigmav_lo68"] < r1["sigmav_best"] < r1["sigmav_hi68"]
    assert r1["ts_vs_zero"] > r0["ts_vs_zero"]


@pytest.mark.skipif(not (_CELLS.exists() and dmspectra.PPPC_PATH.exists()), reason="caches not built")
def test_scale_matches_totani():
    """Absolute scale: Totani reports <sigma v> ~ 0.5-0.9e-24 cm^3/s for
    0.5-0.8 TeV DM from a ~6,000-count 20 GeV halo (NFW-rho^2). Fitting his
    configuration (his model, his residual templates, shells) with b b-bar at
    600 GeV must land within a factor of a few of 1e-24 -- a factor 1e6 off
    is what the missing 1/dE gave."""
    from reduce_totani import TOTANI_EDGES_GEV
    from totani import ranking, selection
    edges = TOTANI_EDGES_GEV[:10]
    sh = ranking.shared_cells(edges)
    st = selection.base_stack("SL_Z6_R30_T150_C2", "gas1_ics1", edges, sh)
    for c in ("bubble_totani", "loopI_shells"):
        st, _, _ = selection.add_candidate(st, c, "SL_Z6_R30_T150_C2", "gas1_ics1", edges, sh)
    r = dmspectra.joint_fit(st, sh["J"], "b", 600.0, edges, n_scan=21)
    assert 1e-25 < r["sigmav_best"] < 1e-23, r["sigmav_best"]
    assert r["ts_vs_zero"] > 25


@have_pppc
def test_no_yield_below_threshold_is_reported_not_fitted():
    edges = [1.162, 1.962, 3.315, 5.598, 9.456, 15.971, 26.976, 45.564, 76.959, 129.985]
    assert dmspectra.bin_yields("W", 50.0, edges).sum() == 0.0
