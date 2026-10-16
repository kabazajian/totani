"""The per-bin fitter reaches the convex minimum: the linear-space polish
(likelihood._polish_linear) closes the plateau failures found 2026-09-22.

Both cases are real-reduction fits that the theta-only fitter got wrong:
adding two free templates RAISED -2lnL (impossible at a convex minimum), and
the residual-construction pair on a five-ring row sat tens of units above
the restart-verified minimum in single bins.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import likelihood as LK, ranking, selection  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
have_run = pytest.mark.skipif(not (ranking.RUN / "ccube.fits").exists(), reason="reduction not present")


def _bins(stack):
    keep = halo_cell_mask()
    sr = LK.restrict_stack(stack, keep)
    names = sorted(sr["templates"])
    for k in range(sr["n_bins"]):
        yield k, np.stack([sr["templates"][n][k] for n in names]), sr["counts"][k]


@have_run
def test_adding_free_templates_never_raises_the_minimum():
    """Five-ring SL_Z10 + the A2020 dust pair: bin 0 was 207 ABOVE the bare fit."""
    sh = ranking.shared_cells(EDGES)
    base = selection.base_stack("SL_Z10_R30_T150_C5", "gas5_ics3", EDGES, sh)
    plus = selection.with_templates(base, {"dust_pos": sh["a2020_dust_pos"], "dust_neg": sh["a2020_dust_neg"]})
    a = selection.bg_neg2lnl(base)
    b = selection.bg_neg2lnl(plus)
    assert np.all(b <= a + 1e-6), (a - b)


@have_run
def test_residual_pair_bins_reach_the_restart_verified_minimum():
    """SS_Z10 five-ring + residual pair (+ sign): bins 1 and 3 were 42 and 22 high.
    A random-restart search must not beat the fitter by more than the
    tolerance slack of the flattest five-ring directions."""
    sh = ranking.shared_cells(EDGES)
    m, v = "SS_Z10_R30_T150_C5", "gas5_ics3"
    res = selection.totani_residual_cells(m, v, EDGES)
    st = selection.with_templates(selection.base_stack(m, v, EDGES, sh),
                                  {"bubble_totani_pos": res["bubble_totani_pos"],
                                   "bubble_totani_neg": res["bubble_totani_neg"]})
    rng = np.random.default_rng(0)
    for k, T, n in _bins(st):
        if k not in (0, 1, 2, 3):
            continue
        val, amp, r = LK.fit_bin(T, n)
        best = val
        for _ in range(5):
            x0 = np.clip(r.x + rng.uniform(-3, 3, r.x.size), *LK.THETA_BOUNDS)
            best = min(best, LK.fit_bin(T, n, x0=x0)[0])
        assert val - best < 3.0, (k, val - best)


@have_run
def test_polish_is_a_no_op_on_a_clean_fit():
    """A one-gas row: the theta fit was already at the minimum to 1e-2."""
    sh = ranking.shared_cells(EDGES)
    st = selection.base_stack("SL_Z6_R30_T150_C2", "gas1_ics1", EDGES, sh)
    for k, T, n in _bins(st):
        val, amp, r = LK.fit_bin(T, n)
        v_theta = LK.neg2lnl(r.x, T / np.abs(T).sum(axis=(1, 2))[:, None, None], n)[0]
        assert val == pytest.approx(v_theta, abs=1e-6)
