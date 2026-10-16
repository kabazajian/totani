"""Poisson cell likelihood -- background components only, no halo template.

Three kinds of check: the analytic gradient against finite differences (a wrong
gradient would still "converge", just to the wrong answer, and nothing else
here would catch it); recovery of known amplitudes from synthetic Poisson data;
and, where the real reduction exists, that the assembled background stack fits
without pathology.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from totani.likelihood import (  # noqa: E402
    THETA_BOUNDS,
    expected_counts,
    fit_all_bins,
    fit_bin,
    neg2lnl,
)

RUN = Path(__file__).resolve().parents[1] / "runs/fff2de5bc7e0"
have_run = pytest.mark.skipif(not (RUN / "expcube.fits").exists(),
                              reason="reduction not present")
have_catalog = pytest.mark.skipif(
    not Path("/Users/aba/data/fermi/catalog/gll_psc_v35.fit").exists(),
    reason="4FGL-DR4 not present")


def _synthetic(seed=0, n_comp=5, n_cell=12):
    rng = np.random.default_rng(seed)
    T = rng.uniform(0.1, 10.0, size=(n_comp, n_cell, n_cell))
    true_amp = rng.uniform(1.0, 500.0, size=n_comp)
    n_obs = rng.poisson(expected_counts(np.log10(true_amp), T)).astype(float)
    return T, n_obs, true_amp


def test_analytic_gradient_matches_finite_differences():
    """🚨 THE LOAD-BEARING TEST. `neg2lnl`'s gradient formula was reused from
    `a2020_stage2_fit.py`, but its tensordot axes were generalised for this
    module's tensor shapes -- exactly the kind of adaptation that can silently
    contract the wrong axes and still produce a plausible-looking converged
    fit at the WRONG optimum. Checked against central finite differences."""
    T, n_obs, _ = _synthetic()
    theta = np.array([0.3, 1.7, -0.5, 2.1, 0.9])
    _, grad = neg2lnl(theta, T, n_obs)
    h = 1e-6
    numeric = np.empty_like(theta)
    for i in range(len(theta)):
        tp, tm = theta.copy(), theta.copy()
        tp[i] += h
        tm[i] -= h
        vp, _ = neg2lnl(tp, T, n_obs)
        vm, _ = neg2lnl(tm, T, n_obs)
        numeric[i] = (vp - vm) / (2 * h)
    assert np.allclose(grad, numeric, rtol=1e-4, atol=1e-3), (grad, numeric)


def test_recovers_known_amplitudes_and_conserves_total_flux():
    """Total predicted counts must match total observed almost exactly (the
    global normalization is always well constrained); individual amplitudes
    are recovered well when the templates are not near-degenerate."""
    T, n_obs, true_amp = _synthetic(seed=1)
    val, amp, r = fit_bin(T, n_obs)
    assert r.status == 0, "did not converge"
    mu = expected_counts(np.log10(amp), T)
    assert mu.sum() == pytest.approx(n_obs.sum(), rel=1e-3)
    # median ratio across components should be close to 1, even if any single
    # weakly-identified component is not
    assert np.median(amp / true_amp) == pytest.approx(1.0, abs=0.5)


def test_amplitudes_are_never_negative():
    """The whole point of the 10**theta parameterisation: no bounds needed on
    the coefficient itself, and it cannot be violated regardless of what the
    optimizer does with theta."""
    T, n_obs, _ = _synthetic(seed=2)
    _, amp, _ = fit_bin(T, n_obs)
    assert (amp > 0).all()
    assert (amp >= 10.0 ** THETA_BOUNDS[0]).all()
    assert (amp <= 10.0 ** THETA_BOUNDS[1]).all()


def test_fit_is_deterministic_given_the_same_start():
    T, n_obs, _ = _synthetic(seed=3)
    _, amp1, _ = fit_bin(T, n_obs)
    _, amp2, _ = fit_bin(T, n_obs)
    assert np.array_equal(amp1, amp2)


@have_run
@have_catalog
def test_background_only_fit_on_the_real_reduction():
    """🚨 NO HALO TEMPLATE HERE -- deliberately. This checks the PIPELINE
    (convergence, positive amplitudes, sane total-flux agreement, no absurd
    per-cell blowup); it is NOT a goodness-of-fit gate and makes no claim about
    the 20 GeV excess. Measured 2026-09-04: total flux agreement is excellent
    (relative ~1e-7) while per-cell residuals are NOT small (halo mean ~7
    sigma, max ~22 sigma; disk max ~39 sigma) -- recorded in
    reference_values.yaml as an open finding, not asserted as either a pass or
    a discovery here. What IS asserted: the halo residual stays in the
    'notable but not obviously broken' range, not in the 'something is
    catastrophically wrong' range (hundreds of sigma, which WOULD mean a units
    or exposure bug rather than model inadequacy or real structure)."""
    from reduce_totani import TOTANI_EDGES_GEV
    from totani.stack import build_background_stack

    stack = build_background_stack(RUN, TOTANI_EDGES_GEV[:10])
    result = fit_all_bins(stack)
    assert len(result["per_bin"]) == 9
    assert len(result["component_order"]) == 19
    for row in result["per_bin"]:
        assert row["converged"], f"bin {row['bin']} did not converge"
        assert all(a > 0 for a in row["amplitudes"].values())
        ratio = row["predicted_total"] / row["observed_total"]
        assert 0.2 < ratio < 5.0, (row["bin"], ratio)
        assert row["max_abs_sigma_halo"] < 100.0, (
            f"bin {row['bin']}: {row['max_abs_sigma_halo']} sigma in the halo "
            "is large enough to suspect a units/exposure bug, not just an "
            "imperfect background model")
