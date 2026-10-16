"""The signed-signal fit: background (non-negative) + one signal amplitude
that Totani explicitly allows to be negative (sec 2.3).

🚨 `test_joint_lbfgsb_conditioning_bug_is_the_reason_for_alternation` documents
the bug that motivated the whole design: a single L-BFGS-B call over both
parameter blocks recovered essentially ZERO signal for injections of +80, -40
AND 0 alike, "converged" via the ftol stopping rule with a genuine leftover
gradient. It is not re-tested directly (that code no longer exists), but the
recovery test below is what would have caught it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from totani.likelihood import (  # noqa: E402
    _fit_signal_given_bg,
    _mu_bg,
    expected_counts,
    fit_all_bins_with_signal,
    fit_bin_with_signal,
    neg2lnl_with_signal,
    profile_signal_uncertainty,
)


def _synthetic(seed=0, n_bg=5, n_cell=12, true_sig=80.0):
    rng = np.random.default_rng(seed)
    T_bg = rng.uniform(0.1, 10.0, size=(n_bg, n_cell, n_cell))
    T_sig = rng.uniform(0.1, 5.0, size=(n_cell, n_cell))
    true_bg = rng.uniform(10.0, 500.0, size=n_bg)
    mu_true = expected_counts(np.log10(true_bg), T_bg) + true_sig * T_sig
    assert (mu_true > 0).all(), "test fixture produced negative mu; adjust true_sig"
    n_obs = rng.poisson(mu_true).astype(float)
    return T_bg, T_sig, n_obs, true_bg, true_sig


@pytest.mark.parametrize("true_sig", [80.0, -40.0, 0.0])
def test_recovers_positive_negative_and_zero_signal(true_sig):
    """THE test the conditioning bug failed. A signal roughly this size against
    this many counts should be recovered to within a modest factor of Poisson
    noise -- not silently absorbed to exactly zero regardless of the true
    injected value, which is what the joint L-BFGS-B version did."""
    T_bg, T_sig, n_obs, true_bg, _ = _synthetic(true_sig=true_sig)
    val, bg, sig, r = fit_bin_with_signal(T_bg, n_obs, T_sig)
    # loose tolerance: this is ONE noisy realization on an unstructured toy
    # template set, not a precision check -- the point is "moved substantially
    # toward the truth", not "matches it exactly".
    assert abs(sig - true_sig) < 30.0, (true_sig, sig)
    assert (bg > 0).all()


def test_fit_bin_with_signal_val_matches_neg2lnl_with_signal():
    """Two independently-written objective evaluations (the alternation loop's
    running total, and the standalone checker) must agree at the reported
    optimum -- if they disagree, one of them has a units/formula bug."""
    T_bg, T_sig, n_obs, _, _ = _synthetic(seed=1)
    val, bg, sig, r = fit_bin_with_signal(T_bg, n_obs, T_sig)
    params = np.concatenate([np.log10(bg), [sig]])
    val2 = neg2lnl_with_signal(params, T_bg, T_sig, n_obs)
    assert val == pytest.approx(val2, rel=1e-6)


def test_1d_newton_solve_reaches_a_first_order_stationary_point():
    """`_fit_signal_given_bg` is convex in the signal amplitude; the analytic
    first derivative at its answer must be ~0 (allowing for the clip at
    SIGNAL_BOUND, which is not reached in this fixture)."""
    T_bg, T_sig, n_obs, true_bg, _ = _synthetic(seed=2)
    bg_contribution = _mu_bg(np.log10(true_bg), T_bg)
    a = _fit_signal_given_bg(T_sig, n_obs, bg_contribution, a0=0.0)
    mu = np.maximum(bg_contribution + a * T_sig, 1e-30)
    grad = 2.0 * np.sum(T_sig * (1.0 - n_obs / mu))
    assert abs(grad) < 1e-4, grad


def test_profile_interval_brackets_the_best_fit_and_has_the_right_delta():
    T_bg, T_sig, n_obs, _, _ = _synthetic(seed=3)
    val, bg, sig, r = fit_bin_with_signal(T_bg, n_obs, T_sig)
    lo, hi = profile_signal_uncertainty(T_bg, n_obs, T_sig, sig, bg)
    assert lo < sig < hi

    def neg2lnl_at(a):
        from totani.likelihood import _fit_bg_given_offset

        theta, _ = _fit_bg_given_offset(T_bg, n_obs, a * T_sig, np.log10(bg))
        mu = _mu_bg(theta, T_bg) + a * T_sig
        from scipy.special import xlogy
        return 2.0 * (mu.sum() - xlogy(n_obs, mu + 1e-30).sum())

    val0 = neg2lnl_at(sig)
    assert neg2lnl_at(lo) - val0 == pytest.approx(1.0, abs=0.05)
    assert neg2lnl_at(hi) - val0 == pytest.approx(1.0, abs=0.05)


def test_fit_all_bins_with_signal_shapes_and_keys():
    from totani.cells import halo_cell_mask

    rng = np.random.default_rng(4)
    n_bins, n_comp = 3, 4
    names = [f"comp{i}" for i in range(n_comp)]
    templates = {n: rng.uniform(1, 10, size=(n_bins, 12, 12)) for n in names}
    true_amp = {n: rng.uniform(10, 100) for n in names}
    signal = rng.uniform(0.5, 5, size=(n_bins, 12, 12))
    mu = sum(true_amp[n] * templates[n] for n in names) + 30.0 * signal
    counts = rng.poisson(mu).astype(float)
    stack = {"templates": templates, "counts": counts, "n_bins": n_bins}

    result = fit_all_bins_with_signal(stack, signal, component_order=names,
                                      with_uncertainty=True)
    assert len(result["per_bin"]) == n_bins
    for row in result["per_bin"]:
        assert set(row["bg_amplitudes"]) == set(names)
        assert "signal_amplitude" in row
        assert "signal_amplitude_lo" in row and "signal_amplitude_hi" in row
        assert row["signal_amplitude_lo"] <= row["signal_amplitude"] <= row["signal_amplitude_hi"]


def test_signal_amplitude_can_be_negative():
    """The whole point: unlike every background component, this parameter is
    not constrained to be positive."""
    T_bg, T_sig, n_obs, _, _ = _synthetic(seed=5, true_sig=-40.0)
    _, _, sig, _ = fit_bin_with_signal(T_bg, n_obs, T_sig)
    assert sig < 0


# ---------------------------------------------------------------------------
# The profile estimator and its acceptance gate (added 2026-09-04 after the
# audit found block alternation biased low on real templates).
# ---------------------------------------------------------------------------

def test_injection_slope_is_scale_free():
    """🚨 REGRESSION FOR A BROKEN GATE. `injection_slope` first used
    `lstsq` on a [injections, ones] design matrix. Real amplitudes are ~1e-33,
    so that matrix has a singular-value ratio ~1e-33 -- below lstsq's rcond
    cutoff, which truncates the slope direction and returns exactly 0.0. The
    gate reported slope = 0.000 for an estimator independently known to be at
    0.35: it would have condemned a working estimator and could equally have
    passed a broken one. The covariance form is scale-free."""
    from totani.likelihood import injection_slope

    rng = np.random.default_rng(11)
    n_bins = 2
    templates = {"a": rng.uniform(1, 5, size=(n_bins, 12, 12)),
                "b": rng.uniform(1, 5, size=(n_bins, 12, 12))}
    signal = rng.uniform(0.5, 2.0, size=(n_bins, 12, 12))
    counts = (40.0 * templates["a"] + 25.0 * templates["b"]).round()
    stack = {"templates": templates, "counts": counts, "n_bins": n_bins}

    # A tiny-amplitude scale is exactly what broke the lstsq version.
    for scale in (1.0, 1e-33):
        sig = signal / scale          # so amplitudes come out ~scale
        out = injection_slope(stack, sig, [0.0, 10.0 * scale, 30.0 * scale])
        assert all(abs(s - 1.0) < 0.1 for s in out["slopes"]), (scale, out["slopes"])


def test_profile_estimator_recovers_injections_at_unit_slope():
    """THE ACCEPTANCE GATE. A correct linear-model MLE recovers an injected
    component at slope ~1 no matter how collinear the templates are --
    collinearity inflates the VARIANCE, not the bias. Measured on the real
    reduction 2026-09-04: profile 0.998-1.012, block alternation 0.355-0.636.
    """
    from totani.likelihood import fit_all_bins_profile, injection_slope

    rng = np.random.default_rng(12)
    n_bins = 2
    # deliberately near-degenerate templates: the signal is mostly, but not
    # entirely, inside the background span
    base = rng.uniform(1, 5, size=(n_bins, 12, 12))
    templates = {"a": base, "b": base * 1.05 + rng.uniform(0, 0.4, size=base.shape)}
    signal = 0.8 * base + rng.uniform(0, 0.5, size=base.shape)
    counts = (60.0 * templates["a"] + 30.0 * templates["b"]).round()
    stack = {"templates": templates, "counts": counts, "n_bins": n_bins}

    out = injection_slope(stack, signal, [0.0, 20.0, 60.0], fitter=fit_all_bins_profile)
    for s in out["slopes"]:
        assert abs(s - 1.0) < 0.15, out["slopes"]


def test_profile_and_alternation_agree_when_well_conditioned():
    """Where the templates are NOT degenerate the two estimators must give the
    same answer -- establishing that the profile is not simply a different
    estimator, but the same MLE reached reliably."""
    from totani.likelihood import fit_bin_profile

    T_bg, T_sig, n_obs, _, true_sig = _synthetic(seed=13, true_sig=80.0)
    _, _, sig_alt, _ = fit_bin_with_signal(T_bg, n_obs, T_sig)
    _, _, sig_prof, _ = fit_bin_profile(T_bg, n_obs, T_sig)
    assert abs(sig_prof - sig_alt) / max(abs(sig_alt), 1.0) < 0.1, (sig_alt, sig_prof)


def test_profile_curve_is_returned_and_brackets_the_minimum():
    from totani.likelihood import fit_bin_profile

    T_bg, T_sig, n_obs, _, _ = _synthetic(seed=14)
    val, bg, sig, curve = fit_bin_profile(T_bg, n_obs, T_sig)
    assert len(curve["grid"]) == len(curve["neg2lnl"])
    assert min(curve["grid"]) < sig < max(curve["grid"])
    assert val <= min(curve["neg2lnl"]) + 1e-6
