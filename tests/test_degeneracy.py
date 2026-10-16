"""Leverage and the latitude lever arm — the linear algebra, on synthetic bases.

These are fast and use no data on purpose: the claim being defended is a
property of least squares, not of this dataset, and it should fail loudly if
someone changes the convention (amplitude ratio vs variance ratio) or the
partition rule.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from totani.degeneracy import lever_arm_share, leverage, residual_power  # noqa: E402


def test_template_inside_the_span_has_no_leverage():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(50, 4))
    h = A @ np.array([1.0, -2.0, 0.5, 3.0])
    assert leverage(h, A) == pytest.approx(0.0, abs=1e-10)


def test_template_orthogonal_to_the_span_has_full_leverage():
    A = np.zeros((6, 2))
    A[0, 0] = A[1, 1] = 1.0
    h = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
    assert leverage(h, A) == pytest.approx(1.0)


def test_leverage_is_the_amplitude_ratio_not_the_power_ratio():
    """Half the power orthogonal is leverage 1/sqrt(2), not 1/2."""
    A = np.array([[1.0], [0.0]])
    h = np.array([1.0, 1.0])
    num, den = residual_power(h, A)
    assert num / den == pytest.approx(0.5)
    assert leverage(h, A) == pytest.approx(np.sqrt(0.5))


def test_residual_is_orthogonal_to_every_basis_vector():
    rng = np.random.default_rng(1)
    A = rng.normal(size=(40, 5))
    h = rng.normal(size=40)
    coef, *_ = np.linalg.lstsq(A, h, rcond=None)
    r = h - A @ coef
    assert np.allclose(A.T @ r, 0.0, atol=1e-10)


def test_splitting_a_fit_can_only_reduce_the_residual():
    """The identity the lever arm rests on: share is in [0, 1)."""
    rng = np.random.default_rng(2)
    for _ in range(20):
        A = rng.normal(size=(60, 6))
        h = rng.normal(size=60)
        m = np.zeros(60, bool)
        m[: rng.integers(10, 50)] = True
        out = lever_arm_share(h, A, [m, ~m])
        assert 0.0 <= out["share"] < 1.0
        assert out["leverage_separate"] <= out["leverage_common"] + 1e-12


def test_lever_arm_is_zero_when_regions_cannot_disagree():
    """Two regions that are exact copies impose no extra constraint."""
    rng = np.random.default_rng(3)
    A0 = rng.normal(size=(25, 3))
    h0 = rng.normal(size=25)
    A = np.vstack([A0, A0])
    h = np.concatenate([h0, h0])
    m = np.zeros(50, bool)
    m[:25] = True
    assert lever_arm_share(h, A, [m, ~m])["share"] == pytest.approx(0.0, abs=1e-10)


def test_lever_arm_is_large_when_regions_demand_opposite_coefficients():
    """Same basis vector, opposite sign relative to h in the two halves."""
    n = 20
    A = np.zeros((2 * n, 1))
    A[:n, 0] = 1.0
    A[n:, 0] = 1.0
    h = np.concatenate([np.ones(n), -np.ones(n)])
    m = np.zeros(2 * n, bool)
    m[:n] = True
    out = lever_arm_share(h, A, [m, ~m])
    # separately, each half is perfectly absorbed; together, neither is
    assert out["leverage_separate"] == pytest.approx(0.0, abs=1e-10)
    assert out["leverage_common"] == pytest.approx(1.0)
    assert out["share"] == pytest.approx(1.0)


def test_masks_must_partition():
    A = np.eye(4)
    h = np.ones(4)
    overlap = [np.array([1, 1, 0, 0], bool), np.array([0, 1, 1, 1], bool)]
    gap = [np.array([1, 0, 0, 0], bool), np.array([0, 1, 0, 0], bool)]
    for bad in (overlap, gap):
        with pytest.raises(ValueError, match="partition"):
            lever_arm_share(h, A, bad)


def test_weights_change_the_answer_and_are_applied_to_both_sides():
    rng = np.random.default_rng(4)
    A = rng.normal(size=(30, 3))
    h = rng.normal(size=30)
    w = rng.uniform(0.5, 2.0, size=30)
    # weighting by a CONSTANT cannot change a ratio
    assert leverage(h, A, np.full(30, 3.0)) == pytest.approx(leverage(h, A))
    # a non-constant weight generally does
    assert leverage(h, A, w) != pytest.approx(leverage(h, A))
