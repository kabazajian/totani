"""How much of the halo template the background basis cannot absorb.

    from totani.degeneracy import leverage, lever_arm_share

The halo amplitude in this analysis is not determined by the halo template. It
is determined by the part of the halo template that is ORTHOGONAL to the span
of the background components, because everything inside that span is absorbed
by the free per-bin normalizations at no cost to the likelihood. `leverage`
measures that part:

    leverage = ||h - A c|| / ||h||,   c = argmin ||h - A c||

so leverage = 1 for a template no background can imitate, and leverage -> 0 as
the halo becomes a linear combination of backgrounds.

🚨 IT IS AN AMPLITUDE RATIO, NOT A VARIANCE RATIO. The variance of the fitted
signal amplitude scales as 1/||h_perp||^2, so relative to an orthogonal signal
the variance inflation is 1/leverage^2 -- at the measured 0.266 that is 14x,
which is the number quoted in the Galp21 proposal. Reporting leverage^2 (0.071)
as "7% orthogonal" would describe the same geometry and mislead about the
error bar, so the convention here is the amplitude ratio throughout.

WHY `lever_arm_share` EXISTS. Leverage computed on the whole ROI is larger
than leverage computed on any sub-region of it, and that is not a paradox: the
fit uses ONE coefficient vector for all latitudes, so a background combination
that imitates the halo at 10-30 deg must ALSO imitate it at 30-60 deg with the
same normalizations. It usually cannot do both. That constraint -- not the
halo's brightness anywhere in particular -- is where the constraining power
comes from, and `lever_arm_share` measures what fraction of the total
orthogonal power it supplies.

The identity that makes the comparison fair is that a common-coefficient least
squares residual can never be smaller than the sum of independently-fitted
per-region residuals (the separate fit has strictly more freedom), so the share
is bounded in [0, 1) by construction. `tests/test_degeneracy.py` asserts it.
"""
from __future__ import annotations

import numpy as np

__all__ = ["residual_power", "leverage", "lever_arm_share"]


def _apply_weights(h, A, weights):
    if weights is None:
        return np.asarray(h, float), np.asarray(A, float)
    w = np.asarray(weights, float)
    return np.asarray(h, float) * w, np.asarray(A, float) * w[:, None]


def residual_power(h, A, weights=None) -> tuple:
    """`(||h_perp||^2, ||h||^2)` for template `h` against basis `A` (n_cell, m).

    `weights` multiplies both sides elementwise -- pass `1/sqrt(counts)` for the
    Poisson-weighted version, or None for the plain least-squares one the
    proposal quotes.
    """
    hw, Aw = _apply_weights(h, A, weights)
    if Aw.ndim != 2 or Aw.shape[0] != hw.shape[0]:
        raise ValueError(f"basis {Aw.shape} does not match template {hw.shape}")
    coef, *_ = np.linalg.lstsq(Aw, hw, rcond=None)
    r = hw - Aw @ coef
    return float(r @ r), float(hw @ hw)


def leverage(h, A, weights=None) -> float:
    """`||h_perp|| / ||h||` -- the fraction of the template's AMPLITUDE that no
    combination of the basis can reproduce. See the module docstring."""
    num, den = residual_power(h, A, weights)
    if den <= 0:
        raise ValueError("template has zero norm; leverage is undefined")
    return float(np.sqrt(num / den))


def lever_arm_share(h, A, region_masks, weights=None) -> dict:
    """What fraction of the orthogonal power comes from sharing coefficients.

    `region_masks` is a sequence of boolean arrays over the same cells as `h`;
    they must PARTITION those cells, because the quantity being computed is a
    comparison between one fit over the union and independent fits over the
    parts. A gap or an overlap would silently compare different amounts of sky.

    Returns the full-ROI and per-region leverages, the leverage the same data
    would have if each region were fitted with its own normalizations, and

        share = (||r_common||^2 - sum_i ||r_i||^2) / ||r_common||^2

    the part of the constraint that exists ONLY because one set of
    normalizations has to serve every region at once.
    """
    masks = [np.asarray(m, bool) for m in region_masks]
    if len(masks) < 2:
        raise ValueError("need at least two regions to have a lever arm")
    stacked = np.stack(masks)
    if not np.all(stacked.sum(axis=0) == 1):
        raise ValueError("region masks must partition the cells exactly "
                         "(every cell in exactly one region)")

    r_full, h_full = residual_power(h, A, weights)
    parts = []
    for m in masks:
        w = None if weights is None else np.asarray(weights, float)[m]
        parts.append(residual_power(np.asarray(h)[m], np.asarray(A)[m], w))
    r_sep = sum(p[0] for p in parts)
    h_sep = sum(p[1] for p in parts)
    return {
        "leverage_common": float(np.sqrt(r_full / h_full)),
        "leverage_separate": float(np.sqrt(r_sep / h_sep)),
        "leverage_per_region": [float(np.sqrt(a / b)) for a, b in parts],
        "share": float((r_full - r_sep) / r_full),
    }
