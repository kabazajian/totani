"""Oscar's numerical checks of the Macias et al. 2018 mixture (note of 2026-09-29)."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scipy.stats import chi2

from totani import stats


def test_k1_is_half_chi2_1():
    for q in (0.5, 3.0, 25.0):
        assert abs(math.exp(stats.log_p_mixture(q, 1)) - 0.5 * chi2.sf(q, 1)) < 1e-15


def test_thresholds():
    assert abs(stats.threshold(4.0, 19) - 38.422) < 1e-3
    assert abs(stats.threshold(4.0, 9) - 27.616) < 1e-3
    assert abs(stats.threshold(4.0, 9, "one") - 29.235) < 1e-3


def test_q_plus_negative_bin_is_zero():
    assert stats.q_plus(12.0, -3.0) == 0.0 and stats.q_plus(12.0, 3.0) == 12.0


def test_large_q_no_underflow():
    lp = stats.log_p_mixture(2000.0, 9)
    assert lp < -400 and math.isfinite(stats.z_2018(lp))


def test_asymptotic_matches_scipy_where_both_work():
    import numpy as np
    for q, d in ((600.0, 3), (900.0, 9), (1200.0, 1)):
        ref = chi2.logsf(q, d)
        a, x = d / 2, q / 2
        from scipy.special import gammaln
        term, ssum = 1.0, 1.0
        for k in range(1, 30):
            term *= (a - k) / x; ssum += term
        asym = (a - 1) * np.log(x) - x + np.log(ssum) - gammaln(a)
        assert abs(asym - ref) < 1e-8 * abs(ref)
