"""Template-detection statistics in the convention Oscar set out on 2026-09-29
(Macias et al. 2018, Nature Astronomy 2, 387, Methods "Testing Extended
Emission Templates" and Supplementary Sec. 2, Eq. 7; Cowan et al. 2011).

For one fixed spatial template with a NON-NEGATIVE amplitude in each of K
energy bins, with every nuisance parameter refitted in both hypotheses,

    q_e^+ = 2 [max_{a_e >= 0} ln L_e - ln L_e(a_e = 0)],      Q^+ = sum_e q_e^+,

and under independent bins and regular nuisances the null distribution is the
binomial chi-square mixture

    f_K(q) = 2^-K delta(q) + 2^-K sum_{j=1..K} C(K, j) f_{chi2_j}(q),

with upper-tail probability p_K(q) = 2^-K sum_j C(K, j) S_{chi2_j}(q) for q > 0.
For a concave amplitude likelihood, a bin whose signed optimum is negative
contributes q_e^+ = 0. Two significance conventions:

    Z_2018 = sqrt(F^-1_{chi2_1}(1 - p)) = Phi^-1(1 - p/2)     (two-sided; the 2018 paper)
    Z_one  = Phi^-1(1 - p)                                   (one-sided)

Everything is evaluated in log space (log-survival functions, log-sum-exp)
so that large Q^+ do not underflow. A signed-amplitude sum of improvements
tests departures in either direction; its reference is chi2_K, not f_K.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import gammaln, logsumexp, ndtri_exp
from scipy.stats import chi2

LN2 = math.log(2.0)


def q_plus(ts_signed: float, a_signed: float) -> float:
    """Positive-emission contribution of one bin from a signed fit: the signed
    TS if the signed optimum is >= 0, else 0 (concave amplitude likelihood)."""
    return float(max(ts_signed, 0.0)) if a_signed >= 0 else 0.0


def _chi2_logsf(q: float, dof) -> np.ndarray:
    """ln S_{chi2_dof}(q), with an asymptotic continuation where scipy underflows:
    ln Gamma(a, x) = (a-1) ln x - x + ln sum_k (a-1)...(a-k) / x^k, a = dof/2, x = q/2."""
    dof = np.atleast_1d(np.asarray(dof, float))
    out = chi2.logsf(q, dof)
    bad = ~np.isfinite(out)
    if bad.any():
        a = dof[bad] / 2.0; x = q / 2.0
        term = np.ones_like(a); ssum = np.ones_like(a)
        for k in range(1, 30):
            term = term * (a - k) / x
            ssum = ssum + term
            if np.all(np.abs(term) < 1e-17 * np.abs(ssum)):
                break
        out[bad] = (a - 1) * np.log(x) - x + np.log(ssum) - gammaln(a)
    return out


def log_p_mixture(q: float, K: int) -> float:
    """ln p_K(q) for the binomial chi-square mixture; q > 0. Returns 0 (p = 1) at q <= 0."""
    if q <= 0:
        return 0.0
    j = np.arange(1, K + 1)
    log_binom = gammaln(K + 1) - gammaln(j + 1) - gammaln(K - j + 1)
    return float(-K * LN2 + logsumexp(log_binom + _chi2_logsf(q, j)))


def log_p_chi2(q: float, dof: int) -> float:
    return float(_chi2_logsf(max(q, 0.0), dof)[0])


def z_2018(log_p: float) -> float:
    """Phi^-1(1 - p/2): the two-sided / chi2_1-equivalent convention of Macias et al. 2018, Eq. 3."""
    return float(-ndtri_exp(log_p - LN2))


def z_one(log_p: float) -> float:
    """Phi^-1(1 - p): the one-sided convention."""
    return float(-ndtri_exp(log_p))


def threshold(z: float, K: int, convention: str = "2018") -> float:
    """The Q^+ at which the mixture reaches significance z (bisection in log space)."""
    zf = z_2018 if convention == "2018" else z_one
    lo, hi = 0.0, 10.0
    while zf(log_p_mixture(hi, K)) < z:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if zf(log_p_mixture(mid, K)) < z:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def summarize(ts_signed, a_signed, bins=None) -> dict:
    """Per-bin q^+, Q^+, K, p_K and both Z for the bins used (default all);
    also the signed sum with its chi2_K reference."""
    ts = np.asarray(ts_signed, float); a = np.asarray(a_signed, float)
    idx = list(range(len(ts))) if bins is None else list(bins)
    qp = [q_plus(ts[i], a[i]) for i in idx]
    Q = float(sum(qp)); K = len(idx)
    lp = log_p_mixture(Q, K)
    signed = float(sum(ts[i] for i in idx))
    lps = log_p_chi2(signed, K)
    return {"K": K, "bins": idx, "q_plus": qp, "Q_plus": Q, "log10_p": lp / math.log(10),
            "Z_2018": z_2018(lp) if Q > 0 else 0.0, "Z_one": z_one(lp) if Q > 0 else 0.0,
            "signed_sum": signed, "signed_log10_p_chi2K": lps / math.log(10)}


def self_check() -> dict:
    """Oscar's numerical checks: K=1 gives p = S_chi2_1 / 2; 2018 4-sigma thresholds
    38.422 (K=19), 27.616 (K=9); one-sided 4-sigma for K=9 is 29.235."""
    q = 7.3
    k1 = abs(math.exp(log_p_mixture(q, 1)) - 0.5 * chi2.sf(q, 1))
    return {"K1_residual": k1, "thr_2018_K19": threshold(4.0, 19), "thr_2018_K9": threshold(4.0, 9),
            "thr_one_K9": threshold(4.0, 9, "one")}
