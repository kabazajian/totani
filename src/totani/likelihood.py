"""Poisson likelihood over 12x12 cells, one energy bin at a time — Totani sec 2.1/2.2.

    from totani.stack import build_background_stack
    from totani.likelihood import fit_bin, fit_all_bins

    stack = build_background_stack(run_dir, edges_gev)
    result = fit_all_bins(stack)

🚨 BACKGROUND ONLY. NO HALO TEMPLATE ENTERS THIS MODULE, deliberately. This is
the pipeline's self-check -- do the 19 named background components alone
describe the data sensibly? -- and it stays that way until a separate, explicit
step adds the halo component. Keeping the two apart means this module can be
imported, run and committed without ever being the thing that crosses the rule
6 boundary that `scripts/hooks/rule6_halo_guard.py` enforces.

THE PARAMETERISATION IS `a_l = 10**theta_l`, NOT A BOUNDED LINEAR AMPLITUDE.
This is gce-fisher's own established convention
(`scripts/a2020_stage2_fit.py::neg2lnl`), reused rather than reinvented, for
three reasons that transfer directly to this project:

  * non-negativity by CONSTRUCTION -- no bounds needed on the coefficients,
    which matters here because every one of the 19 background templates is a
    "known" component with f_l > 0 (Totani sec 2.3), so this is exactly the
    right shape of constraint;
  * scale invariance -- a template's absolute units (recall `totani.stack`'s
    own note: gas/dust carry gas COLUMN DENSITY, not gamma-ray flux; ICS has no
    stated units at all) never has to be normalised away, because theta simply
    sits at a different value;
  * well-conditioned gradients, a decade in amplitude being one unit in theta.

MATCHING THE SAME LIKELIHOOD FORM. `xlogy(n, mu)` is exactly zero where the
data is zero (empty cells contribute nothing without being dropped from the
sum), and the small epsilon floor on `mu` prevents a literal log(0) at the
boundary of the search. Constants in n (`ln(n!)`) are dropped throughout,
matching `a2020_stage2_fit.py`'s own stated convention ("TS is the anchor,
absolute -lnL is not") -- this module's `-2lnL` is for COMPARISON, never an
absolute number to quote on its own.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import xlogy

__all__ = ["expected_counts", "neg2lnl", "fit_bin", "fit_all_bins",
           "MAX_RESTARTS", "THETA_BOUNDS"]

#: Same bounds as the established A2020 fit: generous enough that a real
#: solution never sits at the edge, tight enough to prevent overflow in
#: `10**theta`.
THETA_BOUNDS = (-8.0, 14.0)

#: Re-enter L-BFGS-B until it actually reports convergence, not once blindly.
#: `a2020_stage2_fit.py`'s own note applies unchanged: gtol tight enough to be
#: meaningful is not reachable in one pass for a function of this magnitude.
MAX_RESTARTS = 8

EPS = 1e-30


def _poisson_val_grad_mu(mu, n):
    """(-2lnL, d/dmu) with a CONSISTENT clamp for non-positive mu.

    🚨 THREE WRONG VERSIONS PRECEDED THIS ONE, and the third was upstream's
    wrong version too. A signed template (dust_neg) times a large trial
    amplitude can drive mu negative during a LINE SEARCH even though the true
    optimum never does. What happens next depends on how the objective treats
    that region, and each naive treatment fails differently:

      1. `xlogy(n, mu + EPS)` -- log of a negative number is NaN; L-BFGS-B
         aborts ABNORMAL/nit=0 the moment a trial step crosses. (First
         version here; bins 0-6 of the real fit all died this way.)
      2. clamping ONLY the log term, keeping raw `mu.sum()` -- now the
         optimizer is REWARDED for negative mu (the sum term plunges while
         the log penalty is frozen at EPS): measured 2026-09-04, one L-BFGS-B
         step jumped to the bounds with -2lnL = -2e14 against a true optimum
         of -1.07e8. Unbounded below. (Second version here.)
      3. clamping the VALUE while leaving the GRADIENT computed from raw mu --
         upstream's own attempt, recorded in `a2020_stage2_fit.py`: "neg2lnl's
         clip keeps the value finite while making the GRADIENT WRONG, and
         L-BFGS-B ran away to -2lnL = 1.5e18".

    The consistent treatment: clamp mu in BOTH terms, and zero the gradient in
    clamped cells (the exact subgradient of the clamped objective). A deeply
    infeasible configuration then scores 2*sum(n)*|ln EPS| -- enormous and BAD
    -- with no spurious pull further in, and every feasible cell still
    contributes its correct Poisson gradient to steer the search back out.
    """
    mus = np.maximum(mu, EPS)
    val = 2.0 * (mus.sum() - xlogy(n, mus).sum())
    dval_dmu = np.where(mu > EPS, 2.0 * (1.0 - n / mus), 0.0)
    return val, dval_dmu


def expected_counts(theta, T) -> np.ndarray:
    """mu_i = sum_l 10**theta_l * T_li, flattened over cells."""
    return np.tensordot(10.0 ** np.asarray(theta), T, axes=(0, 0))


def neg2lnl(theta, T, n):
    """-2lnL and its analytic gradient, amplitudes as 10**theta.

    `T` is `(n_comp, n_cell, n_cell)` or already flattened to `(n_comp,
    n_cell**2)`; `n` matches its trailing shape. Returns `(value, grad)` for
    `scipy.optimize.minimize(..., jac=True)`.
    """
    a = 10.0 ** np.asarray(theta)
    mu = np.tensordot(a, T, axes=(0, 0))
    val, dval_dmu = _poisson_val_grad_mu(mu, n)
    axes = (tuple(range(1, T.ndim)), tuple(range(n.ndim)))
    grad_a = np.tensordot(T, dval_dmu, axes=axes)
    grad = np.log(10.0) * a * grad_a
    return val, grad


#: Linear-space polish rounds after the theta fit (see `_polish_linear`).
POLISH_ROUNDS = 2


def _polish_linear(T_scaled, n_bin, a0, offset=0.0):
    """L-BFGS-B in LINEAR amplitude space, bounds [0, 1e14], from `a0`.

    🚨 WHY THE THETA FIT ALONE IS NOT ENOUGH (found 2026-09-22). In the
    log10 parameterisation the gradient is d/dtheta = ln10 * a * d/da: once
    the optimizer has driven a component small, its gradient VANISHES even
    where d/da is large, and L-BFGS-B reports convergence ("relative
    reduction of F") on a plateau. Measured on the real reduction: a 5-ring
    row with the A2020 dust templates added sat 207-296 ABOVE its own fit
    without them (adding free templates cannot raise a convex minimum); the
    residual-construction pair on 5-ring rows sat 17-62 above the minimum in
    single 1.5-7 GeV bins, in one of the two sign choices; 1-gas rows were
    clean to 0.01. The problem is convex in `a` (Poisson, mu linear in a), so
    a descent in a-space from the theta solution has no vanishing-gradient
    trap; it recovered every case to the restart-verified minimum. Same
    objective, same clamp (`_poisson_val_grad_mu`), `offset` as in
    `_fit_bg_given_offset`.
    """
    T_flat = T_scaled.reshape(T_scaled.shape[0], -1)
    n_flat = np.asarray(n_bin).reshape(-1)
    off = np.asarray(offset, float).reshape(-1) if np.ndim(offset) else offset

    def obj(a):
        mu = a @ T_flat + off
        val, dmu = _poisson_val_grad_mu(mu, n_flat)
        return val, T_flat @ dmu

    kw = dict(jac=True, method="L-BFGS-B", bounds=[(0.0, 1e14)] * len(a0),
              options={"maxiter": 4000, "ftol": 1e-15, "gtol": 1e-8})
    r = minimize(obj, np.clip(a0, 0.0, 1e14), **kw)
    for _ in range(MAX_RESTARTS):
        if r.status == 0:
            break
        r = minimize(obj, r.x, **kw)
    return r.fun, r.x, r


def _fit_theta_then_polish(objective_theta, T_scaled, n_bin, x0, kw, offset=0.0):
    """The theta fit, then `POLISH_ROUNDS` of (linear polish, theta refit),
    keeping the lowest -2lnL seen. Returns `(val, a, r, nfev)` with `r.x`,
    `r.fun` set to the returned solution."""
    r = minimize(objective_theta, x0, **kw)
    nfev = r.nfev
    for _ in range(MAX_RESTARTS):
        if r.status == 0:
            break
        r = minimize(objective_theta, r.x, **kw)
        nfev += r.nfev
    best_val, best_a = r.fun, 10.0 ** r.x
    for _ in range(POLISH_ROUNDS):
        v2, a2, r2 = _polish_linear(T_scaled, n_bin, best_a, offset)
        nfev += r2.nfev
        if v2 < best_val - 1e-9:
            best_val, best_a = v2, a2
        else:
            break
        x1 = np.clip(np.log10(np.maximum(best_a, 1e-300)), *THETA_BOUNDS)
        r3 = minimize(objective_theta, x1, **kw)
        nfev += r3.nfev
        if r3.fun < best_val - 1e-9:
            best_val, best_a = r3.fun, 10.0 ** r3.x
    # the polish may park a component at exactly 0; report it at the theta
    # lower bound instead (1e-8 of a unit-L1 template: no change in -2lnL at
    # this precision), so amplitudes stay positive for every log10 downstream
    best_a = np.maximum(best_a, 10.0 ** THETA_BOUNDS[0])
    r.x = np.clip(np.log10(best_a), *THETA_BOUNDS)
    r.fun = best_val
    r.nfev = nfev
    return best_val, best_a, r, nfev


def fit_bin(T_bin, n_bin, x0=None):
    """Poisson fit of one energy bin. Returns `(neg2lnl_value, amplitudes, result)`.

    `T_bin`: `(n_comp, n_cell, n_cell)`. `n_bin`: `(n_cell, n_cell)`. `amplitudes`
    are relative to `T_bin` AS PASSED IN (`expected_counts(log10(amp), T_bin)`
    reproduces the fit) -- the internal rescaling below is invisible to callers.

    🚨 PER-TEMPLATE RESCALING IS NOT OPTIONAL HERE, unlike in
    `a2020_stage2_fit.py`, whose templates already share a common
    (gcepy-normalised) scale. This project's 19 background components span
    THIRTEEN ORDERS OF MAGNITUDE in basis sum -- measured 2026-09-04 on the
    real reduction, bin 0: point sources ~7.7e5, ics_ring6 ~6.8e5, loopI_shell1
    ~2.2e22. The fixed `THETA_BOUNDS = (-8, 14)` copied from that established
    convention is calibrated for ITS templates' scale, not this project's, and
    the first version of this function used it unscaled: bin 0 failed to
    converge, pinned against the lower bound, because loopI_shell1's true
    optimum needed an amplitude far below 1e-8 to avoid a ~100x over-prediction.

    Each template is divided by its own L1 norm (`sum(abs(T_bin[i]))`, the same
    measure already used for `x0`) before fitting, so every rescaled template
    has the SAME O(1) scale regardless of its original units -- exactly
    restoring the "a decade in amplitude is one unit in theta" regime the log10
    parameterisation is designed for. The recovered amplitude is divided back
    by the same factor before being returned, so nothing about this is visible
    to a caller multiplying `amp` against the ORIGINAL `T_bin`.
    """
    m = T_bin.shape[0]
    scale = np.maximum(np.array([np.abs(T_bin[i]).sum() for i in range(m)]), 1e-300)
    T_scaled = T_bin / scale[:, None, None]
    if x0 is None:
        # Every rescaled template has L1 norm 1, so an equal-share start is the
        # SAME theta for every component regardless of the original units --
        # this is what the rescaling buys, not just a numerical nicety.
        x0 = np.full(m, np.log10(max(n_bin.sum() / m, 1e-300)))
    kw = dict(args=(T_scaled, n_bin), method="L-BFGS-B", jac=True,
             bounds=[THETA_BOUNDS] * m,
             options={"maxiter": 4000, "ftol": 1e-15, "gtol": 1e-5})
    val, a, r, _ = _fit_theta_then_polish(lambda th: neg2lnl(th, T_scaled, n_bin), T_scaled, n_bin, x0,
                                          dict(kw, args=()))
    return val, a / scale, r


def fit_all_bins(stack, component_order=None, halo_mask=None) -> dict:
    """Fit every bin in `stack` independently (Totani sec 2.2: the maximum
    likelihood fitting IS independent per bin, not an approximation of a joint
    fit). Returns per-bin amplitudes, convergence, predicted-vs-observed
    totals, AND a per-cell goodness-of-fit split by disk/halo -- the
    background-only self-check this module exists for.

    🚨 TOTAL-FLUX AGREEMENT IS NOT GOODNESS OF FIT, AND THE FIRST VERSION OF
    THIS FUNCTION REPORTED ONLY THE FORMER. Measured 2026-09-04 on the real
    reduction, bin 0 (19 components, 144 cells): total predicted matched total
    observed to a relative 3e-8 -- which reads as an excellent fit -- while the
    per-cell residual (n - mu)/sqrt(mu) reached 38.6 SIGMA in the disk and 21.8
    sigma even restricted to the halo (mean 6.9 sigma there). A Poisson MLE with
    enough free normalizations gets the AGGREGATE right far more easily than
    the SHAPE; reporting only the aggregate would have looked like a clean pass
    while hiding exactly the kind of large, unexplained residual that needs
    chasing before this pipeline can support any claim about the data.
    `halo_mask` defaults to `cells.halo_cell_mask()` so the disk/halo split is
    always available without extra plumbing at the call site.
    """
    from .cells import halo_cell_mask

    names = component_order or sorted(stack["templates"])
    T = np.stack([stack["templates"][n] for n in names])   # (n_comp, n_bins, 12, 12)
    counts = stack["counts"]                                # (n_bins, 12, 12)
    n_bins = stack["n_bins"]
    halo = halo_cell_mask() if halo_mask is None else halo_mask

    # ⚠ NO CROSS-BIN WARM START. `fit_bin` rescales each template internally
    # (L1-norm 1) before fitting theta, and that rescaling factor is DIFFERENT
    # per bin (the same component's spatial map is regridded onto a different
    # exposure/spectral weighting per bin). A warm start built from bin k's
    # AMPLITUDES (original units) would need to be re-expressed in bin k+1's
    # own rescaled theta space to mean anything -- simpler and equally fast to
    # let each bin use its own default (equal-share) start, which the
    # convergence tests already show works reliably on its own.
    per_bin = []
    for k in range(n_bins):
        val, amp, r = fit_bin(T[:, k], counts[k])
        mu = expected_counts(np.log10(np.clip(amp, 1e-300, None)), T[:, k])
        sig = (counts[k] - mu) / np.sqrt(np.maximum(mu, 1.0))
        per_bin.append({
            "bin": k, "neg2lnl": val, "converged": bool(r.status == 0),
            "n_iter": int(r.nit), "amplitudes": dict(zip(names, amp)),
            "observed_total": float(counts[k].sum()),
            "predicted_total": float(mu.sum()),
            "max_abs_sigma_disk": float(np.abs(sig[~halo]).max()) if (~halo).any() else None,
            "mean_abs_sigma_halo": float(np.abs(sig[halo]).mean()),
            "max_abs_sigma_halo": float(np.abs(sig[halo]).max()),
        })
    return {"component_order": names, "per_bin": per_bin}


# ---------------------------------------------------------------------------
# Signal fit: 19 non-negative background amplitudes PLUS one SIGNED signal
# amplitude, jointly. This is what actually searches for the halo.
# ---------------------------------------------------------------------------

#: Symmetric bound on the signal amplitude, on a template pre-scaled to unit
#: L1 norm (matching the background rescaling in `fit_bin`).
SIGNAL_BOUND = 10.0 ** THETA_BOUNDS[1]

#: Block-alternation stopping rule and iteration cap. See `fit_bin_with_signal`
#: for why a single JOINT L-BFGS-B call over both parameter blocks is not used.
ALTERNATION_TOL = 1e-8
MAX_ALTERNATIONS = 30


def _mu_bg(theta_bg, T_bg):
    return np.tensordot(10.0 ** theta_bg, T_bg, axes=(0, 0))


def _fit_bg_given_offset(T_bg, n_bin, offset, x0):
    """Background-only Poisson MLE with a FIXED additive offset (the current
    signal contribution). Same objective as `neg2lnl`, plus `offset` -- the
    identical device `a2020_stage2_fit.py::neg2lnl` uses for "fixed model flux
    that is not a fitted component" (there, catalogued sources held at their
    4FGL prediction; here, the other block in an alternating solve)."""
    def obj(theta_bg):
        a = 10.0 ** theta_bg
        mu = np.tensordot(a, T_bg, axes=(0, 0)) + offset
        val, dval_dmu = _poisson_val_grad_mu(mu, n_bin)
        axes = (tuple(range(1, T_bg.ndim)), tuple(range(n_bin.ndim)))
        grad_a = np.tensordot(T_bg, dval_dmu, axes=axes)
        return val, np.log(10.0) * a * grad_a

    m = T_bg.shape[0]
    kw = dict(args=(), method="L-BFGS-B", jac=True, bounds=[THETA_BOUNDS] * m,
             options={"maxiter": 2000, "ftol": 1e-15, "gtol": 1e-7})
    _, _, r, _ = _fit_theta_then_polish(obj, T_bg, n_bin, x0, kw, offset=offset)
    return r.x, r


def _fit_signal_given_bg(T_sig, n_bin, bg_contribution, a0):
    """1-D Newton solve for the signed signal amplitude, background FIXED.

    `-2lnL(a) = 2 sum[(bg_i + a*T_i) - n_i ln(bg_i + a*T_i)]` is CONVEX in `a`
    wherever `mu > 0` (its second derivative `2 sum T_i^2 n_i / mu_i^2` is a sum
    of squares), so Newton's method converges fast and reliably from any start
    -- unlike the joint L-BFGS-B call this function exists to avoid, there is
    no scale-mismatch here: it is one parameter, solved on its own.
    """
    a = a0
    for _ in range(50):
        mu = bg_contribution + a * T_sig
        mu_safe = np.maximum(mu, EPS)
        g = 2.0 * np.sum(T_sig * (1.0 - n_bin / mu_safe))
        h = 2.0 * np.sum(T_sig ** 2 * n_bin / mu_safe ** 2)
        if h <= 0:
            break
        step = g / h
        a_new = a - step
        # Do not let a single Newton step push any cell's mu <= 0 -- halve the
        # step instead of jumping into the unphysical region the likelihood
        # itself would reject anyway.
        # 🚨 TEST ONLY CELLS THE SIGNAL CAN AFFECT (T_sig != 0). In a
        # cell-restricted stack (`restrict_stack`) the EXCLUDED cells have
        # mu = 0 by construction -- inert, not infeasible -- and the first
        # version's `np.any(mu <= 0)` read them as a violation and halved
        # EVERY Newton step ~50 times into nothing: steps of ~7e-16 toward a
        # true optimum of ~5e4, "converged" at zero each outer iteration.
        # Every disk-excluded result produced before 2026-09-04 -- the
        # all-zero quadrant fits, the frozen "upper bounds" identical across
        # configurations, the injection test recovering 0.00 of an injected
        # 200 -- was this one guard. Caught because the injection test
        # existed; none of it was visible from the real-data fits alone.
        affected = T_sig != 0
        while (np.any((bg_contribution + a_new * T_sig)[affected] <= 0)
               and abs(step) > 1e-15):
            step *= 0.5
            a_new = a - step
        if abs(a_new - a) < 1e-10 * max(abs(a), 1.0):
            a = a_new
            break
        a = a_new
    return float(np.clip(a, -SIGNAL_BOUND, SIGNAL_BOUND))


def neg2lnl_with_signal(params, T_bg, T_sig, n):
    """-2lnL for `mu = sum_l 10**theta_l * T_bg_l + a_sig * T_sig`, value only
    (no gradient -- kept for `profile_signal_uncertainty` and tests; the actual
    fit uses block alternation, not a joint gradient-based solve)."""
    m = T_bg.shape[0]
    theta_bg, a_sig = params[:m], params[m]
    mu = _mu_bg(theta_bg, T_bg) + a_sig * T_sig
    return _poisson_val_grad_mu(mu, n)[0]


def fit_bin_with_signal(T_bg, n_bin, T_sig, x0_bg=None, x0_sig=0.0):
    """Poisson fit of one bin: m non-negative background amplitudes plus one
    signed signal amplitude. Returns `(neg2lnl, bg_amplitudes, signal_amp, result)`,
    all relative to `T_bg`/`T_sig` AS PASSED IN.

    🚨 BLOCK ALTERNATION, NOT A SINGLE JOINT L-BFGS-B CALL. A first version put
    theta_bg (log-scale, non-negative) and a_sig (linear, signed) into ONE
    parameter vector for one `minimize` call. Measured 2026-09-04 on synthetic
    data with an injected signal of +80 (clearly detectable in principle): the
    fit recovered a_sig = 8e-6 -- essentially zero -- for injections of +80,
    -40 AND 0 alike. `r.status == 0` ("converged") but `r.message` was
    "RELATIVE REDUCTION OF F <= FACTR*EPSMCH" (the ftol stopping rule), with a
    LEFTOVER gradient of -0.019 in the signal direction -- L-BFGS-B's implicit
    identity-matrix initial Hessian guess is a poor match for two parameter
    blocks with very different local curvature scales (exponential vs linear).

    The fix sidesteps the conditioning problem rather than fighting it:
    alternate a background-only fit (`_fit_bg_given_offset`, holding the
    current signal fixed as an additive offset) against a 1-D Newton solve for
    the signal alone (`_fit_signal_given_bg`, convex, holding background
    fixed). Each sub-problem is on ONE scale, so neither optimizer ever has to
    reconcile the two.

    🚨 A SECOND, INDEPENDENT SCALE PROBLEM -- CAUGHT ON REAL DATA, NOT SYNTHETIC.
    Even with alternation, the FIRST version of this function fit `T_bg`/`T_sig`
    UNSCALED, and every bin came back `converged=False` with `sig_amp ~ 1e-16`
    (numerical noise) against the real reduction. `fit_bin` already carries the
    fix for exactly this (this project's templates span 13+ orders of
    magnitude in basis sum) -- but that fix was NOT carried over when this
    function was written, and the halo template makes it WORSE, not just
    present: measured 2026-09-04, `nfw_rho2_on_target`'s cell sums run
    ~5.8e42 to 3.7e44, twenty orders of magnitude past even the largest
    background component (loopI_shell1, ~2.2e22). `THETA_BOUNDS = (-8, 14)`
    cannot express an amplitude anywhere near what an UNSCALED halo template
    needs. Fixed the same way as `fit_bin`: every array -- background AND
    signal -- rescaled to unit L1 norm before fitting, unscaled after.
    """
    m = T_bg.shape[0]
    scale_bg = np.maximum(np.array([np.abs(T_bg[i]).sum() for i in range(m)]), 1e-300)
    scale_sig = max(np.abs(T_sig).sum(), 1e-300)
    Tb = T_bg / scale_bg[:, None, None]
    Ts = T_sig / scale_sig

    theta_bg = (np.full(m, np.log10(max(n_bin.sum() / (m + 1), 1e-300)))
               if x0_bg is None
               else np.log10(np.clip(np.asarray(x0_bg) * scale_bg, 1e-300, None)))
    a_sig = x0_sig * scale_sig
    prev_val = np.inf
    r = None
    for _ in range(MAX_ALTERNATIONS):
        offset = a_sig * Ts
        theta_bg, r = _fit_bg_given_offset(Tb, n_bin, offset, theta_bg)
        bg_contribution = _mu_bg(theta_bg, Tb)
        a_sig = _fit_signal_given_bg(Ts, n_bin, bg_contribution, a_sig)
        mu = bg_contribution + a_sig * Ts
        val = _poisson_val_grad_mu(mu, n_bin)[0]
        if abs(prev_val - val) < ALTERNATION_TOL * max(abs(val), 1.0):
            break
        prev_val = val
    return val, 10.0 ** theta_bg / scale_bg, a_sig / scale_sig, r


def profile_signal_uncertainty(T_bg, n_bin, T_sig, sig_best, bg_best,
                               delta: float = 1.0, max_iter: int = 40):
    """1-sigma interval on the signal amplitude via Delta(-2lnL) = `delta`,
    profiling (re-fitting) the background at each trial signal value, HELD
    FIXED at that value -- not re-optimized.

    Matches A2020's own established convention for exactly this kind of
    profile (`scripts/a2020_stage2_fit.py::_dm_scan`: "the source normalisation
    is held at a series of values and every other free parameter is refitted
    at each one"). `delta=1.0` is the standard one-parameter 68% interval.

    `sig_best`/`bg_best` are in ORIGINAL (unscaled) `T_sig`/`T_bg` units, as
    returned by `fit_bin_with_signal` -- this function rescales internally the
    SAME way `fit_bin_with_signal` does, for the SAME reason. Returns
    `(lo, hi)` in ORIGINAL units too.
    """
    m = T_bg.shape[0]
    scale_bg = np.maximum(np.array([np.abs(T_bg[i]).sum() for i in range(m)]), 1e-300)
    scale_sig = max(np.abs(T_sig).sum(), 1e-300)
    Tb = T_bg / scale_bg[:, None, None]
    Ts = T_sig / scale_sig
    theta0 = np.log10(np.clip(bg_best * scale_bg, 1e-300, None))
    a_best = sig_best * scale_sig

    def neg2lnl_at(a_sig):
        theta, _ = _fit_bg_given_offset(Tb, n_bin, a_sig * Ts, theta0)
        mu = _mu_bg(theta, Tb) + a_sig * Ts
        return _poisson_val_grad_mu(mu, n_bin)[0]

    val0 = neg2lnl_at(a_best)

    def _bracket_and_bisect(direction):
        step = max(abs(a_best), 1e-3) * 0.5
        a = a_best
        for _ in range(max_iter):
            a_next = a + direction * step
            if neg2lnl_at(a_next) - val0 >= delta:
                break
            a = a_next
            step *= 1.5
        else:
            # Never crossed delta: the profile is FLAT in this direction as
            # far as the bracket reached. Returning the last trial point would
            # fabricate a bound set by the step schedule, not the likelihood
            # (measured: identical fake "bounds" to 4 sig figs across
            # different configurations and injections). Report unbounded.
            return float("-inf") if direction < 0 else float("inf")
        lo, hi = (a, a_next) if direction > 0 else (a_next, a)
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            if (neg2lnl_at(mid) - val0) < delta:
                lo, hi = (mid, hi) if direction > 0 else (lo, mid)
            else:
                lo, hi = (lo, mid) if direction > 0 else (mid, hi)
        return 0.5 * (lo + hi) / scale_sig

    return _bracket_and_bisect(-1), _bracket_and_bisect(+1)


def fit_all_bins_with_signal(background_stack, signal_template, halo_mask=None,
                             component_order=None, with_uncertainty: bool = True) -> dict:
    """The halo search: background stack + one signal template, fit jointly
    per bin. `signal_template` is `(n_bins, 12, 12)` cell-space, matching the
    background templates' shape (see `stack.build_background_stack` and
    `stack.halo_cell_template`).
    """
    from .cells import halo_cell_mask

    names = component_order or sorted(background_stack["templates"])
    T = np.stack([background_stack["templates"][n] for n in names])
    counts = background_stack["counts"]
    n_bins = background_stack["n_bins"]
    halo = halo_cell_mask() if halo_mask is None else halo_mask

    per_bin = []
    for k in range(n_bins):
        val, bg_amp, sig_amp, r = fit_bin_with_signal(T[:, k], counts[k], signal_template[k])
        mu = expected_counts(np.log10(np.clip(bg_amp, 1e-300, None)), T[:, k]) \
            + sig_amp * signal_template[k]
        sig_resid = (counts[k] - mu) / np.sqrt(np.maximum(mu, 1.0))
        row = {
            "bin": k, "neg2lnl": val, "converged": bool(r.status == 0),
            "n_iter": int(r.nit), "bg_amplitudes": dict(zip(names, bg_amp)),
            "signal_amplitude": float(sig_amp),
            "observed_total": float(counts[k].sum()),
            "predicted_total": float(mu.sum()),
            "mean_abs_sigma_halo": float(np.abs(sig_resid[halo]).mean()),
            "max_abs_sigma_halo": float(np.abs(sig_resid[halo]).max()),
        }
        if with_uncertainty and r.status == 0:
            lo, hi = profile_signal_uncertainty(T[:, k], counts[k], signal_template[k],
                                                sig_amp, bg_amp)
            row["signal_amplitude_lo"] = float(lo)
            row["signal_amplitude_hi"] = float(hi)
        per_bin.append(row)
    return {"component_order": names, "per_bin": per_bin}


def restrict_stack(stack, keep, signal_template=None):
    """A copy of `stack` (and optionally the signal template) with every cell
    OUTSIDE the boolean `keep` mask excluded from the fit.

    Exclusion is implemented by zeroing both counts AND every template in the
    dropped cells: with n = 0 and mu = 0 there, a cell contributes nothing to
    the Poisson sum or any gradient (`xlogy(0, EPS) = 0`), which is exactly
    what "excluded from the analysis" means in this cell-space likelihood.

    🚨 THE BASELINE FIT NEEDED THIS AND DID NOT HAVE IT. Totani's halo-search
    fits EXCLUDE the disk (sec 3.2: "From here, the fits will be performed
    excluding the Galactic disk region") -- only his bubble-template
    construction includes it. The first-pass and v2 fits here fit ALL 144
    cells, disk included: 65% of the counts and the worst residuals (35-40
    sigma) sat in cells his analysis never uses, steering every background
    amplitude. The disk-excluded fit is not a variant; it is the procedure.
    His quadrant systematic (sec 3.4.1) is the same mechanism with four
    smaller masks.
    """
    keep = np.asarray(keep, bool)[None, :, :]
    out = {"counts": stack["counts"] * keep,
          "templates": {n: t * keep for n, t in stack["templates"].items()},
          "n_bins": stack["n_bins"]}
    if signal_template is None:
        return out
    return out, signal_template * keep


def quadrant_masks(b_min: float = 10.0):
    """Totani sec 3.4.1's four Galactic-coordinate quadrants, as cell masks
    (each already excluding the disk). Keys name the (l, b) signs."""
    from .cells import cell_centers, halo_cell_mask

    Lc, Bc = cell_centers()
    halo = halo_cell_mask(b_min=b_min)
    return {"l+b+": halo & (Lc > 0) & (Bc > 0),
            "l-b+": halo & (Lc < 0) & (Bc > 0),
            "l+b-": halo & (Lc > 0) & (Bc < 0),
            "l-b-": halo & (Lc < 0) & (Bc < 0)}


def fit_bin_profile(T_bg, n_bin, T_sig, n_scan: int = 25, n_refine: int = 3):
    """Signal amplitude by 1-D PROFILE over a_sig -- the estimator that does
    NOT alternate. Returns `(neg2lnl, bg_amplitudes, signal_amp, curve)`.

    🚨 REPLACES `fit_bin_with_signal` FOR ANY REPORTED NUMBER. Block
    alternation is correct in principle and stalls in practice: the halo
    template is only ~26% orthogonal to the background span, and in that soft
    valley each half-step moves the other block most of the way back, so the
    iteration converges by its own tolerance long before reaching the optimum.
    Measured 2026-09-04, injected-signal recovery SLOPE (a working estimator
    gives ~1 regardless of collinearity): 0.35-0.57 disk-excluded, 0.67-0.98
    disk-included. That bias is invisible in a real-data fit -- it just returns
    a plausible smaller number -- which is why the acceptance gate for this
    function is an INJECTION TEST, not a convergence flag.

    ⚠ THE REPORTED AMPLITUDE IS QUANTISED, and the lattice is the SAME for
    every model. `a_full` depends only on the bin's total counts and the
    signal template's sum -- both identical across rows -- and the nested
    `linspace` refinement lands `a_best` on a grid of spacing
    `a_full / 768` (verified exactly 2026-09-21: 43,195 counts in the
    20.8 GeV bin / 768 = 56.243490, the observed spacing among 257 scored
    rows, which take only 34 distinct values). That is 0.9% of a typical
    6,200-count halo, so it is a rounding rather than a bias and the TS error
    is second order -- but DIFFERENCES BETWEEN ROWS SMALLER THAN ~56 COUNTS
    ARE NOT REAL, and medians over rows can come out bit-identical between
    genuinely different runs for this reason. Raise `n_refine` if a finer
    point estimate is ever needed: each extra round divides the spacing by 4
    at the cost of 9 convex solves per bin.

    The profile has no such failure mode: `a_sig` is held FIXED on a scan grid
    and ONLY the background is refit at each point (`_fit_bg_given_offset`, a
    single well-conditioned convex solve). The reported MLE is the scan
    minimum, refined by golden-section bracketing. This is also exactly the
    construction `profile_signal_uncertainty` already uses for the interval, so
    the point estimate and its error bars now come from ONE curve rather than
    two different procedures that need not agree.
    """
    m = T_bg.shape[0]
    scale_bg = np.maximum(np.array([np.abs(T_bg[i]).sum() for i in range(m)]), 1e-300)
    scale_sig = max(np.abs(T_sig).sum(), 1e-300)
    Tb = T_bg / scale_bg[:, None, None]
    Ts = T_sig / scale_sig
    theta_seed = np.full(m, np.log10(max(n_bin.sum() / (m + 1), 1e-300)))

    cache = {}

    def profile(a):
        key = float(a)
        if key not in cache:
            theta, _ = _fit_bg_given_offset(Tb, n_bin, a * Ts, theta_seed)
            mu = _mu_bg(theta, Tb) + a * Ts
            cache[key] = (_poisson_val_grad_mu(mu, n_bin)[0], theta)
        return cache[key][0]

    # Scan symmetrically about zero on a scale set by the DATA: the amplitude
    # that would make the signal alone explain the whole bin is the natural
    # upper end, and the profile is smooth on that scale.
    a_full = n_bin.sum() / max(Ts.sum(), 1e-300)
    grid = np.linspace(-0.5 * a_full, 1.5 * a_full, n_scan)
    vals = np.array([profile(a) for a in grid])
    i = int(np.argmin(vals))

    lo = grid[max(i - 1, 0)]
    hi = grid[min(i + 1, len(grid) - 1)]
    for _ in range(n_refine):
        sub = np.linspace(lo, hi, 9)
        sv = np.array([profile(a) for a in sub])
        j = int(np.argmin(sv))
        lo, hi = sub[max(j - 1, 0)], sub[min(j + 1, len(sub) - 1)]
    a_best = 0.5 * (lo + hi)

    val = profile(a_best)
    theta = cache[float(a_best)][1]
    curve = {"grid": (grid / scale_sig).tolist(), "neg2lnl": vals.tolist()}
    return val, 10.0 ** theta / scale_bg, a_best / scale_sig, curve


def fit_all_bins_profile(background_stack, signal_template, halo_mask=None,
                         component_order=None, with_uncertainty: bool = True,
                         n_refine: int = 3) -> dict:
    """The halo search using the profile estimator (`fit_bin_profile`).

    Same output shape as `fit_all_bins_with_signal`, but every signal
    amplitude comes from the 1-D profile rather than block alternation. Use
    THIS for anything reported; see `fit_bin_profile` for why.
    """
    from .cells import halo_cell_mask

    names = component_order or sorted(background_stack["templates"])
    T = np.stack([background_stack["templates"][n] for n in names])
    counts = background_stack["counts"]
    n_bins = background_stack["n_bins"]
    halo = halo_cell_mask() if halo_mask is None else halo_mask

    per_bin = []
    for k in range(n_bins):
        val, bg_amp, sig_amp, curve = fit_bin_profile(T[:, k], counts[k],
                                                      signal_template[k], n_refine=n_refine)
        mu = expected_counts(np.log10(np.clip(bg_amp, 1e-300, None)), T[:, k]) \
            + sig_amp * signal_template[k]
        sig_resid = (counts[k] - mu) / np.sqrt(np.maximum(mu, 1.0))
        row = {"bin": k, "neg2lnl": val, "bg_amplitudes": dict(zip(names, bg_amp)),
              "signal_amplitude": float(sig_amp),
              "observed_total": float(counts[k].sum()),
              "predicted_total": float(mu.sum()),
              "mean_abs_sigma_halo": float(np.abs(sig_resid[halo]).mean()),
              "max_abs_sigma_halo": float(np.abs(sig_resid[halo]).max())}
        if with_uncertainty:
            lo, hi = profile_signal_uncertainty(T[:, k], counts[k],
                                                signal_template[k], sig_amp, bg_amp)
            row["signal_amplitude_lo"] = float(lo)
            row["signal_amplitude_hi"] = float(hi)
        per_bin.append(row)
    return {"component_order": names, "per_bin": per_bin}


def injection_slope(background_stack, signal_template, injections,
                    halo_mask=None, component_order=None, fitter=None):
    """Recovered-vs-injected SLOPE per bin -- the acceptance gate for any
    signal estimator used here.

    A correct linear-model MLE recovers an injected component at slope ~1
    however collinear the templates are: collinearity inflates the VARIANCE,
    not the bias. A slope materially below 1 therefore measures the OPTIMIZER,
    not the physics, and that is precisely the failure block alternation hid
    behind a converged flag.
    """
    fit = fitter or fit_all_bins_profile
    names = component_order or sorted(background_stack["templates"])
    recovered = []
    for a in injections:
        st = {"counts": background_stack["counts"] + a * signal_template,
             "templates": background_stack["templates"],
             "n_bins": background_stack["n_bins"]}
        r = fit(st, signal_template, halo_mask=halo_mask,
               component_order=names, with_uncertainty=False)
        recovered.append([row["signal_amplitude"] for row in r["per_bin"]])
    rec = np.array(recovered)
    inj = np.asarray(injections, float)
    # 🚨 COVARIANCE FORM, NOT `lstsq` ON [inj, 1]. The amplitudes here are of
    # order 1e-33, so a design matrix pairing that column with a column of
    # ONES has a singular-value ratio of ~1e-33 -- below `lstsq`'s rcond
    # cutoff, which then TRUNCATES the slope direction and returns exactly
    # 0.0. Measured 2026-09-04: this gate reported slope = 0.000 for BOTH
    # estimators, including one already known from a hand check to be 0.35 --
    # i.e. the gate would have condemned a working estimator and could equally
    # have passed a broken one. The covariance ratio is scale-free.
    dinj = inj - inj.mean()
    denom = float(np.sum(dinj ** 2))
    if denom <= 0:
        raise ValueError("injections must span more than one value")
    slopes = [float(np.sum(dinj * (rec[:, k] - rec[:, k].mean())) / denom)
             for k in range(rec.shape[1])]
    return {"injections": list(inj), "recovered": rec.tolist(), "slopes": slopes}
