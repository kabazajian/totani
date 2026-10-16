"""DM annihilation templates and J-factors, following Abazajian et al. 2020.

Task 4 / Task 7. gcepy ships **one** fixed DM morphology (gNFW γ = 1.2), so
varying γ — or using a cored profile at all — requires generating templates
ourselves. This module does that.

**Array-module agnostic.** Every function takes `xp`, defaulting to `numpy`.
Pass `jax.numpy` and the whole chain is differentiable, which is what Task 4
needs to put γ, r_c, ρ⊙ and the axis ratios into the Fisher matrix. Nothing
here imports jax, so it remains safe to use from the `fermi` environment too.

PROFILES (arXiv:2003.10416)
---------------------------
**gNFW**, their Eq. 2 — normalized so ρ(R⊙) = ρ⊙ exactly, for any (γ, r_s):

    ρ(r) = ρ⊙ (R⊙/r)^γ [(r_s + R⊙)/(r_s + r)]^(3−γ)

**Cored ("Read" et al. 2016)**, their Eq. 3 — defined through the *enclosed
mass*, not the density:

    M_c(r) = M_NFW(r) · tanh(r/r_c),     M_NFW is the γ = 1 NFW mass

so the density follows from ρ_c = (1/4πr²) dM_c/dr, giving

    ρ_c(r) = ρ_NFW(r)·tanh(r/r_c) + M_NFW(r)·sech²(r/r_c) / (4π r² r_c)

Both terms are finite as r → 0 (the first because ρ_NFW ∝ 1/r cancels the
tanh's linear rise), which is what makes the profile cored. A2020 fixes
r_c = 1 kpc.

**Ellipsoidal**, their Eq. B1 — substitute r → r′ with

    r′² = x² + y²/(b/a)² + z²/(c/a)²

A2020 uses c/a = 0.7 (oblate, short axis perpendicular to the disc), b/a = 1.

CONVENTIONS
-----------
Lengths in kpc, densities in GeV cm⁻³. ``j_map`` returns ∫ρ² ds in
**GeV² cm⁻⁵ sr⁻¹**; ``j_factor_roi`` returns ∫dΩ∫ds ρ² in **GeV² cm⁻⁵**, which
is A2020's convention (note it carries the solid-angle integral — Stenhouse et
al. use a *per-steradian* pole value instead, and the two must never be
compared without conversion; see `prompts/task10_positioning.md` §3).
"""

from __future__ import annotations

import numpy as _np

__all__ = ["rho_gnfw", "mass_nfw", "rho_read_cored", "j_map", "dm_template",
           "j_factor_roi", "KPC_CM", "A2020_DEFAULTS", "A2020_CANONICAL",
           "A2020_MARGINALISATION_CENTRES",
           "GCEPY_EXCESS_KING_SIGMA", "GCEPY_EXCESS_KING_GAMMA",
           "gcepy_convention_template"]

KPC_CM = 3.0856775814913673e21          # cm per kpc

#: Effective PSF of **gcepy's own excess templates**, fitted 2026-07-30 directly
#: from its shipped gNFW γ=1.2 `dm` template. Exposure was divided out exactly by
#: forming ``dm / isotropic`` — the isotropic template is a uniform sky times the
#: same exposure and PSF, so the ratio is the PSF-convolved DM morphology alone.
#: Agreement of our generator ⊗ this kernel with gcepy's template is ≤ 10.2%
#: over 0.15–15° in every bin.
#:
#: ⚠ **RE-FITTED 2026-07-30 after a bug fix in the line-of-sight integration.**
#: The first fit reached 8.4% — against a generator whose ∫ρ²ds was numerically
#: wrong at small θ. The two-parameter kernel had partly **absorbed** that error,
#: which is worth remembering: fitting a nuisance function to data can hide a
#: bug in the model and still look like a good fit. With the integration fixed
#: the kernel narrows (σ₀ 0.5135 → 0.5060, σ₁₃ 0.0833 → 0.0750) and the residual
#: is slightly *larger*, 10.2%. That residual is the King form's inability to
#: match the true effective kernel, not a model error — which is exactly why
#: this is called an effective kernel and not "the PSF".
#:
#: ⚠ These are NOT `gcefisher.cholis_ingest.KING_SIGMA`. That kernel additionally
#: absorbs the Cholis release's 0.25° pixel boxcar, which gcepy's excess templates
#: do not carry — hence it is ~21% too wide in bin 0. Use THESE when generating an
#: excess template meant to sit alongside gcepy's own.
GCEPY_EXCESS_KING_SIGMA = _np.array(
    [0.5060, 0.3903, 0.2917, 0.2167, 0.1917, 0.1500, 0.1256,
     0.1003, 0.0917, 0.0833, 0.0667, 0.0667, 0.0750, 0.0750])
GCEPY_EXCESS_KING_GAMMA = _np.array(
    [1.725, 1.657, 1.657, 1.536, 1.661, 1.622, 1.673,
     1.648, 1.748, 1.789, 1.633, 1.655, 2.138, 2.315])

#: Canonical MW values (ρ⊙ = 0.30, R⊙ = 8.25, r_s = 20), offered by Kev as the
#: set possibly behind A2020's published J-factors.
#:
#: 🚨 **That hypothesis is retired, and so is the "0.3% agreement" once claimed
#: here.** The agreement was a COINCIDENCE — it compared our *oblate* (q_c = 0.7)
#: J against A2020's *spherical* `nomask` file, and two factors cancelled:
#: ours_oblate/ours_spherical = 0.743 × ours_spherical/theirs = 1.342 = 0.997.
#: Like for like, this set gives **1.29–1.34× their mean**.
#:
#: The deeper point is that the question was ill-posed: A2020's files are
#: **marginalised distributions** (their own loader: "the prior on the J-factor,
#: calculated from a convolution of priors on rho_local, gamma, and the scale
#: radius"), so no fixed parameter set underlies them at all. Kept for reference
#: as standard canonical values; NOT an A2020 anchor.
A2020_CANONICAL = dict(rho_sun=0.30, R_sun=8.25, rs=20.0)

#: 🚨 **These are the ones to use against A2020.** The centres of the priors they
#: marginalise over — ρ⊙ = 0.28 ± 0.08 (Zhang et al.), r_s log-normal about
#: 26 kpc with width 0.14 dex, quoted in their §III. Marginalising over exactly
#: these reproduces their published distributions; see `A2020_J_STATUS`.
#:
#: ⚠ R⊙ = 8.18 here is *our* choice, not A2020's — the paper never states one.
#: It was scanned over 7.9–8.6 kpc and does not resolve the residual either way.
A2020_MARGINALISATION_CENTRES = dict(rho_sun=0.28, R_sun=8.18, rs=26.0)

#: ✅ RESOLVED 2026-07-30, once the spherical default was fixed. A2020's files
#: are named `J_MC_*` — they are **Monte-Carlo distributions** obtained by
#: marginalising over the §III priors, so their median sits at the prior
#: CENTRES, not at any separate "canonical" set.
#:
#: Reproducing that marginalisation ourselves (ρ⊙ ~ N(0.28, 0.08), r_s
#: log-normal about 26 kpc at 0.14 dex, spherical) matches their **distribution
#: shape**: their mean/median is 1.086 at every γ, and our median-ratio and
#: mean-ratio agree with each other to <0.5% everywhere. So the §III priors are
#: right and only a multiplicative piece is left.
#:
#: 🚨 That piece is **not** R⊙, and it is **not** uniform. Scanned R⊙ over
#: 7.9–8.6 kpc (`scripts/a2020_j_normalisation.py`): it rescales every γ almost
#: identically, so it cannot remove what is actually a **γ-dependent** residual —
#: ratio to their median 1.000 / 1.007 / 1.035 / 1.127 at γ = 1.0 / 1.1 / 1.2 /
#: 1.3. Our own integral is excluded as the cause (converged to 0.2% at all
#: three γ over θ_min = 1e-4…1e-8, n_θ = 300…1200).
#:
#: Monotonic, accelerating toward the γ = 1.5 divergence — the signature of an
#: under-resolved ROI quadrature on *their* side, the same bug we found in
#: ourselves (`j_factor_bugs.bug_2`). Hypothesis only: their J generator is not
#: in the public release. Question for Oscar and Kev.
#:
#: The "canonical fixed values" hypothesis is retired: their files are
#: marginalised distributions, so no fixed parameter set underlies them.
A2020_J_STATUS = ("gamma<=1.2 agrees to 5-10%; residual TRENDS with gamma "
                  "(19% at 1.3) and is not R_sun -- do not quote J at gamma>=1.3")

#: ⚠ STILL OPEN: our ELLIPSOIDAL J falls (0.743×) where A2020's published
#: "triaxial" J rises (1.467×).
#:
#: 🚨 An earlier note here claimed compression "must REDUCE" J. **That is wrong**
#: (Kev, 2026-07-30): J is a density-SQUARED integral, so compressing a halo at
#: fixed MASS concentrates the density and RAISES J. The sign depends entirely
#: on what is held fixed.
#:
#: What we implement is the bare substitution r → r′ at fixed ρ(R⊙). Since the
#: Sun sits in the plane, r′ = R⊙ there and the normalisation is untouched,
#: while off-plane r′ > r lowers the density — so mass is REMOVED rather than
#: redistributed, and J falls. A mass-conserving flattening, or a compression
#: along a different axis (e.g. along the line of sight rather than the disc
#: normal), would raise it instead.
#:
#: Which of those A2020 did is not yet determined, and our implementation may
#: simply not match their intent. Do not quote an absolute ellipsoidal J.

#: ⚠ **q_c = 1.0 (SPHERICAL) is the default, corrected 2026-07-30.** A2020's
#: BASELINE profile is the spherical gNFW (their Fig. 3c); the c/a = 0.7
#: ellipsoid (Fig. 3d) is one of several morphology VARIANTS they test, not the
#: fiducial. An earlier version made 0.7 the default, which silently entered
#: every bare call and produced a 34% error in J that masqueraded as a 0.3%
#: agreement. Pass `q_c=0.7` explicitly to get the ellipsoidal variant.
A2020_DEFAULTS = dict(gamma=1.2, rc=1.0, q_b=1.0, q_c=1.0,
                      **A2020_MARGINALISATION_CENTRES)


def _rho_s_rs3(rs, rho_sun, R_sun):
    """ρ_s r_s³ for the γ=1 NFW matching ρ(R⊙) = ρ⊙, used by the mass profile."""
    return rho_sun * R_sun * (rs + R_sun) ** 2


def rho_gnfw(r, gamma=1.2, rs=26.0, rho_sun=0.28, R_sun=8.18, xp=_np):
    """Generalized NFW, A2020 Eq. 2. ρ(R⊙) = ρ⊙ exactly for any (γ, r_s)."""
    r = xp.maximum(xp.asarray(r, dtype=float), 1e-8)
    return rho_sun * (R_sun / r) ** gamma * ((rs + R_sun) / (rs + r)) ** (3.0 - gamma)


def mass_nfw(r, rs=26.0, rho_sun=0.28, R_sun=8.18, xp=_np):
    """Enclosed mass of the γ = 1 NFW, in GeV cm⁻³ kpc³."""
    r = xp.maximum(xp.asarray(r, dtype=float), 1e-8)
    x = r / rs
    # M(r) = 4π ρ_s r_s³ [ln(1+x) − x/(1+x)], and matching ρ(R⊙)=ρ⊙ for γ=1
    # gives ρ_s r_s³ = ρ⊙ R⊙ (r_s + R⊙)². No further r_s factor — an earlier
    # version divided by r_s² here, which the M'(r) = 4πr²ρ test caught.
    return 4.0 * xp.pi * _rho_s_rs3(rs, rho_sun, R_sun) * (
        xp.log1p(x) - x / (1.0 + x))


def rho_read_cored(r, rc=1.0, rs=26.0, rho_sun=0.28, R_sun=8.18, xp=_np):
    """Read et al. cored profile as used by A2020 (their Eq. 3).

    Derived from M_c = M_NFW·tanh(r/r_c) rather than written down directly, so
    that the mass really is the tanh-suppressed NFW mass — the defining
    property. Both terms stay finite at r → 0.
    """
    r = xp.maximum(xp.asarray(r, dtype=float), 1e-8)
    t = xp.tanh(r / rc)
    rho_n = rho_gnfw(r, 1.0, rs, rho_sun, R_sun, xp=xp)
    m_n = mass_nfw(r, rs, rho_sun, R_sun, xp=xp)
    sech2 = 1.0 - t ** 2
    return rho_n * t + m_n * sech2 / (4.0 * xp.pi * r ** 2 * rc)


def _los_nodes(grid_l, grid_b, R_sun, s_max, n_s, q_b, q_c, xp=_np):
    """Line-of-sight sample points, concentrated at the point of closest approach.

    🚨 A UNIFORM grid in s is unusable for a cusped profile. The integrand of
    ∫ρ² ds peaks over a width ~p about the tangent point, where p = R⊙ sin ψ is
    the impact parameter — and p falls below the step size ds for θ ≲ 1°, so the
    peak lands *between* samples and the integral becomes meaningless exactly
    where a cusped J-factor gets most of its value. Measured before the fix:
    J over a 4°×4° ROI at γ=1.4 changed by ×8.9 between 0.4° and 0.02° pixels,
    i.e. it never converged, while the cored profile (no cusp, no sharp peak)
    converged to 0.06%.

    The fix is the substitution ``s = s₀ + p·sinh(u)``, for which a *spherical*
    profile gives ``r = p·cosh(u)`` exactly — the integrand becomes smooth in u
    and a uniform u grid resolves it at every θ. For the ellipsoidal case the
    substitution is no longer exact, but it still concentrates samples where the
    integrand peaks, and convergence is verified numerically either way.

    Returns ``(r, jacobian, t)`` with shapes broadcasting to (nb, nl, n_s).
    """
    L = xp.deg2rad(xp.asarray(grid_l, dtype=float))[None, :]
    B = xp.deg2rad(xp.asarray(grid_b, dtype=float))[:, None]

    cos_psi = xp.cos(B) * xp.cos(L)
    s0 = R_sun * cos_psi                                     # tangent point
    # impact parameter, floored so a sight line straight at the GC stays finite
    p_imp = R_sun * xp.sqrt(xp.maximum(1.0 - cos_psi ** 2, 1e-12))

    u_lo = xp.arcsinh((0.0 - s0) / p_imp)
    u_hi = xp.arcsinh((s_max - s0) / p_imp)
    t = xp.linspace(0.0, 1.0, n_s)
    u = u_lo[..., None] + (u_hi - u_lo)[..., None] * t
    S = s0[..., None] + p_imp[..., None] * xp.sinh(u)
    jac = p_imp[..., None] * xp.cosh(u) * (u_hi - u_lo)[..., None]

    Lc, Bc = L[..., None], B[..., None]
    x = R_sun - S * xp.cos(Bc) * xp.cos(Lc)
    y = -S * xp.cos(Bc) * xp.sin(Lc)
    z = S * xp.sin(Bc)
    r = xp.sqrt(x ** 2 + (y / q_b) ** 2 + (z / q_c) ** 2)
    return r, jac, t


def j_map(grid_l, grid_b, profile="gnfw", s_max=100.0, n_s=800, xp=_np, **kw):
    """∫ρ² ds along each sight line, in GeV² cm⁻⁵ sr⁻¹.

    `profile` is "gnfw" or "cored"; remaining keywords go to the profile
    (γ, r_s, ρ⊙, R⊙, r_c) or to the geometry (q_b, q_c).
    """
    p = dict(A2020_DEFAULTS)
    # 🚨 REJECT UNKNOWN KEYS. `p.update(kw)` alone ADDS them, so a mis-cased or
    # mistyped profile parameter is silently discarded and the DEFAULT profile is
    # returned with no error. Caught 2026-08-09: a scan over `r_sun` (the field
    # is `R_sun`) and `rs` produced four BYTE-IDENTICAL maps, which read as a
    # physics result and was the parameters never arriving. A caller that
    # believes it varied the halo and did not is exactly the silent-wrong-number
    # failure rule 1 exists for.
    unknown = set(kw) - set(p) - {"profile"}
    if unknown:
        raise TypeError(
            f"unknown profile parameter(s) {sorted(unknown)}; expected any of "
            f"{sorted(p)}. Refusing to silently return the default profile.")
    p.update(kw)
    r, jac, t = _los_nodes(grid_l, grid_b, p["R_sun"], s_max, n_s,
                           p["q_b"], p["q_c"], xp=xp)

    if profile == "gnfw":
        rho = rho_gnfw(r, p["gamma"], p["rs"], p["rho_sun"], p["R_sun"], xp=xp)
    elif profile == "cored":
        rho = rho_read_cored(r, p["rc"], p["rs"], p["rho_sun"], p["R_sun"], xp=xp)
    else:
        raise ValueError(f"unknown profile {profile!r}; use 'gnfw' or 'cored'")

    return xp.trapezoid(rho ** 2 * jac, t, axis=-1) * KPC_CM


def dm_template(grid_l, grid_b, normalize=True, xp=_np, **kw):
    """J map on a grid, optionally normalized to unit sum.

    A2020's Fig. 3 maps are "normalized to unit flux"; every template in this
    project carries a free per-bin normalization anyway, so the scale is
    convention only — but the normalization must be applied consistently, and
    `j_factor_roi` is the right entry point when the absolute J is wanted.
    """
    j = j_map(grid_l, grid_b, xp=xp, **kw)
    return j / xp.sum(j) if normalize else j


def j_factor_roi(half_width_deg=20.0, pixel_deg=None, n_theta=600, n_phi=360,
                 theta_min_deg=1e-5, xp=_np, **kw):
    """J = ∫dΩ ∫ds ρ² over a square ROI, in GeV² cm⁻⁵ — A2020's convention.

    🚨 Integrated in POLAR coordinates with a LOG-SPACED θ grid, not on the
    Cartesian analysis grid. For a cusped profile J(θ) ∝ θ^(1−2γ), so the
    solid-angle integrand goes as θ^(2−2γ) and a midpoint rule on square pixels
    converges as δ^(2−2γ) — hopelessly slow for steep γ. Measured on a 4°×4° ROI
    with the Cartesian scheme: γ=1.4 moved 62% between 0.4° and 0.02° pixels and
    was still climbing, while γ=1.0 moved 5%. The log-polar grid resolves the
    cusp directly and converges at every γ.

    Exact spherical geometry: a point at angular distance θ and position angle φ
    from the GC sits at b = arcsin(sin θ sin φ), ℓ = arctan2(sin θ cos φ, cos θ),
    with dΩ = sin θ dθ dφ. Points outside the square ROI are dropped.

    ``pixel_deg`` is accepted and ignored, so existing callers keep working; it
    no longer has any meaning here.

    ⚠ **Inner truncation matters more as γ steepens.** The solid-angle integrand
    goes as θ^(2−2γ), so the fraction lost below ``theta_min_deg`` scales as
    θ_min^(3−2γ) — exponent 1.0, 0.6, 0.2 at γ = 1.0, 1.2, 1.4. Measured: γ=1.0
    and 1.2 are converged by θ_min = 1e-4, but **γ = 1.4 is still climbing at
    1e-6** and should not be quoted without settling the inner treatment. (The
    integral formally diverges at γ ≥ 1.5.)

    ⚠ This is an ROI-integrated quantity. Stenhouse et al. (arXiv:2607.08552)
    instead quote a per-steradian line-of-sight integral toward the Galactic
    poles. Converting between them requires recomputing from a common profile;
    do not rescale one into the other.
    """
    th = xp.exp(xp.linspace(xp.log(theta_min_deg),
                            xp.log(half_width_deg * 2.0 ** 0.5), n_theta))
    ph = xp.linspace(0.0, 2.0 * xp.pi, n_phi, endpoint=False)

    T = xp.deg2rad(th)[:, None]
    P = ph[None, :]
    b = xp.arcsin(xp.sin(T) * xp.sin(P))
    l = xp.arctan2(xp.sin(T) * xp.cos(P), xp.cos(T))
    ld, bd = xp.rad2deg(l), xp.rad2deg(b)

    inside = (xp.abs(ld) <= half_width_deg) & (xp.abs(bd) <= half_width_deg)

    # j_map expects 1-D grids; evaluate pointwise on the flattened polar mesh
    p = dict(A2020_DEFAULTS)
    # 🚨 REJECT UNKNOWN KEYS. `p.update(kw)` alone ADDS them, so a mis-cased or
    # mistyped profile parameter is silently discarded and the DEFAULT profile is
    # returned with no error. Caught 2026-08-09: a scan over `r_sun` (the field
    # is `R_sun`) and `rs` produced four BYTE-IDENTICAL maps, which read as a
    # physics result and was the parameters never arriving. A caller that
    # believes it varied the halo and did not is exactly the silent-wrong-number
    # failure rule 1 exists for.
    unknown = set(kw) - set(p) - {"profile"}
    if unknown:
        raise TypeError(
            f"unknown profile parameter(s) {sorted(unknown)}; expected any of "
            f"{sorted(p)}. Refusing to silently return the default profile.")
    p.update(kw)
    profile = p.pop("profile", kw.get("profile", "gnfw"))
    r, jac, t = _los_nodes(ld.ravel(), xp.zeros(1), p["R_sun"], 100.0, 600,
                           p["q_b"], p["q_c"], xp=xp)  # placeholder, replaced below

    # build radii directly for the (theta, phi) points
    from_l = xp.deg2rad(ld).ravel()
    from_b = xp.deg2rad(bd).ravel()
    cos_psi = xp.cos(from_b) * xp.cos(from_l)
    s0 = p["R_sun"] * cos_psi
    p_imp = p["R_sun"] * xp.sqrt(xp.maximum(1.0 - cos_psi ** 2, 1e-12))
    u_lo = xp.arcsinh((0.0 - s0) / p_imp)
    u_hi = xp.arcsinh((100.0 - s0) / p_imp)
    tt = xp.linspace(0.0, 1.0, 600)
    u = u_lo[:, None] + (u_hi - u_lo)[:, None] * tt
    S = s0[:, None] + p_imp[:, None] * xp.sinh(u)
    JJ = p_imp[:, None] * xp.cosh(u) * (u_hi - u_lo)[:, None]
    x = p["R_sun"] - S * xp.cos(from_b)[:, None] * xp.cos(from_l)[:, None]
    y = -S * xp.cos(from_b)[:, None] * xp.sin(from_l)[:, None]
    z = S * xp.sin(from_b)[:, None]
    rr = xp.sqrt(x ** 2 + (y / p["q_b"]) ** 2 + (z / p["q_c"]) ** 2)

    if profile == "gnfw":
        rho = rho_gnfw(rr, p["gamma"], p["rs"], p["rho_sun"], p["R_sun"], xp=xp)
    elif profile == "cored":
        rho = rho_read_cored(rr, p["rc"], p["rs"], p["rho_sun"], p["R_sun"], xp=xp)
    else:
        raise ValueError(f"unknown profile {profile!r}; use 'gnfw' or 'cored'")

    Jvals = (xp.trapezoid(rho ** 2 * JJ, tt, axis=-1) * KPC_CM).reshape(ld.shape)
    Jvals = xp.where(inside, Jvals, 0.0)

    dphi = 2.0 * xp.pi / n_phi
    return xp.trapezoid(xp.sum(Jvals, axis=1) * dphi * xp.sin(xp.deg2rad(th)),
                        xp.deg2rad(th))

def gcepy_convention_template(half_width_deg=20.0, pixel_deg=0.1, normalize=True,
                              **kw):
    """A (14, N, N) excess template in **gcepy's own convention**.

    This is the Task 4 deliverable: gcepy ships a single gNFW γ=1.2 morphology
    and no way to vary it, so DM-profile marginalization needs templates
    generated here that are interchangeable with gcepy's. Applies
    `GCEPY_EXCESS_KING_*` per energy bin — the kernel fitted from gcepy's own
    `dm` template — so the output can sit alongside the shipped ones.

    ⚠ Returned in the **generator's** (FITS/WCS) longitude orientation. Apply
    ``[:, :, ::-1]`` before use with gcepy arrays, exactly as elsewhere in this
    project. Note that a *spherical* profile is symmetric in ℓ so the flip is a
    no-op there — but it is NOT for the ellipsoidal case (q_b ≠ 1), which is
    precisely when a silent mistake would matter.
    """
    from gcereduce.forward import convolve_map, king_kernel

    n = int(round(2 * half_width_deg / pixel_deg))
    c = (_np.arange(n) + 0.5 - n / 2.0) * pixel_deg
    raw = j_map(c, c, **kw)
    out = _np.stack([
        convolve_map(raw, king_kernel(GCEPY_EXCESS_KING_SIGMA[b],
                                      GCEPY_EXCESS_KING_GAMMA[b], pixel_deg))
        for b in range(len(GCEPY_EXCESS_KING_SIGMA))])
    if normalize:
        out = out / out.sum(axis=(1, 2), keepdims=True)
    return out
