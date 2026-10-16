"""PSF convolution and exposure application — the rest of the template forward model.

Task 7, downstream of `regrid.py`. Jax-free.

⚠ **The King parameters here are Track A's EFFECTIVE kernel, not the LAT PSF.**
`gcefisher.cholis_ingest` fitted them to map 0.25° Cholis release maps onto
gcepy's 0.1° counts, so they absorb the PSF *and* that release's 0.25° pixel
boxcar *and* the in-bin energy integration *and* the inclination average. The
Pohl maps are 0.5°, so the boxcar term is wrong for them. Reused here because:

- for **cost measurement** the kernel's provenance is irrelevant — a King
  convolution of a given width costs what it costs; and
- for the **exposure-free ratio validation** (see `validate_ratio`) the kernel
  is common to numerator and denominator and largely divides out.

For anything where the absolute kernel matters, use `gtsrcmaps` with the real
IRF (Task 7 proper). Recorded so nobody promotes this to "the PSF".
"""

from __future__ import annotations

import numpy as np

__all__ = ["king_kernel", "convolve_map", "apply_forward_model"]


def king_kernel(sigma_deg: float, gamma: float, pixel_deg: float,
                n: int | None = None) -> np.ndarray:
    """Normalized King profile on a square grid, matching Track A's form.

    k(r) ∝ (1 + r² / (2 γ σ²))^(−γ)
    """
    if n is None:                        # cover ~12 sigma, odd size, capped
        n = int(min(401, max(15, 2 * int(12 * sigma_deg / pixel_deg) + 1)))
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    r2 = ((yy - c) ** 2 + (xx - c) ** 2) * pixel_deg ** 2
    k = (1.0 + r2 / (2.0 * gamma * sigma_deg ** 2)) ** (-gamma)
    return k / k.sum()


def convolve_map(a: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """FFT convolution with edge padding.

    ⚠ Padded by the kernel half-width rather than wrapped. A periodic
    convolution would wrap flux from one edge of the ROI to the other, which on
    a 40°×40° Galactic map means the northern edge leaking into the southern —
    a small effect, but a silent and entirely avoidable one.
    """
    from numpy.fft import irfft2, rfft2

    pad = kernel.shape[0] // 2
    ap = np.pad(a, pad, mode="edge")
    sh = ap.shape
    kp = np.zeros(sh)
    kc = kernel.shape[0]
    kp[:kc, :kc] = kernel
    kp = np.roll(kp, (-(kc // 2), -(kc // 2)), axis=(0, 1))
    out = np.real(irfft2(rfft2(ap) * rfft2(kp), s=sh))
    return out[pad:pad + a.shape[0], pad:pad + a.shape[1]]


def apply_forward_model(template, sigma_deg, gamma, pixel_deg, exposure=None):
    """template → PSF-convolved (× exposure if given).

    `template` is 2-D (pure morphology, one map for all bins) or 3-D with a
    leading energy axis. `sigma_deg`/`gamma` are scalars or per-bin arrays.
    """
    t = np.asarray(template, float)
    sig = np.atleast_1d(sigma_deg)
    gam = np.atleast_1d(gamma)
    nb = len(sig)

    if t.ndim == 2:
        planes = [convolve_map(t, king_kernel(sig[b], gam[b], pixel_deg))
                  for b in range(nb)]
    else:
        if t.shape[0] != nb:
            raise ValueError(f"template has {t.shape[0]} planes, {nb} kernels given")
        planes = [convolve_map(t[b], king_kernel(sig[b], gam[b], pixel_deg))
                  for b in range(nb)]
    out = np.stack(planes)
    if exposure is not None:
        e = np.asarray(exposure, float)
        out = out * (e if e.ndim == 3 else e[None])
    return out


def validate_ratio(ours_a, ours_b, ref_a, ref_b, mask=None):
    """Compare a template RATIO against a reference — exposure-free by construction.

    Counts are ``T ⊗ PSF × exposure``, so for two templates in the same energy
    bin the exposure cancels exactly in the ratio and the PSF very largely does.
    That makes this a validation we can run against Song's gcepy-format arrays
    **without knowing their exposure**, which is what made the earlier direct
    comparison inconclusive.

    Returns (median ratio-of-ratios, robust scatter).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        r_ours = np.asarray(ours_a, float) / np.asarray(ours_b, float)
        r_ref = np.asarray(ref_a, float) / np.asarray(ref_b, float)
        q = r_ours / r_ref
    good = np.isfinite(q) & (q > 0)
    if mask is not None:
        good &= mask
    if not good.any():
        return np.nan, np.nan
    v = q[good]
    return float(np.median(v)), float(np.subtract(*np.percentile(v, [84, 16])) / 2)
