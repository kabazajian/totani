"""The LAT point-spread function folded into the pixel templates (2026-09-27).

Two tables, chosen by FINE_TAG (read when the table is first used):

* FINE_TAG=psfuc (2026-09-30, the consistent one): gtpsf for P8R3_ULTRACLEAN_V3,
  evtype 3, the Galactic Center direction, with THIS run's livetime cube --
  the data's own class and exposure. `runs/fff2de5bc7e0/psf_P8R3_ULTRACLEAN_V3_evtype3.fits`
  (+ .json provenance).
* FINE_TAG=psf (2026-09-27 to 09-29 results): the gtpsf table of gce-fisher run
  1bbcd8ac5ddb (P8R3 CLEAN, evtype 3, that run's livetime cube), used as a
  proxy. `results/psf_class_check_2026-09-30.md` measures the difference:
  R68 within 4% at >= 5.6 GeV, the CLEAN kernel 3-9% wider at 1.5-4.3 GeV.

Per analysis bin the table's PSF(theta) is averaged over the bin's energies
with weight E^-2.4 x exposure (the table's own column), then represented by a
sum of three Gaussians fitted to its containment curve. A template's pixel map
(plate carree, 0.125 deg) is convolved separably per Gaussian: along latitude
with sigma/0.125 pixels, then along longitude row by row with
sigma/(0.125 cos b) pixels, which is exact on the sphere to first order for
the sub-degree to few-degree widths here. Edges are handled by reflection
('nearest'); the maps include the disk, so disk emission spills into the
|b| >= 10 cells as it does on the sky. `check_point_source` verifies the
containment radii and flux conservation of the implementation.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

PSF_TABLE_CLEAN = Path("/Users/aba/work/gce-fisher/runs/1bbcd8ac5ddb/psf.fits")
PSF_TABLE_ULTRACLEAN = Path(__file__).resolve().parents[2] / "runs" / "fff2de5bc7e0" / "psf_P8R3_ULTRACLEAN_V3_evtype3.fits"


def table_path() -> Path:
    """The PSF table for the current FINE_TAG: ULTRACLEAN for 'psfuc', CLEAN otherwise."""
    import os
    return PSF_TABLE_ULTRACLEAN if os.environ.get("FINE_TAG", "") == "psfuc" else PSF_TABLE_CLEAN


PSF_TABLE = PSF_TABLE_CLEAN          # kept for older imports; use table_path()
PIX = 0.125
N_GAUSS = 3


@lru_cache(maxsize=None)
def _table():
    from astropy.io import fits
    with fits.open(table_path()) as h:
        E = np.asarray(h["PSF"].data["Energy"], float)
        X = np.asarray(h["PSF"].data["Exposure"], float)
        P = np.asarray(h["PSF"].data["Psf"], float)
        th = np.asarray(h["THETA"].data["Theta"], float)
    return E, X, P, th


def bin_profile(e_lo_gev, e_hi_gev):
    """(theta deg, PSF sr^-1) averaged over the bin with weight E^-2.4 x exposure."""
    E, X, P, th = _table()
    Eg = np.logspace(np.log10(e_lo_gev * 1e3), np.log10(e_hi_gev * 1e3), 16)
    w = Eg ** -2.4 * np.interp(np.log(Eg), np.log(E), X) * Eg          # per d lnE
    logP = np.log(np.maximum(P, 1e-300))
    Pg = np.array([np.exp(np.array([np.interp(np.log(e), np.log(E), logP[:, j]) for j in range(len(th))])) for e in Eg])
    prof = (w[:, None] * Pg).sum(0) / w.sum()
    t = np.radians(th)
    prof /= np.trapezoid(prof * 2 * np.pi * np.sin(t), t)
    return th, prof


def containment(th, prof):
    t = np.radians(th)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (prof[1:] * np.sin(t[1:]) + prof[:-1] * np.sin(t[:-1])) * np.diff(t))]) * 2 * np.pi
    return c / c[-1]


@lru_cache(maxsize=None)
def gaussians(e_lo_gev, e_hi_gev):
    """((w1, w2, w3), (s1, s2, s3) deg): three Gaussians matching the bin's containment curve."""
    from scipy.optimize import least_squares
    th, prof = bin_profile(e_lo_gev, e_hi_gev)
    c = containment(th, prof)
    r68 = float(np.interp(0.68, c, th))
    x = th[th <= 12.0]; y = c[th <= 12.0]

    def model(p):
        w = np.exp(p[:N_GAUSS]); w = w / w.sum(); s = np.exp(p[N_GAUSS:])
        return (w[:, None] * (1 - np.exp(-np.radians(x)[None] ** 2 / (2 * np.radians(s)[:, None] ** 2)))).sum(0)

    p0 = np.concatenate([np.log([0.6, 0.3, 0.1]), np.log([0.5 * r68, 1.5 * r68, 5 * r68])])
    r = least_squares(lambda p: model(p) - y, p0)
    w = np.exp(r.x[:N_GAUSS]); w = w / w.sum(); s = np.exp(r.x[N_GAUSS:])
    o = np.argsort(s)
    return tuple(float(v) for v in w[o]), tuple(float(v) for v in s[o])


def convolve(img, B, e_lo_gev, e_hi_gev):
    """PSF-convolve one pixel map (ny, nx) on the plate-carree grid with latitudes B (deg)."""
    from scipy.ndimage import gaussian_filter1d
    w, s = gaussians(float(e_lo_gev), float(e_hi_gev))
    img = np.asarray(img, float)
    b_row = np.asarray(B, float)[:, 0]
    cosb = np.maximum(np.cos(np.radians(b_row)), 0.05)
    out = np.zeros_like(img)
    for wi, si in zip(w, s):
        tmp = gaussian_filter1d(img, si / PIX, axis=0, mode="nearest")
        for r in range(img.shape[0]):
            tmp[r] = gaussian_filter1d(tmp[r], si / (PIX * cosb[r]), mode="nearest")
        out += wi * tmp
    return out


def convolve_cube(cube, B, edges_gev):
    return np.stack([convolve(cube[k], B, edges_gev[k], edges_gev[k + 1]) for k in range(cube.shape[0])])


def check_point_source(L, B, edges_gev, l0=10.0, b0=30.0):
    """Convolve a delta at (l0, b0) in each bin; return the measured R68/R95 on the
    sky against the table's, and the flux ratio."""
    out = []
    i = int(np.argmin(np.abs(B[:, 0] - b0))); j = int(np.argmin(np.abs(L[0] - l0)))
    lr, br = np.radians(L), np.radians(B)
    lc, bc = np.radians(L[i, j]), np.radians(B[i, j])
    dist = np.degrees(np.arccos(np.clip(np.sin(br) * np.sin(bc) + np.cos(br) * np.cos(bc) * np.cos(lr - lc), -1, 1)))
    for k in range(len(edges_gev) - 1):
        img = np.zeros(L.shape); img[i, j] = 1.0
        cimg = convolve(img, B, edges_gev[k], edges_gev[k + 1])
        o = np.argsort(dist.ravel()); c = np.cumsum(cimg.ravel()[o]); c /= c[-1]; d = dist.ravel()[o]
        th, prof = bin_profile(edges_gev[k], edges_gev[k + 1]); ct = containment(th, prof)
        out.append({"bin": k, "flux_ratio": float(cimg.sum()), "R68": float(np.interp(0.68, c, d)), "R68_table": float(np.interp(0.68, ct, th)),
                    "R95": float(np.interp(0.95, c, d)), "R95_table": float(np.interp(0.95, ct, th))})
    return out


def point_sources_exact(l_deg, b_deg, counts, L, B, edges_gev, radius_deg: float = 8.0, inner_deg: float = 1.5, oversample: int = 4):
    """Point sources placed at their exact positions (2026-10-01): each source's
    counts per bin are spread over the pixels within `radius_deg` with the bin's
    three-Gaussian kernel evaluated at the true angular separation of every
    pixel centre from the source (within `inner_deg` the kernel is integrated
    over each pixel on an `oversample` x `oversample` sub-grid, since above
    10 GeV the core is narrower than a pixel), weighted by the pixel solid
    angle and renormalised so the counts on the grid are exactly those given.
    This replaces depositing a source in the pixel containing it and
    convolving, which displaces it by up to half a pixel (0.06 deg): the
    gtsrcmaps benchmark (scripts/srcmaps_benchmark.py) showed that to move
    1-14% of a source's counts across 1- and 2-deg cell boundaries.
    `counts` is (n_sources, n_bins)."""
    l = (np.asarray(l_deg, float) + 180.0) % 360.0 - 180.0; b = np.asarray(b_deg, float)
    counts = np.asarray(counts, float)
    ny, nx = L.shape
    lrow = np.radians(L[0]); brow = np.radians(B[:, 0])
    dl = np.radians(L[0, 1] - L[0, 0]); db = np.radians(B[1, 0] - B[0, 0])
    sub = (np.arange(oversample) + 0.5) / oversample - 0.5          # sub-pixel offsets in pixel units
    out = np.zeros((counts.shape[1],) + L.shape)
    kern = [gaussians(float(edges_gev[k]), float(edges_gev[k + 1])) for k in range(counts.shape[1])]

    def kernel(k, br, lr, bi, li):
        cosd = np.sin(br) * np.sin(bi) + np.cos(br) * np.cos(bi) * np.cos(lr - li)
        d2 = np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0))) ** 2
        w, s = kern[k]
        return sum(wk / (2 * np.pi * sk ** 2) * np.exp(-d2 / (2 * sk ** 2)) for wk, sk in zip(w, s)) * np.cos(br)

    for i in range(len(l)):
        bi = np.radians(b[i]); li = np.radians(l[i])
        r0 = int(np.argmin(np.abs(np.degrees(brow) - b[i]))); c0 = int(np.argmin(np.abs(np.degrees(lrow) - l[i])))
        cosb = max(np.cos(bi), 0.05)
        dr = int(np.ceil(radius_deg / PIX)); dc = int(np.ceil(radius_deg / (PIX * cosb)))
        rs, re = max(r0 - dr, 0), min(r0 + dr + 1, ny); cs, ce = max(c0 - dc, 0), min(c0 + dc + 1, nx)
        ir = int(np.ceil(inner_deg / PIX)); ic = int(np.ceil(inner_deg / (PIX * cosb)))
        irs, ire = max(r0 - ir, 0), min(r0 + ir + 1, ny); ics, ice = max(c0 - ic, 0), min(c0 + ic + 1, nx)
        # oversampled inner box: every pixel split into oversample^2 sub-pixels
        bsub = (brow[irs:ire][:, None] + sub[None] * db).ravel(); lsub = (lrow[ics:ice][:, None] + sub[None] * dl).ravel()
        for k in range(counts.shape[1]):
            ker = kernel(k, brow[rs:re][:, None], lrow[cs:ce][None], bi, li)
            fine = kernel(k, bsub[:, None], lsub[None], bi, li).reshape(ire - irs, oversample, ice - ics, oversample).mean(axis=(1, 3))
            ker[irs - rs:ire - rs, ics - cs:ice - cs] = fine
            tot = ker.sum()
            if tot > 0:
                out[k, rs:re, cs:ce] += counts[i, k] * ker / tot
    return out
