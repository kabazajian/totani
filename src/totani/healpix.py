"""HEALPix RING-scheme pixelisation, without healpy.

The environment has no healpy, and two external maps are HEALPix: the Haslam
408 MHz map (NSIDE 512) and the GALPROP skymaps from Ackermann+ 2012 (NSIDE
128). This is the standard algorithm (Gorski et al. 2005), transcribed from
the reference implementation, and tested by an exact round trip over every
pixel at several NSIDE including 512 (`tests/test_haslam.py`).

Sampling elsewhere is nearest-pixel via `ang2pix_ring`. That is exact for a
map whose pixels are finer than the target grid (Haslam: 6.9' vs 7.5') and a
step function for one whose pixels are coarser (GALPROP: 27' vs 7.5'); the
fit works on 10-degree cells, so neither moves any flux that matters, but
the coarser case is said where it is used.
"""
from __future__ import annotations

import numpy as np

__all__ = ["ang2pix_ring", "pix2ang_ring", "npix"]


def npix(nside: int) -> int:
    return 12 * int(nside) * int(nside)


def ang2pix_ring(nside: int, l_deg, b_deg) -> np.ndarray:
    """RING pixel index of Galactic (l, b) in degrees; vectorised.

    Transcribed from the reference algorithm (Gorski+ 2005; healpy's
    `ang2pix_ring`). theta is colatitude = 90 - b; phi = l in [0, 360).
    """
    nside = int(nside)
    theta = np.radians(90.0 - np.asarray(b_deg, float))
    phi = np.radians(np.asarray(l_deg, float) % 360.0)
    z = np.cos(theta)
    za = np.abs(z)
    tt = (phi / (0.5 * np.pi)) % 4.0            # in [0, 4)
    npix = 12 * nside * nside
    ncap = 2 * nside * (nside - 1)
    pix = np.empty(z.shape, dtype=np.int64)

    eq = za <= 2.0 / 3.0
    if eq.any():
        t1 = nside * (0.5 + tt[eq])
        t2 = nside * z[eq] * 0.75
        jp = np.floor(t1 - t2).astype(np.int64)      # ascending edge line
        jm = np.floor(t1 + t2).astype(np.int64)      # descending edge line
        ir = nside + 1 + jp - jm                     # ring number, 1..2nside+1
        kshift = 1 - (ir & 1)
        ip = (jp + jm - nside + kshift + 1) // 2
        ip = ip % (4 * nside)
        pix[eq] = ncap + (ir - 1) * 4 * nside + ip

    po = ~eq
    if po.any():
        ttp = tt[po]
        tp = ttp - np.floor(ttp)
        tmp = nside * np.sqrt(3.0 * (1.0 - za[po]))
        jp = np.floor(tp * tmp).astype(np.int64)
        jm = np.floor((1.0 - tp) * tmp).astype(np.int64)
        ir = jp + jm + 1                             # ring number, 1..nside-1
        ip = np.floor(ttp * ir).astype(np.int64) % (4 * ir)
        north = z[po] > 0
        out = np.where(north, 2 * ir * (ir - 1) + ip,
                       npix - 2 * ir * (ir + 1) + ip)
        pix[po] = out
    return pix


def pix2ang_ring(nside: int, pix) -> tuple:
    """Galactic (l, b) in degrees of RING pixel centres; vectorised."""
    nside = int(nside)
    pix = np.asarray(pix, dtype=np.int64)
    npix = 12 * nside * nside
    ncap = 2 * nside * (nside - 1)
    z = np.empty(pix.shape, float)
    phi = np.empty(pix.shape, float)

    nc = pix < ncap
    if nc.any():
        p = pix[nc]
        iring = ((1 + np.floor(np.sqrt(1 + 2 * p + 0.0))) // 2).astype(np.int64)
        # guard the integer sqrt at ring boundaries
        iring = np.where(2 * iring * (iring - 1) > p, iring - 1, iring)
        iring = np.where(2 * (iring + 1) * iring <= p, iring + 1, iring)
        iphi = p + 1 - 2 * iring * (iring - 1)
        z[nc] = 1.0 - iring ** 2 / (3.0 * nside ** 2)
        phi[nc] = (iphi - 0.5) * np.pi / (2.0 * iring)

    eq = (pix >= ncap) & (pix < npix - ncap)
    if eq.any():
        ip = pix[eq] - ncap
        iring = ip // (4 * nside) + nside
        iphi = ip % (4 * nside) + 1
        fodd = 0.5 * (1 + ((iring + nside) & 1))
        z[eq] = (2 * nside - iring) / (1.5 * nside)
        phi[eq] = (iphi - fodd) * np.pi / (2.0 * nside)

    sc = pix >= npix - ncap
    if sc.any():
        ip = npix - pix[sc]
        iring = ((1 + np.floor(np.sqrt(2 * ip - 1 + 0.0))) // 2).astype(np.int64)
        iring = np.where(2 * iring * (iring - 1) >= ip, iring - 1, iring)
        iring = np.where(2 * (iring + 1) * iring < ip, iring + 1, iring)
        iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1))
        z[sc] = -1.0 + iring ** 2 / (3.0 * nside ** 2)
        phi[sc] = (iphi - 0.5) * np.pi / (2.0 * iring)

    b = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    l = (np.degrees(phi) + 180.0) % 360.0 - 180.0
    return l, b


