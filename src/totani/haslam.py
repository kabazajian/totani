"""The Haslam 408 MHz map as an EMPIRICAL Loop I template -- the alternative to
the two geometric Wolleben shells in `loopi.py`.

Kev's decision, 2026-09-18: after the structured bubble map and the geometric
Loop I shells were found to disagree about the low-energy shape of the
northern l < 0 sky (`results/rebaseline_candidate_noLoopI_2026-09-18.md`),
try a radio-synchrotron template for Loop I alongside the shells. The map is
the Remazeilles et al. (2015, MNRAS 451, 4311; arXiv:1411.3268)
destriped/desourced reprocessing of Haslam et al. (1982), from LAMBDA:

    https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/haslam408_dsds_Remazeilles2014.fits
    sha256 a1d36df6305fa39fd461de462990bdc19b1e9bbbe438125b0eef5956b07c8e3e
    HEALPix NSIDE 512 RING, Galactic, brightness temperature in K, 56' beam

⚠ WHAT THIS IS AND IS NOT. 408 MHz brightness traces synchrotron: the
cosmic-ray electron density times the magnetic field, integrated along the
line of sight. The gamma rays from Loop I are inverse Compton off the same
electrons times the radiation field. Related morphology, not the same one --
which is exactly why the LAT team (Ackermann+ 2014, 2017) and Totani used a
geometric model instead. This template is a MODELLING CHOICE, scored as a
variant; it is not a drop-in replacement for the shells.

Two shapes are offered:
  * `haslam_on_target`           -- the raw map over the whole target grid;
  * `haslam_nps_on_target`       -- the map restricted to the disc subtended
                                    by Wolleben shell 2 (the North Polar Spur
                                    shell; see `nps_disc` for why not both).
The restricted form is the one that answers "alternate Loop I"; the raw form
is a synchrotron template competing with the IC rings over the whole ROI.

🚨 NO healpy HERE. The environment does not have it, so the RING-scheme
pixelisation is implemented below (the standard HEALPix algorithm, Gorski et
al. 2005) and tested by an exact round trip over every pixel of the map in
`tests/test_haslam.py`. Sampling is nearest-pixel: the map's 6.9' pixels are
finer than the 7.5' target grid and far finer than the 56' beam, so no flux
is moved by sampling.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

#: LAMBDA delivery, gitignored (12.6 MB); `haslam_path` raises with the URL if
#: it is absent rather than letting a caller fit a zero map.
_HASLAM = Path(__file__).resolve().parents[2] / "external" / "haslam" \
    / "haslam408_dsds_Remazeilles2014.fits"
HASLAM_URL = ("https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/"
              "haslam408_dsds_Remazeilles2014.fits")
HASLAM_SHA256 = "a1d36df6305fa39fd461de462990bdc19b1e9bbbe438125b0eef5956b07c8e3e"


def haslam_path(path=None) -> Path:
    if path is not None:
        return Path(path)
    if not _HASLAM.exists():
        raise FileNotFoundError(
            f"{_HASLAM} is missing. It is the Remazeilles+ 2015 destriped/"
            f"desourced Haslam 408 MHz map, gitignored; download it from "
            f"{HASLAM_URL} (sha256 {HASLAM_SHA256}).")
    return _HASLAM


from .healpix import ang2pix_ring, pix2ang_ring  # noqa: F401  (re-exported; tests import them from here)


# ---- the map ------------------------------------------------------------------

def load_haslam(path=None) -> tuple:
    """(temperature_K per RING pixel, nside). Asserts the header is what the
    sampler assumes: RING, Galactic, full sky."""
    with fits.open(haslam_path(path)) as h:
        hd = h[1].header
        if hd.get("ORDERING") != "RING" or hd.get("COORDSYS") != "GALACTIC":
            raise ValueError(f"expected RING/GALACTIC, got "
                             f"{hd.get('ORDERING')}/{hd.get('COORDSYS')}")
        nside = int(hd["NSIDE"])
        data = np.asarray(h[1].data[hd["TTYPE1"]], float).ravel()
    if data.size != 12 * nside * nside:
        raise ValueError(f"{data.size} pixels is not full sky at NSIDE {nside}")
    if not np.isfinite(data).all() or (data <= 0).any():
        raise ValueError("Haslam map has non-finite or non-positive temperatures")
    return data, nside


def haslam_on_target(target, path=None) -> np.ndarray:
    """408 MHz brightness temperature (K) at every target pixel centre."""
    from .templates import grid_lb

    data, nside = load_haslam(path)
    L, B = grid_lb(target)
    return data[ang2pix_ring(nside, L, B)]


def nps_disc(L, B) -> np.ndarray:
    """The sky disc subtended by Wolleben shell 2 -- the North Polar Spur
    shell -- as a boolean map.

    ⚠ WHY SHELL 2 AND NOT "THE SHELLS' FOOTPRINT". Shell 1 has the Sun INSIDE
    its outer sphere (d = 78 pc < r_out = 81 pc), so its path length is
    non-zero along every sight line: its footprint is the whole sky and cannot
    restrict anything (found by `tests/test_haslam.py`, 2026-09-18). Shell 2
    (centre l = 332, b = 37 deg, d = 95 pc, r_out = 82 pc) has the Sun outside
    it, so it subtends a disc of angular radius arcsin(82/95) = 59.7 deg. That
    disc is the geometric model's own statement of where Loop I's bright
    northern structure is.
    """
    from .loopi import SHELLS  # noqa: F401  -- documents the source of the numbers

    l0, b0, d, r_out = 332.0, 37.0, 95.0, 82.0
    rho = np.degrees(np.arcsin(r_out / d))
    cos_sep = (np.sin(np.radians(B)) * np.sin(np.radians(b0))
               + np.cos(np.radians(B)) * np.cos(np.radians(b0))
               * np.cos(np.radians(L - l0)))
    return np.degrees(np.arccos(np.clip(cos_sep, -1, 1))) <= rho


def haslam_nps_on_target(target, path=None) -> np.ndarray:
    """The Haslam map inside shell 2's disc, zero outside -- an empirical
    SHAPE for the geometric model's own Loop I REGION. One template where the
    shells were two."""
    from .templates import grid_lb

    L, B = grid_lb(target)
    return haslam_on_target(target, path) * nps_disc(L, B)
