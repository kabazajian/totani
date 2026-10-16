"""The structure library for the ranking programme (plan §4, agreed by Kev
2026-09-19): every large structure in or crossing the Totani ROI, as a
builder `f(target, edges_gev)` returning either a 2-D spatial shape (goes
through the flat-spectrum exposure route, one free norm per bin) or an
(n_bins, ny, nx) per-bin intensity cube (already spectrally integrated; goes
through the single-energy route, like the IC rings).

`STRUCTURES` is the registry the ranking's Stage C selects from. Model-
independent entries are cached by `ranking.structure_cells`; the Fermi-bubble
flat/structured and Loop I shell templates come from `ranking.shared_cells`
because they are part of the standard stack. The Totani-residual bubble
templates are built per diffuse model by `bubble_totani` and are registered
separately in Stage C (Kev: a legitimate candidate, carried with its
injection recovery curve).
"""
from __future__ import annotations

import numpy as np
from astropy.io import fits

from ._paths import GCE_FISHER_ROOT
from .galprop import integrate_bins
from .haslam import haslam_nps_on_target, haslam_on_target, nps_disc
from .templates import grid_lb


# ---- eROSITA bubbles (geometric outer footprint; the brightness map is Oscar's, pending)

def erosita_outer(target, edges_gev):
    from .erosita import erosita_fraction_on_target
    from .templates import bubble_on_target
    return erosita_fraction_on_target(target) * (~(bubble_on_target(target) > 0))


def erosita_outer_north(target, edges_gev):
    L, B = grid_lb(target)
    return erosita_outer(target, edges_gev) * (B > 0)


def erosita_outer_south(target, edges_gev):
    L, B = grid_lb(target)
    return erosita_outer(target, edges_gev) * (B < 0)


# ---- eROSITA X-ray brightness lobes (Oscar's full-sky package, 2026-09-21)
#
# Fitted with INDEPENDENT NORTH AND SOUTH normalisations per bin, which is what
# Oscar recommends and what our own Stage C had already preferred for the
# geometric footprint (`erosita_outer_NS` beat `erosita_outer`). The `unif_*`
# entries are his uniform controls on the identical support: brightness versus
# geometry, everything else held fixed. 🚨 Roughly 45% of the in-ROI weight of
# the brightness maps is the RGB-reconstructed Galactic East -- see
# `erosita.EAST_IS_RECONSTRUCTED` -- so the control is not optional.

def _brightness(key):
    def f(target, edges_gev, _k=key):
        from .erosita import brightness_on_target
        return brightness_on_target(target, _k)
    f.__name__ = f"erosita_{key}"
    return f


# ---- Loop I alternatives (Haslam 408 MHz)

def loopI_haslam_nps(target, edges_gev):
    return haslam_nps_on_target(target)


def loopI_haslam_raw(target, edges_gev):
    return haslam_on_target(target)


# ---- radio loops II-IV (Berkhuijsen 1971 geometry as tabulated by Vidal+ 2015)

#: (l, b, angular radius) in degrees, the CONTINUUM columns of Vidal, Dickinson,
#: Davies & Leahy 2015 (MNRAS 452, 656) Table 1, read from page 5 of
#: `docs/Vidal_2015_radio_loops_1410.4438.pdf` on 2026-09-19 (rule 1a: the
#: PDF, not memory). Loop I is the Wolleben two-shell model in `loopi.py` and
#: is not repeated here. Loops V and VI (127.5, 18.0, 67.2) and (120.5, 30.8,
#: 72.4) are centred far outside the ROI and only their edges could enter;
#: they are listed for completeness and not registered.
#:
#: The same table lists arcs visible only in WMAP polarisation, several INSIDE
#: this ROI: GCS (344.0, 4.8, 18.5), VIIb (0.7, -23.3, 45.9), IX (332.0, 16,
#: 46.5, "inside the NPS"), VIIa (345.0, 0.0, 65). Not registered -- no
#: continuum counterpart, so no gamma-ray expectation -- but flagged for Kev.
RADIO_LOOPS: dict = {
    "loopII": (100.0, -32.5, 45.5),      # Cetus Arc; Large, Quigley & Haslam 1962
    "loopIII": (124.0, 15.5, 32.5),      # Quigley & Haslam 1965
    "loopIV": (315.0, 48.5, 19.8),       # Large, Quigley & Haslam 1966; l = -45, in the northern l < 0 ROI
    "loopV": (127.5, 18.0, 67.2),
    "loopVI": (120.5, 30.8, 72.4),
}


def _thin_shell_projection(L, B, l0, b0, radius_deg, thickness=0.15):
    """A uniform-emissivity spherical shell of angular radius R seen from
    outside, via `loopi.shell_map` with the Sun placed at d = 100 pc so that
    the shell's outer radius subtends R; inner radius (1 - thickness) x outer.
    Path length through the shell: limb-brightened, zero outside R."""
    from .loopi import shell_map
    d = 100.0
    r_out = d * np.sin(np.radians(radius_deg))
    return shell_map(L, B, l0, b0, d, (1.0 - thickness) * r_out, r_out)


def radio_loop(name):
    def build(target, edges_gev):
        l0, b0, R = RADIO_LOOPS[name]
        L, B = grid_lb(target)
        return _thin_shell_projection(L, B, l0, b0, R)
    return build


# ---- solar (and lunar) inverse-Compton halo

#: Part of the A2020 delivery (Pohl et al. 2022), under the data root with the rest of it.
_SUN_IC = GCE_FISHER_ROOT / "data" / "external" / "Pohl_etal_2022" / "Sun_and_Moon_Template_IC.fits"


def sun_moon_ic(target, edges_gev, path=_SUN_IC):
    """The A2020-delivery Sun+Moon IC template: all-sky RA/Dec CAR 0.25 deg,
    46 log-energy planes 31.6 MeV - 1 TeV, time-averaged over 2008-2016 (its
    GTI table). Sampled nearest-pixel onto the Galactic target grid through
    astropy's ICRS <-> Galactic transform, then integrated onto the analysis
    bins (log-log, as the GALPROP cubes). ⚠ Our data run 2008-2023; the
    ecliptic-band morphology is set by the exposure history, so this is the
    right shape with a period caveat; one free norm per bin absorbs the level."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    with fits.open(path) as h:
        cube = h[0].data.astype(float)                      # (nE, n_dec, n_ra)
        hd = h[0].header
        E = np.asarray(h["ENERGIES"].data["Energy"], float)  # MeV
    L, B = grid_lb(target)
    eq = SkyCoord(l=L.ravel() * u.deg, b=B.ravel() * u.deg, frame="galactic").icrs
    ra, dec = eq.ra.deg, eq.dec.deg
    col = np.round((ra - hd["CRVAL1"]) / hd["CDELT1"] + hd["CRPIX1"] - 1).astype(int) % hd["NAXIS1"]
    row = np.clip(np.round((dec - hd["CRVAL2"]) / hd["CDELT2"] + hd["CRPIX2"] - 1).astype(int), 0, hd["NAXIS2"] - 1)
    spec = cube[:, row, col].T                              # (npix, nE)
    spec = np.where(spec < 0, 0.0, spec)
    binned = integrate_bins(spec, E, edges_gev)             # (npix, n_bins)
    return np.transpose(binned.reshape(L.shape + (-1,)), (2, 0, 1))


# ---- gll_iem_v07 (reference model, not a structure; registered here so the cache machinery serves it)

_GIEM = None  # set by ranking once the file is on disk: external/fssc/gll_iem_v07.fits


def gll_iem_v07(target, edges_gev, path=None):
    """The LAT team's interpolated model, all-sky CAR 0.125 deg, 28 planes,
    dN/dE in MeV^-1 cm^-2 s^-1 sr^-1; sampled by exact pixel-centre match onto
    the target (same pixel size and phase; asserted) and bin-integrated."""
    from pathlib import Path
    p = Path(path) if path else Path(__file__).resolve().parents[2] / "external" / "fssc" / "gll_iem_v07.fits"
    with fits.open(p) as h:
        cube = h[0].data.astype(float)
        hd = h[0].header
        E = np.asarray(h["ENERGIES"].data["Energy"], float)
    L, B = grid_lb(target)
    x = (L - hd["CRVAL1"]) / hd["CDELT1"] + hd["CRPIX1"] - 1
    y = (B - hd["CRVAL2"]) / hd["CDELT2"] + hd["CRPIX2"] - 1
    # Longitude: same 0.125 deg grid and half-integer CRPIX1 -> exact match.
    xi = np.rint(x).astype(int) % hd["NAXIS1"]
    if float(np.abs(x - np.rint(x)).max()) > 1e-3:
        raise ValueError("gll_iem_v07 longitude phase differs from the target grid")
    # Latitude: CRPIX2 is an INTEGER (721.0), so its pixel centres sit at
    # b = k x 0.125 while ours sit at (k + 1/2) x 0.125 -- every target centre
    # is exactly halfway between two model rows. Average the two (linear
    # interpolation at the midpoint); asserted rather than assumed.
    frac = y - np.floor(y)
    if not np.allclose(frac, 0.5, atol=1e-3):
        raise ValueError("gll_iem_v07 latitude phase is neither aligned nor half-pixel offset; re-derive")
    y0 = np.clip(np.floor(y).astype(int), 0, hd["NAXIS2"] - 2)
    spec = 0.5 * (cube[:, y0.ravel(), xi.ravel()] + cube[:, (y0 + 1).ravel(), xi.ravel()]).T
    binned = integrate_bins(spec, E, edges_gev)
    return np.transpose(binned.reshape(L.shape + (-1,)), (2, 0, 1))


STRUCTURES = {
    "erosita_bright_N": _brightness("whole_N"),
    "erosita_bright_S": _brightness("whole_S"),
    "erosita_bright_outer_N": _brightness("outer_N"),
    "erosita_bright_outer_S": _brightness("outer_S"),
    "erosita_unif_N": _brightness("unif_whole_N"),
    "erosita_unif_S": _brightness("unif_whole_S"),
    "erosita_unif_outer_N": _brightness("unif_outer_N"),
    "erosita_unif_outer_S": _brightness("unif_outer_S"),
    "erosita_outer": erosita_outer,
    "erosita_outer_N": erosita_outer_north,
    "erosita_outer_S": erosita_outer_south,
    "loopI_haslam_nps": loopI_haslam_nps,
    "loopI_haslam_raw": loopI_haslam_raw,
    "loopII": radio_loop("loopII"),
    "loopIII": radio_loop("loopIII"),
    "loopIV": radio_loop("loopIV"),
    "sun_moon_ic": sun_moon_ic,
    "gll_iem_v07": gll_iem_v07,
}
