"""Every background template, on the 960x960 Totani grid.

    tgt = target_header()
    gas = gas_dust_on_target(tgt)            # 4 HI + 4 H2 + 2 dust residuals
    ics = ics_on_target(tgt, EDGES_GEV)      # 6 rings, integrated onto the bins
    bub = bubble_on_target(tgt)              # Ackermann+ 2014 flat template
    loo = loopi_on_target(tgt)               # 2 shells, generated here

🚨 THE FILE NAMING, THE VARIANT CHOICE AND THE HEADER FIXES ARE NOT REDERIVED
HERE. They took a delivery from Oscar Macias and two rounds of correspondence to
pin down, and every one of them is a place where a plausible-looking wrong file
would produce a plausible-looking wrong answer. So this module IMPORTS
`a2020_templates_official` (CONFIGS, load_map, ICS_4RING_GROUPS) and
`a2020_stage2_diffuse.ics_on_bins` -- the VENDORED, pinned copies of the
gce-fisher scripts in `totani/_vendor/` (see `_paths.py`) -- rather than
restating them. What is new here is only the TARGET GRID.

The three things that would silently break if they were rederived:

1. **The baseline variant is `pohl` gas + `mag5` dust at T_S = 150 K** (Oscar,
   2026-08-11). Not guessable: in the normalized delivery the UNTAGGED files are
   the 170 K set and the baseline 150 K set carries the `T150K` tag, which is
   backwards from what the names suggest.
2. **`H2_pohl_4annuli.fits` has a half-pixel latitude error in its header**
   (`CRVAL2 = -90.0` against -89.75 everywhere else) and we apply the
   correction. Settled by his 2026-08-13 delivery, whose H2 is bit-identical in
   pixel data on the common grid -- the header was the defect. Every headline
   A2020 number uses this fix.
3. **The gas and dust maps declare bare `GLON`/`GLAT` with no projection code**,
   which `regrid_car` refuses by design. `load_map` sets `GLON-CAR`, guarded by
   a full-sky span check so a cutout can never be assumed to be plate carree.
   Reading these files with `fits.getdata` instead would hit that refusal -- or
   worse, be "fixed" by someone setting the code without the span guard.

⚠ **THE SOLID-ANGLE MEASURE, RE-EXAMINED FOR THIS LATITUDE RANGE.** `regrid_car`
defaults to `solid_angle=True` and the validated A2020 path uses the default, so
we match it. The regridder's own argument for why the measure barely matters --
that flux is divided by the target pixel weight in the SAME measure, so only its
second-order variation WITHIN one source pixel survives -- was made over
|b| <= 20, where cos b varies by 6%. Here cos b varies by a factor of two. The
argument still holds, because what survives is the curvature of sin(b) across a
single 0.5 deg source pixel, not across the ROI. Recorded because inheriting the
conclusion without re-checking its premise is how a latitude-dependent bias gets
in.
"""
from __future__ import annotations

import numpy as np
from pathlib import Path
from astropy.io import fits

from . import ensure_scripts_on_path
from .loopi import loopi_maps

ensure_scripts_on_path()

from gcereduce.regrid import car_edges, regrid_car          # noqa: E402
from a2020_templates_official import (           # noqa: E402
    CONFIGS,
    ICS_4RING_GROUPS,
    load_map,
)
from a2020_stage2_diffuse import ics_on_bins     # noqa: E402

__all__ = ["target_header", "grid_lb", "axis_coords", "pixel_solid_angle",
           "delivered_paths", "gas_dust_on_target",
           "ics_on_target", "ics_four_ring", "bubble_on_target", "bubble_macias_on_target",
           "loopi_on_target", "isotropic_on_target", "CONFIGS",
           "ICS_4RING_GROUPS"]

#: Where A2020's own delivery lives, relative to the gce-fisher checkout.
_DELIVERY = "data/external/Templates_arxiv2003.10416"
#: The Ackermann+ 2014 bubble template, Mollweide copy -- see `bubble_on_target`.
_BUBBLE = "data/external/FermiBubblesPublicFiles/bubbleTemp_wcs.fits"

#: Oscar Macias' structured all-sky bubble map, delivered 2026-09-06 in reply
#: to the Galp21 proposal. Lives in THIS repo's `external/`, not gce-fisher's.
_BUBBLE_MACIAS = "external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits"


def target_header(npix: int = 960, pixel_deg: float = 0.125) -> fits.Header:
    """The analysis grid: Galactic CAR, centred on (0, 0), longitude to the LEFT.

    ⚠ `CDELT1` is NEGATIVE, matching what `gtbin` writes and every other map
    this pipeline feeds. gcepy's arrays run the other way; this project has met
    four longitude-handedness traps, so never index these arrays by hand --
    go through the WCS, which `grid_lb` does.
    """
    h = fits.Header()
    h["NAXIS"], h["NAXIS1"], h["NAXIS2"] = 2, npix, npix
    h["CTYPE1"], h["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    h["CRVAL1"], h["CRVAL2"] = 0.0, 0.0
    h["CDELT1"], h["CDELT2"] = -pixel_deg, pixel_deg
    h["CRPIX1"] = h["CRPIX2"] = (npix + 1) / 2.0
    return h


def axis_coords(header) -> tuple:
    """1-D (l_c, b_c) pixel-centre coordinate VECTORS -- the length-nx and
    length-ny axes separately, NOT the full (ny, nx) meshgrid `grid_lb` builds.

    🚨 THIS IS THE SHAPE `gcereduce.dm_templates.j_map` ACTUALLY WANTS.
    `_los_nodes` does `L = asarray(grid_l)[None, :]`, `B = asarray(grid_b)[:,
    None]` and broadcasts them together -- i.e. it builds the 2-D grid
    INTERNALLY from two 1-D vectors, exactly like `np.meshgrid`. Passing the
    full 2-D `grid_lb()` output instead is not merely inefficient: each
    already-2-D array picks up ANOTHER broadcast axis inside `_los_nodes`, so
    an (ny, nx) input silently becomes an (ny, nx, nx)-ish intermediate.
    Measured 2026-09-04: this took RSS to 13.85 GB and climbing on what should
    have been a small row-chunk, on a MACHINE SHARED WITH ANOTHER SESSION'S
    30+ GB queue -- caught and killed before it repeated the 2026-08-20
    OOM/reboot incident, not after.
    """
    l_edges = car_edges(header, 1, int(header["NAXIS1"]))
    b_edges = car_edges(header, 2, int(header["NAXIS2"]))
    l_c = (l_edges[:-1] + l_edges[1:]) / 2.0
    b_c = (b_edges[:-1] + b_edges[1:]) / 2.0
    return (l_c + 180.0) % 360.0 - 180.0, b_c


def grid_lb(header) -> tuple:
    """(l, b) of every pixel centre. Longitudes wrapped to [-180, 180).

    🚨 COMPUTED FROM THE WCS KEYWORDS DIRECTLY, VIA `gcereduce.regrid.car_edges`,
    AND NOT THROUGH `astropy.wcs`. That is a deliberate reversal of the first
    version of this function, which used astropy and was silently wrong.

    A2020's delivered gas and dust maps carry `CRVAL2 = -89.75` with
    `CRPIX2 = 1`, and span 0-360 deg in longitude. Once `load_map` supplies the
    missing `-CAR` projection code, FITS-WCS reads that as a plate carree whose
    reference point is 89.75 deg off the equator and whose native longitudes run
    past 180 -- both outside the projection's valid domain. astropy then returns
    **NaN**, and measured on 2026-09-02 it returned NaN for every latitude in the
    HI cube, so an ROI mask selected ZERO pixels and the flux check reported
    "source flux = 0" rather than an error. Moving the reference to the equator
    fixed half of it; the longitude wrap still produced NaN for 129,240 of
    259,200 pixels.

    `car_edges` is the project's own validated reading of these headers -- it is
    what `regrid_car` uses, it is tested against a brute-force reference, and it
    handles both longitude handedness conventions through one path. Using it
    here means the coordinates and the regrid can never disagree about what a
    header means, which is the failure this whole class of bug comes from.

    ⚠ This still honours "never index without going through the WCS": the sign
    of CDELT, the reference pixel and the reference value are all read from the
    header. What is dropped is astropy's spherical rotation machinery, which
    these plate-carree maps do not need and which their headers do not satisfy.
    """
    ny, nx = int(header["NAXIS2"]), int(header["NAXIS1"])
    l_edges = car_edges(header, 1, nx)
    b_edges = car_edges(header, 2, ny)
    l_c = 0.5 * (l_edges[:-1] + l_edges[1:])
    b_c = 0.5 * (b_edges[:-1] + b_edges[1:])
    L, B = np.meshgrid(l_c, b_c)
    L = (L + 180.0) % 360.0 - 180.0
    if not (np.isfinite(L).all() and np.isfinite(B).all()):
        raise ValueError("non-finite pixel coordinates from the header")
    if np.abs(B).max() > 90.0 + 1e-9:
        raise ValueError(
            f"pixel latitudes reach {np.abs(B).max():.4f} deg, past the pole. "
            "The header is not the plate carree this assumes.")
    return L, B


def delivered_paths(root, config: str = "baseline", ts: str = "150K") -> dict:
    """A2020's delivered gas/dust files for one variant.

    The pattern is `a2020_templates_official.main`'s, kept in one place so the
    two cannot drift: `HI_{gas}_4annuli_TS{ts}`, `H2_{gas}_4annuli`,
    `posres_{gas}_{mag}_{ts}`, `negres_{gas}_{mag}_{ts}`.
    """
    if config not in CONFIGS:
        raise ValueError(f"unknown config {config!r}; expected {sorted(CONFIGS)}")
    if ts not in ("150K", "170K"):
        raise ValueError(f"T_S must be '150K' or '170K', not {ts!r} -- there is "
                         "no default; A2020's baseline is 150K (Oscar, "
                         "2026-08-11) but it must be said explicitly")
    gas, mag = CONFIGS[config]
    d = root / _DELIVERY
    return {
        "HI": d / f"HI_{gas}_4annuli_TS{ts}.fits",
        "H2": d / f"H2_{gas}_4annuli.fits",
        "dust_pos": d / f"posres_{gas}_{mag}_{ts}.fits",
        "dust_neg": d / f"negres_{gas}_{mag}_{ts}.fits",
    }


def gas_dust_on_target(target, root=None, config: str = "baseline",
                       ts: str = "150K") -> dict:
    """The 4 HI annuli, 4 H2 annuli and 2 dust residuals, on the target grid.

    Energy-independent spatial templates: each carries one free normalization
    per energy bin in the fit, so no spectrum is applied here.
    """
    from . import GCE_FISHER_ROOT

    root = GCE_FISHER_ROOT if root is None else root
    paths = delivered_paths(root, config, ts)
    out = {}
    for key in ("HI", "H2"):
        # `fix_h2` corrects H2_pohl's half-pixel CRVAL2; it is a no-op on any
        # other file and asserts rather than shifting an already-corrected one.
        data, hdr, _ = _load(paths[key], fix_h2=True)
        if data.ndim != 3 or data.shape[0] != 4:
            raise ValueError(f"{paths[key].name} is {data.shape}, expected 4 annuli")
        for k in range(4):
            out[f"{key}_ring{k + 1}"] = regrid_car(data[k], hdr, target)
    for key in ("dust_pos", "dust_neg"):
        data, hdr, _ = _load(paths[key], fix_h2=True)
        out[key] = regrid_car(data if data.ndim == 2 else data[0], hdr, target)
    return out


def _load(path, fix_h2: bool):
    """`load_map` with a stable 3-tuple return, whatever it hands back."""
    got = load_map(path, fix_h2=fix_h2)
    if isinstance(got, tuple):
        if len(got) == 3:
            return got
        if len(got) == 2:
            return got[0], got[1], {}
    raise TypeError(f"load_map returned {type(got)} with unexpected shape; "
                    "the upstream signature changed and this wrapper must "
                    "be updated rather than guessed at")


def ics_on_target(target, edges_gev, root=None) -> list:
    """The six IC rings, regridded and integrated onto the analysis bins.

    Delegates to `a2020_stage2_diffuse.ics_on_bins`, which owns the log10-E
    integration. ⚠ That integration assumes the third axis samples a density in
    log10 E; no BUNIT confirms it, and the assumption is recorded upstream. Each
    ring carries a free per-bin normalization, so it affects only the in-bin
    weighting, not the scale.
    """
    from . import GCE_FISHER_ROOT

    root = GCE_FISHER_ROOT if root is None else root
    return ics_on_bins(np.asarray(edges_gev, float), target,
                       src_dir=root / _DELIVERY)


def ics_four_ring(six: list) -> list:
    """Collapse Pohl's six IC rings to A2020's four: (1+2+3), 4, 5, 6.

    CONFIRMED by Oscar 2026-08-11 -- it had been our inference. Valid only
    because the IC maps are in PHYSICAL units rather than individually
    normalized; the gas maps in the public release ARE normalized and summing
    those would silently misweight the inner Galaxy.
    """
    if len(six) != 6:
        raise ValueError(f"expected 6 IC rings, got {len(six)}")
    return [sum(six[i - 1] for i in group) for group in ICS_4RING_GROUPS]


def bubble_on_target(target, root=None) -> np.ndarray:
    """The Ackermann+ 2014 flat bubble template, sampled onto the target grid.

    🚨 SAMPLED FROM THE MOLLWEIDE COPY, NOT THE HEALPIX ONE, AND DELIBERATELY.
    The HEALPix file (nside=64) carries NO `ORDERING` keyword, so RING vs NESTED
    is unstated; reading it under the wrong one gives a scattered map that a free
    normalization would happily absorb into the fit. `healpy` is also absent
    here. The Mollweide copy has an unambiguous FITS-WCS that astropy inverts,
    so this route avoids the guess entirely.

    Nearest-neighbour sampling is correct HERE because the template is a binary
    0/1 mask -- Ackermann+ 2017 Table 1 calls it a "Flat template" -- so there is
    no flux to conserve, only a boundary. It would be wrong for any of the
    intensity maps, which go through `regrid_car` instead.
    """
    from astropy.wcs import WCS

    from . import GCE_FISHER_ROOT

    root = GCE_FISHER_ROOT if root is None else root
    with fits.open(root / _BUBBLE) as h:
        data, hdr = h[0].data.astype(float), h[0].header
    values = np.unique(data[np.isfinite(data)])
    if not np.all(np.isin(values, (0.0, 1.0))):
        raise ValueError(
            f"bubble template is not a 0/1 mask (values {values[:5]}...); "
            "nearest-neighbour sampling would no longer be flux-correct")

    L, B = grid_lb(target)
    x, y = WCS(hdr).wcs_world2pix(L.ravel(), B.ravel(), 0)
    xi, yi = np.rint(x).astype(int), np.rint(y).astype(int)
    ok = ((xi >= 0) & (xi < data.shape[1]) & (yi >= 0) & (yi < data.shape[0]))
    if not ok.all():
        raise ValueError(f"{(~ok).sum()} target pixels fall outside the "
                         "bubble map's projection; it should be all-sky")
    out = np.nan_to_num(data[yi, xi]).reshape(L.shape)
    return out


def bubble_macias_on_target(target, path=None) -> np.ndarray:
    """Oscar Macias' structured Fermi bubble map on the target grid.

    🚨 PROVENANCE CORRECTED 2026-09-18 (Oscar, by email): this is Figure 12 of
    **Ackermann et al. 2017 (arXiv:1704.03910)**, NOT the 2014 bubbles paper,
    which every note before 2026-09-18 says. Its morphology was derived from
    **1-10 GeV** data by a 3-component spectral decomposition (bubble-like,
    softer background, MSP-like), documented in sec. 2.5 and footnote 4 of
    Macias et al. 2019. So it should be at its BEST at low energy, and above
    ~10 GeV the morphology is an extrapolation -- the opposite of the
    "derived above a few GeV" reading the low-energy caveat had assumed.

    The STRUCTURED template, not the flat one, with the point-source holes
    inpainted by Oscar (HDU 1 `INPAINT` flags the filled
    pixels: 1 inside the inner 40 deg, 2 elsewhere). All-sky CAR, 0.5 deg,
    CDELT1 < 0 like every fermitools product; no BUNIT, which is fine because
    the fit gives it one free normalization per bin.

    🚨 THIS IS AN INTENSITY, SO IT GOES THROUGH `regrid_car`, NOT THE
    NEAREST-NEIGHBOUR ROUTE `bubble_on_target` USES. That route is correct for
    a 0/1 mask and refuses anything else by design; a structured map must be
    regridded conserving flux or the 0.5 -> 0.125 deg resampling would
    misweight every edge.

    Tier 1 (2026-09-07, `scripts/bubble_macias.py`): full ROI coverage, no
    NaN, handedness consistent with the Ackermann mask (pixel correlation
    0.375 as-is vs 0.329 l-flipped; same sign of l-asymmetry in both
    hemispheres). 94% of the flat footprint is non-zero in this map, and 51%
    of this map's flux lies OUTSIDE that footprint -- Oscar's point about the
    shells not capturing the structure, measured.
    """
    from gcereduce.regrid import regrid_car

    root = Path(__file__).resolve().parents[2]
    path = root / _BUBBLE_MACIAS if path is None else Path(path)
    with fits.open(path) as h:
        data, hdr = h[0].data.astype(float), h[0].header
    if not np.isfinite(data).all() or (data < 0).any():
        raise ValueError(f"{path.name}: non-finite or negative values; not an intensity map")
    out = regrid_car(data, hdr, target)
    if np.isnan(out).any():
        raise ValueError("regrid produced NaN inside the ROI; the map is not all-sky")
    return out


def loopi_on_target(target) -> dict:
    """Both Loop I shells on the target grid; two independent templates."""
    L, B = grid_lb(target)
    return loopi_maps(L, B)


def isotropic_on_target(target) -> np.ndarray:
    """The isotropic component's spatial shape: uniform, sec 2.3.

    A flat array of ones. Genuinely dimensionless -- the absolute isotropic
    intensity is exactly what f_l fits per bin (Totani initializes it at
    E^2 dN/dE = 1e-4 MeV cm^-2 s^-1 sr^-1, but that is an MCMC STARTING POINT,
    not a spatial template value).
    """
    ny, nx = int(target["NAXIS2"]), int(target["NAXIS1"])
    return np.ones((ny, nx))


def pixel_solid_angle(header) -> np.ndarray:
    """Solid angle of every pixel, in steradians -- (ny, nx), constant along rows.

    Exact for CAR: dOmega = dl_rad * (sin(b_hi) - sin(b_lo)) per pixel, which is
    dl*db*cos(b) to second order. Ranges over a factor of TWO across this ROI
    (cos 60 = 0.5).

    🚨 THE FACTOR THE FIRST TEMPLATE STACK OMITTED ENTIRELY -- found by the
    validation step, 2026-09-04, not by inspection. Every diffuse template
    (gas, dust, ICS, bubble, Loop I, isotropic, halo J) is an INTENSITY-like
    map: a value per sky DIRECTION, whose counts contribution per pixel is
    value x exposure x THIS. Omitting it is a latitude-dependent SHAPE error of
    up to 2x across |b| <= 60 -- exactly the class of error the module
    docstrings said would "not be harmless", while wrongly asserting
    `regrid_car` had it covered (`regrid_car` conserves flux while REGRIDDING;
    it never converts intensity to counts -- and the A2020 path never faced
    this because gtsrcmaps applies solid angle internally). The tell was the
    fitted isotropic amplitude jumping five orders of magnitude between
    adjacent bins while point sources -- the one component that needs NO
    solid angle -- fit at f = 0.64-0.98.
    """
    l_edges = np.deg2rad(car_edges(header, 1, int(header["NAXIS1"])))
    b_edges = np.deg2rad(car_edges(header, 2, int(header["NAXIS2"])))
    dl = np.abs(np.diff(l_edges))
    domega_b = np.abs(np.diff(np.sin(b_edges)))
    return np.outer(domega_b, dl)
