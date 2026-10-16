"""The eROSITA Bubble boundary as a region definition, and matched-theta weights.

The boundary is Oscar Macias' reconstruction of the solid outline in Zhang et
al. (arXiv:2408.06312v2) Extended Data Figure 5(a), delivered 2026-09-08 in
`external/erosita_validation/`. It is a GEOMETRIC FOOTPRINT, not an X-ray
brightness map: it says where the edge runs and nothing about what is inside.
So it never enters the fit as a template. It defines two regions, and this
module builds the weights that make comparing them a test.

🚨 THE NULL IS AGAIN A SYMMETRY, NOT A MODEL -- but a weaker one than the
reflection test in `scripts/patch_test.py` uses. There, reflecting a patch in
l -> -l left theta_GC exactly invariant pixel by pixel. Here the two regions
are genuinely different shapes, so they are matched in theta_GC only in
DISTRIBUTION: inside and outside are reweighted, annulus by annulus, to carry
the same theta_GC distribution. Any spherically symmetric intensity then has
the same response-corrected expectation on both sides, and a lobe bounded by
the eROSITA edge need not. `spherical_contrast` measures how well that
cancellation actually holds and is asserted in the tests.

⚠ THE TEST HAS NO SUPPORT INSIDE theta_GC ~ 34.5 deg. Inside that radius every
pixel of the halo search region is inside the eROSITA footprint, so there is
no outside comparison at the same angular distance and those annuli drop out
entirely. Measured here 2026-09-08 (34.6 deg north, 34.4 deg south); Oscar's
independent quadrature in `matched_theta_geometric_support.csv` gives 34.62 and
34.38. That is most of the halo's leverage gone -- `results/
band_forecast_2026-09-06.json` puts 54-61% of it in the 10-30 vs 30-60 deg
latitude contrast -- so this bounds a lobe's morphology, and is a poor
constraint on the halo's amplitude. Say so wherever the numbers are quoted.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

#: Oscar's delivery. Gitignored (80 MB); `boundary_path` raises with the
#: provenance if it is absent rather than letting a caller fit a zero map.
_BOUNDARY = Path(__file__).resolve().parents[2] / "external" / "erosita_validation" \
    / "eROSITA_boundary_Macias.fits"

#: sha256 of the ED Figure 5(a) JPEG the reconstruction was measured from,
#: verified 2026-09-08 byte-for-byte against the copy embedded on page 28 of
#: our own `docs/arXiv_2408.06312v2.pdf`. Checked against the file's SRCSHA so
#: a re-delivery built from a different raster cannot pass silently.
SOURCE_SHA256 = "2f3f5accf707449c00a6e41ceca979692f07f9cd795d7bbdc7872b6276f459b6"

#: Area-weighted mean longitude of each lobe over the halo search region, in
#: degrees, as delivered. The two have OPPOSITE SIGN -- the northern lobe leans
#: to negative l and the southern to positive -- which is the cheapest decisive
#: guard against an l-flipped delivery, the trap `CONVENTIONS.md` opens with.
#: A longitude flip swaps both signs; nothing else in the pipeline would.
LOBE_MEAN_LON_DEG = {"north": -9.461, "south": +3.752}


def boundary_path(path=None) -> Path:
    if path is not None:
        return Path(path)
    if not _BOUNDARY.exists():
        raise FileNotFoundError(
            f"{_BOUNDARY} is missing. It is Oscar Macias' eROSITA boundary "
            "reconstruction (delivered 2026-09-08), gitignored at 80 MB, so a "
            "fresh clone does not have it. Ask him, or restore it from the "
            "machine it was delivered to; do not substitute a traced curve.")
    return _BOUNDARY


def erosita_fraction_on_target(target, path=None) -> np.ndarray:
    """Per-pixel fraction of solid angle inside the eROSITA boundary, on `target`.

    🚨 THIS IS AN EXACT PIXEL-FOR-PIXEL TRANSFER, NOT A RESAMPLING, AND IT IS
    CHECKED. Oscar's map is all-sky CAR at 0.125 deg with pixel centres on the
    same half-pixel phase as the Totani grid, so every target pixel centre
    lands exactly on a source pixel centre (asserted: max |x - round(x)| = 0).
    Nearest-neighbour is therefore not an approximation here and no flux is
    moved. If a future delivery changes CDELT or CRPIX this raises instead of
    quietly interpolating -- which for a partial-coverage FRACTION map would
    bias the boundary by up to half a pixel.

    The FRACTION extension is used, not the 0/1 PRIMARY: 5,209 pixels straddle
    the boundary and their partial coverage is the difference between an edge
    good to 0.125 deg and one good to the stroke width.
    """
    from astropy.wcs import WCS

    from .templates import grid_lb

    with fits.open(boundary_path(path)) as h:
        if h[0].header.get("SRCSHA") != SOURCE_SHA256:
            raise ValueError(
                f"SRCSHA is {h[0].header.get('SRCSHA')!r}, expected "
                f"{SOURCE_SHA256!r} -- this file was reconstructed from a "
                "different raster than the one verified against our own copy "
                "of arXiv:2408.06312v2 page 28. Re-verify before using it.")
        frac = h["FRACTION"].data.astype(float)
        hdr = h[0].header

    if not (0.0 <= frac.min() and frac.max() <= 1.0):
        raise ValueError(f"FRACTION is not in [0, 1]: {frac.min()}..{frac.max()}")

    L, B = grid_lb(target)
    x, y = WCS(hdr).wcs_world2pix(L.ravel(), B.ravel(), 0)
    off = max(float(np.abs(x - np.rint(x)).max()), float(np.abs(y - np.rint(y)).max()))
    if off > 1e-6:
        raise ValueError(
            f"target pixel centres miss the boundary map's centres by {off:.4f} "
            "pixels; the grids no longer share a phase, so nearest-neighbour "
            "sampling of a partial-coverage map would bias the edge")
    xi, yi = np.rint(x).astype(int), np.rint(y).astype(int)
    if not ((xi >= 0) & (xi < frac.shape[1]) & (yi >= 0) & (yi < frac.shape[0])).all():
        raise ValueError("target pixels fall outside the boundary map; it should be all-sky")
    return frac[yi, xi].reshape(L.shape)


def check_handedness(L, B, inside, omega, atol_deg=0.5) -> dict:
    """Assert the delivered footprint has the published l-asymmetry.

    Returns the measured area-weighted mean longitude per lobe. Raises if
    either differs from `LOBE_MEAN_LON_DEG` by more than `atol_deg`, which an
    l-flip cannot survive: the two lobes lean opposite ways, by 9.5 and 3.8
    deg, so a flip moves both by twice that and changes both signs.
    """
    out = {}
    for lobe, sign in (("north", +1), ("south", -1)):
        m = (sign * B >= 10.0)
        w = omega * inside * m
        out[lobe] = float((w * L).sum() / w.sum())
    bad = {k: v for k, v in out.items()
           if abs(v - LOBE_MEAN_LON_DEG[k]) > atol_deg}
    if bad:
        raise ValueError(
            f"eROSITA footprint fails the handedness check: measured {out}, "
            f"expected {LOBE_MEAN_LON_DEG} within {atol_deg} deg. Both signs "
            "flipping means the delivered map is mirrored in longitude.")
    return out


def theta_gc(L, B) -> np.ndarray:
    """Angular distance from the Galactic Centre, degrees."""
    return np.degrees(np.arccos(np.clip(
        np.cos(np.radians(B)) * np.cos(np.radians(L)), -1.0, 1.0)))


def matched_theta_weights(L, B, omega, inside, keep, response, *,
                          theta_step_deg=0.25, node_weight="min"):
    """Per-pixel weights making inside and outside share a theta_GC distribution.

    🚨 THE MATCHING IS RESPONSE-WEIGHTED, NOT GEOMETRIC, AND THE DIFFERENCE IS
    NOT COSMETIC. The statistic these weights build sums COUNTS. A spherical
    component of profile f contributes `A f(theta_i) R_i` counts to pixel i,
    where `R_i = exposure_i x omega_i` is `flat_sr` -- the same response every
    template in `stack.build_background_stack` is built through. Cancelling it
    needs `sum_in c_i R_i = sum_out d_i R_i` per annulus, so the annulus
    normalisation is by RESPONSE-weighted area, not solid angle.

    Matching on solid angle instead was the first version of this function, on
    2026-09-08, and it left the fitted NFW halo predicting 329 more counts
    outside than inside at 20.8 GeV against a signal of 960 -- a 34% leak of
    exactly the component the test is supposed to be blind to. It was caught by
    watching `halo_predicted_diff`, which is now ~0.1 counts and is written to
    the JSON every run for that reason. `response` is therefore required, not
    optional, and it is per energy bin because the exposure is.

    In each annulus of width `theta_step_deg` the two sides are separately
    normalised and then given the SAME annulus weight, so a spherical intensity
    contributes equally to both sums and cancels in the difference whatever its
    radial profile. Annuli where either side is empty carry no weight: they
    cannot be matched, so they are dropped rather than compared. That is the
    whole of the theta_GC > 34.5 deg restriction -- it is not a cut anyone chose.

    `node_weight` sets the annulus weight. "min" is min(R_in, R_out), the
    response-weighted analogue of Oscar's choice in
    `matched_theta_geometric_support.csv`. "harmonic" is
    2 R_in R_out / (R_in + R_out), the variance-optimal weighting for a lobe of
    uniform surface brightness; it lies within a factor two of "min" everywhere,
    so the choice moves the significance by a few percent, not more.

    Returns `(w_in, w_out, info)`, both dimensionless per-pixel multipliers on
    counts. A pixel that straddles the boundary carries weight on BOTH sides.
    """
    keep = np.asarray(keep).astype(bool)
    return matched_theta_weights_regions(
        L, B, omega, inside * keep, (1.0 - inside) * keep, response,
        theta_step_deg=theta_step_deg, node_weight=node_weight)


def matched_theta_weights_regions(L, B, omega, f_in, f_out, response, *,
                                  theta_step_deg=0.25, node_weight="min"):
    """`matched_theta_weights` for two ARBITRARY fractional-coverage maps.

    `f_in` and `f_out` are per-pixel coverage fractions in [0, 1] of the two
    regions to be matched in theta_GC. `matched_theta_weights` is the special
    case `f_in = inside * keep`, `f_out = (1 - inside) * keep`. This general
    form exists for the spatial breakdown (`scripts/boundary_breakdown.py`,
    2026-09-18): a SUB-REGION of the inside is matched against the SAME full
    outside, so every sub-region's statistic refers to one common reference
    and the sub-regions are comparable with each other. The two maps may
    overlap on boundary pixels; that is handled as before (a pixel then
    carries weight on both sides).

    The matching, the normalisation and the `info` dict are unchanged; see
    `matched_theta_weights` for the response-weighting argument, which applies
    verbatim.
    """
    if node_weight not in ("min", "harmonic"):
        raise ValueError(f"node_weight must be 'min' or 'harmonic', got {node_weight!r}")

    th = theta_gc(L, B)
    f_in = np.asarray(f_in, float)
    f_out = np.asarray(f_out, float)
    # Response-weighted area per side: counts per unit spherical intensity.
    r_in = f_in * response
    r_out = f_out * response

    node = np.floor(th / theta_step_deg).astype(int)
    n_node = int(node.max()) + 1
    R_in = np.bincount(node.ravel(), r_in.ravel(), minlength=n_node)
    R_out = np.bincount(node.ravel(), r_out.ravel(), minlength=n_node)

    supported = (R_in > 0) & (R_out > 0)
    if not supported.any():
        raise ValueError("no theta annulus has both an inside and an outside part")
    if node_weight == "min":
        W = np.where(supported, np.minimum(R_in, R_out), 0.0)
    else:
        W = np.where(supported, 2 * R_in * R_out / np.maximum(R_in + R_out, 1e-300), 0.0)
    W = W / W.sum()

    sup_pix = supported[node]
    # Scale so that `sum_in c_i R_i` equals the inside region's OWN total
    # response over the supported annuli. The statistic is then in counts on
    # the inside region's scale -- directly comparable to the fitted halo
    # counts inside -- rather than in the units of 1/(exposure x solid angle)
    # that the annulus normalisation leaves behind. Ratios and significances
    # are unaffected by this choice; only the readability of the counts is.
    scale = float(r_in[sup_pix].sum())

    with np.errstate(divide="ignore", invalid="ignore"):
        w_in = scale * W[node] * np.where(R_in[node] > 0, f_in / R_in[node], 0.0)
        w_out = scale * W[node] * np.where(R_out[node] > 0, f_out / R_out[node], 0.0)

    omega_in = omega * f_in
    omega_out = omega * f_out
    th_sup = np.where(supported)[0] * theta_step_deg
    info = {
        "theta_step_deg": theta_step_deg,
        "node_weight": node_weight,
        "supported_nodes": int(supported.sum()),
        "theta_min_deg": float(th_sup.min()),
        "theta_max_deg": float(th_sup.max() + theta_step_deg),
        "inside_area_sr": float(omega_in[sup_pix].sum()),
        "outside_area_sr": float(omega_out[sup_pix].sum()),
        "inside_area_dropped_sr": float(omega_in[~sup_pix].sum()),
        "outside_area_dropped_sr": float(omega_out[~sup_pix].sum()),
        "response_matched_in": float((w_in * response).sum()),
        "response_matched_out": float((w_out * response).sum()),
        "max_weight": float(max(w_in.max(), w_out.max())),
    }
    return w_in, w_out, info


#: Spherical intensity profiles used to verify the cancellation. Deliberately
#: unlike each other and unlike NFW: if the matching only worked for profiles
#: resembling the halo it would prove nothing about a lobe.
SPHERICAL_PROFILES = {
    "constant": lambda t: np.ones_like(t),
    "exponential_20deg": lambda t: np.exp(-t / 20.0),
    "cored_inverse_square_8deg": lambda t: 1.0 / (1.0 + (t / 8.0) ** 2),
    "gaussian_25deg": lambda t: np.exp(-0.5 * (t / 25.0) ** 2),
    "smooth_oscillatory": lambda t: 1.0 + 0.5 * np.sin(np.radians(3.0 * t)),
}


def spherical_contrast(w_in, w_out, theta, response, profiles=None) -> dict:
    """Relative inside-minus-outside contrast of each spherical profile, in counts.

    `response` must be the same `flat_sr` the weights were built with: the
    profiles are compared as they would enter the statistic, in expected
    COUNTS, not as intensities on the sky. Measuring the contrast without it
    would report a cancellation the statistic does not actually have -- which
    is how the solid-angle version of `matched_theta_weights` looked correct
    while leaking a third of the halo.

    A constant profile cancels to rounding: that part is exact. A VARYING
    profile leaves the curvature across one annulus, of order 1e-5 at the
    0.25 deg default (see `tests/test_erosita.py`, which also shows 0.25 deg is
    the optimum step). Compare against `unmatched_contrast`, the same quantity
    with the two regions taken as they are: that is of order 100%, because
    inside and outside sit at very different mean theta_GC. The matching is not
    a refinement, it is the whole test.
    """
    profiles = SPHERICAL_PROFILES if profiles is None else profiles
    out = {}
    for name, f in profiles.items():
        v = f(theta) * response
        a = float((w_in * v).sum())
        b = float((w_out * v).sum())
        out[name] = 2 * (a - b) / (a + b) if (a + b) else 0.0
    return out


def unmatched_contrast(inside, keep, omega, theta, response, profiles=None) -> dict:
    """`spherical_contrast` for the two regions with no theta matching at all."""
    profiles = SPHERICAL_PROFILES if profiles is None else profiles
    f_in = inside * keep
    f_out = (1.0 - inside) * keep
    out = {}
    for name, f in profiles.items():
        v = f(theta) * response
        i = float((f_in * v).sum() / f_in.sum())
        o = float((f_out * v).sum() / f_out.sum())
        out[name] = 2 * (i - o) / (i + o)
    return out


# ---- Oscar's full-sky X-ray BRIGHTNESS templates (delivered 2026-09-21) --------

#: The full-sky brightness package. Unlike `_BOUNDARY`, which is a geometric
#: footprint, these carry X-ray surface brightness and so enter the fit as
#: templates. Gitignored (1.4 GB).
_BRIGHTNESS = Path(__file__).resolve().parents[2] / "external" \
    / "eROSITA_fullsky_lat_handoff" / "erosita_brightness_fullsky_normalized.fits"

#: 🚨 HALF OF THIS MAP IS NOT DATA. The Galactic West (-180 < l < 0) is the
#: German consortium's native 0.6-1.0 keV survey. The Galactic East (0 < l <
#: 180) belongs to the Russian consortium and is not public, so Oscar
#: reconstructed it from the published MPE RGB image, calibrated on the West.
#: Our ROI is |l| <= 60, split almost evenly, so ~45% of the in-ROI weight of
#: every one of these templates rests on that reconstruction. The package
#: header says so itself: QUALSTAT = RECONSTRUCTED_APPROXIMATION.
#:
#: Measured here 2026-09-21 (`scripts/erosita_brightness_checks.py`): the
#: East/West mean surface brightness ratio over the halo search region is
#: 0.911 north and 0.819 south, and it is STABLE at the few-per-cent level
#: across every coherent alternative Oscar ships -- calibration percentiles,
#: smoothing, registration shifts, residual-shape transport -- including
#: `FULL_ED1`, which is registered and calibrated from a DIFFERENT published
#: display (0.953 / 0.797). So the contrast is not an artefact of one
#: calibration choice. What none of those variants can bound is a common-mode
#: error in the RGB -> intensity inversion shared by all of them; ED1 is the
#: only partial handle on it. Always fit the uniform controls alongside.
EAST_IS_RECONSTRUCTED = True

#: Extension names. `T_*` whole footprint, `O_T_*` outer only, `U_*`/`O_U_*`
#: the matching uniform (geometry-only) controls, which are Oscar's own
#: brightness-versus-geometry test.
BRIGHTNESS_EXTENSIONS = {
    "whole_N": "T_N", "whole_S": "T_S",
    "outer_N": "O_T_N", "outer_S": "O_T_S",
    "unif_whole_N": "U_N", "unif_whole_S": "U_S",
    "unif_outer_N": "O_U_N", "unif_outer_S": "O_U_S",
}

#: Area-weighted mean longitude of the UNIFORM controls over the halo search
#: region, degrees. These are pure geometry, so they must reproduce the
#: independently delivered `LOBE_MEAN_LON_DEG` of the 2026-09-08 boundary --
#: measured -9.670 / +3.748 against -9.461 / +3.752, i.e. the same Zhang et
#: al. outline. That agreement is the handedness guard for this delivery: a
#: longitude flip would swap both signs and nothing else would.
UNIFORM_MEAN_LON_DEG = {"north": -9.670, "south": +3.748}


def brightness_path(path=None) -> Path:
    if path is not None:
        return Path(path)
    if not _BRIGHTNESS.exists():
        raise FileNotFoundError(
            f"{_BRIGHTNESS} is missing. It is Oscar Macias' full-sky eROSITA "
            "brightness package (delivered 2026-09-21), gitignored at 1.4 GB. "
            "Ask him for the handoff archive; do not substitute the geometric "
            "boundary, which carries no brightness information.")
    return _BRIGHTNESS


def brightness_on_target(target, key: str, path=None) -> np.ndarray:
    """One lobe of Oscar's brightness package on `target`, in sr^-1.

    `key` is a name in `BRIGHTNESS_EXTENSIONS`. As for the boundary map this
    is an exact pixel-for-pixel transfer, asserted, not a resampling: the
    delivery is CAR 0.125 deg on CRPIX 1440.5/720.5 and the Totani grid is CAR
    0.125 deg on CRPIX 480.5, so both axes share a half-pixel phase and every
    target centre lands on a source centre.

    Each lobe carries unit solid-angle integral over its OWN support before any
    ROI cut, so the returned array integrates to less than one on our grid.
    That is correct and must not be renormalised: the fit gives each lobe a
    free normalisation per bin anyway, and renormalising after an ROI change
    would silently redefine the flux convention the README warns about.
    """
    from astropy.wcs import WCS

    from .templates import grid_lb

    if key not in BRIGHTNESS_EXTENSIONS:
        raise KeyError(f"{key!r} is not one of {sorted(BRIGHTNESS_EXTENSIONS)}")
    ext = BRIGHTNESS_EXTENSIONS[key]
    with fits.open(brightness_path(path), memmap=True) as h:
        hdr = h[ext].header
        if hdr.get("NATIVE", False):
            raise ValueError(
                f"{ext} claims NATIVE=T; this loader documents and expects the "
                "reconstructed-approximation delivery. Re-read the README.")
        img = np.asarray(h[ext].data, float)

    L, B = grid_lb(target)
    x, y = WCS(hdr).wcs_world2pix(L.ravel(), B.ravel(), 0)
    off = max(float(np.abs(x - np.rint(x)).max()), float(np.abs(y - np.rint(y)).max()))
    if off > 1e-6:
        raise ValueError(
            f"target pixel centres miss the brightness map's centres by "
            f"{off:.4f} pixels; the grids no longer share a phase, so "
            "nearest-neighbour sampling would move flux across the lobe edge")
    xi, yi = np.rint(x).astype(int), np.rint(y).astype(int)
    if not ((xi >= 0) & (xi < img.shape[1]) & (yi >= 0) & (yi < img.shape[0])).all():
        raise ValueError("target pixels fall outside the brightness map; it should be all-sky")
    out = img[yi, xi].reshape(L.shape)
    if not np.isfinite(out).all() or out.min() < 0:
        raise ValueError(f"{ext} sampled to non-finite or negative values on the target")
    return out
