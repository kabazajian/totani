"""Flux-conserving regridding of CAR sky maps onto the analysis grid.

Task 7. Built for the Pohl et al. 2022 ring-based templates, whose files span
four pixel scales (1′, 0.1°, 0.2°, 0.5°), seven footprints, one offset centre
and two *mutually contradictory* longitude conventions — see CONVENTIONS.md
"Pohl et al. 2022 originals". Nothing here imports jax or gcepy.

WHY A HAND-ROLLED REGRIDDER RATHER THAN `reproject`
---------------------------------------------------
Both source and target are CAR (plate carrée) in the *same* frame with no
rotation, so pixel boundaries are lines of constant ℓ and constant b. The
overlap of two pixels therefore **factorizes exactly** into a longitude part and
a latitude part, and the whole regrid is two small matrix products. That is
exact — not an approximation of a spherical polygon intersection — and it is
fast enough to run over the whole template library without a dependency. The
separability is verified against a brute-force reference in the tests.

THE ONE THING THIS MODULE WILL NOT DO SILENTLY
----------------------------------------------
It refuses non-Galactic frames. `Sun_and_Moon_Template_*.fits` are RA/DEC, whose
pixel boundaries are *not* aligned with Galactic ℓ/b, so the separable identity
does not hold and using it would be quietly wrong. Those need genuine spherical
reprojection and are out of scope here (and are in any case tied to Pohl's own
GTIs, not our data selection).

THE SOLID-ANGLE CONVENTION — SETTLED 2026-07-30, AND IT BARELY MATTERS
----------------------------------------------------------------------
Every Pohl map satisfies ``sum(map) × Δℓ × Δb == 1.000000`` using **flat**
coordinate area, across all four pixel scales; true solid angle instead gives
0.99993 (H2) and 0.99756 (boxy bulge). So their *normalization* was done flat.

⚠ **An earlier draft of this docstring claimed the two measures differ by a
≈6% cos b gradient across the ROI that a free normalization could not absorb.
That was wrong, and the retraction is kept here because the wrong version was
plausible.** The 6% is the spread of cos b over |b| ≤ 20°, which governs how a
*total* is computed — but the regrid divides the accumulated flux by the target
pixel weight in the *same* measure, so the measure cancels to first order. Only
its second-order variation *within one source pixel* survives, and sin b ≈ b to
high accuracy over half a degree.

Measured difference between ``solid_angle=True`` and ``False``:

    0.5° → 0.1° grid (nested, aligned)      5e-16   (machine precision)
    0.2° → 0.1° grid (nested, aligned)      5e-16
    Nuclear Bulge, 1′ with OFFSET centre    2e-05
    0.5° → 0.2° grid (2.5×, non-nested)     8e-05
    0.5° → 0.3° grid (non-nested)           2e-04

It cancels **exactly** whenever each target pixel falls inside a single source
pixel — which is the case for the whole Pohl stack on the 0.1° analysis grid —
and even the worst non-nested case is 0.017%, far below the adopted anchor
tolerances (0.1% total, 0.3% per bin). The flag is retained because a Phase-3
grid point with 0.2° pixels *is* non-nested against the 0.5° sources, and
because stating the size of an effect is better than assuming it away.
"""

from __future__ import annotations

import numpy as np

__all__ = ["car_edges", "overlap_1d", "regrid_car", "GRID_GAL_CAR"]

GRID_GAL_CAR = ("GLON-CAR", "GLAT-CAR")


def car_edges(header, axis: int, n: int) -> np.ndarray:
    """World coordinates of the ``n + 1`` pixel edges along ``axis`` (1 or 2).

    FITS pixel *p* (1-based) is centred at ``CRVAL + (p - CRPIX)·CDELT`` and
    spans ``p ± 0.5``, so for 0-based array index ``i`` the edges are at
    ``CRVAL + (i + 0.5 - CRPIX)·CDELT``.

    Returned **as-is, including descending order** when ``CDELT < 0`` — the
    caller must not assume ascending. Handling both orders here is exactly what
    makes the Pohl gas maps (ℓ decreasing) and ICS maps (ℓ increasing) work
    through one code path.
    """
    crval = float(header[f"CRVAL{axis}"])
    crpix = float(header[f"CRPIX{axis}"])
    cdelt = float(header[f"CDELT{axis}"])
    return crval + (np.arange(n + 1) + 0.5 - crpix) * cdelt


def _intervals(edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pixel (lo, hi) pairs from an edge array of either orientation."""
    lo = np.minimum(edges[:-1], edges[1:])
    hi = np.maximum(edges[:-1], edges[1:])
    return lo, hi


def overlap_1d(src_edges: np.ndarray, tgt_edges: np.ndarray,
               wrap: bool = False, measure=None) -> np.ndarray:
    """Overlap matrix ``(n_tgt, n_src)`` in the given measure.

    ``measure`` maps a world coordinate to the conserved variable — identity for
    longitude, ``sin`` (in radians) for latitude under true solid angle. ``wrap``
    adds the ±360° images, so an all-sky source meets a target straddling ℓ = 0.
    """
    slo, shi = _intervals(src_edges)
    tlo, thi = _intervals(tgt_edges)
    if measure is not None:
        slo, shi = measure(slo), measure(shi)
        tlo, thi = measure(tlo), measure(thi)
        slo, shi = np.minimum(slo, shi), np.maximum(slo, shi)
        tlo, thi = np.minimum(tlo, thi), np.maximum(tlo, thi)

    shifts = (0.0,) if (wrap is False or measure is not None) else (-360.0, 0.0, 360.0)
    out = np.zeros((len(tlo), len(slo)))
    for s in shifts:
        lo = np.maximum(tlo[:, None], slo[None, :] + s)
        hi = np.minimum(thi[:, None], shi[None, :] + s)
        out += np.clip(hi - lo, 0.0, None)
    return out


def regrid_car(data: np.ndarray, header, target_header,
               solid_angle: bool = True) -> np.ndarray:
    """Regrid a 2-D or 3-D CAR map onto the target CAR grid, conserving flux.

    For a 3-D cube the leading axis is treated as energy and each plane is
    regridded independently.

    The returned array is in the **same units as the input** (an intensity-like
    value per pixel): total flux ``Σ V·w`` is conserved, where ``w`` is the pixel
    weight in the chosen measure. See the module docstring on `solid_angle`.
    """
    for h, what in ((header, "source"), (target_header, "target")):
        ct = (str(h.get("CTYPE1", "")).strip(), str(h.get("CTYPE2", "")).strip())
        if ct != GRID_GAL_CAR:
            raise ValueError(
                f"{what} grid is {ct}, not {GRID_GAL_CAR}. This regridder is exact "
                "only for axis-aligned Galactic CAR; RA/DEC maps (e.g. Sun/Moon) "
                "need genuine spherical reprojection and are refused rather than "
                "silently mishandled."
            )

    plane = data if data.ndim == 2 else data[0]
    nb, nl = plane.shape
    tnb = int(target_header["NAXIS2"])
    tnl = int(target_header["NAXIS1"])

    meas = (lambda x: np.sin(np.deg2rad(x))) if solid_angle else None
    Ol = overlap_1d(car_edges(header, 1, nl), car_edges(target_header, 1, tnl), wrap=True)
    Ob = overlap_1d(car_edges(header, 2, nb), car_edges(target_header, 2, tnb), measure=meas)

    # target pixel weights in the same measure, for converting flux back to a value
    tl = car_edges(target_header, 1, tnl)
    tb = car_edges(target_header, 2, tnb)
    wl = np.abs(np.diff(tl))
    wb = np.abs(np.diff(meas(tb) if meas else tb))
    denom = wb[:, None] * wl[None, :]
    denom[denom == 0] = np.inf

    def one(p):
        return (Ob @ p @ Ol.T) / denom

    if data.ndim == 2:
        return one(data)
    return np.stack([one(p) for p in data])
