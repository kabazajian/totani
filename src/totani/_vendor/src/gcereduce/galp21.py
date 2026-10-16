"""Galp21 ingest — Youyou Li's GALPROP v57 SA/ISRF maps onto the analysis grid.

Seam-roll + cutout + six-ring sums, then a `MapCubeFunction`-ready FITS per
(component × ring group). Tier 1 (`scripts/galp21_tier1.py`) validated the
delivery; this turns it into templates.

**Role, RULED (Kev, PI, 2026-08-15):** λ directions for the **final KMT limit**,
not a grid family — they touch neither the Metric A pass matrix, nor the
16(ii)/A25 family choice, nor the 16(iii) proxies, so the blind selection is
undisturbed.

🚨 THE SEAM. `CRVAL1 = 180` at `CRPIX1 = 1800.5` with `CDELT1 = +0.1`, so
longitude runs 0.05° → 359.95° across the array and **the Galactic Centre sits
on the wrap, between the last column and the first**. Native column `c` carries

    l(c) = 180 + (c + 1 - 1800.5) * 0.1        [degrees, 0-based c]

so |l| <= L is columns `[3600 - 10L, 3600)` followed by `[0, 10L)`. This module
takes that concatenation DIRECTLY rather than `np.roll`-ing the full cube: a
roll materialises all 380 MB per file, and the six-ring sums touch 17 of them.

🚨 THE HANDEDNESS, AND IT GOES THE OTHER WAY FROM EVERY OTHER PRODUCT HERE.
GALPROP writes `CDELT1 = +0.1` — longitude INCREASING with column, which is
gcepy's convention and the opposite of every fermitools product. Tier 1 asserts
it per file rather than inheriting it from one model, because this project has
met four longitude-handedness traps and the one delivery needing no flip is not
a reason to stop looking.

**On write we flip to `CDELT1 = -0.1`**, matching `build_grid_extras.write_j_fits`
and every other map this pipeline feeds `gtsrcmaps`. Both describe the same sky;
what makes it safe is that the WCS is written explicitly and read back in the
test. ⚠ Never index these arrays without going through the WCS.

RING SUMS ARE EXACT MAP ADDITION -- no interpolation, no GALPROP re-run. The six
requested boundaries fall on native ring edges, so the sum is addition on
identical grids. That is what takes Youyou's holiday off the critical path.

🚨 THE CO/H2 RING GRID GOVERNS ALL THREE GAS PHASES, INCLUDING THE ONES LABELLED
"HI" AND "HII". `gen_skymaps.cc` reads `HIR` then `COR` into the same shared
`galaxy.R_bins`, and COR is read second, so the HI-fit edges are discarded
(SETUP_PARAMETERS.md §8 — Youyou's finding, and exactly the provenance detail
that is impossible to recover later). The HI grid differs from CO/H2 by up to
~30 pc in the inner rings; we use CO/H2 because that is what the output used.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

#: 🚨 THE CORRECTED BATCH, 2026-08-20 (Kev). The previous delivery had the H2 gas
#: phase identically zero in all six models -- a GALPROP bug in the CAR/mapcube
#: path that Youyou Li diagnosed and fixed
#: (`validation/galp21_h2_missing_2026-08-15.md`). These are the re-run maps.
#:
#: The path changed with the batch: `galprop-youyou/` is gone and the models sit
#: directly under `maps/` rather than under `maps/extracted/`. The old location
#: is NOT a fallback -- it holds the defective batch, and silently reading it
#: would reintroduce the empty-H2 defect that Tier 1 check 6 exists to catch.
#:
#: ⚠ `/Volumes/auto/...` is the automounter spelling and is what resolves in a
#: shell here; the same share is `/Volumes/Data2/...` when mounted directly.
BASE = Path("/Volumes/auto/Data2/gce-fisher/GCE_Galprop_maps_2026/maps")

MODELS = ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12")

NATIVE_SHAPE = (44, 600, 3600)
NATIVE_PIXEL_DEG = 0.1
NATIVE_CRVAL1, NATIVE_CRPIX1, NATIVE_CDELT1 = 180.0, 1800.5, 0.1
N_RINGS = 17

#: The CO/H2 ring outer edges, kpc (SETUP_PARAMETERS.md §8). Ring `i` (1-based)
#: spans `(RING_EDGES_KPC[i-1], RING_EDGES_KPC[i])`.
RING_EDGES_KPC = (0.000, 1.485, 1.993, 2.494, 3.003, 3.500, 3.999, 4.496,
                  5.004, 5.499, 6.505, 7.000, 7.991, 10.00, 11.50, 16.50,
                  19.00, 50.00)

#: The six groups, as 1-based native ring indices. Boundaries land on native
#: edges exactly -- asserted by `check_ring_grouping()`, not assumed.
SIX_RING_GROUPS = {
    1: (1,),                       # 0.000 - 1.485
    2: (2, 3),                     # 1.485 - 2.494
    3: (4, 5),                     # 2.494 - 3.500
    4: (6, 7, 8, 9, 10, 11, 12),   # 3.500 - 7.991
    5: (13,),                      # 7.991 - 10.00
    6: (14, 15, 16, 17),           # 10.00 - 50.00
}

#: A2020's six-ring boundaries, which the grouping reproduces to <1%. Recorded
#: so the ~1% offsets are a stated fact absorbed by the free per-ring
#: normalizations, not a silent approximation (delivery doc Q2: retain the
#: native boundaries, never interpolate).
A2020_RING_EDGES_KPC = (0.0, 1.5, 2.5, 3.5, 8.0, 10.0, 50.0)

GAS_PHASES = ("H2R", "HIR", "HII")
N_ISRF_COMPONENTS = 3

#: 🚨 MEASURED 2026-08-15: the H2 (molecular hydrogen) gas phase is IDENTICALLY
#: ZERO in this delivery -- every ring, both pion decay and bremsstrahlung, all
#: six models. And it is NOT an export defect: the 17 per-ring maps sum to the
#: delivered `*_mapcube` total at ratio 1.000000, so the totals carry no H2
#: either. The GALDEFs are configured for it (`nH2_model = 9`, `COR_filename`
#: set, 17 non-zero per-ring `X_CO_values`, `gamma_rays = 2` = "compute HI,H2
#: skymaps separately"), so this is a defect in the runs, not a modelling choice.
#:
#: Molecular gas is a major pi0 source in the inner Galaxy -- the CMZ and the
#: molecular ring -- which is exactly our ROI. A background family without it is
#: not usable, so this is REFUSED rather than silently summed to zero.
#: Record: `validation/galp21_h2_missing_2026-08-15.md`. Owner: Youyou Li.
#:
#: ✅ RESOLVED 2026-08-20 -- the guard is retired for the CORRECTED delivery.
#: Youyou's root cause: GALPROP's CAR/mapcube path read the CO/H2 cube's
#: reference coordinate with the H I convention (centre vs first pixel).
#: The corrected batch at GCE_Galprop_maps_2026/maps passes Tier 1 check 6
#: (H2 non-zero, all six models, `runs/galp21_tier1.json`), and BASE points
#: there. The old batch (galprop-youyou/) keeps the defect -- if BASE is ever
#: pointed back at it, check 6 is the tripwire, not this constant.
EMPTY_PHASES = ()


def check_ring_grouping() -> dict:
    """The six groups tile rings 1..17 exactly and land on native edges."""
    seen = [r for g in SIX_RING_GROUPS.values() for r in g]
    if sorted(seen) != list(range(1, N_RINGS + 1)):
        raise ValueError(f"the six groups do not tile rings 1..{N_RINGS}: {seen}")
    edges, offsets = [RING_EDGES_KPC[0]], []
    for g in sorted(SIX_RING_GROUPS):
        rings = SIX_RING_GROUPS[g]
        if tuple(rings) != tuple(range(rings[0], rings[-1] + 1)):
            raise ValueError(f"group {g} is not contiguous: {rings}")
        edges.append(RING_EDGES_KPC[rings[-1]])
    for got, want in zip(edges, A2020_RING_EDGES_KPC):
        offsets.append(0.0 if want == 0 else abs(got - want) / want)
    if max(offsets) > 0.011:
        raise ValueError(
            f"a group boundary is {max(offsets):.1%} from A2020's, beyond the "
            f"~1% the delivery documents; the native ring grid has changed")
    return {"edges_kpc": edges, "a2020_edges_kpc": list(A2020_RING_EDGES_KPC),
            "max_fractional_offset": max(offsets)}


def component_path(model: str, kind: str, *, ring: int | None = None,
                   phase: str | None = None, isrf: int | None = None) -> Path:
    """Path to one native mapcube.

    `kind` is `pion_decay`, `bremss` or `ics_anisotropic`.

    ⚠ **Per-ring inverse Compton exists only for the ANISOTROPIC calculation.**
    `ics_isotropic_mapcube` ships as a total with no ring split. The two are
    alternative treatments of the same emission, not additive components, and
    the anisotropic one is both the better calculation and the only one that can
    be decomposed -- so it is what the ring sums use. Stated because summing the
    two would double the IC.
    """
    if kind in ("pion_decay", "bremss"):
        if phase not in GAS_PHASES:
            raise ValueError(f"{kind} needs phase in {GAS_PHASES}, got {phase!r}")
        if not 1 <= (ring or 0) <= N_RINGS:
            raise ValueError(f"ring must be 1..{N_RINGS}, got {ring}")
        stem = f"{kind}_{phase}_mapcube_comp_{ring - 1}"      # comp index is 0-based
    elif kind == "ics_anisotropic":
        if not 1 <= (ring or 0) <= N_RINGS:
            raise ValueError(f"ring must be 1..{N_RINGS}, got {ring}")
        if not 0 <= (isrf if isrf is not None else -1) < N_ISRF_COMPONENTS:
            raise ValueError(f"isrf must be 0..{N_ISRF_COMPONENTS - 1}, got {isrf}")
        stem = f"ics_anisotropic_comp_ring_{ring}_mapcube_comp_{isrf}"
    else:
        raise ValueError(f"unknown kind {kind!r}")
    return BASE / model / f"{stem}_57_{model}"


def cut_columns(half_deg: float) -> tuple[slice, slice]:
    """The two native column slices whose concatenation is |l| <= `half_deg`.

    Returned as a pair because the GC is on the array's wrap: the high-longitude
    half comes FIRST (l just below 360, i.e. just below 0) and the low half
    second, which puts l increasing through zero across the join.
    """
    n = int(round(half_deg / NATIVE_PIXEL_DEG))
    if n <= 0 or 2 * n > NATIVE_SHAPE[2]:
        raise ValueError(f"half_deg {half_deg} is outside the delivered grid")
    return slice(NATIVE_SHAPE[2] - n, NATIVE_SHAPE[2]), slice(0, n)


def cut_rows(half_deg: float) -> slice:
    """Rows for |b| <= `half_deg` on the 600-row (|b| <= 30) native grid."""
    n = int(round(half_deg / NATIVE_PIXEL_DEG))
    mid = NATIVE_SHAPE[1] // 2
    if n <= 0 or n > mid:
        raise ValueError(
            f"half_deg {half_deg} exceeds the delivered |b| <= "
            f"{mid * NATIVE_PIXEL_DEG}; the maps do not extend further and "
            f"padding them would be inventing template data (rule 3)")
    return slice(mid - n, mid + n)


def cutout(path: Path, *, half_l: float, half_b: float) -> np.ndarray:
    """`(44, 2*half_b*10, 2*half_l*10)` float64 cut about the GC.

    Reads through `memmap` and slices, so only the cut is materialised -- the
    six-ring sums touch 17 of these per component and a full roll would move
    380 MB each time. Cast to float64 on the way out: the cubes are float32 and
    this project has already had a float32 accumulator overflow to `inf`.
    """
    from astropy.io import fits

    lo, hi = cut_columns(half_l)
    rows = cut_rows(half_b)
    with fits.open(path, memmap=True) as h:
        d = h[0].data
        if tuple(d.shape) != NATIVE_SHAPE:
            raise ValueError(f"{path.name}: shape {d.shape} != {NATIVE_SHAPE}")
        left = np.asarray(d[:, rows, lo], dtype=np.float64)
        right = np.asarray(d[:, rows, hi], dtype=np.float64)
    return np.concatenate([left, right], axis=-1)


def energies_mev(model: str) -> np.ndarray:
    """The 44 plane energies, read from a delivered ENERGIES extension."""
    from astropy.io import fits

    p = component_path(model, "pion_decay", ring=1, phase="H2R")
    with fits.open(p, memmap=True) as h:
        col = h["ENERGIES"].data
        return np.asarray(col[col.names[0]], float)


def ring_group_sum(model: str, kind: str, group: int, *, half_l: float,
                   half_b: float, phase: str | None = None,
                   isrf: str | int = "sum") -> np.ndarray:
    """One six-ring group, summed over its native rings — exact map addition.

    For `ics_anisotropic`, `isrf="sum"` adds the three ISRF components (the
    delivery-doc recommendation: downstream there is one IC amplitude per ring,
    matching the Pohl six-ring convention we compare against). An integer keeps
    one component instead.
    """
    if group not in SIX_RING_GROUPS:
        raise ValueError(f"group must be 1..6, got {group}")
    if phase in EMPTY_PHASES:
        raise ValueError(
            f"phase {phase!r} is IDENTICALLY ZERO in this delivery -- every "
            f"ring, both pion decay and bremsstrahlung, all six models -- and "
            f"the per-ring maps sum to the delivered total at ratio 1.000000, "
            f"so the totals carry no H2 either. Summing it would produce a gas "
            f"template silently missing the molecular hydrogen that dominates "
            f"inner-Galaxy pi0 emission. See "
            f"validation/galp21_h2_missing_2026-08-15.md; this is Youyou's to "
            f"diagnose and must not be worked around here (CLAUDE.md rule 3).")
    acc = None
    for ring in SIX_RING_GROUPS[group]:
        if kind == "ics_anisotropic":
            comps = (range(N_ISRF_COMPONENTS) if isrf == "sum" else (int(isrf),))
            for c in comps:
                a = cutout(component_path(model, kind, ring=ring, isrf=c),
                           half_l=half_l, half_b=half_b)
                acc = a if acc is None else acc + a
        else:
            a = cutout(component_path(model, kind, ring=ring, phase=phase),
                       half_l=half_l, half_b=half_b)
            acc = a if acc is None else acc + a
    return acc


def write_mapcube(out: Path, cube: np.ndarray, energies_mev_arr,
                  *, pixel_deg: float = NATIVE_PIXEL_DEG,
                  history: str | None = None) -> Path:
    """Write a `MapCubeFunction`-ready Galactic-CAR FITS.

    🚨 **Flipped in longitude on the way out.** The cut arrives l-INCREASING
    with column (GALPROP's `CDELT1 = +0.1`, which is gcepy's convention); this
    writes `CDELT1 = -pixel_deg` with the array reversed, matching
    `build_grid_extras.write_j_fits` and every other map this pipeline hands to
    `gtsrcmaps`. Same sky, one convention. The WCS is explicit and the test
    reads it back rather than trusting this docstring.
    """
    from astropy.io import fits

    cube = np.asarray(cube, float)
    if cube.ndim != 3:
        raise ValueError(f"cube must be 3-D (energy, b, l), got {cube.shape}")
    e = np.asarray(energies_mev_arr, float)
    if e.size != cube.shape[0]:
        raise ValueError(f"{e.size} energies for {cube.shape[0]} planes")
    n_b, n_l = cube.shape[1], cube.shape[2]

    hdr = fits.Header()
    hdr["CTYPE1"], hdr["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    hdr["CRVAL1"], hdr["CRVAL2"] = 0.0, 0.0
    hdr["CDELT1"], hdr["CDELT2"] = -pixel_deg, pixel_deg
    hdr["CRPIX1"], hdr["CRPIX2"] = (n_l + 1) / 2.0, (n_b + 1) / 2.0
    hdr["CTYPE3"] = "photon energy"
    hdr["CUNIT3"] = "MeV"
    # ⚠ NO BUNIT is written because the delivery carries none. The setup doc
    # states cm^-2 sr^-1 s^-1 MeV^-1 and Tier 1's counts check corroborates it
    # to 85% (right with point sources excluded), but that is not the same as a
    # keyword we were given -- inventing one would be inventing provenance.
    hdr["HISTORY"] = "Galp21 ingest: seam cut about the GC, six-ring sum"
    hdr["HISTORY"] = "l FLIPPED to CDELT1<0 to match this pipeline's other maps"
    if history:
        hdr["HISTORY"] = history[:70]

    primary = fits.PrimaryHDU(cube[:, :, ::-1].astype("float32"), header=hdr)
    ecol = fits.Column(name="Energy", format="D", unit="MeV", array=e)
    ehdu = fits.BinTableHDU.from_columns([ecol], name="ENERGIES")
    out.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([primary, ehdu]).writeto(out, overwrite=True)
    return out


def longitude_of_column(header, col: int) -> float:
    """l for a 0-based column, from the header's own WCS. Never index directly."""
    return float(header["CRVAL1"]
                 + (col + 1 - header["CRPIX1"]) * header["CDELT1"])
