#!/usr/bin/env python
"""Build A2020's gas / dust / IC templates from the authors' OWN delivered maps.

    .venv/bin/python scripts/a2020_templates_official.py --config baseline --ts 150K

WHAT THESE MAPS ARE
-------------------
`data/external/Templates_arxiv2003.10416/` is the template set Oscar Macias sent
on 2026-08-11, in reply to our asking what A2020's diffuse models were built
from. It is the set used for arXiv:2003.10416 itself, and it replaces our earlier
reconstruction from the public Pohl et al. (2022) release for the gas, dust and
inverse-Compton components. His description of the naming, quoted so the mapping
below can be audited without the email:

  * hydrodynamic gas maps carry "pohl"; the standard gas maps are "interpolated";
  * pick ONE set of HI maps at a fixed spin temperature T_S;
  * each HI set has matching positive (posres) and negative (negres) dust
    residual maps at the SAME T_S, each at two magnitude cuts;
  * exactly one posres and one negres go into a given fit;
  * "The six-ring IC maps are also included. To reproduce the four-ring IC setup
    used in the analysis, combine the first three rings, so the four components
    are (1+2+3),4,5,6."

That last line CONFIRMS the grouping we had inferred and were using. It is no
longer an assumption.

WHAT THIS SCRIPT DOES NOT COVER
-------------------------------
Only the components in the delivery: HI, H2, the two dust residuals, and IC.
The Fermi bubbles, Loop I, the bulges, the isotropic and the point sources are
untouched by this and still come from where `a2020_stage2_diffuse.py` gets them.
Sun and Moon remain omitted, as recorded in `results/a2020_stage2_fit_spec.md`.

🚨 TWO CHOICES ARE NOT DETERMINED BY THE DELIVERY AND ARE NOT GUESSED HERE
--------------------------------------------------------------------------
`--ts` is REQUIRED and has no default. Oscar said to pick one fixed T_S but did
not say which A2020 used, and both 150K and 170K are shipped. A silent default
would be exactly the plausible-looking wrong choice rule 1 exists to prevent.
Asked in `validation/email_ts_and_h2_header_to_oscar_2026-08-11.md`.

The magnitude cut IS inferred, and the inference is recorded rather than hidden:
the upstream repo's own `GDE_Dust_corrected/Readme.txt` says that variant
"assumed 2mag cut for the dust maps", so the baseline is taken to be mag5 and
`GDE_Dust` to be mag2. If that is wrong, `CONFIGS` below is the single place to
fix it.

🚨 WHAT FERMITOOLS FORCES US TO CHANGE, AND WHAT WE LEAVE ALONE
---------------------------------------------------------------
Oscar asked (2026-08-11) that the files be used "exactly as delivered,
including their WCS information". That is not literally executable, and each
departure below is forced by a specific Fermitools error, reproduced here so
nobody has to rediscover them:

  1. CUBE -> PLANES. A <SpatialMap> cannot address a plane of a cube, and HI/H2
     arrive as one 3-D file of four annuli.
  2. PROJECTION CODE. gtsrcmaps refuses the delivered bare GLON/GLAT outright:
     "CTYPE1 must be more than 7 characters". So GLON-CAR/GLAT-CAR is set.
  3. CRVAL2 = 0. With a projection code declared and CRVAL2 = -89.75, FITS-WCS
     reads an OBLIQUE (rotated) CAR: gtsrcmaps then exited 0 having convolved
     all ten components to EXACTLY ZERO. The grid is re-expressed with
     CRVAL2 = 0 and CRPIX2 carrying the offset, which is the same sky mapping
     -- asserted, not assumed -- and is what every map in this model that
     convolves correctly already uses.
  4. CUTOUT. An all-sky <SpatialMap> throws "St12out_of_range". Every working
     one is a cutout, so we take the 87x87 GC window the public gas maps use.

Note that the public Pohl gas maps are ALREADY in exactly that processed form
-- 87x87, GLON-CAR, CRVAL2 = 0 -- so this conversion is not novel.

⚠ THE H2 HALF-PIXEL IS *NOT* IN THAT LIST, and is left as delivered.
H2_pohl_4annuli.fits declares CRVAL2 = -90.0 where the rest declare -89.75.
Our measurement says the pixels are not shifted: the emission centroid in ROW
INDEX -- which uses no header -- is 179.4243 for H2_pohl against 179.4375 for
H2_interpolated, a 0.007 deg difference against the 0.25 deg the headers claim.
Taken literally the header puts molecular gas at b = -0.29 in one map and -0.03
in its twin.

Oscar, who built and cross-checked these maps, is "not so sure about changing
the headers". He is the authority on his own data, so the delivered value is
kept and our measurement is reported to him rather than silently applied. The
consequence is visible rather than hidden: H2 is materialised half a pixel south
of HI, and the run report records it. `--fix-h2` applies our correction for a
sensitivity check.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

SRC = ROOT / "data" / "external" / "Templates_arxiv2003.10416"

#: Oscar's SECOND delivery, 2026-08-13. See `--set normalized` below and the
#: full record in `validation/oscar_delivery_2026-08-13.md`.
SRC_NORMALIZED = ROOT / "data" / "external" / "Abazajian2020_Maps"

#: A2020 GDE variant -> (gas family, dust magnitude cut). See the header note on
#: where the magnitude-cut inference comes from.
CONFIGS = {
    "baseline": ("pohl", "mag5"),
    "GDE_Gas":  ("interpolated", "mag5"),
    "GDE_Dust": ("pohl", "mag2"),
}

#: The four gas annuli, which the normalized delivery ships one file per ring.
NORMALIZED_RINGS = (2, 3, 4, 5)

#: 🚨 (config, T_S) -> filename stems in the NORMALIZED delivery.
#:
#: EVERY ENTRY HERE WAS ESTABLISHED BY CORRELATING PIXEL DATA AGAINST THE FIRST
#: DELIVERY, NOT BY READING THE FILENAMES, because the filenames mislead: the
#: UNTAGGED files are the 170 K set, and the baseline T_S = 150 K set is the one
#: carrying the `T150K` tag. Each mapping below is an exact match, rel L1 <=
#: 1.7e-14 (`validation/oscar_delivery_2026-08-13.md` section 3).
#:
#: H2 correctly has no T_S variants -- T_S is the HI spin temperature.
NORMALIZED_FILES = {
    ("baseline", "150K"): {
        "HI": [f"HI_pohl_T150K_ring_{k}_Normalized" for k in NORMALIZED_RINGS],
        "H2": [f"H2_pohl_ring_{k}_Normalized" for k in NORMALIZED_RINGS],
        "dust_pos": "DustPositive_T150K_Normalized",
        "dust_neg": "DustNegative_T150K_Normalized",
    },
    ("baseline", "170K"): {
        "HI": [f"HI_pohl_ring_{k}_Normalized" for k in NORMALIZED_RINGS],
        "H2": [f"H2_pohl_ring_{k}_Normalized" for k in NORMALIZED_RINGS],
        "dust_pos": "DustPositive_Normalized",
        "dust_neg": "DustNegative_Normalized",
    },
}

#: Why GDE_Gas and GDE_Dust are NOT in the table above, rather than being mapped
#: to the nearest-looking file. The normalized delivery does ship candidates --
#: `hi_/h2_ring_*` (lowercase), `*_galprop_*` and `*_2magcut_*` -- but NONE of
#: them is an exact match to any map in the first delivery: the lowercase gas
#: differs by up to 1.8e-3 and the galprop/2magcut dust only reaches r = 0.96-0.98.
#: They are a different generation, so which A2020 variant they correspond to is
#: NOT determined by the data we have. Guessing would be exactly the
#: plausible-looking wrong choice rule 1 exists to prevent. Ask Oscar.
NORMALIZED_UNSUPPORTED = (
    "The normalized delivery (2026-08-13) contains an exact, verified match for "
    "the BASELINE gas/dust only (pohl, mag5, at both 150K and 170K). Its "
    "candidates for {config} -- the lowercase `hi_/h2_ring_*` gas and the "
    "`_galprop_` / `_2magcut_` dust -- do NOT match any map in the first "
    "delivery exactly (r = 0.96-0.98, rel L1 up to 1.8e-3), so which A2020 "
    "variant they are is not determined by the data. Not guessed here. "
    "Use `--set delivered` for {config}, or ask Oscar which file is which. "
    "See validation/oscar_delivery_2026-08-13.md section 3."
)

#: 6-ring -> 4-ring IC grouping, CONFIRMED by Oscar 2026-08-11 (was inferred).
ICS_4RING_GROUPS = ((1, 2, 3), (4,), (5,), (6,))

#: The latitude origin every other map in the delivery uses.
COMMON_CRVAL2 = -89.75


def load_map(path: Path, fix_h2: bool = False):
    """Open a delivered map, applying the recorded header normalisations.

    Both are deliberate, guarded and reported. Neither touches pixel data.
    """
    with fits.open(path) as h:
        data = np.asarray(h[0].data, float)
        hd = h[0].header.copy()
    fixes = []

    # (1) The gas and dust maps declare bare GLON/GLAT with no projection code,
    # which `regrid_car` refuses rather than assume -- correctly, since that
    # guard is what stops an RA/DEC map being silently mishandled. These maps
    # ARE axis-aligned Galactic CAR, and that is checkable rather than assumed:
    # 720 x 360 at 0.5 deg spans l over exactly 360 deg and b over exactly
    # 180 deg linearly, which is plate carree by construction, and the IC maps
    # in the SAME delivery carry identical geometry while declaring GLON-CAR.
    if (hd.get("CTYPE1"), hd.get("CTYPE2")) == ("GLON", "GLAT"):
        span_l = abs(float(hd["CDELT1"])) * int(hd["NAXIS1"])
        span_b = abs(float(hd["CDELT2"])) * int(hd["NAXIS2"])
        assert abs(span_l - 360.0) < 1e-6 and abs(span_b - 180.0) < 1e-6, (
            f"{path.name} declares bare GLON/GLAT but spans {span_l} x {span_b} "
            f"deg, not the full sky, so it cannot be assumed to be CAR.")
        hd["CTYPE1"], hd["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
        fixes.append({"key": "CTYPE1/CTYPE2", "from": "GLON/GLAT",
                      "to": "GLON-CAR/GLAT-CAR",
                      "reason": f"missing projection code; grid spans "
                                f"{span_l:.0f}x{span_b:.0f} deg linearly = CAR"})

    # (2) The half-pixel latitude discrepancy, evidenced in the module
    # docstring. 🚨 OFF BY DEFAULT since 2026-08-11: Oscar, who built and
    # cross-checked these maps, is "not so sure about changing the headers" and
    # asked that the delivery be used as-is. He is the authority on his own
    # data, so our measurement is recorded and reported to him rather than
    # silently applied. --fix-h2 re-enables it for a sensitivity check.
    if fix_h2 and path.name == "H2_pohl_4annuli.fits":
        got = float(hd["CRVAL2"])
        assert abs(got - (-90.0)) < 1e-6, (
            f"H2_pohl CRVAL2 is {got}, not the -90.0 this override was written "
            f"for. The file has changed; re-check the offset before trusting it."
        )
        hd["CRVAL2"] = COMMON_CRVAL2
        fixes.append({"key": "CRVAL2", "from": got, "to": COMMON_CRVAL2,
                      "reason": "half-pixel latitude error; data are on the "
                                "common grid (row-centroid evidence in "
                                "the module docstring)"})
    return data, hd, fixes


#: Flat (non-solid-angle) area of a 0.5 deg pixel, in steradians. The normalized
#: delivery is normalised with THIS, not with cos(b) * dOmega -- measured, see
#: `validation/oscar_delivery_2026-08-13.md` section 4. The two differ by up to
#: 17% depending on how much of a map sits at high latitude, so they are not
#: interchangeable.
FLAT_PIXEL_SR = (0.5 * np.pi / 180.0) ** 2

#: Tolerance on that normalisation. The pohl set hits 1.000000 to ~1e-6; the
#: lowercase interpolated set drifts to 6.1e-4. Set with headroom over the
#: measured worst case, tight enough that a genuinely unnormalised map fails.
NORMALIZED_SUM_TOL = 2e-3


def load_map_normalized(paths: list[Path]):
    """Open map(s) from the NORMALIZED delivery, asserting what we used to fix.

    🚨 THIS FUNCTION DELIBERATELY REPAIRS NOTHING. Three of the four conversions
    `load_map` is forced to apply -- projection code, CRVAL2 = 0, and the signed
    negres -- are already correct in this delivery, so here they become
    PRECONDITIONS. If a future re-delivery regresses on any of them we find out
    by a loud failure instead of by a silent repair that hides the change. The
    fourth (the 87x87 cutout) is geometry gtsrcmaps forces on us regardless and
    is still applied downstream by `cutout_like_public`.

    HI and H2 arrive one file per ring. They are stacked into the same
    (4, 360, 720) array the first delivery shipped as a cube, so everything
    downstream -- `check_registration`, `materialise`, the 4-annuli assertion --
    is shared between the two paths rather than duplicated.
    """
    data, ref_hd = [], None
    for p in paths:
        with fits.open(p) as h:
            d = np.asarray(h[0].data, float)
            hd = h[0].header.copy()

        assert (hd.get("CTYPE1"), hd.get("CTYPE2")) == ("GLON-CAR", "GLAT-CAR"), (
            f"{p.name} declares {hd.get('CTYPE1')}/{hd.get('CTYPE2')}, not "
            f"GLON-CAR/GLAT-CAR. The normalized delivery is supposed to carry "
            f"the projection code; gtsrcmaps refuses a bare GLON/GLAT.")
        assert float(hd["CRVAL2"]) == 0.0, (
            f"{p.name} has CRVAL2 = {hd['CRVAL2']}, not 0. With a projection "
            f"code declared that is an OBLIQUE CAR, which gtsrcmaps convolves "
            f"to EXACTLY ZERO while exiting 0.")
        assert d.ndim == 2, f"{p.name} is {d.ndim}-D; expected one plane per file"

        # 🚨 A <SpatialMap> MAY NOT CONTAIN NEGATIVE PIXELS. Oscar states the
        # normalized set already implements the positive-flux convention, and it
        # does -- but that is checked, not taken on trust. Negative zero passes:
        # it is equal to zero, and DustNegative carries it from the sign flip.
        assert d.min() >= 0.0, (
            f"{p.name} has {(d < 0).sum()} negative pixels (min {d.min():.4e}); "
            f"the normalized delivery is supposed to be positive-valued and "
            f"gtsrcmaps rejects negatives outright.")

        got = float(d.sum() * FLAT_PIXEL_SR)
        assert abs(got - 1.0) < NORMALIZED_SUM_TOL, (
            f"{p.name} integrates to {got:.6f}, not 1, under the flat-pixel "
            f"convention this delivery uses. Either it is not from the "
            f"normalized set or the convention has changed.")

        keys = ("NAXIS1", "NAXIS2", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2",
                "CRPIX1", "CRPIX2", "CDELT1", "CDELT2")
        this = tuple(str(hd.get(k)) for k in keys)
        if ref_hd is None:
            ref_hd = (this, hd)
        assert this == ref_hd[0], (
            f"{p.name} is on a different grid from {paths[0].name}: "
            f"{dict(zip(keys, this))} vs {dict(zip(keys, ref_hd[0]))}. "
            f"Stacking these rings would misregister them.")
        data.append(d)

    stacked = np.stack(data) if len(data) > 1 else data[0]
    return stacked, ref_hd[1]


#: The window every WORKING <SpatialMap> in this model uses: 87 x 87 at 0.5 deg
#: centred on the Galactic Centre, i.e. the public Pohl gas/dust cutouts.
CUTOUT_N, CUTOUT_CRPIX = 87, 43.5


def cutout_like_public(plane, h):
    """Cut an all-sky map down to the GC window the public gas maps use.

    🚨 gtsrcmaps CANNOT TAKE AN ALL-SKY <SpatialMap> HERE. Handed the full
    720 x 360 map it threw `std::out_of_range` (2026-08-11). Every component in
    this model that convolves correctly is a cutout -- the public gas and dust
    are 87 x 87 about the GC, the bubbles 80 x 80 -- so A2020's own pipeline was
    evidently fed cutouts too. We reproduce that window exactly.

    The cut is a pure integer slice: source and target share the same 0.5 deg
    grid, so no interpolation happens and no pixel value changes. The resulting
    (l, b) of the corner pixels are asserted against the target WCS rather than
    trusted.
    """
    import numpy as _np
    cd1, cd2 = float(h["CDELT1"]), float(h["CDELT2"])
    # target: CRVAL (0,0), CRPIX (43.5, 43.5), same CDELTs as the source
    def tgt_l(i): return 0.0 + (i - CUTOUT_CRPIX) * cd1
    def tgt_b(j): return 0.0 + (j - CUTOUT_CRPIX) * cd2
    def src_i(l): return (l - float(h["CRVAL1"])) / cd1 + float(h["CRPIX1"])
    def src_j(b): return (b - float(h["CRVAL2"])) / cd2 + float(h["CRPIX2"])

    # The window is snapped to THIS map's own pixel grid rather than forced
    # onto the public one. A map whose latitude origin differs by half a pixel
    # -- H2_pohl does, and Oscar has asked us to keep it -- then keeps its own
    # registration, and the header below states honestly where the slice
    # actually is. Forcing CRPIX to the public value would have MOVED that map
    # by a quarter degree while claiming not to touch it.
    i0, j0 = int(round(src_i(tgt_l(1)))), int(round(src_j(tgt_b(1))))
    sl = plane[j0 - 1:j0 - 1 + CUTOUT_N, i0 - 1:i0 - 1 + CUTOUT_N]
    assert sl.shape == (CUTOUT_N, CUTOUT_N), (
        f"cutout ran off the map: got {sl.shape}, need {(CUTOUT_N, CUTOUT_N)}")
    # true sky position of the slice's first pixel, from the source WCS
    l0 = float(h["CRVAL1"]) + (i0 - float(h["CRPIX1"])) * cd1
    b0 = float(h["CRVAL2"]) + (j0 - float(h["CRPIX2"])) * cd2
    h = h.copy()
    h["CRVAL1"], h["CRVAL2"] = 0.0, 0.0
    h["CRPIX1"] = 1.0 - l0 / cd1
    h["CRPIX2"] = 1.0 - b0 / cd2
    h["NAXIS1"] = h["NAXIS2"] = CUTOUT_N
    # the header must describe where the pixels REALLY are
    assert abs((0.0 + (1 - h["CRPIX1"]) * cd1) - l0) < 1e-9
    assert abs((0.0 + (1 - h["CRPIX2"]) * cd2) - b0) < 1e-9
    return _np.ascontiguousarray(sl), h


def materialise(loaded: dict, out_dir: Path, files: dict | None = None,
                as_delivered: bool = False) -> dict:
    """Write the delivered maps as corrected SINGLE-PLANE FITS on disk.

    🚨 THIS EXISTS BECAUSE `gtsrcmaps` READS FITS FROM DISK. The header fixes in
    `load_map` are applied to an in-memory copy, which is enough for our own
    regridder and reaches Fermitools not at all. Two things therefore have to be
    made real on disk before the srcmaps path can use these maps:

      * the corrected headers (CTYPE projection code; H2_pohl's CRVAL2), and
      * one file per annulus -- HI and H2 arrive as a single 3-D file with four
        planes, while the diffuse XML declares each component as a 2-D
        <SpatialMap>, which cannot address a plane of a cube.

    The delivered files are never modified; corrected copies go under `runs/`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    written, sign_flipped = {}, set()
    for name, (d, hd) in loaded.items():
        if as_delivered:
            # 🚨 OSCAR'S INSTRUCTION, 2026-08-11: "keep the FITS files
            # exactly as delivered, including their WCS information, and process
            # them through the same Fermitools source-map/model pipeline used in
            # the analysis. Fermitools uses the WCS to place and convolve the
            # spatial templates, so this procedure would produce the actual maps
            # that went into our analysis."
            #
            # So the header is taken VERBATIM from the delivered file: no
            # projection code added, no CRVAL2 corrected, no cutout. The ONLY
            # unavoidable change is dropping the third-axis keywords, because a
            # <SpatialMap> cannot address a plane of a cube -- the annuli must
            # be separate files whatever else we do.
            hd = fits.getheader(files[name])
        planes = [d[k] for k in range(d.shape[0])] if d.ndim == 3 else [d]
        for k, plane in enumerate(planes):
            h = hd.copy()
            for key in ("NAXIS3", "CTYPE3", "CRVAL3", "CDELT3", "CRPIX3",
                        "CUNIT3"):
                h.remove(key, ignore_missing=True)
            h["NAXIS"] = 2
            # 🚨 THE LATITUDE WCS IS REWRITTEN INTO CAR-STANDARD FORM, and this
            # is not cosmetic -- it is the difference between a template and a
            # field of zeros. In FITS-WCS the CAR projection is the simple
            # linear grid only when the reference point lies on the native
            # equator, i.e. CRVAL2 = 0. The delivery declares CRVAL2 = -89.75
            # with CRPIX2 = 1, which is legal but means an OBLIQUE (rotated)
            # CAR. Our own regridder read it loosely and got the right answer;
            # Fermitools reads it strictly, mapped the ROI off the sky, and
            # convolved every redirected component to EXACTLY ZERO while
            # exiting 0 and writing a full set of maps (2026-08-11).
            #
            # Every map in this analysis that convolves correctly -- the public
            # gas maps, the bubbles, and the delivered IC rings -- has
            # CRVAL2 = 0. So the pixel-for-pixel identical grid is re-expressed
            # that way: b(i) is preserved exactly, only its encoding changes.
            # 🚨 A <SpatialMap> MAY NOT CONTAIN NEGATIVE PIXELS.
            # gtsrcmaps: "Image pixel value less than zero." The delivered
            # negres maps are signed -- non-positive everywhere, 78% of pixels
            # strictly negative -- while the PUBLIC release ships the same
            # component as its MAGNITUDE (min -0.0, max +22.5). So A2020's own
            # chain evidently carried the sign in the fitted amplitude, not in
            # the template, and we match that: the magnitude is written and the
            # sign convention stays where the public run already puts it.
            #
            # Guarded to the non-positive case ONLY. A mixed-sign map cannot be
            # represented as a magnitude without losing information, so it fails
            # loudly rather than being silently rectified.
            if float(plane.min()) < 0.0:
                assert float(plane.max()) <= 1e-12, (
                    f"{name} has BOTH signs (min {plane.min():.3e}, max "
                    f"{plane.max():.3e}); |x| would destroy information")
                plane = -plane
                sign_flipped.add(name)
            if as_delivered:
                written[f"{name}_ring{k + 1}" if len(planes) > 1 else name] = None
                p = out_dir / ((f"{name}_ring{k + 1}" if len(planes) > 1
                                else name) + ".fits")
                fits.PrimaryHDU(np.asarray(plane, "float32"),
                                header=h).writeto(p, overwrite=True)
                written[p.stem] = p.name
                continue
            b1 = float(h["CRVAL2"]) + (1 - float(h["CRPIX2"])) * float(h["CDELT2"])
            h["CRVAL2"] = 0.0
            h["CRPIX2"] = 1.0 - b1 / float(h["CDELT2"])
            # the rewrite must be an identity on the sky, not a shift
            b1_new = 0.0 + (1 - float(h["CRPIX2"])) * float(h["CDELT2"])
            assert abs(b1_new - b1) < 1e-9, (
                f"CAR rewrite moved row 1 from b={b1} to b={b1_new}")
            plane, h = cutout_like_public(plane, h)
            key = f"{name}_ring{k + 1}" if len(planes) > 1 else name
            p = out_dir / f"{key}.fits"
            fits.PrimaryHDU(np.asarray(plane, "float32"), header=h).writeto(
                p, overwrite=True)
            written[key] = p.name
    if sign_flipped:
        print(f"  \U0001F6A8 written as MAGNITUDE (gtsrcmaps rejects negative "
              f"pixels): {sorted(sign_flipped)}")
    return written


def row_centroid(plane: np.ndarray) -> float:
    """Emission centroid in ROW INDEX -- deliberately header-independent."""
    p = np.abs(plane).sum(axis=1)
    i = np.arange(len(p))
    m = (i > len(p) // 2 - 10) & (i < len(p) // 2 + 10)
    return float((p[m] * i[m]).sum() / p[m].sum())


#: Maps whose emission is dominated by the Galactic plane, so that "where is
#: b = 0" is a meaningful check on the header. Dust residual maps are NOT in
#: this set: they are the positive and negative parts of a gas-model residual,
#: structured rather than plane-peaked, so their latitude centroid carries no
#: registration information and gating on it would be a false alarm. (It fired
#: on the first run for exactly that reason.)
PLANE_DOMINATED = ("HI", "H2")


def check_registration(maps: dict) -> dict:
    """Every map on the same pixel grid; the gas maps also agreeing on b = 0.

    Two independent checks, because either alone can pass while the set is still
    misregistered: the WCS keywords can agree while the data are shifted, and the
    data can agree while a header lies -- which is exactly the H2_pohl case, and
    is what the second check would have caught had it been in place.
    """
    # SHAPE AND SPACING are asserted -- a map on a different pixel size or size
    # is a different product and nothing downstream would be comparable.
    #
    # THE LATITUDE ORIGIN IS REPORTED, NOT ASSERTED. H2_pohl declares a CRVAL2
    # half a pixel from the rest, and since 2026-08-11 we keep that as delivered
    # at Oscar's request. Asserting equality here would make honouring his
    # instruction impossible; what the check owes the reader is to SAY that the
    # set is not uniform, not to refuse to proceed.
    ref, rows, origins = None, {}, {}
    for name, (d, hd) in maps.items():
        key = (hd["NAXIS1"], hd["NAXIS2"], float(hd["CDELT1"]), float(hd["CDELT2"]))
        if ref is None:
            ref = key
        assert key == ref, f"{name} has a different shape/spacing: {key} vs {ref}"
        origins[name] = (float(hd["CRVAL2"])
                         + (1 - float(hd["CRPIX2"])) * float(hd["CDELT2"]))
        plane = d[0] if d.ndim == 3 else d
        rows[name] = row_centroid(plane)

    gas = {k: v for k, v in rows.items() if k in PLANE_DOMINATED}
    spread_deg = (max(gas.values()) - min(gas.values())) * abs(ref[3])
    origin_spread = max(origins.values()) - min(origins.values())
    report = {"grid": {"naxis1": ref[0], "naxis2": ref[1], "cdelt1": ref[2],
                       "cdelt2": ref[3]},
              "row_centroids": rows,
              "declared_b_of_row1": origins,
              "declared_origin_spread_deg": origin_spread,
              "gas_centroid_spread_deg": spread_deg,
              "centroid_gate_applies_to": list(PLANE_DOMINATED)}
    assert spread_deg < 0.15, (
        f"the gas maps disagree on the plane by {spread_deg:.3f} deg "
        f"({gas}) -- they are not co-registered")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", choices=sorted(CONFIGS), default="baseline")
    ap.add_argument("--set", choices=("delivered", "normalized"),
                    default="delivered",
                    help="which of Oscar's two deliveries to build from. "
                         "'delivered' = Templates_arxiv2003.10416 (2026-08-11), "
                         "the set the existing runs used. 'normalized' = "
                         "Abazajian2020_Maps (2026-08-13), the corrected set he "
                         "says actually entered A2020. Default stays "
                         "'delivered' so existing runs keep reproducing.")
    ap.add_argument("--ts", choices=("150K", "170K"), required=True,
                    help="HI spin temperature. REQUIRED: A2020's choice is not "
                         "stated in the delivery and is not guessed here.")
    ap.add_argument("--no-regrid", action="store_true",
                    help="validate and collapse on the native 0.5 deg grid "
                         "without touching the analysis grid or the counts cube")
    ap.add_argument("--fix-h2", action="store_true",
                    help="apply our half-pixel CRVAL2 correction to "
                         "H2_pohl_4annuli.fits. OFF by default: Oscar asked "
                         "that the delivered headers be kept (2026-08-11).")
    ap.add_argument("--as-delivered", action="store_true",
                    help="materialise with the DELIVERED WCS verbatim -- no "
                         "projection code, no CRVAL2 fix, no cutout (Oscar, "
                         "2026-08-11). Only the cube split is applied.")
    ap.add_argument("--materialise", action="store_true",
                    help="write corrected SINGLE-PLANE FITS under runs/ for "
                         "gtsrcmaps, which reads from disk and so cannot see "
                         "the in-memory header fixes")
    a = ap.parse_args()

    normalized = a.set == "normalized"

    # 🚨 ONE SUFFIX, USED BY EVERY OUTPUT PATH THIS SCRIPT WRITES. It previously
    # tagged only the `official_maps` cutouts, so a `--set normalized` run
    # silently overwrote the first delivery's regridded `templates_official`
    # product with a different generation of the same-named files -- no error,
    # no diff, nothing in the output saying which set was on disk. The two
    # deliveries must never share an output directory.
    set_suffix = ("_norm" if normalized else "") + ("_fixh2" if a.fix_h2 else "")

    if normalized and a.fix_h2:
        ap.error(
            "--fix-h2 is meaningless with --set normalized and is refused "
            "rather than ignored. That delivery already carries H2 on the "
            "common grid (CRVAL2 = 0, row 1 at b = -89.75) with pixel data "
            "IDENTICAL to the first delivery -- r = 1.00000000, rel L1 = 0 at "
            "zero shift on all four rings. Accepting the flag here would put a "
            "'_fixh2' tag on a run in which no correction was applied. See "
            "validation/oscar_delivery_2026-08-13.md section 2.")
    if normalized and a.as_delivered:
        ap.error("--as-delivered is a first-delivery workaround (Oscar asked "
                 "that its verbatim WCS be tried); the normalized set needs no "
                 "such treatment and the flag would only skip the cutout.")

    gas, mag = CONFIGS[a.config]
    src_dir = SRC_NORMALIZED if normalized else SRC
    if normalized:
        reg_files = NORMALIZED_FILES.get((a.config, a.ts))
        if reg_files is None:
            print("🚫 " + NORMALIZED_UNSUPPORTED.format(config=a.config),
                  file=sys.stderr)
            return 2
        files = {k: ([src_dir / f"{s}.fits" for s in v]
                     if isinstance(v, list) else src_dir / f"{v}.fits")
                 for k, v in reg_files.items()}
    else:
        files = {
            "HI": SRC / f"HI_{gas}_4annuli_TS{a.ts}.fits",
            "H2": SRC / f"H2_{gas}_4annuli.fits",
            "dust_pos": SRC / f"posres_{gas}_{mag}_{a.ts}.fits",
            "dust_neg": SRC / f"negres_{gas}_{mag}_{a.ts}.fits",
        }
    print(f"A2020 OFFICIAL TEMPLATES -- config '{a.config}' "
          f"(gas={gas}, dust={mag}, T_S={a.ts}) "
          f"[{'NORMALIZED set, 2026-08-13' if normalized else 'first delivery, 2026-08-11'}]")
    flat = [p for v in files.values()
            for p in (v if isinstance(v, list) else [v])]
    missing = [str(p) for p in flat if not p.exists()]
    if missing:
        print("🚫 missing delivered maps:\n  " + "\n  ".join(missing),
              file=sys.stderr)
        return 2

    loaded, overrides = {}, {}
    for name, p in files.items():
        if normalized:
            ps = p if isinstance(p, list) else [p]
            d, hd = load_map_normalized(ps)
            shown = (f"{ps[0].name} (+{len(ps) - 1} rings)" if len(ps) > 1
                     else ps[0].name)
        else:
            d, hd, fixes = load_map(p, fix_h2=a.fix_h2)
            shown = p.name
            if fixes:
                overrides[p.name] = fixes
                for f in fixes:
                    print(f"  🚨 header normalised on {p.name}: "
                          f"{f['key']} {f['from']} -> {f['to']}")
        loaded[name] = (d, hd)
        print(f"  {name:<9} {shown:<44} {str(d.shape):>16}")
    if normalized:
        print("  ✅ preconditions hold on all inputs: GLON-CAR/GLAT-CAR, "
              "CRVAL2 = 0, no negative pixels, flat-pixel integral = 1")

    reg = check_registration(loaded)
    print(f"\n  co-registration OK: all four on {reg['grid']['naxis1']}x"
          f"{reg['grid']['naxis2']} at {abs(reg['grid']['cdelt1'])} deg; "
          f"HI and H2 agree on the plane to "
          f"{reg['gas_centroid_spread_deg']:.3f} deg")
    if reg["declared_origin_spread_deg"] > 1e-9:
        print(f"  \u26a0 the delivery is NOT uniform in declared latitude origin "
              f"(spread {reg['declared_origin_spread_deg']:.2f} deg):")
        for k, v in reg["declared_b_of_row1"].items():
            print(f"      {k:<9} header puts row 1 at b = {v:+.3f}")
        print(f"      kept as delivered per Oscar 2026-08-11; --fix-h2 overrides")
    print(f"  (dust residual centroids reported, not gated: "
          + ", ".join(f"{k} {v:.2f}" for k, v in reg['row_centroids'].items()
                      if k not in PLANE_DOMINATED) + " in row index)")

    if a.materialise:
        from a2020_stage2_reduce import a2020_config
        mdir = (ROOT / "runs" / a2020_config().hash / "official_maps"
                / f"{a.config}_{a.ts}{set_suffix}")
        written = materialise(loaded, mdir, files=files,
                              as_delivered=a.as_delivered)
        print(f"\n  materialised {len(written)} single-plane FITS -> "
              f"{mdir.relative_to(ROOT)}")
        for k in sorted(written):
            print(f"    {k:<12} {written[k]}")
        # 🚨 VERIFIED ON DISK against the geometry gtsrcmaps actually
        # accepts, not against what we intended. Every invariant here is one a
        # failed run taught us: the projection code (regrid_car refuses without
        # it), CRVAL2 = 0 (an oblique CAR convolved to all zeros), 2-D (a
        # SpatialMap cannot index a cube), and the 87x87 GC window (an all-sky
        # SpatialMap threw std::out_of_range).
        chk = fits.getheader(mdir / "H2_ring1.fits")
        if a.as_delivered:
            raw = fits.getheader(files["H2"])
            for k in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1",
                      "CRPIX2", "CDELT1", "CDELT2", "NAXIS1", "NAXIS2"):
                assert str(chk.get(k)) == str(raw.get(k)), (
                    f"as-delivered materialise altered {k}: "
                    f"{raw.get(k)} -> {chk.get(k)}")
            assert chk["NAXIS"] == 2, "materialised map is still a cube"
            print(f"  verified on disk: WCS is VERBATIM from the delivery "
                  f"({chk['CTYPE1']}/{chk['CTYPE2']}, CRVAL2={chk['CRVAL2']}, "
                  f"{chk['NAXIS1']}x{chk['NAXIS2']}); only the cube split applied")
            return 0
        assert chk["CTYPE1"] == "GLON-CAR" and chk["CTYPE2"] == "GLAT-CAR", \
            "materialised map lost its projection code on the way to disk"
        assert float(chk["CRVAL2"]) == 0.0, (
            f"CRVAL2 is {chk['CRVAL2']}, not 0 -- that is an OBLIQUE CAR and "
            f"gtsrcmaps will convolve it to zero")
        assert chk["NAXIS"] == 2, "materialised map is still a cube"
        assert float(chk["CRPIX2"]) == CUTOUT_CRPIX or a.fix_h2 is False, \
            "unexpected H2 registration"
        assert chk["NAXIS1"] == chk["NAXIS2"] == CUTOUT_N, (
            f"materialised map is {chk['NAXIS1']}x{chk['NAXIS2']}, not the "
            f"{CUTOUT_N}x{CUTOUT_N} GC window gtsrcmaps accepts")
        print(f"  verified on disk: {chk['CTYPE1']}/{chk['CTYPE2']}, "
              f"CRVAL2={chk['CRVAL2']}, NAXIS={chk['NAXIS']}, "
              f"{chk['NAXIS1']}x{chk['NAXIS2']}")

    for name in ("HI", "H2"):
        n = loaded[name][0].shape[0]
        assert n == 4, f"{name} has {n} annuli, expected 4"
    print(f"  HI and H2 both carry 4 annuli, as A2020's Table II footnote says")

    report = {"config": a.config, "gas_family": gas, "dust_mag_cut": mag,
              "spin_temperature": a.ts,
              "template_set": a.set,
              "source_dir": str(src_dir.relative_to(ROOT)),
              "provenance": ("Oscar Macias, email 2026-08-13 (normalized set)"
                             if normalized else
                             "Oscar Macias, email 2026-08-11"),
              "ics_4ring_groups": [list(g) for g in ICS_4RING_GROUPS],
              "ics_grouping_status": "CONFIRMED by Oscar 2026-08-11",
              "header_overrides": overrides,
              "registration": reg,
              "components_not_covered": ["bubbles", "loopI", "boxy_bulge",
                                         "nuclear_bulge", "isotropic",
                                         "point_sources", "sun_and_moon"],
              "files": {k: ([p.name for p in v] if isinstance(v, list)
                            else v.name) for k, v in files.items()}}

    if a.no_regrid:
        out = (ROOT / "runs"
               / f"a2020_official_{a.config}_{a.ts}{set_suffix}_native.json")
        out.write_text(json.dumps(report, indent=2))
        print(f"\n  --no-regrid: validation only, wrote {out}")
        return 0

    from a2020_stage2_reduce import a2020_config
    from a2020_stage2_diffuse import ics_on_bins, target_header
    from gcereduce.regrid import regrid_car

    cfg = a2020_config()
    ccube = ROOT / "runs" / cfg.hash / "ccube.fits"
    if not ccube.exists():
        print(f"🚫 no counts cube at {ccube}", file=sys.stderr)
        return 2
    tgt = target_header(ccube)
    edges = np.array(cfg.energy_edges_gev, float)
    out_dir = (ROOT / "runs" / cfg.hash / "templates_official"
               / f"{a.config}_{a.ts}{set_suffix}")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n  target : {cfg.npix} x {cfg.npix} at {cfg.pixel_deg} deg, "
          f"{cfg.n_ebins} bins")
    print(f"  out    : {out_dir.relative_to(ROOT)}\n")

    t0 = time.time()
    frac = {}
    for name, (d, hd) in loaded.items():
        planes = [d[k] for k in range(d.shape[0])] if d.ndim == 3 else [d]
        for k, plane in enumerate(planes):
            g = regrid_car(plane, hd, tgt)
            key = f"{name}_ring{k + 1}" if len(planes) > 1 else name
            np.save(out_dir / f"tmpl_{key}.npy", g.astype("float32"))
            # regrid_car conserves FLUX, not the sum of pixel values; divide by
            # the pixel-area ratio for the physical in-ROI fraction. Same
            # correction, and same reason, as in a2020_stage2_diffuse.py.
            ar = (abs(float(hd["CDELT1"])) / cfg.pixel_deg) ** 2
            f = float(g.sum() / plane.sum() / ar) if plane.sum() else float("nan")
            frac[key] = f
            print(f"  {key:<12} -> {g.shape}   flux in ROI {f:6.3f}")

    print()
    # 🚨 The (1+2+3) collapse is only valid if the ring maps carry PHYSICAL units.
    # Had they been individually normalised -- as the public release's gas maps
    # are, every one summing to the same value -- adding them would have silently
    # misweighted the inner Galaxy. Checked here rather than assumed.
    ring_totals = [float(np.asarray(fits.getdata(SRC / f"ics_Ring_{i}.fits"),
                                    float).sum()) for i in range(1, 7)]
    spread = max(ring_totals) / min(ring_totals)
    assert spread > 1.05, (
        f"IC ring totals agree to {spread:.4f}x -- they look individually "
        f"normalised, and summing (1+2+3) would misweight the inner Galaxy")
    print(f"  IC rings carry physical units (totals span {spread:.0f}x), "
          f"so the (1+2+3) collapse is valid")
    report["ics_ring_totals"] = ring_totals
    report["ics_total_spread"] = spread

    rings6 = ics_on_bins(edges, tgt, src_dir=SRC)
    for j, grp in enumerate(ICS_4RING_GROUPS, start=1):
        stack = sum(rings6[i - 1] for i in grp)
        np.save(out_dir / f"tmpl_ics_ring{j}.npy", stack.astype("float32"))
        print(f"  ics_ring{j} <- 6-ring {list(grp)}   {stack.shape}  "
              f"sum {stack.sum():.4e}")

    report["in_roi_fraction"] = frac
    report["seconds"] = time.time() - t0
    report["out_dir"] = str(out_dir.relative_to(ROOT))
    (out_dir / "report.json").write_text(json.dumps(report, indent=2))
    n_written = len(frac) + len(ICS_4RING_GROUPS)
    report["n_templates"] = n_written
    print(f"\n  wrote {n_written} templates + report.json in "
          f"{report['seconds']:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
