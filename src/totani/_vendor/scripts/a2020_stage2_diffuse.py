#!/usr/bin/env python
"""A2020 redo, STAGE 2b-ii: their diffuse templates on our analysis grid.

    .venv/bin/python scripts/a2020_stage2_diffuse.py

Regrids the public Pohl et al. (2022) maps -- which ARE A2020's Table I template
library -- onto the 400x400, 0.1 deg Galactic CAR grid of the Stage 2a reduction,
and puts the inverse-Compton cube on A2020's 15 energy bins.

WHAT THIS DOES AND DELIBERATELY DOES NOT DO
--------------------------------------------
It produces **regridded, energy-binned spatial templates**. It does NOT apply the
PSF, and that omission is the point:

🚨 `gcereduce.forward.king_kernel` IS NOT THE LAT PSF. Its own module docstring
says so: the King parameters were fitted to map 0.25 deg Cholis maps onto gcepy's
0.1 deg counts, so they absorb the PSF *and* a 0.25 deg pixel boxcar *and* the
in-bin energy integration. "The Pohl maps are 0.5 deg, so the boxcar term is
wrong for them... For anything where the absolute kernel matters, use gtsrcmaps
with the real IRF. Recorded so nobody promotes this to 'the PSF'."

Reproducing Table II's TS is exactly a case where the absolute kernel matters, so
the convolution is left to `gtsrcmaps` with A2020's own
`P8R3_ULTRACLEANVETO_V2` -- the same path the point sources take.

THE FOUR-RING ICS IS BUILT BY SUMMING, WHICH IS EXACT
-----------------------------------------------------
A2020's baseline uses a 3D IC map in FOUR rings (0-3.5, 3.5-8, 8-10, 10-50 kpc),
matching the gas rings. Pohl ships SIX (0-1.5, 1.5-2.5, 2.5-3.5, 3.5-8, 8-10,
10-50), so

    4-ring #1 = 6-ring #1 + #2 + #3      (0-1.5, 1.5-2.5, 2.5-3.5 -> 0-3.5)
    4-ring #2,3,4 = 6-ring #4,5,6        (identical intervals)

Summing is valid ONLY because the ICS maps are in physical units rather than
individually normalised -- verified: the six ring sums span 1.7e-5 to 7.7e-4, so
they are not each scaled to unity the way the gas maps are (`*_Normalized`,
every one summing to 13131.2254). Had they been normalised, adding them would
have silently weighted the inner Galaxy wrongly.

⚠ SUN AND MOON ARE NOT REGRIDDED. They are the only Table I components in
RA---CAR/DEC--CAR rather than Galactic CAR, and `regrid_car` refuses non-Galactic
grids by design rather than silently mis-rotating them. They are a real part of
A2020's baseline and their absence is reported, not hidden.

⚠ THE ICS ENERGY REBINNING ASSUMES the third axis samples a density in log10 E
(CTYPE3 'photon energy', CRVAL3 = log10(500 MeV), CDELT3 = 0.0969 dex). No BUNIT
is present to confirm it. The templates carry a FREE NORMALISATION PER BIN, so
this affects only the in-bin weighting, not the overall scale -- but it is an
assumption, and it is recorded rather than buried.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

POHL = ROOT / "data" / "external" / "Pohl_etal_2022"

#: Table I components that are plain 2-D spatial maps, fitted with one free
#: normalisation per energy bin ("the spectra needs to be determined with a
#: bin-by-bin likelihood analysis" -- Pohl Readme).
SPATIAL = {
    "HI_ring1": "HI_pohl_Ts_bestTexc_ring_1_Normalized.fits",
    "HI_ring2": "HI_pohl_Ts_bestTexc_ring_2_Normalized.fits",
    "HI_ring3": "HI_pohl_Ts_bestTexc_ring_3_Normalized.fits",
    "HI_ring4": "HI_pohl_Ts_bestTexc_ring_4_Normalized.fits",
    "H2_ring1": "H2_pohl_ring_1_Normalized.fits",
    "H2_ring2": "H2_pohl_ring_2_Normalized.fits",
    "H2_ring3": "H2_pohl_ring_3_Normalized.fits",
    "H2_ring4": "H2_pohl_ring_4_Normalized.fits",
    "dust_pos": "Dust_Positive_Residuals_Normalized.fits",
    "dust_neg": "Dust_Negative_Residuals_Normalized.fits",
    "bubbles": "Fermi_Bubbles_Normalized.fits",
    "loopI": "LoopI_Normalized.fits",
    "boxy_bulge": "Boxy_Bulge_Normalized.fits",
    "nuclear_bulge": "Nuclear_Bulge_Normalized.fits",
}

#: 6-ring -> 4-ring grouping. A2020's baseline is the four-ring IC map.
ICS_4RING_GROUPS = ((1, 2, 3), (4,), (5,), (6,))


def target_header(ccube: Path):
    with fits.open(ccube) as h:
        hd = h[0].header.copy()
    for k in ("NAXIS3", "CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3"):
        hd.remove(k, ignore_missing=True)
    hd["NAXIS"] = 2
    return hd


def ics_on_bins(edges_gev, tgt, src_dir=None):
    """The six ICS rings regridded and integrated onto the analysis bins.

    `src_dir` selects which copy of the ring maps to read; it defaults to the
    public Pohl release. `a2020_templates_official.py` passes A2020's own
    delivered maps instead, and reuses this function rather than restating the
    log10-E integration, which is the fiddly part.
    """
    from gcereduce.regrid import regrid_car
    src_dir = POHL if src_dir is None else Path(src_dir)
    out = []
    for i in range(1, 7):
        with fits.open(src_dir / f"ics_Ring_{i}.fits") as h:
            cube, hd = h[0].data.astype(float), h[0].header
            e = np.array([r[0] for r in h["ENERGIES"].data], float)   # MeV
        lg = np.log10(e)
        planes = np.stack([regrid_car(cube[k], hd, tgt) for k in range(cube.shape[0])])
        binned = np.empty((len(edges_gev) - 1,) + planes.shape[1:])
        for b in range(len(edges_gev) - 1):
            lo, hi = np.log10(edges_gev[b] * 1e3), np.log10(edges_gev[b + 1] * 1e3)
            # integrate the log10-E density over the bin, on a fine sub-grid so
            # bins narrower than the 0.0969 dex slice spacing are still sampled
            xs = np.linspace(lo, hi, 25)
            # Interpolate the whole map at each sub-energy and accumulate, rather
            # than looping pixels: 25 map-sized operations per bin instead of
            # 160000 scalar interpolations.
            acc = np.zeros(planes.shape[1:])
            w = np.gradient(xs)
            for x, wi in zip(xs, w):
                k = np.clip(np.searchsorted(lg, x) - 1, 0, len(lg) - 2)
                f = (x - lg[k]) / (lg[k + 1] - lg[k])
                acc += wi * (planes[k] * (1 - f) + planes[k + 1] * f)
            binned[b] = acc
        out.append(binned)
    return out


def main() -> int:
    from a2020_stage2_reduce import a2020_config
    from gcereduce.regrid import regrid_car

    cfg = a2020_config()
    workdir = ROOT / "runs" / cfg.hash
    ccube = workdir / "ccube.fits"
    if not ccube.exists():
        print(f"🚫 no counts cube at {ccube}", file=sys.stderr)
        return 2
    tgt = target_header(ccube)
    edges = np.array(cfg.energy_edges_gev, float)
    out_dir = workdir / "templates"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("A2020 STAGE 2b-ii -- diffuse templates onto the analysis grid")
    print(f"  target : {cfg.npix} x {cfg.npix} at {cfg.pixel_deg} deg, "
          f"{cfg.n_ebins} bins\n")

    t0 = time.time()
    report = {}
    for name, fn in SPATIAL.items():
        with fits.open(POHL / fn) as h:
            d, hd = h[0].data.astype(float), h[0].header
        g = regrid_car(d, hd, tgt)
        np.save(out_dir / f"tmpl_{name}.npy", g.astype("float32"))
        # 🚨 regrid_car conserves FLUX (sum of V*w), not the sum of pixel VALUES.
        # Going 0.5 deg -> 0.1 deg there are 25x more pixels each holding a
        # similar intensity, so a naive sum ratio comes out ~25 and looks like a
        # fraction above 1. Divide by the pixel-area ratio to get the physical
        # fraction of all-sky flux landing in the ROI. Checked against the table
        # recorded in CONVENTIONS.md: 0.952 / 0.950 / 0.041 / 1.000 / 0.997,
        # reproduced exactly.
        area_ratio = (abs(float(hd["CDELT1"])) / cfg.pixel_deg) ** 2
        frac = float(g.sum() / d.sum() / area_ratio) if d.sum() else float("nan")
        report[name] = {"in_roi_fraction": frac, "shape": list(g.shape),
                        "native_pixel_deg": abs(float(hd["CDELT1"]))}
        print(f"  {name:<16} {str(d.shape):>14} -> {g.shape}   "
              f"flux in ROI {frac:6.3f}")

    print()
    rings6 = ics_on_bins(edges, tgt)
    for j, grp in enumerate(ICS_4RING_GROUPS, start=1):
        stack = sum(rings6[i - 1] for i in grp)
        np.save(out_dir / f"tmpl_ics_ring{j}.npy", stack.astype("float32"))
        report[f"ics_ring{j}"] = {"from_6ring": list(grp),
                                  "shape": list(stack.shape)}
        print(f"  ics_ring{j} <- 6-ring {list(grp)}   {stack.shape}  "
              f"sum {stack.sum():.4e}")

    dt = time.time() - t0
    report["_meta"] = {
        "seconds": dt, "n_templates": len(SPATIAL) + 4,
        "psf_applied": False,
        "psf_note": "left to gtsrcmaps with P8R3_ULTRACLEANVETO_V2; "
                    "forward.king_kernel is NOT the LAT PSF",
        "missing_from_a2020_baseline": ["Sun_and_Moon (RA/DEC-CAR, not regriddable "
                                        "by regrid_car)", "isotropic (added at fit "
                                        "time from iso_P8R3_ULTRACLEANVETO_V2_v1.txt)"],
    }
    (ROOT / "runs" / "a2020_stage2_diffuse.json").write_text(
        json.dumps(report, indent=1))
    print(f"\n  {len(SPATIAL) + 4} templates in {dt:.1f} s -> {out_dir}")
    print("  ⚠ PSF NOT applied here (see module docstring); Sun+Moon not regridded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
