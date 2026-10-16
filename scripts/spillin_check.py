#!/usr/bin/env python
"""Emission entering the |l|,|b| <= 60 grid from beyond it (Oscar, 2026-10-01,
item 2): refold every template on a grid padded to +-62 deg and measure what
changes in the fitted cells and in the halo fit.

    FINE_TAG=psfuc python scripts/spillin_check.py --cell 2
    FINE_TAG=psfuc python scripts/spillin_check.py --cell 1

The production libraries fold the PSF on the 960x960 grid with the edge value
replicated beyond it ('nearest'). Here the pixel maps are rebuilt on a 992x992
grid (+-62 deg; `runs/spillin/expcube_pad992.fits` is gtexpcube2 on that grid
with the same livetime cube and IRF, identical to the production cube on its
central crop), folded there, and binned to the same cells; pixels beyond
+-60 deg drop out of the binning as always. The diffuse cubes (all-sky
HEALPix), the catalog sources (all-sky), isotropic, the bubble maps, Loop I
and the Haslam map are evaluated on the padded grid; a template whose source
map does not extend beyond +-60 deg keeps the edge replication and is listed.
Totani's construction is built from counts and is unchanged.

Reported per template and bin: the change of the summed expected counts in
the fitted cells, and the largest relative change of any cell in the outer
ring (cells touching +-60 deg) and in the interior. Then the three fixed
configurations of oscar_conditional.py are refitted with the padded library
for the representative rows: halo at 20.8 GeV, Q+ (K = 6 and 9), and the
flat - structured difference with the halo in, against the production values.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

DATE = "2026-10-01"
PAD_DEG = 2.0
PAD_EXPCUBE = ROOT / "runs" / "spillin" / "expcube_pad992.fits"
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1")]


def padded_grid(edges_gev):
    from astropy.io import fits
    from totani.exposure import flat_spectrum_bin_exposure, single_energy_bin_exposure
    from totani.pointsources import default_catalog_path, extended_source_mask
    from totani.stack import _strip_to_2d
    from totani.templates import grid_lb, pixel_solid_angle, target_header
    npix = 960 + int(round(2 * PAD_DEG / 0.125))
    target = target_header(npix)
    with fits.open(PAD_EXPCUBE) as h:
        cube = h[0].data.astype(float); emev = h["ENERGIES"].data["Energy"].astype(float); hdr = _strip_to_2d(h[0].header)
    assert int(hdr["NAXIS1"]) == npix and abs(float(hdr["CRPIX1"]) - (npix + 1) / 2) < 1e-6
    L, B = grid_lb(target); om = pixel_solid_angle(target)
    with fits.open(default_catalog_path()) as h:
        mask = extended_source_mask(L, B, h["ExtendedSources"].data)
    return {"target": target, "L": L, "B": B, "mask": mask, "exposure_cube": cube, "exposure_mev": emev, "exp_header": hdr,
            "basis": flat_spectrum_bin_exposure(cube, emev, edges_gev) * om[None],
            "single_sr": single_energy_bin_exposure(cube, emev, edges_gev) * om[None]}


def fold_and_bin(cube_pix, g, cell, edges_gev):
    from totani import psf
    from totani.cells import bin_to_cells
    conv = psf.convolve_cube(np.asarray(cube_pix, float), g["B"], edges_gev)
    return np.stack([bin_to_cells(conv[k], g["L"], g["B"], cell_deg=cell, mask=g["mask"]) for k in range(conv.shape[0])])


def padded_templates(g, cell, edges_gev, names):
    """The model-independent templates and J on the padded grid; returns (dict, notes)."""
    from astropy.io import fits
    from totani import structures
    from totani.halo import nfw_rho2_on_target
    from totani.pointsources import default_catalog_path, exposure_pixel, extended_source_mask, source_expected_counts
    from totani.templates import bubble_macias_on_target, bubble_on_target, isotropic_on_target, loopi_on_target
    target = g["target"]; T = {}; notes = {}
    shapes = {"bubble_flat": lambda: (bubble_on_target(target) > 0).astype(float), "bubble_structured": lambda: bubble_macias_on_target(target),
              "isotropic": lambda: isotropic_on_target(target)}
    loops = {}
    for n in names:
        t0 = time.time()
        try:
            if n in shapes:
                cube = shapes[n]()[None] * g["basis"]
            elif n.startswith("loopI_shell"):
                if not loops:
                    loops = loopi_on_target(target)
                cube = loops[n][None] * g["basis"]
            elif n == "point_sources":
                with fits.open(default_catalog_path()) as h:
                    cat = h["LAT_Point_Source_Catalog"].data; ext = h["ExtendedSources"].data
                cl = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"]); cb = cat["GLAT"]
                half = 60.0 + PAD_DEG
                keep = (np.abs(cl) <= half) & (np.abs(cb) <= half) & ~extended_source_mask(cl, cb, ext)
                sc = source_expected_counts(cat[keep], cl[keep], cb[keep], edges_gev, g["exposure_cube"], g["exposure_mev"], g["exp_header"])
                # the production placement (fine.point_sources_cells, exact=True): exact positions, kernel integrated over pixels
                from totani import psf as _psf
                from totani.cells import bin_to_cells as _b2c
                conv = _psf.point_sources_exact(cl[keep], cb[keep], sc, g["L"], g["B"], edges_gev)
                T[n] = np.stack([_b2c(conv[k], g["L"], g["B"], cell_deg=cell, mask=g["mask"]) for k in range(conv.shape[0])])
                notes[n] = f"{int(keep.sum())} sources within +-{half:g} deg, {int(keep.sum()) - 2411} of them in the padding; exact-position placement as in production"
                print(f"  padded {n} ({time.time()-t0:.0f}s)", flush=True)
                continue
            elif n == "J":
                cube = nfw_rho2_on_target(target)[None] * g["basis"]
            else:
                arr = structures.STRUCTURES[n](target, edges_gev)
                cube = arr[None] * g["basis"] if arr.ndim == 2 else arr * g["single_sr"]
            if not np.isfinite(cube).all():
                raise ValueError("non-finite values on the padded grid")
            T[n] = fold_and_bin(cube, g, cell, edges_gev)
            print(f"  padded {n} ({time.time()-t0:.0f}s)", flush=True)
        except Exception as e:  # noqa: BLE001 -- the source map does not cover the padding
            notes[n] = f"kept edge replication: {type(e).__name__}: {str(e)[:80]}"
            print(f"  {n}: {notes[n]}", flush=True)
    return T, notes


def compare(name, new, old, keep, cell):
    """Per-bin change: total over fitted cells, max relative change in the outer ring and interior."""
    n = old.shape[-1]
    idx = np.arange(n)
    outer = (idx[:, None] == 0) | (idx[:, None] == n - 1) | (idx[None, :] == 0) | (idx[None, :] == n - 1)
    out = []
    for k in range(old.shape[0]):
        o, w = old[k], new[k]; d = w - o
        ref = np.maximum(np.abs(o), 1e-3 * np.abs(o[keep]).max())
        rel = np.where(keep, np.abs(d) / ref, 0.0)
        tot = o[keep].sum()
        out.append({"total_rel": float(d[keep].sum() / tot) if tot else 0.0,
                    "max_rel_outer": float(rel[outer & keep].max()), "max_rel_inner": float(rel[~outer & keep].max()),
                    "abs_change_outer": float(np.abs(d)[outer & keep].sum()), "abs_change_inner": float(np.abs(d)[~outer & keep].sum())})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, default=2.0); a = ap.parse_args()
    assert os.environ.get("FINE_TAG") == "psfuc", "run with FINE_TAG=psfuc"
    from reduce_totani import TOTANI_EDGES_GEV
    from totani import fine, ranking
    import oscar_conditional as OC
    edges = TOTANI_EDGES_GEV[:10]
    cell = a.cell
    lib = fine.load_library(cell); keep = lib["keep"]
    g = padded_grid(edges)
    names = list(lib["templates"]) + ["J"]
    t0 = time.time()
    T, notes = padded_templates(g, cell, edges, names)
    res = {"cell": cell, "pad_deg": PAD_DEG, "notes": notes, "templates": {}, "rows": {}}
    for n in names:
        if n in T:
            old = lib["J"] if n == "J" else lib["templates"][n]
            res["templates"][n] = compare(n, T[n], old, keep, cell)
    # the rows' diffuse cubes on the padded grid, then the fixed configurations refitted
    lib_pad = {**lib, "templates": {**lib["templates"], **{n: T[n] for n in T if n != "J"}}, "J": T.get("J", lib["J"])}
    OC._W["cell"] = cell
    for code, variant in ROWS:
        t1 = time.time()
        cubes = ranking.galprop_cubes(code, variant, g["target"], edges)
        diffuse = {f"galprop_{n}": fold_and_bin(c * g["single_sr"], g, cell, edges) for n, c in cubes.items()}
        row = fine.load_row(cell, code, variant)
        for n in diffuse:
            res["templates"][f"{code} {variant} {n}"] = compare(n, diffuse[n], row["diffuse"][n], keep, cell)
        print(f"  {code} {variant} diffuse on the padded grid ({time.time()-t1:.0f}s)", flush=True)
        # monkeypatch load_row so OC.stack_for sees the padded diffuse templates
        orig = fine.load_row
        fine.load_row = lambda c, *k, _r={"diffuse": diffuse, "residual": row["residual"]}: _r
        out = {}
        OC._W["lib"] = lib_pad
        for config in OC.CONFIGS:
            out[(config, "padded")] = OC.evaluate(lib_pad["counts"], (code, variant), config, with_intervals=False)
        fine.load_row = orig
        prod = json.loads((ROOT / "results" / "fits" / f"oscar_conditional_cells{cell:g}_psfuc_2026-09-29.json").read_text())["results"]
        for config in OC.CONFIGS:
            p = next(r for r in prod if tuple(r["row"]) == (code, variant) and r["config"] == config)
            out[(config, "production")] = p
        rr = {}
        for config in OC.CONFIGS:
            n_, p_ = out[(config, "padded")], out[(config, "production")]
            rr[config] = {"halo20": [n_["counts"][5], p_["counts"][5]], "Q6": [n_["resp"]["Q_plus"], p_["resp"]["Q_plus"]], "Q9": [n_["all"]["Q_plus"], p_["all"]["Q_plus"]],
                          "Z6": [n_["resp"]["Z_2018"], p_["resp"]["Z_2018"]], "v1": [n_["v1"], p_["v1"]], "v0": [n_["v0"], p_["v0"]], "counts": [n_["counts"], p_["counts"]]}
        res["rows"][f"{code} {variant}"] = rr
        print(f"  {code} {variant} refitted ({time.time()-t1:.0f}s)", flush=True)
    (ROOT / "results" / "fits" / f"spillin_check_cells{cell:g}_psfuc_{DATE}.json").write_text(json.dumps(res, indent=1, default=float))
    render(res)
    print(f"done ({time.time()-t0:.0f}s)")
    return 0


def render(res):
    cell = res["cell"]; E = [1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0]
    L = [f"# Emission from beyond ±60°: templates refolded on a grid padded to ±{60 + res['pad_deg']:g}°, {cell:g}° cells, PSF-folded ({DATE})", "",
         "Production libraries replicate the edge value beyond ±60° when folding the PSF. Here every template is rebuilt on the padded grid "
         "(exposure from gtexpcube2 on that grid, same livetime cube and IRF; identical to the production cube on its central crop), folded there, and binned to the same cells. "
         "'total' = change of the summed expected counts over the fitted cells; 'outer' / 'inner' = largest relative change of any fitted cell in the ring of cells touching ±60° / elsewhere.", ""]
    if res["notes"]:
        L += ["Notes:"] + [f"- {k}: {v}" for k, v in res["notes"].items()] + [""]
    L += ["## Templates: change in the fitted cells (per bin, 1.5 → 100 GeV)", "", "| template | quantity | " + " | ".join(f"{e:g}" for e in E) + " |", "|---|---|" + "---|" * 9]
    for n, rows in res["templates"].items():
        L.append(f"| {n} | total | " + " | ".join(f"{r['total_rel']*100:+.3f}%" for r in rows) + " |")
        L.append(f"| | max outer | " + " | ".join(f"{r['max_rel_outer']*100:.2f}%" for r in rows) + " |")
        L.append(f"| | max inner | " + " | ".join(f"{r['max_rel_inner']*100:.2f}%" for r in rows) + " |")
    L += ["", "## Fixed configurations refitted with the padded library (padded / production)", "",
          "| row | bubbles | halo 20.8 GeV | Q⁺ K=6 | Z₂₀₁₈ K=6 | Q⁺ K=9 |", "|---|---|---|---|---|---|"]
    RESP = [3, 4, 5, 6, 7, 8]
    for row, rr in res["rows"].items():
        for c, x in rr.items():
            L.append(f"| {row} | {c} | {x['halo20'][0]:,.0f} / {x['halo20'][1]:,.0f} | {x['Q6'][0]:.1f} / {x['Q6'][1]:.1f} | {x['Z6'][0]:.1f} / {x['Z6'][1]:.1f} | {x['Q9'][0]:.1f} / {x['Q9'][1]:.1f} |")
    L += ["", "| row | flat − structured, H1, K=6 (padded / production) | all bins | construction − flat, H1, K=6 | all bins |", "|---|---|---|---|---|"]
    for row, rr in res["rows"].items():
        d = lambda a, b, i, bins: sum(rr[a]["v1"][i][k] for k in bins) - sum(rr[b]["v1"][i][k] for k in bins)
        L.append(f"| {row} | {d('flat','structured',0,RESP):+.0f} / {d('flat','structured',1,RESP):+.0f} | {d('flat','structured',0,range(9)):+.0f} / {d('flat','structured',1,range(9)):+.0f} | "
                 f"{d('construction','flat',0,RESP):+.0f} / {d('construction','flat',1,RESP):+.0f} | {d('construction','flat',0,range(9)):+.0f} / {d('construction','flat',1,range(9)):+.0f} |")
    L += ["", "| row | bubbles | halo counts per bin, padded − production |", "|---|---|---|"]
    for row, rr in res["rows"].items():
        for c, x in rr.items():
            L.append(f"| {row} | {c} | " + " ".join(f"{a-b:+.0f}" for a, b in zip(x["counts"][0], x["counts"][1])) + " |")
    out = ROOT / "results" / f"spillin_check_cells{cell:g}_psfuc_{DATE}.md"
    out.write_text("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    raise SystemExit(main())
