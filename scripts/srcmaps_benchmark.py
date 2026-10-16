#!/usr/bin/env python
"""Benchmark of this pipeline's PSF folding against standard Fermitools source
maps (Oscar, 2026-10-01, item 2): gtsrcmaps in 20x20 deg sub-regions for a few
4FGL point sources, the structured bubble map and the NFW J-map, compared cell
by cell with the folded templates at 2 and 1 deg.

    python scripts/srcmaps_benchmark.py prepare      # gcepy env: cropped ccubes, maps, XML, shell script
    bash runs/srcmaps/run_fermitools.sh               # fermi-2.5.3 env (activated inside): gtexpcube2 + gtsrcmaps
    FINE_TAG=psfuc python scripts/srcmaps_benchmark.py compare   # gcepy env

Sub-regions (centers on the 2-deg cell lattice so cells coincide with the
production cells): (l, b) = (0, +16), (0, +40), (-40, -40). Each sub-region
is a 160x160 crop of the production counts cube (the same 0.125 deg plate
carree grid, 13 energy planes), so gtsrcmaps sees the data's own binning;
the binned exposure is gtexpcube2 on a 30x30 deg grid around the center with
the production livetime cube and IRF (P8R3_ULTRACLEAN_V3, evtype 3). The
spatial maps are written on the +-62 deg padded grid, in intensity per sr.

Fermitools convolves each energy plane with the PSF at that energy (point
sources: at the source position; diffuse maps: the mean PSF over the region)
and the likelihood integrates across the bin; this pipeline averages the PSF
over the bin (E^-2.4 x exposure), fits it with three Gaussians and folds the
bin-integrated template separably on the plate carree grid. Both are compared
as bin-integrated expected counts (flat spectrum within the bin), normalized
per bin over the interior of the sub-region (the outer 2 deg excluded, where
the source map lacks emission from beyond its edge), as cell-by-cell ratios.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

DATE = "2026-10-01"
RUN = ROOT / "runs" / "fff2de5bc7e0"
OUT = ROOT / "runs" / "srcmaps"
LT = ROOT / "runs" / "_ltcubes" / "ltcube_788ab00a4d3b.fits"
IRF = "P8R3_ULTRACLEAN_V3"
CENTERS = {"b16": (0.0, 16.0), "b40": (0.0, 40.0), "corner": (-40.0, -40.0)}
HALF = 10.0            # sub-region half-width, deg
MARGIN = 2.0           # excluded from the comparison at the sub-region edge
PIX = 0.125
N_SRC = 3
FERMI_ENV = "fermi-2.5.3"


def crop_ccube(l0, b0, path):
    """The production counts cube cropped to the sub-region, header shifted."""
    from astropy.io import fits
    from totani.templates import grid_lb
    from totani.stack import _strip_to_2d
    with fits.open(RUN / "ccube.fits") as h:
        hdr = h[0].header.copy(); data = h[0].data; eb = h["EBOUNDS"].copy(); gti = h["GTI"].copy()
        L, B = grid_lb(_strip_to_2d(h[0].header))
    n = int(round(2 * HALF / PIX))
    # longitude runs to the LEFT (CDELT1 < 0): find the pixel range from the coordinate arrays
    cols = np.where(np.abs(L[0] - l0) <= HALF - PIX / 2 + 1e-9)[0]; rows = np.where(np.abs(B[:, 0] - b0) <= HALF - PIX / 2 + 1e-9)[0]
    assert len(cols) == n and len(rows) == n, (len(cols), len(rows), n)
    c0, r0 = cols.min(), rows.min()
    sub = data[:, r0:r0 + n, c0:c0 + n]
    hdr["NAXIS1"] = hdr["NAXIS2"] = n
    hdr["CRPIX1"] = float(hdr["CRPIX1"]) - c0; hdr["CRPIX2"] = float(hdr["CRPIX2"]) - r0
    fits.HDUList([fits.PrimaryHDU(sub, hdr), eb, gti]).writeto(path, overwrite=True)
    return (int(r0), int(c0), int(n))


def write_map(arr, target, path):
    from astropy.io import fits
    hdr = fits.Header()
    for k in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CRPIX1", "CRPIX2"):
        hdr[k] = target[k]
    hdr["CUNIT1"] = hdr["CUNIT2"] = "deg"; hdr["BUNIT"] = "sr-1"
    fits.PrimaryHDU(np.asarray(arr, "f4"), hdr).writeto(path, overwrite=True)


def pick_sources(l0, b0):
    from astropy.io import fits
    from totani.pointsources import default_catalog_path, extended_source_mask
    with fits.open(default_catalog_path()) as h:
        cat = h["LAT_Point_Source_Catalog"].data; ext = h["ExtendedSources"].data
    cl = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"]); cb = cat["GLAT"]
    inside = (np.abs(cl - l0) <= HALF - MARGIN - 1.0) & (np.abs(cb - b0) <= HALF - MARGIN - 1.0) & ~extended_source_mask(cl, cb, ext)
    idx = np.where(inside)[0]; idx = idx[np.argsort(-cat["Flux1000"][idx])][:N_SRC]
    return [{"name": str(cat["Source_Name"][i]).strip(), "ra": float(cat["RAJ2000"][i]), "dec": float(cat["DEJ2000"][i]),
             "l": float(cl[i]), "b": float(cb[i])} for i in idx]


def xml_model(sources, maps):
    L = ['<?xml version="1.0" ?>', '<source_library title="srcmaps benchmark">']
    for s in sources:
        L += [f'  <source name="{s["name"].replace(" ", "_")}" type="PointSource">',
              '    <spectrum type="PowerLaw"><parameter free="1" max="1e3" min="1e-5" name="Prefactor" scale="1e-9" value="1"/>'
              '<parameter free="0" max="5" min="0" name="Index" scale="-1" value="2"/><parameter free="0" max="2e5" min="30" name="Scale" scale="1" value="1000"/></spectrum>',
              f'    <spatialModel type="SkyDirFunction"><parameter free="0" max="360" min="-360" name="RA" scale="1" value="{s["ra"]:.5f}"/>'
              f'<parameter free="0" max="90" min="-90" name="DEC" scale="1" value="{s["dec"]:.5f}"/></spatialModel>', '  </source>']
    for name, path in maps.items():
        L += [f'  <source name="{name}" type="DiffuseSource">',
              '    <spectrum type="PowerLaw"><parameter free="1" max="1e3" min="1e-5" name="Prefactor" scale="1e-9" value="1"/>'
              '<parameter free="0" max="5" min="0" name="Index" scale="-1" value="2"/><parameter free="0" max="2e5" min="30" name="Scale" scale="1" value="1000"/></spectrum>',
              f'    <spatialModel file="{path}" type="SpatialMap"><parameter free="0" max="1000" min="0.001" name="Prefactor" scale="1" value="1"/></spatialModel>', '  </source>']
    L.append("</source_library>")
    return "\n".join(L) + "\n"


def prepare():
    from totani.halo import nfw_rho2_on_target
    from totani.templates import bubble_macias_on_target, target_header
    OUT.mkdir(exist_ok=True)
    target = target_header(992)
    maps = {"bubble_structured": OUT / "map_bubble_structured.fits", "nfw": OUT / "map_nfw.fits"}
    write_map(bubble_macias_on_target(target), target, maps["bubble_structured"])
    write_map(nfw_rho2_on_target(target), target, maps["nfw"])
    meta = {"centers": CENTERS, "half": HALF, "margin": MARGIN, "regions": {}}
    sh = ["#!/bin/bash", "set -e", f"source /opt/homebrew/Caskroom/miniforge/base/bin/activate {FERMI_ENV}", f"cd {OUT}"]
    for key, (l0, b0) in CENTERS.items():
        crop = crop_ccube(l0, b0, OUT / f"ccube_{key}.fits")
        srcs = pick_sources(l0, b0)
        (OUT / f"model_{key}.xml").write_text(xml_model(srcs, {k: str(v) for k, v in maps.items()}))
        meta["regions"][key] = {"l0": l0, "b0": b0, "crop": crop, "sources": srcs}
        nb = int(round(2 * (HALF + 5.0) / PIX))
        sh += [f"echo '== {key}: gtexpcube2'",
               f"gtexpcube2 infile={LT} cmap=none outfile=bexpmap_{key}.fits irfs={IRF} evtype=3 nxpix={nb} nypix={nb} binsz={PIX} "
               f"coordsys=GAL xref={l0} yref={b0} axisrot=0 proj=CAR ebinalg=FILE ebinfile={RUN}/ebounds.fits chatter=2",
               f"echo '== {key}: gtsrcmaps'",
               f"gtsrcmaps scfile=none expcube={LT} cmap=ccube_{key}.fits srcmdl=model_{key}.xml bexpmap=bexpmap_{key}.fits outfile=srcmaps_{key}.fits "
               f"irfs={IRF} evtype=3 ptsrc=yes psfcorr=yes emapbnds=no copyall=no chatter=2"]
    sh.append("echo done")
    (OUT / "run_fermitools.sh").write_text("\n".join(sh) + "\n")
    (OUT / "meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1)); print(f"wrote {OUT / 'run_fermitools.sh'}")


# ---- comparison ----------------------------------------------------------------

def bin_integrate(planes, e_mev, edges_gev, n_sub=32):
    """Flat-spectrum integral of per-energy maps across each bin, log-linear in E (as exposure.flat_spectrum_bin_exposure)."""
    out = []
    le = np.log(e_mev)
    for lo, hi in zip(edges_gev[:-1], edges_gev[1:]):
        E = np.exp(np.linspace(np.log(lo * 1e3), np.log(hi * 1e3), n_sub + 1))
        acc = np.zeros(planes.shape[1:])
        vals = []
        for e in E:
            j = np.clip(np.searchsorted(le, np.log(e)) - 1, 0, len(le) - 2)
            w = (np.log(e) - le[j]) / (le[j + 1] - le[j])
            vals.append(planes[j] * (1 - w) + planes[j + 1] * w)
        vals = np.array(vals)
        acc = np.trapezoid(vals, E, axis=0)
        out.append(acc)
    return np.array(out)


def ours_pixels(meta, edges):
    """This pipeline's PSF-folded pixel cubes on the production grid for the same sources and maps."""
    from astropy.io import fits
    from totani import fine, psf
    from totani.halo import nfw_rho2_on_target
    from totani.pointsources import default_catalog_path, exposure_pixel, source_expected_counts, dnde
    from totani.templates import bubble_macias_on_target
    g = fine._psf_grid(edges); target = g["target"]
    cubes = {"bubble_structured": bubble_macias_on_target(target)[None] * g["basis"], "nfw": nfw_rho2_on_target(target)[None] * g["basis"]}
    with fits.open(default_catalog_path()) as h:
        cat = h["LAT_Point_Source_Catalog"].data
    names = [str(n).strip() for n in cat["Source_Name"]]
    out = {n: psf.convolve_cube(c, g["B"], edges) for n, c in cubes.items()}
    exact = os.environ.get("SRC_EXACT", "1") == "1"
    for key, reg in meta["regions"].items():
        for s in reg["sources"]:
            i = names.index(s["name"])
            # a flat spectrum in the bin, as the source map comparison is spatial: unit intensity x exposure at the source pixel
            row, col = exposure_pixel(np.array([s["l"]]), np.array([s["b"]]), target)
            cnt = np.array([[g["basis"][k, row[0], col[0]] for k in range(len(edges) - 1)]])
            if exact:
                out[s["name"]] = psf.point_sources_exact(np.array([s["l"]]), np.array([s["b"]]), cnt, g["L"], g["B"], edges)
            else:
                c = np.zeros((len(edges) - 1,) + g["L"].shape)
                for k in range(len(edges) - 1):
                    c[k, row[0], col[0]] = cnt[0, k]
                out[s["name"]] = psf.convolve_cube(c, g["B"], edges)
    return out, g


def compare():
    from astropy.io import fits
    from reduce_totani import TOTANI_EDGES_GEV
    from totani.cells import bin_to_cells
    edges = TOTANI_EDGES_GEV[:10]
    meta = json.loads((OUT / "meta.json").read_text())
    ours, g = ours_pixels(meta, edges)
    E = np.sqrt(np.array(edges[:-1]) * np.array(edges[1:]))
    L = [f"# This pipeline's PSF folding against Fermitools gtsrcmaps, {DATE}" + ("" if os.environ.get("SRC_EXACT", "1") == "1" else " (point sources deposited in their pixel, the pre-10-01 placement)"), "",
         f"Sub-regions of {2*HALF:g}°; the outer {MARGIN:g}° excluded; both sides integrated over each bin for a flat spectrum and normalized per bin over the interior. "
         "Ratio = ours / gtsrcmaps per cell: median, and the 5th–95th percentile range over interior cells weighted by the gtsrcmaps counts (cells below 1% of the brightest are dropped for the point sources). "
         "'moved' = the fraction of the template's counts that this pipeline places in a different cell than gtsrcmaps (half the L1 distance of the normalized cell maps).", ""]
    res = {}
    for key, reg in meta["regions"].items():
        p = OUT / f"srcmaps_{key}.fits"
        if not p.exists():
            L += [f"## {key}: not available", ""]; continue
        r0, c0, n = reg["crop"]
        with fits.open(p) as h:
            names_h = [x.name for x in h]
            if "ENERGIES" in names_h:
                e_mev = h["ENERGIES"].data["Energy"].astype(float)
            else:
                eb = h["EBOUNDS"].data; e_mev = np.append(eb["E_MIN"].astype(float), float(eb["E_MAX"][-1])) / 1e3   # EBOUNDS are in keV
            sm = {x.name: x.data.astype(float) for x in h if x.name not in ("PRIMARY", "EBOUNDS", "GTI", "ENERGIES")}
        Lsub, Bsub = g["L"][r0:r0 + n, c0:c0 + n], g["B"][r0:r0 + n, c0:c0 + n]
        interior = (np.abs(Lsub - reg["l0"]) <= HALF - MARGIN) & (np.abs(Bsub - reg["b0"]) <= HALF - MARGIN)
        L += [f"## {key}: center ({reg['l0']:g}, {reg['b0']:g})", "", "| template | cells | " + " | ".join(f"{e:.3g}" for e in E) + " |", "|---|---|" + "---|" * 9]
        for name, planes in sm.items():
            ours_name = next((k for k in ours if k.replace(" ", "_") == name), None)
            if ours_name is None:
                continue
            fermi = bin_integrate(planes[:len(edges)], e_mev[:len(edges)], edges)
            mine = ours[ours_name][:, r0:r0 + n, c0:c0 + n]
            res[(key, name)] = {}
            for cell in (2.0, 1.0):
                m = bin_to_cells(np.where(interior, 1.0, 0.0), Lsub, Bsub, cell_deg=cell) > 0.5 * (cell / PIX) ** 2
                a = np.stack([bin_to_cells(np.where(interior, fermi[k], 0), Lsub, Bsub, cell_deg=cell) for k in range(9)])
                b = np.stack([bin_to_cells(np.where(interior, mine[k], 0), Lsub, Bsub, cell_deg=cell) for k in range(9)])
                cells = []
                for k in range(9):
                    A, Bv = a[k][m], b[k][m]
                    A, Bv = A / A.sum(), Bv / Bv.sum()
                    ok = A > 0.01 * A.max() if name not in ("bubble_structured", "nfw") else A > 0
                    r = Bv[ok] / A[ok]; w = A[ok]
                    order = np.argsort(r); cw = np.cumsum(w[order]) / w.sum()
                    q = lambda x: float(r[order][np.searchsorted(cw, x)])
                    cells.append((q(0.5), q(0.05), q(0.95), int(ok.sum()), float(0.5 * np.abs(Bv - A).sum())))
                res[(key, name)][cell] = cells
                L.append(f"| {name} | {cell:g}° ({cells[5][3]}) | " + " | ".join(f"{c[0]:.3f} [{c[1]:.2f}, {c[2]:.2f}]" for c in cells) + " |")
                L.append(f"| | moved | " + " | ".join(f"{100*c[4]:.1f}%" for c in cells) + " |")
        L.append("")
    out = ROOT / "results" / f"srcmaps_benchmark{'' if os.environ.get('SRC_EXACT', '1') == '1' else '_pixelsrc'}_{DATE}.md"
    out.write_text("\n".join(L) + "\n"); print("\n".join(L))
    (ROOT / "results" / "fits" / f"srcmaps_benchmark_{DATE}.json").write_text(json.dumps({f"{k[0]} {k[1]}": {str(c): v for c, v in d.items()} for k, d in res.items()}, indent=1))


if __name__ == "__main__":
    {"prepare": prepare, "compare": compare}[sys.argv[1]]()
