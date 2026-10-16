#!/usr/bin/env python
"""Oscar's item 3: the decisive comparisons at finer spatial resolution.

    python scripts/finer_cells_test.py                 # ~1 h, one process, ~10 GB

Every fit in the paper is on 10 x 10 deg cells. Here representative rows are
rebuilt from the pixel maps on 2 and 1 deg cells (and on 10 deg through the
same path, for a like-for-like reference) and the conditional tests are
repeated: H0/H1 TS, the 20.8 GeV halo with its profile interval, the halo
spectrum, the quadrant gate, and the flat-vs-structured preference with the
halo allowed. Totani's residual construction is still built on his 144
cells of 10 deg (that is his procedure) and only projected onto the finer
cells.

Caveat stated in the output: the expected-count maps are exposure-weighted
templates without PSF convolution; the LAT 68% containment is < 0.5 deg
above 4 GeV and ~1 deg at 1.5 GeV, so at 1 deg cells the lowest bin is the
one to read with care. Likelihood values are compared only within one
binning.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking, structures  # noqa: E402
from totani.bubble_totani import bubble_image, fit_all_cells, pixel_model_for_bin  # noqa: E402
from totani.cells import bin_to_cells, cell_centers, halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402
from totani.stack import build_background_stack, cell_template_from_shape, halo_cell_template  # noqa: E402
from totani.templates import bubble_macias_on_target  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-24"
K20 = 5
RUN = ranking.RUN
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("SO_Z6_R30_T150_C5", "gas1_ics1"),
        ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-"), ("G21_SA0_R12+iso", "gas1_ics1")]
CELLS = (10.0, 2.0, 1.0)
PANELS = {"residual": ("bubble_totani_pos", "bubble_totani_neg", "loopI_shell1", "loopI_shell2", "erosita_unif_N", "erosita_unif_S"),
          "flat": ("bubble_flat", "loopI_shell1", "loopI_shell2", "erosita_unif_N", "erosita_unif_S"),
          "structured": ("bubble_structured", "loopI_shell1", "loopI_shell2", "erosita_unif_N", "erosita_unif_S")}
NORES = {"gas1_ics1": ("bubble_structured", "erosita_unif_N", "erosita_unif_S", "loopI_haslam_raw", "sun_moon_ic"),
         "gas5_ics3": ("bubble_flat", "erosita_unif_outer_N", "erosita_unif_outer_S", "loopI_haslam_raw", "sun_moon_ic")}
_PIX: dict = {}


def pixel_library(cell_deg):
    """Everything model-independent at this cell size, built once from the pixel maps."""
    if cell_deg in _PIX:
        return _PIX[cell_deg]
    bg = build_background_stack(RUN, EDGES, bubble="ackermann", cell_deg=cell_deg)
    target, L, B, single_sr = ranking._grid_and_basis(EDGES); mask = ranking._ext_mask(L, B)
    lib = {"counts": bg["counts"], "n_bins": 9, "templates": {}, "a2020": {}}
    for n in ("bubble", "loopI_shell1", "loopI_shell2", "isotropic", "point_sources"):
        lib["templates"]["bubble_flat" if n == "bubble" else n] = bg["templates"][n]
    for n in ("HI_ring1", "HI_ring2", "HI_ring3", "HI_ring4", "H2_ring1", "H2_ring2", "H2_ring3", "H2_ring4",
              "dust_pos", "dust_neg", "ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6"):
        lib["a2020"][n] = bg["templates"][n]
    lib["templates"]["bubble_structured"] = cell_template_from_shape(RUN, EDGES, bubble_macias_on_target(target), cell_deg=cell_deg)
    for n in ("erosita_unif_N", "erosita_unif_S", "erosita_unif_outer_N", "erosita_unif_outer_S", "loopI_haslam_raw", "sun_moon_ic"):
        arr = structures.STRUCTURES[n](target, EDGES)
        if arr.ndim == 2:
            lib["templates"][n] = cell_template_from_shape(RUN, EDGES, arr, cell_deg=cell_deg)
        else:
            pix = arr * single_sr
            lib["templates"][n] = np.stack([bin_to_cells(pix[k], L, B, cell_deg=cell_deg, mask=mask) for k in range(9)])
    lib["J"] = halo_cell_template(RUN, EDGES, cell_deg=cell_deg)
    lib["keep"] = halo_cell_mask(cell_deg=cell_deg)
    Lc, Bc = cell_centers(cell_deg=cell_deg)
    lib["quads"] = {q: lib["keep"] & m for q, m in {"l+b+": (Lc > 0) & (Bc > 0), "l-b+": (Lc < 0) & (Bc > 0),
                                                    "l+b-": (Lc > 0) & (Bc < 0), "l-b-": (Lc < 0) & (Bc < 0)}.items()}
    lib["grid"] = (target, L, B, single_sr, mask)
    _PIX[cell_deg] = lib
    return lib


def row_templates(row, cell_deg, lib):
    """The diffuse templates and the residual pair for one row at this cell size."""
    target, L, B, single_sr, mask = lib["grid"]
    code, variant = row
    if code == "a2020":
        diffuse = dict(lib["a2020"])
        pm = pixel_model_for_bin(RUN, EDGES, 2, diffuse="a2020", with_gce=True)
    else:
        cubes = ranking.galprop_cubes(code, variant, target, EDGES)
        diffuse = {f"galprop_{n}": np.stack([bin_to_cells((c * single_sr)[k], L, B, cell_deg=cell_deg, mask=mask) for k in range(9)])
                   for n, c in cubes.items()}
        pm = pixel_model_for_bin(RUN, EDGES, 2, diffuse="galprop_totani", with_gce=True, galprop_cubes=cubes)
    img = bubble_image(pm, fit_all_cells(pm))                      # his 144-cell construction, unchanged
    res = {"bubble_totani_pos": cell_template_from_shape(RUN, EDGES, img["positive"], cell_deg=cell_deg),
           "bubble_totani_neg": cell_template_from_shape(RUN, EDGES, img["negative"], cell_deg=cell_deg)}
    return diffuse, res


def evaluate(counts, templates, J, keep, quads, sign_free=None):
    st = {"counts": counts, "n_bins": 9, "templates": templates}
    if sign_free is not None:                                   # the negative residual template: better sign per bin
        best = None
        for s in (1, -1):
            t2 = dict(templates); t2[sign_free] = s * templates[sign_free]
            v = sum(b["neg2lnl"] for b in fit_all_bins(restrict_stack({**st, "templates": t2}, keep), halo_mask=keep)["per_bin"])
            if best is None or v < best[0]:
                best = (v, t2)
        st["templates"] = best[1]
    sr, Jr = restrict_stack(st, keep, J)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=True)["per_bin"]
    roi = [float(Jr[k][keep].sum()) for k in range(9)]
    qa = {}
    for q, m in quads.items():
        stq, Jq = restrict_stack(st, m, J)
        one = {"counts": stq["counts"][K20:K20 + 1], "n_bins": 1, "templates": {n: t[K20:K20 + 1] for n, t in stq["templates"].items()}}
        qa[q] = fit_all_bins_profile(one, Jq[K20:K20 + 1], halo_mask=m, with_uncertainty=False)["per_bin"][0]["signal_amplitude"] * float(Jq[K20][m].sum())
    ts20 = bg[K20]["neg2lnl"] - pf[K20]["neg2lnl"]
    return {"neg2lnl_H0": float(sum(b["neg2lnl"] for b in bg)), "neg2lnl_H1": float(sum(p["neg2lnl"] for p in pf)),
            "ts20": ts20, "ts_sum": float(sum(bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(9))),
            "halo20": pf[K20]["signal_amplitude"] * roi[K20], "lo": pf[K20]["signal_amplitude_lo"] * roi[K20], "hi": pf[K20]["signal_amplitude_hi"] * roi[K20],
            "halo_bins": [pf[k]["signal_amplitude"] * roi[k] for k in range(9)], "quadrants": qa,
            "gate": bool(ts20 > 4 and pf[K20]["signal_amplitude"] > 0 and all(v > 0 for v in qa.values()))}


def main() -> int:
    t0 = time.time()
    out = {"date": DATE, "cells": list(CELLS), "results": []}
    lines = [f"# The decisive comparisons at 10, 2 and 1 degree cells ({DATE})", "",
             "Exposure-weighted templates without PSF convolution (LAT 68% containment < 0.5° above 4 GeV, ~1° at 1.5 GeV). "
             "Likelihoods compared only within one binning. 'Δ flat−struct.' is −2lnL(H1, flat) − −2lnL(H1, structured) with the halo allowed and the same companions.", "",
             "| row | cells | set | TS@20.8 / TS sum | halo@20.8 [68%] | gate | Δ flat−struct. (halo in) | halo spectrum (counts) |",
             "|---|---|---|---|---|---|---|---|"]
    for cell_deg in CELLS:
        lib = pixel_library(cell_deg)
        print(f"library at {cell_deg} deg ({time.time()-t0:.0f}s)", flush=True)
        for row in ROWS:
            try:
                diffuse, res = row_templates(row, cell_deg, lib)
            except (FileNotFoundError, KeyError) as e:
                print(f"  skip {row}: {e}"); continue
            base = {**diffuse, "isotropic": lib["templates"]["isotropic"], "point_sources": lib["templates"]["point_sources"]}
            sets = dict(PANELS)
            v = row[1] if row[1] != "-" else "gas1_ics1"
            sets["no-residual selection"] = NORES[v]
            r_flat = r_struct = None
            for sname, names in sets.items():
                t = dict(base)
                for n in names:
                    t[n] = res[n] if n in res else lib["templates"][n]
                r = evaluate(lib["counts"], t, lib["J"], lib["keep"], lib["quads"], sign_free="bubble_totani_neg" if "bubble_totani_neg" in names else None)
                if sname == "flat": r_flat = r
                if sname == "structured": r_struct = r
                rec = {"model": row[0], "variant": row[1], "cell_deg": cell_deg, "set": sname, **r}
                out["results"].append(rec)
                dfs = (r_flat["neg2lnl_H1"] - r_struct["neg2lnl_H1"]) if (r_flat and r_struct and sname == "structured") else None
                lines.append(f"| {row[0]} {row[1]} | {cell_deg:g}° | {sname} | {r['ts20']:.1f} / {r['ts_sum']:.0f} | {r['halo20']:,.0f} [{r['lo']:,.0f}, {r['hi']:,.0f}] | "
                             f"{'Y' if r['gate'] else 'n'} | {'' if dfs is None else f'{dfs:+.0f}'} | " + " ".join(f"{x:,.0f}" for x in r["halo_bins"]) + " |")
                print(lines[-1] + f"   ({time.time()-t0:.0f}s)", flush=True)
            (ROOT / "results" / "fits" / f"finer_cells_test_{DATE}.json").write_text(json.dumps(out, indent=1))
            (ROOT / "results" / f"finer_cells_test_{DATE}.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
