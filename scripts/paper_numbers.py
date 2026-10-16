#!/usr/bin/env python
"""Numbers quoted in the paper that no analysis script wrote out (audit
2026-10-15, items C): each recomputed here from the cached inputs and recorded,
with its location in the paper and its definition.

    python scripts/paper_numbers.py            # FINE_TAG not needed; ~10 min, < 4 GB

Writes results/paper_numbers_2026-10-15.{md,json}. Items:

  catalog_share       Sec. II B   catalog-source share of the observed counts, |b| >= 10, 2-deg cells
  erosita_east_share  Sec. II C   share of the brightness-weighted eROSITA template from the reconstructed East
  chz_cells           Sec. IV C   cells of the 2014 outside region touched by the enlarged CHZ boundaries
  psf_containment     App. A      R68 / R95 of the three-Gaussian kernels against the gtpsf table
  ps_redistribution   App. A      fraction of each template's counts moved between 1- and 2-deg cells by the PSF
  pion_spectrum       App. B      E^2 dN/dE peak energy and 1-10 GeV slope of the 2012 pi0 map toward the GC
  galp21_gas_ratio    App. B      Galp21 / 2012 (SL_Z6_R30_T150_C2) gas at 1.5 GeV
  galp21_swap         App. B      background-only -2lnL deficits with gas and IC swapped between Galp21 and 2012
  unit_injection      App. D      counts of a unit constant-coefficient NFW injection, 1.5 and 100 GeV
  spectral_slopes     App. D      E^2 dN/dE log-log slopes of the null-sky and data spectra
  erosita_ew          App. E      East/West brightness ratios, nominal and the second published display (ED1)
  hestia_shapes       App. E      HESTIA DM axis ratio c/a at 1.4-30 kpc and minor-axis alignment with the disk
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402

DATE = "2026-10-15"
EDGES = TOTANI_EDGES_GEV[:10]
E_CENTER = np.sqrt(np.asarray(EDGES[:-1]) * np.asarray(EDGES[1:]))
CACHE = ROOT / "results" / "fits" / "ranking_cells"
HESTIA = Path(os.environ.get("HESTIA_ROOT", Path.home() / "Data2" / "hestia_muru"))
OUT: dict = {}


def record(key, where, definition, values, summary):
    OUT[key] = {"paper": where, "definition": definition, "values": values, "summary": summary}
    print(f"{key}: {summary}", flush=True)


def catalog_share():
    z = np.load(CACHE / "cells2_psfuc" / "library.npz")
    keep = z["keep"]; T = z["T__point_sources"]; n = z["counts"]
    f = [float(T[k][keep].sum() / n[k][keep].sum()) for k in range(9)]
    record("catalog_share", "Sec. II B (catalog sources predict ... of the observed counts)",
           "point-source template at its catalog normalization over the observed counts, fitted 2-deg cells "
           "(|b| >= 10, extended sources masked), ULTRACLEAN response, per energy bin",
           {"E_gev": E_CENTER.tolist(), "fraction": f},
           f"{100*min(f):.1f}%-{100*max(f):.1f}% over 1.16-130 GeV; {100*f[5]:.1f}% at 20.8 GeV")


def erosita_east_share():
    from totani.erosita import brightness_on_target
    from totani.templates import grid_lb, pixel_solid_angle, target_header
    tgt = target_header(); L, B = grid_lb(tgt); om = pixel_solid_angle(tgt)
    roi = (np.abs(L) <= 60) & (np.abs(B) >= 10) & (np.abs(B) <= 60)
    out = {}
    for lab, key in (("north", "whole_N"), ("south", "whole_S")):
        w = brightness_on_target(tgt, key) * om
        out[lab] = float(w[roi & (L > 0)].sum() / w[roi].sum())
    both = None
    record("erosita_east_share", "Sec. II C (this reconstruction supplies approximately ...% of the template weight)",
           "fraction of each brightness-weighted lobe's weight (template x pixel solid angle) in the fitted region "
           "|l| <= 60, 10 <= |b| <= 60 that lies at l > 0 (the reconstructed East)",
           out, f"north {100*out['north']:.1f}%, south {100*out['south']:.1f}%")


def chz_cells():
    z = np.load(CACHE / "cells2" / "library.npz")
    keep = z["keep"]; out14 = keep & (np.abs(z["T__bubble_flat"]).sum(0) == 0)
    fp = np.load(CACHE / "cells2_psfuc" / "cholis2026.npz")["chz_flatE__footprint"]
    J = np.load(CACHE / "cells2_psfuc" / "library.npz")["J"][5]
    touched = out14 & fp
    record("chz_cells", "Sec. IV C (the enlarged boundaries ... intersect N of the 2190 cells, containing x% of the NFW counts)",
           "2-deg cells of the 2014 outside region touched by any pixel of the CHZ edges-E footprint; their share of "
           "the response-folded NFW counts at 20.8 GeV in that region",
           {"outside_cells": int(out14.sum()), "touched_cells": int(touched.sum()),
            "nfw_share": float(J[touched].sum() / J[out14].sum())},
           f"{int(touched.sum())} of {int(out14.sum())} cells, {100*J[touched].sum()/J[out14].sum():.1f}% of the NFW counts")


def psf_containment():
    from totani import psf, ranking
    _, L, B, _ = ranking._grid_and_basis(EDGES)
    os.environ.setdefault("FINE_TAG", "psfuc")
    rows = psf.check_point_source(L, B, EDGES, l0=10.0, b0=30.0)
    d68 = [r["R68"] / r["R68_table"] - 1 for r in rows]; d95 = [r["R95"] / r["R95_table"] - 1 for r in rows]
    vals = {"E_gev": E_CENTER.tolist(), "R68_deg": [r["R68"] for r in rows], "R68_table_deg": [r["R68_table"] for r in rows],
            "R95_deg": [r["R95"] for r in rows], "R95_table_deg": [r["R95_table"] for r in rows],
            "R68_dev": d68, "R95_dev": d95, "flux_ratio": [r["flux_ratio"] for r in rows]}
    below_pixel = [float(E_CENTER[k]) for k in range(9) if rows[k]["R68_table"] < 0.125]
    record("psf_containment", "App. A (kernels reproduce the 68% and 95% containment radii ...)",
           "a delta at (l, b) = (10, 30) deg convolved with each bin's three-Gaussian kernel; measured R68, R95 against "
           "the gtpsf table for P8R3_ULTRACLEAN_V3 evtype 3",
           vals, "R68 deviations " + ", ".join(f"{100*x:+.1f}%" for x in d68[:5]) + " at 1.5-12.3 GeV; R95 deviations "
           + ", ".join(f"{100*x:+.1f}%" for x in d95[:5]) + f"; table R68 < 0.125 deg from {min(below_pixel):.1f} GeV")


def ps_redistribution():
    out = {}
    for cell in ("cells2", "cells1"):
        plain = np.load(CACHE / cell / "library.npz")
        keep = plain["keep"]
        for lib_name in ("library_pixelsrc", "library"):
            path = CACHE / f"{cell}_psfuc" / f"{lib_name}.npz"
            if not path.exists():
                continue
            z = np.load(path)
            moved = {}
            for t in ("point_sources", "bubble_flat", "bubble_structured", "isotropic"):
                if f"T__{t}" not in z or f"T__{t}" not in plain:
                    continue
                a, b = plain[f"T__{t}"], z[f"T__{t}"]
                moved[t] = [float(0.5 * np.abs(a[k][keep] / a[k][keep].sum() - b[k][keep] / b[k][keep].sum()).sum()) for k in range(9)]
            out[f"{cell}/{lib_name}"] = moved
        rp, rq = (np.load(CACHE / c / "row__SL_Z6_R30_T150_C2__gas1_ics1.npz") for c in (cell, f"{cell}_psfuc"))
        out[f"{cell}/diffuse_SLZ6R30T150C2"] = {t: [float(0.5 * np.abs(rp[f"D__{t}"][k][keep] / rp[f"D__{t}"][k][keep].sum()
                                                       - rq[f"D__{t}"][k][keep] / rq[f"D__{t}"][k][keep].sum()).sum()) for k in range(9)]
                                                for t in ("galprop_gas", "galprop_ics")}
    def at100(key):
        d = out.get(key, {}); return ", ".join(f"{t} {100*v[8]:.1f}%" for t, v in d.items())
    record("ps_redistribution", "App. A (point sources show the largest redistribution between cells ...)",
           "half the L1 distance between the normalized fitted-cell maps of a template without and with the PSF "
           "folding (1- and 2-deg cells): the fraction of its counts that the response moves to other cells. "
           "'library_pixelsrc' = sources at pixel centres (before 2026-10-01), 'library' = exact positions",
           out, "at 100 GeV -- " + "; ".join(f"{k}: {at100(k)}" for k in out))


def pion_spectrum():
    from totani.galprop import load_skymap
    from totani.healpix import ang2pix_ring
    spec, E, nside = load_skymap("pi0_decay")
    gc = ang2pix_ring(nside, 0.0, 0.0)
    e2 = spec[gc] * E ** 2
    i = int(np.argmax(e2))
    x = np.log(E[i - 1:i + 2]); c = np.polyfit(x, np.log(e2[i - 1:i + 2]), 2)
    e_peak = float(np.exp(-c[1] / (2 * c[0])))           # parabola in (ln E, ln E^2F) through the maximum
    w = (E >= 1000.0) & (E <= 10000.0)
    slope = float(np.polyfit(np.log(E[w]), np.log(spec[gc, w]), 1)[0])
    record("pion_spectrum", "App. B (the pion-decay E^2 dN/dE spectrum peaks near ... and has a slope of ...)",
           "2012 grid SL_Z6_R30_T150_C2 pi0 HEALPix map, pixel at the Galactic Centre: energy of the E^2 dN/dE maximum, "
           "from a parabola in ln E through the tabulated maximum and its neighbours, and the least-squares dN/dE "
           "slope over the tabulated energies 1-10 GeV",
           {"E_peak_mev": e_peak, "E_tabulated_max_mev": float(E[i]), "slope_dnde_1_10gev": slope, "n_energies": int(w.sum())},
           f"peak near {e_peak:.0f} MeV (tabulated maximum {E[i]:.0f} MeV); dN/dE slope {slope:.2f} over 1-10 GeV")


def galp21_gas_ratio():
    from totani import galprop
    from totani.templates import target_header
    from totani.cells import halo_cell_mask
    tgt = target_header()
    ref = galprop.galprop_on_target(tgt, EDGES)["gas"][0]
    from totani.templates import grid_lb
    L, B = grid_lb(tgt)
    roi = (np.abs(L) <= 60) & (np.abs(B) >= 10) & (np.abs(B) <= 60)
    out = {}
    for m in ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12"):
        g = np.load(CACHE / f"galp21_cubes__{m}.npz")["gas_total"][0]
        out[m] = {"total_ratio": float(g[roi].sum() / ref[roi].sum()),
                  "pixel_ratio_median": float(np.median(g[roi] / ref[roi]))}
    t = [v["total_ratio"] for v in out.values()]
    record("galp21_gas_ratio", "App. B (Galp21 gas templates agree with the 2012 baseline to factors of ...)",
           "Galp21 gas (pi0 + bremsstrahlung, all rings) over the 2012 SL_Z6_R30_T150_C2 gas, summed over "
           "|l| <= 60, 10 <= |b| <= 60, in the 1.16-1.96 GeV bin",
           out, f"total ratios {min(t):.3f}-{max(t):.3f} over the six models")


def galp21_swap():
    from totani import ranking, selection
    sh = ranking.shared_cells(EDGES)
    v = "gas1_ics1"
    ref = {c: float(selection.bg_neg2lnl(selection.base_stack(c, v, EDGES, sh)).sum())
           for c in ranking.MODEL_CODES if ranking.cache_path("model", f"{c}__{v}").exists()}
    best = min(ref, key=ref.get)
    b_st = selection.base_stack(best, v, EDGES, sh)
    out = {"best_2012": best, "best_2012_neg2lnl": ref[best], "models": {}}
    for m in ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12"):
        an = selection.base_stack(f"G21_{m}", v, EDGES, sh)
        iso = selection.base_stack(f"G21_{m}+iso", v, EDGES, sh)
        def fit(gas_from, ics_from):
            st = {**b_st, "templates": dict(b_st["templates"])}
            st["templates"]["galprop_gas"] = gas_from["templates"]["galprop_gas"]
            st["templates"]["galprop_ics"] = ics_from["templates"]["galprop_ics"]
            return float(selection.bg_neg2lnl(st).sum()) - ref[best]
        out["models"][m] = {"galp21_aniso": float(selection.bg_neg2lnl(an).sum()) - ref[best],
                            "gas2012_icsG21aniso": fit(b_st, an), "gasG21_ics2012": fit(an, b_st),
                            "gas2012_icsG21iso": fit(b_st, iso), "gasG21_icsG21iso": fit(an, iso)}
        print(m, out["models"][m], flush=True)
    def rng(k):
        x = [d[k] for d in out["models"].values()]; return f"{min(x):,.0f}-{max(x):,.0f}"
    record("galp21_swap", "App. B (2012 gas with Galp21 anisotropic inverse Compton gives a deficit of ...)",
           "background-only -2lnL (10-deg cells, |b| >= 10, nine bins) of diffuse + isotropic + point sources, single-gas "
           "decomposition, relative to the best 2012 model, with the gas and inverse Compton templates taken from the "
           "models named", out,
           f"best 2012 = {best}; Galp21 as delivered {rng('galp21_aniso')}; 2012 gas + G21 aniso IC {rng('gas2012_icsG21aniso')}; "
           f"G21 gas + 2012 IC {rng('gasG21_ics2012')}; 2012 gas + G21 iso IC {rng('gas2012_icsG21iso')}; "
           f"G21 gas + G21 iso IC {rng('gasG21_icsG21iso')}")


def unit_injection():
    d = json.loads((ROOT / "results" / "fits" / "injection_refined_2026-09-29.json").read_text())
    counts = [float(d["a_ref"] * j) for j in d["roi_J"]]
    record("unit_injection", "App. D (the unit injection gives ... counts at 1.5 GeV and ... at 100 GeV)",
           "expected counts of the unit injection (constant coefficient a_ref of injection_refined_2026-09-29.json) "
           "times the response-folded NFW J summed over the fitted region, per energy bin",
           {"E_gev": E_CENTER.tolist(), "a_ref": d["a_ref"], "counts": counts},
           f"{counts[0]:,.0f} counts at 1.5 GeV and {counts[8]:,.0f} at 100 GeV per unit coefficient")


def spectral_slopes():
    d = json.loads((ROOT / "results" / "fits" / "halo_spectra_compare_2026-09-24.json").read_text())
    E = np.asarray(d["E_gev"]); out = {}
    for name, s in d["series"].items():
        y = np.asarray(s.get("e2f_median", s.get("e2f")), float); ok = y > 0
        if ok.sum() >= 3:
            out[name] = {"slope": float(np.polyfit(np.log(E[ok]), np.log(y[ok]), 1)[0]), "bins_used": int(ok.sum())}
    record("spectral_slopes", "App. D (coarse null spectra ... slope ..., compared with ... for the flat-template data)",
           "least-squares slope of ln(E^2 dN/dE) against ln E over the bins with positive E^2 dN/dE, 10-deg cells "
           "(results/fits/halo_spectra_compare_2026-09-24.json)", out,
           "; ".join(f"{k} {v['slope']:+.2f} ({v['bins_used']} bins)" for k, v in out.items()))


def erosita_ew():
    from astropy.io import fits
    from totani.healpix import pix2ang_ring
    d = fits.open(ROOT / "external" / "eROSITA_fullsky_lat_handoff" / "systematics"
                  / "erosita_brightness_uncertainties.fits")["UNCERTAINTY"].data
    nside = 256
    l, b = pix2ang_ring(nside, np.arange(12 * nside * nside)); l = np.where(l > 180, l - 360, l)
    roi = (np.abs(l) <= 60) & (np.abs(b) >= 10) & (np.abs(b) <= 60)
    def ratio(col, valid):
        v = np.asarray(d[col], float); m = np.asarray(d[valid]); o = {}
        for lab, sel in (("north", b > 0), ("south", b < 0)):
            k = m & sel & (v > 0) & roi
            o[lab] = float(v[k & (l > 0)].mean() / v[k & (l < 0)].mean())
        return o
    out = {"nominal": ratio("FULL_NOMINAL", "VALID_FULL"), "ED1": ratio("FULL_ED1", "VALID_ED1")}
    record("erosita_ew", "App. E (East/West ratios ...; another published display gives ...)",
           "mean East (l > 0) over mean West surface brightness in the fitted region, per lobe, HEALPix 256 "
           "(as tests/test_erosita_brightness.py)", out,
           f"nominal {out['nominal']['north']:.3f} / {out['nominal']['south']:.3f}; ED1 {out['ED1']['north']:.3f} / {out['ED1']['south']:.3f}")


def hestia_shapes():
    """AHF radius is in kpc/h; negative radii mark bins inside the convergence radius
    (Power et al. 2003) and are dropped, as in gce-hestia-comment/code/ahf_frames.py."""
    h = 0.677
    R_MIN = 8.2 * np.sin(np.radians(10.0))   # smallest Galactocentric radius on any |b| >= 10 deg sightline
    frames = json.loads((HESTIA / "muru_frames.json").read_text())
    name_of = {"G12": "G1.1", "G11": "G1.2", "G22": "G2.1", "G21": "G2.2", "G32": "G3.1", "G31": "G3.2"}
    out = {}
    for prof, halo in name_of.items():
        a = np.loadtxt(HESTIA / "profiles" / f"profile_{prof}_AHF.txt")
        r_kpc = a[:, 0] / h
        sel = (a[:, 0] > 0) & (r_kpc >= R_MIN) & (r_kpc <= 30.0)
        c = a[sel, 12]; Ec = a[sel, 19:22]
        z = np.asarray(frames[halo]["z_hat"], float); z = z / np.linalg.norm(z)
        ang = [float(np.degrees(np.arccos(min(1.0, abs(e @ z) / np.linalg.norm(e))))) for e in Ec]
        out[halo] = {"profile": prof, "n_bins": int(sel.sum()), "r_kpc_range": [float(r_kpc[sel].min()), float(r_kpc[sel].max())],
                     "c_over_a_min": float(c.min()), "c_over_a_max": float(c.max()),
                     "minor_axis_vs_disk_normal_deg_max": max(ang)}
    cmin = min(v["c_over_a_min"] for v in out.values()); cmax = max(v["c_over_a_max"] for v in out.values())
    amax = max(v["minor_axis_vs_disk_normal_deg_max"] for v in out.values())
    record("hestia_shapes", "App. E (the HESTIA halos have c/a = ... within 30 kpc and minor axes within ... of the disk normal)",
           "AHF DM profiles of the six HESTIA analogues (Muru et al. 2025), converged radial bins between R0 sin 10 deg "
           "= 1.42 kpc (the closest approach to the Galactic Center of any sightline at |b| >= 10 deg) and 30 kpc: "
           "minor-to-major axis ratio c/a, and the largest angle between the minor axis Ec and the stellar disk "
           "normal z_hat (muru_frames.json)", out,
           f"c/a {cmin:.2f}-{cmax:.2f}; minor axis within {amax:.1f} deg of the disk normal")


def main() -> int:
    for f in (catalog_share, erosita_east_share, chz_cells, psf_containment, ps_redistribution, pion_spectrum,
              galp21_gas_ratio, unit_injection, spectral_slopes, erosita_ew, hestia_shapes, galp21_swap):
        try:
            f()
        except Exception as e:                       # record the failure rather than hide it
            OUT[f.__name__] = {"error": repr(e)}; print(f"{f.__name__}: ERROR {e!r}", flush=True)
    (ROOT / "results" / "fits" / f"paper_numbers_{DATE}.json").write_text(json.dumps(OUT, indent=1))
    L = [f"# Numbers quoted in the paper, recomputed and recorded ({DATE})", "",
         "Each number the audit of 2026-10-15 found without a script-written source, recomputed from the cached "
         "inputs by `scripts/paper_numbers.py`. Values and definitions in `results/fits/paper_numbers_{DATE}.json`.", "",
         "| item | paper | result | definition |", "|---|---|---|---|"]
    for k, v in OUT.items():
        if "error" in v:
            L.append(f"| {k} | | ERROR: {v['error']} | |")
        else:
            L.append(f"| {k} | {v['paper']} | {v['summary']} | {v['definition']} |")
    (ROOT / "results" / f"paper_numbers_{DATE}.md").write_text("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
