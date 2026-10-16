#!/usr/bin/env python
"""Kev 2026-10-11: test the Fermi-bubble model of Cholis, Hu & Zhong
(arXiv:2610.10695, Zenodo 23199004) inside our likelihood. Their final model
replaces the halo-like template with three bubble components -- flat bubbles
with enlarged edges ('E'), 'Inner Bubbles' (zero at |b| < 21 deg and outside
the edges), and a tilted pair of 'Cocoons' (southern 10x the northern) -- but
they never fit an NFW component on top of that model, and they run no
injection or null simulations. Here, on the PSF-folded 2-deg cells:

  data    For the six representative rows, with the standard companions
          (eROSITA uniform N/S, Haslam Loop I, solar IC), each bubble
          description -- 2014 flat, structured, theirs (flat-E alone,
          flat-E + cocoons, flat-E + inner + cocoons) -- gives the H0 and
          signed NFW profile fits per bin, A(20.8) with its 68% interval,
          Q+ / Z_2018 over 5.6-130 GeV, and -2lnL for the conditional
          comparisons. The same fits OUTSIDE the bubbles: our region (cells
          with no 2014 flat-mask coverage) and that region with their
          edges-E footprint also removed. For the three quadrant-calibration
          data rows, the standard and 45-deg-rotated Q6 with their model.
  sims    Truth = the real-data fit of the true row (SO_Z6_R30_T150_C5,
          1-gas) with their full model: (i) background only, Poisson skies
          refitted with their model + NFW (spurious NFW) and with the 2014
          flat mask + NFW (bubble-substructure mismatch); (ii) the same
          background + Totani's spectrum (6412 counts at 20.8 GeV) refitted
          with their model + NFW (recovery); (iii) background + the fitted
          common halo, refitted for the quadrant statistics (null
          calibration of Q6 under their model).

    FINE_TAG=psfuc python scripts/cholis_bubbles.py --stage data --workers 6
    FINE_TAG=psfuc python scripts/cholis_bubbles.py --stage sims --skies 50 --workers 6
    FINE_TAG=psfuc python scripts/cholis_bubbles.py --render

Their FITS files carry no WCS. Orientation (FITS axis 1 = l increasing with
column, axis 2 = b, 0.25-deg pixels with centers at -69.875 ... +69.875) was
established from the map content (overlap with the 2014 template, the North
Polar Spur in their Loop I map); pixel edges are taken at -70 + 0.25 i
(+-0.125 deg ambiguity, far below the 2-deg cells). The spatial shape of each
component is energy independent (plane ratios constant to 1e-6), so one plane
gives the shape and each component has a free normalization per bin.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import oscar_conditional as OC  # noqa: E402
from totani import ranking  # noqa: E402
from ns_halo_spectra import joint  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine, simsky, stats  # noqa: E402
from totani.likelihood import expected_counts, fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

DATE = "2026-10-11"
CELL = 2.0
EDGES = TOTANI_EDGES_GEV[:10]
K20 = 5
RESP = OC.RESP_BINS
ZEN = ROOT / "external" / "cholis2026"
PIECES = {"chz_flatE": "Bubbles_Edges_Alternative_E_*.fits", "chz_inner": "BubblesInner_with_Boundaries_*.fits",
          "chz_cocoons": "Cocoons_*.fits", "chz_loop_bright": "LoopI_Bright_Shell_*.fits", "chz_loop_dim": "LoopI_Dim_Shell_*.fits"}
CHZ_LOOP = ("chz_loop_bright", "chz_loop_dim")
CONFIGS = {"flat": None, "structured": None,                       # ours, through oscar_conditional
           "chz_flatE": ("chz_flatE",), "chz_flatE_cocoons": ("chz_flatE", "chz_cocoons"),
           "chz_full": ("chz_flatE", "chz_inner", "chz_cocoons"),
           # their two Loop I shells (free per bin) in place of our Haslam Loop I template
           "chz_full_chzloop": ("chz_flatE", "chz_inner", "chz_cocoons") + CHZ_LOOP}
ROWS = OC.ROWS
QROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("a2020", "-")]
TRUE_ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
_W: dict = {}


def shapes_on_target(g):
    """Their three components as energy-independent pixel shapes on our 0.125-deg grid."""
    from astropy.io import fits
    L, B = g["L"], g["B"]
    j = np.floor((L + 70.0) / 0.25).astype(int); i = np.floor((B + 70.0) / 0.25).astype(int)
    assert j.min() >= 0 and j.max() < 560 and i.min() >= 0 and i.max() < 560
    out = {}
    for name, pat in PIECES.items():
        (f,) = glob.glob(str(ZEN / pat))
        with fits.open(f, memmap=True) as h:
            plane = np.array(h[0].data[15], float)          # (b, l); shape is the same in every plane
        out[name] = plane[i, j]
    return out


def cache_path(cell=None):
    return fine.cache_dir(CELL if cell is None else cell) / "cholis2026.npz"


def build():
    """Fold every component in PIECES onto this cell size; components already cached are kept."""
    p = cache_path()
    out = dict(np.load(p)) if p.exists() else {}
    if all(k in out for k in PIECES):
        return p
    from totani.cells import bin_to_cells
    g = fine._psf_grid(EDGES)
    sh = shapes_on_target(g)
    for name, s in sh.items():
        if name in out:
            continue
        out[name] = fine._cells_psf(s[None] * g["basis"], g, CELL, EDGES)
        out[f"{name}__footprint"] = bin_to_cells((s > 0).astype(float), g["L"], g["B"], cell_deg=CELL, mask=g["mask"]) > 0
        print(f"built {name} at {CELL:g} deg", flush=True)
    np.savez(p, **out)
    return p


def _init(cell=2.0):
    global CELL
    CELL = cell
    OC._init(CELL); _W["lib"] = OC._W["lib"]; lib = _W["lib"]; keep = lib["keep"]
    z = np.load(cache_path()); _W["chz"] = {k: z[k] for k in PIECES}
    # Outside regions are defined once on the 2-deg cells and carried to finer cells (each 2-deg
    # cell is exactly f x f finer cells), so both binnings fit the same sky, as in the paper.
    plain = np.load(ranking.CACHE / "cells2" / "library.npz")
    out2 = np.abs(plain["T__bubble_flat"]).sum(0) == 0
    fp2 = np.load(cache_path(2.0))["chz_flatE__footprint"]
    f = int(round(2.0 / CELL)); up = lambda m: np.kron(m, np.ones((f, f), bool))
    out14 = keep & up(out2)
    _W["regions"] = {"all": keep, "outside": out14, "outside_union": out14 & up(~fp2)}
    J = lib["J"][K20]
    _W["J_frac"] = {r: float(J[m].sum() / J[keep].sum()) for r, m in _W["regions"].items()}
    q = np.load(fine.cache_dir(CELL) / "halo_quadrants.npz")
    _W["quads"] = {d: [q[k] for k in q.files if k.startswith(d + "__")] for d in ("std", "rot")}
    _W["amp"] = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")


def stack_for(counts, row_key, config):
    if CONFIGS[config] is None:
        return OC.stack_for(counts, row_key, config)
    lib = _W["lib"]; row = fine.load_row(CELL, *row_key)
    st = fine.base_stack(lib, row, counts)
    for n in OC.COMP:
        st["templates"][n] = lib["templates"][n]
    for n in CONFIGS[config]:
        st["templates"][n] = _W["chz"][n]
    if any(n in CHZ_LOOP for n in CONFIGS[config]):
        st["templates"].pop("loopI_haslam_raw")
    return st


def fit(counts, row_key, config, region="all", with_intervals=True):
    """H0 and signed NFW profile per bin in a region; the oscar_conditional summaries."""
    lib = _W["lib"]; J = lib["J"]; m = _W["regions"][region]
    st = stack_for(counts, row_key, config)
    sr, Jr = restrict_stack(st, m, J)
    bg = fit_all_bins(sr, halo_mask=m)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=m, with_uncertainty=with_intervals, n_refine=OC.N_REFINE)["per_bin"]
    roi = [float(Jr[k][m].sum()) for k in range(9)]
    a = [pf[k]["signal_amplitude"] for k in range(9)]
    r = {"v0": [bg[k]["neg2lnl"] for k in range(9)], "v1": [pf[k]["neg2lnl"] for k in range(9)],
         "a": a, "counts": [a[k] * roi[k] for k in range(9)], "roiJ": roi}
    r["ts"] = [r["v0"][k] - r["v1"][k] for k in range(9)]
    r["resp"] = stats.summarize(r["ts"], a, RESP); r["all"] = stats.summarize(r["ts"], a)
    if with_intervals:
        r["lo"] = [pf[k].get("signal_amplitude_lo", np.nan) * roi[k] for k in range(9)]
        r["hi"] = [pf[k].get("signal_amplitude_hi", np.nan) * roi[k] for k in range(9)]
    return r


def qquad(counts, row_key, config):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    st = stack_for(counts, row_key, config)
    sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    res = {"std": [], "rot": []}
    for k in range(9):
        T = np.stack([sr["templates"][n][k] for n in names]); n = sr["counts"][k]
        v1, _ = joint(T, n, [Jr[k]])
        for d in ("std", "rot"):
            v4, _ = joint(T, n, [q[k] * keep for q in _W["quads"][d]])
            res[d].append(float(v1 - v4))
    return res


def truth(config, with_halo, row=TRUE_ROW):
    """Expected counts of `row` fitted to the data with `config`: background
    only (with_halo False) or background + the fitted common halo."""
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    st = stack_for(lib["counts"], row, config)
    sr, Jr = restrict_stack(st, keep, J)
    names = sorted(sr["templates"]); full = np.stack([st["templates"][n] for n in names])
    mu = np.zeros_like(lib["counts"], dtype=float)
    if with_halo:
        pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False, n_refine=OC.N_REFINE)
        order = pf["component_order"]; pb = pf["per_bin"]
        for k in range(9):
            bga = np.array([pb[k]["bg_amplitudes"][n] for n in order])
            mu[k] = expected_counts(np.log10(np.clip(bga, 1e-300, None)), np.stack([st["templates"][n][k] for n in order])) + pb[k]["signal_amplitude"] * J[k]
    else:
        fb = fit_all_bins(sr, halo_mask=keep)
        order = fb["component_order"]; pb = fb["per_bin"]
        for k in range(9):
            bga = np.array([pb[k]["amplitudes"][n] for n in order])
            mu[k] = expected_counts(np.log10(np.clip(bga, 1e-300, None)), np.stack([st["templates"][n][k] for n in order]))
    return mu


def task(t):
    kind = t[0]; lib = _W["lib"]; keep = lib["keep"]
    if kind == "data":
        _, row_key, config, region = t
        return {"kind": kind, "row": list(row_key), "config": config, "region": region,
                **fit(lib["counts"], row_key, config, region, True), "J_frac": _W["J_frac"][region]}
    if kind == "dataq":
        _, row_key, config = t
        return {"kind": kind, "row": list(row_key), "config": config, **qquad(lib["counts"], row_key, config)}
    if kind == "shape":
        _, row_key, config = t
        return {"kind": kind, "row": list(row_key), "config": config, **shape_test(row_key, config)}
    if kind == "qnull_row":
        _, i, mu, row_key = t
        rng = np.random.default_rng([31, i, 20261011, [r for r in QROWS].index(row_key)])
        counts = lib["counts"].astype(float); counts[:, keep] = rng.poisson(np.clip(mu[:, keep], 0, None))   # fitted cells only
        return {"kind": kind, "sim": i, "row": list(row_key), "chz_full": qquad(counts, row_key, "chz_full")}
    _, i, mu, fit_configs = t
    rng = np.random.default_rng([29, i, 20261011, ["null", "inj", "qnull"].index(kind)])
    counts = lib["counts"].astype(float); counts[:, keep] = rng.poisson(np.clip(mu[:, keep], 0, None))   # fitted cells only
    out = {"kind": kind, "sim": i}
    for c in fit_configs:
        if kind == "qnull":
            out[c] = qquad(counts, TRUE_ROW, c)
        else:
            r = fit(counts, TRUE_ROW, c, "all", False)
            out[c] = {"counts": r["counts"], "ts": r["ts"], "Z_2018": r["resp"]["Z_2018"], "Q_plus": r["resp"]["Q_plus"]}
    return out


def shape_test(row_key, config):
    """The Sec. IV B spectral-shape test (dm_shape_test.py) with this configuration's stack."""
    import dm_shape_test as DS
    from scipy.stats import chi2
    from totani import dmspectra
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    st = stack_for(lib["counts"], row_key, config)
    out = {}
    for rname, bins in DS.RANGES.items():
        s, Js, edges = DS.sub_stack(st, J, bins)
        sr, Jr = restrict_stack(s, keep, Js)
        pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False, n_refine=OC.N_REFINE)["per_bin"]
        v_free = float(sum(p["neg2lnl"] for p in pf))
        rec = {"v_free": v_free, "K": len(bins), "spectra": {}}
        for ch, m in DS.SPECTRA:
            jf = dmspectra.joint_fit(s, Js, ch, m, edges, halo_mask=keep, n_scan=31, with_intervals=False)
            d = jf["neg2lnl_best"] - v_free
            rec["spectra"][f"{ch} {m:g}"] = {"sigmav_best": jf["sigmav_best"], "ts_vs_zero": jf["ts_vs_zero"], "delta": d,
                                             "p_shape": float(chi2.sf(max(d, 0.0), len(bins) - 1))}
        out[rname] = rec
    return out


def run(tasks, workers):
    res = []; t0 = time.time()
    with Pool(workers, initializer=_init, initargs=(CELL,)) as pool:
        for k, r in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(r)
            if k % 10 == 0 or k == len(tasks):
                print(f"[{k}/{len(tasks)}] ({time.time()-t0:.0f}s)", flush=True)
    return res


def main() -> int:
    global CELL
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", choices=("data", "sims", "loop", "shape", "qnull_rows"))
    ap.add_argument("--workers", type=int, default=6); ap.add_argument("--cell", type=float, default=2.0)
    ap.add_argument("--skies", type=int, default=50); ap.add_argument("--render", action="store_true"); a = ap.parse_args()
    CELL = a.cell
    assert fine.tag() == "_psfuc", "run with FINE_TAG=psfuc (the production response)"
    fits_dir = ROOT / "results" / "fits"
    if a.render:
        import cholis_bubbles_render
        cholis_bubbles_render.render(); return 0
    build()
    sfx = "" if CELL == 2.0 else f"_cells{CELL:g}"
    if a.stage == "data" and CELL != 2.0:
        tasks = [("data", r, c, "all") for r in ROWS for c in ("flat", "structured", "chz_full")]
        tasks += [("data", r, c, g) for r in ROWS for c in ("flat", "chz_full") for g in ("outside", "outside_union")]
        tasks += [("dataq", r, c) for r in QROWS for c in ("flat", "chz_full")]
        res = run(tasks, a.workers)
        (fits_dir / f"cholis_bubbles_data{sfx}_{DATE}.json").write_text(json.dumps({"date": DATE, "cell": CELL, "results": res}, indent=1, default=float))
        return 0
    if a.stage == "loop":
        tasks = [("data", r, "chz_full_chzloop", g) for r in ROWS for g in ("all", "outside_union")]
        tasks += [("dataq", r, "chz_full_chzloop") for r in QROWS]
        res = run(tasks, a.workers)
        (fits_dir / f"cholis_bubbles_loop_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1, default=float))
        return 0
    if a.stage == "shape":
        tasks = [("shape", r, c) for r in ROWS for c in ("chz_full", "chz_full_chzloop")]
        res = run(tasks, a.workers)
        (fits_dir / f"cholis_bubbles_shape_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1, default=float))
        return 0
    if a.stage == "qnull_rows":
        _init(CELL)
        tasks = []
        for r in QROWS:
            if r == TRUE_ROW:
                continue                                # already calibrated in the 'sims' stage
            mu = truth("chz_full", True, r)
            tasks += [("qnull_row", i, mu, r) for i in range(4 * a.skies)]
        res = run(tasks, a.workers)
        (fits_dir / f"cholis_bubbles_qnullrows_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1, default=float))
        return 0
    if a.stage == "data":
        tasks = [("data", r, c, "all") for r in ROWS for c in CONFIGS]
        tasks += [("data", r, c, g) for r in ROWS for c in ("flat", "structured", "chz_full") for g in ("outside", "outside_union")]
        tasks += [("dataq", r, c) for r in QROWS for c in ("flat", "chz_full")]
        res = run(tasks, a.workers)
        (fits_dir / f"cholis_bubbles_data_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1, default=float))
    else:
        _init()
        mu0 = truth("chz_full", False); mu1 = truth("chz_full", True)
        inj = mu0 + _W["amp"][:, None, None] * _W["lib"]["J"]
        tasks = [("null", i, mu0, ("chz_full", "flat")) for i in range(a.skies)]
        tasks += [("inj", i, inj, ("chz_full",)) for i in range(a.skies)]
        tasks += [("qnull", i, mu1, ("chz_full",)) for i in range(4 * a.skies)]
        res = run(tasks, a.workers)
        injected = (_W["amp"] * np.array([float(_W["lib"]["J"][k][_W["lib"]["keep"]].sum()) for k in range(9)])).tolist()
        (fits_dir / f"cholis_bubbles_sims_{DATE}.json").write_text(json.dumps(
            {"date": DATE, "skies": a.skies, "true_row": list(TRUE_ROW), "injected_counts": injected, "results": res}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
