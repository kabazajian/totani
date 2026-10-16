#!/usr/bin/env python
"""Rank every cached diffuse configuration at a finer cell size under a fixed
companion set (Kev, 2026-10-15): the finer-binning, response-folded version of
the 10-deg ranking figure -- the signed 20.8 GeV significance against the
BIC-penalized background fit quality, one panel per Fermi bubbles description.

    FINE_TAG=psfuc python scripts/rank_fine.py --cell 2 --workers 6     # score rows not yet scored
    FINE_TAG=psfuc python scripts/rank_fine.py --cell 2 --summarize     # collect into one json + md

Rows: every row__<code>__<variant>.npz in the cell cache (the 2012 grid, A2020,
Galp21). Companions, as in the 10-deg figure (`plot_halo_significance.py`):
the bubbles template -- Totani's residual pair, the flat 2014 mask, or the
structured 2017 map -- plus the two Loop I shells, the uniform north/south
eROSITA lobes, isotropic emission and point sources. Bins 5.6-130 GeV (K = 6),
as in the final analysis. Scoring as `stage_select._score` without the
<sigma v> fit. Incremental: one JSON per (bubbles, row) under
results/fits/rank_fine/cells<c><tag>/<bubbles>/, skipped when present, so the
scoring can run while the row cache is still being built.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_select as ss  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine, stats  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-10-15"
BINS = ss.BINSETS["resp"]
BUBBLES = {"residual": "bubble_totani", "flat": "bubble_flat", "structured": "bubble_structured"}
COMPANIONS = ("loopI_shells", "erosita_unif_NS")


def out_dir(cell):
    d = ROOT / "results" / "fits" / "rank_fine" / f"cells{cell:g}{fine.tag()}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cached_rows(cell):
    return [tuple(p.stem[len("row__"):].split("__")) for p in sorted(fine.cache_dir(cell).glob("row__*.npz"))]


def _one(args):
    r, bub = args
    t0 = time.time()
    lib = ss._W["lib"]; keep, J = lib["keep"], lib["J"]; b = ss._W["bins"]
    st = ss._stack_with(r, [BUBBLES[bub], *COMPANIONS]); sub = ss._sub(st); Jb = J[b]
    sr, Jr = restrict_stack(sub, keep, Jb)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=True, n_refine=6)["per_bin"]
    ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(len(b))]; a = [pf[k]["signal_amplitude"] for k in range(len(b))]
    k20 = b.index(5); roi20 = float(Jr[k20][keep].sum())
    qa = {}
    for q, m in lib["quads"].items():
        stq, Jq = restrict_stack(st, m, J)
        one = {"counts": stq["counts"][5:6], "n_bins": 1, "templates": {n: t[5:6] for n, t in stq["templates"].items()}}
        qa[q] = fit_all_bins_profile(one, Jq[5:6], halo_mask=m, with_uncertainty=False)["per_bin"][0]["signal_amplitude"] * float(Jq[5][m].sum())
    return {"model": r[0], "variant": r[1], "bubbles": bub, "bins": list(b), "ts": ts, "a": a, "Q": stats.summarize(ts, a),
            "halo20": a[k20] * roi20, "halo20_lo": pf[k20]["signal_amplitude_lo"] * roi20, "halo20_hi": pf[k20]["signal_amplitude_hi"] * roi20,
            "quadrants20": qa, "favor": bool(ts[k20] > 4 and a[k20] > 0 and all(v > 0 for v in qa.values())),
            "neg2lnl_H0": float(sum(x["neg2lnl"] for x in bg)), "n_free_per_bin": len(st["templates"]),
            "n_cells": int(keep.sum()), "seconds": time.time() - t0}


def _default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(repr(o))


def score(cell, workers, limit):
    od = out_dir(cell)
    jobs = []
    for bub in BUBBLES:
        (od / bub).mkdir(exist_ok=True)
        for r in cached_rows(cell):
            if not (od / bub / f"{r[0]}__{r[1]}.json").exists():
                jobs.append((r, bub))
    if limit:
        jobs = jobs[:limit]
    print(f"{len(jobs)} (row, bubbles) jobs to score at {cell:g} deg{fine.tag()}; {workers} workers", flush=True)
    t0 = time.time()
    with Pool(workers, initializer=ss._init, initargs=(cell, BINS)) as pool:
        for i, d in enumerate(pool.imap_unordered(_one, jobs), 1):
            (od / d["bubbles"] / f"{d['model']}__{d['variant']}.json").write_text(json.dumps(d, default=_default))
            print(f"[{i}/{len(jobs)}] {d['bubbles']} {d['model']} {d['variant']}: halo20 {d['halo20']:.0f}, TS20 {d['ts'][d['bins'].index(5)]:.1f}, "
                  f"favor {d['favor']} ({d['seconds']:.0f}s; {time.time()-t0:.0f}s total)", flush=True)


def summarize(cell):
    od = out_dir(cell); K = len(BINS)
    out = {"date": DATE, "cell_deg": cell, "psf": fine.tag(), "bins": BINS, "companions": list(COMPANIONS), "panels": {}}
    L = [f"# Finer-binning ranking, {cell:g} deg cells{fine.tag()}, bins {BINS} (K = {K}) ({DATE})", "",
         "Companions: bubbles template + Loop I shells + uniform eROSITA N/S + isotropic + point sources. "
         "Margin = −(1/2)[(−2lnL_H0 + p ln n) − min over the panel], p = templates × K, n = cells × K.", ""]
    for bub in BUBBLES:
        rows = [json.loads(p.read_text()) for p in sorted((od / bub).glob("*.json"))]
        for r in rows:
            r["bic"] = r["neg2lnl_H0"] + r["n_free_per_bin"] * K * np.log(r["n_cells"] * K)
        best = min(r["bic"] for r in rows) if rows else 0.0
        for r in rows:
            r["margin"] = -0.5 * (r["bic"] - best)
            r["sig20"] = float(np.sign(r["a"][r["bins"].index(5)]) * np.sqrt(max(r["ts"][r["bins"].index(5)], 0.0)))
        rows.sort(key=lambda r: -r["margin"])
        out["panels"][bub] = rows
        fav = sum(r["favor"] for r in rows)
        L += [f"## {bub}: {len(rows)} rows, {fav} favor a halo", "",
              "| model | variant | margin | sig 20.8 | halo 20.8 [68%] | Z2018 | favor |", "|---|---|---|---|---|---|---|"]
        for r in rows:
            z = r["Q"]["Z_2018"]
            L.append(f"| {r['model']} | {r['variant']} | {r['margin']:.0f} | {r['sig20']:+.1f} | {r['halo20']:,.0f} [{r['halo20_lo']:,.0f}, {r['halo20_hi']:,.0f}] | {z:.1f} | {'Y' if r['favor'] else 'n'} |")
        L.append("")
    tag = f"cells{cell:g}{fine.tag()}_{DATE}"
    (ROOT / "results" / "fits" / f"rank_fine_{tag}.json").write_text(json.dumps(out, default=_default))
    (ROOT / "results" / f"rank_fine_{tag}.md").write_text("\n".join(L) + "\n")
    print("\n".join(l for l in L if l.startswith("## ")))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, required=True)
    ap.add_argument("--workers", type=int, default=6); ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--summarize", action="store_true")
    a = ap.parse_args()
    if a.summarize:
        summarize(a.cell)
    else:
        score(a.cell, a.workers, a.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
