#!/usr/bin/env python
"""Oscar's note of 2026-09-29, item 2: where the halo information comes from
at finer resolution -- the halo fitted only in cells OUTSIDE the flat-bubble
footprint and only INSIDE it, under each bubble description, with an injection
sensitivity check in each region.

    FINE_TAG=psf python scripts/inout_bubbles.py --cell 2 --workers 5
    FINE_TAG=psfuc python scripts/inout_bubbles.py --cell 1 --common   # Oscar 2026-10-01: one mask at both resolutions

Regions (fixed by geometry, not by result): 'outside' = fitted cells
(|b| >= 10, extended sources masked) whose UNCONVOLVED flat-mask template is
zero in every bin; 'inside' = the other fitted cells. With --common the
region is defined once, on the 2-deg cells, and carried to the finer cells
(each 2-deg cell is exactly four 1-deg cells), so both resolutions fit the
same pixels and the same fraction of the NFW template (Oscar, 2026-10-01).
Configurations: structured
map, flat mask, construction, each with the 2-deg companions (eROSITA uniform
N/S, Haslam, solar IC), fitted in the region only (backgrounds refitted there).
Per region: the halo counts at 20.8 GeV with 68% interval (in counts of the
halo template within that region), per-bin TS, Q^+ (K = 9 and the
response-valid K = 6) with p and Z. Injection: a spherical halo on Totani's
spectrum at 0, 1, 2 times his amplitude added to the counts; recovered/injected
slope per bin with the templates held fixed (the construction's own absorption
of an injected halo is item 3, injection_refined.py).
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

import oscar_conditional as OC  # noqa: E402
from totani import fine, ranking, simsky, stats  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

DATE = "2026-09-29"
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("a2020", "-"),
        ("SS_Z10_R30_T150_C5", "gas5_ics3")]      # representative five-ring row added 2026-10-15 (audit)
_W: dict = {}


def _init(cell, common=False):
    OC._init(cell); _W["cell"] = cell; _W["lib"] = OC._W["lib"]
    keep = _W["lib"]["keep"]
    if common:
        plain = np.load(ranking.CACHE / "cells2" / "library.npz")
        out2 = np.abs(plain["T__bubble_flat"]).sum(0) == 0
        f = int(round(2.0 / cell))
        out = keep & np.kron(out2, np.ones((f, f), bool))
    else:
        plain = np.load(ranking.CACHE / f"cells{cell:g}" / "library.npz")
        out = keep & (np.abs(plain["T__bubble_flat"]).sum(0) == 0)
    _W["regions"] = {"outside": out, "inside": keep & ~out}
    J = _W["lib"]["J"][OC.fine.K20]
    _W["J_frac"] = {g: float(J[m].sum() / J[keep].sum()) for g, m in _W["regions"].items()}
    _W["amp"] = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")


def fit_region(counts, row_key, config, region, with_intervals):
    lib = _W["lib"]; J = lib["J"]; m = _W["regions"][region]
    st = OC.stack_for(counts, row_key, config)
    sr, Jr = restrict_stack(st, m, J)
    bg = fit_all_bins(sr, halo_mask=m)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=m, with_uncertainty=with_intervals, n_refine=OC.N_REFINE)["per_bin"]
    roi = [float(Jr[k][m].sum()) for k in range(9)]
    ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(9)]; a = [pf[k]["signal_amplitude"] for k in range(9)]
    r = {"a": a, "counts": [a[k] * roi[k] for k in range(9)], "ts": ts, "roiJ": roi, "all": stats.summarize(ts, a), "resp": stats.summarize(ts, a, OC.RESP_BINS)}
    if with_intervals:
        r["lo20"] = pf[5]["signal_amplitude_lo"] * roi[5]; r["hi20"] = pf[5]["signal_amplitude_hi"] * roi[5]
    return r


def task(t):
    row_key, config, region = t
    lib = _W["lib"]
    data = fit_region(lib["counts"], row_key, config, region, True)
    inj = []
    for lev in (1.0, 2.0):
        counts = lib["counts"] + lev * _W["amp"][:, None, None] * lib["J"]
        inj.append(fit_region(counts, row_key, config, region, False)["a"])
    rec = np.array([data["a"]] + inj); injected = np.array([[lev * _W["amp"][k] for k in range(9)] for lev in (0.0, 1.0, 2.0)])
    slopes = []
    for k in range(9):
        d = injected[:, k] - injected[:, k].mean()
        slopes.append(float((d * (rec[:, k] - rec[:, k].mean())).sum() / (d ** 2).sum()) if np.any(d != 0) else None)
    return {"row": list(row_key), "config": config, "region": region, **data, "injection_slopes": slopes,
            "n_cells": int(_W["regions"][region].sum()), "J_frac": _W["J_frac"][region]}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, default=2.0); ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--common", action="store_true", help="region defined on the 2-deg cells and carried to this cell size"); a = ap.parse_args()
    global DATE
    if a.common:
        DATE = "2026-10-01"
    t0 = time.time()
    tasks = [(r, c, g) for r in ROWS for c in OC.CONFIGS for g in ("outside", "inside")]
    with Pool(a.workers, initializer=_init, initargs=(a.cell, a.common)) as pool:
        res = []
        for k, x in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(x); print(f"[{k}/{len(tasks)}] {x['row']} {x['config']} {x['region']}: {x['counts'][5]:,.0f} ({time.time()-t0:.0f}s)", flush=True)
    tag = fine.tag() + ("_common" if a.common else "")
    (ROOT / "results" / "fits" / f"inout_bubbles_cells{a.cell:g}{tag}_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1, default=float))
    n_out = res[0]["n_cells"] if res[0]["region"] == "outside" else next(x["n_cells"] for x in res if x["region"] == "outside")
    n_in = next(x["n_cells"] for x in res if x["region"] == "inside")
    jf = {g: next(x["J_frac"] for x in res if x["region"] == g) for g in ("outside", "inside")}
    lines = [f"# Halo inside vs outside the flat-bubble footprint, {a.cell:g}° cells{' (PSF-folded)' if fine.tag() else ''}{' — common 2° mask' if a.common else ''} ({DATE})", "",
             (f"Outside: {n_out} fitted cells ({jf['outside']*100:.1f}% of the NFW template at 20.8 GeV)" + (" inside the 2° cells with zero unconvolved flat-mask template" if a.common else " with zero unconvolved flat-mask template")
              + f"; inside: {n_in} ({jf['inside']*100:.1f}%). Counts are the halo template's own counts in the region. "
              "Injection slopes: recovered/injected per bin for a spherical halo on Totani's spectrum (templates fixed); '—' where nothing was injected."), "",
             "| row | bubbles | region | halo 20.8 GeV [68%] | TS 20.8 | Q⁺ K=9 (Z₂₀₁₈) | Q⁺ K=6 (Z₂₀₁₈) | injection slope 4.3 → 100 GeV |", "|---|---|---|---|---|---|---|---|"]
    for r in ROWS:
        for c in OC.CONFIGS:
            for g in ("outside", "inside"):
                x = next(y for y in res if tuple(y["row"]) == r and y["config"] == c and y["region"] == g)
                sl = " ".join("—" if s is None else f"{s:.2f}" for s in x["injection_slopes"][2:])
                lines.append(f"| {r[0]} {r[1]} | {c} | {g} | {x['counts'][5]:,.0f} [{x['lo20']:,.0f}, {x['hi20']:,.0f}] | {x['ts'][5]:.1f} | "
                             f"{x['all']['Q_plus']:.1f} ({x['all']['Z_2018']:.1f}) | {x['resp']['Q_plus']:.1f} ({x['resp']['Z_2018']:.1f}) | {sl} |")
    out = ROOT / "results" / f"inout_bubbles_cells{a.cell:g}{tag}_{DATE}.md"
    out.write_text("\n".join(lines) + "\n"); print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
