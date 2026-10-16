#!/usr/bin/env python
"""The 2-deg simulated skies for the vote that flipped: does a sky whose
bubbles ARE the structured map choose the flat mask at 2 deg, as the data do?

    python scripts/stage_e_2deg.py --workers 6

Truths (2-deg cells): the true row SO_Z6_R30_T150_C5 (single-gas) with the 2-deg
no-residual companions (eROSITA uniform N/S, Haslam raw, solar IC) and EITHER
the structured map OR the flat mask, fitted to the data background-only and
Poisson-sampled on the |b|>=10 cells (disk cells keep the real counts); each
with no halo and with a halo on Totani's spectrum (6412 counts at 20.8 GeV).
The two truths differ only in the bubble template.

Per sky, each of the 33 fit rows (16 stratified 2012 models x 2 splits +
A2020, truth model excluded) fits the flat-mask set and the structured-map set
with the same companions: -2lnL(H0) of each (the vote) and the 20.8 GeV halo
(profile, one bin) with its TS under each. 'Sky 0' is the data itself.
"""
from __future__ import annotations

import argparse
import json
import sys
import zlib
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from totani import fine, simsky  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

CELL = 2.0
DATE = "2026-09-26"
K20 = 5
TRUE_ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
COMP = ("erosita_unif_N", "erosita_unif_S", "loopI_haslam_raw", "sun_moon_ic")
N_SKIES = 20
_W: dict = {}


def _init():
    _W["lib"] = fine.load_library(CELL)


def voters():
    return [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()] + [("a2020", "-")]


def stack_for(lib, row, bubble, counts):
    st = fine.base_stack(lib, row, counts)
    for n in (bubble,) + COMP:
        st["templates"][n] = lib["templates"][n]
    return st


def halo20(st, lib):
    keep, J = lib["keep"], lib["J"]
    sr, Jr = restrict_stack(st, keep, J)
    one = {"counts": sr["counts"][K20:K20 + 1], "n_bins": 1, "templates": {n: t[K20:K20 + 1] for n, t in sr["templates"].items()}}
    b = fit_all_bins(one, halo_mask=keep)["per_bin"][0]["neg2lnl"]
    p = fit_all_bins_profile(one, Jr[K20:K20 + 1], halo_mask=keep, with_uncertainty=False)["per_bin"][0]
    return float(p["signal_amplitude"] * Jr[K20][keep].sum()), float(b - p["neg2lnl"])


def truth_cells(bubble, with_halo):
    lib = fine.load_library(CELL); row = fine.load_row(CELL, *TRUE_ROW)
    st = stack_for(lib, row, bubble, lib["counts"])
    fit = fit_all_bins(restrict_stack(st, lib["keep"]), halo_mask=lib["keep"])
    names = fit["component_order"]
    amps = np.array([[fit["per_bin"][k]["amplitudes"][n] for n in names] for k in range(9)])
    mu = np.einsum("kn,nkij->kij", amps, np.stack([st["templates"][n] for n in names]))
    if with_halo:
        amp = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")
        mu = mu + amp[:, None, None] * lib["J"]
    return mu


def run_sky(task):
    truth, i, mu = task
    lib = _W["lib"]; keep = lib["keep"]
    if mu is None:
        counts = lib["counts"]
    else:
        rng = np.random.default_rng([zlib.crc32(truth.encode()), i, 20260926])
        counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    rows = []
    for r in voters():
        row = fine.load_row(CELL, *r)
        rec = {"model": r[0], "variant": r[1]}
        for b in ("bubble_flat", "bubble_structured"):
            st = stack_for(lib, row, b, counts)
            h0 = float(fine.bg_neg2lnl(st, keep).sum()); h, ts = halo20(st, lib)
            rec[b] = {"neg2lnl_H0": h0, "halo20": h, "ts20": ts}
        rec["flat_minus_structured"] = rec["bubble_flat"]["neg2lnl_H0"] - rec["bubble_structured"]["neg2lnl_H0"]
        rows.append(rec)
    return {"truth": truth, "sim": i, "rows": rows}


def summarise(out):
    lines = [f"# 2-degree simulated skies: the flat-vs-structured vote ({DATE})", "",
             f"Truth row {TRUE_ROW[0]} (1-gas), companions {', '.join(COMP)}; 33 fit rows vote by −2lnL(H0) with the same companions. "
             "'flat wins' = fraction of rows whose flat-mask fit is better; halo at 20.8 GeV under the chosen template of each row.", "",
             "| sky | n skies | flat wins (median over skies, range) | halo@20.8 under chosen template (median) | TS (median) | halo under flat / structured (median) |",
             "|---|---|---|---|---|---|"]
    by = {}
    for s in out:
        by.setdefault(s["truth"], []).append(s)
    for t, S in by.items():
        fw = [np.mean([r["flat_minus_structured"] < 0 for r in s["rows"]]) for s in S]
        chosen = [np.median([r["bubble_flat" if r["flat_minus_structured"] < 0 else "bubble_structured"]["halo20"] for r in s["rows"]]) for s in S]
        tsc = [np.median([r["bubble_flat" if r["flat_minus_structured"] < 0 else "bubble_structured"]["ts20"] for r in s["rows"]]) for s in S]
        hf = [np.median([r["bubble_flat"]["halo20"] for r in s["rows"]]) for s in S]
        hs = [np.median([r["bubble_structured"]["halo20"] for r in s["rows"]]) for s in S]
        lines.append(f"| {t} | {len(S)} | {np.median(fw):.2f} [{min(fw):.2f}, {max(fw):.2f}] | {np.median(chosen):,.0f} | {np.median(tsc):.0f} | {np.median(hf):,.0f} / {np.median(hs):,.0f} |")
    (ROOT / "results" / f"stage_e_2deg{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    t0 = time.time()
    tasks = [("data", 0, None)]
    for b, lab in (("bubble_structured", "structured"), ("bubble_flat", "flat")):
        for h in (False, True):
            mu = truth_cells(b, h)
            tasks += [(f"{'halo' if h else 'null'}_{lab}", i, mu) for i in range(N_SKIES)]
    print(f"{len(tasks)} skies ({time.time()-t0:.0f}s)", flush=True)
    out = []
    with Pool(a.workers, initializer=_init) as pool:
        for k, s in enumerate(pool.imap_unordered(run_sky, tasks), 1):
            out.append(s)
            fw = np.mean([r["flat_minus_structured"] < 0 for r in s["rows"]])
            print(f"[{k}/{len(tasks)}] {s['truth']} {s['sim']}: flat wins {fw:.2f} ({time.time()-t0:.0f}s)", flush=True)
            if k % 10 == 0 or k == len(tasks):
                (ROOT / "results" / "fits" / f"stage_e_2deg{fine.tag()}_{DATE}.json").write_text(json.dumps(out))
    summarise(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
