#!/usr/bin/env python
"""Does a lobe mismatch make the data's south-heavy (east-south) halo?
Skies whose TRUE eROSITA lobes are the best-fitting flexible brightness
lobes (12 pieces, fitted to the data background-only with the true row
SO_Z6_R30_T150_C5 1-gas, flat mask, Haslam, solar IC), with no halo and with
a spherical halo on Totani's spectrum; each fitted under every lobe model of
erosita_flex_ns.py (the 2-deg default is uniform N/S). PSF-folded 2-deg cells.

    FINE_TAG=psf python scripts/erosita_flex_truth.py --workers 6
"""
from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import erosita_flex_ns as EF  # noqa: E402
from stage_e_2deg import TRUE_ROW  # noqa: E402
from totani import fine, simsky  # noqa: E402
from totani.likelihood import fit_all_bins, restrict_stack  # noqa: E402

DATE = EF.DATE


def truth_mu(with_halo):
    EF._init(); lib = EF._W["lib"]
    st = EF.stack(lib["counts"], TRUE_ROW, EF._W["models"]["bright_flex"])
    fit = fit_all_bins(restrict_stack(st, lib["keep"]), halo_mask=lib["keep"])
    names = fit["component_order"]
    amps = np.array([[fit["per_bin"][k]["amplitudes"][n] for n in names] for k in range(9)])
    mu = np.einsum("kn,nkij->kij", amps, np.stack([st["templates"][n] for n in names]))
    if with_halo:
        amp = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")
        mu = mu + amp[:, None, None] * lib["J"]
    return mu


def task(t):
    kind, i, mu = t
    lib = EF._W["lib"]; keep = lib["keep"]
    rng = np.random.default_rng([11, i, 20260927])
    counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    return {"kind": kind, "sim": i, "models": EF.analyse(counts, TRUE_ROW, with_h0=False)}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    tasks = []
    for h in (False, True):
        mu = truth_mu(h); tasks += [(f"{'halo' if h else 'null'}_brightflex", i, mu) for i in range(20)]
    with Pool(a.workers, initializer=EF._init) as pool:
        res = pool.map(task, tasks)
    (ROOT / "results" / "fits" / f"erosita_flex_truth{fine.tag()}_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1))
    lines = [f"# Skies with the data's best-fitting flexible brightness lobes as truth ({DATE}, PSF-folded 2°)", "",
             "True row SO_Z6_R30_T150_C5 1-gas; halo at 20.8 GeV in four quadrant pieces under each fitted lobe model. Medians [16–84%] over 20 skies.", "",
             "| sky | fitted lobes | halo total | S/N | S/N West | S/N East | WN / WS / EN / ES |", "|---|---|---|---|---|---|---|"]
    for kind in ("null_brightflex", "halo_brightflex"):
        S = [r for r in res if r["kind"] == kind]
        for m in EF.MODEL_ORDER:
            H = [r["models"][m]["halo"] for r in S]; rt = np.array([EF.ratios(h) for h in H]); tot = np.array([sum(h.values()) for h in H])
            f = lambda c: f"{np.nanmedian(rt[:,c]):.2f} [{np.nanpercentile(rt[:,c],16):.2f}, {np.nanpercentile(rt[:,c],84):.2f}]"
            q = " / ".join(f"{np.median([h[k] for h in H]):,.0f}" for k in ("WN", "WS", "EN", "ES"))
            lines.append(f"| {kind} | {m} | {np.median(tot):,.0f} | {f(0)} | {f(1)} | {f(2)} | {q} |")
    (ROOT / "results" / f"erosita_flex_truth{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
