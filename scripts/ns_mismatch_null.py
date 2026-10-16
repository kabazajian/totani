#!/usr/bin/env python
"""Item 3 (Kev, 2026-09-27): does any tested foreground mismatch, on a sky with
NO halo, make the flat-mask halo fit return the data's south-heavy excess
(S/N ~ 1.4 at ~6000 counts at 20.8 GeV)? PSF-folded 2-deg cells.

    FINE_TAG=psf python scripts/ns_mismatch_null.py --workers 6

Truths (all halo-free, true row SO_Z6_R30_T150_C5 1-gas, background-only fits
to the data, Poisson-sampled on |b|>=10 cells):
  flat          flat mask + companions (the reference: fitted correctly)
  structured    structured map + companions
  construction  Totani's residual pair of the true row + companions (the
                'third morphology': bubbles shaped like his construction)
  flat_NS       flat mask with independent N/S normalizations + companions
                (bubbles brighter in one hemisphere, as the data prefer)
Fits (all with a single-normalization flat mask, halo split N/S at b = 0,
20.8 GeV):
  correct       the true row and companions
  wrong_row     eight other diffuse rows of the Stage E subset
  no_erosita    the true row without the eROSITA lobes
  shells        the true row with Loop I shells in place of Haslam
Reported per (truth, fit): the halo's N and S counts, S/N, and the fraction of
skies with total >= 4000 counts and S/N >= 1.3. 20 skies per truth.
"""
from __future__ import annotations

import argparse
import json
import sys
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from ns_halo_spectra import joint  # noqa: E402
from stage_e_2deg import COMP, TRUE_ROW  # noqa: E402
from totani import fine  # noqa: E402
from totani.cells import cell_centers  # noqa: E402
from totani.likelihood import fit_all_bins, restrict_stack  # noqa: E402

CELL = 2.0
DATE = "2026-09-27"
K20 = 5
WRONG_ROWS = [(c, "gas1_ics1") for c in stage_e.fit_codes()[::2]]
_W: dict = {}


def _init():
    _W["lib"] = fine.load_library(CELL)
    _, Bc = cell_centers(cell_deg=CELL); _W["north"] = Bc > 0


def templates_for(lib, row, bubble_names, companions):
    t = {"isotropic": lib["templates"]["isotropic"], "point_sources": lib["templates"]["point_sources"], **row["diffuse"]}
    for n in bubble_names:
        t[n] = row["residual"][n] if n.startswith("bubble_totani_") else lib["templates"][n]
    for n in companions:
        t[n] = lib["templates"][n]
    return t


def truth_mu(kind):
    lib = fine.load_library(CELL); row = fine.load_row(CELL, *TRUE_ROW)
    _, Bc = cell_centers(cell_deg=CELL); north = (Bc > 0)[None]
    if kind == "flat":
        t = templates_for(lib, row, ("bubble_flat",), COMP)
    elif kind == "structured":
        t = templates_for(lib, row, ("bubble_structured",), COMP)
    elif kind == "construction":
        t = templates_for(lib, row, ("bubble_totani_pos", "bubble_totani_neg"), COMP)
        # the negative image's sign per bin as the data choose it
        st = {"counts": lib["counts"], "n_bins": 9, "templates": t}
        vp = fine.bg_neg2lnl(st, lib["keep"]); t2 = dict(t); t2["bubble_totani_neg"] = -t["bubble_totani_neg"]
        vm = fine.bg_neg2lnl({**st, "templates": t2}, lib["keep"]); sg = np.where(vp <= vm, 1, -1)
        t["bubble_totani_neg"] = t["bubble_totani_neg"] * sg[:, None, None]
    elif kind == "flat_NS":
        t = templates_for(lib, row, (), COMP); fb = lib["templates"]["bubble_flat"]
        t["bubble_flat_N"] = fb * north; t["bubble_flat_S"] = fb * ~north
    st = {"counts": lib["counts"], "n_bins": 9, "templates": t}
    fit = fit_all_bins(restrict_stack(st, lib["keep"]), halo_mask=lib["keep"])
    names = fit["component_order"]; amps = np.array([[fit["per_bin"][k]["amplitudes"][n] for n in names] for k in range(9)])
    return np.einsum("kn,nkij->kij", amps, np.stack([t[n] for n in names]))


def halo_ns(counts, row_key, companions):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]; north = _W["north"]
    t = templates_for(lib, fine.load_row(CELL, *row_key), ("bubble_flat",), companions)
    st = {"counts": counts, "n_bins": 9, "templates": t}
    sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    T = np.stack([sr["templates"][n][K20] for n in names]); n = sr["counts"][K20]
    hN, hS = Jr[K20] * north, Jr[K20] * ~north
    _, a = joint(T, n, [hN, hS])
    return float(a[0] * hN[keep].sum()), float(a[1] * hS[keep].sum())


FITS = {"correct": ([TRUE_ROW], COMP), "wrong_row": (WRONG_ROWS, COMP),
        "no_erosita": ([TRUE_ROW], ("loopI_haslam_raw", "sun_moon_ic")),
        "shells": ([TRUE_ROW], ("erosita_unif_N", "erosita_unif_S", "loopI_shell1", "loopI_shell2", "sun_moon_ic"))}


def task(t):
    truth, i, mu = t
    lib = _W["lib"]; keep = lib["keep"]
    if mu is None:
        counts = lib["counts"]
    else:
        rng = np.random.default_rng([zlib.crc32(truth.encode()), i, 20260927])
        counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    rec = {"truth": truth, "sim": i}
    for fname, (rows, comp) in FITS.items():
        rec[fname] = [dict(zip(("N", "S"), halo_ns(counts, r, comp)), row=list(r)) for r in rows]
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    tasks = [("data", 0, None)]
    for kind in ("flat", "structured", "construction", "flat_NS"):
        mu = truth_mu(kind); tasks += [(kind, i, mu) for i in range(20)]
    with Pool(a.workers, initializer=_init) as pool:
        res = pool.map(task, tasks)
    (ROOT / "results" / "fits" / f"ns_mismatch_null{fine.tag()}_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1))
    lines = [f"# Which halo-free foreground mismatch makes a south-heavy flat-mask halo? 2° cells{' (PSF-folded)' if fine.tag() else ''} ({DATE})", "",
             "Halo at 20.8 GeV split at b = 0 under a single-normalization flat mask. 'data-like' = total ≥ 4000 counts and S/N ≥ 1.3.", "",
             "| truth → fit | N (median) | S (median) | total (median) | S/N (median) | data-like fraction |", "|---|---|---|---|---|---|"]
    for truth in ("data", "flat", "structured", "construction", "flat_NS"):
        S = [r for r in res if r["truth"] == truth]
        for fname in FITS:
            vals = [(x["N"], x["S"]) for r in S for x in r[fname]]
            Nn = np.array([v[0] for v in vals]); Ss = np.array([v[1] for v in vals]); tot = Nn + Ss
            ratio = np.where(np.abs(Nn) > 50, Ss / np.where(np.abs(Nn) > 50, Nn, 1), np.nan)
            like = float(np.mean((tot >= 4000) & (ratio >= 1.3)))
            lines.append(f"| {truth} → {fname} | {np.median(Nn):,.0f} | {np.median(Ss):,.0f} | {np.median(tot):,.0f} | {np.nanmedian(ratio):.2f} | {like:.2f} |")
    (ROOT / "results" / f"ns_mismatch_null{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
