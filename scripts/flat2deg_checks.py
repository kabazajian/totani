#!/usr/bin/env python
"""Under the flat mask at 2 deg -- the template the data choose there -- is
the excess Totani's halo? Two tests, on the data and on flat-bubble skies
with and without an injected halo (the stage_e_2deg truths):

  (i)  the quadrant likelihood test at 20.8 GeV: one halo amplitude vs four,
       same background freedoms (standard division), Delta(-2lnL) ~ chi^2_3;
  (ii) the halo spectrum, 9 bins, exposure-corrected E^2 dN/dE (the
       conversion of halo_spectra_compare, which is cell-size independent).

Data rows: the true row, Totani's baseline (both splits), SS_Z10 5-ring,
A2020. Sims: the true row, 20 skies per truth.

    python scripts/flat2deg_checks.py --workers 6
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

from quadrant_lr_test import fit_joint  # noqa: E402
from stage_e_2deg import COMP, TRUE_ROW, stack_for, truth_cells  # noqa: E402
from totani import fine  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, fit_bin_profile, restrict_stack  # noqa: E402

CELL = 2.0
DATE = "2026-09-26"
K20 = 5
DATA_ROWS = [TRUE_ROW, ("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-")]
_W: dict = {}


def _init():
    _W["lib"] = fine.load_library(CELL)
    from totani.cells import cell_centers
    Lc, Bc = cell_centers(cell_deg=CELL)
    _W["lab"] = (Lc > 0).astype(int) + 2 * (Bc > 0).astype(int)


SPLIT_NS = "--ns" in sys.argv


def analyse(row_key, counts):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    row = fine.load_row(CELL, *row_key)
    st = stack_for(lib, row, "bubble_flat", counts)
    if SPLIT_NS:                                  # the flat bubbles with independent N and S normalizations
        from totani.cells import cell_centers
        _, Bc = cell_centers(cell_deg=CELL); north = (Bc > 0)[None]
        fb = st["templates"].pop("bubble_flat")
        st["templates"]["bubble_flat_N"] = fb * north; st["templates"]["bubble_flat_S"] = fb * ~north
    sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False)["per_bin"]
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    spec = [pf[k]["signal_amplitude"] * float(Jr[k][keep].sum()) for k in range(9)]
    ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(9)]
    T = np.stack([sr["templates"][n][K20] for n in names]); n = sr["counts"][K20]
    v_c, bg_amp, a_c, _ = fit_bin_profile(T, n, Jr[K20]); v_cl, bgc, ac = fit_joint(T, n, [Jr[K20]], bg_amp, [a_c]); v_c = min(v_c, v_cl)
    pieces = [Jr[K20] * (_W["lab"] == q) for q in range(4)]
    v_q, _, aq = fit_joint(T, n, pieces, bgc, [ac[0]] * 4)
    return {"spectrum": spec, "ts": ts, "lr_delta20": float(v_c - v_q), "quadrants20": [float(a * p[keep].sum()) for a, p in zip(aq, pieces)]}


def task(t):
    kind, i, mu, row_key = t
    lib = _W["lib"]; keep = lib["keep"]
    if mu is None:
        counts = lib["counts"]
    else:
        rng = np.random.default_rng([zlib.crc32(kind.encode()), i, 20260926])
        counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    return {"kind": kind, "sim": i, "row": list(row_key), **analyse(row_key, counts)}


def main() -> int:
    from scipy.stats import chi2
    sys.path.insert(0, str(ROOT / "scripts"))
    from halo_spectra_compare import conversion, e2f
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); ap.add_argument("--ns", action="store_true"); a = ap.parse_args()
    tasks = [("data", 0, None, r) for r in DATA_ROWS]
    for h in (False, True):
        mu = truth_cells("bubble_flat", h)
        tasks += [(f"{'halo' if h else 'null'}_flat", i, mu, TRUE_ROW) for i in range(20)]
    with Pool(a.workers, initializer=_init) as pool:
        res = pool.map(task, tasks)
    conv = conversion()
    out = {"date": DATE, "cell_deg": CELL, "results": res, "flux_per_count": conv.tolist()}
    (ROOT / "results" / "fits" / f"flat2deg_checks{fine.tag()}{'_ns' if SPLIT_NS else ''}_{DATE}.json").write_text(json.dumps(out, indent=1))
    E = [1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100]
    lines = [f"# Under the flat mask{' with N/S bubble normalizations' if SPLIT_NS else ''} at 2 degrees: sphericity and spectrum ({DATE})", "",
             f"Companions: {', '.join(COMP)}. Quadrant test: Δ(−2lnL) common vs four amplitudes at 20.8 GeV (χ²₃: median 2.4). "
             "Spectrum: exposure-corrected E²dN/dE (GeV cm⁻² s⁻¹) of the halo template, 9 bins.", "",
             "| sky | row | Δ quadrant @20.8 (p) | halo@20.8 | quadrants W-S E-S W-N E-N | E²dN/dE by bin (×10⁻⁷) |", "|---|---|---|---|---|---|"]
    for r in res:
        if r["kind"] == "data":
            sp = e2f(r["spectrum"], conv) * 1e7
            lines.append(f"| data | {r['row'][0]} {r['row'][1]} | {r['lr_delta20']:.1f} ({chi2.sf(max(r['lr_delta20'],0),3):.2g}) | {r['spectrum'][K20]:,.0f} | "
                         + " ".join(f"{x:,.0f}" for x in r["quadrants20"]) + " | " + " ".join(f"{x:+.2f}" for x in sp) + " |")
    for kind in ("null_flat", "halo_flat"):
        S = [r for r in res if r["kind"] == kind]
        d = np.array([r["lr_delta20"] for r in S]); sp = e2f(np.median([r["spectrum"] for r in S], axis=0), conv) * 1e7
        qm = np.median([r["quadrants20"] for r in S], axis=0)
        lines.append(f"| {kind} (20) | true row, median | {np.median(d):.1f} [{np.percentile(d,16):.1f}, {np.percentile(d,84):.1f}] | {np.median([r['spectrum'][K20] for r in S]):,.0f} | "
                     + " ".join(f"{x:,.0f}" for x in qm) + " | " + " ".join(f"{x:+.2f}" for x in sp) + " |")
    inj = e2f(np.clip(json.loads((ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json").read_text())["a0_halo_counts"], 0, None), conv) * 1e7
    lines += ["", "Injected (Totani's spectrum) E²dN/dE ×10⁻⁷: " + " ".join(f"{x:+.2f}" for x in inj)]
    (ROOT / "results" / f"flat2deg_checks{fine.tag()}{'_ns' if SPLIT_NS else ''}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
