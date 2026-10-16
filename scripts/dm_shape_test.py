#!/usr/bin/env python
"""Does an annihilation spectrum describe the halo template's bin amplitudes?
(Oscar, 2026-10-01: the 420 GeV W+W- benchmark of Yamashita's model, and the
500 GeV b-bbar fit used throughout, tested as SHAPES, not only as positive
amplitudes.)

    FINE_TAG=psfuc python scripts/dm_shape_test.py --cell 2 --workers 6

For each representative row and fixed bubble description (the configurations
of oscar_conditional.py, 2-deg PSF-folded), over the response-valid bins
(5.6-130 GeV, K = 6) and over all nine bins: the joint one-parameter fit of a
PPPC4DMID spectrum (<sigma v> signed, every background normalization free per
bin) against the fit with a free signed halo amplitude in every bin. The
difference Delta = -2lnL(spectrum) - sum_e -2lnL_e(free) is a likelihood-ratio
statistic for the spectral SHAPE with K - 1 degrees of freedom when the shape
is right (the one-parameter family is nested in the K-amplitude one), reported
with its chi-square p-value. Both fits use the same halo template, cells and
backgrounds; the free-amplitude fit is the profile of likelihood.fit_bin_profile
with n_refine = 6, as in the tables.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import oscar_conditional as OC  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import dmspectra, fine  # noqa: E402
from totani.likelihood import fit_all_bins_profile, restrict_stack  # noqa: E402

DATE = "2026-10-01"
EDGES = TOTANI_EDGES_GEV[:10]
SPECTRA = [("b", 500.0), ("W", 420.0), ("tau", 500.0)]
RANGES = {"resp": OC.RESP_BINS, "all": list(range(9))}


def sub_stack(st, J, bins):
    s = {"counts": st["counts"][bins], "n_bins": len(bins), "templates": {n: t[bins] for n, t in st["templates"].items()}}
    return s, J[bins], [EDGES[b] for b in bins] + [EDGES[bins[-1] + 1]]


def task(t):
    row_key, config = t
    lib = OC._W["lib"]; keep, J = lib["keep"], lib["J"]
    st = OC.stack_for(lib["counts"], row_key, config)
    out = {"row": list(row_key), "config": config}
    for rname, bins in RANGES.items():
        assert bins == list(range(bins[0], bins[-1] + 1))
        s, Js, edges = sub_stack(st, J, bins)
        sr, Jr = restrict_stack(s, keep, Js)
        pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False, n_refine=OC.N_REFINE)["per_bin"]
        v_free = float(sum(p["neg2lnl"] for p in pf))
        rec = {"v_free": v_free, "K": len(bins), "spectra": {}}
        for ch, m in SPECTRA:
            jf = dmspectra.joint_fit(s, Js, ch, m, edges, halo_mask=keep, n_scan=31, with_intervals=False)
            d = jf["neg2lnl_best"] - v_free
            rec["spectra"][f"{ch} {m:g}"] = {"sigmav_best": jf["sigmav_best"], "ts_vs_zero": jf["ts_vs_zero"], "delta": d,
                                             "p_shape": float(chi2.sf(max(d, 0.0), len(bins) - 1))}
        out[rname] = rec
    return out


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, default=2.0); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    t0 = time.time()
    tasks = [(r, c) for r in OC.ROWS for c in OC.CONFIGS]
    with Pool(a.workers, initializer=OC._init, initargs=(a.cell,)) as pool:
        res = []
        for k, x in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(x); print(f"[{k}/{len(tasks)}] {x['row']} {x['config']} ({time.time()-t0:.0f}s)", flush=True)
    tag = fine.tag()
    (ROOT / "results" / "fits" / f"dm_shape_test_cells{a.cell:g}{tag}_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res}, indent=1, default=float))
    L = [f"# Annihilation spectra as shapes: joint one-parameter fit against free bin amplitudes, {a.cell:g}° cells, PSF-folded ({DATE})", "",
         "Δ = −2lnL(one spectrum, ⟨σv⟩ signed, backgrounds free per bin) − Σ −2lnL(free signed halo amplitude per bin); χ²(K−1) p-value when the shape is right. "
         "⟨σv⟩ in 10⁻²⁵ cm³ s⁻¹; TS₀ = improvement of the spectrum over ⟨σv⟩ = 0. K = 6: 5.6–130 GeV; K = 9: all bins.", ""]
    for rname, label in (("resp", "5.6–130 GeV (K = 6, 5 dof)"), ("all", "all nine bins (K = 9, 8 dof)")):
        L += [f"## {label}", "", "| row | bubbles | " + " | ".join(f"{ch} {m:g} GeV: ⟨σv⟩, TS₀, Δ (p)" for ch, m in SPECTRA) + " |", "|---|---|" + "---|" * len(SPECTRA)]
        for r in OC.ROWS:
            for c in OC.CONFIGS:
                x = next(y for y in res if tuple(y["row"]) == r and y["config"] == c)[rname]
                cells = [f"{s['sigmav_best']*1e25:+.1f}, {s['ts_vs_zero']:.0f}, {s['delta']:.1f} ({s['p_shape']:.2g})" for s in x["spectra"].values()]
                L.append(f"| {r[0]} {r[1]} | {c} | " + " | ".join(cells) + " |")
        L.append("")
    out = ROOT / "results" / f"dm_shape_test_cells{a.cell:g}{tag}_{DATE}.md"
    out.write_text("\n".join(L) + "\n"); print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
