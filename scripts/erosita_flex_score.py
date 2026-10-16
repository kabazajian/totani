#!/usr/bin/env python
"""fine.score (TS, 20.8 GeV halo, quadrant gate, joint b-bbar 500 GeV <sigma v>) for the
34 data rows under the uniform N/S lobes and the flexible 12-piece brightness lobes,
flat mask, Haslam, solar IC. PSF-folded 2-deg cells.

    FINE_TAG=psf python scripts/erosita_flex_score.py --workers 4
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
import stage_e  # noqa: E402
from totani import fine  # noqa: E402

MODELS = ("unif", "bright_flex")


def task(r):
    lib = EF._W["lib"]; out = {"row": list(r)}
    for m in MODELS:
        d = fine.score(EF.stack(lib["counts"], r, EF._W["models"][m]), lib, with_sigmav=True, edges_gev=EF.EDGES)
        out[m] = {"ts": d["ts"], "halo20": d["halo_counts"][5], "favor": d["favor"], "sigmav": d["sigmav"]["sigmav_best"],
                  "sigmav_ts": d["sigmav"]["ts_vs_zero"]}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4); a = ap.parse_args()
    rows = [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()] + [("a2020", "-"), ("SL_Z6_R30_T150_C2", "gas1_ics1")]
    with Pool(a.workers, initializer=EF._init) as pool:
        res = pool.map(task, rows)
    (ROOT / "results" / "fits" / f"erosita_flex_score{fine.tag()}_{EF.DATE}.json").write_text(json.dumps(res, indent=1, default=float))
    lines = [f"# Scoring under uniform vs flexible brightness lobes, PSF-folded 2° ({EF.DATE})", "",
             "| lobes | split | favor | TS (sum) median | halo20 median [range] | <σv> b b̄ 500 GeV ×1e-25 median [range] | 2σ+ |", "|---|---|---|---|---|---|---|"]
    for m in MODELS:
        for split in ("gas1_ics1", "gas5_ics3", "all"):
            S = [r[m] for r in res if split == "all" or r["row"][1] == split]
            ts = [np.sum(x["ts"]) if np.ndim(x["ts"]) else x["ts"] for x in S]; h = [x["halo20"] for x in S]; sv = [x["sigmav"] * 1e25 for x in S]
            lines.append(f"| {m} | {split} | {sum(bool(x['favor']) for x in S)}/{len(S)} | {np.median(ts):.0f} | {np.median(h):,.0f} [{min(h):,.0f}–{max(h):,.0f}] | "
                         f"{np.median(sv):+.1f} [{min(sv):+.1f}, {max(sv):+.1f}] | {sum(x['sigmav_ts'] > 4 and x['sigmav'] > 0 for x in S)}/{len(S)} |")
    b = next(r for r in res if r["row"][0] == "SL_Z6_R30_T150_C2")
    lines += ["", "Totani's baseline row: " + "; ".join(f"{m}: halo20 {b[m]['halo20']:,.0f}, <σv> {b[m]['sigmav']*1e25:+.1f}e-25, favor {b[m]['favor']}" for m in MODELS)]
    (ROOT / "results" / f"erosita_flex_score{fine.tag()}_{EF.DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
