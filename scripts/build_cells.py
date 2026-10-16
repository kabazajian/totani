#!/usr/bin/env python
"""Build a cell library and row caches at any cell size, unconvolved and
PSF-folded (Oscar, 2026-09-29, item 1: the 2 deg vs 1 deg comparison with a
response check). One process, <= 10 GB; never alongside another stack build.

    python scripts/build_cells.py --cell 1              # unconvolved library + rows, then PSF-folded
    python scripts/build_cells.py --cell 1 --rows rep   # only the representative rows
    python scripts/build_cells.py --cell 2 --rows grid  # the whole 2012 grid (~1 min per row, existing rows skipped)

Rows: 'rep' = Totani's baseline (both splits), the simulation truth rows,
A2020 and Galp21 SA0_R12 (1-gas); 'subset' = rep + the 33-row stratified
2012 subset (the selection ensemble). Galp21 beyond SA0_R12 waits on
Youyou's input checks (Oscar, 09-29).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
REP = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("SO_Z6_R30_T150_C5", "gas1_ics1"),
       ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-"), ("G21_SA0_R12", "gas1_ics1")]


def rows(which, cell=None):
    import stage_e
    if which == "existing":                       # every row in the unconvolved cache at this cell size
        from totani import ranking
        d = ranking.CACHE / f"cells{cell:g}"
        return [tuple(p.stem[len("row__"):].split("__")) for p in sorted(d.glob("row__*.npz"))]
    from totani import ranking
    out = [r for r in REP if ranking.available(*r)]       # the Galp21 row only when its cubes are ingested
    if which == "subset":
        out += [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()]
    if which == "grid":                           # the whole 2012 grid (Kev, 2026-10-15: the finer-binning ranking figure)
        from totani import ranking
        out += [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()]
        out += [(c, v) for c in ranking.MODEL_CODES for v in ("gas1_ics1", "gas5_ics3")]
    if which == "galp21":                         # the 24 Galp21 rows (six models, both IC treatments, both splits)
        out += [(f"G21_{m}{ic}", v) for m in ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12")
                for ic in ("", "+iso") for v in ("gas1_ics1", "gas5_ics3")]
    return list(dict.fromkeys(out))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, required=True)
    ap.add_argument("--rows", default="subset", choices=("rep", "subset", "grid", "galp21", "existing")); ap.add_argument("--psf-only", action="store_true")
    ap.add_argument("--_psf_child", action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args()
    t0 = time.time()
    if not a._psf_child and not a.psf_only:
        from totani import fine
        fine.build_library(a.cell, EDGES); print(f"library {a.cell:g} deg ({time.time()-t0:.0f}s)", flush=True)
        R = rows(a.rows, a.cell)
        for i, r in enumerate(R, 1):
            fine.build_row(a.cell, r[0], r[1], EDGES); print(f"[{i}/{len(R)}] row {r} ({time.time()-t0:.0f}s)", flush=True)
        # PSF pass in a fresh process (FINE_TAG is read at call time, but keep the pixel caches' memory separate)
        env = dict(os.environ, FINE_TAG=os.environ.get("PSF_TAG", "psf"))
        return subprocess.call([sys.executable, __file__, "--cell", str(a.cell), "--rows", a.rows, "--_psf_child"], env=env)
    if a.psf_only:
        os.environ["FINE_TAG"] = os.environ.get("PSF_TAG", "psf")
    from totani import fine
    assert fine.tag().startswith("_psf")
    fine.build_library_psf(a.cell, EDGES); print(f"PSF library {a.cell:g} deg ({time.time()-t0:.0f}s)", flush=True)
    R = rows(a.rows, a.cell)
    for i, r in enumerate(R, 1):
        fine.build_row_psf(a.cell, r[0], r[1], EDGES); print(f"[{i}/{len(R)}] PSF row {r} ({time.time()-t0:.0f}s)", flush=True)
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
