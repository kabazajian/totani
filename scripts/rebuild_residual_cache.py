#!/usr/bin/env python
"""Rebuild every cached Totani-residual template with the current fitter.

    python scripts/rebuild_residual_cache.py        # ~75 min, one process, ~8 GB

Needed after any change to `likelihood.fit_bin` (2026-09-22: the linear-space
polish moved the 5-ring rows' templates by 1-4%). One row at a time: the
pixel-level build reads the GALPROP cubes for each row.
"""
from __future__ import annotations
import sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from reduce_totani import TOTANI_EDGES_GEV
from totani import ranking, selection
EDGES = TOTANI_EDGES_GEV[:10]
rows = sorted({tuple(p.stem.split("__")[1:]) for p in ranking.CACHE.glob("totani_residual__*.npz")})
t0 = time.time()
for i, (code, v) in enumerate(rows, 1):
    try:
        selection.totani_residual_cells(code, v, EDGES, force=True)
        print(f"[{i}/{len(rows)}] {code} {v} ({time.time()-t0:.0f}s)", flush=True)
    except Exception as e:
        print(f"[{i}/{len(rows)}] {code} {v} FAILED: {e}", flush=True)
