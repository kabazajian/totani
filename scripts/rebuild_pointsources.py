#!/usr/bin/env python
"""Replace the point-source template of an existing PSF-folded library by the
exact-position version (fine.point_sources_cells, 2026-10-01), keeping the old
library as library_pixelsrc.npz.

    FINE_TAG=psfuc python scripts/rebuild_pointsources.py --cell 2
"""
from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, required=True); a = ap.parse_args()
    from reduce_totani import TOTANI_EDGES_GEV
    from totani import fine
    assert fine.tag().startswith("_psf")
    edges = TOTANI_EDGES_GEV[:10]
    p = fine.cache_dir(a.cell) / "library.npz"; bak = p.with_name("library_pixelsrc.npz")
    if not bak.exists():
        shutil.copy2(p, bak)
    t0 = time.time()
    g = fine._psf_grid(edges)
    new = fine.point_sources_cells(g, a.cell, edges, exact=True)
    z = dict(np.load(p))
    old = z["T__point_sources"]
    print(f"{a.cell:g} deg: totals new/old per bin", np.round(new.reshape(9, -1).sum(1) / old.reshape(9, -1).sum(1), 5))
    keep = z["keep"].astype(bool)
    for k in range(9):
        d = np.abs(new[k] - old[k])[keep]
        print(f"  bin {k}: moved fraction {0.5 * d.sum() / old[k][keep].sum():.4f}, max cell change {d.max():,.0f} counts")
    z["T__point_sources"] = new
    np.savez(p, **z)
    print(f"written {p} ({time.time()-t0:.0f}s); previous library kept as {bak.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
