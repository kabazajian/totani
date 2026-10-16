#!/usr/bin/env python
"""Stream the Galp21 |b| <= 60 archives onto the Totani target grid.

    python scripts/galp21_ingest.py SA0_R12 SA100_F98        # ~1 h per model, one at a time
    python scripts/galp21_ingest.py SA0_R12 --dry-run 4      # parse the first 4 members only

Writes `results/fits/ranking_cells/galp21_cubes__<model>.npz` (+ `.json`
manifest). One process, one model at a time: ~5 GB peak (eleven (44, 960,
960) float64 spectrum accumulators), xz on 4 threads. See `totani.galp21`.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import galp21, ranking  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("models", nargs="+")
    ap.add_argument("--dry-run", type=int, default=None, metavar="N")
    ap.add_argument("--xz-threads", type=int, default=4)
    args = ap.parse_args()
    _, L, B, _ = ranking._grid_and_basis(EDGES)
    for m in args.models:
        out = ranking.CACHE / f"galp21_cubes__{m}.npz"
        if out.exists() and args.dry_run is None:
            print(f"{m}: {out.name} exists"); continue
        r = galp21.ingest_model(m, L, B, EDGES, None if args.dry_run else out,
                                xz_threads=args.xz_threads, limit=args.dry_run,
                                log=lambda s: print(s, flush=True))
        print(json.dumps({k: v for k, v in r.items() if k != "manifest"}, indent=1)[:3000], flush=True)
        print("peak RSS %.2f GB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
