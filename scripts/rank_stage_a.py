#!/usr/bin/env python
"""Stage A: build every cell-template cache the ranking needs.

    python scripts/rank_stage_a.py                 # everything available now
    python scripts/rank_stage_a.py --models-only   # just the Paper II models
    python scripts/rank_stage_a.py --force         # rebuild

Models are taken from `external/galprop_paperII/` (extracted) or from the
tarballs gce-fisher is fetching to NFS; models with neither are listed and
skipped, so this can be re-run as the download proceeds. Structures whose
inputs are not yet on disk (radio-loop geometry, gll_iem_v07, Oscar's
eROSITA brightness map) are skipped with a reason.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking  # noqa: E402
from totani.structures import STRUCTURES  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models-only", action="store_true")
    ap.add_argument("--structures-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--variants", default="gas1_ics1,gas5_ics3")
    args = ap.parse_args()
    t0 = time.time()

    if not args.structures_only:
        print("shared cells ...", flush=True)
        ranking.shared_cells(EDGES, force=args.force)
        print(f"  done ({time.time()-t0:.0f}s)", flush=True)

    if not args.models_only:
        for name in STRUCTURES:
            try:
                ranking.structure_cells(name, EDGES, force=args.force)
                print(f"structure {name:18s} ok ({time.time()-t0:.0f}s)", flush=True)
            except (KeyError, FileNotFoundError, ValueError) as e:
                print(f"structure {name:18s} SKIPPED: {e}", flush=True)

    if not args.structures_only:
        variants = args.variants.split(",")
        missing, built = [], 0
        for code in ranking.MODEL_CODES:
            if not ranking.paperii_available(code):
                missing.append(code); continue
            for v in variants:
                ranking.model_cells(code, v, EDGES, force=args.force)
            built += 1
            print(f"model {code:22s} {'/'.join(variants)} ok ({time.time()-t0:.0f}s)", flush=True)
        print(f"\n{built} models cached, {len(missing)} not yet available: {missing[:6]}{' ...' if len(missing) > 6 else ''}")
    print(f"total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
