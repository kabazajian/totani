#!/usr/bin/env python
"""A2020 redo, STAGE 2a: reduce the data at A2020's own selection.

    .venv/bin/python scripts/a2020_stage2_reduce.py [--plan]

Stage 1 validated the limit chain in isolation and PASSED (0.15%,
`validation/a2020_stage1_result_2026-08-08.md`). Stage 2 validates everything
upstream of it -- reduction, templates, likelihood -- against Table II's TS
values. This script is the first half: the reduction.

EVERY FIELD IS FROM THE PAPER'S APPENDIX A OR ITS TABLE I, NOT INFERRED
-----------------------------------------------------------------------
    eight years, 2008-08-04 to 2016-08-02   weeks 9-426 (see below)
    P8R3 ULTRACLEANVETO, all event types    evclass=1024, evtype=3
    0.667-158 GeV, 15 bins                  their Ebands.dat, exactly
    40 x 40 deg about the GC                roi_half_deg=20
    zenith < 90 deg                         zmax_deg=90
    (DATA_QUAL>0)&&(LAT_CONFIG==1)          their filter, NOT our Cholis default
    IRF P8R3_ULTRACLEANVETO_V2              irf_override -- our default is V3

⚠ THREE FIDELITY DIFFERENCES, RECORDED BEFORE RUNNING RATHER THAN DISCOVERED
AFTER:

1. **Week granularity.** A2020 cut on dates; we select whole weekly files. Week
   426 spans MET 491358610-491962466 and their stop (2016-08-02, MET 491788800)
   falls inside it, so including week 426 takes ~2 days more data than they did
   -- 0.07% of eight years. `week_min=9` needs no such caveat: the mission itself
   begins 2008-08-04T15:43:36, so week 9 IS their start.
2. ~~Pixel size is UNCONFIRMED.~~ ✅ **RESOLVED 2026-08-08: 0.1 deg, confirmed by
   Kev, an author of A2020.** The paper's Appendix A does not state it, so it
   could not come from the PDF; rule 1a's "a human who checked it" is what makes
   it anchorable. It matches what was assumed, but was only recorded once
   confirmed. So this is now a MATCH, not a difference -- our grid and theirs are
   the same, and the Pohl templates interpolate onto it.
3. ~~4FGL-DR4, not DR1.~~ ✅ **RESOLVED 2026-08-08: we use DR1, their catalog**
   (Kev, correcting an earlier choice of mine). This was not the small difference
   my "note it" framing implied: DR1 has **473** sources in the 40x40 ROI against
   DR4's **715**, so DR4 would have added **256 point sources A2020 never had**.
   The DR1 version is not recoverable from the paper, which cites only "4FGL" --
   but it does not matter, because v21 and v22 have IDENTICAL in-ROI source lists.
   v21 is used, v22 postdating submission.

None of these is expected to move TS materially, which is why TS is the anchor
and absolute -log L is informational (rule: `a2020_redo_plan_2026-08-08.md` §3a).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

#: A2020's 15 analysis bins, GeV. Read from their own `Ebands.dat` (keV), which
#: ships with the published pipeline -- a machine-readable source, anchorable
#: under rule 1a. The top bin 37.49-158.1 GeV is their single "macroenergy" bin.
A2020_EDGES_GEV = (
    0.6667607, 0.8891397, 1.18569, 1.581139, 2.108482, 2.811707, 3.749471,
    5.0, 6.667607, 8.891397, 11.85687, 15.81139, 21.08483, 28.11707, 37.49471,
    158.1139,
)


def a2020_config():
    from gcereduce.reduce import ReductionConfig
    return ReductionConfig(
        week_min=9,
        week_max=426,
        energy_edges_gev=A2020_EDGES_GEV,
        roi_half_deg=20.0,
        pixel_deg=0.1,
        evclass=1024,                 # P8R3 ULTRACLEANVETO
        evtype=3,                     # FRONT+BACK, i.e. "all event types"
        zmax_deg=90.0,
        gtmktime_filter="(DATA_QUAL>0)&&(LAT_CONFIG==1)",
        irf_override="P8R3_ULTRACLEANVETO_V2",
        note="A2020 (2003.10416) Appendix A selection, for the Stage 2 validation",
    )


def main() -> int:
    from gcereduce.reduce import reduce as run_reduce

    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", action="store_true", help="print the config, run nothing")
    a = ap.parse_args()

    cfg = a2020_config()
    print("A2020 STAGE 2a -- reduction at their selection")
    print(f"  config hash : {cfg.hash}")
    print(f"  grid        : {cfg.n_ebins} x {cfg.npix} x {cfg.npix}")
    print(f"  ltcube key  : {cfg.ltcube_key}")
    print(f"  weeks       : {cfg.week_min}-{cfg.week_max}")
    print(f"  IRF         : {cfg.irf_override}  (our default for this evclass is V3)")
    print(f"  filter      : {cfg.gtmktime_filter}")
    if a.plan:
        print("\n  --plan: nothing run")
        return 0

    t0 = time.time()
    wd = run_reduce(cfg)
    dt = time.time() - t0
    print(f"\n  reduction  : {wd}")
    print(f"  wall clock : {dt/60:.1f} min")
    out = ROOT / "runs" / "a2020_stage2_reduction.json"
    out.write_text(json.dumps(
        {"hash": cfg.hash, "workdir": str(wd), "reduction_seconds": dt,
         "config": cfg.to_dict()}, indent=1, default=str))
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
