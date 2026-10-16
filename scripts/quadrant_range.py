#!/usr/bin/env python
"""How much of the quadrant-gate result is the companion set? (Kev, 2026-09-21)

Every "x of N rows favor a halo" number we quote is conditional on ONE
companion set. Measured 2026-09-21, holding the 128 models fixed and changing
only the companions, Oscar's structured panel moves between 15/257 and
133/257 -- 6% to 52% -- while the halo amplitude and its significance barely
move. So the gate responds to how much large-scale freedom the rest of the
model carries, and a single number presents a choice as a measurement.

This is the cheap, non-hand-picked version of the fix: score each row under
the companion set THAT ROW'S OWN greedy path selected (Stage C already
records them in `own_paths`; there are 13 distinct sets across the 257 rows).
One fit per row, so it costs one Stage B pass, not one per set.

    python scripts/quadrant_range.py                       # the residual panel
    python scripts/quadrant_range.py --stage-c <json>

Writes `results/quadrant_range_<date>.{json,md}`. Compare its headline with
the consensus-set number in the matching Stage C run: the difference is the
width the figure's subtitle is currently hiding.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rank_stage_c import _score, EDGES  # noqa: E402  (same scorer, same conventions)

DATE = "2026-09-21"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-c", default=str(ROOT / "results" / "fits" / f"rank_stage_c_{DATE}.json"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    t0 = time.time()

    src = json.loads(Path(args.stage_c).read_text())
    consensus = tuple(src["consensus_selected"])
    paths = {(p["model"], p["variant"]): tuple(p["chosen"])
             for p in src["own_paths"] if "chosen" in p}
    jobs = [(m, v, s) for (m, v), s in paths.items()]
    distinct = collections.Counter(s for s in paths.values())
    print(f"{len(jobs)} rows, {len(distinct)} distinct self-selected companion sets; "
          f"consensus was {list(consensus)}", flush=True)

    with Pool(args.workers) as pool:
        rows = pool.map(_score, jobs)
    for r, (m, v, s) in zip(rows, jobs):
        r["chosen"] = list(s)
        r["neg2lnl_total"] = float(sum(r["neg2lnl_bg_only"]))
    print(f"scored in {time.time()-t0:.0f}s", flush=True)

    # the consensus-set scoring of the same rows, for the side-by-side
    by_row = {(r["model"], r["variant"]): r for r in src["scored"]}

    def gate(r):
        return (r["ts"][5] > 4 and r["halo_counts"][5] > 0
                and not r["quadrant_any_negative"][5])

    phys = [r for r in rows if r["model"] != "gll_iem_v07"]
    own_fav = sum(1 for r in phys if gate(r))
    con = [by_row[(r["model"], r["variant"])] for r in phys if (r["model"], r["variant"]) in by_row]
    con_fav = sum(1 for r in con if gate(r))

    # per-row agreement: does the verdict depend on which set that row is given?
    flip = [(r["model"], r["variant"], gate(by_row[(r["model"], r["variant"])]), gate(r))
            for r in phys if (r["model"], r["variant"]) in by_row]
    both = sum(1 for _, _, a, b in flip if a and b)
    neither = sum(1 for _, _, a, b in flip if not a and not b)
    only_con = sum(1 for _, _, a, b in flip if a and not b)
    only_own = sum(1 for _, _, a, b in flip if not a and b)

    h_own = [r["halo_counts"][5] for r in phys]
    h_con = [r["halo_counts"][5] for r in con]

    out = {"date": DATE, "stage_c": str(args.stage_c), "consensus": list(consensus),
           "n_rows": len(phys), "n_distinct_sets": len(distinct),
           "favor_consensus_set": con_fav, "favor_own_set": own_fav,
           "agree_both": both, "agree_neither": neither,
           "only_consensus": only_con, "only_own": only_own,
           "median_halo_consensus": float(np.median(h_con)),
           "median_halo_own": float(np.median(h_own)),
           "distinct_sets": [{"set": list(s), "n_rows": n} for s, n in distinct.most_common()],
           "rows": rows}
    tag = f"_{args.tag}" if args.tag else ""
    jp = ROOT / "results" / "fits" / f"quadrant_range_{DATE}{tag}.json"
    jp.write_text(json.dumps(out, indent=2))

    pct = lambda k, n: f"{k}/{n} ({100*k/n:.0f}%)"
    lines = [
        f"# How much of the quadrant gate is the companion set? — {DATE}", "",
        f"Source: `{Path(args.stage_c).name}`. {len(phys)} physical rows, "
        f"{len(distinct)} distinct self-selected companion sets.", "",
        "| scoring | >2σ AND quadrant-consistent | median halo @ 20.8 |",
        "|---|---|---|",
        f"| one consensus set for every row: `{', '.join(consensus)}` | {pct(con_fav, len(con))} | {np.median(h_con):,.0f} |",
        f"| each row under its OWN selected set | {pct(own_fav, len(phys))} | {np.median(h_own):,.0f} |",
        "",
        "Per-row agreement between the two scorings:", "",
        f"- passes under both: **{both}**", f"- passes under neither: **{neither}**",
        f"- passes only with the consensus set: **{only_con}**",
        f"- passes only with its own set: **{only_own}**", "",
        f"**{only_con + only_own} of {len(flip)} rows ({100*(only_con+only_own)/len(flip):.0f}%) change verdict "
        "on the companion set alone.** That is the width a single subtitle number hides.", "",
        "## The self-selected sets", "",
        "| rows | set |", "|---|---|",
    ]
    for s, n in distinct.most_common():
        lines.append(f"| {n} | `{', '.join(s)}` |")
    mp = ROOT / "results" / f"quadrant_range_{DATE}{tag}.md"
    mp.write_text("\n".join(lines) + "\n")

    print(f"\n  consensus set for every row : {pct(con_fav, len(con))} quadrant-consistent at >2σ")
    print(f"  each row its own set        : {pct(own_fav, len(phys))}")
    print(f"  change verdict on companions: {only_con + only_own}/{len(flip)} "
          f"({100*(only_con+only_own)/len(flip):.0f}%)")
    print(f"\nsaved {jp}, {mp}; total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
