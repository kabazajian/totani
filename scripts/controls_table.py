#!/usr/bin/env python
"""tab:controls (App. C): the original quadrant diagnostic for the structured
Fermi bubbles template with Loop I shells and each set of accompanying foregrounds.

    python scripts/controls_table.py

Reads the Stage B outputs (results/fits/rank_stage_b_<date>_structured_shells_<set>.json),
all made with the fitter as corrected on 2026-09-22, and writes
results/controls_table_2026-10-15.md. Over the 257 configurations (gll_iem_v07 excluded):

  pass     TS_20.8 > 4, positive 20.8 GeV NFW counts, and no negative quadrant amplitude
  >2sigma  TS_20.8 > 4
  counts   median 20.8 GeV NFW counts
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FITS = ROOT / "results" / "fits"
SETS = [
    ("solar IC", "rank_stage_b_2026-09-24_structured_shells_SOLARIC_ctrl_rerun2026-10-15.json"),
    ("geometric eROSITA outer N/S", "rank_stage_b_2026-09-24_structured_shells_GEOMNS_ctrl_rerun2026-10-15.json"),
    ("uniform eROSITA outer N/S", "rank_stage_b_2026-09-24_structured_shells_UNIFOUTER_ctrl_rerun2026-10-15.json"),
    ("uniform eROSITA whole N/S", "rank_stage_b_2026-09-24_structured_shells_sun.json"),
]


def summarize(path):
    d = json.loads(path.read_text())
    rows = d.get("rows") or d.get("results") or d
    if isinstance(rows, dict):
        rows = list(rows.values())
    rows = [r for r in rows if r.get("model") != "gll_iem_v07"]
    ts = np.array([r["ts"][5] for r in rows])
    bad = np.array([bool(r["quadrant_any_negative"][5]) or r["halo_counts"][5] <= 0 for r in rows])
    counts = np.array([r["halo_counts_20p8"] for r in rows])
    return len(rows), int(((ts > 4) & ~bad).sum()), int((ts > 4).sum()), float(np.median(counts))


def main() -> int:
    L = ["# tab:controls: accompanying-foreground tests, structured bubbles + Loop I shells, 10-deg cells (2026-10-15)", "",
         "| companions | pass | >2sigma | median counts at 20.8 GeV | source |", "|---|---|---|---|---|"]
    for name, f in SETS:
        n, npass, n2, med = summarize(FITS / f)
        L.append(f"| {name} | {npass}/{n} | {n2} | {med:.0f} | `{f}` |")
        print(L[-1], flush=True)
    (ROOT / "results" / "controls_table_2026-10-15.md").write_text("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
