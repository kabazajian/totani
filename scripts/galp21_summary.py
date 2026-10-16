#!/usr/bin/env python
"""The six Galp21 models in one table, both IC treatments, from the first-look
JSONs: fit rank, and the 20.8 GeV halo / gate / <sigma v> under the residual
selection, the no-residual selection, and each row's own no-residual path.

    python scripts/galp21_summary.py results/fits/galp21_first_look_*_2026-09-23.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-09-23"


def main(paths) -> int:
    rows = []
    for p in paths:
        rows += json.loads(Path(p).read_text())["rows"]
    rows.sort(key=lambda r: (r["model"].endswith("+iso"), r["variant"], r["model"]))
    lines = [f"# Galp21 in the Totani ROI — six models, two IC treatments ({DATE})", "",
             "Ranks are of 130 (128 2012 models + A2020 + the row). Halo at 20.8 GeV in counts; 'gate' = >2σ and "
             "quadrant-consistent; ⟨σv⟩ b b̄ 500 GeV in 10⁻²⁵ cm³ s⁻¹ with its TS.", "",
             "| model | IC | split | bare rank | rank, residual sel. | halo / gate / ⟨σv⟩ (TS), residual sel. | "
             "rank, no-residual sel. | halo / gate / ⟨σv⟩, no-residual sel. | own no-residual path: bubble → halo / gate / ⟨σv⟩ |",
             "|---|---|---|---|---|---|---|---|---|"]
    agg = {}
    for r in rows:
        m = r["model"][4:]; ic = "iso" if m.endswith("+iso") else "aniso"; m = m.replace("+iso", "")
        cr, cn, own = r["sets"]["C_residual"], r["sets"]["C_no_residual"], r["own_no_residual"]
        bub = next((c for c in own["chosen"] if c.startswith("bubble_")), "none").replace("bubble_", "")
        f = lambda s: f"{s['halo_20']:,.0f} / {'Y' if s['favor'] else 'n'} / {s['sigmav']['sigmav_best']*1e25:+.1f} ({s['sigmav']['ts_vs_zero']:.0f})"
        lines.append(f"| {m} | {ic} | {r['variant']} | {r['bare']['rank']} | {cr['rank']} | {f(cr)} | {cn['rank']} | {f(cn)} | "
                     f"{bub} → {own['halo_20']:,.0f} / {'Y' if own['favor'] else 'n'} / {own['sigmav_best']*1e25:+.1f} ({own['sigmav_ts']:.0f}) |")
        agg.setdefault((ic, r["variant"]), []).append((bub, own["favor"], own["sigmav_best"], cr["favor"], cr["sigmav"]["sigmav_best"], r["bare"]["rank"]))
    lines += ["", "## Summary by IC treatment and split (6 models each)", "",
              "| IC | split | bare rank range | residual sel.: gate passed / ⟨σv⟩ median | own no-residual: bubble chosen | gate passed | ⟨σv⟩ median |",
              "|---|---|---|---|---|---|---|"]
    for (ic, v), a in sorted(agg.items()):
        bubs = {b: sum(x[0] == b for x in a) for b in sorted({x[0] for x in a})}
        lines.append(f"| {ic} | {v} | {min(x[5] for x in a)}–{max(x[5] for x in a)} | {sum(x[3] for x in a)}/6 / {np.median([x[4] for x in a])*1e25:+.1f} | "
                     + ", ".join(f"{b} {n}" for b, n in bubs.items()) + f" | {sum(x[1] for x in a)}/6 | {np.median([x[2] for x in a])*1e25:+.1f} |")
    out = ROOT / "results" / f"galp21_summary_{DATE}.md"
    out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
