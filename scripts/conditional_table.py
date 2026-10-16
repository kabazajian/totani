#!/usr/bin/env python
"""Oscar's item 1: conditional likelihood comparisons for the principal
complete models.

    python scripts/conditional_table.py

For each representative row and each complete model (the split's 09-23
selections, with the residual construction allowed and excluded):
  forward   H0 = diffuse + structures;  H1 = H0 + halo (signed, per bin):
            TS per bin and summed, the 20.8 GeV amplitude with its 68%
            profile interval, the joint b-bbar 500 GeV <sigma v> with interval;
  reverse   with the halo ALLOWED, swap the bubble template among {flat,
            structured, residual construction} keeping the other structures:
            -2lnL of each (H1 minimum summed over bins) relative to the
            selected one; and add/remove the eROSITA pair.
Every comparison refits all shared normalizations. Cell space, 10 deg.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import dmspectra, ranking, selection  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-24"
FITS = ROOT / "results" / "fits"
K20 = 5
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"),
        ("SL_Z10_R30_T150_C5", "gas1_ics1"), ("SL_Z10_R30_T150_C5", "gas5_ics3"),
        ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SS_Z10_R30_T150_C5", "gas5_ics3"),
        ("a2020", "-"), ("G21_SA0_R12+iso", "gas1_ics1"), ("G21_SA0_R12+iso", "gas5_ics3"),
        ("G21_SA100_F98+iso", "gas1_ics1"), ("G21_SA100_F98+iso", "gas5_ics3")]
BUBBLES = ("bubble_flat", "bubble_structured", "bubble_totani")
EROSITA = tuple(selection.GROUPS["erosita"])


def sel(variant, residual):
    p = FITS / f"rank_stage_c_2026-09-23_{variant}{'' if residual else '_no_residual'}.json"
    return tuple(json.loads(p.read_text())["consensus_selected"])


def stack(row, cands, sh):
    st = selection.base_stack(row[0], row[1], EDGES, sh)
    for c in cands:
        st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
    return st


def h0_h1(st, J):
    keep = halo_cell_mask(); sr, Jr = restrict_stack(st, keep, J)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=True)["per_bin"]
    roi = float(Jr[K20][keep].sum())
    ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(9)]
    return {"neg2lnl_H0": float(sum(b["neg2lnl"] for b in bg)), "neg2lnl_H1": float(sum(p["neg2lnl"] for p in pf)),
            "ts_bins": ts, "ts_sum": float(sum(ts)), "halo20": pf[K20]["signal_amplitude"] * roi,
            "halo20_lo": pf[K20]["signal_amplitude_lo"] * roi, "halo20_hi": pf[K20]["signal_amplitude_hi"] * roi,
            "halo_bins": [pf[k]["signal_amplitude"] * float(Jr[k][keep].sum()) for k in range(9)]}


def main() -> int:
    sh = ranking.shared_cells(EDGES); J = sh["J"]
    out = {"date": DATE, "rows": []}
    lines = [f"# Conditional likelihood comparisons ({DATE})", "",
             "H0 = diffuse + selected structures; H1 = H0 + halo (signed amplitude per bin). TS summed over 9 bins and at 20.8 GeV; "
             "halo at 20.8 GeV with 68% profile interval; joint b b̄ 500 GeV ⟨σv⟩ (10⁻²⁵) with 68% interval and TS. "
             "Reverse: with the halo allowed, −2lnL(H1) of the alternative bubble templates and of the eROSITA pair relative to the selected model "
             "(negative = alternative preferred; +p = parameters added per bin).", "",
             "| row | model (split, residual) | selected | TS sum / TS@20.8 | halo@20.8 [68%] | ⟨σv⟩ [68%] (TS) | reverse: flat / structured / residual | reverse: eROSITA pair |",
             "|---|---|---|---|---|---|---|---|"]
    for row in ROWS:
        v = row[1] if row[1] != "-" else "gas1_ics1"
        for residual in (True, False):
            cands = sel(v, residual)
            if row[0] == "gll_iem_v07":
                continue
            try:
                st = stack(row, cands, sh)
            except (FileNotFoundError, KeyError):
                continue
            fwd = h0_h1(st, J)
            sv = dmspectra.joint_fit(st, J, "b", 500, EDGES, n_scan=31)
            # reverse: bubble alternatives with the halo allowed
            bub = next((c for c in cands if c in BUBBLES), None)
            others = tuple(c for c in cands if c not in BUBBLES)
            rev = {}
            for b in BUBBLES:
                try:
                    h = h0_h1(stack(row, others + (b,), sh), J)
                    rev[b] = h["neg2lnl_H1"] - fwd["neg2lnl_H1"]
                except (FileNotFoundError, KeyError):
                    rev[b] = None
            has_ero = [c for c in cands if c in EROSITA]
            if has_ero:
                h = h0_h1(stack(row, tuple(c for c in cands if c not in EROSITA), sh), J)
                ero = {"removed": h["neg2lnl_H1"] - fwd["neg2lnl_H1"], "which": has_ero[0]}
            else:
                best = None
                for e in ("erosita_unif_NS", "erosita_bright_NS"):
                    h = h0_h1(stack(row, cands + (e,), sh), J)
                    d = h["neg2lnl_H1"] - fwd["neg2lnl_H1"]
                    if best is None or d < best[1]:
                        best = (e, d)
                ero = {"added": best[1], "which": best[0]}
            rec = {"model": row[0], "variant": row[1], "split": v, "residual_allowed": residual, "selected": list(cands),
                   **fwd, "sigmav": {k: sv[k] for k in ("sigmav_best", "sigmav_lo68", "sigmav_hi68", "ts_vs_zero")},
                   "reverse_bubbles": rev, "reverse_erosita": ero}
            out["rows"].append(rec)
            rv = " / ".join("—" if rev[b] is None else f"{rev[b]:+.0f}" for b in BUBBLES)
            er = (f"remove {ero['which']}: {ero['removed']:+.0f}" if "removed" in ero else f"add {ero['which']}: {ero['added']:+.0f}")
            lines.append(f"| {row[0]} {row[1]} | {v}, {'with' if residual else 'without'} | {' → '.join(c.replace('bubble_', '').replace('_', ' ') for c in cands)} | "
                         f"{fwd['ts_sum']:.0f} / {fwd['ts_bins'][K20]:.1f} | {fwd['halo20']:,.0f} [{fwd['halo20_lo']:,.0f}, {fwd['halo20_hi']:,.0f}] | "
                         f"{sv['sigmav_best']*1e25:+.1f} [{sv['sigmav_lo68']*1e25:+.1f}, {sv['sigmav_hi68']*1e25:+.1f}] ({sv['ts_vs_zero']:.0f}) | {rv} | {er} |")
            print(lines[-1], flush=True)
    (FITS / f"conditional_table_{DATE}.json").write_text(json.dumps(out, indent=1))
    (ROOT / "results" / f"conditional_table_{DATE}.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
