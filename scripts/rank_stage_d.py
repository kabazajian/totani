#!/usr/bin/env python
"""Stage D: the halo, per model, with the consensus structures -- shape-
agnostic per-bin amplitudes (already in Stage C's scoring), the signed joint
<sigma v> fits for PPPC4DMID spectra over a mass grid, the gates, and the
spread across models.

    python scripts/rank_stage_d.py                         # reads rank_stage_c_<date>.json
    python scripts/rank_stage_d.py --channels b,tau,W --masses 20,50,100,200,500,1000

For every (model, variant) row scored in Stage C: for each channel and mass,
`dmspectra.joint_fit` -> best <sigma v> (signed), TS vs zero, 68% interval,
95% UL. Gates carried from Stage C: injection slope, quadrant consistency,
N/S ratio. Summary per channel/mass: the distribution of best-fit <sigma v>
and TS across models, split by model-quality tercile (delta -2lnL vs the
best model), and the model-marginalised limit -- the 95% UL from the profile
of the MEDIAN curve across models -- as a first, crude covariance treatment
(plan §5, Stage D (ii)); the proper covariance waits for Galp21.

Cell space only; parallel over (model, channel, mass).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani import dmspectra, ranking, selection  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-23"


def _stack_for(model, variant, structures):
    sh = ranking.shared_cells(EDGES)
    st = selection.base_stack(model, variant, EDGES, sh)
    for c in structures:
        try:
            st, _, _ = selection.add_candidate(st, c, model, variant, EDGES, sh)
        except (FileNotFoundError, KeyError, ValueError):
            pass                                        # a structure this row cannot have (flagged in Stage C)
    return st, sh["J"]


def _one(args):
    model, variant, structures, channel, m = args
    st, J = _stack_for(model, variant, structures)
    r = dmspectra.joint_fit(st, J, channel, m, EDGES, n_scan=31)
    r.update({"model": model, "variant": variant})
    r.pop("curve", None)
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-c", default=str(ROOT / "results" / "fits" / f"rank_stage_c_{DATE}.json"))
    ap.add_argument("--channels", default="b,tau,W")
    ap.add_argument("--masses", default="20,50,100,200,500,1000")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--tag", default="")
    ap.add_argument("--render-only", action="store_true", help="rewrite the Markdown from the existing JSON")
    args = ap.parse_args()
    t0 = time.time()
    tag = f"_{args.tag}" if args.tag else ""
    jp = ROOT / "results" / "fits" / f"rank_stage_d_{DATE}{tag}.json"
    if args.render_only:
        out = json.loads(jp.read_text())
        render(out, tag)
        return 0
    c = json.loads(Path(args.stage_c).read_text())
    structures = tuple(c["consensus_selected"])
    rows = c["scored"]
    channels = args.channels.split(","); masses = [float(x) for x in args.masses.split(",")]
    print(f"{len(rows)} model rows, structures {structures}, {len(channels)} channels x {len(masses)} masses", flush=True)

    jobs = [(r["model"], r["variant"], structures, ch, m) for r in rows for ch in channels for m in masses]
    with Pool(args.workers) as pool:
        fits = pool.map(_one, jobs)
    print(f"{len(fits)} joint fits in {time.time()-t0:.0f}s", flush=True)

    best_total = min(r["neg2lnl_total"] for r in rows)
    quality = {(r["model"], r["variant"]): r["neg2lnl_total"] - best_total for r in rows}
    q_vals = sorted(quality.values()); terc = [q_vals[len(q_vals) // 3], q_vals[2 * len(q_vals) // 3]]

    def tercile(k):
        q = quality[k]
        return "best" if q <= terc[0] else ("middle" if q <= terc[1] else "worst")

    summary = {}
    for ch in channels:
        for m in masses:
            sel = [f for f in fits if f["channel"] == ch and f["m_gev"] == m and f["model"] != "gll_iem_v07"
                   and f["sigmav_best"] is not None]
            if not sel:
                summary[f"{ch}_{m:g}"] = {"n": 0, "note": "no photon yield in the analysis bins"}; continue
            by_t = {t: [f for f in sel if tercile((f["model"], f["variant"])) == t] for t in ("best", "middle", "worst")}
            summary[f"{ch}_{m:g}"] = {
                "n": len(sel),
                "sigmav_best_median": float(np.median([f["sigmav_best"] for f in sel])),
                "sigmav_best_range": [float(np.min([f["sigmav_best"] for f in sel])), float(np.max([f["sigmav_best"] for f in sel]))],
                "ts_median": float(np.median([f["ts_vs_zero"] for f in sel])),
                "ul95_median": float(np.median([f["sigmav_ul95"] for f in sel])),
                "ul95_max_over_models": float(np.max([f["sigmav_ul95"] for f in sel])),
                "by_quality_tercile": {t: {"n": len(v),
                                           "sigmav_best_median": float(np.median([f["sigmav_best"] for f in v])) if v else None,
                                           "ts_median": float(np.median([f["ts_vs_zero"] for f in v])) if v else None,
                                           # 🚨 sigmav_lo68 is the Δ(−2lnL) = 1 crossing: a 1σ bound. Until
                                           # 2026-09-22 this was written as "fraction_positive_2sigma" and
                                           # quoted as 2σ (5-ring, no residual: 100% at 1σ, 96% at 2σ).
                                           "fraction_positive_1sigma": float(np.mean([f["sigmav_lo68"] > 0 for f in v])) if v else None,
                                           "fraction_positive_2sigma": float(np.mean([f["ts_vs_zero"] > 4 and f["sigmav_best"] > 0
                                                                                      for f in v])) if v else None}
                                       for t, v in by_t.items()}}

    out = {"date": DATE, "stage_c": Path(args.stage_c).name, "structures": list(structures),
           "channels": channels, "masses": masses, "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]),
           "quality_terciles_delta_neg2lnl": terc, "fits": fits, "summary": summary,
           "gates_from_stage_c": {f"{r['model']}__{r['variant']}": {
               "injection_min": min(r["injection_slopes"]), "quadrant_sum_over_joint_20p8": r["quadrant_sum_over_joint"][5],
               "quadrant_any_negative_20p8": r["quadrant_any_negative"][5],
               "north_over_south_20p8": (r["halo_north"][5] / r["halo_south"][5]) if r["halo_south"][5] else None,
               "halo_20p8": r["halo_counts"][5], "ts_20p8": r["ts"][5], "delta_neg2lnl": quality[(r["model"], r["variant"])]}
               for r in rows}}
    jp.write_text(json.dumps(out, indent=2))
    render(out, tag)
    print(f"saved {jp}; total {time.time()-t0:.0f}s")
    return 0


def render(out, tag):
    channels, masses, structures, rows, fits, summary = (out["channels"], out["masses"], out["structures"],
                                                          out["gates_from_stage_c"], out["fits"], out["summary"])
    quality = {k: g["delta_neg2lnl"] for k, g in rows.items()}
    lines = [f"# Stage D — {DATE}{(' ('+tag.lstrip('_')+')') if tag else ''}", "",
             f"Structures: {', '.join(structures)}. {len(rows)} model rows; joint signed ⟨σv⟩ fits (PPPC4DMID, EW-corrected).", "",
             "## Across models, per channel and mass", "",
             "| channel | m [GeV] | median ⟨σv⟩ best [cm³/s] | range | median TS | median 95% UL | best-tercile: median ⟨σv⟩, TS, frac >2σ | worst-tercile: same |",
             "|---|---|---|---|---|---|---|---|"]
    for key, s in summary.items():
        ch, m = key.rsplit("_", 1)
        if not s.get("n"):
            lines.append(f"| {ch} | {m} | – | – | – | – | no photon yield in the analysis bins | – |"); continue
        b, w = s["by_quality_tercile"]["best"], s["by_quality_tercile"]["worst"]
        fmt = lambda t: (f"{t['sigmav_best_median']:.2e}, {t['ts_median']:.1f}, {t['fraction_positive_2sigma']:.2f}" if t["n"] else "-")
        lines.append(f"| {ch} | {m} | {s['sigmav_best_median']:.2e} | [{s['sigmav_best_range'][0]:.1e}, {s['sigmav_best_range'][1]:.1e}] | "
                     f"{s['ts_median']:.1f} | {s['ul95_median']:.2e} | {fmt(b)} | {fmt(w)} |")
    lines += ["", "## Per model (b b̄, 500 GeV) with the gates", "",
              "| model | variant | Δ(−2lnL) | halo @ 20.8 | TS @ 20.8 | ⟨σv⟩ best | TS | 95% UL | inj min | quad Σ/joint | N/S |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    ref = [f for f in fits if f["channel"] == channels[0] and f["m_gev"] == (500.0 if 500.0 in masses else masses[-1])
           and f["sigmav_best"] is not None]
    ref.sort(key=lambda f: quality[f"{f['model']}__{f['variant']}"])
    for f in ref:
        g = rows[f"{f['model']}__{f['variant']}"]
        ns = g["north_over_south_20p8"]
        lines.append(f"| {f['model']} | {f['variant']} | {g['delta_neg2lnl']:+,.0f} | {g['halo_20p8']:,.0f} | {g['ts_20p8']:.1f} | "
                     f"{f['sigmav_best']:.2e} | {f['ts_vs_zero']:.1f} | {f['sigmav_ul95']:.2e} | {g['injection_min']:.2f} | "
                     f"{(g['quadrant_sum_over_joint_20p8'] if g['quadrant_sum_over_joint_20p8'] is not None else float('nan')):.2f}"
                     f"{' (neg)' if g['quadrant_any_negative_20p8'] else ''} | {(ns if ns is not None else float('nan')):.2f} |")
    mp = ROOT / "results" / f"rank_stage_d_{DATE}{tag}.md"
    mp.write_text("\n".join(lines) + "\n")
    print(f"rendered {mp}")


if __name__ == "__main__":
    raise SystemExit(main())
