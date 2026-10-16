#!/usr/bin/env python
"""Galp21 in the Totani ROI: how the new models fit, and what the 20 GeV
excess does with them -- every number set beside the 2012-grid rows it is
compared with.

    python scripts/galp21_first_look.py                  # the ingested Galp21 models
    python scripts/galp21_first_look.py --models SA0_R12

For each Galp21 model x component split (Totani's 1-gas and 5-ring):

  A. FIT QUALITY. Background-only -2lnL (|b| >= 10, 9 bins) of the bare
     diffuse model (diffuse + isotropic + point sources), and under the three
     Figure-1 structure sets and the split's two Stage C selections; each
     placed in the ranking of the 2012 rows (+ A2020) under the same set. The
     fitted gas / IC normalizations (1 = Galp21's absolute prediction) and
     the mean |residual| per halo cell.
  B. THE EXCESS. Under each set: the 20.8 GeV halo amplitude with 68%
     interval, TS, the quadrant gate ('favor'), N/S, and the joint b-bbar
     500 GeV <sigma v> (TS > 4 = 2 sigma; NOT lo68 > 0).
  C. ITS OWN SELECTION: each Galp21 row's greedy path, with and without
     Totani's residual construction (rebuilt from the Galp21 model).

Single process; the residual-template builds read the full cubes (~3 GB).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import dmspectra, galp21, ranking, selection  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-23"
FITS = ROOT / "results" / "fits"
K20 = 5

#: set name -> (source JSON, row key, how the set is specified)
SETS = {
    "fig1_residual": (FITS / "rank_stage_c_2026-09-21.json", "scored", "candidates"),
    "fig1_flat": (FITS / "rank_stage_b_2026-09-21_flat_shells_sun.json", "rows", "structures"),
    "fig1_structured": (FITS / "rank_stage_b_2026-09-21_structured_shells_sun.json", "rows", "structures"),
    # 09-23 = the corrected-fitter selections (Loop I shells + NPS in the residual
    # branch). ⚠ Every first-look run up to and including the six-model set of
    # 2026-09-23 used the 09-21 files here (shells without NPS); the six are
    # mutually comparable, and the halo numbers barely respond to that template.
    "C_residual": (FITS / "rank_stage_c_2026-09-23_{v}.json", "scored", "candidates"),
    "C_no_residual": (FITS / "rank_stage_c_2026-09-23_{v}_no_residual.json", "scored", "candidates"),
}


def gate(d) -> bool:
    return d["ts"][K20] > 4 and d["halo_counts"][K20] > 0 and not d["quadrant_any_negative"][K20]


def set_spec(name, variant):
    path, key, kind = SETS[name]
    d = json.loads(Path(str(path).format(v=variant)).read_text())
    spec = d["consensus_selected"] if kind == "candidates" else d["structures"]
    rows = [r for r in d[key] if r["model"] != "gll_iem_v07" and r["variant"] in (variant, "-")]
    return kind, tuple(spec), rows


def stack_for(model, variant, kind, spec, sh):
    if kind == "structures":
        return ranking.assemble(model, variant, spec, EDGES)
    st = selection.base_stack(model, variant, EDGES, sh)
    for c in spec:
        st, _, _ = selection.add_candidate(st, c, model, variant, EDGES, sh)
    return st


def diagnostics(st):
    keep = halo_cell_mask()
    fit = fit_all_bins(restrict_stack(st, keep), halo_mask=keep)["per_bin"]
    norms = {n: [fit[k]["amplitudes"][n] for k in range(st["n_bins"])]
             for n in st["templates"] if n.startswith("galprop_")}
    return {"neg2lnl": [b["neg2lnl"] for b in fit], "norms": norms,
            "mean_abs_sigma_halo": [b["mean_abs_sigma_halo"] for b in fit]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default=None, help="comma list; default every ingested model")
    ap.add_argument("--tag", default="", help="output name suffix (the iso run of 09-23 overwrote the aniso file without one)")
    args = ap.parse_args()
    out_path = FITS / f"galp21_first_look_{DATE}{('_' + args.tag) if args.tag else ''}.json"
    models = args.models.split(",") if args.models else sorted(
        p.stem.split("__")[1] for p in ranking.CACHE.glob("galp21_cubes__*.npz"))
    missing = [m for m in models if not galp21.available(m)]
    if missing or not models:
        sys.exit(f"Galp21 cubes not ingested for {missing or 'any model'}: run scripts/galp21_ingest.py first "
                 f"(the archives are available from the authors on request)")
    t0 = time.time()
    sh = ranking.shared_cells(EDGES)
    out = {"date": DATE, "models": models, "rows": []}
    for m in models:
        code = galp21.code(m)
        for v in ranking.VARIANTS:
            row = {"model": code, "variant": v, "sets": {}}
            base = selection.base_stack(code, v, EDGES, sh)
            row["bare"] = diagnostics(base)
            # the 2012 rows' bare -2lnL, for the ranking
            ref = [float(selection.bg_neg2lnl(selection.base_stack(c, v, EDGES, sh)).sum())
                   for c in ranking.MODEL_CODES]
            ref.append(float(selection.bg_neg2lnl(selection.base_stack("a2020", "-", EDGES, sh)).sum()))
            tot = float(sum(row["bare"]["neg2lnl"]))
            row["bare"].update({"neg2lnl_total": tot, "rank": 1 + sum(x < tot for x in ref), "of": len(ref) + 1,
                                "best_2012_or_a2020": min(ref), "median_2012": float(np.median(ref[:-1]))})
            print(f"{code} {v}: bare -2lnL rank {row['bare']['rank']}/{row['bare']['of']} ({time.time()-t0:.0f}s)",
                  flush=True)
            for sname in SETS:
                kind, spec, refrows = set_spec(sname, v)
                st = stack_for(code, v, kind, spec, sh)
                d = ranking.score(st, sh["J"])
                sv = dmspectra.joint_fit(st, sh["J"], "b", 500, EDGES, n_scan=31)
                tot = float(sum(d["neg2lnl_bg_only"]))
                ref_tot = [r["neg2lnl_total"] for r in refrows]
                ref_halo = [r["halo_counts"][K20] for r in refrows]
                row["sets"][sname] = {
                    "spec": list(spec), "neg2lnl_total": tot, "rank": 1 + sum(x < tot for x in ref_tot),
                    "of": len(ref_tot) + 1, "delta_vs_best_ref": tot - min(ref_tot),
                    "delta_vs_median_ref": tot - float(np.median(ref_tot)),
                    "halo_20": d["halo_counts"][K20], "halo_20_lo": d["halo_lo"][K20], "halo_20_hi": d["halo_hi"][K20],
                    "ts_20": d["ts"][K20], "halo_counts": d["halo_counts"], "ts": d["ts"],
                    "quadrant_counts_20": {q: c[K20] for q, c in d["quadrant_counts"].items()},
                    "favor": bool(gate(d)), "north_20": d["halo_north"][K20], "south_20": d["halo_south"][K20],
                    "injection_slope_min": min(d["injection_slopes"]),
                    "ref_halo_20_median": float(np.median(ref_halo)),
                    "ref_halo_20_range": [float(min(ref_halo)), float(max(ref_halo))],
                    "ref_favor": sum(bool(r.get("favor", gate(r))) for r in refrows), "ref_n": len(refrows),
                    "sigmav": {k: sv[k] for k in ("sigmav_best", "sigmav_lo68", "sigmav_hi68", "sigmav_ul95",
                                                  "ts_vs_zero")},
                    "diag": diagnostics(st)}
                s = row["sets"][sname]
                print(f"   {sname:16s} rank {s['rank']:3d}/{s['of']}  halo20 {s['halo_20']:7.0f} "
                      f"[{s['halo_20_lo']:.0f}, {s['halo_20_hi']:.0f}] TS {s['ts_20']:6.1f} favor {s['favor']}  "
                      f"sv {sv['sigmav_best']:+.2e} (TS {sv['ts_vs_zero']:.1f})  "
                      f"| 2012 halo median {s['ref_halo_20_median']:.0f} ({time.time()-t0:.0f}s)", flush=True)
            for label, cands in (("with_residual", None),
                                 ("no_residual", [c for c in selection.CANDIDATES if c != "bubble_totani"])):
                g = selection.greedy_for_model(code, v, EDGES, sh, candidates=cands)
                st = stack_for(code, v, "candidates", tuple(g["chosen"]), sh)
                d = ranking.score(st, sh["J"])
                sv = dmspectra.joint_fit(st, sh["J"], "b", 500, EDGES, n_scan=31)
                row[f"own_{label}"] = {"chosen": g["chosen"], "final_neg2lnl": g["final_neg2lnl"],
                                       "steps": [{k: s[k] for k in ("best", "delta", "threshold", "accepted")}
                                                 for s in g["steps"]],
                                       "halo_20": d["halo_counts"][K20], "halo_20_lo": d["halo_lo"][K20],
                                       "halo_20_hi": d["halo_hi"][K20], "ts_20": d["ts"][K20],
                                       "favor": bool(gate(d)), "north_20": d["halo_north"][K20],
                                       "south_20": d["halo_south"][K20],
                                       "sigmav_best": sv["sigmav_best"], "sigmav_ts": sv["ts_vs_zero"],
                                       "sigmav_lo68": sv["sigmav_lo68"], "sigmav_hi68": sv["sigmav_hi68"]}
                o = row[f"own_{label}"]
                print(f"   own path {label:13s}: {' > '.join(g['chosen'])}  halo20 {o['halo_20']:.0f} TS {o['ts_20']:.1f} "
                      f"favor {o['favor']}  sv {o['sigmav_best']:+.2e} (TS {o['sigmav_ts']:.1f})", flush=True)
            out["rows"].append(row)
            out_path.write_text(json.dumps(out, indent=1))
    print(f"done {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
