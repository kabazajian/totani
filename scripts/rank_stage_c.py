#!/usr/bin/env python
"""Stage C: greedy structure selection, then re-rank and score with the
selected set.

    python scripts/rank_stage_c.py                    # all cached models, both variants
    python scripts/rank_stage_c.py --variant gas1_ics1 --workers 6

Steps (`totani.selection`):
  0. build the per-model Totani-residual bubble templates (pixel-level, one
     at a time, cached);
  1. the CONSENSUS path: at each step every remaining candidate's Δ(−2lnL)
     is computed for every model (parallel over models). 🚨 SELECTION IS BY
     FRACTION OF ROWS, NOT BY MEDIAN (Kev, 2026-09-21): a candidate is
     eligible if it clears its own BIC threshold in >= FRACTION_CUTOFF of the
     rows, and the eligible candidate that is the largest improvement in the
     most rows is added, its morphology group retired; stop when nothing is
     eligible. The median across rows is recorded as a diagnostic only -- on
     this ensemble, bimodal by component split, it fell in the gap between
     the two populations (report_2026-09-21_explained §4). Run this per
     split with `--variant`. `gll_iem_v07` is scored along but excluded from
     every fraction (it is not a physical model);
  2. each model's OWN greedy path (parallel over models);
  3. the full scoring (halo, TS, intervals, injection, quadrants, N/S) of every
     model with the consensus set, ranked by −2lnL.

Outputs `results/fits/rank_stage_c_<date>.json`, `results/rank_stage_c_<date>.md`.
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
from totani import ranking, selection  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-23"
#: A candidate enters the consensus only if it clears its own BIC threshold in
#: at least this fraction of rows (Kev, 2026-09-21). Run per component split
#: (`--variant`): the two splits are different analyses and disagree.
FRACTION_CUTOFF = 0.8


def cached_models(variant):
    out = []
    for p in sorted(ranking.CACHE.glob("model__*.npz")):
        _, code, v = p.stem.split("__")
        # the 2012 grid only: Galp21 rows (G21_*) share this cache since 2026-09-22
        if code in ranking.MODEL_CODES and variant in ("both", v):
            out.append((code, v))
    return out


def _deltas_for_model(args):
    """One model: Δ per candidate given the current selected set."""
    model, variant, selected, cands = args
    sh = ranking.shared_cells(EDGES)
    st = selection.base_stack(model, variant, EDGES, sh)
    for c in selected:
        try:
            st, _, _ = selection.add_candidate(st, c, model, variant, EDGES, sh)
        except (FileNotFoundError, KeyError, ValueError):
            # a selected structure this model cannot have (gll_iem_v07 has no
            # residual construction): evaluate the rest without it, flagged
            pass
    cur = selection.bg_neg2lnl(st)
    out = {}
    for c in cands:
        try:
            _, v, _ = selection.add_candidate(st, c, model, variant, EDGES, sh)
            out[c] = float((cur - v).sum())
        except (FileNotFoundError, KeyError, ValueError) as e:
            out[c] = None
    return (model, variant), out


def _own_path(args):
    # 🚨 Before 2026-09-22 this took no candidate list, so the --skip-residual
    # runs' "own paths" all contained bubble_totani (every 09-21 no_residual
    # JSON: 129/129 rows). The consensus was never affected.
    model, variant, cands = args
    sh = ranking.shared_cells(EDGES)
    try:
        return selection.greedy_for_model(model, variant, EDGES, sh, candidates=list(cands))
    except Exception as e:  # a model without every candidate still gets a path
        return {"model": model, "variant": variant, "error": str(e)}


def _score(args):
    model, variant, structures = args
    sh = ranking.shared_cells(EDGES)
    st = selection.base_stack(model, variant, EDGES, sh)
    signs, unavailable = {}, []
    for c in structures:
        try:
            st, _, s = selection.add_candidate(st, c, model, variant, EDGES, sh)
        except (FileNotFoundError, KeyError, ValueError):
            unavailable.append(c); continue
        if s is not None:
            signs[c] = s
    d = ranking.score(st, sh["J"])
    d.update({"model": model, "variant": variant, "structures": [c for c in structures if c not in unavailable],
              "structures_unavailable": unavailable,
              "neg_template_signs": signs, "neg2lnl_total": float(sum(d["neg2lnl_bg_only"]))})
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=("gas1_ics1", "gas5_ics3", "both"), default="both")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--tag", default="")
    ap.add_argument("--skip-residual", action="store_true", help="leave bubble_totani out")
    ap.add_argument("--render-only", action="store_true", help="rewrite the Markdown from the existing JSON")
    args = ap.parse_args()
    t0 = time.time()
    if args.render_only:
        tag = f"_{args.tag}" if args.tag else ""
        render(json.loads((ROOT / "results" / "fits" / f"rank_stage_c_{DATE}{tag}.json").read_text()), tag)
        return 0
    models = cached_models(args.variant) + [("a2020", "-")]
    refs = [("gll_iem_v07", "-")]

    # 0. per-model Totani-residual templates (pixel level; one at a time)
    cands = list(selection.CANDIDATES)
    if args.skip_residual:
        cands.remove("bubble_totani")
    else:
        for m, v in models:
            try:
                selection.totani_residual_cells(m, v, EDGES)
                print(f"  residual templates {m:22s} {v:9s} ok ({time.time()-t0:.0f}s)", flush=True)
            except Exception as e:
                print(f"  residual templates {m:22s} {v:9s} FAILED: {e}", flush=True)

    # 1. consensus path
    selected, trace = [], []
    remaining = list(cands)
    with Pool(args.workers) as pool:
        while remaining:
            rows = dict(pool.map(_deltas_for_model, [(m, v, tuple(selected), tuple(remaining)) for m, v in models + refs]))
            phys = [rows[k] for k in rows if k[0] != "gll_iem_v07"]
            usable = [c for c in remaining if all(r.get(c) is not None for r in phys)]
            dropped = [c for c in remaining if c not in usable]
            if dropped:
                # 🚨 a candidate unavailable for ANY physical model would be silently
                # excluded from the consensus (this hid bubble_totani on the first
                # run, 2026-09-19, when the a2020 residual build failed). Say so.
                print(f"  WARNING: candidates {dropped} unavailable for some model; excluded from this step", flush=True)
            # 🚨 FRACTION-OF-ROWS RULE (Kev, 2026-09-21), replacing the median.
            # The rows are competing hypotheses about ONE sky, not independent
            # measurements, so no average of their Δ is a likelihood of
            # anything; and this ensemble is bimodal by construction (two
            # component splits), so the median fell in the gap between the two
            # populations and was pinned there -- see report_2026-09-21_explained
            # §4. Each row now answers two yes/no questions for itself:
            #   clears(c) : does c beat ITS OWN BIC threshold in this row?
            #   best(c)   : among the candidates most rows need, is c this
            #               row's largest improvement?
            # A candidate is ELIGIBLE if clears(c) >= FRACTION_CUTOFF; the one
            # selected is the eligible candidate with the largest best(c).
            # Stop when nothing is eligible. The median is still recorded, as a
            # diagnostic only, so the old and new rules can be compared.
            rule = selection.fraction_rule_step(phys, usable, FRACTION_CUTOFF)
            thr_c, clears, eligible = rule["threshold_per_candidate"], rule["clears_fraction"], rule["eligible"]
            best_frac, best, accepted = rule["best_fraction_among_eligible"], rule["best"], rule["accepted"]
            med = {c: float(np.median([r[c] for r in phys])) for c in usable}
            step = {"step": len(selected) + 1, "selected_so_far": list(selected), "dropped_candidates": dropped,
                    "rule": "fraction-of-rows", "fraction_cutoff": FRACTION_CUTOFF,
                    "clears_fraction": clears, "best_fraction_among_eligible": best_frac,
                    "eligible": eligible, "threshold_per_candidate": thr_c,
                    "median_delta_diagnostic_only": med,
                    "min_delta": {c: float(np.min([r[c] for r in phys])) for c in usable},
                    "max_delta": {c: float(np.max([r[c] for r in phys])) for c in usable},
                    "per_model": {f"{k[0]}__{k[1]}": v for k, v in rows.items()},
                    "best": best, "accepted": accepted}
            trace.append(step)
            others = sorted(usable, key=lambda c: (clears[c], best_frac.get(c, 0.0)), reverse=True)
            print(f"step {step['step']}: {'ACCEPT' if accepted else 'stop'} {best}"
                  f" (clears {100*clears[best]:.0f}% of rows, best in {100*best_frac.get(best, 0.0):.0f}% of eligible-row contests;"
                  f" median Δ {med[best]:,.0f} for reference); others: "
                  + ", ".join(f"{c} clears {100*clears[c]:.0f}%/best {100*best_frac.get(c, 0.0):.0f}%"
                              for c in others if c != best)[:400]
                  + f" ({time.time()-t0:.0f}s)", flush=True)
            if not accepted:
                break
            selected.append(best)
            remaining = [c for c in remaining if c != best and not any(
                best in g and c in g for g in selection.GROUPS.values())]

        # 2. per-model own paths
        own = pool.map(_own_path, [(m, v, tuple(cands)) for m, v in models])
        print(f"own paths done ({time.time()-t0:.0f}s)", flush=True)

        # 3. full scoring with the consensus set
        scored = pool.map(_score, [(m, v, tuple(selected)) for m, v in models + refs])
    base = next((r for r in scored if (r["model"], r["variant"]) == ("SL_Z6_R30_T150_C2", "gas1_ics1")), None)
    for r in scored:
        r["delta_neg2lnl_vs_totani_baseline"] = (r["neg2lnl_total"] - base["neg2lnl_total"]) if base else None
    scored.sort(key=lambda r: r["neg2lnl_total"])

    out = {"date": DATE, "run": "fff2de5bc7e0", "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]),
           "candidates": cands, "groups": selection.GROUPS, "threshold": "BIC: p ln(1080)",
           "selection_rule": f"fraction-of-rows: eligible if clears own BIC threshold in >= {FRACTION_CUTOFF:.0%} of rows; "
                             "select the eligible candidate that is the largest improvement in the most rows",
           "variant": args.variant,
           "consensus_selected": selected, "consensus_trace": trace, "own_paths": own, "scored": scored}
    tag = f"_{args.tag}" if args.tag else ""
    jp = ROOT / "results" / "fits" / f"rank_stage_c_{DATE}{tag}.json"
    jp.write_text(json.dumps(out, indent=2))
    render(out, tag)
    print(f"\nconsensus: {selected}\nsaved {jp}; total {time.time()-t0:.0f}s")
    return 0


def render(out: dict, tag: str) -> None:
    """The Markdown from the JSON (also `--render-only`)."""
    selected, trace, own, scored = out["consensus_selected"], out["consensus_trace"], out["own_paths"], out["scored"]
    n_rows = len([r for r in scored if r["model"] != "gll_iem_v07"])
    lines = [f"# Stage C — {DATE}{(' ('+tag.lstrip('_')+')') if tag else ''}", "",
             f"Consensus selection across {n_rows} model rows ({out.get('variant', 'both')}) by the "
             f"**fraction-of-rows rule**: a candidate is eligible if it clears its own BIC threshold (p·ln 1080 = 63 per "
             f"template) in ≥ {100*trace[0].get('fraction_cutoff', 0.8):.0f}% of rows; the eligible candidate that is the "
             f"largest improvement in the most rows is selected. Median Δ is shown for reference only. "
             f"Selected: **{' → '.join(selected) if selected else 'nothing'}**", "",
             "| step | chosen | clears BIC in | best among eligible in | median Δ (ref.) | runner-up (clears / best) |",
             "|---|---|---|---|---|---|"]
    for st in trace:
        cl, bf, med = st["clears_fraction"], st["best_fraction_among_eligible"], st["median_delta_diagnostic_only"]
        b = st["best"]
        others = sorted((c for c in cl if c != b), key=lambda c: (cl[c], bf.get(c, 0.0)), reverse=True)
        ru = others[0] if others else "-"
        lines.append(f"| {st['step']} | {b if st['accepted'] else '(stop) ' + b} | {100*cl[b]:.0f}% | "
                     f"{100*bf.get(b, 0.0):.0f}% | {med[b]:,.0f} | "
                     f"{ru} ({100*cl.get(ru, 0):.0f}% / {100*bf.get(ru, 0.0):.0f}%) |")
    lines += ["", "## Each model's own path", "", "| model | variant | path |", "|---|---|---|"]
    for o in own:
        lines.append(f"| {o['model']} | {o['variant']} | {' → '.join(o.get('chosen', [])) or o.get('error', '-')} |")
    lines += ["", f"## Ranking with the consensus set ({', '.join(selected)})", "",
              "| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(scored, 1):
        q = r["quadrant_sum_over_joint"][5]
        dtb = r["delta_neg2lnl_vs_totani_baseline"]        # None in a five-ring-only run: the baseline is a single-gas row
        lines.append(f"| {i} | {r['model']} | {r['variant']} | {r['n_free_per_bin']} | {(f'{dtb:+,.0f}' if dtb is not None else '-')} | "
                     f"{r['halo_counts'][5]:,.0f} [{r['halo_lo'][5]:,.0f}, {r['halo_hi'][5]:,.0f}] | {r['ts'][5]:.1f} | "
                     f"{(q if q is not None else float('nan')):.2f} | {'yes' if r['quadrant_any_negative'][5] else 'no'} | "
                     f"{min(r['injection_slopes']):.2f} | {r['halo_north'][5]:,.0f} / {r['halo_south'][5]:,.0f} |")
    mp = ROOT / "results" / f"rank_stage_c_{DATE}{tag}.md"
    mp.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
