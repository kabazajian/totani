#!/usr/bin/env python
"""The quadrant gate as a RANGE over a pre-registered family of companion sets
(Kev, 2026-09-21, item 2).

Every "x of N rows favor a halo" number is conditional on what else is in the
fit. Measured 2026-09-21 with the models held fixed, Oscar's structured panel
moved from 15/257 to 133/257 on the companion set alone, while the halo
amplitude and its significance barely moved (report_2026-09-21_explained §6.3).
So the number is reported as a range over a family fixed in advance:

    bubble   in {Totani residual, flat 2014, Oscar structured 2017}   (the 3 panels)
    Loop I   in {two shells, two shells + Haslam NPS}
    eROSITA  in {none, Oscar's uniform N/S, Oscar's brightness N/S}
    solar IC in {out, in}

12 companion sets per panel, 36 in all, every one scored on all 258 rows with
Stage C's scorer unchanged. Resumable: each set is written to its own JSON as
it completes and skipped on a re-run.

    python scripts/quadrant_family.py                 # everything (~12 h at 6 workers)
    python scripts/quadrant_family.py --panel flat    # one panel
    python scripts/quadrant_family.py --summarise     # rebuild the summary from what exists

Outputs `results/fits/quadrant_family/<date>__<panel>__<set>.json`, a
compact per-row record; and `results/quadrant_family_<date>.{json,md}` with
the range per panel.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rank_stage_c import _score, cached_models  # noqa: E402

DATE = "2026-09-23"
OUT = ROOT / "results" / "fits" / "quadrant_family"

PANELS = {
    "residual": "bubble_totani",
    "flat": "bubble_flat",
    "structured": "bubble_structured",
}
FAMILY = {
    "loopI": [("shells", "loopI_shells"), ("shells+nps", "loopI_shells+nps")],
    "erosita": [("none", None), ("unif", "erosita_unif_NS"), ("bright", "erosita_bright_NS")],
    "sun": [("nosun", None), ("sun", "sun_moon_ic")],
}


def companion_sets():
    for (ln, lc), (en, ec), (sn, sc) in itertools.product(*FAMILY.values()):
        name = f"{ln}__{en}__{sn}"
        yield name, tuple(c for c in (lc, ec, sc) if c is not None)


def gate(r):
    return r["ts"][5] > 4 and r["halo_counts"][5] > 0 and not r["quadrant_any_negative"][5]


def fav(r):
    """Sets scored before 2026-09-22 wrote the key as "favour"; read both."""
    return r["favor"] if "favor" in r else r["favour"]


def compact(r):
    return {"model": r["model"], "variant": r["variant"],
            "halo_counts": r["halo_counts"], "halo_lo": r.get("halo_lo"), "halo_hi": r.get("halo_hi"),
            "ts": r["ts"], "quadrant_any_negative": r["quadrant_any_negative"],
            "neg2lnl_total": r["neg2lnl_total"], "n_free_per_bin": r["n_free_per_bin"],
            "structures": r["structures"], "structures_unavailable": r.get("structures_unavailable", []),
            "favor": bool(gate(r))}


def run_set(panel, set_name, structures, workers):
    p = OUT / f"{DATE}__{panel}__{set_name}.json"
    if p.exists():
        return json.loads(p.read_text())
    t0 = time.time()
    models = cached_models("both") + [("a2020", "-"), ("gll_iem_v07", "-")]
    with Pool(workers) as pool:
        rows = pool.map(_score, [(m, v, structures) for m, v in models])
    out = {"date": DATE, "panel": panel, "set": set_name, "structures": list(structures),
           "rows": [compact(r) for r in rows], "seconds": time.time() - t0}
    OUT.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out))
    return out


def summarise():
    res = {}
    for p in sorted(OUT.glob(f"{DATE}__*.json")):
        d = json.loads(p.read_text())
        phys = [r for r in d["rows"] if r["model"] != "gll_iem_v07"]
        res.setdefault(d["panel"], {})[d["set"]] = {
            "structures": d["structures"], "n": len(phys),
            "favor": sum(fav(r) for r in phys),
            "sig2": sum(1 for r in phys if r["ts"][5] > 4 and r["halo_counts"][5] > 0),
            "median_halo": float(np.median([r["halo_counts"][5] for r in phys])),
            "by_variant": {v: {"n": sum(1 for r in phys if r["variant"] == v),
                               "favor": sum(fav(r) for r in phys if r["variant"] == v)}
                           for v in ("gas1_ics1", "gas5_ics3")},
            # per-row verdicts, so the figure can mark rows that pass under all / some / none
            "favor_by_row": {f"{r['model']}__{r['variant']}": fav(r) for r in phys},
        }
    summary = {"date": DATE, "family": {k: [n for n, _ in v] for k, v in FAMILY.items()}, "panels": {}}
    lines = [f"# The quadrant gate as a range over companion sets — {DATE}", "",
             "Family (pre-registered): Loop I ∈ {shells, shells+NPS} × eROSITA ∈ {none, Oscar uniform N/S, Oscar brightness N/S} "
             "× solar IC ∈ {out, in}; isotropic and point sources always in. 'favor' = >2σ AND quadrant-consistent at 20.8 GeV.", ""]
    for panel in PANELS:
        sets = res.get(panel, {})
        if not sets:
            continue
        favs = [s["favor"] for s in sets.values()]; n = next(iter(sets.values()))["n"]
        rows_all = set.intersection(*[set(k for k, f in s["favor_by_row"].items() if f) for s in sets.values()])
        rows_any = set.union(*[set(k for k, f in s["favor_by_row"].items() if f) for s in sets.values()])
        summary["panels"][panel] = {"n_sets": len(sets), "n_rows": n, "favor_min": min(favs), "favor_max": max(favs),
                                    "rows_pass_under_all": len(rows_all), "rows_pass_under_any": len(rows_any),
                                    "sets": sets}
        lines += [f"## {panel}: favor a halo in **{min(favs)}–{max(favs)} of {n}** rows "
                  f"({100*min(favs)/n:.0f}–{100*max(favs)/n:.0f}%) across {len(sets)} companion sets; "
                  f"{len(rows_all)} rows pass under every set, {len(rows_any)} under at least one", "",
                  "| companion set | templates | favor | >2σ | median halo | 1-gas | 5-ring |", "|---|---|---|---|---|---|---|"]
        for name, s in sorted(sets.items(), key=lambda kv: -kv[1]["favor"]):
            bv = s["by_variant"]
            lines.append(f"| {name} | {len(s['structures'])} | {s['favor']}/{s['n']} | {s['sig2']} | {s['median_halo']:,.0f} | "
                         f"{bv['gas1_ics1']['favor']}/{bv['gas1_ics1']['n']} | {bv['gas5_ics3']['favor']}/{bv['gas5_ics3']['n']} |")
        lines.append("")
    (ROOT / "results" / "fits" / f"quadrant_family_{DATE}.json").write_text(json.dumps(summary, indent=1))
    (ROOT / "results" / f"quadrant_family_{DATE}.md").write_text("\n".join(lines) + "\n")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--summarise", action="store_true")
    args = ap.parse_args()
    if not args.summarise:
        t0 = time.time()
        panels = [args.panel] if args.panel else list(PANELS)
        todo = [(pn, sn, (PANELS[pn],) + st) for pn in panels for sn, st in companion_sets()]
        print(f"{len(todo)} (panel, set) scorings, {args.workers} workers", flush=True)
        for i, (pn, sn, st) in enumerate(todo, 1):
            d = run_set(pn, sn, st, args.workers)
            phys = [r for r in d["rows"] if r["model"] != "gll_iem_v07"]
            print(f"[{i}/{len(todo)}] {pn:10s} {sn:22s} favor {sum(fav(r) for r in phys):3d}/{len(phys)}"
                  f"  ({d['seconds']:.0f}s; total {time.time()-t0:.0f}s)", flush=True)
    s = summarise()
    for pn, v in s["panels"].items():
        print(f"  {pn:10s}: {v['favor_min']}–{v['favor_max']} of {v['n_rows']} across {v['n_sets']} sets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
