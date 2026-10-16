#!/usr/bin/env python
"""The 10-deg programme rerun on 2-deg cells (Kev, 2026-09-26), for the
33-row stratified subset (16 models x 2 splits + A2020, the Stage E fit
ensemble) and the 24 Galp21 rows (6 models x 2 splits x 2 IC treatments).

    python scripts/stage_2deg.py --build              # library + 57 rows, one process, ~1 h, <= 10 GB
    python scripts/stage_2deg.py --select --workers 6 # Stage C per split, with/without residual; own paths; scoring

Stage C: the fraction-of-rows rule of rank_stage_c (eligible if a candidate
clears its own BIC threshold -- p ln n with n = 27,000 at 2 deg -- in >= 80%
of the 2012-subset rows of the split; the eligible candidate that is the
largest improvement in the most rows is selected). Galp21 rows do not vote;
they get their own greedy paths and are scored under the consensus sets.
Scoring (fine.score): TS, 20.8 GeV amplitude with interval, spectrum,
quadrant gate, joint b-bbar 500 GeV <sigma v>.
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
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine, ranking, selection  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
CELL = 2.0
DATE = "2026-09-26" if not fine.tag() else "2026-09-27"
CUTOFF = 0.8
G21 = [f"G21_{m}{ic}" for m in ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12") for ic in ("", "+iso")]
G21 = [g for g in G21 if ranking.available(g)]             # Galp21 rows only when their cubes are ingested


def rows_2012(variant):
    return [(c, variant) for c in stage_e.fit_codes()] + [("a2020", "-")]


def all_rows():
    out = [("a2020", "-")] + [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()]
    out += [(g, v) for g in G21 for v in ("gas1_ics1", "gas5_ics3")]
    out += [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"),      # Totani's baseline
            ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SS_Z10_R30_T150_C5", "gas5_ics3")]     # the simulation truths
    return list(dict.fromkeys(out))


_W: dict = {}


def _init():
    _W["lib"] = fine.load_library(CELL)


def _row(r):
    if r not in _W:
        _W[r] = fine.load_row(CELL, *r)
    return _W[r]


def _stack_with(r, selected):
    lib = _W["lib"]; row = _row(r); st = fine.base_stack(lib, row)
    for c in selected:
        st, _ = fine.add_candidate(st, c, lib, row, lib["keep"])
    return st


def _deltas(args):
    r, selected, cands = args
    lib = _W["lib"]; row = _row(r)
    st = _stack_with(r, selected); cur = fine.bg_neg2lnl(st, lib["keep"])
    out = {}
    for c in cands:
        _, v = fine.add_candidate(st, c, lib, row, lib["keep"])
        out[c] = float((cur - v).sum())
    return r, out


def _own(args):
    r, cands = args
    lib = _W["lib"]; row = _row(r); st = fine.base_stack(lib, row); cur = fine.bg_neg2lnl(st, lib["keep"])
    chosen, steps, rem = [], [], list(cands)
    while rem:
        d = {}
        for c in rem:
            _, v = fine.add_candidate(st, c, lib, row, lib["keep"]); d[c] = float((cur - v).sum())
        best = max(d, key=d.get); thr = fine.bic_threshold(9 * len(selection.CANDIDATES[best]), lib["n_data"])
        steps.append({"best": best, "delta": d[best], "threshold": thr, "accepted": d[best] > thr})
        if d[best] <= thr:
            break
        st, cur = fine.add_candidate(st, best, lib, row, lib["keep"]); chosen.append(best)
        rem = [c for c in rem if c != best and not any(best in g and c in g for g in selection.GROUPS.values())]
    return {"model": r[0], "variant": r[1], "chosen": chosen, "steps": steps}


def _score(args):
    r, selected = args
    d = fine.score(_stack_with(r, selected), _W["lib"], with_sigmav=True, edges_gev=EDGES)
    d.update({"model": r[0], "variant": r[1], "structures": list(selected)})
    return d


def build():
    t0 = time.time()
    psf_mode = bool(fine.tag())
    (fine.build_library_psf if psf_mode else fine.build_library)(CELL, EDGES); print(f"library{fine.tag()} ({time.time()-t0:.0f}s)", flush=True)
    for i, r in enumerate(all_rows(), 1):
        (fine.build_row_psf if psf_mode else fine.build_row)(CELL, r[0], r[1], EDGES)
        print(f"[{i}/{len(all_rows())}] {r} ({time.time()-t0:.0f}s)", flush=True)


def select(workers):
    t0 = time.time(); out = {"date": DATE, "cell_deg": CELL, "cutoff": CUTOFF, "splits": {}}
    lib = fine.load_library(CELL)
    with Pool(workers, initializer=_init) as pool:
        for v in ("gas1_ics1", "gas5_ics3"):
            voters = rows_2012(v)
            for branch in ("residual", "no_residual"):
                cands = [c for c in selection.CANDIDATES if branch == "residual" or c != "bubble_totani"]
                selected, trace, rem = [], [], list(cands)
                while rem:
                    rows = dict(pool.map(_deltas, [(r, tuple(selected), tuple(rem)) for r in voters]))
                    phys = list(rows.values())
                    thr = {c: fine.bic_threshold(9 * len(selection.CANDIDATES[c]), lib["n_data"]) for c in rem}
                    clears = {c: float(np.mean([p[c] > thr[c] for p in phys])) for c in rem}
                    elig = [c for c in rem if clears[c] >= CUTOFF]
                    if elig:
                        tops = [max(elig, key=lambda c: p[c]) for p in phys]
                        bf = {c: float(np.mean([t == c for t in tops])) for c in elig}; best, acc = max(bf, key=bf.get), True
                    else:
                        bf = {}; best, acc = max(clears, key=clears.get), False
                    trace.append({"best": best, "accepted": acc, "clears": clears, "best_fraction": bf,
                                  "per_row": {f"{r[0]}__{r[1]}": d for r, d in rows.items()}})
                    print(f"{v} {branch} step {len(trace)}: {'ACCEPT' if acc else 'stop'} {best} (clears {100*clears[best]:.0f}%, "
                          f"best {100*bf.get(best, 0):.0f}%) ({time.time()-t0:.0f}s)", flush=True)
                    if not acc:
                        break
                    selected.append(best)
                    rem = [c for c in rem if c != best and not any(best in g and c in g for g in selection.GROUPS.values())]
                g21_rows = [(g, v) for g in G21]
                own = pool.map(_own, [(r, tuple(cands)) for r in voters + g21_rows])
                scored = pool.map(_score, [(r, tuple(selected)) for r in voters + g21_rows])
                out["splits"].setdefault(v, {})[branch] = {"selected": selected, "trace": trace, "own_paths": own, "scored": scored}
                (ROOT / "results" / "fits" / f"stage_2deg{fine.tag()}_{DATE}.json").write_text(json.dumps(out, indent=1))
                print(f"{v} {branch}: selected {selected}; own paths and scoring done ({time.time()-t0:.0f}s)", flush=True)
    render(out)


def render(out):
    lines = [f"# The 10-degree programme on 2-degree cells ({DATE})", "",
             "Fit ensemble: 16 stratified 2012 models per split + A2020 (the Stage E subset); Galp21 rows (6 models × 2 IC) scored and given own paths, not voting. "
             "BIC threshold p ln n with n = 27,000. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive.", ""]
    for v, br in out["splits"].items():
        for b, d in br.items():
            sc = d["scored"]; s12 = [x for x in sc if not x["model"].startswith("G21")]; sg = [x for x in sc if x["model"].startswith("G21")]
            def summ(S):
                if not S:
                    return "-"
                h = [x["halo_counts"][5] for x in S]; sv = [x["sigmav"]["sigmav_best"] * 1e25 for x in S]
                return (f"favor {sum(x['favor'] for x in S)}/{len(S)}; halo20 median {np.median(h):,.0f} [{min(h):,.0f}–{max(h):,.0f}]; "
                        f"⟨σv⟩ median {np.median(sv):+.1f}, 2σ+ {sum(x['sigmav']['ts_vs_zero'] > 4 and x['sigmav']['sigmav_best'] > 0 for x in S)}/{len(S)}")
            own_b = {}
            for o in d["own_paths"]:
                bub = next((c for c in o["chosen"] if c.startswith("bubble_")), "none")
                key = "G21" if o["model"].startswith("G21") else "2012"
                own_b.setdefault(key, {}); own_b[key][bub] = own_b[key].get(bub, 0) + 1
            lines += [f"## {v}, {b.replace('_', ' ')}", "", f"Selected: **{' → '.join(d['selected']) or 'nothing'}**", "",
                      f"- 2012 subset: {summ(s12)}", f"- Galp21: {summ(sg)}", f"- own paths, bubble chosen: {own_b}", ""]
    (ROOT / "results" / f"stage_2deg{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true"); ap.add_argument("--select", action="store_true")
    ap.add_argument("--render", action="store_true"); ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.build:
        build()
    if a.select:
        select(a.workers)
    if a.render:
        render(json.loads((ROOT / "results" / "fits" / f"stage_2deg{fine.tag()}_{DATE}.json").read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
