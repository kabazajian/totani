#!/usr/bin/env python
"""Stage B: rank the diffuse models with a FIXED reference structure set.

    python scripts/rank_stage_b.py                     # every cached model, both variants
    python scripts/rank_stage_b.py --variant gas1_ics1 --workers 8

Reference structures (plan §5): the flat 2014 bubble mask and the two Wolleben
Loop I shells -- the Ackermann+ 2012 / Totani choice -- so the ranking is
about the diffuse model, not the structures. Isotropic and point sources are
always in. Each (model, variant) is scored the standard way in cell space
(`ranking.score`): background-only -2lnL per bin, per-cell mean |sigma|, the
NFW-rho^2 halo per bin with TS and 68% interval, the injection gate, the
quadrant diagnostic and the N/S split. The a2020 ring model and gll_iem_v07
are scored as reference rows.

Ranking metric: -2lnL summed over the 9 bins, relative to Totani's baseline
model (SL_Z6_R30_T150_C2, gas1_ics1). Everything is written to JSON and a
Markdown table; a Foster-style scatter (halo TS improvement vs diffuse-model
-2lnL relative to Totani's, colored by halo counts fraction) is saved as a
figure, with our axes.

Cell-space only: each score is ~40 s and ~100 MB; `--workers` runs them in
parallel. Nothing here touches the 960x960 grids.
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
from totani import ranking  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
REFERENCE_STRUCTURES = ("bubble_flat", "loopI_shell1", "loopI_shell2")
BASELINE = ("SL_Z6_R30_T150_C2", "gas1_ics1")
DATE = "2026-09-24"


def cached_models():
    out = []
    for p in sorted(ranking.CACHE.glob("model__*.npz")):
        _, code, variant = p.stem.split("__")
        if code in ranking.MODEL_CODES:          # the 2012 grid only: Galp21 rows (G21_*) share this cache
            out.append((code, variant))
    return out


def _score_one(args):
    model, variant, structures = args
    sh = ranking.shared_cells(EDGES)
    st = ranking.assemble(model, variant, structures, EDGES)
    d = ranking.score(st, sh["J"])
    d.update({"model": model, "variant": variant, "structures": list(structures),
              "neg2lnl_total": float(sum(d["neg2lnl_bg_only"]))})
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=("gas1_ics1", "gas5_ics3", "both"), default="both")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--structures", default=",".join(REFERENCE_STRUCTURES))
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    structures = tuple(args.structures.split(","))
    t0 = time.time()

    jobs = [(m, v, structures) for m, v in cached_models() if args.variant in ("both", v)]
    jobs += [("a2020", "-", structures), ("gll_iem_v07", "-", ("bubble_flat", "loopI_shell1", "loopI_shell2"))]
    print(f"{len(jobs)} (model, variant) rows with structures {structures}; {args.workers} workers", flush=True)
    with Pool(args.workers) as pool:
        rows = pool.map(_score_one, jobs)
    print(f"scored in {time.time()-t0:.0f}s", flush=True)

    base = next((r for r in rows if (r["model"], r["variant"]) == BASELINE), None)
    base_total = base["neg2lnl_total"] if base else None
    for r in rows:
        r["delta_neg2lnl_vs_totani_baseline"] = (r["neg2lnl_total"] - base_total) if base_total is not None else None
        r["halo_counts_20p8"] = r["halo_counts"][5]
        r["ts_20p8"] = r["ts"][5]
        r["ts_sum"] = float(sum(r["ts"]))
        tot = sum(sum(v) for v in [st for st in [r["halo_counts"]]])  # noqa: F841
    rows.sort(key=lambda r: r["neg2lnl_total"])

    out = {"date": DATE, "run": "fff2de5bc7e0", "structures": list(structures),
           "baseline": list(BASELINE), "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]), "rows": rows}
    tag = f"_{args.tag}" if args.tag else ""
    jp = ROOT / "results" / "fits" / f"rank_stage_b_{DATE}{tag}.json"
    jp.write_text(json.dumps(out, indent=2))

    lines = ["# Stage B ranking — " + DATE + (f" ({args.tag})" if args.tag else ""),
             "", f"Reference structures: {', '.join(structures)} + isotropic + point sources. "
             f"Δ(−2lnL) relative to Totani's baseline model ({BASELINE[0]}, {BASELINE[1]}), summed over 9 bins.", "",
             "| rank | model | variant | n | Δ(−2lnL) | ⟨|σ|⟩ 1.5 GeV | halo @ 20.8 | TS @ 20.8 | ΣTS | quad Σ/joint @ 20.8 | any neg | inj min | N / S @ 20.8 |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1):
        q = r["quadrant_sum_over_joint"][5]
        lines.append(f"| {i} | {r['model']} | {r['variant']} | {r['n_free_per_bin']} | "
                     f"{r['delta_neg2lnl_vs_totani_baseline']:+,.0f} | {r['mean_abs_sigma_halo'][0]:.2f} | "
                     f"{r['halo_counts_20p8']:,.0f} [{r['halo_lo'][5]:,.0f}, {r['halo_hi'][5]:,.0f}] | {r['ts_20p8']:.1f} | {r['ts_sum']:.0f} | "
                     f"{(q if q is not None else float('nan')):.2f} | {'yes' if r['quadrant_any_negative'][5] else 'no'} | "
                     f"{min(r['injection_slopes']):.2f} | {r['halo_north'][5]:,.0f} / {r['halo_south'][5]:,.0f} |")
    mp = ROOT / "results" / f"rank_stage_b_{DATE}{tag}.md"
    mp.write_text("\n".join(lines) + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5.5))
        pr = [r for r in rows if r["variant"] in ("gas1_ics1", "gas5_ics3")]
        for v, mk in (("gas1_ics1", "o"), ("gas5_ics3", "^")):
            rr = [r for r in pr if r["variant"] == v]
            if not rr:
                continue
            x = [r["delta_neg2lnl_vs_totani_baseline"] for r in rr]
            y = [r["ts_sum"] for r in rr]
            c = [r["halo_counts_20p8"] for r in rr]
            sc = ax.scatter(x, y, c=c, marker=mk, cmap="coolwarm", edgecolor="k", s=45, label=v)
        for r in rows:
            if r["variant"] == "-":
                ax.scatter([r["delta_neg2lnl_vs_totani_baseline"]], [r["ts_sum"]], marker="*", s=160,
                           color="k" if r["model"] == "a2020" else "gray", label=r["model"], zorder=3)
        ax.axvline(0, ls="--", color="r", lw=1)
        ax.set_xlabel("diffuse model −2lnL (bg-only, 9 bins) relative to Totani's baseline")
        ax.set_ylabel("halo TS summed over bins")
        ax.legend(fontsize=8)
        if pr:
            plt.colorbar(sc, ax=ax, label="halo counts @ 20.8 GeV")
        fig.tight_layout(); fig.savefig(ROOT / "results" / f"rank_stage_b_{DATE}{tag}.png", dpi=110)
    except Exception as e:  # matplotlib optional
        print(f"(no figure: {e})")

    print(f"\ntop 10 by -2lnL:")
    for r in rows[:10]:
        print(f"  {r['model']:22s} {r['variant']:9s} d(-2lnL) {r['delta_neg2lnl_vs_totani_baseline']:+9,.0f}  "
              f"halo@20.8 {r['halo_counts_20p8']:7,.0f} TS {r['ts_20p8']:6.1f}  quad {r['quadrant_sum_over_joint'][5] if r['quadrant_sum_over_joint'][5] is not None else float('nan'):+.2f}")
    print(f"\nsaved {jp}, {mp}; total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
