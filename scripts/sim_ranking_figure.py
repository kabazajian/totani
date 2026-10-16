#!/usr/bin/env python
"""Oscar's item 7: Figure 1 for simulated skies, with and without a halo, on
a common model subset and the same ranking definition as the data.

    python scripts/sim_ranking_figure.py

From the Stage E skies (null_structured, halo_structured; 16 models x 2
decompositions + A2020 per sky): for the residual path and the structured
(no-residual) path, every (sky, row) point at margin = -(1/2) x BIC-penalized
(-2lnL(H0) - best in that sky) against sign x sqrt(TS) at 20.8 GeV, coloured
by decomposition. The same quantities for the DATA rows of the same 33-model
subset are drawn alongside. Spearman correlations of significance with margin
are printed per decomposition, for data, null skies and injected skies: if a
genuine halo also produces the anticorrelation, the trend does not
discriminate; if it does not, it does.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from totani.selection import CANDIDATES  # noqa: E402

DATE = "2026-09-24"
FITS = ROOT / "results" / "fits"
K20 = 5
N_DATA = 1080
N_DIFF = {"gas1_ics1": 2, "gas5_ics3": 8, "-": 14}


def points(rows, selected):
    n_set = sum(len(CANDIDATES[c]) for c in selected)
    pts = []
    for r in rows:
        p = N_DIFF[r["variant"]] + 2 + n_set
        pts.append({"variant": r["variant"], "score": r["neg2lnl_total"] + 9 * p * np.log(N_DATA),
                    "sig": np.sign(r["halo_counts"][K20]) * np.sqrt(max(r["ts"][K20], 0.0)), "halo": r["halo_counts"][K20]})
    best = min(p["score"] for p in pts)
    for p in pts:
        p["margin"] = -0.5 * (p["score"] - best)
    return pts


def data_points(path, key, subset):
    d = json.loads(Path(path).read_text())
    rows = [r for r in d[key] if (r["model"], r["variant"]) in subset]
    pts = []
    for r in rows:
        pts.append({"variant": r["variant"], "score": r["neg2lnl_total"] + 9 * r["n_free_per_bin"] * np.log(N_DATA),
                    "sig": np.sign(r["halo_counts"][K20]) * np.sqrt(max(r["ts"][K20], 0.0)), "halo": r["halo_counts"][K20]})
    best = min(p["score"] for p in pts)
    for p in pts:
        p["margin"] = -0.5 * (p["score"] - best)
    return pts


def main() -> int:
    from scipy.stats import spearmanr
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rc("text", usetex=True); plt.rc("text.latex", preamble=r"\usepackage{amsmath}"); plt.rc("font", family="serif", size=11); plt.rc("legend", fontsize=9, frameon=False)
    subset = {(c, v) for c in stage_e.fit_codes() for v in ("gas1_ics1", "gas5_ics3")} | {("a2020", "-")}
    data = {"residual": data_points(FITS / "quadrant_family" / "2026-09-23__residual__shells__unif__nosun.json", "rows", subset),
            "structured": data_points(FITS / "rank_stage_b_2026-09-24_structured_shells_sun.json", "rows", subset)}
    sims = {}
    for truth in ("null_structured", "halo_structured"):
        S = [json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(FITS / "stage_e" / f"{truth}__*.json")))]
        for branch in ("residual", "no_residual"):
            pts = []
            for s in S:
                for v in ("gas1_ics1", "gas5_ics3"):
                    b = s["splits"][v][branch]
                    pts += points(b["rows"], b["selected"])
            sims[(truth, branch)] = pts
    # Drawn at its printed size (Kev, 2026-10-02), as Fig. 1: \textwidth = 510 pt,
    # captions 9 pt, so every label is >= 10 pt and no tight bbox rescales it.
    from matplotlib.ticker import MaxNLocator
    from matplotlib.lines import Line2D
    FS = 10
    plt.rc("font", size=FS); plt.rc("legend", fontsize=FS); plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)
    W = 510 / 72.27
    fig, axes = plt.subplots(2, 3, figsize=(W, 5.6), sharey="row")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.955, bottom=0.215, hspace=0.36, wspace=0.10)
    cols = {"gas1_ics1": "tab:green", "gas5_ics3": "tab:purple", "-": "k"}
    summary = {}
    panels = [("data", "residual", "data, residual", data["residual"]),
              ("null_structured", "residual", "no halo, residual", sims[("null_structured", "residual")]),
              ("halo_structured", "residual", "halo, residual", sims[("halo_structured", "residual")]),
              ("data", "structured", "data, structured", data["structured"]),
              ("null_structured", "no_residual", "no halo, structured", sims[("null_structured", "no_residual")]),
              ("halo_structured", "no_residual", "halo, structured", sims[("halo_structured", "no_residual")])]
    # the Spearman coefficients go in an empty corner of each panel
    corners = [(0.03, 0.04, "left", "bottom"), (0.03, 0.96, "left", "top"), (0.03, 0.04, "left", "bottom"),
               (0.97, 0.96, "right", "top"), (0.03, 0.96, "left", "top"), (0.03, 0.04, "left", "bottom")]
    for ax, (src, br, title, pts), cn in zip(axes.ravel(), panels, corners):
        for v, c in cols.items():
            sub = [p for p in pts if p["variant"] == v]
            if sub:
                ax.scatter([p["margin"] for p in sub], [p["sig"] for p in sub], s=3 if src != "data" else 16, c=c,
                           alpha=0.35 if src != "data" else 0.9, edgecolor="none" if src != "data" else "k", linewidth=0.4)
        ax.axhline(0, color="0.4", lw=0.8); ax.axhline(2, color="0.6", lw=0.7, ls=":"); ax.grid(alpha=0.3)
        ax.set_title(title, fontsize=FS + 0.5)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        rho = {}
        for v in ("gas1_ics1", "gas5_ics3"):
            sub = [p for p in pts if p["variant"] == v]
            if len(sub) > 5:
                rho[v] = float(spearmanr([p["margin"] for p in sub], [p["sig"] for p in sub]).correlation)
        summary[f"{src}__{br}"] = rho
        ax.text(cn[0], cn[1], r"$\rho$: " + r" / ".join(f"${r:+.2f}$" for r in rho.values()), transform=ax.transAxes,
                ha=cn[2], va=cn[3], fontsize=FS, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.2))
    fig.text(0.535, 0.125, r"$\left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm row} - \left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm best}$"
             r"\quad(BIC-penalized fit quality)", ha="center", va="center", fontsize=FS)
    fig.text(0.012, 0.585, r"significance at 20.8\,GeV, sign$\,\times\sqrt{\mathrm{TS}}$", rotation=90, ha="center", va="center", fontsize=FS)
    handles = [Line2D([], [], marker="o", ls="none", mfc=cols[v], mec="k", mew=0.4, ms=4.5) for v in cols]
    labels = [r"1 gas $+$ 1 ICS", r"5 rings $+$ 3 ICS", r"A2020"]
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.535, 0.0), handletextpad=0.3, columnspacing=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "results" / f"sim_ranking_figure_{DATE}.{ext}", dpi=200)
    (FITS / f"sim_ranking_figure_{DATE}.json").write_text(json.dumps({"date": DATE, "spearman": summary}, indent=1))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
