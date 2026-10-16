#!/usr/bin/env python
"""The finer-binning ranking figure (Kev, 2026-10-15): the signed 20.8 GeV
significance of the NFW component against the BIC-penalized background fit
quality of the diffuse configuration, response-folded, one panel per Fermi
bubbles description -- the 10-deg figure of `plot_halo_significance.py`
repeated on finer cells over 5.6-130 GeV.

    python scripts/plot_rank_fine.py --cell 2            # reads results/fits/rank_fine_cells2_psfuc_<date>.json
    python scripts/plot_rank_fine.py --cell 1 --no-title

Same conventions as the 10-deg figure: x is the margin, -(1/2) of the
BIC-penalized -2lnL relative to the best configuration in the panel (p =
templates x K, n = cells x K); y is sign(A) sqrt(TS) at 20.8 GeV; colour is
the 20.8 GeV NFW counts; hollow markers fail the quadrant test; rings mark
Totani's configurations, the star A2020, diamonds Galp21 (isotropic IC; the
anisotropic rows are drawn as hollow-edged diamonds in a second colour).
Written at the printed width, labels >= 10 pt, no tight bbox.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-10-15"
TOTANI_BASELINE = ("SL_Z6_R30_T150_C2", "gas1_ics1")
TOTANI_VARIANTS = {("SS_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z10_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R20_T150_C2", "gas1_ics1"),
                   ("SL_Z6_R30_T100000_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C5", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3")}
TITLES = {"residual": r"residual (4.3\,GeV) bubbles", "flat": r"flat 2014 mask", "structured": r"structured 2017 map"}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, required=True)
    ap.add_argument("--tag", default="psfuc"); ap.add_argument("--date", default=DATE)
    ap.add_argument("--no-title", action="store_true"); ap.add_argument("--vmax", type=float, default=9000)
    args = ap.parse_args()
    src = ROOT / "results" / "fits" / f"rank_fine_cells{args.cell:g}_{args.tag}_{args.date}.json"
    d = json.loads(src.read_text())
    K = len(d["bins"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    plt.rc("text", usetex=True); plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
    FS = 10
    plt.rc("font", family="serif", size=FS); plt.rc("legend", fontsize=FS, frameon=False)
    plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)

    W = 510 / 72.27
    fig, axes = plt.subplots(1, 3, figsize=(W, 4.1), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.865, top=0.93, bottom=0.36, wspace=0.12)
    vmin, vmax = 0, args.vmax
    summary = []
    for ax, bub in zip(axes, ("residual", "flat", "structured")):
        rows = d["panels"].get(bub, [])
        if not rows:
            ax.set_title(f"{TITLES[bub]} (missing)"); continue
        grid = [r for r in rows if not r["model"].startswith("G21_") and r["model"] != "a2020"]
        g21 = [r for r in rows if r["model"].startswith("G21_")]
        a20 = [r for r in rows if r["model"] == "a2020"]
        phys = grid + a20
        n_fav = sum(1 for r in phys if r["favor"]); n_any = sum(1 for r in phys if r["sig20"] > 2.0)
        summary.append((bub, n_fav, n_any, len(phys), len(g21)))
        for variant, mk, lab in (("gas1_ics1", "o", r"GALPROP, 1 gas $+$ 1 ICS"), ("gas5_ics3", "^", r"GALPROP, 5 rings $+$ 3 ICS")):
            for filled in (True, False):
                sub = [r for r in grid if r["variant"] == variant and r["favor"] == filled]
                if not sub:
                    continue
                xs = [r["margin"] for r in sub]; ys = [r["sig20"] for r in sub]
                if filled:
                    ax.scatter(xs, ys, c=[r["halo20"] for r in sub], cmap="viridis", vmin=vmin, vmax=vmax,
                               marker=mk, s=16, edgecolor="k", linewidth=0.4, label=lab)
                else:
                    ax.scatter(xs, ys, facecolors="none", marker=mk, s=16, edgecolor="k", linewidth=0.6,
                               label=r"hollow: quadrant test fails")
        for r in grid:
            mv = (r["model"], r["variant"])
            if mv == TOTANI_BASELINE:
                ax.scatter([r["margin"]], [r["sig20"]], marker="o", s=95, facecolors="none", edgecolor="crimson",
                           linewidth=1.5, zorder=5, label=r"Totani's baseline model")
            elif mv in TOTANI_VARIANTS:
                ax.scatter([r["margin"]], [r["sig20"]], marker="o", s=70, facecolors="none", edgecolor="crimson",
                           linewidth=0.8, zorder=5, label=r"Totani's \S3.4.4 variants")
        for r in a20:
            ax.scatter([r["margin"]], [r["sig20"]], marker="*", s=80, facecolors="k" if r["favor"] else "none",
                       edgecolor="k", linewidth=0.7, zorder=6, label=r"A2020 ring model" if r["favor"] else r"hollow: quadrant test fails")
        for iso, color, lab in ((True, "darkorange", r"Galp21 (v57, isotropic IC)"), (False, "firebrick", r"Galp21 (v57, anisotropic IC)")):
            for filled in (True, False):
                sub = [r for r in g21 if r["model"].endswith("+iso") == iso and r["favor"] == filled]
                if not sub:
                    continue
                xs = [r["margin"] for r in sub]; ys = [r["sig20"] for r in sub]
                if filled:
                    ax.scatter(xs, ys, c=[r["halo20"] for r in sub], cmap="viridis", vmin=vmin, vmax=vmax, marker="D", s=20,
                               edgecolor=color, linewidth=0.8, zorder=4, label=lab)
                else:
                    ax.scatter(xs, ys, facecolors="none", marker="D", s=20, edgecolor=color, linewidth=0.8, zorder=4,
                               label=r"hollow: quadrant test fails")
        xmin = min(r["margin"] for r in rows)
        ax.set_xlim(xmin * 1.10, -xmin * 0.08)
        ax.axhline(0, color="0.4", lw=0.9); ax.axhline(2, color="0.6", lw=0.7, ls=":")
        ax.set_title(TITLES[bub], fontsize=FS + 0.5)
        corner = (0.04, 0.10, "left", "bottom") if bub == "flat" else (0.97, 0.96, "right", "top")
        ax.text(corner[0], corner[1], rf"{n_fav}/{len(phys)} favor a halo", transform=ax.transAxes, ha=corner[2], va=corner[3],
                fontsize=FS, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
        ax.xaxis.set_major_locator(MaxNLocator(3, prune="upper" if ax is not axes[-1] else None))
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"significance at 20.8\,GeV, sign$\,\times\sqrt{\mathrm{TS}}$", fontsize=FS)
    handles, labels = [], []
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h); labels.append(l)
    keep = [(h, l) for h, l in zip(handles, labels) if not l.startswith("hollow")]
    keep.append((Line2D([], [], marker="o", ls="none", mfc="none", mec="k", ms=4.5, mew=0.6), r"hollow: quadrant test fails"))
    order = [r"GALPROP, 1 gas $+$ 1 ICS", r"GALPROP, 5 rings $+$ 3 ICS", r"A2020 ring model", r"Galp21 (v57, isotropic IC)",
             r"Galp21 (v57, anisotropic IC)", r"Totani's baseline model", r"Totani's \S3.4.4 variants", r"hollow: quadrant test fails"]
    keep.sort(key=lambda x: order.index(x[1]) if x[1] in order else 99)
    fig.legend([h for h, _ in keep], [l for _, l in keep], loc="upper center", ncol=3, bbox_to_anchor=(0.47, 0.175),
               handletextpad=0.3, columnspacing=1.0, borderaxespad=0.2)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cax = fig.add_axes([0.878, 0.36, 0.015, 0.57]); cb = fig.colorbar(sm, cax=cax)
    cb.set_label(r"halo counts at 20.8\,GeV", fontsize=FS); cb.ax.tick_params(labelsize=FS)
    fig.text(0.47, 0.25, r"$\left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm model} - \left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm best}$"
             r"\quad(BIC-penalized fit quality, 5.6--130\,GeV)", ha="center", va="center", fontsize=FS)
    if not args.no_title:
        n_models = len({r["model"] for r in d["panels"]["flat"] if not r["model"].startswith("G21_") and r["model"] != "a2020"})
        fig.suptitle(rf"Halo significance at 20.8\,GeV against diffuse-model fit quality: ${args.cell:g}^\circ$ cells with the PSF, "
                     rf"{K} bins 5.6--130\,GeV, {n_models} GALPROP models $\times$ 2 splits, A2020, Galp21", fontsize=FS, y=0.995)
    tag = "_notitle" if args.no_title else ""
    out = ROOT / "results" / f"rank_fine_figure_cells{args.cell:g}_{args.tag}_{args.date}{tag}.png"
    fig.savefig(out, dpi=200); fig.savefig(out.with_suffix(".pdf"))
    print(f"saved {out} (+ .pdf)")
    for bub, nf, na, n, ng in summary:
        print(f"{bub}: {nf}/{n} favor a halo, {na} above 2 sigma; {ng} Galp21 rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
