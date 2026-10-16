#!/usr/bin/env python
"""The nine-panel ranking figure for Appendix C (Kev, 2026-10-15): the signed
20.8 GeV significance of the NFW component against the BIC-penalized background
fit quality of the diffuse configuration, for the three Fermi bubbles
descriptions (columns) at the original 10-deg binning over nine bins and at
2 and 1 deg with the response folded in over 5.6-130 GeV (rows).

    python scripts/plot_rank_all.py              # results/rank_all_figure_<date>.{png,pdf}
    python scripts/plot_rank_all.py --no-title

Row 1 is `plot_halo_significance.py` (its loaders are reused); rows 2 and 3
read results/fits/rank_fine_cells{2,1}_psfuc_<date>.json from `rank_fine.py`.
Conventions as in the 10-deg figure: x is the margin, -(1/2) of the
BIC-penalized -2lnL relative to the best 2012-grid or A2020 configuration in
the panel; y is sign(A) sqrt(TS) at 20.8 GeV; colour is the 20.8 GeV NFW
counts; hollow markers fail the quadrant test; rings mark Totani's
configurations, the star A2020, diamonds the Galp21 isotropic-IC rows. Drawn
at the printed width with labels >= 10 pt and no tight bbox.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import plot_halo_significance as p10  # noqa: E402

DATE = "2026-10-15"
TITLES = (r"residual (4.3\,GeV) bubbles", r"flat 2014 mask", r"structured 2017 map")
ROWS = ((r"$10^\circ$ cells, nine bins, no PSF", None),
        (r"$2^\circ$ cells with PSF, 5.6--130\,GeV", 2),
        (r"$1^\circ$ cells with PSF, 5.6--130\,GeV", 1))


def rows_10deg(panel):
    """The 10-deg rows in the common format: grid (2012 + A2020) and Galp21 isotropic."""
    title, path, key = p10.PANELS[panel]
    rows, _ = p10.load(path, key, True)
    grid = [dict(r, favor=(r["sig"] > 2.0 and not r["quad_fail"])) for r in rows if r["model"] != "gll_iem_v07"]
    g21 = [dict(r, favor=(r["sig"] > 2.0 and not r["quad_fail"])) for r in p10.load_galp21(panel, True)]
    return grid, g21


def rows_fine(cell, bub):
    d = json.loads((ROOT / "results" / "fits" / f"rank_fine_cells{cell:g}_psfuc_{DATE}.json").read_text())
    out = []
    for r in d["panels"][bub]:
        out.append({"model": r["model"], "variant": r["variant"], "score": r["bic"], "sig": r["sig20"], "halo": r["halo20"],
                    "quad_fail": any(v <= 0 for v in r["quadrants20"].values()), "favor": r["favor"]})
    grid = [r for r in out if not r["model"].startswith("G21_")]
    g21 = [r for r in out if r["model"].endswith("+iso")]
    return grid, g21


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--no-title", action="store_true"); ap.add_argument("--vmax", type=float, default=9000)
    ap.add_argument("--draft", action="store_true", help="watermark: the finer rows are the cached subset, not the full grid")
    args = ap.parse_args()

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
    # 8.0 in tall: with its caption it fills one revtex4-2 page (text height 9.25 in) at \textwidth
    fig, axes = plt.subplots(3, 3, figsize=(W, 7.6), sharey="row")
    fig.subplots_adjust(left=0.085, right=0.875, top=0.955 if args.no_title else 0.945, bottom=0.175, wspace=0.10, hspace=0.27)
    vmin, vmax = 0, args.vmax
    counts = {}
    for i, (rlabel, cell) in enumerate(ROWS):
        for j, title in enumerate(TITLES):
            ax = axes[i, j]
            bub = ("residual", "flat", "structured")[j]
            grid, g21 = rows_10deg(j) if cell is None else rows_fine(cell, bub)
            best = min(r["score"] for r in grid)
            for r in grid + g21:
                r["margin"] = -0.5 * (r["score"] - best)
            n_fav = sum(r["favor"] for r in grid); counts[(i, j)] = (n_fav, len(grid), sum(r["favor"] for r in g21), len(g21))
            for variant, mk, lab in (("gas1_ics1", "o", r"GALPROP, 1 gas $+$ 1 ICS"), ("gas5_ics3", "^", r"GALPROP, 5 rings $+$ 3 ICS")):
                for filled in (True, False):
                    sub = [r for r in grid if r["variant"] == variant and r["favor"] == filled]
                    if not sub:
                        continue
                    xs = [r["margin"] for r in sub]; ys = [r["sig"] for r in sub]
                    if filled:
                        ax.scatter(xs, ys, c=[r["halo"] for r in sub], cmap="viridis", vmin=vmin, vmax=vmax, marker=mk, s=14,
                                   edgecolor="k", linewidth=0.4, label=lab)
                    else:
                        ax.scatter(xs, ys, facecolors="none", marker=mk, s=14, edgecolor="k", linewidth=0.6, label="hollow")
            for r in grid:
                mv = (r["model"], r["variant"])
                if mv == p10.TOTANI_BASELINE:
                    ax.scatter([r["margin"]], [r["sig"]], marker="o", s=85, facecolors="none", edgecolor="crimson", linewidth=1.4,
                               zorder=5, label=r"Totani's baseline model")
                elif mv in p10.TOTANI_VARIANTS:
                    ax.scatter([r["margin"]], [r["sig"]], marker="o", s=60, facecolors="none", edgecolor="crimson", linewidth=0.8,
                               zorder=5, label=r"Totani's \S3.4.4 variants")
                elif r["model"] == "a2020":
                    ax.scatter([r["margin"]], [r["sig"]], marker="*", s=75, facecolors="k" if r["favor"] else "none", edgecolor="k",
                               linewidth=0.7, zorder=6, label=r"A2020 ring model" if r["favor"] else "hollow")
            for filled in (True, False):
                sub = [r for r in g21 if r["favor"] == filled]
                if not sub:
                    continue
                xs = [r["margin"] for r in sub]; ys = [r["sig"] for r in sub]
                if filled:
                    ax.scatter(xs, ys, c=[r["halo"] for r in sub], cmap="viridis", vmin=vmin, vmax=vmax, marker="D", s=18,
                               edgecolor="darkorange", linewidth=0.8, zorder=4, label=r"Galp21 (v57, isotropic IC)")
                else:
                    ax.scatter(xs, ys, facecolors="none", marker="D", s=18, edgecolor="darkorange", linewidth=0.8, zorder=4, label="hollow")
            xmin = min(r["margin"] for r in grid + g21)
            ax.set_xlim(xmin * 1.10, -xmin * 0.08)
            ax.axhline(0, color="0.4", lw=0.9); ax.axhline(2, color="0.6", lw=0.7, ls=":")
            if i == 0:
                ax.set_title(title, fontsize=FS + 0.5)
            # the empty corner: lower left for the residual and flat panels, upper right for the structured ones
            corner = (0.97, 0.95, "right", "top") if j == 2 else (0.04, 0.08, "left", "bottom")
            ax.text(corner[0], corner[1], rf"{n_fav}/{len(grid)} favor a halo", transform=ax.transAxes, ha=corner[2], va=corner[3],
                    fontsize=FS - 0.5, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.2))
            ax.xaxis.set_major_locator(MaxNLocator(3, prune="upper" if j < 2 else None))
            ax.grid(alpha=0.3)
        axes[i, 0].set_ylabel(rlabel + "\n" + r"sign$\,\times\sqrt{\mathrm{TS}}$ at 20.8\,GeV", fontsize=FS)
    handles, labels = [], []
    for ax in axes.flat:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels and l != "hollow":
                handles.append(h); labels.append(l)
    keep = list(zip(handles, labels))
    keep.append((Line2D([], [], marker="o", ls="none", mfc="none", mec="k", ms=4.5, mew=0.6), r"hollow: quadrant test fails"))
    order = [r"GALPROP, 1 gas $+$ 1 ICS", r"GALPROP, 5 rings $+$ 3 ICS", r"A2020 ring model", r"Galp21 (v57, isotropic IC)",
             r"Totani's baseline model", r"Totani's \S3.4.4 variants", r"hollow: quadrant test fails"]
    keep.sort(key=lambda x: order.index(x[1]) if x[1] in order else 99)
    fig.legend([h for h, _ in keep], [l for _, l in keep], loc="upper center", ncol=3, bbox_to_anchor=(0.48, 0.085),
               handletextpad=0.3, columnspacing=1.0, borderaxespad=0.2)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cax = fig.add_axes([0.888, 0.165, 0.014, 0.79]); cb = fig.colorbar(sm, cax=cax)
    cb.set_label(r"NFW counts at 20.8\,GeV", fontsize=FS); cb.ax.tick_params(labelsize=FS)
    fig.text(0.48, 0.115, r"$\left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm model} - \left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm best}$"
             r"\quad(BIC-penalized fit quality)", ha="center", va="center", fontsize=FS)
    if not args.no_title:
        fig.suptitle(r"Halo significance at 20.8\,GeV against diffuse-model fit quality, three bubbles templates and three binnings", fontsize=FS, y=0.995)
    if args.draft:
        fig.text(0.48, 0.56, "[draft]", ha="center", va="center", fontsize=110, color="crimson", alpha=0.18, rotation=30, zorder=100)
    tag = ("_notitle" if args.no_title else "") + ("_draft" if args.draft else "")
    out = ROOT / "results" / f"rank_all_figure_{DATE}{tag}.png"
    fig.savefig(out, dpi=200); fig.savefig(out.with_suffix(".pdf"))
    print(f"saved {out} (+ .pdf)")
    for (i, j), (nf, n, gf, ng) in sorted(counts.items()):
        print(f"row {i} {['residual','flat','structured'][j]:10s}: {nf}/{n} grid+A2020 favor a halo; Galp21 iso {gf}/{ng}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
