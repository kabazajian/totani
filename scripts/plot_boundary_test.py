#!/usr/bin/env python
"""The boundary-test figure for the September 8 report.

    python scripts/plot_boundary_test.py

Two panels sharing the energy axis: the lobe amplitude with its error bar, and
the significance. Both carry the same two series -- the flat 2014 mask and
Oscar's structured map -- because the comparison between them IS the result.

Identity is never color alone: each series also has its own marker and dash,
which is what makes the figure survive a grayscale print and color-vision
deficiency. Colours are slots 1 and 2 of the validated categorical palette, in
fixed order.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "fits" / "boundary_test_2026-09-08.json"
OUT = ROOT / "results" / "boundary_test_2026-09-08_fig.pdf"

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8983"
SERIES = {"ackermann": "#2a78d6", "macias": "#eb6834"}
LABEL = {"ackermann": "flat 2014 mask", "macias": "structured map (baseline)"}
MARKER = {"ackermann": "o", "macias": "s"}
DASH = {"ackermann": (0, (5, 2)), "macias": "solid"}


def main() -> int:
    d = json.loads(SRC.read_text())
    E = np.array(d["E_centers_gev"])

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.4, 6.0), sharex=True,
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.12})

    for ax in (ax1, ax2):
        ax.set_facecolor("white")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
            ax.spines[s].set_linewidth(0.8)
        ax.tick_params(colors=INK_2, labelsize=9, width=0.8)
        ax.grid(True, which="major", color="#e8e7e3", lw=0.7, zorder=0)
        ax.set_axisbelow(True)

    # ---- (a) amplitude -----------------------------------------------------
    for bub in ("ackermann", "macias"):
        rows = d["runs"][bub]["hemispheres"]["combined"]
        y = np.array([r["lobe"] for r in rows])
        e = np.array([r["sigma_excess"] for r in rows])
        pos = y > 0
        ax1.errorbar(E[pos], y[pos], yerr=e[pos], color=SERIES[bub],
                     marker=MARKER[bub], ms=5.5, mew=0, lw=1.8,
                     ls=DASH[bub], capsize=2.5, elinewidth=1.1, zorder=3,
                     label=LABEL[bub])
        # A point that went negative cannot sit on a log axis. Show it as an
        # open marker on the floor rather than dropping it silently; the value
        # is in the caption.
        if (~pos).any():
            ax1.plot(E[~pos], np.full((~pos).sum(), 40.0), marker="v", ms=5.5,
                     ls="none", mfc="white", mec=SERIES[bub], mew=1.4, zorder=3)
    ax1.set_yscale("log")
    ax1.set_ylim(32, 7e4)
    ax1.set_ylabel("lobe  (counts)", fontsize=10, color=INK)
    ax1.legend(frameon=False, fontsize=9, loc="lower left",
               labelcolor=INK_2, handlelength=2.6)
    ax1.set_title("Emission bounded by the eROSITA edge that the bubble\n"
                  "templates do not capture",
                  fontsize=11, color=INK, loc="left", pad=10)

    # ---- (b) significance --------------------------------------------------
    ax2.axhspan(-2, 2, color="#f0efeb", zorder=0)
    ax2.axhline(0, color=MUTED, lw=0.8, zorder=1)
    for bub in ("ackermann", "macias"):
        rows = d["runs"][bub]["hemispheres"]["combined"]
        s = np.array([r["lobe_sigma"] for r in rows])
        ax2.plot(E, s, color=SERIES[bub], marker=MARKER[bub], ms=5.5, mew=0,
                 lw=1.8, ls=DASH[bub], zorder=3)
    ax2.set_ylabel("significance ($\\sigma$)", fontsize=10, color=INK)
    ax2.set_xlabel("photon energy (GeV)", fontsize=10, color=INK)
    ax2.set_xscale("log")
    ax2.set_xlim(1.0, 130)
    ax2.set_ylim(-2.8, 20.5)
    ax2.set_xticks([1, 2, 5, 10, 20, 50, 100])
    ax2.set_xticklabels(["1", "2", "5", "10", "20", "50", "100"])
    ax2.annotate("consistent with zero", xy=(1.12, 0.2), fontsize=7.5,
                 color=MUTED, va="bottom")

    # The bin the whole analysis turns on.
    for ax in (ax1, ax2):
        ax.axvline(E[5], color=MUTED, lw=0.8, ls=(0, (1, 3)), zorder=1)
    ax1.annotate("20.8 GeV", xy=(E[5], 7e4), fontsize=8, color=INK_2,
                 ha="center", va="bottom", xytext=(0, 3),
                 textcoords="offset points")

    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    print(f"saved {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
