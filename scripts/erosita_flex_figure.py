#!/usr/bin/env python
"""Figure: the halo's south/north ratio at 20.8 GeV per eROSITA lobe model,
West (native eROSITA) and East (reconstructed) halves, data vs spherical-halo
reference skies. PSF-folded 2-deg cells, flat mask.

    python scripts/erosita_flex_figure.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAG = "_" + os.environ.get("FINE_TAG", "psf")
SUFFIX = "" if TAG == "_psf" else TAG          # the CLEAN-table figures keep their names
DATE = "2026-09-27"
ORDER = ["none", "unif", "bright", "unif_outer", "bright_outer", "unif_EW", "bright_EW",
         "unif_flex", "bright_flex", "unif_mirrorW", "bright_mirrorW", "unif_Wonly", "bright_Wonly"]
LABEL = {"none": "none", "unif": "uniform N/S", "bright": "brightness N/S", "unif_outer": "uniform, outer",
         "bright_outer": "brightness, outer", "unif_EW": "uniform N/S×E/W", "bright_EW": "brightness N/S×E/W",
         "unif_flex": "uniform, 12 pieces", "bright_flex": "brightness, 12 pieces", "unif_mirrorW": "uniform, West mirrored",
         "bright_mirrorW": "brightness, West mirrored", "unif_Wonly": "uniform, West only", "bright_Wonly": "brightness, West only"}


def ratios(h):
    r = lambda s, n: s / n if abs(n) > 50 else np.nan
    return r(h["WS"] + h["ES"], h["WN"] + h["EN"]), r(h["WS"], h["WN"]), r(h["ES"], h["EN"])


def main() -> int:
    res = json.loads((ROOT / "results" / "fits" / f"erosita_flex_ns{TAG}_{DATE}.json").read_text())["results"]
    D = [r for r in res if r["kind"] == "data"]; S = [r for r in res if r["kind"] == "halo_flat"]
    # Drawn at its printed size (Kev, 2026-10-03): the manuscript includes it at
    # 0.9\textwidth = 459 pt (revtex4-2 prd), captions are 9 pt, so every label is
    # >= 10 pt and nothing rescales it. The key to the markers is in the caption.
    FS = 10
    plt.rc("text", usetex=True); plt.rc("font", family="serif", size=FS)
    plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)
    W = 0.9 * 510 / 72.27
    fig, axes = plt.subplots(1, 3, figsize=(W, 4.0), sharey=True)
    fig.subplots_adjust(left=0.30, right=0.985, top=0.93, bottom=0.14, wspace=0.08)
    y = np.arange(len(ORDER))[::-1]
    for ax, c, title in zip(axes, (0, 1, 2), ("whole region", r"West (native)", r"East (reconstr.)")):
        for yy, m in zip(y, ORDER):
            dr = np.array([ratios(r["models"][m]["halo"])[c] for r in D]); sr = np.array([ratios(r["models"][m]["halo"])[c] for r in S])
            ax.fill_betweenx([yy - 0.3, yy + 0.3], np.nanpercentile(sr, 16), np.nanpercentile(sr, 84), color="0.82", lw=0)
            ax.plot([np.nanmedian(sr)] * 2, [yy - 0.3, yy + 0.3], color="0.5", lw=0.9)
            ax.errorbar(np.nanmedian(dr), yy, xerr=[[np.nanmedian(dr) - np.nanpercentile(dr, 16)], [np.nanpercentile(dr, 84) - np.nanmedian(dr)]],
                        fmt="o", ms=4, color="C3" if m == "bright_flex" else "k", capsize=1.8, lw=0.9)
        ax.axvline(1.0, color="0.3", ls=":", lw=0.8)
        ax.set_xlim(0.5, 2.1); ax.set_xticks([1.0, 1.5, 2.0]); ax.set_xticklabels(["1", "1.5", "2"])
        ax.set_title(title, fontsize=FS)
    axes[0].set_yticks(y); axes[0].set_yticklabels([LABEL[m].replace("×", r"$\times$") for m in ORDER])
    fig.text(0.6425, 0.03, r"halo counts at 20.8\,GeV, south/north", ha="center", va="bottom", fontsize=FS)
    out = ROOT / "results" / f"erosita_flex_figure_{DATE}{SUFFIX}.pdf"
    fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=200)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
