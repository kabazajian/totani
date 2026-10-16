#!/usr/bin/env python
"""Figures for the reply to Oscar's note of 2026-09-29.

    python scripts/oscar_figures.py

(1) recovery_refined: item 3 -- recovered/injected slope per bin, templates
    fixed and rebuilt, refined estimator (n_refine 7) with the 09-19 lattice
    (n_refine 3) for comparison; and the -2lnL profiles of the 1.5, 2.6 and
    4.3 GeV bins at the original injection strength.
(2) spectra_vs_sims: item 2 -- the halo spectrum on 2-deg PSF-folded cells for
    the simulation's true row under each bubble description, with 68% profile
    intervals, against the simulated skies (halo-free and injected, the
    construction rebuilt per sky): medians and 16-84% bands.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

TAG = "_" + os.environ.get("FINE_TAG", "psf")
SUFFIX = "" if TAG == "_psf" else TAG          # the CLEAN-table figures keep their names
DATE = "2026-09-29"
E = np.array([1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0])
plt.rcParams.update({"font.size": 8, "font.family": "serif"})
INK = "#222222"; MUTED = "#8a8a8a"
C_FIXED = "#2b6cb0"; C_REBUILT = "#c05621"


def recovery():
    d = json.loads((ROOT / "results" / "fits" / f"injection_refined_{DATE}.json").read_text())
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.3), gridspec_kw={"width_ratios": [1.6, 1, 1, 1]})
    ax = axes[0]
    for sname, ls, lab in (("flat", "-", "flat spectrum"), ("totani_like", "--", "Totani-like")):
        s = d["spectra"][sname]["slopes"]
        for kind, col in (("fixed", C_FIXED), ("rebuilt", C_REBUILT)):
            y = np.array([np.nan if v is None else v for v in s[kind]["7"]])
            ax.plot(E, y, ls, color=col, marker="o", ms=3, lw=1.2, label=f"{kind}, {lab}")
            y3 = np.array([np.nan if v is None else v for v in s[kind]["3"]])
            ax.plot(E, y3, ls, color=col, lw=0.6, alpha=0.45)
    ax.axhline(1, color=MUTED, lw=0.6, ls=":"); ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xscale("log"); ax.set_ylim(-1.2, 1.2); ax.set_xlabel("energy [GeV]"); ax.set_ylabel("recovered / injected")
    ax.text(1.35, -1.1, "1.5 GeV rebuilt: −4.7", fontsize=6, color=C_REBUILT)
    ax.legend(fontsize=5, frameon=False, loc="center right"); ax.set_title("item 3: refined estimator (thin: old grid)", fontsize=7)
    for ax, k in zip(axes[1:], ("0", "1", "2")):
        c = d["spectra"]["flat"]["profiles"][k]
        for kind, col in (("fixed", C_FIXED), ("rebuilt", C_REBUILT)):
            x = np.array(c[kind]["counts"]); v = np.array(c[kind]["neg2lnl"]); v = v - v.min()
            ax.plot(x / 1e3, v, color=col, lw=1.2)
        inj = c["injected"] * d["roi_J"][int(k)]
        base = d["a0"]["7"]["a"][int(k)] * d["roi_J"][int(k)]
        ax.axvline((base + inj) / 1e3, color=INK, lw=0.6, ls=":")
        ax.set_ylim(0, 25); ax.set_title(f"{E[int(k)]:g} GeV profile", fontsize=7); ax.set_xlabel("halo counts [10³]")
    axes[1].set_ylabel("Δ(−2 ln L)")
    fig.tight_layout()
    out = ROOT / "results" / f"oscar_recovery_refined_{DATE}{SUFFIX}.pdf"; fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=160); print(out)


def spectra():
    from halo_spectra_compare import conversion
    conv = conversion(); e2 = lambda c: np.asarray(c, float) * conv * (E * 1e3) ** 2 / 1e3 * 1e7
    sims = json.loads((ROOT / "results" / "fits" / f"sims_calibration_cells2{TAG}_{DATE}.json").read_text())["results"]
    data = json.loads((ROOT / "results" / "fits" / f"oscar_conditional_cells2{TAG}_{DATE}.json").read_text())["results"]
    # Drawn at its printed size (Kev, 2026-10-03): the manuscript includes it at
    # 0.95\textwidth = 484.5 pt (revtex4-2 prd), captions 9 pt, every label >= 10 pt.
    FS = 10
    plt.rc("text", usetex=True); plt.rc("font", family="serif", size=FS)
    plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)
    rows = {("SO_Z6_R30_T150_C5", "gas1_ics1"): ("data, representative single-gas row", INK, "o"),
            ("SL_Z6_R30_T150_C2", "gas1_ics1"): ("data, Totani's baseline", "#6b46c1", "s"),
            ("a2020", "-"): ("data, A2020", "#2f855a", "^")}
    configs = (("structured", "structured map"), ("flat", "flat mask"), ("construction", "residual construction"))
    W = 0.95 * 510 / 72.27
    fig, axes = plt.subplots(1, 3, figsize=(W, 4.3), sharey=True)
    fig.subplots_adjust(left=0.085, right=0.99, top=0.935, bottom=0.30, wspace=0.07)
    for ax, (c, title) in zip(axes, configs):
        for kind, col, lab in (("null_structured", "#c05621", "no halo, structured bubbles"), ("null_flat", "#b7791f", "no halo, flat bubbles"),
                               ("halo_flat", "#2b6cb0", r"halo injected, flat bubbles")):
            M = np.array([e2(r[c]["counts"]) for r in sims if r["kind"] == kind])
            ax.fill_between(E, np.percentile(M, 16, 0), np.percentile(M, 84, 0), color=col, alpha=0.18, lw=0)
            ax.plot(E, np.median(M, 0), color=col, lw=1.2, label=lab)
        for j, (rk, (lab, col, mk)) in enumerate(rows.items()):
            x = next(r for r in data if tuple(r["row"]) == rk and r["config"] == c)
            y = e2(x["counts"]); lo = e2(x["lo"]); hi = e2(x["hi"])
            xs = E * (1 + 0.04 * (j - 1))
            ax.errorbar(xs, y, yerr=[y - lo, hi - y], fmt=mk, ms=3.5, color=col, lw=0.9, capsize=0, label=lab)
        ax.axhline(0, color=MUTED, lw=0.6); ax.set_xscale("log"); ax.set_title(title, fontsize=FS)
        ax.set_xlabel(r"$E$ [GeV]", fontsize=FS)
    axes[0].set_ylabel(r"$E^2\,dN/dE$ [$10^{-7}$ GeV cm$^{-2}$ s$^{-1}$]", fontsize=FS); axes[0].set_ylim(-4, 12); axes[0].set_yticks([0, 5, 10])
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, fontsize=FS, frameon=False, bbox_to_anchor=(0.54, 0.0),
               handletextpad=0.4, columnspacing=1.5)
    out = ROOT / "results" / f"oscar_spectra_vs_sims_{DATE}{SUFFIX}.pdf"; fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=200); print(out)

if __name__ == "__main__":
    recovery(); spectra()
