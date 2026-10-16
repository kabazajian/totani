#!/usr/bin/env python
"""Component-plus-residual maps for every large component (Kev, 2026-10-03;
replaces the residual-map figure Oscar asked for on 2026-09-29).

    FINE_TAG=psfuc python scripts/component_residual_maps.py

Totani's baseline row (single-gas) on 2-deg PSF-folded cells, each bubble
description with the 2-deg companions (eROSITA uniform N/S lobes, Haslam
Loop I, solar IC), the NFW halo fitted (signed amplitude per bin). Rows: the
three bubble descriptions. Columns: the residual of the full model, then for
the halo, the bubbles, Loop I and the eROSITA lobes the component plus the
residual, i.e. data minus every OTHER fitted component -- the emission the
template is fitted to, in the template's place. Each map is summed over
12.3-35 GeV (bins 4-6, each bin's own fit) in 4-deg blocks and shown as
(data - other components) / sqrt(full model), in sigma per block, so the
columns share one scale. Cells outside the fit in gray; the flat-bubble
footprint outlined.
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
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import oscar_conditional as OC  # noqa: E402
from totani import fine, ranking  # noqa: E402
from totani.cells import cell_centers  # noqa: E402
from totani.likelihood import fit_all_bins_profile, restrict_stack  # noqa: E402

DATE = "2026-10-03"
CELL = 2.0
ROW = ("SL_Z6_R30_T150_C2", "gas1_ics1")
BINS = (4, 5, 6)
AGG = 2
ROWS = {"structured": "structured map", "flat": "flat mask", "construction": "residual\nconstruction"}
COLS = [("residual", "residual"), ("halo", "halo"), ("bubbles", "bubbles"), ("loopI", r"Loop\,I"), ("erosita", "eROSITA")]


def group(name):
    if name.startswith("bubble_"):
        return "bubbles"
    if name.startswith("loopI"):
        return "loopI"
    if name.startswith("erosita"):
        return "erosita"
    return None


def main() -> int:
    OC._init(CELL); lib = OC._W["lib"]; keep, J = lib["keep"], lib["J"]
    Lc, Bc = cell_centers(cell_deg=CELL)
    plain = np.load(ranking.CACHE / f"cells{CELL:g}" / "library.npz"); foot = np.abs(plain["T__bubble_flat"]).sum(0) > 0
    agg = lambda x: x.reshape(x.shape[0] // AGG, AGG, x.shape[1] // AGG, AGG).sum((1, 3))
    nk = agg(keep.astype(float))
    maps, summary = {}, {}
    for c in OC.CONFIGS:
        st = OC.stack_for(lib["counts"], ROW, c)
        sr, Jr = restrict_stack(st, keep, J)
        pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False, n_refine=OC.N_REFINE)
        names = pf["component_order"]
        N = np.zeros(keep.shape); M = np.zeros(keep.shape); comp = {g: np.zeros(keep.shape) for g, _ in COLS[1:]}
        for k in BINS:
            amps = pf["per_bin"][k]["bg_amplitudes"]
            for n in names:
                m = amps[n] * st["templates"][n][k]
                M += np.where(keep, m, 0)
                g = group(n)
                if g:
                    comp[g] += np.where(keep, m, 0)
            h = pf["per_bin"][k]["signal_amplitude"] * J[k]
            M += np.where(keep, h, 0); comp["halo"] += np.where(keep, h, 0)
            N += np.where(keep, lib["counts"][k], 0)
        sd = np.sqrt(np.maximum(agg(M), 1))
        res = agg(N) - agg(M)
        maps[c] = {"residual": np.where(nk > 0, res / sd, np.nan)}
        for g, _ in COLS[1:]:
            maps[c][g] = np.where(nk > 0, (res + agg(comp[g])) / sd, np.nan)
        summary[c] = {g: float(np.nansum(comp[g])) for g, _ in COLS[1:]}
        summary[c]["fitted_names"] = names
        print(c, {g: round(v) for g, v in summary[c].items() if g != "fitted_names"},
              {g: (round(float(np.nanmin(m)), 1), round(float(np.nanmax(m)), 1)) for g, m in maps[c].items()}, flush=True)

    # Drawn at its printed size: figure* at \textwidth = 510 pt (revtex4-2 prd),
    # captions 9 pt, every label >= 10 pt.
    from matplotlib.colors import TwoSlopeNorm
    FS = 10
    plt.rc("text", usetex=True); plt.rc("font", family="serif", size=FS)
    plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)
    W = 510 / 72.27
    # square panels: 5 across the width, so the height follows from the width
    left, right, wsp = 0.135, 0.875, 0.04
    pw = (right - left) * W / (5 + 4 * wsp)                 # panel width, inches
    H = 3 * pw * (1 + 0.06) + 0.95
    fig, axes = plt.subplots(3, 5, figsize=(W, H))
    fig.subplots_adjust(left=left, right=right, top=1 - 0.30 / H, bottom=0.62 / H, wspace=wsp, hspace=0.06)
    ext = [Lc.max() + CELL / 2, Lc.min() - CELL / 2, Bc.min() - CELL / 2, Bc.max() + CELL / 2]
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("0.85")
    norm = TwoSlopeNorm(vmin=-4, vcenter=0, vmax=VMAX)
    for i, c in enumerate(OC.CONFIGS):
        for j, (g, title) in enumerate(COLS):
            ax = axes[i, j]
            im = ax.imshow(maps[c][g], origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="nearest", aspect="equal")
            ax.contour(Lc, Bc, foot.astype(float), levels=[0.5], colors="k", linewidths=0.5)
            ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([40, 0, -40]); ax.set_yticks([-40, 0, 40])
            if j > 0:
                ax.set_yticklabels([])
            if i < 2:
                ax.set_xticklabels([])
            if i == 0:
                ax.set_title(title, fontsize=FS)
        axes[i, 0].set_ylabel(r"$b$ [deg]", fontsize=FS, labelpad=1)
        pos = axes[i, 0].get_position()
        fig.text(0.012, 0.5 * (pos.y0 + pos.y1), ROWS[c], rotation=90, ha="left", va="center", fontsize=FS, multialignment="center")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\ell$ [deg]", fontsize=FS)
    p0, p1 = axes[-1, -1].get_position(), axes[0, -1].get_position()
    cax = fig.add_axes([0.89, p0.y0, 0.017, p1.y1 - p0.y0])
    cb = fig.colorbar(im, cax=cax); cb.set_label(r"12--35\,GeV, $\sigma$ per $4^\circ$ block", fontsize=FS)
    out = ROOT / "results" / f"component_residual_maps_cells2{fine.tag()}_{DATE}.pdf"
    fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=200)
    (ROOT / "results" / "fits" / f"component_residual_maps_cells2{fine.tag()}_{DATE}.json").write_text(json.dumps(summary, indent=1))
    print(out)
    return 0


VMAX = 8.0

if __name__ == "__main__":
    raise SystemExit(main())
