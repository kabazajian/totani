#!/usr/bin/env python
"""Oscar's note of 2026-09-29, item 2: residual maps at finer resolution --
where the halo-like counts sit under each bubble description.

    FINE_TAG=psf python scripts/residual_maps.py --cell 2

Totani's baseline row (single-gas), each bubble description with the 2-deg
companions. Top row: (data - background-only model) / sqrt(model), summed over
12.3-35 GeV (bins 4-6) with each bin's own background fit. Bottom row: the
same after the halo is added (H1). Cells outside the fit (|b| < 10, masked
sources) in gray; the flat-bubble footprint outlined.
"""
from __future__ import annotations

import argparse
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
from totani.likelihood import expected_counts, fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

DATE = "2026-09-29"
ROW = ("SL_Z6_R30_T150_C2", "gas1_ics1")
BINS = (4, 5, 6)
AGG = 2                     # 2 x 2 cells: 4-deg blocks at 2 deg
TITLES = {"structured": "structured map", "flat": "flat mask", "construction": "residual construction"}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, default=2.0); a = ap.parse_args()
    OC._init(a.cell); lib = OC._W["lib"]; keep, J = lib["keep"], lib["J"]
    Lc, Bc = cell_centers(cell_deg=a.cell)
    plain = np.load(ranking.CACHE / f"cells{a.cell:g}" / "library.npz"); foot = np.abs(plain["T__bubble_flat"]).sum(0) > 0
    maps = {}
    for c in OC.CONFIGS:
        st = OC.stack_for(lib["counts"], ROW, c)
        sr, Jr = restrict_stack(st, keep, J)
        bg = fit_all_bins(sr, halo_mask=keep); pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False, n_refine=OC.N_REFINE)
        names = bg["component_order"]; T = np.stack([st["templates"][n] for n in names])
        N = np.zeros(keep.shape); M0 = np.zeros(keep.shape); M1 = np.zeros(keep.shape)
        for k in BINS:
            b0 = np.array([bg["per_bin"][k]["amplitudes"][n] for n in names])
            mu0 = np.tensordot(b0, T[:, k], axes=(0, 0))
            b1 = np.array([pf["per_bin"][k]["bg_amplitudes"][n] for n in pf["component_order"]])
            T1 = np.stack([st["templates"][n] for n in pf["component_order"]])
            mu1 = np.tensordot(b1, T1[:, k], axes=(0, 0)) + pf["per_bin"][k]["signal_amplitude"] * J[k]
            N += np.where(keep, lib["counts"][k], 0); M0 += np.where(keep, mu0, 0); M1 += np.where(keep, mu1, 0)
        # aggregate AGG x AGG cells (fitted cells only), significance of the summed residual
        agg = lambda x: x.reshape(x.shape[0] // AGG, AGG, x.shape[1] // AGG, AGG).sum((1, 3))
        nk = agg(keep.astype(float))
        sig = lambda M: np.where(nk > 0, (agg(N) - agg(M)) / np.sqrt(np.maximum(agg(M), 1)), np.nan)
        maps[c] = (sig(M0), sig(M1))
    plt.rcParams.update({"font.size": 8, "font.family": "serif"})
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 5.0), sharex=True, sharey=True)
    ext = [Lc.max() + a.cell / 2, Lc.min() - a.cell / 2, Bc.min() - a.cell / 2, Bc.max() + a.cell / 2]
    vmax = 4.0
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("0.85")
    for j, c in enumerate(OC.CONFIGS):
        for i in range(2):
            ax = axes[i, j]
            im = ax.imshow(maps[c][i], origin="lower", extent=ext, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation="nearest", aspect="equal")
            ax.contour(Lc, Bc, foot.astype(float), levels=[0.5], colors="k", linewidths=0.6)
            if i == 0:
                ax.set_title(TITLES[c], fontsize=8)
        axes[0, 0].set_ylabel("no halo\n$b$ [deg]"); axes[1, 0].set_ylabel("with halo\n$b$ [deg]")
    for ax in axes[1]:
        ax.set_xlabel(r"$\ell$ [deg]")
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02); cb.set_label(rf"residual, 12–35 GeV [$\sigma$ per {AGG * a.cell:g}° block]")
    fig.suptitle(f"Totani's baseline diffuse model, {a.cell:g}° cells{', PSF-folded' if fine.tag() else ''}", fontsize=8)
    out = ROOT / "results" / f"residual_maps_cells{a.cell:g}{fine.tag()}_{DATE}.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.with_suffix(".png"), dpi=150, bbox_inches="tight")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
