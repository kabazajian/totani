#!/usr/bin/env python
"""The band-forecast leverage table with Oscar's bubble in the basis.

    python scripts/band_forecast_bubble.py

`band_forecast.py` computed the halo's band leverage under Galp21 substitution
with the Ackermann flat mask as the bubble column, because that was the
background on 2026-09-06. Oscar's structured map arrived that evening and is
itself more halo-like (correlation 0.74 vs 0.52 in the halo cells). This
re-emits the same table for BOTH bubble choices from the two caches --
`band_forecast_cells_2026-09-06.npz` (Galp21 in band cells) and
`bubble_macias_cells_2026-09-07.npz` (our stack, Oscar's template, halo) --
so no stack is rebuilt and the two notes can share one baseline.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from band_forecast import KINDS, _phases  # noqa: E402
from totani.cells import cell_centers, halo_cell_mask  # noqa: E402
from totani.degeneracy import leverage  # noqa: E402

GALP = ROOT / "results" / "fits" / "band_forecast_cells_2026-09-06.npz"
STACK = ROOT / "results" / "fits" / "bubble_macias_cells_2026-09-07.npz"
OUT = ROOT / "results" / "fits" / "band_forecast_bubble_2026-09-07.json"
MODELS = ("SA0_F98", "SA0_R12", "SA50_F98", "SA50_R12", "SA100_F98", "SA100_R12")
K20 = 23   # native plane nearest 20.76 GeV, as in band_forecast.py


def main() -> int:
    g = np.load(GALP)
    z = np.load(STACK)
    names = list(z["names"])
    T, J, mac = z["T"], z["J"], z["macias"]
    keep = halo_cell_mask()
    band = keep & (np.abs(cell_centers()[1]) < 30.0)
    col = {n: names.index(n) for n in names}
    ic_cols = [col[n] for n in ("ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6")]
    h2_cols = [col[f"H2_ring{i}"] for i in range(1, 5)]
    hi_cols = [col[f"HI_ring{i}"] for i in range(1, 5)]
    h = J[5][band]
    expo = T[col["isotropic"], 5][band]

    def cells(model, kind, phase):
        return g[f"{model}__{kind}__{phase or 'all'}"][:, K20]

    def basis(bubble):
        A = np.stack([T[i, 5] for i in range(len(names))])[:, band].T.copy()
        if bubble == "macias":
            A[:, col["bubble"]] = mac[5][band]
        return A

    def substituted(A, model, swap_ic, swap_gas):
        A_m = A.copy()
        if swap_ic:
            ic = cells(model, "ics_anisotropic", None)
            for j, c in enumerate(ic_cols):
                A_m[:, c] = ic[j][band] * expo
        if swap_gas:
            h2 = cells(model, "pion_decay", "H2R") + cells(model, "bremss", "H2R")
            hi = sum(cells(model, k, p) for k in ("pion_decay", "bremss") for p in ("HIR", "HII"))
            for j, c in enumerate(h2_cols):
                A_m[:, c] = h2[j][band] * expo
            for j, c in enumerate(hi_cols):
                A_m[:, c] = hi[j][band] * expo
        return leverage(h, A_m)

    results = {"date": "2026-09-07", "plane": K20, "n_band_cells": int(band.sum()), "bubble": {}}
    for bubble in ("ackermann", "macias"):
        A = basis(bubble)
        r = {"pohl_baseline": leverage(h, A), "per_model": {}}
        for m in MODELS:
            r["per_model"][m] = {"ic_only": substituted(A, m, True, False),
                                 "gas_only": substituted(A, m, False, True),
                                 "self_consistent": substituted(A, m, True, True)}
        results["bubble"][bubble] = r
    OUT.write_text(json.dumps(results, indent=2))

    print(f"{'':14s}{'-- Ackermann flat --':^32s}{'-- Oscar structured --':^32s}")
    print(f"{'basis':14s}{'IC':>10s}{'gas':>10s}{'both':>10s}  {'IC':>10s}{'gas':>10s}{'both':>10s}")
    a, o = results["bubble"]["ackermann"], results["bubble"]["macias"]
    print(f"{'Pohl':14s}{a['pohl_baseline']:10.3f}{a['pohl_baseline']:10.3f}{a['pohl_baseline']:10.3f}  "
          f"{o['pohl_baseline']:10.3f}{o['pohl_baseline']:10.3f}{o['pohl_baseline']:10.3f}")
    for m in MODELS:
        pa, po = a["per_model"][m], o["per_model"][m]
        print(f"{m:14s}{pa['ic_only']:10.3f}{pa['gas_only']:10.3f}{pa['self_consistent']:10.3f}  "
              f"{po['ic_only']:10.3f}{po['gas_only']:10.3f}{po['self_consistent']:10.3f}")
    sc_a = [v["self_consistent"] for v in a["per_model"].values()]
    sc_o = [v["self_consistent"] for v in o["per_model"].values()]
    print(f"\nself-consistent range: flat {min(sc_a):.3f}-{max(sc_a):.3f} ({max(sc_a)/min(sc_a):.2f}x); "
          f"Oscar {min(sc_o):.3f}-{max(sc_o):.3f} ({max(sc_o)/min(sc_o):.2f}x)")
    print(f"Pohl/self-consistent: flat {a['pohl_baseline']/max(sc_a):.1f}-{a['pohl_baseline']/min(sc_a):.1f}x; "
          f"Oscar {o['pohl_baseline']/max(sc_o):.1f}-{o['pohl_baseline']/min(sc_o):.1f}x")
    print(f"saved {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
