#!/usr/bin/env python
"""Youyou's overlap-band forecast, RUN rather than described -- all three components.

    python scripts/band_forecast.py

Her test 1 (2026-09-05), verbatim: "the spread across the six models in the IC
prediction restricted to 10 < |b| < 30, per energy bin and specifically near
20 GeV -- both absolute and after renormalizing each model to a common total in
that band. The renormalized number is the one that matters, since a pure level
difference is absorbed by the free per-bin amplitudes and cannot become a
usable direction."

The band is inside the |b| <= 30 maps that already exist, so this runs it --
for inverse Compton as she asked, and for pion decay and bremsstrahlung as Kev
asked (2026-09-06), because the extension is needed for all three and the
hadronic question is not the IC question.

🚨 RENORMALIZING TO A COMMON TOTAL REMOVES ONE DIRECTION. THE FIT ABSORBS 19.
The statistic that forecasts whether the six models will move the HALO is the
leverage of the model-to-model difference against the full background basis
(`totani.degeneracy`), not its raw size. Our own worked example: `ics_ring123`
correlates with the halo at r = 0.94 and dropping it moves the answer < 2%.

🚨 THE SELF-CONSISTENT BASIS CLOSES A CAVEAT THE IC-ONLY RUN LEFT OPEN. The
first run (IC only) found each Galp21 model's IC 2.2-2.8x more halo-degenerate
than Pohl's, but in a HYBRID basis -- Galp21 IC beside Pohl gas -- which can
span more than either family alone. With pion decay and bremsstrahlung read
too, the gas columns can be swapped as well, giving a basis in which every
diffuse component comes from ONE model. That is the fair comparison.

GAS PHASES -> OUR COLUMNS. Pohl's four gas annuli (0-3.5, 3.5-8, 8-10,
10-50 kpc; Readme.txt) are exactly the (1+2+3), 4, 5, 6 collapse of the six
ring groups, the same collapse A2020 uses for IC, so one FOUR_RING serves
both. H2R -> H2_ring1..4. HIR -> HI_ring1..4, WITH HII FOLDED IN: we carry no
ionized-gas template, and the alternative -- adding four columns -- would give
the Galp21 basis more freedom than Pohl's and bias the leverage comparison.
Pion decay and bremsstrahlung are summed per phase, as `galp21_ingest.py` does,
because both trace gas x cosmic rays and that sum is what a "gas emission"
template is.

⚠ THE MAP PATH IN `gcereduce.galp21.BASE` IS STALE. The automounter map moved
to `/private/nfs/Data2` (`/etc/auto_nucosmo`). Overridden here rather than
edited there, because that module belongs to the frozen gce-fisher analysis.

⚠ EXPOSURE IS APPLIED AT CELL LEVEL to the Galp21 columns, via the per-cell
exposure x solid angle our own `isotropic` template carries. Adequate for a
FORECAST statistic on 10 deg cells; not for a fitted amplitude. In the IC-only
run this approximation pushed AGAINST the observed leverage drop, so it cannot
be the cause of it.

THE CELL ARRAYS ARE SAVED. 918 native files, ~90 min over the share, reduce to
a few MB of cell sums; `CELLS` caches them so any re-analysis is free. Delete
it to force a re-read.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import totani  # noqa: E402,F401  -- puts the vendored gcereduce on sys.path
import gcereduce.galp21 as g21  # noqa: E402

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani.cells import bin_to_cells, cell_centers, halo_cell_mask  # noqa: E402
from totani.degeneracy import leverage  # noqa: E402
from totani.stack import build_background_stack, halo_cell_template  # noqa: E402

g21.BASE = Path(os.environ.get("GALP21_NATIVE_MAPS", "/private/nfs/Data2/gce-fisher/GCE_Galprop_maps_2026/maps"))

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
HALF_L, HALF_B = 60.0, 30.0
#: (1+2+3), 4, 5, 6 -- Pohl's four gas annuli AND A2020's four-ring IC collapse.
FOUR_RING = ((1, 2, 3), (4,), (5,), (6,))
KINDS = ("ics_anisotropic", "pion_decay", "bremss")
CELLS = ROOT / "results" / "fits" / "band_forecast_cells_2026-09-06.npz"
OUT = ROOT / "results" / "fits" / "band_forecast_2026-09-06.json"


def native_grid():
    """(L, B) pixel centres of the cutout, degrees, matching `galp21.cutout`."""
    c_left = np.arange(3600 - int(10 * HALF_L), 3600)
    c_right = np.arange(0, int(10 * HALF_L))
    l = np.concatenate([g21.NATIVE_CRVAL1 + (c_left + 1 - g21.NATIVE_CRPIX1) * 0.1 - 360.0,
                        g21.NATIVE_CRVAL1 + (c_right + 1 - g21.NATIVE_CRPIX1) * 0.1])
    rows = g21.cut_rows(HALF_B)
    r = np.arange(rows.start, rows.stop)
    b = 0.0 + (r + 1 - 300.5) * 0.1
    return np.meshgrid(l, b)


def _files(model, kind, phase):
    """Every native file that feeds one (model, kind, phase), by ring group."""
    out = {}
    for group in range(1, 7):
        fs = []
        for ring in g21.SIX_RING_GROUPS[group]:
            if kind == "ics_anisotropic":
                fs += [g21.component_path(model, kind, ring=ring, isrf=i)
                       for i in range(g21.N_ISRF_COMPONENTS)]
            else:
                fs.append(g21.component_path(model, kind, ring=ring, phase=phase))
        out[group] = fs
    return out


def _phases(kind):
    return (None,) if kind == "ics_anisotropic" else g21.GAS_PHASES


def preflight():
    """Fail before the first read if any of the 918 files is missing."""
    missing = [p for m in g21.MODELS for k in KINDS for ph in _phases(k)
               for fs in _files(m, k, ph).values() for p in fs if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} files missing, first: {missing[0]}")
    n = sum(len(fs) for m in g21.MODELS for k in KINDS for ph in _phases(k)
            for fs in _files(m, k, ph).values())
    print(f"preflight: all {n} native files present", flush=True)


def kind_cells(model, kind, phase, L, B, omega):
    """(4, 44, 12, 12): one (model, kind, phase), four-ring collapsed, in cells.

    Solid angle is applied before binning: the cubes are intensities and
    `bin_to_cells` sums, so an unweighted sum would misweight latitude.
    """
    six = np.empty((6, g21.NATIVE_SHAPE[0], 12, 12))
    for group, fs in _files(model, kind, phase).items():
        acc = None
        for p in fs:
            a = g21.cutout(p, half_l=HALF_L, half_b=HALF_B)
            acc = a if acc is None else acc + a
        six[group - 1] = bin_to_cells(acc * omega[None, :, :], L, B)
        del acc
    return np.stack([six[[i - 1 for i in grp]].sum(axis=0) for grp in FOUR_RING])


def read_all():
    if CELLS.exists():
        z = np.load(CELLS)
        print(f"loaded cached cells from {CELLS.name}", flush=True)
        return {k: z[k] for k in z.files}
    preflight()
    t0 = time.time()
    L, B = native_grid()
    omega = np.cos(np.radians(B))
    cells = {}
    for model in g21.MODELS:
        for kind in KINDS:
            for phase in _phases(kind):
                key = f"{model}__{kind}__{phase or 'all'}"
                cells[key] = kind_cells(model, kind, phase, L, B, omega)
                print(f"  {key} ({time.time()-t0:.0f}s)", flush=True)
    np.savez_compressed(CELLS, **cells)
    print(f"saved {CELLS} ({CELLS.stat().st_size/1e6:.1f} MB)", flush=True)
    return cells


def spread_stats(maps, band):
    """maps: list over models of (12, 12) cell maps at one plane."""
    tot = np.array([m[band].sum() for m in maps])
    shapes = np.array([m[band] / m[band].sum() for m in maps])
    return {"absolute_max_over_min": float(tot.max() / tot.min()),
            "renormalized_worst_cell_ratio": float((shapes.max(0) / shapes.min(0)).max()),
            "renormalized_rms_frac_dev": float(
                np.sqrt(((shapes - shapes.mean(0)) ** 2).mean(0)).sum() / shapes.mean(0).sum())}


def main() -> int:
    import resource

    cells = read_all()
    models = list(g21.MODELS)
    energies = g21.energies_mev(models[0])
    k20 = int(np.argmin(np.abs(energies - TOTANI_CENTERS_GEV[5] * 1000.0)))
    band = halo_cell_mask() & (np.abs(cell_centers()[1]) < 30.0)

    def total(model, kind):
        """Sum over phases and the four rings -> (44, 12, 12)."""
        return sum(cells[f"{model}__{kind}__{ph or 'all'}"].sum(axis=0)
                   for ph in _phases(kind))

    results = {"date": "2026-09-06", "models": models,
               "energies_mev": energies.tolist(), "plane_near_20gev": k20,
               "n_band_cells": int(band.sum()), "kinds": {}}

    # ---- test 1 per component, plus the gas sum ---------------------------
    kinds_plus = list(KINDS) + ["gas_pi0_plus_bremss"]

    def total_any(model, kind):
        if kind == "gas_pi0_plus_bremss":
            return total(model, "pion_decay") + total(model, "bremss")
        return total(model, kind)

    for kind in kinds_plus:
        tots = [total_any(m, kind) for m in models]
        r = spread_stats([t[k20] for t in tots], band)
        r["absolute_band_total_20gev"] = [float(t[k20][band].sum()) for t in tots]
        r["per_plane"] = [dict(E_gev=energies[k] / 1000.0,
                               **spread_stats([t[k] for t in tots], band))
                          for k in range(len(energies)) if 1e3 <= energies[k] <= 2e5]
        # F98/R12 pairs: pion decay is ISRF-independent, so they should coincide
        r["isrf_pair_max_frac_diff_20gev"] = float(max(
            np.abs(tots[i][k20][band] - tots[i + 1][k20][band]).sum()
            / tots[i][k20][band].sum() for i in (0, 2, 4)))
        results["kinds"][kind] = r

    # ---- the basis -----------------------------------------------------------
    bg = build_background_stack(RUN, EDGES)
    J = halo_cell_template(RUN, EDGES)
    names = sorted(bg["templates"])
    A = np.stack([bg["templates"][n][5] for n in names])[:, band].T      # (48, 19)
    expo = bg["templates"]["isotropic"][5][band]
    h = J[5][band]
    col = {n: names.index(n) for n in names}

    # (c) what survives the 19-component basis, per component, vs SA0_F98
    for kind in kinds_plus:
        ref = total_any(models[0], kind)[k20][band]
        ref = ref / ref.sum()
        rows = []
        for m in models[1:]:
            s = total_any(m, kind)[k20][band]
            d = (s / s.sum() - ref) * expo
            size = float(np.abs(s / s.sum() - ref).sum())
            rows.append({"model": m, "fractional_size_vs_reference": size,
                         "leverage_vs_19_components":
                         leverage(d, A) if size > 1e-9 else None})
        results["kinds"][kind]["basis_projection"] = rows

    # (d) halo leverage under substitution: IC only, gas only, both
    ic_cols = [col[n] for n in ("ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6")]
    h2_cols = [col[f"H2_ring{i}"] for i in range(1, 5)]
    hi_cols = [col[f"HI_ring{i}"] for i in range(1, 5)]

    def galp_gas(model):
        """(H2 four rings, HI+HII four rings) at plane k20, pi0 + bremss."""
        def ph(kind, phase):
            return cells[f"{model}__{kind}__{phase}"][:, k20]
        h2 = ph("pion_decay", "H2R") + ph("bremss", "H2R")
        hi = (ph("pion_decay", "HIR") + ph("bremss", "HIR")
              + ph("pion_decay", "HII") + ph("bremss", "HII"))
        return h2, hi

    def substituted(model, swap_ic, swap_gas):
        A_m = A.copy()
        if swap_ic:
            ic = cells[f"{model}__ics_anisotropic__all"][:, k20]
            for j, c in enumerate(ic_cols):
                A_m[:, c] = ic[j][band] * expo
        if swap_gas:
            h2, hi = galp_gas(model)
            for j, c in enumerate(h2_cols):
                A_m[:, c] = h2[j][band] * expo
            for j, c in enumerate(hi_cols):
                A_m[:, c] = hi[j][band] * expo
        return leverage(h, A_m)

    results["halo_leverage_pohl_baseline"] = leverage(h, A)
    results["halo_leverage_substituted"] = {
        m: {"ic_only": substituted(m, True, False),
            "gas_only": substituted(m, False, True),
            "self_consistent_ic_and_gas": substituted(m, True, True)}
        for m in models}

    OUT.write_text(json.dumps(results, indent=2))

    print(f"\nplane {k20} = {energies[k20]/1000:.2f} GeV, {band.sum()} band cells")
    print(f"{'component':22s}{'abs max/min':>12s}{'renorm cell':>12s}{'RMS':>8s}{'ISRF-pair':>11s}")
    for kind in kinds_plus:
        r = results["kinds"][kind]
        print(f"{kind:22s}{r['absolute_max_over_min']:12.3f}"
              f"{r['renormalized_worst_cell_ratio']:12.3f}"
              f"{100*r['renormalized_rms_frac_dev']:7.2f}%"
              f"{100*r['isrf_pair_max_frac_diff_20gev']:10.3f}%")
    print("\nleverage of model difference vs the 19-component basis (ref SA0_F98):")
    for kind in kinds_plus:
        vals = [r["leverage_vs_19_components"] for r in results["kinds"][kind]["basis_projection"]]
        print(f"  {kind:22s} " + "  ".join("  --  " if v is None else f"{v:.3f}" for v in vals))
    print(f"\nhalo leverage in the band; Pohl baseline {results['halo_leverage_pohl_baseline']:.4f}")
    print(f"  {'model':12s}{'IC only':>10s}{'gas only':>10s}{'both':>10s}")
    for m, r in results["halo_leverage_substituted"].items():
        print(f"  {m:12s}{r['ic_only']:10.4f}{r['gas_only']:10.4f}"
              f"{r['self_consistent_ic_and_gas']:10.4f}")
    print(f"saved {OUT}")
    print("peak RSS: %.2f GB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
