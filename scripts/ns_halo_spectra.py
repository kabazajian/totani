#!/usr/bin/env python
"""Item 2 (Kev, 2026-09-27): fit the northern and southern halo amplitudes
separately and compare their spectra, on PSF-folded 2-deg cells.

    FINE_TAG=psf python scripts/ns_halo_spectra.py --workers 6

The NFW-rho^2 template is split at b = 0 into halo_N and halo_S, each with a
signed amplitude per bin, fitted jointly with the background (linear space,
background >= 0) under the flat mask and the 2-deg companions -- and again
with independent N/S flat-bubble normalizations. Per bin: the two amplitudes
with 68% intervals (profile over each, the other refitted); exposure-corrected
E^2 dN/dE per hemisphere (J-weighted conversion computed per hemisphere); and
a shape test: Delta(-2lnL) between S = r x N in every bin (one ratio r) and
free S, N per bin (8 extra parameters).

Data: the 33-row subset (16 models x 2 splits + A2020) and Totani's baseline.
Sims (true row): flat-bubble skies with a spherical halo injected on Totani's
spectrum (the halves must agree), and halo-free (both ~0). 20 skies each.
If both halves are peaked with a common shape, the excess is an asymmetric
halo-like signal; if the southern excess has a different spectrum, it is a
southern foreground.
"""
from __future__ import annotations

import argparse
import json
import sys
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from stage_e_2deg import COMP, TRUE_ROW, stack_for, truth_cells  # noqa: E402
from totani import fine  # noqa: E402
from totani import likelihood as LK  # noqa: E402
from totani.cells import cell_centers  # noqa: E402
from totani.likelihood import restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
CELL = 2.0
DATE = "2026-09-27"
_W: dict = {}


def _init():
    _W["lib"] = fine.load_library(CELL)
    _, Bc = cell_centers(cell_deg=CELL); _W["north"] = (Bc > 0)


def joint(T_bg, n, sigs, fixed=None):
    """-2lnL with background >= 0 and signed signals; `fixed` = {index: value} holds signals fixed."""
    m = T_bg.shape[0]; s = len(sigs)
    sb = np.maximum(np.abs(T_bg).reshape(m, -1).sum(1), 1e-300); ss = np.maximum(np.array([np.abs(t).sum() for t in sigs]), 1e-300)
    Tb = (T_bg / sb[:, None, None]).reshape(m, -1); Ts = np.stack([t / c for t, c in zip(sigs, ss)]).reshape(s, -1); nf = n.reshape(-1)
    def obj(x):
        mu = x[:m] @ Tb + x[m:] @ Ts; v, dmu = LK._poisson_val_grad_mu(mu, nf)
        return v, np.concatenate([Tb @ dmu, Ts @ dmu])
    x0 = np.concatenate([np.full(m, nf.sum() / (m + s)), np.zeros(s)])
    bounds = [(0.0, 1e14)] * m + [(-1e14, 1e14)] * s
    if fixed:
        for i, val in fixed.items():
            x0[m + i] = val * ss[i]; bounds[m + i] = (val * ss[i], val * ss[i])
    r = minimize(obj, x0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 8000, "ftol": 1e-15, "gtol": 1e-8})
    for _ in range(6):
        if r.status == 0:
            break
        r = minimize(obj, r.x, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 8000, "ftol": 1e-15, "gtol": 1e-8})
    return float(r.fun), r.x[m:] / ss


def analyse(counts, row_key, split_bubbles, with_errors):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]; north = _W["north"]
    st = stack_for(lib, fine.load_row(CELL, *row_key), "bubble_flat", counts)
    if split_bubbles:
        fb = st["templates"].pop("bubble_flat"); st["templates"]["bubble_flat_N"] = fb * north[None]; st["templates"]["bubble_flat_S"] = fb * ~north[None]
    sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    out = {"N": [], "S": [], "N_err": [], "S_err": [], "v_free": [], "counts_N_per_amp": [], "counts_S_per_amp": []}
    for k in range(9):
        T = np.stack([sr["templates"][nm][k] for nm in names]); n = sr["counts"][k]
        hN = Jr[k] * north; hS = Jr[k] * ~north
        v, a = joint(T, n, [hN, hS])
        out["N"].append(float(a[0] * hN[keep].sum())); out["S"].append(float(a[1] * hS[keep].sum())); out["v_free"].append(v)
        out["counts_N_per_amp"].append(float(hN[keep].sum())); out["counts_S_per_amp"].append(float(hS[keep].sum()))
        if with_errors:
            errs = []
            for i, h in ((0, hN), (1, hS)):
                d = []
                for f in (0.8, 1.2):
                    step = max(abs(a[i]) * (f - 1), 0.05 * max(abs(a[0]), abs(a[1]), 1e-40))
                    vv, _ = joint(T, n, [hN, hS], fixed={i: a[i] + (step if f > 1 else -step)}); d.append((step, vv - v))
                curv = np.mean([dv / st_ ** 2 for st_, dv in d])
                errs.append(float(np.sqrt(1.0 / max(curv, 1e-300)) * h[keep].sum()))
            out["N_err"].append(errs[0]); out["S_err"].append(errs[1])
    # shape test: S = r N in every bin (one ratio) vs free
    def v_ratio(r):
        tot = 0.0
        for k in range(9):
            T = np.stack([sr["templates"][nm][k] for nm in names]); n = sr["counts"][k]
            h = Jr[k] * north + r * Jr[k] * ~north
            v, _ = joint(T, n, [h]); tot += v
        return tot
    rs = np.linspace(0.8, 2.2, 15); vr = [v_ratio(r) for r in rs]
    out["shape_delta"] = float(min(vr) - sum(out["v_free"])); out["best_ratio"] = float(rs[int(np.argmin(vr))])
    return out


def task(t):
    kind, i, mu, row_key, split = t
    lib = _W["lib"]; keep = lib["keep"]
    if mu is None:
        counts = lib["counts"]
    else:
        rng = np.random.default_rng([zlib.crc32(kind.encode()), i, 20260927])
        counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    return {"kind": kind, "sim": i, "row": list(row_key), "split": split, **analyse(counts, row_key, split, with_errors=(kind == "data"))}


def hemi_conversion():
    """Per-hemisphere J-weighted counts -> dN/dE conversion (as halo_spectra_compare, per hemisphere)."""
    from totani import ranking, simsky
    from totani.halo import nfw_rho2_on_target
    from totani.templates import pixel_solid_angle
    basis, L, B, mask = simsky._shape_basis(EDGES)
    target, _, _, _ = ranking._grid_and_basis(EDGES); Jp = nfw_rho2_on_target(target); om = pixel_solid_angle(target)
    roi = (np.abs(B) >= 10) & (np.abs(B) <= 60) & (np.abs(L) <= 60) & ~mask
    out = {}
    for h, m in (("N", roi & (B > 0)), ("S", roi & (B < 0))):
        out[h] = np.array([float((Jp * om)[m].sum()) / float((Jp * basis[k])[m].sum()) for k in range(9)])
    return out


def main() -> int:
    from scipy.stats import chi2
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    data_rows = [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()] + [("a2020", "-"), ("SL_Z6_R30_T150_C2", "gas1_ics1")]
    tasks = [("data", 0, None, r, s) for r in data_rows for s in (False, True)]
    for h in (False, True):
        mu = truth_cells("bubble_flat", h)
        tasks += [(f"{'halo' if h else 'null'}_flat", i, mu, TRUE_ROW, False) for i in range(20)]
    with Pool(a.workers, initializer=_init) as pool:
        res = pool.map(task, tasks)
    conv = hemi_conversion(); E = np.array([1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0])
    e2 = lambda c, h: np.asarray(c) * conv[h] * (E * 1e3) ** 2 / 1e3 * 1e7          # x1e-7 GeV cm^-2 s^-1
    out = {"date": DATE, "psf": fine.tag(), "results": res, "conversion": {k: v.tolist() for k, v in conv.items()}}
    (ROOT / "results" / "fits" / f"ns_halo_spectra{fine.tag()}_{DATE}.json").write_text(json.dumps(out, indent=1))
    lines = [f"# Northern and southern halo amplitudes and spectra, flat mask, 2° cells{' (PSF-folded)' if fine.tag() else ''} ({DATE})", "",
             f"Companions: {', '.join(COMP)}. E²dN/dE per hemisphere (×10⁻⁷ GeV cm⁻² s⁻¹, J-weighted per hemisphere). "
             "Shape test: Δ(−2lnL) of S = r·N in every bin against free N, S (8 extra parameters; χ²₈ median 7.3).", ""]
    for split in (False, True):
        D = [r for r in res if r["kind"] == "data" and r["split"] == split]
        lines += [f"## Data, {'N/S flat-bubble normalizations' if split else 'one flat-bubble normalization'}", "",
                  "| row | shape Δ (p) | best S/N ratio | S/N at 20.8 | N E²dN/dE by bin | S E²dN/dE by bin |", "|---|---|---|---|---|---|"]
        for r in D:
            lines.append(f"| {r['row'][0]} {r['row'][1]} | {r['shape_delta']:.1f} ({chi2.sf(max(r['shape_delta'],0),8):.2g}) | {r['best_ratio']:.2f} | "
                         f"{r['S'][5]/max(r['N'][5],1e-9):.2f} | " + " ".join(f"{x:+.2f}" for x in e2(r['N'], 'N')) + " | " + " ".join(f"{x:+.2f}" for x in e2(r['S'], 'S')) + " |")
        sd = np.array([r["shape_delta"] for r in D]); rat = np.array([r["S"][5] / max(r["N"][5], 1e-9) for r in D])
        lines += ["", f"Median over rows: shape Δ {np.median(sd):.1f}, S/N at 20.8 GeV {np.median(rat):.2f}; median N spectrum "
                  + " ".join(f"{x:+.2f}" for x in e2(np.median([r['N'] for r in D], 0), 'N')) + "; median S spectrum "
                  + " ".join(f"{x:+.2f}" for x in e2(np.median([r['S'] for r in D], 0), 'S')), ""]
    lines += ["## Simulated flat-bubble skies (true row)", "", "| sky | shape Δ median [16–84%] | S/N at 20.8 median | N / S E²dN/dE medians by bin |", "|---|---|---|---|"]
    for kind in ("halo_flat", "null_flat"):
        S = [r for r in res if r["kind"] == kind]
        sd = np.array([r["shape_delta"] for r in S]); rat = np.median([r["S"][5] / max(r["N"][5], 1e-9) for r in S]) if kind == "halo_flat" else float("nan")
        lines.append(f"| {kind} | {np.median(sd):.1f} [{np.percentile(sd,16):.1f}, {np.percentile(sd,84):.1f}] | {rat:.2f} | "
                     + " ".join(f"{x:+.2f}" for x in e2(np.median([r['N'] for r in S], 0), 'N')) + " / " + " ".join(f"{x:+.2f}" for x in e2(np.median([r['S'] for r in S], 0), 'S')) + " |")
    (ROOT / "results" / f"ns_halo_spectra{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
