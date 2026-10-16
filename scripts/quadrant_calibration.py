#!/usr/bin/env python
"""Oscar's note of 2026-09-29, 'The quadrant comparison tests a different
hypothesis': Q_quad = sum_e [-2lnL(common) - (-2lnL(four quadrant amplitudes))]
with the halo partitioned in SKY INTENSITY before PSF folding (equal quadrant
amplitudes reproduce the folded common halo, spillover included), identical
background freedoms refitted in both, and its null distribution CALIBRATED on
simulated skies drawn from the fitted common-halo model -- not assumed chi2_3K.

    FINE_TAG=psf python scripts/quadrant_calibration.py --cell 2 --skies 200 --workers 5

Divisions: standard (l<>0 x b<>0) and rotated 45 deg (N/S/E/W wedges); the
maximum over the two is also calibrated (the rotation choice repeated in every
sky). Statistics: the 20.8 GeV bin alone (chi2_3 reference), the sum over all
nine bins (chi2_27) and over the response-valid bins 5.6-130 GeV (chi2_18).
Configurations: flat mask and structured map with the 2-deg companions, and
the flat mask with the flexible (12-piece brightness) eROSITA lobes. Rows for
the data: Totani's baseline, the simulation truth row, A2020. Null skies: the
truth row SO_Z6_R30_T150_C5 1-gas, background + common halo fitted to the data
per bin (signed), Poisson-sampled on the fitted cells; each sky refitted with
the same row and configuration.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import beta, chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from ns_halo_spectra import joint  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine  # noqa: E402
from totani.likelihood import restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-29"
K20 = 5
RESP = [3, 4, 5, 6, 7, 8]
COMP = ("erosita_unif_N", "erosita_unif_S", "loopI_haslam_raw", "sun_moon_ic")
CONFIGS = ("flat", "structured", "flat_flexlobes")
DATA_ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("a2020", "-")]
TRUE_ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
_W: dict = {}


def quad_path(cell):
    return fine.cache_dir(cell) / "halo_quadrants.npz"


def build_quadrants(cell):
    """NFW rho^2 split in sky intensity (pixel level), each piece exposure-weighted and PSF-folded."""
    p = quad_path(cell)
    if p.exists():
        return p
    from totani.halo import nfw_rho2_on_target
    g = fine._psf_grid(EDGES); L, B = g["L"], g["B"]
    Jp = nfw_rho2_on_target(g["target"])
    div = {"std": {"WN": (L < 0) & (B > 0), "WS": (L < 0) & (B < 0), "EN": (L >= 0) & (B > 0), "ES": (L >= 0) & (B < 0)},
           "rot": {"N": B > np.abs(L), "S": B < -np.abs(L), "E": (L >= 0) & (L >= np.abs(B)), "W": (L < 0) & (L <= -np.abs(B))}}
    # rot: assign the diagonal ties (|l| == |b|) to E/W; the four masks partition the grid
    out = {}
    for d, masks in div.items():
        tot = sum(m.astype(int) for m in masks.values())
        assert np.all(tot == 1), f"{d} masks do not partition the grid"
        for q, m in masks.items():
            out[f"{d}__{q}"] = fine._cells_psf((Jp * m)[None] * g["basis"], g, cell, EDGES)
    lib = fine.load_library(cell)
    for d in div:
        s = sum(out[f"{d}__{q}"] for q in div[d])
        assert np.allclose(s, lib["J"], rtol=1e-6, atol=1e-6 * np.abs(lib["J"]).max()), f"{d}: quadrants do not sum to the folded halo"
    np.savez(p, **out)
    return p


def _init(cell):
    import erosita_flex_ns as EF
    _W["cell"] = cell; _W["lib"] = fine.load_library(cell)
    z = np.load(quad_path(cell)); _W["quads"] = {d: [z[k] for k in z.files if k.startswith(d + "__")] for d in ("std", "rot")}
    if cell == EF.CELL:
        EF._init(); _W["flex"] = EF._W["models"]["bright_flex"]
    else:
        _W["flex"] = None


def stack(counts, row_key, config):
    lib = _W["lib"]; row = fine.load_row(_W["cell"], *row_key)
    st = fine.base_stack(lib, row, counts)
    st["templates"]["bubble_structured" if config == "structured" else "bubble_flat"] = lib["templates"]["bubble_structured" if config == "structured" else "bubble_flat"]
    if config == "flat_flexlobes":
        st["templates"].update(_W["flex"]); st["templates"]["loopI_haslam_raw"] = lib["templates"]["loopI_haslam_raw"]
        st["templates"]["sun_moon_ic"] = lib["templates"]["sun_moon_ic"]
    else:
        for n in COMP:
            st["templates"][n] = lib["templates"][n]
    return st


def qquad(counts, row_key, config):
    """Per-bin Delta(-2lnL) common vs four, both divisions; and the common-halo fit."""
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    st = stack(counts, row_key, config)
    sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    res = {"std": [], "rot": [], "common_amp": [], "v_common": []}
    for k in range(9):
        T = np.stack([sr["templates"][n][k] for n in names]); n = sr["counts"][k]
        v1, a1 = joint(T, n, [Jr[k]])
        res["common_amp"].append(float(a1[0])); res["v_common"].append(v1)
        for d in ("std", "rot"):
            v4, _ = joint(T, n, [q[k] * keep for q in _W["quads"][d]])
            res[d].append(float(v1 - v4))
    return res


def truth_mu(cell, config):
    """Background + common halo fitted to the data (true row), per bin: the calibration null."""
    from totani.likelihood import expected_counts
    _init(cell)
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    st = stack(lib["counts"], TRUE_ROW, config)
    sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    mu = np.zeros_like(lib["counts"], dtype=float)
    from scipy.optimize import minimize
    from totani import likelihood as LK
    for k in range(9):
        T = np.stack([sr["templates"][n][k] for n in names]); nk = sr["counts"][k]
        m = T.shape[0]; sb = np.maximum(np.abs(T).reshape(m, -1).sum(1), 1e-300); ss = max(np.abs(Jr[k]).sum(), 1e-300)
        Tb = (T / sb[:, None, None]).reshape(m, -1); Ts = (Jr[k] / ss).reshape(-1); nf = nk.reshape(-1)
        def obj(x):
            mu_ = x[:m] @ Tb + x[m] * Ts; v, d = LK._poisson_val_grad_mu(mu_, nf)
            return v, np.concatenate([Tb @ d, [Ts @ d]])
        x0 = np.concatenate([np.full(m, nf.sum() / (m + 1)), [0.0]])
        r = minimize(obj, x0, jac=True, method="L-BFGS-B", bounds=[(0, 1e14)] * m + [(-1e14, 1e14)], options={"maxiter": 8000, "ftol": 1e-15, "gtol": 1e-8})
        full = np.stack([st["templates"][n][k] for n in names])
        mu[k] = np.tensordot(r.x[:m] / sb, full, axes=(0, 0)) + (r.x[m] / ss) * J[k]
    return mu


def task(t):
    kind, i, payload = t
    lib = _W["lib"]; keep = lib["keep"]
    if kind == "data":
        row_key, config = payload
        return {"kind": "data", "row": list(row_key), "config": config, **qquad(lib["counts"], row_key, config)}
    config, mu = payload
    rng = np.random.default_rng([13, i, 20260929, CONFIGS.index(config)])
    counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    return {"kind": "null", "sim": i, "config": config, **qquad(counts, TRUE_ROW, config)}


def summaries(r):
    s = {}
    for d in ("std", "rot"):
        x = np.array(r[d])
        s[d] = {"q20": float(x[K20]), "sum9": float(x.sum()), "sum6": float(x[RESP].sum())}
    s["max"] = {k: max(s["std"][k], s["rot"][k]) for k in ("q20", "sum9", "sum6")}
    return s


def emp_p(obs, null):
    null = np.asarray(null); n = len(null); k = int((null >= obs).sum())
    lo = beta.ppf(0.16, k, n - k + 1) if k > 0 else 0.0; hi = beta.ppf(0.84, k + 1, n - k) if k < n else 1.0
    return k, n, (k + 1) / (n + 1), float(lo), float(hi)


def render(res, cell, n_skies):
    ref = {"q20": 3, "sum9": 27, "sum6": 18}
    lines = [f"# Quadrant comparison, calibrated — {cell:g}° cells, PSF-folded ({DATE})", "",
             "Q_quad = −2lnL(common halo) − (−2lnL(four quadrant amplitudes)), halo partitioned in sky intensity before PSF folding, "
             "identical backgrounds refitted. 'q20' = 20.8 GeV bin alone (the 12–19 quoted on 09-26 were this, on 10° cells); "
             "'sum9' = all nine bins; 'sum6' = response-valid bins 5.6–130 GeV. 'max' = the larger of the standard and 45°-rotated divisions, "
             "the maximization repeated in every simulated sky. Null: skies drawn from the fitted common-halo model of the true row "
             f"(SO_Z6_R30_T150_C5, 1-gas), {n_skies} per configuration, refitted with the same row and set. Empirical p = (k+1)/(n+1) "
             "with the 68% Clopper–Pearson range of k/n; χ² p for reference only.", ""]
    for config in CONFIGS:
        N = [summaries(r) for r in res if r["kind"] == "null" and r["config"] == config]
        if not N:
            continue
        lines += [f"## {config}", "", "| statistic | null median [16–84%] | null 95% / 99% | χ² dof, median |", "|---|---|---|---|"]
        for d in ("std", "rot", "max"):
            for k in ("q20", "sum9", "sum6"):
                v = np.array([n[d][k] for n in N])
                lines.append(f"| {d} {k} | {np.median(v):.1f} [{np.percentile(v,16):.1f}, {np.percentile(v,84):.1f}] | {np.percentile(v,95):.1f} / {np.percentile(v,99):.1f} | "
                             f"{ref[k]} ({chi2.median(ref[k]):.1f}) |")
        lines += ["", "| data row | statistic | observed | empirical p [68%] | k/n | χ² p |", "|---|---|---|---|---|---|"]
        for r in [x for x in res if x["kind"] == "data" and x["config"] == config]:
            s = summaries(r)
            for d in ("std", "max"):
                for k in ("q20", "sum9", "sum6"):
                    kk, nn, p, lo, hi = emp_p(s[d][k], [n[d][k] for n in N])
                    lines.append(f"| {r['row'][0]} {r['row'][1]} | {d} {k} | {s[d][k]:.1f} | {p:.3g} [{lo:.3g}, {hi:.3g}] | {kk}/{nn} | {chi2.sf(max(s[d][k], 0), ref[k]):.2g} |")
        lines.append("")
    out = ROOT / "results" / f"quadrant_calibration_cells{cell:g}{fine.tag()}_{DATE}.md"
    out.write_text("\n".join(lines) + "\n"); print("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, default=2.0); ap.add_argument("--skies", type=int, default=200)
    ap.add_argument("--workers", type=int, default=5); ap.add_argument("--render", action="store_true"); a = ap.parse_args()
    assert fine.tag().startswith("_psf"), "run with FINE_TAG=psf or psfuc"
    js = ROOT / "results" / "fits" / f"quadrant_calibration_cells{a.cell:g}{fine.tag()}_{DATE}.json"
    if a.render:
        d = json.loads(js.read_text()); render(d["results"], a.cell, d["skies"]); return 0
    t0 = time.time()
    build_quadrants(a.cell); print(f"quadrant templates ({time.time()-t0:.0f}s)", flush=True)
    configs = CONFIGS if a.cell == 2.0 else ("flat", "structured")
    tasks = [("data", 0, (r, c)) for r in DATA_ROWS for c in configs]
    for c in configs:
        mu = truth_mu(a.cell, c); tasks += [("null", i, (c, mu)) for i in range(a.skies)]
    print(f"{len(tasks)} tasks ({time.time()-t0:.0f}s)", flush=True)
    with Pool(a.workers, initializer=_init, initargs=(a.cell,)) as pool:
        res = []
        for k, r in enumerate(pool.imap_unordered(task, tasks, chunksize=4), 1):
            res.append(r)
            if k % 50 == 0 or k == len(tasks):
                print(f"[{k}/{len(tasks)}] ({time.time()-t0:.0f}s)", flush=True)
    js.write_text(json.dumps({"date": DATE, "cell": a.cell, "skies": a.skies, "results": res}, indent=1))
    render(res, a.cell, a.skies)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
