#!/usr/bin/env python
"""Oscar's note of 2026-09-29, item 3: the injection test of Totani's
construction (bubble_totani_injection.py, 09-19) redone with the REFINED
amplitude estimator at the ORIGINAL injection strengths, templates fixed and
rebuilt side by side, with the profile curves of the previously discrepant
bins and a stability check of the refinement.

    python scripts/injection_refined.py      # 10-deg cells, one process, ~6 GB

Estimator: likelihood.fit_bin_profile with n_refine = 3 (the 09-19 run,
lattice a_full/768 -- 2,423 counts at 1.5 GeV), 7 (a_full/49,152 -- ~38 counts)
and 8 (stability). Injections: 0, 1, 2 x the fitted 20.8 GeV amplitude,
flat and Totani-like spectra, exactly as on 09-19. Per injection and bin:
recovered amplitude (fixed and rebuilt templates), TS, and min expected count
(signed fits must keep mu > 0). For bins 1.5, 2.6, 4.3 GeV at injection level
1, the full -2lnL profile over the halo amplitude (61 points) for fixed and
rebuilt templates.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from bubble_totani import with_residual_templates  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani.bubble_totani import bubble_image, fit_all_cells, pixel_model_for_bin  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.halo import nfw_rho2_on_target  # noqa: E402
from totani.likelihood import (_fit_bg_given_offset, _mu_bg, _poisson_val_grad_mu, expected_counts,  # noqa: E402
                               fit_all_bins, fit_bin_profile, restrict_stack)
from totani.stack import build_background_stack, cell_template_from_shape, halo_cell_template  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9
K_BUILD = 2
DATE = "2026-09-29"
REFINES = (3, 7, 8)


def slope(inj, rec):
    inj = np.asarray(inj, float); rec = np.asarray(rec, float); d = inj - inj.mean()
    return float(np.sum(d * (rec - rec.mean())) / np.sum(d ** 2))


def main() -> int:
    import resource
    t0 = time.time()
    keep = halo_cell_mask()
    J = halo_cell_template(RUN, EDGES)
    pm = pixel_model_for_bin(RUN, EDGES, K_BUILD, diffuse="galprop_totani", with_gce=True)
    J_pix = nfw_rho2_on_target(pm["target"]) * pm["flat_sr"]
    bg = build_background_stack(RUN, EDGES, bubble="macias", diffuse="galprop_totani")
    roi_J = [float(J[k][keep].sum()) for k in range(N_BINS)]
    print(f"stack ({time.time()-t0:.0f}s)", flush=True)

    def templates_for(counts_pix):
        pm2 = dict(pm); pm2["counts"] = counts_pix
        img = bubble_image(pm2, fit_all_cells(pm2))
        return cell_template_from_shape(RUN, EDGES, img["positive"]), cell_template_from_shape(RUN, EDGES, img["negative"])

    def prepared(counts_cells, pos, neg):
        st = {"counts": counts_cells, "n_bins": N_BINS, "templates": bg["templates"]}
        st, signs = with_residual_templates(st, pos, neg)
        sr, Jr = restrict_stack(st, keep, J)
        names = sorted(sr["templates"])
        return sr, Jr, names

    def fit(sr, Jr, names, n_refine):
        T = np.stack([sr["templates"][n] for n in names])
        v0 = [b["neg2lnl"] for b in fit_all_bins(sr, component_order=names, halo_mask=keep)["per_bin"]]
        out = {"a": [], "ts": [], "min_mu": []}
        for k in range(N_BINS):
            v1, bga, a, _ = fit_bin_profile(T[:, k], sr["counts"][k], Jr[k], n_refine=n_refine)
            mu = expected_counts(np.log10(np.clip(bga, 1e-300, None)), T[:, k]) + a * Jr[k]
            out["a"].append(float(a)); out["ts"].append(float(v0[k] - v1)); out["min_mu"].append(float(mu[keep].min()))
        return out

    def profile_curve(sr, Jr, names, k, a_center, half_width):
        T = np.stack([sr["templates"][n] for n in names])[:, k]
        m = T.shape[0]
        sb = np.maximum(np.array([np.abs(T[i]).sum() for i in range(m)]), 1e-300); ss = max(np.abs(Jr[k]).sum(), 1e-300)
        Tb = T / sb[:, None, None]; Ts = Jr[k] / ss
        seed = np.full(m, np.log10(max(sr["counts"][k].sum() / (m + 1), 1e-300)))
        grid = np.linspace(a_center - half_width, a_center + half_width, 61)
        vals = []
        for a in grid:
            th, _ = _fit_bg_given_offset(Tb, sr["counts"][k], a * ss * Ts, seed)
            vals.append(float(_poisson_val_grad_mu(_mu_bg(th, Tb) + a * ss * Ts, sr["counts"][k])[0]))
        return {"amplitude": grid.tolist(), "counts": (grid * roi_J[k]).tolist(), "neg2lnl": vals}

    pos0, neg0 = templates_for(pm["counts"])
    sr0, Jr0, names0 = prepared(bg["counts"], pos0, neg0)
    base = {n: fit(sr0, Jr0, names0, n) for n in REFINES}
    a_ref = base[7]["a"][5]
    spectra = {"flat": np.ones(N_BINS), "totani_like": np.clip(np.array(base[7]["a"]) / a_ref, 0, None)}
    # the 09-19 run defined the Totani-like shape from its own n_refine=3 fit; record both
    old = json.loads((ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json").read_text())
    spectra_old_shape = np.array(old["spectra"]["totani_like"]["shape"])
    out = {"date": DATE, "a_ref": a_ref, "roi_J": roi_J, "refines": REFINES, "a0": {str(n): base[n] for n in REFINES},
           "totani_like_shape": spectra["totani_like"].tolist(), "totani_like_shape_0919": spectra_old_shape.tolist(), "spectra": {}}
    for sname, shape in (("flat", spectra["flat"]), ("totani_like", spectra_old_shape)):
        rec = {"fixed": {n: [] for n in REFINES}, "rebuilt": {n: [] for n in REFINES}}
        curves = {}
        for lev in (0.0, 1.0, 2.0):
            a_bins = lev * a_ref * shape
            counts_cells = bg["counts"] + a_bins[:, None, None] * J
            pos, neg = templates_for(pm["counts"] + a_bins[K_BUILD] * J_pix) if lev > 0 else (pos0, neg0)
            srF, JrF, nF = prepared(counts_cells, pos0, neg0)
            srR, JrR, nR = prepared(counts_cells, pos, neg)
            for n in REFINES:
                rec["fixed"][n].append(fit(srF, JrF, nF, n)); rec["rebuilt"][n].append(fit(srR, JrR, nR, n))
            if lev == 1.0:
                for k in (0, 1, 2):
                    aF = rec["fixed"][7][-1]["a"][k]; aR = rec["rebuilt"][7][-1]["a"][k]
                    hw = 3 * max(abs(aF - aR), abs(a_bins[k]), 0.02 * abs(a_ref))
                    curves[str(k)] = {"injected": float(a_bins[k]), "fixed": profile_curve(srF, JrF, nF, k, 0.5 * (aF + aR), hw),
                                      "rebuilt": profile_curve(srR, JrR, nR, k, 0.5 * (aF + aR), hw)}
            print(f"[{sname}] level {lev}: rebuilt n7 counts " + " ".join(f"{a*j:8.0f}" for a, j in zip(rec['rebuilt'][7][-1]['a'], roi_J))
                  + f" ({time.time()-t0:.0f}s)", flush=True)
        inj = np.array([[lev * a_ref * shape[k] for k in range(N_BINS)] for lev in (0.0, 1.0, 2.0)])
        sl = {kind: {str(n): [slope(inj[:, k], [r["a"][k] for r in rec[kind][n]]) if inj[:, k].std() > 0 else None for k in range(N_BINS)]
                     for n in REFINES} for kind in ("fixed", "rebuilt")}
        out["spectra"][sname] = {"shape": shape.tolist(), "injected": inj.tolist(),
                                 "recovered": {kind: {str(n): rec[kind][n] for n in REFINES} for kind in rec}, "slopes": sl, "profiles": curves}
        for kind in ("fixed", "rebuilt"):
            for n in REFINES:
                print(f"  [{sname}] {kind:7s} n_refine {n}: " + " ".join("  n/a" if s is None else f"{s:6.2f}" for s in sl[kind][str(n)]), flush=True)
    (ROOT / "results" / "fits" / f"injection_refined_{DATE}.json").write_text(json.dumps(out, indent=1))
    print(f"total {time.time()-t0:.0f}s; peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e9:.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
