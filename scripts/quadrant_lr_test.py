#!/usr/bin/env python
"""Oscar's item 4b: replace 'four positive quadrants' by a likelihood comparison
between ONE halo normalization and FOUR independent quadrant normalizations,
with the same background freedoms in both; and repeat with the quadrant
division rotated by 45 degrees.

    python scripts/quadrant_lr_test.py

H_common: background + one signed halo amplitude (per bin).
H_quad:   background + four signed amplitudes, one per quadrant of the halo
          template (3 extra parameters per bin). Delta(-2lnL) ~ chi^2_3 if the
          halo really is one spherical component.
Both are fitted jointly in linear amplitude space (background >= 0, halo
amplitudes signed), the convex problem, from the common-halo solution.
Representative rows x the three Fig.-1 panel sets, all nine bins; the 20.8
GeV bin is reported and the sum over bins.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking, selection  # noqa: E402
from totani import likelihood as LK  # noqa: E402
from totani.cells import cell_centers, halo_cell_mask  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-24"
K20 = 5
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SL_Z10_R30_T150_C5", "gas1_ics1"),
        ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-"),
        ("G21_SA0_R12+iso", "gas1_ics1"), ("G21_SA100_F98+iso", "gas5_ics3")]
PANELS = {"residual": ("bubble_totani", "loopI_shells", "erosita_unif_NS"),
          "flat": ("bubble_flat", "loopI_shells", "erosita_unif_NS"),
          "structured": ("bubble_structured", "loopI_shells", "erosita_unif_NS")}


def quadrant_labels(rotated: bool):
    Lc, Bc = cell_centers()
    if not rotated:
        return (Lc > 0).astype(int) + 2 * (Bc > 0).astype(int)          # E/W x N/S
    # rotated by 45 deg: |l| > |b| -> east/west wedges, |l| < |b| -> north/south wedges
    ew = np.abs(Lc) >= np.abs(Bc)
    return np.where(ew, (Lc > 0).astype(int), 2 + (Bc > 0).astype(int))


def fit_joint(T_bg, n, T_sig_list, x0_bg, a0_sig):
    """Linear-space L-BFGS-B: background >= 0, signals signed. Returns -2lnL and amplitudes."""
    m = T_bg.shape[0]; s = len(T_sig_list)
    scale_bg = np.maximum(np.array([np.abs(T_bg[i]).sum() for i in range(m)]), 1e-300)
    scale_sg = np.maximum(np.array([np.abs(t).sum() for t in T_sig_list]), 1e-300)
    Tb = (T_bg / scale_bg[:, None, None]).reshape(m, -1); Ts = np.stack([t / sc for t, sc in zip(T_sig_list, scale_sg)]).reshape(s, -1)
    nf = n.reshape(-1)
    def obj(x):
        mu = x[:m] @ Tb + x[m:] @ Ts
        val, dmu = LK._poisson_val_grad_mu(mu, nf)
        return val, np.concatenate([Tb @ dmu, Ts @ dmu])
    x0 = np.concatenate([np.asarray(x0_bg) * scale_bg, np.asarray(a0_sig) * scale_sg])
    bounds = [(0.0, 1e14)] * m + [(-1e14, 1e14)] * s
    r = minimize(obj, x0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 6000, "ftol": 1e-15, "gtol": 1e-8})
    for _ in range(6):
        if r.status == 0:
            break
        r = minimize(obj, r.x, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 6000, "ftol": 1e-15, "gtol": 1e-8})
    return r.fun, r.x[:m] / scale_bg, r.x[m:] / scale_sg


def main() -> int:
    from scipy.stats import chi2
    sh = ranking.shared_cells(EDGES); J = sh["J"]; keep = halo_cell_mask()
    out = {"date": DATE, "rows": []}
    lines = [f"# Quadrant likelihood test: one halo normalization vs four ({DATE})", "",
             "Δ(−2lnL) = −2lnL(common) − −2lnL(four quadrant amplitudes), 3 extra parameters per bin; p from χ²₃ at 20.8 GeV. "
             "Quadrant amplitudes in counts at 20.8 GeV (E−S, E+ S, E− N, E+ N for the standard division; W, E, S, N wedges for the rotated one).", "",
             "| row | panel | division | Δ(−2lnL) @20.8 (p) | Δ summed 9 bins | common halo@20.8 | quadrant amplitudes @20.8 | all four positive |",
             "|---|---|---|---|---|---|---|---|"]
    for row in ROWS:
        for pname, cands in PANELS.items():
            try:
                st = selection.base_stack(row[0], row[1], EDGES, sh)
                for c in cands:
                    st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
            except (FileNotFoundError, KeyError):
                continue
            sr, Jr = LK.restrict_stack(st, keep, J); names = sorted(sr["templates"])
            for rotated in (False, True):
                lab = quadrant_labels(rotated)
                rec = {"model": row[0], "variant": row[1], "panel": pname, "rotated": rotated, "bins": []}
                for k in range(9):
                    T = np.stack([sr["templates"][n][k] for n in names]); n = sr["counts"][k]
                    v_common, bg_amp, a_common, _ = LK.fit_bin_profile(T, n, Jr[k])
                    v_common_lin, bgc, ac = fit_joint(T, n, [Jr[k]], bg_amp, [a_common])
                    v_common = min(v_common, v_common_lin)
                    pieces = [Jr[k] * (lab == q) for q in range(4)]
                    v_quad, bgq, aq = fit_joint(T, n, pieces, bgc, [ac[0]] * 4)
                    roi = [float(p[keep].sum()) for p in pieces]
                    rec["bins"].append({"delta": float(v_common - v_quad), "common": float(a_common * Jr[k][keep].sum()),
                                        "quadrants": [float(a * r) for a, r in zip(aq, roi)]})
                d20 = rec["bins"][K20]["delta"]; dsum = sum(b["delta"] for b in rec["bins"])
                q20 = rec["bins"][K20]["quadrants"]
                rec.update({"delta20": d20, "p20": float(chi2.sf(max(d20, 0), 3)), "delta_sum": dsum})
                out["rows"].append(rec)
                lines.append(f"| {row[0]} {row[1]} | {pname} | {'rotated 45°' if rotated else 'standard'} | {d20:.1f} ({rec['p20']:.2g}) | {dsum:.0f} | "
                             f"{rec['bins'][K20]['common']:,.0f} | " + ", ".join(f"{x:,.0f}" for x in q20) + f" | {'Y' if all(x > 0 for x in q20) else 'n'} |")
                print(lines[-1], flush=True)
    (ROOT / "results" / "fits" / f"quadrant_lr_test_{DATE}.json").write_text(json.dumps(out, indent=1))
    (ROOT / "results" / f"quadrant_lr_test_{DATE}.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
