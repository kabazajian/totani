#!/usr/bin/env python
"""Quadrant likelihood test on skies with (a) a spherical halo injected on
Totani's spectrum and (b) the same halo injected in the NORTHERN half only,
both on the structured-bubble null truth, fitted with the true row and the
residual construction (rebuilt from the injected pixels) or the structured
map. If the construction returns a common-amplitude fit for the asymmetric
injection, it launders asymmetry and the data's sphericity under it is not
informative; if the asymmetry survives, the data's sphericity counts.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import stage_e
from quadrant_lr_test import fit_joint, quadrant_labels, K20, EDGES
from totani import likelihood as LK, ranking, selection, simsky
from totani.cells import cell_centers, halo_cell_mask
from totani.halo import nfw_rho2_on_target
ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
SETS = {"construct": ("bubble_totani", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
        "structured": ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic")}
sh = ranking.shared_cells(EDGES); J = sh["J"]; keep = halo_cell_mask(); shared = simsky.load_shared()
truth = dict(np.load(stage_e.truth_path("null_structured"))); rowpix = simsky.load_row(*ROW); lab = quadrant_labels(False)
amp = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")
target, L, B, _ = ranking._grid_and_basis(EDGES); Jpix = nfw_rho2_on_target(target) * np.asarray(shared["flat_sr"])
_, Bc = cell_centers(); north_cells = (Bc > 0)[None]; north_pix = B > 0
halo_cells = amp[:, None, None] * J; halo_pix = amp[2] * Jpix
out = {}
for inj in ("spherical", "north_only"):
    hc = halo_cells * (north_cells if inj == "north_only" else 1); hp = halo_pix * (north_pix if inj == "north_only" else 1)
    out[inj] = {s: [] for s in SETS}
    for i in range(20):
        rng = np.random.default_rng([stage_e.TRUTH_SEED["null_structured"], i, 20260922])
        t2 = dict(truth); t2["halo_cells"] = hc; t2["halo_pix2"] = hp
        cells, pix = simsky.draw(t2, sh["counts"], shared, rng)
        shr = {**sh, "counts": cells, **simsky.residual_cells(pix, rowpix, shared, EDGES)}
        for sname, cands in SETS.items():
            st = selection.base_stack(ROW[0], ROW[1], EDGES, shr)
            for c in cands:
                st, _, _ = selection.add_candidate(st, c, ROW[0], ROW[1], EDGES, shr)
            sr, Jr = LK.restrict_stack(st, keep, J); names = sorted(sr["templates"])
            T = np.stack([sr["templates"][n][K20] for n in names]); n = sr["counts"][K20]
            v_c, bg_amp, a_c, _ = LK.fit_bin_profile(T, n, Jr[K20]); v_cl, bgc, ac = fit_joint(T, n, [Jr[K20]], bg_amp, [a_c]); v_c = min(v_c, v_cl)
            pieces = [Jr[K20] * (lab == q) for q in range(4)]; v_q, _, aq = fit_joint(T, n, pieces, bgc, [ac[0]] * 4)
            out[inj][sname].append({"delta": float(v_c - v_q), "common": float(a_c * Jr[K20][keep].sum()),
                                    "quadrants": [float(a * p[keep].sum()) for a, p in zip(aq, pieces)]})
        print(f"[{inj} {i+1}/20] " + "  ".join(f"{s}: Δ={out[inj][s][-1]['delta']:.1f} q={[round(x) for x in out[inj][s][-1]['quadrants']]}" for s in SETS), flush=True)
summ = {}
for inj, d in out.items():
    for s, L_ in d.items():
        dd = np.array([x["delta"] for x in L_]); q = np.array([x["quadrants"] for x in L_])
        summ[f"{inj}__{s}"] = {"delta_median": float(np.median(dd)), "delta_p16": float(np.percentile(dd, 16)), "delta_p84": float(np.percentile(dd, 84)),
                               "common_median": float(np.median([x["common"] for x in L_])), "quadrants_median": np.median(q, 0).tolist()}
        print(f"{inj:10s} {s:10s}", json.dumps(summ[f'{inj}__{s}']))
(ROOT / "results" / "fits" / "quadrant_lr_injected_2026-09-24.json").write_text(json.dumps({"date": "2026-09-24", "summary": summ, "skies": out}, indent=1))
