#!/usr/bin/env python
"""The quadrant likelihood test (quadrant_lr_test.py) on the halo-free
structured-bubble skies: does the residual construction make a FALSE
component look spherical too? True row, own companions, bubble template
replaced by the construction or the flat mask; 20 skies; standard division,
20.8 GeV bin.
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
from totani.cells import halo_cell_mask
ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
SETS = {"construct": ("bubble_totani", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
        "flat_truth": ("bubble_flat", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
        "truth": ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic")}
sh = ranking.shared_cells(EDGES); J = sh["J"]; keep = halo_cell_mask(); shared = simsky.load_shared()
truth = dict(np.load(stage_e.truth_path("null_structured"))); rowpix = simsky.load_row(*ROW); lab = quadrant_labels(False)
out = {s: [] for s in SETS}
for i in range(20):
    rng = np.random.default_rng([stage_e.TRUTH_SEED["null_structured"], i, 20260922])
    cells, pix = simsky.draw(truth, sh["counts"], shared, rng)
    shr = {**sh, "counts": cells, **simsky.residual_cells(pix, rowpix, shared, EDGES)}
    for sname, cands in SETS.items():
        st = selection.base_stack(ROW[0], ROW[1], EDGES, shr)
        for c in cands:
            st, _, _ = selection.add_candidate(st, c, ROW[0], ROW[1], EDGES, shr)
        sr, Jr = LK.restrict_stack(st, keep, J); names = sorted(sr["templates"])
        T = np.stack([sr["templates"][n][K20] for n in names]); n = sr["counts"][K20]
        v_c, bg_amp, a_c, _ = LK.fit_bin_profile(T, n, Jr[K20]); v_cl, bgc, ac = fit_joint(T, n, [Jr[K20]], bg_amp, [a_c]); v_c = min(v_c, v_cl)
        pieces = [Jr[K20] * (lab == q) for q in range(4)]; v_q, _, aq = fit_joint(T, n, pieces, bgc, [ac[0]] * 4)
        out[sname].append({"delta": float(v_c - v_q), "common": float(a_c * Jr[K20][keep].sum()),
                           "quadrants": [float(a * p[keep].sum()) for a, p in zip(aq, pieces)]})
    print(f"[{i+1}/20] " + "  ".join(f"{s}: Δ={out[s][-1]['delta']:.1f}" for s in SETS), flush=True)
from scipy.stats import chi2
summ = {}
for s, L in out.items():
    d = np.array([x["delta"] for x in L]); summ[s] = {"delta_median": float(np.median(d)), "delta_p16": float(np.percentile(d, 16)), "delta_p84": float(np.percentile(d, 84)),
                                                      "delta_max": float(d.max()), "fraction_p_below_0.01": float(np.mean(chi2.sf(np.clip(d, 0, None), 3) < 0.01)),
                                                      "fraction_all_positive": float(np.mean([all(v > 0 for v in x["quadrants"]) for x in L]))}
    print(s, json.dumps(summ[s]))
(ROOT / "results" / "fits" / "quadrant_lr_null_2026-09-24.json").write_text(json.dumps({"date": "2026-09-24", "summary": summ, "skies": out}, indent=1))
