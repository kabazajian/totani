#!/usr/bin/env python
"""Injection test of Totani's bubble construction: does a planted halo survive
being absorbed into the 4.3 GeV residual templates?

    python scripts/bubble_totani_injection.py     # ~6 min

The standard injection gate (`likelihood.injection_slope`) holds the
templates FIXED and asks whether the fitter recovers a planted halo: it does
(slope 0.99-1.00, `bubble_totani_2026-09-19`). But Totani's positive and
negative templates are BUILT FROM THE DATA at 4.3 GeV -- a real halo present
at 4.3 GeV goes into them by construction, and thereafter a template with
one free normalisation per bin carries the halo's morphology. The question
this answers: if a halo with a given spectrum is planted in the data and the
templates are rebuilt as he builds them, what fraction of it does the halo
fit recover, per energy bin?

Method: for injection amplitudes a in {0, a1, a2} (units of the NFW-rho^2
template, the same as `injection_slope`; a1 ~ the fitted 20.8 GeV amplitude
scale), add a x J(pixel) to the 4.3 GeV counts, redo the 144-cell
construction fit (flat mask + his model + GCE), rebuild the positive and
negative templates from the new residual, then add a x J(cells) to every
bin's counts and run the halo fit at |b| >= 10 with the rebuilt templates.
Recovered amplitude vs a, per bin, in covariance form -> slope. Alongside:
the same injections with templates FIXED at a = 0 (the standard gate). The
difference between the two slopes is what the construction absorbs.

Two injected spectra: FLAT in amplitude (the same a in every bin -- the
halo's counts then follow the exposure) and TOTANI-LIKE (a scaled per bin to
his fitted spectrum, `halo_intervals` of the a = 0 fit, so the planted halo
peaks at 20 GeV and is zero-ish below 3 GeV). The second is the case that
matters: does his construction let his own signal through?
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from bubble_totani import with_residual_templates  # noqa: E402
from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani.bubble_totani import bubble_image, fit_all_cells, pixel_model_for_bin  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.halo import nfw_rho2_on_target  # noqa: E402
from totani.likelihood import fit_all_bins_profile, restrict_stack  # noqa: E402
from totani.stack import build_background_stack, cell_template_from_shape, halo_cell_template  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9
K_BUILD = 2                      # 4.3 GeV, the bin the templates are built from
DATE = "2026-09-19"


def slope(inj, rec):
    inj = np.asarray(inj, float); rec = np.asarray(rec, float)
    d = inj - inj.mean()
    return float(np.sum(d * (rec - rec.mean())) / np.sum(d ** 2))


def main() -> int:
    import resource
    t0 = time.time()
    keep = halo_cell_mask()
    J = halo_cell_template(RUN, EDGES)                        # (n_bins, 12, 12) counts per unit amplitude
    pm = pixel_model_for_bin(RUN, EDGES, K_BUILD, diffuse="galprop_totani", with_gce=True)
    J_pix = nfw_rho2_on_target(pm["target"]) * pm["flat_sr"]  # pixel counts per unit amplitude at 4.3 GeV
    bg = build_background_stack(RUN, EDGES, bubble="macias", diffuse="galprop_totani")
    print(f"stack, pixel model, J ({time.time()-t0:.0f}s)", flush=True)

    def templates_for(counts_pix):
        pm2 = dict(pm); pm2["counts"] = counts_pix
        fit = fit_all_cells(pm2)
        img = bubble_image(pm2, fit)
        return (cell_template_from_shape(RUN, EDGES, img["positive"]),
                cell_template_from_shape(RUN, EDGES, img["negative"]), fit)

    def halo_fit(counts_cells, pos, neg):
        st = {"counts": counts_cells, "n_bins": N_BINS, "templates": bg["templates"]}
        st, signs = with_residual_templates(st, pos, neg)
        sr, Jr = restrict_stack(st, keep, J)
        fit = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False)["per_bin"]
        return [row["signal_amplitude"] for row in fit], signs

    # a = 0: the templates as built from the data, and the fitted spectrum
    pos0, neg0, fit0 = templates_for(pm["counts"])
    amp0, signs0 = halo_fit(bg["counts"], pos0, neg0)
    roi_J = [float(J[k][keep].sum()) for k in range(N_BINS)]
    print("a = 0 halo counts: " + " ".join(f"{a * j:8.0f}" for a, j in zip(amp0, roi_J)), flush=True)
    a_ref = amp0[5]                                            # the fitted 20.8 GeV amplitude
    inj_levels = [0.0, 1.0, 2.0]                               # in units of a_ref

    results = {"date": DATE, "run": "fff2de5bc7e0", "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]),
               "a_ref": a_ref, "roi_J": roi_J, "a0_halo_counts": [a * j for a, j in zip(amp0, roi_J)],
               "spectra": {}}
    spectra = {
        "flat": np.ones(N_BINS),
        # Totani-like: his fitted spectrum, positive part, normalised to 1 at 20.8 GeV
        "totani_like": np.clip(np.array(amp0) / a_ref, 0, None),
    }
    for sname, shape in spectra.items():
        rec_rebuilt, rec_fixed, signs_all = [], [], []
        for lev in inj_levels:
            a_bins = lev * a_ref * shape                        # amplitude per bin
            counts_cells = bg["counts"] + a_bins[:, None, None] * J
            counts_pix = pm["counts"] + a_bins[K_BUILD] * J_pix
            pos, neg, fit = templates_for(counts_pix) if lev > 0 else (pos0, neg0, fit0)
            amp_r, signs = halo_fit(counts_cells, pos, neg)
            amp_f, _ = halo_fit(counts_cells, pos0, neg0)
            rec_rebuilt.append(amp_r); rec_fixed.append(amp_f); signs_all.append(signs)
            print(f"  [{sname}] level {lev:.1f}: rebuilt-template halo counts "
                  + " ".join(f"{a * j:8.0f}" for a, j in zip(amp_r, roi_J)) + f"  ({time.time()-t0:.0f}s)", flush=True)
        rec_rebuilt = np.array(rec_rebuilt); rec_fixed = np.array(rec_fixed)
        inj = np.array([[lev * a_ref * shape[k] for k in range(N_BINS)] for lev in inj_levels])
        s_rebuilt = [slope(inj[:, k], rec_rebuilt[:, k]) if inj[:, k].std() > 0 else None for k in range(N_BINS)]
        s_fixed = [slope(inj[:, k], rec_fixed[:, k]) if inj[:, k].std() > 0 else None for k in range(N_BINS)]
        results["spectra"][sname] = {"shape": shape.tolist(), "levels": inj_levels,
                                     "injected_amplitude": inj.tolist(),
                                     "recovered_rebuilt": rec_rebuilt.tolist(), "recovered_fixed": rec_fixed.tolist(),
                                     "slope_rebuilt": s_rebuilt, "slope_fixed": s_fixed,
                                     "neg_signs": signs_all}
        print(f"  [{sname}] slope, templates REBUILT: " + " ".join("   n/a " if s is None else f"{s:7.2f}" for s in s_rebuilt))
        print(f"  [{sname}] slope, templates FIXED  : " + " ".join("   n/a " if s is None else f"{s:7.2f}" for s in s_fixed), flush=True)

    out = ROOT / "results" / "fits" / f"bubble_totani_injection_{DATE}.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out}\ntotal {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
