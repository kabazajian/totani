#!/usr/bin/env python
"""The halo search with TOTANI'S OWN diffuse model: GALPROP SLZ6R30T150C2.

    python scripts/totani_diffuse.py        # ~4 min; two stack builds (~3 GB each)

Kev, 2026-09-18: "I want to use the same diffuse model as Totani." His model
(arXiv:2507.07209 sec 2.3) is GALPROP v54 with the Ackermann+ 2012 parameter
set SLZ6R30T150C2 -- the pre-computed output of which the LAT team published
(`totani.galprop` has the provenance) -- used as TWO templates: `gas` = pi0 +
bremsstrahlung with the model's ratio, and `ics`, each free per bin. His
sec 3.4.4 variant "gas rings & 3 ICS" splits gas into five Galactocentric
groups (Ackermann+ 2012 Table 2 annuli 1 | 2-5 | 6-12 | 13 | 14-17) and ICS
into its three ISRF components; both are built here.

⚠ RUN 1 (2026-09-18 evening) double-counted ionized bremsstrahlung (5% of
bremss): `bremss_healpix` already includes it. Fixed in `totani.galprop`;
this file's JSON is from the corrected run.

What this does: build the stack with `diffuse="galprop_totani"` in place of
the 14 A2020 ring templates, keep everything else as in our baseline
(structured bubbles, two Loop I shells, isotropic, 4FGL-DR4 point sources),
and score it the standard way. Then the same with the flat bubble mask, and
with the 3-component ICS. The A2020 baseline is re-scored alongside for a
like-for-like comparison on the same data, cells and point sources.

⚠ WHAT DIFFERS FROM TOTANI EVEN SO: his bubble template is his own
construction (sec 3.1, from the 4.3 GeV residual map), his point-source
treatment and his fit machinery are his; and his ROI cells are the same 10
degree cells (sec 2.2). So this is "our pipeline with his diffuse model", not
a reproduction of his fit. It answers: does the A2020-vs-GALPROP diffuse
choice move the halo, and does it change the quadrant test?
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from systematics import bubble_cell_template, quadrant_diagnostic, replace_component  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, injection_slope, restrict_stack  # noqa: E402
from totani.stack import build_background_stack, halo_cell_template  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9


def score(stack, J):
    d = quadrant_diagnostic(stack, J)
    keep = halo_cell_mask()
    sr, Jr = restrict_stack(stack, keep, J)
    d["injection_slopes"] = injection_slope(sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
    fit = fit_all_bins_profile(sr, Jr, halo_mask=keep)["per_bin"]
    roi_J = [float(Jr[k][keep].sum()) for k in range(N_BINS)]
    d["halo_intervals"] = [{"bin": k, "E_gev": TOTANI_CENTERS_GEV[k],
                            "counts": fit[k]["signal_amplitude"] * roi_J[k],
                            "lo": fit[k]["signal_amplitude_lo"] * roi_J[k],
                            "hi": fit[k]["signal_amplitude_hi"] * roi_J[k]} for k in range(N_BINS)]
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    d["neg2lnl_bg_only"] = [bg[k]["neg2lnl"] for k in range(N_BINS)]
    d["mean_abs_sigma_bg_only"] = [bg[k].get("mean_abs_sigma_halo") for k in range(N_BINS)]
    d["component_counts_bg_only"] = {
        n: [bg[k]["amplitudes"][n] * float(sr["templates"][n][k][keep].sum()) for k in range(N_BINS)]
        for n in sr["templates"]}
    d["components"] = sorted(stack["templates"])
    d["n_free_per_bin"] = len(stack["templates"])
    return d


def main() -> int:
    import resource
    t0 = time.time()
    J = halo_cell_template(RUN, EDGES)
    flat = bubble_cell_template(RUN, EDGES, "ackermann")

    stacks = {}
    a2020 = build_background_stack(RUN, EDGES, bubble="macias")
    stacks["a2020_structured"] = a2020
    stacks["a2020_flat"] = replace_component(a2020, "bubble", flat)
    print(f"A2020 stack ({time.time()-t0:.0f}s)", flush=True)
    del a2020

    tot = build_background_stack(RUN, EDGES, bubble="macias", diffuse="galprop_totani")
    stacks["totani_structured"] = tot
    stacks["totani_flat"] = replace_component(tot, "bubble", flat)
    print(f"Totani-model stack ({time.time()-t0:.0f}s)", flush=True)
    tot3 = build_background_stack(RUN, EDGES, bubble="macias", diffuse="galprop_totani", ics_split=True)
    stacks["totani_ics3_structured"] = tot3
    stacks["totani_ics3_flat"] = replace_component(tot3, "bubble", flat)
    print(f"Totani-model 3-ICS stack ({time.time()-t0:.0f}s)", flush=True)
    # His sec 3.4.4 "gas rings & 3 ICS": five Galactocentric gas groups + 3 ICS.
    totr = build_background_stack(RUN, EDGES, bubble="macias", diffuse="galprop_totani",
                                  ics_split=True, gas_rings=True)
    stacks["totani_rings3ics_structured"] = totr
    stacks["totani_rings3ics_flat"] = replace_component(totr, "bubble", flat)
    print(f"Totani-model rings+3-ICS stack ({time.time()-t0:.0f}s)", flush=True)

    results = {"date": "2026-09-18", "run": "fff2de5bc7e0",
               "galprop_model": "SLZ6R30T150C2 (Ackermann+ 2012), Totani 2507.07209 sec 2.3",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]), "variants": {}}
    for name, st in stacks.items():
        d = score(st, J)
        results["variants"][name] = d
        h5 = d["halo_intervals"][5]; q5 = d["summary"]["bin5"]
        print(f"  {name:24s} n={d['n_free_per_bin']:2d} halo@20.8 {h5['counts']:6.0f} "
              f"[{h5['lo']:.0f}, {h5['hi']:.0f}] TS {d['ts'][5]:5.1f} quad {q5['sum_over_baseline']:+.2f} "
              f"neg={q5['any_negative']} slope {min(d['injection_slopes']):.2f}-{max(d['injection_slopes']):.2f} "
              f"({time.time()-t0:.0f}s)", flush=True)
        print("      halo/bin   " + " ".join(f"{h['counts']:8.0f}" for h in d["halo_intervals"]))
        print("      TS         " + " ".join(f"{t:8.1f}" for t in d["ts"]))
        print("      <|sigma|>  " + " ".join(f"{(s if s is not None else float('nan')):8.2f}" for s in d["mean_abs_sigma_bg_only"]))
        print("      quad bin5  " + str([round(q) for q in q5["quadrants"]]))
    out = ROOT / "results" / "fits" / "totani_diffuse_2026-09-18.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out}\ntotal {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
