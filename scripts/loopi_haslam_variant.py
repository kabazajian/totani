#!/usr/bin/env python
"""Loop I from the Haslam 408 MHz map instead of the two geometric shells.

    python scripts/loopi_haslam_variant.py      # ~3 min, one 7 GB stack build

Kev's decision, 2026-09-18: an empirical alternative for Loop I, as a variant
beside the shells. Two questions:

(a) Does the structured map's low-energy disagreement with the shells
    (structured - flat = +243 at 1.5 GeV WITH the shells, -215 WITHOUT them;
    `results/bubble_lowE_check_2026-09-18.md`) go away with an empirical
    Loop I? Reported as delta(-2lnL) structured - flat per bin for each
    Loop I model.
(b) What does each Loop I model do to the halo, the quadrants and the
    injection gate? Standard scoring (`systematics.py`'s).

Loop I models:
    shells        two Wolleben shells (the baseline)
    haslam_nps    Haslam inside the shell-2 disc (`haslam.nps_disc`), one norm
    haslam_raw    Haslam over the whole ROI, one norm -- a synchrotron template
                  competing with the IC rings everywhere, not only at Loop I
    shells+nps    both: does the empirical shape carry something the shells miss?
    none          for reference (the tested-and-rejected candidate)

⚠ 408 MHz synchrotron is electrons x B-field; the gamma rays are electrons x
ISRF. Related, not identical. This is a modelling choice, scored, not a
replacement adopted.
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

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from systematics import LOOPI, drop, quadrant_diagnostic, replace_component  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.haslam import haslam_nps_on_target, haslam_on_target  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, injection_slope, restrict_stack  # noqa: E402
from totani.stack import _strip_to_2d, build_background_stack, cell_template_from_shape, halo_cell_template  # noqa: E402
from totani.templates import bubble_on_target  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9


def add(stack, new):
    t = dict(stack["templates"]); t.update(new)
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def bg_only(stack):
    keep = halo_cell_mask()
    sr = restrict_stack(stack, keep)
    fit = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    return fit, sr


def score(stack, J):
    d = quadrant_diagnostic(stack, J)
    keep = halo_cell_mask()
    sr, Jr = restrict_stack(stack, keep, J)
    d["injection_slopes"] = injection_slope(sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
    fit = fit_all_bins_profile(sr, Jr, halo_mask=keep)["per_bin"]
    roi_J = [float(Jr[k][keep].sum()) for k in range(N_BINS)]
    d["halo_intervals"] = [{"bin": k, "counts": fit[k]["signal_amplitude"] * roi_J[k],
                            "lo": fit[k]["signal_amplitude_lo"] * roi_J[k],
                            "hi": fit[k]["signal_amplitude_hi"] * roi_J[k]} for k in range(N_BINS)]
    bg, srb = bg_only(stack)
    d["neg2lnl_bg_only"] = [bg[k]["neg2lnl"] for k in range(N_BINS)]
    d["loopi_counts_bg_only"] = {
        n: [bg[k]["amplitudes"][n] * float(srb["templates"][n][k][keep].sum()) for k in range(N_BINS)]
        for n in stack["templates"] if n.startswith("loopI") or n.startswith("haslam")}
    d["components"] = sorted(stack["templates"])
    return d


def main() -> int:
    import resource
    t0 = time.time()
    bg = build_background_stack(RUN, EDGES, bubble="macias")
    J = halo_cell_template(RUN, EDGES)
    with fits.open(RUN / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    flat_cells = cell_template_from_shape(RUN, EDGES, (bubble_on_target(target) > 0).astype(float))
    nps = cell_template_from_shape(RUN, EDGES, haslam_nps_on_target(target))
    raw = cell_template_from_shape(RUN, EDGES, haslam_on_target(target))
    print(f"stack + templates ({time.time()-t0:.0f}s)", flush=True)

    no_shells = drop(bg, set(LOOPI))
    models = {
        "shells": bg,
        "haslam_nps": add(no_shells, {"haslam_nps": nps}),
        "haslam_raw": add(no_shells, {"haslam_raw": raw}),
        "shells+nps": add(bg, {"haslam_nps": nps}),
        "none": no_shells,
    }
    results = {"date": "2026-09-18", "run": "fff2de5bc7e0", "bubble_baseline": "macias",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]), "models": {}}
    for name, st in models.items():
        d = score(st, J)
        flat = replace_component(st, "bubble", flat_cells)
        fb, _ = bg_only(flat)
        d["neg2lnl_bg_only_flat_bubble"] = [fb[k]["neg2lnl"] for k in range(N_BINS)]
        d["delta_structured_minus_flat"] = [a - b for a, b in zip(d["neg2lnl_bg_only"], d["neg2lnl_bg_only_flat_bubble"])]
        results["models"][name] = d
        h5 = d["halo_intervals"][5]; q5 = d["summary"]["bin5"]
        print(f"  {name:11s} halo@20.8 {h5['counts']:6.0f} [{h5['lo']:.0f}, {h5['hi']:.0f}] TS {d['ts'][5]:5.1f} "
              f"quad {q5['sum_over_baseline']:+.2f} neg={q5['any_negative']} "
              f"slope {min(d['injection_slopes']):.2f}-{max(d['injection_slopes']):.2f} ({time.time()-t0:.0f}s)", flush=True)
        print("      halo/bin        " + " ".join(f"{h['counts']:8.0f}" for h in d["halo_intervals"]))
        print("      -2lnL bg-only   " + " ".join(f"{v:8.0f}" for v in d["neg2lnl_bg_only"]))
        print("      struct - flat   " + " ".join(f"{v:+8.0f}" for v in d["delta_structured_minus_flat"]))
    base = results["models"]["shells"]["neg2lnl_bg_only"]
    print("\n  -2lnL (bg-only, structured) relative to the shells, per bin:")
    for name, d in results["models"].items():
        print(f"  {name:11s} " + " ".join(f"{a - b:+8.0f}" for a, b in zip(d["neg2lnl_bg_only"], base))
              + f"   params: {d['n_bubble_params'] if 'n_bubble_params' in d else ''}"
              f"{len([n for n in d['components'] if n.startswith('loopI') or n.startswith('haslam')])} Loop I")
    out = ROOT / "results" / "fits" / "loopi_haslam_variant_2026-09-18.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out}\ntotal {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
