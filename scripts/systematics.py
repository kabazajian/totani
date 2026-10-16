#!/usr/bin/env python
"""Totani's sec 3.4 systematics, plus the variants his quadrant test implicates.

    python scripts/systematics.py                 # structured bubble baseline (2026-09-07 on)
    python scripts/systematics.py --bubble ackermann   # the 2026-09-04 table

🚨 THE BASELINE BUBBLE TEMPLATE CHANGED ON 2026-09-07 (Kev). Oscar's structured
map is the more physical template -- the bubbles as the data show them rather
than a flat geometric mask -- and it fits better with the same parameter count
at 4-7 and 35-60 GeV (`results/bubble_macias_2026-09-07.md`). It is now the
baseline; the flat mask stays in the table as the variant `bubble_flat`. The
output file is dated so the earlier table remains the record of what the flat
baseline gave. The one place the structured map fits WORSE, 1.5-2.6 GeV, is
where the background is independently known to be inadequate.

🚨 EVERY VARIANT IS SCORED BY THE QUADRANT DIAGNOSTIC, NOT BY ITS HALO
AMPLITUDE. The v3 baseline already produces a positive halo-shaped amplitude
in every bin; what refutes reading it as a halo is that the four quadrants
disagree by ~3x, one goes NEGATIVE, and they sum to roughly half the baseline
(sec 3.4.1 -- a spherical halo cannot do that). So the question each variant
answers is: does REMOVING or CHANGING this component make the quadrants agree?
Whatever does is the component the "halo" was absorbing.
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
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import (  # noqa: E402
    fit_all_bins,
    fit_all_bins_profile,
    injection_slope,
    quadrant_masks,
    restrict_stack,
)
from totani.pointsources import psf95_radius_deg  # noqa: E402
from totani.stack import build_background_stack, halo_cell_template  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
OUT = {"macias": ROOT / "results" / "fits" / "systematics_2026-09-07.json",
       "ackermann": ROOT / "results" / "fits" / "systematics_2026-09-04.json"}
#: ⚠ `--no-loopI` IS A TESTED-AND-NOT-ADOPTED CANDIDATE (2026-09-18). Kev asked
#: for structured bubbles with no Loop I shells as the baseline, on the
#: reading that Oscar's 1-10 GeV residual map (Ackermann+ 2017 Fig. 12) fights
#: the shells at low energy (`results/bubble_lowE_check_2026-09-18.md`). The
#: run showed the shells are NOT double-counted by the map -- they fit only
#: 8-15% lower beside it and buy delta(-2lnL) = -2,986 at 1.5 GeV -- and that
#: dropping them raises the 1.5 GeV halo 23%, doubles the eROSITA residual and
#: worsens the quadrant test (`results/rebaseline_candidate_noLoopI_2026-09-18.md`).
#: The flag drops the shells from the baseline and every variant, turns the
#: old `no_loopI` variant into `with_loopI` (the 09-07 baseline), and writes
#: its own dated file. THE 09-07 BASELINE STANDS.
LOOPI = ("loopI_shell1", "loopI_shell2")
OUT_NO_LOOPI = ROOT / "results" / "fits" / "systematics_2026-09-18.json"


def bubble_cell_template(run_dir, edges_gev, which: str):
    """One bubble template through the diffuse pipeline (exposure x solid
    angle, extended-source mask, cells) WITHOUT rebuilding the stack -- so the
    other template can be swapped in as a variant at cell-space cost."""
    from astropy.io import fits

    from totani.cells import bin_to_cells
    from totani.exposure import flat_spectrum_bin_exposure
    from totani.pointsources import default_catalog_path, extended_source_mask
    from totani.stack import _strip_to_2d
    from totani.templates import (bubble_macias_on_target, bubble_on_target,
                                  grid_lb, pixel_solid_angle)

    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure = h[0].data.astype(float)
        emev = h["ENERGIES"].data["Energy"].astype(float)
    L, B = grid_lb(target)
    with fits.open(default_catalog_path()) as h:
        mask = extended_source_mask(L, B, h["ExtendedSources"].data)
    shape = bubble_on_target(target) if which == "ackermann" else bubble_macias_on_target(target)
    cube = (shape[None, :, :] * flat_spectrum_bin_exposure(exposure, emev, edges_gev)
            * pixel_solid_angle(target)[None, :, :])
    return bin_to_cells(cube, L, B, mask=mask)


def replace_component(stack, name, arr):
    t = dict(stack["templates"]); t[name] = arr
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def split_ns(stack, name="bubble"):
    from totani.cells import cell_centers
    _, Bc = cell_centers()
    t = dict(stack["templates"]); b = t.pop(name)
    t[f"{name}_N"] = b * (Bc > 0)[None, :, :]
    t[f"{name}_S"] = b * (Bc < 0)[None, :, :]
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def drop(stack, names):
    """Stack with components removed -- a background-model variant."""
    return {"counts": stack["counts"], "n_bins": stack["n_bins"],
            "templates": {k: v for k, v in stack["templates"].items()
                          if k not in names}}


def quadrant_diagnostic(stack, J, bins=(0, 2, 5)):
    """Per-variant: baseline amplitude, the four quadrant amplitudes, and the
    three numbers that decide whether the component behaves like a halo."""
    out = {}
    base_keep = halo_cell_mask()
    st, Jr = restrict_stack(stack, base_keep, J)
    base = fit_all_bins_profile(st, Jr, halo_mask=base_keep, with_uncertainty=False)
    b0 = fit_all_bins(st, halo_mask=base_keep)
    base_counts = [row["signal_amplitude"] * float(Jr[row["bin"]][base_keep].sum())
                   for row in base["per_bin"]]

    quads = {}
    for qn, qm in quadrant_masks().items():
        stq, Jq = restrict_stack(stack, qm, J)
        rq = fit_all_bins_profile(stq, Jq, halo_mask=qm, with_uncertainty=False)
        quads[qn] = [row["signal_amplitude"] * float(Jq[row["bin"]][qm].sum())
                     for row in rq["per_bin"]]

    out["baseline_counts"] = base_counts
    out["quadrant_counts"] = quads
    out["ts"] = [b0["per_bin"][k]["neg2lnl"] - base["per_bin"][k]["neg2lnl"]
                 for k in range(stack["n_bins"])]
    out["mean_abs_sigma"] = [row["mean_abs_sigma_halo"] for row in base["per_bin"]]
    summary = {}
    for k in bins:
        qv = np.array([quads[q][k] for q in quads])
        summary[f"bin{k}"] = {
            "sum_over_baseline": float(qv.sum() / base_counts[k]) if base_counts[k] else None,
            "spread_max_over_min": float(qv.max() / qv.min()) if qv.min() > 0 else None,
            "any_negative": bool((qv < 0).any()),
            "quadrants": qv.tolist()}
    out["summary"] = summary
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bubble", choices=("macias", "ackermann"), default="macias")
    ap.add_argument("--no-loopI", action="store_true",
                    help="2026-09-18 baseline: drop both Loop I shells everywhere")
    args = ap.parse_args()
    bubble, no_loopi = args.bubble, args.no_loopI
    other = "ackermann" if bubble == "macias" else "macias"
    if no_loopi and bubble != "macias":
        raise SystemExit("--no-loopI is defined for the structured baseline only")

    t0 = time.time()
    bg_full = build_background_stack(RUN, EDGES, bubble=bubble)
    bg = drop(bg_full, set(LOOPI)) if no_loopi else bg_full
    J = halo_cell_template(RUN, EDGES)
    print(f"base stack (bubble={bubble}, loopI={'out' if no_loopi else 'in'}): "
          f"{bg_full['n_sources_kept']} sources kept, "
          f"{bg_full['n_sources_dropped_in_mask']} dropped in mask "
          f"({time.time()-t0:.0f}s)", flush=True)

    variants = {
        "baseline": bg,
        f"bubble_{'flat' if other == 'ackermann' else 'structured'}":
            replace_component(bg, "bubble", bubble_cell_template(RUN, EDGES, other)),
        "bubble_NS": split_ns(bg),
        "no_bubble": drop(bg, {"bubble"}),
        "no_dust": drop(bg, {"dust_pos", "dust_neg"}),
    }
    if no_loopi:
        variants["with_loopI"] = bg_full                       # the 09-07 baseline
        variants["with_loopI_bubble_flat"] = replace_component(
            bg_full, "bubble", bubble_cell_template(RUN, EDGES, other))
    else:
        variants["no_loopI"] = drop(bg, set(LOOPI))
        variants["no_loopI_no_bubble"] = drop(bg, set(LOOPI) | {"bubble"})

    # Totani sec 3.4.3: mask sources instead of modelling them. Two radii, so
    # the conclusion is demonstrably not an artefact of the approximate 95%
    # containment scaling (see pointsources.PSF95_OVER_PSF68).
    for label, radii in (("ps_masked_95pct",
                          [psf95_radius_deg(e) for e in TOTANI_CENTERS_GEV[:9]]),
                         ("ps_masked_fixed_1deg", 1.0)):
        print(f"building {label} ...", flush=True)
        variants[label] = build_background_stack(
            RUN, EDGES, ps_mask_radii_deg=radii, bubble=bubble,
            drop_components=LOOPI if no_loopi else ())

    results = {"date": ("2026-09-18" if no_loopi else
                        "2026-09-07" if bubble == "macias" else "2026-09-04"),
               "run": "fff2de5bc7e0", "bubble_baseline": bubble,
               "loopI_in_baseline": not no_loopi,
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]), "variants": {}}
    for name, st in variants.items():
        d = quadrant_diagnostic(st, J)
        keep = halo_cell_mask()
        sr, Jr = restrict_stack(st, keep, J)
        d["injection_slopes"] = injection_slope(
            sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
        d["components"] = sorted(st["templates"])
        # 68% intervals for every bin, so the envelope is quoted with errors
        fit = fit_all_bins_profile(sr, Jr, halo_mask=keep)["per_bin"]
        d["halo_intervals"] = [
            {"bin": k, "counts": fit[k]["signal_amplitude"] * float(Jr[k][keep].sum()),
             "lo": fit[k]["signal_amplitude_lo"] * float(Jr[k][keep].sum()),
             "hi": fit[k]["signal_amplitude_hi"] * float(Jr[k][keep].sum())}
            for k in range(9)]
        results["variants"][name] = d
        s0 = d["summary"]["bin0"]
        print(f"  {name:22s} bin0 sum/base={s0['sum_over_baseline']:.2f} "
              f"neg={s0['any_negative']} slope="
              f"{min(d['injection_slopes']):.2f}-{max(d['injection_slopes']):.2f} "
              f"({time.time()-t0:.0f}s)", flush=True)

    out = OUT_NO_LOOPI if no_loopi else OUT[bubble]
    out.write_text(json.dumps(results, indent=2))
    print(f"saved {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
