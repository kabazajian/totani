#!/usr/bin/env python
"""Oscar's structured bubble template against the Ackermann flat mask.

    python scripts/bubble_macias.py

🚨 WHY THIS MATTERS MORE THAN ONE TEMPLATE SWAP USUALLY WOULD. Removing the
bubble template DOUBLES the fitted halo (7,705 -> 15,242 at 20.8 GeV,
`systematics_2026-09-04`), so the halo is absorbing bubble photons that the
flat 2014 mask does not claim. And giving that mask independent north/south
normalizations removes essentially all matched-theta asymmetry at 20.8 GeV
while moving the halo 0.7% (2026-09-05). So the questions for Oscar's map are
three: does it move the halo amplitude, does it fix the quadrant and
matched-theta failures, and does it do so with ONE normalization rather than
two.

Tier 1 on the delivery is in `templates.bubble_macias_on_target`. The variants
here are cell-space surgery on ONE stack build (7 GB peak, once): the bubble
column is replaced, split, or dropped; nothing else changes, so every
difference between rows is the bubble template and nothing else.
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
from systematics import drop, quadrant_diagnostic  # noqa: E402
from totani.cells import bin_to_cells, cell_centers, halo_cell_mask  # noqa: E402
from totani.degeneracy import leverage  # noqa: E402
from totani.exposure import flat_spectrum_bin_exposure  # noqa: E402
from totani.likelihood import (  # noqa: E402
    expected_counts,
    fit_all_bins,
    fit_all_bins_profile,
    injection_slope,
    restrict_stack,
)
from totani.pointsources import default_catalog_path, extended_source_mask  # noqa: E402
from totani.stack import _strip_to_2d, build_background_stack, halo_cell_template  # noqa: E402
from totani.templates import bubble_macias_on_target, grid_lb, pixel_solid_angle  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
OUT = ROOT / "results" / "fits" / "bubble_macias_2026-09-07.json"
#: Cell-space stack + Oscar's template + halo template, so every re-analysis
#: after the first is free and no second 7 GB build is ever needed.
CELLS = ROOT / "results" / "fits" / "bubble_macias_cells_2026-09-07.npz"


def macias_cell_template(run_dir, edges_gev) -> np.ndarray:
    """Oscar's map through the SAME pipeline as every other diffuse template
    (flat-spectrum exposure x solid angle, extended-source mask, cells) --
    mirrors `stack.halo_cell_template` so the stack need not be rebuilt."""
    run_dir = Path(run_dir)
    with fits.open(run_dir / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(run_dir / "expcube.fits") as h:
        exposure = h[0].data.astype(float)
        emev = h["ENERGIES"].data["Energy"].astype(float)
    L, B = grid_lb(target)
    with fits.open(default_catalog_path()) as h:
        mask = extended_source_mask(L, B, h["ExtendedSources"].data)
    cube = (bubble_macias_on_target(target)[None, :, :]
            * flat_spectrum_bin_exposure(exposure, emev, edges_gev)
            * pixel_solid_angle(target)[None, :, :])
    return bin_to_cells(cube, L, B, mask=mask)


def split_ns(stack, name="bubble"):
    _, Bc = cell_centers()
    t = dict(stack["templates"])
    b = t.pop(name)
    t[f"{name}_N"] = b * (Bc > 0)[None, :, :]
    t[f"{name}_S"] = b * (Bc < 0)[None, :, :]
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def replace(stack, name, arr):
    t = dict(stack["templates"]); t[name] = arr
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def on_off(stack, keep, k):
    """Matched-theta on-axis (|l|<=15) vs off-axis (|l|>=35) residual
    differences by hemisphere and annulus, background-only fit -- the
    2026-09-05 diagnostic that the N/S bubble split was found to cure."""
    Lc, Bc = cell_centers()
    th = np.degrees(np.arccos(np.clip(np.cos(np.radians(Bc)) * np.cos(np.radians(Lc)), -1, 1)))
    st, _ = restrict_stack(stack, keep, None) if False else (restrict_stack(stack, keep), None)
    names = sorted(st["templates"])
    T = np.stack([st["templates"][n] for n in names])
    r = fit_all_bins(st, halo_mask=keep)["per_bin"][k]
    mu = expected_counts(np.log10(np.clip(np.array([r["amplitudes"][n] for n in names]), 1e-300, None)), T[:, k])
    n = stack["counts"][k]
    out = {}
    for hemi, hm in (("N", Bc > 0), ("S", Bc < 0)):
        for a, b in ((25, 40), (40, 55)):
            m = keep & hm & (th >= a) & (th < b)
            on, off = m & (np.abs(Lc) <= 15), m & (np.abs(Lc) >= 35)
            son = (n[on].sum() - mu[on].sum()) / np.sqrt(mu[on].sum())
            sof = (n[off].sum() - mu[off].sum()) / np.sqrt(mu[off].sum())
            out[f"{hemi}_{a}_{b}"] = float(son - sof)
    return out


def main() -> int:
    import resource

    t0 = time.time()
    if CELLS.exists():
        z = np.load(CELLS)
        names0 = list(z["names"])
        base = {"counts": z["counts"], "n_bins": int(z["counts"].shape[0]),
                "templates": {n: z["T"][i] for i, n in enumerate(names0)}}
        J, mac = z["J"], z["macias"]
        print(f"loaded cached cells from {CELLS.name}", flush=True)
    else:
        base = build_background_stack(RUN, EDGES)          # bubble="ackermann"
        J = halo_cell_template(RUN, EDGES)
        mac = macias_cell_template(RUN, EDGES)
        names0 = sorted(base["templates"])
        np.savez_compressed(CELLS, counts=base["counts"], names=np.array(names0),
                            T=np.stack([base["templates"][n] for n in names0]),
                            J=J, macias=mac)
        print(f"stack + Macias template built and cached ({time.time()-t0:.0f}s)", flush=True)

    keep = halo_cell_mask()
    _, Bc = cell_centers()
    band = keep & (np.abs(Bc) < 30.0)
    variants = {
        "ackermann":    base,
        "macias":       replace(base, "bubble", mac),
        "ackermann_NS": split_ns(base),
        "macias_NS":    split_ns(replace(base, "bubble", mac)),
        "no_bubble":    drop(base, {"bubble"}),
    }

    results = {"date": "2026-09-07", "run": "fff2de5bc7e0",
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]), "variants": {}}
    for name, st in variants.items():
        sr, Jr = restrict_stack(st, keep, J)
        fit = fit_all_bins_profile(sr, Jr, halo_mask=keep)["per_bin"]
        halo = [{"bin": k,
                 "counts": fit[k]["signal_amplitude"] * float(Jr[k][keep].sum()),
                 "lo": fit[k]["signal_amplitude_lo"] * float(Jr[k][keep].sum()),
                 "hi": fit[k]["signal_amplitude_hi"] * float(Jr[k][keep].sum()),
                 "mean_abs_sigma": fit[k]["mean_abs_sigma_halo"]} for k in range(9)]
        q = quadrant_diagnostic(st, J)
        sl = injection_slope(sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
        # Is the template a better BACKGROUND, or just a more halo-shaped one?
        # Same data, same free-parameter count (ackermann vs macias), so the
        # background-only -2lnL is directly comparable between those two rows.
        bg_only = fit_all_bins(sr, halo_mask=keep)["per_bin"]
        neg2lnl_bg = [float(bg_only[k]["neg2lnl"]) for k in range(9)]
        ts_halo = [float(bg_only[k]["neg2lnl"] - fit[k]["neg2lnl"]) for k in range(9)]
        # degeneracy of the halo with the bubble column(s) and with the whole basis, bin 5
        names = sorted(sr["templates"])
        A = np.stack([sr["templates"][n][5] for n in names])[:, keep].T
        h = Jr[5][keep]
        bcols = [n for n in names if n.startswith("bubble")]
        corr = {}
        for n in bcols:
            c = sr["templates"][n][5][keep]
            corr[n] = float(np.corrcoef(h, c)[0, 1]) if c.std() > 0 else None
        Ab = np.stack([sr["templates"][n][5] for n in names])[:, band].T
        results["variants"][name] = {
            "components": names, "halo": halo,
            "neg2lnl_background_only": neg2lnl_bg, "ts_halo": ts_halo,
            "halo_leverage_band_bin5": leverage(Jr[5][band], Ab),
            "quadrant_summary": q["summary"], "quadrant_counts": q["quadrant_counts"],
            "injection_slope_min": float(min(sl)), "injection_slope_max": float(max(sl)),
            "on_off_bin3": on_off(st, keep, 3), "on_off_bin5": on_off(st, keep, 5),
            "halo_bubble_corr_bin5": corr,
            "halo_leverage_bin5": leverage(h, A),
            "bubble_fitted_fraction_bin5": {
                n: float(fit[5]["bg_amplitudes"][n] * sr["templates"][n][5][keep].sum()
                         / fit[5]["predicted_total"]) for n in bcols},
        }
        h5 = halo[5]; s5 = q["summary"]["bin5"]
        print(f"  {name:13s} bg-only -2lnL bin5 {neg2lnl_bg[5]:.1f}  TS_halo {ts_halo[5]:.1f}"
              f"  band lev {results['variants'][name]['halo_leverage_band_bin5']:.3f}", flush=True)
        print(f"  {name:13s} bin5 halo {h5['counts']:7.0f} [{h5['lo']:.0f},{h5['hi']:.0f}]"
              f"  quad sum/base {s5['sum_over_baseline']:.2f} spread {s5['spread_max_over_min'] or float('nan'):.2f}"
              f"  slope {min(sl):.3f}-{max(sl):.3f}  lev {results['variants'][name]['halo_leverage_bin5']:.3f}"
              f"  ({time.time()-t0:.0f}s)", flush=True)

    OUT.write_text(json.dumps(results, indent=2))
    print(f"saved {OUT}")
    print("peak RSS: %.2f GB" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
