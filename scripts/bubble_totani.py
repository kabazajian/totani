#!/usr/bin/env python
"""Build Totani's bubble templates (sec 3.1) and run his halo search with them.

    python scripts/bubble_totani.py            # ~8 min; ~7 GB peak (A2020 stack)

Step 1 (`totani.bubble_totani`): at 1.5 and 4.3 GeV, fit the sharp-edged flat
mask with every other component -- his GALPROP diffuse model, the two Loop I
shells, isotropic, point sources, and the GC GeV excess (NFW-rho^2.5) -- over
all 144 cells, disk included. The bubble image is flat template + residual,
per pixel, in flux, smoothed with sigma = 1 deg. Saved as FITS and PNG (the
analogue of his figure 1). The 4.3 GeV image is split into POSITIVE and
NEGATIVE residual templates.

Step 2: the halo search at |b| >= 10 with the two residual templates in place
of any bubble template, halo = NFW-rho^2 as everywhere here (his figure 6).
His rule (sec 3.2): f >= 0 for the positive template, the negative template's
sign is free. This fitter keeps background amplitudes non-negative, so for
each bin the background-only fit is done with the negative template and with
its inverse, and the sign with the better -2lnL is kept for that bin. Run
under his diffuse model (his configuration end to end) and under A2020 (his
bubbles, our diffuse model).

The consistency check, from his sec 3.2: the residual templates absorb
everything at 4.3 GeV EXCEPT the GC excess, whose counts were subtracted in
the construction (it was a known component there). In his figure 5 the
rho^2.5 GCE recovers exactly its construction-fit flux at 4.3 GeV; in the
rho^2 fits (figure 6) the halo template replaces it, so at 4.3 GeV the halo
must be of order the GCE's |b| >= 10 counts from the construction fit (4,914
here), not zero. Measured 2026-09-19: 7,826 [6,002, 9,339] -- the same order,
1.6x, the flatter profile taking some residual as well. Printed and stored.

Not reproduced: his iteration of the flat boundary by eye (the Ackermann+
2014 mask is used as is), his MCMC, his PSF-convolved point sources.
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
from systematics import drop  # noqa: E402
from totani.bubble_totani import bubble_image, fit_all_cells, pixel_model_for_bin  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, restrict_stack  # noqa: E402
from totani.stack import build_background_stack, cell_template_from_shape, halo_cell_template  # noqa: E402
from totani_diffuse import score  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
N_BINS = 9
DATE = "2026-09-19"
OUT_FITS = ROOT / "results" / "fits" / f"bubble_totani_{DATE}.fits"


def build_images() -> dict:
    t0 = time.time()
    images, fits_out = {}, {}
    for k in (0, 2):
        pm = pixel_model_for_bin(RUN, EDGES, k, diffuse="galprop_totani", with_gce=True)
        fit = fit_all_cells(pm)
        img = bubble_image(pm, fit)
        images[k] = img
        fits_out[k] = fit
        v = img["flux_smooth"][img["valid"]]
        print(f"  bin {k} ({TOTANI_CENTERS_GEV[k]:.1f} GeV): 144-cell fit -2lnL {fit['neg2lnl']:.0f}, "
              f"converged {fit['converged']}, {pm['n_sources']} sources; bubble amp {fit['amplitudes']['bubble']:.3g}; "
              f"smoothed flux: pos {float((v > 0).mean()):.2f} of valid pixels, |min| {-v.min():.2e}, max {v.max():.2e} "
              f"({time.time()-t0:.0f}s)", flush=True)
        amps = fit["amplitudes"]
        print("      amplitudes: " + ", ".join(f"{n} {a:.3g}" for n, a in amps.items()))
    return images, fits_out


def save_images(images, fits_out, L, B):
    hdus = [fits.PrimaryHDU()]
    hdus[0].header["DATE"] = DATE
    hdus[0].header["RECIPE"] = "Totani 2507.07209 sec 3.1: flat mask + residual, flux, 1 deg Gaussian"
    hdus[0].header["DIFFUSE"] = "GALPROP SLZ6R30T150C2 gas+ics"
    for k, img in images.items():
        for key in ("flux_smooth", "flux_raw"):
            h = fits.ImageHDU(img[key].astype(np.float32), name=f"{key.upper()}_BIN{k}")
            h.header["EGEV"] = TOTANI_CENTERS_GEV[k]; h.header["BUBAMP"] = fits_out[k]["amplitudes"]["bubble"]
            h.header["BUNIT"] = "ph cm-2 s-1 sr-1 per bin"
            hdus.append(h)
    hdus.append(fits.ImageHDU(images[2]["positive"].astype(np.float32), name="POSITIVE"))
    hdus.append(fits.ImageHDU(images[2]["negative"].astype(np.float32), name="NEGATIVE"))
    hdus.append(fits.ImageHDU(images[2]["valid"].astype(np.uint8), name="VALID"))
    for h in hdus[1:]:
        h.header["CTYPE1"], h.header["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
        h.header["CRVAL1"], h.header["CRVAL2"] = 0.0, 0.0
        h.header["CDELT1"], h.header["CDELT2"] = -0.125, 0.125
        h.header["CRPIX1"], h.header["CRPIX2"] = 480.5, 480.5
    fits.HDUList(hdus).writeto(OUT_FITS, overwrite=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
        for ax, k in zip(axes, (0, 2)):
            s = images[k]["flux_smooth"].copy(); s[np.abs(B) < 10] = np.nan
            v = np.nanpercentile(np.abs(s), 99)
            im = ax.imshow(s, origin="lower", extent=[60, -60, -60, 60], cmap="RdBu_r", vmin=-v, vmax=v)
            ax.set_title(f"Totani-style bubble image, {TOTANI_CENTERS_GEV[k]:.1f} GeV (flat mask + residual, σ = 1°)")
            ax.set_xlabel("l [deg]"); ax.set_ylabel("b [deg]")
            plt.colorbar(im, ax=ax, label="flux per bin [ph cm$^{-2}$ s$^{-1}$ sr$^{-1}$]")
        fig.tight_layout()
        fig.savefig(ROOT / "results" / f"bubble_totani_{DATE}.png", dpi=110)
        print(f"  figure: results/bubble_totani_{DATE}.png")
    except Exception as e:  # matplotlib optional
        print(f"  (no figure: {e})")


def with_residual_templates(stack, pos_cells, neg_cells) -> tuple:
    """Bubble replaced by pos (+ neg with the per-bin sign the background-only
    fit prefers). Returns (stack, signs)."""
    keep = halo_cell_mask()
    base = drop(stack, {"bubble"})
    t = dict(base["templates"]); t["bubble_pos"] = pos_cells
    signs = []
    neg_final = np.empty_like(neg_cells)
    for k in range(N_BINS):
        best = None
        for sign in (+1, -1):
            tt = dict(t); tt["bubble_neg"] = sign * neg_cells
            st = {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": tt}
            sr = restrict_stack(st, keep)
            v = fit_all_bins(sr, halo_mask=keep)["per_bin"][k]["neg2lnl"]
            if best is None or v < best[0]:
                best = (v, sign)
        signs.append(best[1]); neg_final[k] = best[1] * neg_cells[k]
    t["bubble_neg"] = neg_final
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}, signs


def main() -> int:
    import resource
    t0 = time.time()
    print("step 1: bubble images", flush=True)
    images, fits_out = build_images()
    from totani.templates import grid_lb, target_header
    L, B = grid_lb(target_header())
    save_images(images, fits_out, L, B)
    pos_shape, neg_shape = images[2]["positive"], images[2]["negative"]
    from totani.templates import bubble_on_target
    fb = bubble_on_target(target_header()) > 0
    print(f"  positive template covers {float((pos_shape > 0).mean()):.2f} of the grid, "
          f"negative {float((neg_shape > 0).mean()):.2f}; of the flat-mask footprint, "
          f"{float((pos_shape[fb] > 0).mean()):.2f} is positive in the 4.3 GeV image", flush=True)

    print("step 2: halo search with the residual templates", flush=True)
    J = halo_cell_template(RUN, EDGES)
    pos_cells = cell_template_from_shape(RUN, EDGES, pos_shape)
    neg_cells = cell_template_from_shape(RUN, EDGES, neg_shape)
    results = {"date": DATE, "run": "fff2de5bc7e0", "E_centers_gev": list(TOTANI_CENTERS_GEV[:N_BINS]),
               "construction": {str(k): {"neg2lnl_144cells": f["neg2lnl"], "amplitudes": f["amplitudes"]}
                                for k, f in fits_out.items()},
               "variants": {}}
    for tag, kw in (("totani_diffuse", dict(diffuse="galprop_totani")), ("a2020_diffuse", {})):
        bg = build_background_stack(RUN, EDGES, bubble="macias", **kw)
        st, signs = with_residual_templates(bg, pos_cells, neg_cells)
        d = score(st, J); d["neg_template_sign_per_bin"] = signs
        results["variants"][f"{tag}_posneg"] = d
        h = d["halo_intervals"]; q5 = d["summary"]["bin5"]
        print(f"  {tag:15s} + pos/neg: halo@20.8 {h[5]['counts']:6.0f} [{h[5]['lo']:.0f}, {h[5]['hi']:.0f}] "
              f"TS {d['ts'][5]:5.1f} quad {q5['sum_over_baseline']:+.2f} neg={q5['any_negative']} "
              f"slope {min(d['injection_slopes']):.2f}-{max(d['injection_slopes']):.2f}  neg-sign {signs} ({time.time()-t0:.0f}s)", flush=True)
        print("      halo/bin   " + " ".join(f"{x['counts']:8.0f}" for x in h))
        print("      TS         " + " ".join(f"{t:8.1f}" for t in d["ts"]))
        print(f"      4.3 GeV consistency: halo {h[2]['counts']:.0f} [{h[2]['lo']:.0f}, {h[2]['hi']:.0f}] "
              f"(expect ~ the GCE's |b|>=10 counts from the construction fit, ~5e3 under his model)")
        print("      quad bin5  " + str([round(x) for x in q5["quadrants"]]))
        print("      E x counts " + " ".join(f"{x['counts'] * e / 1e3:8.0f}" for x, e in zip(h, TOTANI_CENTERS_GEV)))
        # positive-only variant: his "resid.(-) masked" is a mask, not this; but
        # dropping the negative template shows what it carries
        st2 = drop(st, {"bubble_neg"}); d2 = score(st2, J)
        results["variants"][f"{tag}_pos_only"] = d2
        print(f"  {tag:15s} + pos only: halo@20.8 {d2['halo_intervals'][5]['counts']:6.0f} TS {d2['ts'][5]:5.1f} "
              f"quad {d2['summary']['bin5']['sum_over_baseline']:+.2f}", flush=True)
        del bg, st, st2

    out = ROOT / "results" / "fits" / f"bubble_totani_{DATE}.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out} and {OUT_FITS}\ntotal {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
