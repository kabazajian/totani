#!/usr/bin/env python
"""Does Galp21's 0.1 deg resolution make it fit worse than the 0.46 deg 2012
grid? Two filtered copies of each Galp21 model -- averaged into HEALPix
nside-128 pixels (exactly the 2012 discretization) and Gaussian-smoothed to
0.46 deg FWHM -- refit bare and under the split's no-residual selection.
Plus: how much small-scale structure the Galp21 and 2012 GAS maps carry at
4.3 GeV in the ROI (rms about a 2 deg smoothing), and the cell-to-cell scatter
of their ratio. Kev's question, 2026-09-23.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from reduce_totani import TOTANI_EDGES_GEV
from totani import ranking, selection, simsky
from totani.cells import halo_cell_mask, bin_to_cells
from totani.healpix import ang2pix_ring, npix
EDGES = TOTANI_EDGES_GEV[:10]
MODELS = ("SA0_R12", "SA100_F98", "SA0_F98", "SA50_F98")
SEL = {"gas1_ics1": ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic"),
       "gas5_ics3": ("bubble_flat", "erosita_unif_outer_NS", "loopI_haslam_raw", "sun_moon_ic")}
SIG_PX = 0.46 / 2.3548 / 0.125          # 0.46 deg FWHM in 0.125 deg pixels

def main():
    t0 = time.time()
    target, L, B, _ = ranking._grid_and_basis(EDGES)
    pix = ang2pix_ring(128, L, B).ravel(); n_hp = npix(128); cnt = np.bincount(pix, minlength=n_hp)
    def hp128(plane):
        s = np.bincount(pix, weights=plane.ravel(), minlength=n_hp)
        return (s / np.maximum(cnt, 1))[pix].reshape(plane.shape)
    missing = [m for m in MODELS if not ranking.available(f"G21_{m}")]
    if missing:
        sys.exit(f"Galp21 cubes not ingested for {missing}: run scripts/galp21_ingest.py first (the archives are "
                 f"available from the authors on request)")
    for m in MODELS:
        src = ranking.CACHE / f"galp21_cubes__{m}.npz"
        for tag, f in (("hp128", hp128), ("g046", lambda p: gaussian_filter(p, SIG_PX, mode="nearest"))):
            out = ranking.CACHE / f"galp21_cubes__{m}+{tag}.npz"
            if out.exists(): continue
            z = np.load(src)
            np.savez(out, **{k: np.stack([f(z[k][i].astype(float)) for i in range(z[k].shape[0])]).astype(np.float32) for k in z.files})
            print(f"  wrote {out.name} ({time.time()-t0:.0f}s)", flush=True)
    sh = ranking.shared_cells(EDGES); res = {}
    print("\nfit -2lnL (|b|>=10, 9 bins): native, then hp128 and g046 minus native")
    for m in MODELS:
        for v in ranking.VARIANTS:
            vals = {}
            for tag in ("", "+hp128", "+g046"):
                code = f"G21_{m}{tag}"
                st = selection.base_stack(code, v, EDGES, sh); bare = float(selection.bg_neg2lnl(st).sum())
                for c in SEL[v]: st, _, _ = selection.add_candidate(st, c, code, v, EDGES, sh)
                vals[tag or "native"] = (bare, float(selection.bg_neg2lnl(st).sum()))
            res[f"{m}__{v}"] = vals
            n = vals["native"]
            print(f"{m:10s} {v}: bare {n[0]:.0f}  hp128 {vals['+hp128'][0]-n[0]:+.1f}  g046 {vals['+g046'][0]-n[0]:+.1f}   | "
                  f"structured-set {n[1]:.0f}  hp128 {vals['+hp128'][1]-n[1]:+.1f}  g046 {vals['+g046'][1]-n[1]:+.1f}", flush=True)
    # gas structure at 4.3 GeV, pixel level, from the Stage E precompute (counts per pixel)
    print("\ngas map structure at 4.3 GeV, ROI |b|>=10: rms of (map/smooth2deg - 1), and 10-deg cell scatter of ln(G21/2012)")
    keep = halo_cell_mask(); shared = simsky.load_shared(); mask = np.asarray(shared["mask"])
    roi = (np.abs(B) >= 10) & ~mask
    ref = np.asarray(simsky.load_row("SL_Z6_R30_T150_C2", "gas1_ics1")["galprop_gas"], float)
    def structure(mp):
        sm = gaussian_filter(mp, 2.0 / 0.125, mode="nearest")
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(sm > 0, mp / sm - 1.0, 0.0)
        return float(np.sqrt(np.mean(r[roi] ** 2)))
    cref = bin_to_cells(ref, L, B, mask=mask)
    print(f"  2012 SL_Z6 gas: structure rms {structure(ref):.3f}")
    for m in ("SA0_R12", "SA100_F98"):
        try:
            g = np.asarray(simsky.load_row(f"G21_{m}", "gas1_ics1")["galprop_gas"], float)
        except FileNotFoundError:
            simsky.precompute_row(f"G21_{m}", "gas1_ics1", EDGES); g = np.asarray(simsky.load_row(f"G21_{m}", "gas1_ics1")["galprop_gas"], float)
        cg = bin_to_cells(g, L, B, mask=mask)
        lr = np.log(cg[keep] / cref[keep])
        print(f"  G21 {m:10s} gas: structure rms {structure(g):.3f};  cell ln(G21/2012): mean {lr.mean():+.3f} std {lr.std():.3f} range {lr.min():+.2f}..{lr.max():+.2f}", flush=True)
    (ROOT / "results" / "fits" / "galp21_smoothing_test_2026-09-23.json").write_text(json.dumps(res, indent=1))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
