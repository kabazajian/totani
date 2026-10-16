#!/usr/bin/env python
"""Audit (2026-09-30): is every 4FGL-DR4 source that should be in the fit
actually in the point-source template, at every binning used in the paper?

    python scripts/audit_point_sources.py

Recomputes, independently of the cached templates, each catalog source's
expected counts (catalog spectrum x exposure, integrated per bin) and checks:
  1. the catalog: version, row count, every SpectrumType is one `dnde` handles,
     no source with non-finite or non-positive expected counts;
  2. accounting: every source in |l|,|b| <= 60 is either kept or dropped for a
     stated reason (inside an extended-source mask circle); nothing else;
  3. the cached templates contain exactly the kept sources: per-bin totals of
     the 10-deg shared template, the 2-deg and 1-deg libraries (unconvolved),
     equal the recomputed kept-source sums; the PSF-folded libraries equal them
     minus what the PSF carries out of the |l|,|b| <= 60 grid (measured by
     folding the source map itself);
  4. placement: in the 10-deg template each cell holds exactly the sources
     whose positions fall in it;
  5. what is NOT modeled: sources outside the grid whose PSF spills into it,
     and the flux of the dropped (masked) sources in the fitted |b| >= 10 cells;
  6. every stack used by the fits includes the template (base_stack).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine, ranking  # noqa: E402
from totani.pointsources import (SPECTRUM_TYPES, cell_index, default_catalog_path, exposure_pixel,  # noqa: E402
                                 extended_source_mask, source_expected_counts)
from totani.stack import _strip_to_2d  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
OUT = ROOT / "results" / "audit_point_sources_2026-10-01.md"


def main() -> int:
    import os
    L = []
    say = lambda s="": (L.append(s), print(s))
    p = default_catalog_path()
    with fits.open(p) as h:
        cat = h["LAT_Point_Source_Catalog"].data; ext = h["ExtendedSources"].data
        ver = h[1].header.get("VERSION", "?")
    with fits.open(ranking.RUN / "expcube.fits") as h:
        expo = h[0].data.astype(float); emev = h["ENERGIES"].data["Energy"].astype(float); hdr = _strip_to_2d(h[0].header)
    say(f"# Point-source audit ({p.name}, catalog VERSION {ver}; 2026-10-01, exact-position placement, FINE_TAG={os.environ.get('FINE_TAG')})"); say()
    cl = np.where(cat["GLON"] > 180, cat["GLON"] - 360, cat["GLON"]).astype(float); cb = cat["GLAT"].astype(float)
    types = sorted(set(str(t).strip() for t in cat["SpectrumType"]))
    say(f"## 1. Catalog"); say()
    say(f"- rows: {len(cat)}; extended sources listed: {len(ext)}")
    say(f"- SpectrumType values: {types}; handled: {sorted(SPECTRUM_TYPES)}; unhandled: {[t for t in types if t not in SPECTRUM_TYPES] or 'none'}")
    grid = (np.abs(cl) <= 60) & (np.abs(cb) <= 60)
    in_mask = extended_source_mask(cl, cb, ext)
    keep = grid & ~in_mask
    rows = cat[grid]
    counts = source_expected_counts(rows, cl[grid], cb[grid], EDGES, expo, emev, hdr)
    bad = ~np.isfinite(counts) | (counts <= 0)
    say(f"- sources in |l|,|b| <= 60: {grid.sum()}; with a non-finite or non-positive expected count in any bin: {int(bad.any(1).sum())}")
    if bad.any():
        for i in np.where(bad.any(1))[0][:20]:
            say(f"  - {rows['Source_Name'][i]}: {counts[i]}")
    say(); say("## 2. Accounting"); say()
    fitted = grid & (np.abs(cb) >= 10)
    say(f"- kept: {keep.sum()} (in the fitted |b| >= 10 region: {(keep & fitted).sum()})")
    say(f"- dropped because inside an extended-source mask circle: {(grid & in_mask).sum()} (in |b| >= 10: {(grid & in_mask & fitted).sum()})")
    say(f"- dropped for any other reason: {grid.sum() - keep.sum() - (grid & in_mask).sum()}")
    kc = counts[keep[grid]]
    names_masked_hi = rows["Source_Name"][(in_mask[grid]) & (np.abs(cb[grid]) >= 10)]
    say(f"- masked sources at |b| >= 10: {', '.join(str(n).strip() for n in names_masked_hi)}")
    tot_kept = kc.sum(0)
    say(); say("## 3. Cached templates against the recomputed kept-source totals (per bin, ratio)"); say()
    say("| template | " + " | ".join(f"{e:.3g}" for e in np.sqrt(np.array(EDGES[:-1]) * np.array(EDGES[1:]))) + " |")
    say("|---|" + "---|" * 9)
    sh = ranking.shared_cells(EDGES)
    say("| 10° shared (all 144 cells) | " + " | ".join(f"{a / b:.6f}" for a, b in zip(sh["point_sources"].reshape(9, -1).sum(1), tot_kept)) + " |")
    for cell in (2.0, 1.0):
        plain = np.load(ranking.CACHE / f"cells{cell:g}" / "library.npz")
        t = plain["T__point_sources"].reshape(9, -1).sum(1)
        say(f"| {cell:g}° unconvolved | " + " | ".join(f"{a / b:.6f}" for a, b in zip(t, tot_kept)) + " |")
    # PSF-folded: compare with an independent placement of the kept sources (exact positions, 2026-10-01)
    import os
    os.environ["FINE_TAG"] = os.environ.get("FINE_TAG") or "psfuc"
    g = fine._psf_grid(EDGES)
    from totani import psf
    conv = psf.point_sources_exact(cl[keep], cb[keep], kc, g["L"], g["B"], EDGES)
    kept_after_fold = conv.reshape(9, -1).sum(1)
    say("| exact placement of the kept sources, flux on the grid / before (renormalised per source) | " + " | ".join(f"{a / b:.6f}" for a, b in zip(kept_after_fold, tot_kept)) + " |")
    for cell in (2.0, 1.0):
        z = np.load(fine.cache_dir(cell) / "library.npz")
        from totani.cells import bin_to_cells
        ref = np.stack([bin_to_cells(conv[k], g["L"], g["B"], cell_deg=cell, mask=g["mask"]) for k in range(9)])
        t = z["T__point_sources"]
        say(f"| {cell:g}° PSF-folded library / independent fold | " + " | ".join(f"{a / b:.6f}" for a, b in zip(t.reshape(9, -1).sum(1), ref.reshape(9, -1).sum(1))) + " |")
        say(f"| {cell:g}° PSF-folded, max cell deviation from independent fold | " + " | ".join(f"{np.abs(t[k] - ref[k]).max() / max(ref[k].max(), 1e-30):.2e}" for k in range(9)) + " |")
    say(); say("## 4. Placement in the 10° template"); say()
    r, c = cell_index(cl[keep], cb[keep], 60.0, 10.0)
    rebuilt = np.zeros((9, 12, 12))
    for k in range(9):
        np.add.at(rebuilt[k], (r, c), kc[:, k])
    dev = np.abs(rebuilt - sh["point_sources"]).max() / sh["point_sources"].max()
    say(f"- max |rebuilt − cached| / max cell, over all bins and cells: {dev:.2e}")
    say(); say("## 5. What is not modeled"); say()
    ring = (~grid) & (np.abs(cl) <= 63) & (np.abs(cb) <= 63)
    n_ring = int(ring.sum())
    say(f"- sources within 3° outside the |l|,|b| <= 60 grid (their PSF can spill in at the lowest energies): {n_ring}; "
        "at 1.5 GeV the 95% radius is 2.1°, so only the outermost cells are affected, and those are at the ROI edge, not near the halo")
    keep_fit_frac = counts[(keep & fitted)[grid]].sum(0) / counts[fitted[grid]].sum(0)
    say("- fraction of the catalog flux in |b| >= 10 that is modeled (the rest sits inside the extended-source masks, whose pixels are removed from data and model alike): "
        + ", ".join(f"{x:.4f}" for x in keep_fit_frac))
    say(); say("## 6. The template is in every fit"); say()
    import inspect
    src = inspect.getsource(fine.base_stack)
    say(f"- `fine.base_stack` includes 'point_sources': {'point_sources' in src}")
    from totani import selection
    src10 = inspect.getsource(selection.base_stack)
    say(f"- `selection.base_stack` (10° programme) includes 'point_sources': {'point_sources' in src10}")
    OUT.write_text("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
