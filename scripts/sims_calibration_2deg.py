#!/usr/bin/env python
"""Oscar's note of 2026-09-29, items 2 and 'Reuse the simulations to calibrate
the procedure we actually report': simulated skies at 2 deg with the PSF, the
residual construction REBUILT for every sky from that sky's own 4.3 GeV pixels,
each sky fitted exactly as the data are in oscar_conditional.py (fixed
configurations: structured map, flat mask, construction; the 2-deg companions;
signed profile fits, n_refine = 6). Per sky and configuration: the halo counts
per bin (the simulated spectra), per-bin TS and q^+, Q^+ over K = 9 and over
the response-valid bins (K = 6). The independent units are skies.

    FINE_TAG=psf python scripts/sims_calibration_2deg.py --skies 50 --workers 5

Truths (tested scenarios, not exhaustive coverage): the true row
SO_Z6_R30_T150_C5 1-gas with the 2-deg companions and EITHER the structured map
OR the flat mask, background-only fits to the data at 2 deg (PSF-folded);
with no halo and with a spherical NFW halo on Totani's spectrum. Cells: Poisson
draws of the PSF-folded truth on the |b| >= 10 cells, real counts in the disk;
the 4.3 GeV bin is drawn per PIXEL from the PSF-folded pixel truth (disk pixels
= real counts + a draw of the halo), binned to the 2-deg cells, and the same
pixel map feeds Totani's 144-cell construction (unconvolved, as he builds it).
Not repeated per sky: the structure selection (the configurations are fixed,
as in the conditional comparison), so these calibrate the conditional test.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import oscar_conditional as OC  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import fine, psf, ranking, simsky, stats  # noqa: E402
from totani.cells import bin_to_cells  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402
from totani.stack import cell_template_from_shape  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
CELL = 2.0
DATE = "2026-09-29"
TRUE_ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
K_BUILD = 2
_W: dict = {}


def _init():
    OC._init(CELL)
    _W["lib"] = OC._W["lib"]; _W["shared"] = simsky.load_shared(); _W["rowpix"] = simsky.load_row(*TRUE_ROW)
    g = fine._psf_grid(EDGES); _W["g"] = {k: g[k] for k in ("L", "B", "mask")}
    z = np.load(fine.cache_dir(CELL) / "halo_quadrants.npz")
    _W["quads"] = {d: [z[k] for k in z.files if k.startswith(d + "__")] for d in ("std", "rot")}


def truth(bubble, with_halo):
    """PSF-folded cell truth and PSF-folded 4.3 GeV pixel truth (halo separately)."""
    from totani.halo import nfw_rho2_on_target
    _init()
    lib = _W["lib"]; keep = lib["keep"]; shared = _W["shared"]; rowpix = _W["rowpix"]
    st = fine.base_stack(lib, fine.load_row(CELL, *TRUE_ROW))
    st["templates"][bubble] = lib["templates"][bubble]
    for n in OC.COMP:
        st["templates"][n] = lib["templates"][n]
    fit = fit_all_bins(restrict_stack(st, keep), halo_mask=keep)
    names = fit["component_order"]
    amps = np.array([[fit["per_bin"][k]["amplitudes"][n] for n in names] for k in range(9)])
    cells = np.einsum("kn,nkij->kij", amps, np.stack([st["templates"][n] for n in names]))
    pix = np.zeros(shared["L"].shape)
    for i, n in enumerate(names):
        p = rowpix[n] if n in rowpix else simsky.structure_pixels(n, EDGES, shared)
        pix += amps[K_BUILD, i] * np.asarray(p, float)
    g = _W["g"]
    pix = psf.convolve(pix, g["B"], EDGES[K_BUILD], EDGES[K_BUILD + 1])
    out = {"cells": cells, "pix2": pix, "halo_cells": np.zeros_like(cells), "halo_pix2": np.zeros_like(pix)}
    if with_halo:
        amp = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")
        target, _, _, _ = ranking._grid_and_basis(EDGES)
        out["halo_cells"] = amp[:, None, None] * lib["J"]
        out["halo_pix2"] = psf.convolve(amp[K_BUILD] * nfw_rho2_on_target(target) * np.asarray(shared["flat_sr"], float),
                                        g["B"], EDGES[K_BUILD], EDGES[K_BUILD + 1])
    return out


def evaluate_sky(counts, residual):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    row = fine.load_row(CELL, *TRUE_ROW)
    out = {}
    for config, cand in OC.CONFIGS.items():
        st = fine.base_stack(lib, row, counts)
        for n in OC.COMP:
            st["templates"][n] = lib["templates"][n]
        row_used = {"diffuse": row["diffuse"], "residual": residual}
        st, _ = fine.add_candidate(st, cand, lib, row_used, keep)
        sr, Jr = restrict_stack(st, keep, J)
        bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
        pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=False, n_refine=OC.N_REFINE)["per_bin"]
        roi = [float(Jr[k][keep].sum()) for k in range(9)]
        ts = [bg[k]["neg2lnl"] - pf[k]["neg2lnl"] for k in range(9)]; a = [pf[k]["signal_amplitude"] for k in range(9)]
        from ns_halo_spectra import joint
        names = sorted(sr["templates"]); qd = {"std": [], "rot": []}
        for k in range(9):
            T = np.stack([sr["templates"][n][k] for n in names]); nk = sr["counts"][k]
            v1, _ = joint(T, nk, [Jr[k]])
            for d in qd:
                v4, _ = joint(T, nk, [q[k] * keep for q in _W["quads"][d]]); qd[d].append(float(v1 - v4))
        out[config] = {"counts": [a[k] * roi[k] for k in range(9)], "ts": ts, "quad": qd, "v0": [b["neg2lnl"] for b in bg], "v1": [p["neg2lnl"] for p in pf],
                       "all": stats.summarize(ts, a), "resp": stats.summarize(ts, a, OC.RESP_BINS)}
    return out


TRUTH_DIR = fine.cache_dir(CELL) / "sims_truths"
KINDS = {"null_structured": ("bubble_structured", False), "null_flat": ("bubble_flat", False),
         "halo_structured": ("bubble_structured", True), "halo_flat": ("bubble_flat", True)}


def load_truth(kind):
    if kind not in _W.setdefault("truths", {}):
        z = np.load(TRUTH_DIR / f"{kind}.npz", mmap_mode="r"); _W["truths"][kind] = {k: z[k] for k in z.files}
    return _W["truths"][kind]


def task(t):
    kind, i = t
    tr = load_truth(kind)
    lib = _W["lib"]; keep = lib["keep"]; shared = _W["shared"]; g = _W["g"]
    rng = np.random.default_rng([17, i, 20260929, ("null_structured", "null_flat", "halo_structured", "halo_flat").index(kind)])
    mu = tr["cells"] + tr["halo_cells"]
    cells = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    disk = simsky.disk_pixels(shared); real_pix = np.asarray(shared["counts"], float)
    pix = np.where(disk, real_pix + rng.poisson(np.clip(tr["halo_pix2"], 0, None)),
                   rng.poisson(np.clip(tr["pix2"] + tr["halo_pix2"], 0, None))).astype(float)
    b2 = bin_to_cells(pix, g["L"], g["B"], cell_deg=CELL, mask=g["mask"])
    cells[K_BUILD] = np.where(keep, b2, lib["counts"][K_BUILD])
    # project Totani's images onto 2-deg cells (as the data's residual pair is), not his 10-deg cells
    from totani.bubble_totani import bubble_image, fit_all_cells
    T = {n: shared[n] for n in simsky.SHARED_NAMES}; T.update(_W["rowpix"])
    pm = {"templates": T, "counts": pix, "flat_sr": shared["flat_sr"], "mask": shared["mask"], "L": shared["L"], "B": shared["B"]}
    img = bubble_image(pm, fit_all_cells(pm))
    residual = {"bubble_totani_pos": cell_template_from_shape(ranking.RUN, EDGES, img["positive"], cell_deg=CELL),
                "bubble_totani_neg": cell_template_from_shape(ranking.RUN, EDGES, img["negative"], cell_deg=CELL)}
    return {"kind": kind, "sim": i, **evaluate_sky(cells, residual)}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--skies", type=int, default=50); ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    assert fine.tag().startswith("_psf")
    t0 = time.time()
    TRUTH_DIR.mkdir(exist_ok=True)
    for kind, (bub, h) in KINDS.items():
        if not (TRUTH_DIR / f"{kind}.npz").exists():
            np.savez(TRUTH_DIR / f"{kind}.npz", **truth(bub, h))
    tasks = [(kind, i) for kind in KINDS for i in range(a.skies)]
    print(f"{len(tasks)} skies ({time.time()-t0:.0f}s)", flush=True)
    js = ROOT / "results" / "fits" / f"sims_calibration_cells2{fine.tag()}_{DATE}.json"
    res = []
    with Pool(a.workers, initializer=_init) as pool:
        for k, r in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(r)
            print(f"[{k}/{len(tasks)}] {r['kind']} {r['sim']}: " + "  ".join(f"{c} {r[c]['counts'][5]:,.0f} Q+ {r[c]['all']['Q_plus']:.0f}" for c in OC.CONFIGS)
                  + f" ({time.time()-t0:.0f}s)", flush=True)
            if k % 20 == 0 or k == len(tasks):
                js.write_text(json.dumps({"date": DATE, "results": res}, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
