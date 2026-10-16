#!/usr/bin/env python
"""Is the finer-cell structured-map excess a binning artifact? Halo-free and
injected skies built and fitted at 2 deg: the truth is the 2-deg fit of the
data with the true row (SO_Z6_R30_T150_C5, 1-gas) and the structured set,
Poisson-sampled on the |b|>=10 cells; each sky is refitted with the same
set (own), with the flat mask, and with the residual construction (rebuilt
at 10 deg from the sky's 4.3 GeV pixels, projected to 2 deg), with and
without a halo injected on Totani's spectrum. 10 skies.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import stage_e
from finer_cells_test import pixel_library, row_templates, evaluate, EDGES, K20
from totani import ranking, simsky
from totani.bubble_totani import bubble_image, fit_all_cells
from totani.cells import bin_to_cells
from totani.likelihood import fit_all_bins, restrict_stack
from totani.stack import cell_template_from_shape
CELL = 2.0; ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
OWN = ("bubble_structured", "erosita_unif_N", "erosita_unif_S", "loopI_haslam_raw", "sun_moon_ic")
SETS = {"own (structured)": OWN, "flat": ("bubble_flat",) + OWN[1:], "residual": ("bubble_totani_pos", "bubble_totani_neg") + OWN[1:]}
lib = pixel_library(CELL); diffuse, _ = row_templates(ROW, CELL, lib)
base = {**diffuse, "isotropic": lib["templates"]["isotropic"], "point_sources": lib["templates"]["point_sources"]}
own_t = {**base, **{n: lib["templates"][n] for n in OWN}}
keep = lib["keep"]
# truth at 2 deg: background-only fit of the data with the own set
fit = fit_all_bins(restrict_stack({"counts": lib["counts"], "n_bins": 9, "templates": own_t}, keep), halo_mask=keep)
names = fit["component_order"]; amps = np.array([[fit["per_bin"][k]["amplitudes"][n] for n in names] for k in range(9)])
truth_cells = np.einsum("kn,nkij->kij", amps, np.stack([own_t[n] for n in names]))
amp = simsky.halo_amplitudes_totani(ROOT / "results" / "fits" / "bubble_totani_injection_2026-09-19.json")
halo_cells = amp[:, None, None] * lib["J"]
# pixel truth at 4.3 GeV for the residual construction (same recipe as simsky.build_truth, own set)
shared = simsky.load_shared(); rowpix = simsky.load_row(*ROW)
target, L, B, single_sr, mask = lib["grid"]
pix_truth = np.zeros(L.shape)
for i, n in enumerate(names):
    pix_truth += amps[2, i] * (np.asarray(rowpix[n], float) if n in rowpix else simsky.structure_pixels(n, EDGES, shared))
from totani.halo import nfw_rho2_on_target
halo_pix = amp[2] * nfw_rho2_on_target(target) * np.asarray(shared["flat_sr"])
disk = np.abs(B) < 10.0; real_pix = np.asarray(shared["counts"])
out = {"cell_deg": CELL, "skies": []}
for i in range(10):
    rng = np.random.default_rng([99, i, 20260924])
    for inj in (False, True):
        mu = truth_cells + (halo_cells if inj else 0)
        cells = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
        pix = np.where(disk, real_pix + (rng.poisson(np.clip(halo_pix, 0, None)) if inj else 0), rng.poisson(np.clip(pix_truth + (halo_pix if inj else 0), 0, None))).astype(float)
        # residual construction from this sky's pixels (his 144 cells), projected to 2 deg
        T = {n: shared[n] for n in simsky.SHARED_NAMES}; T.update(rowpix)
        pm = {"templates": T, "counts": pix, "flat_sr": shared["flat_sr"], "mask": shared["mask"], "L": shared["L"], "B": shared["B"]}
        img = bubble_image(pm, fit_all_cells(pm))
        res = {"bubble_totani_pos": cell_template_from_shape(ranking.RUN, EDGES, img["positive"], cell_deg=CELL),
               "bubble_totani_neg": cell_template_from_shape(ranking.RUN, EDGES, img["negative"], cell_deg=CELL)}
        rec = {"sim": i, "injected": inj}
        for sname, cands in SETS.items():
            t = dict(base)
            for n in cands:
                t[n] = res[n] if n in res else lib["templates"][n]
            r = evaluate(cells, t, lib["J"], keep, lib["quads"], sign_free="bubble_totani_neg" if "bubble_totani_neg" in cands else None)
            rec[sname] = {"halo20": r["halo20"], "lo": r["lo"], "hi": r["hi"], "ts20": r["ts20"], "gate": r["gate"]}
        out["skies"].append(rec)
        print(f"[{i+1}/10 inj={inj}] " + "  ".join(f"{s}: {rec[s]['halo20']:,.0f} (TS {rec[s]['ts20']:.0f}, {'Y' if rec[s]['gate'] else 'n'})" for s in SETS), flush=True)
summ = {}
for inj in (False, True):
    for s in SETS:
        h = np.array([r[s]["halo20"] for r in out["skies"] if r["injected"] == inj]); ts = np.array([r[s]["ts20"] for r in out["skies"] if r["injected"] == inj])
        summ[f"{'injected' if inj else 'null'}__{s}"] = {"halo20_median": float(np.median(h)), "halo20_p16": float(np.percentile(h, 16)), "halo20_p84": float(np.percentile(h, 84)),
                                                          "ts20_median": float(np.median(ts)), "gate_fraction": float(np.mean([r[s]["gate"] for r in out["skies"] if r["injected"] == inj]))}
        print(f"{'injected' if inj else 'null':8s} {s:16s}", json.dumps(summ[f"{'injected' if inj else 'null'}__{s}"]))
out["summary"] = summ
(ROOT / "results" / "fits" / "finer_cells_null_2026-09-24.json").write_text(json.dumps(out, indent=1))
