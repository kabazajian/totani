#!/usr/bin/env python
"""Item 2 (Kev, 2026-09-27), second half: the HESTIA J-maps as halo templates
on PSF-folded 2-deg cells at 20.8 GeV, flat mask and the 2-deg companions.

    FINE_TAG=psf python scripts/hestia_halo_fits.py --workers 6

Data: the HESTIA Milky-Way analogues of Muru et al., PRL 135, 161005 (2025);
maps from scripts/hestia_roi60_maps.py (6 halos x 36 solar azimuths, and a
spherically randomized copy of each halo as the no-asymmetry control).

Per view: the template's own south/north count ratio over the kept cells,
relative to NFW's (A = 1: as symmetric as NFW); then on the data, for five
rows, (i) Delta(-2lnL) of the single-amplitude HESTIA halo against the
single-amplitude NFW halo (negative = HESTIA fits better), and (ii) with the
template split at b = 0, the fitted amplitude ratio a_S/a_N (1 = the
template's own asymmetry accounts for the data's; NFW gives ~1.4-1.6).
The sign of b in a simulation is arbitrary, so each view is also scored with
its latitude flipped (a halo brighter in its north is the same halo seen
with b -> -b).
"""
from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from ns_halo_spectra import joint  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from stage_e_2deg import TRUE_ROW, stack_for  # noqa: E402
from totani import fine, psf  # noqa: E402
from totani.cells import bin_to_cells, cell_centers  # noqa: E402
from totani.likelihood import restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
CELL = 2.0
DATE = "2026-09-27"
K20 = 5
MAPS = Path.home() / "Data2" / "hestia_muru" / "totani_roi60"
HALOS = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), TRUE_ROW, ("a2020", "-"),
        (stage_e.fit_codes()[0], "gas1_ics1"), (stage_e.fit_codes()[0], "gas5_ics3")]
_W: dict = {}


def _init():
    lib = fine.load_library(CELL); _W["lib"] = lib
    g = fine._psf_grid(EDGES); _W["g"] = {k: g[k] for k in ("L", "B", "basis", "mask")}
    _, Bc = cell_centers(cell_deg=CELL); _W["north"] = Bc > 0
    _W["stacks"] = {}
    for r in ROWS:
        st = stack_for(lib, fine.load_row(CELL, *r), "bubble_flat", lib["counts"])
        sr, Jr = restrict_stack(st, lib["keep"], lib["J"]); names = sorted(sr["templates"])
        _W["stacks"][r] = (np.stack([sr["templates"][n][K20] for n in names]), sr["counts"][K20])


def cells_k20(j025):
    """0.25-deg ROI map (rows b from -60, columns l from +60 down) -> PSF-folded 2-deg cells at 20.8 GeV."""
    g = _W["g"]
    pix = np.repeat(np.repeat(j025.astype(float), 2, axis=0), 2, axis=1) * g["basis"][K20]
    conv = psf.convolve(pix, g["B"], EDGES[K20], EDGES[K20 + 1])
    return bin_to_cells(conv, g["L"], g["B"], cell_deg=CELL, mask=g["mask"])


def score(h):
    keep = _W["lib"]["keep"]; north = _W["north"]
    h = h * keep
    out = {"template_SN": float(h[~north].sum() / h[north].sum())}
    for r, (T, n) in _W["stacks"].items():
        v1, _ = joint(T, n, [h])
        _, a = joint(T, n, [h * north, h * ~north])
        out["|".join(r)] = {"v": v1, "aS_over_aN": float(a[1] / a[0]) if a[0] != 0 else float("nan")}
    return out


def task(t):
    halo, i, kind = t
    m = np.load(MAPS / f"{halo}.npz")["jmap" if kind == "hestia" else "jmap_spherical"][i]
    h = cells_k20(m)
    return {"halo": halo, "view": i, "kind": kind, "as_is": score(h), "flipped": score(cells_k20(m[::-1]))}


def nfw_ref():
    _init() if not _W else None
    return score(_W["lib"]["J"][K20])


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    tasks = [(h, i, k) for h in HALOS for i in range(36) for k in ("hestia", "spherical")]
    with Pool(a.workers, initializer=_init) as pool:
        ref = pool.apply(nfw_ref)
        res = []
        for k, r in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(r)
            if k % 36 == 0:
                print(f"[{k}/{len(tasks)}]", flush=True)
    (ROOT / "results" / "fits" / f"hestia_halo_fits{fine.tag()}_{DATE}.json").write_text(json.dumps({"date": DATE, "nfw": ref, "results": res}, indent=1))
    render(ref, res)
    return 0


def render(ref, res):
    keys = ["|".join(r) for r in ROWS]
    lines = [f"# HESTIA halos (Muru et al. 2025) as the 20.8 GeV halo template, PSF-folded 2° cells, flat mask ({DATE})", "",
             f"NFW reference: template S/N {ref['template_SN']:.3f}; fitted a_S/a_N on the rows: "
             + ", ".join(f"{k.split('|')[0]} {ref[k]['aS_over_aN']:.2f}" for k in keys), "",
             "A = template S/N relative to NFW's (both orientations of b). Δ = −2lnL(HESTIA) − (−2lnL(NFW)), single amplitude; "
             "a_S/a_N = fitted amplitude ratio with the template split at b = 0 (1 = the template's asymmetry accounts for the data's). Medians [min, max] over 36 views.", "",
             "| halo | kind | A (as is) | Δ best orientation, Totani baseline row | a_S/a_N best orientation, baseline | a_S/a_N, median over the 5 rows |",
             "|---|---|---|---|---|---|"]
    b0 = keys[0]
    for halo in HALOS:
        for kind in ("hestia", "spherical"):
            S = [r for r in res if r["halo"] == halo and r["kind"] == kind]
            A = np.array([r["as_is"]["template_SN"] / ref["template_SN"] for r in S])
            best = [min((r["as_is"], r["flipped"]), key=lambda x: x[b0]["v"]) for r in S]
            d = np.array([x[b0]["v"] - ref[b0]["v"] for x in best]); q = np.array([x[b0]["aS_over_aN"] for x in best])
            q5 = np.array([np.median([x[k]["aS_over_aN"] for k in keys]) for x in best])
            f = lambda v, p=2: f"{np.median(v):.{p}f} [{v.min():.{p}f}, {v.max():.{p}f}]"
            lines.append(f"| {halo} | {kind} | {f(A, 3)} | {f(d, 1)} | {f(q)} | {f(q5)} |")
    allA = np.array([max(r["as_is"]["template_SN"], r["flipped"]["template_SN"]) / ref["template_SN"] for r in res if r["kind"] == "hestia"])
    lines += ["", f"Largest template asymmetry over all 216 views (either orientation): {allA.max():.3f} × NFW's; views with ≥ 1.3: {int((allA >= 1.3).sum())}; ≥ 1.1: {int((allA >= 1.1).sum())}.", ""]
    (ROOT / "results" / f"hestia_halo_fits{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
