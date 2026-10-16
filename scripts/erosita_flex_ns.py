#!/usr/bin/env python
"""Items 1 and 3 (Kev, 2026-09-27), PSF-folded 2-deg cells, flat mask:
does the halo's south/north ratio at 20.8 GeV follow the eROSITA lobes, and
does it track the RGB-reconstructed Galactic East of the eROSITA map?

    FINE_TAG=psf python scripts/erosita_flex_ns.py --workers 6

The NFW-rho^2 halo is split into four quadrant pieces (W = l<0, native
eROSITA; E = l>0, reconstructed; N/S at b = 0), each with a signed amplitude,
fitted at 20.8 GeV jointly with the background (flat mask, one normalization;
Loop I Haslam raw; solar IC) under a series of eROSITA models:

  none            no lobes
  unif, bright    Oscar's uniform / brightness lobes, N and S free (the 2-deg default is unif)
  *_outer         the same restricted outside the bubbles
  *_EW            N/S x E/W: four free pieces
  *_flex          N/S x E/W x |b| bands (0-25, 25-40, 40-90 deg): twelve free pieces
  *_mirrorW       the native West copied onto the East (l -> -l), N and S free:
                  the reconstructed East replaced by the measured West
  *_Wonly         only the western lobes (the East given no lobe template)

Per (row, model): halo counts per quadrant, S/N overall and per half, the
quadrant LR (common amplitude vs 4, chi2_3), the fitted lobe counts N and S
and E and W, and the 9-bin background-only -2lnL (data only; the BIC cost is
9 ln 27000 = 92 per extra template).
Data: the 33-row subset and Totani's baseline. Reference: 20 flat-bubble skies
of the true row with a spherical halo on Totani's spectrum (unif lobes true).
"""
from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from stage_e_2deg import TRUE_ROW, truth_cells  # noqa: E402
from totani import fine  # noqa: E402
from totani import likelihood as LK  # noqa: E402
from totani.cells import cell_centers  # noqa: E402
from totani.likelihood import restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
CELL = 2.0
DATE = "2026-09-27"
K20 = 5
BANDS = ((0, 25), (25, 40), (40, 90))
BASE = ("loopI_haslam_raw", "sun_moon_ic")
_W: dict = {}


def flex_path():
    return fine.cache_dir(CELL) / "erosita_flex.npz"


def build_flex():
    """PSF-folded eROSITA pieces: per map (unif, bright) the 12 N/S x E/W x band
    pieces (sums give EW and whole-lobe N/S) and the West-mirrored N/S lobes."""
    p = flex_path()
    if p.exists():
        return p
    from totani import structures
    g = fine._psf_grid(EDGES); L, B = g["L"], g["B"]
    assert np.allclose(L[:, ::-1], -L)
    out = {}
    for m in ("unif", "bright"):
        for h in ("N", "S"):
            a = structures.STRUCTURES[f"erosita_{m}_{h}"](g["target"], EDGES)
            for side, sm in (("W", L < 0), ("E", L >= 0)):
                for lo, hi in BANDS:
                    piece = a * sm * (np.abs(B) >= lo) * (np.abs(B) < hi)
                    out[f"{m}_{h}_{side}_b{lo}"] = fine._cells_psf(piece[None] * g["basis"], g, CELL, EDGES)
            mir = np.where(L < 0, a, a[:, ::-1])
            out[f"{m}_{h}_mirrorW"] = fine._cells_psf(mir[None] * g["basis"], g, CELL, EDGES)
            print(f"built {m} {h}", flush=True)
    np.savez(p, **out)
    return p


def models(F):
    """{model name: {template name: cells cube}} from the flex pieces."""
    s = lambda keys: sum(F[k] for k in keys)
    out = {"none": {}}
    lib = _W["lib"]
    for m in ("unif", "bright"):
        out[m] = {f"{m}_{h}": s([f"{m}_{h}_{sd}_b{lo}" for sd in "WE" for lo, _ in BANDS]) for h in "NS"}
        out[f"{m}_outer"] = {f"{m}_outer_{h}": lib["templates"][f"erosita_{m}_outer_{h}"] for h in "NS"}
        out[f"{m}_EW"] = {f"{m}_{h}_{sd}": s([f"{m}_{h}_{sd}_b{lo}" for lo, _ in BANDS]) for h in "NS" for sd in "WE"}
        out[f"{m}_flex"] = {f"{m}_{h}_{sd}_b{lo}": F[f"{m}_{h}_{sd}_b{lo}"] for h in "NS" for sd in "WE" for lo, _ in BANDS}
        out[f"{m}_mirrorW"] = {f"{m}_{h}_mirrorW": F[f"{m}_{h}_mirrorW"] for h in "NS"}
        out[f"{m}_Wonly"] = {f"{m}_{h}_W": s([f"{m}_{h}_W_b{lo}" for lo, _ in BANDS]) for h in "NS"}
    return out


MODEL_ORDER = ["none", "unif", "bright", "unif_outer", "bright_outer", "unif_EW", "bright_EW",
               "unif_flex", "bright_flex", "unif_mirrorW", "bright_mirrorW", "unif_Wonly", "bright_Wonly"]


def _init():
    _W["lib"] = fine.load_library(CELL)
    F = dict(np.load(flex_path()))
    # sanity: the rebuilt whole lobes equal the library's (same pixels, same PSF)
    for m in ("unif", "bright"):
        for h in "NS":
            rebuilt = sum(F[f"{m}_{h}_{sd}_b{lo}"] for sd in "WE" for lo, _ in BANDS)
            ref = _W["lib"]["templates"][f"erosita_{m}_{h}"]
            assert np.allclose(rebuilt, ref, rtol=1e-5, atol=1e-6 * np.abs(ref).max()), (m, h)
    _W["models"] = models(F)
    Lc, Bc = cell_centers(cell_deg=CELL); _W["N"] = Bc > 0; _W["E"] = Lc > 0


def joint_full(T_bg, n, sigs):
    """-2lnL; background >= 0, signals signed. Returns (v, background amps, signal amps)."""
    m = T_bg.shape[0]; s = len(sigs)
    sb = np.maximum(np.abs(T_bg).reshape(m, -1).sum(1), 1e-300); ss = np.maximum(np.array([np.abs(t).sum() for t in sigs]), 1e-300)
    Tb = (T_bg / sb.reshape((m,) + (1,) * (T_bg.ndim - 1))).reshape(m, -1); Ts = np.stack([t / c for t, c in zip(sigs, ss)]).reshape(s, -1); nf = n.reshape(-1)
    def obj(x):
        mu = x[:m] @ Tb + x[m:] @ Ts; v, dmu = LK._poisson_val_grad_mu(mu, nf)
        return v, np.concatenate([Tb @ dmu, Ts @ dmu])
    x0 = np.concatenate([np.full(m, nf.sum() / (m + s)), np.zeros(s)])
    bounds = [(0.0, 1e14)] * m + [(-1e14, 1e14)] * s
    r = minimize(obj, x0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 8000, "ftol": 1e-15, "gtol": 1e-8})
    for _ in range(6):
        if r.status == 0:
            break
        r = minimize(obj, r.x, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 8000, "ftol": 1e-15, "gtol": 1e-8})
    return float(r.fun), r.x[:m] / sb, r.x[m:] / ss


def stack(counts, row_key, lobes):
    lib = _W["lib"]; row = fine.load_row(CELL, *row_key)
    t = {"isotropic": lib["templates"]["isotropic"], "point_sources": lib["templates"]["point_sources"], **row["diffuse"],
         "bubble_flat": lib["templates"]["bubble_flat"], **{n: lib["templates"][n] for n in BASE}, **lobes}
    return {"counts": counts, "n_bins": 9, "templates": t}


def analyse(counts, row_key, with_h0):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]; Nm, Em = _W["N"], _W["E"]
    out = {}
    for mname in MODEL_ORDER:
        lobes = _W["models"][mname]
        st = stack(counts, row_key, lobes)
        sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
        T = np.stack([sr["templates"][nm][K20] for nm in names]); n = sr["counts"][K20]; h = Jr[K20]
        quads = {"WN": h * (~Em & Nm), "WS": h * (~Em & ~Nm), "EN": h * (Em & Nm), "ES": h * (Em & ~Nm)}
        v4, bg, a4 = joint_full(T, n, list(quads.values()))
        v1, _, a1 = joint_full(T, n, [h])
        hc = {q: float(a * quads[q][keep].sum()) for q, a in zip(quads, a4)}
        lc = {}
        for nm, amp in zip(names, bg):
            if nm in lobes:
                c = amp * sr["templates"][nm][K20]
                for lab, msk in (("N", Nm), ("S", ~Nm), ("E", Em), ("W", ~Em)):
                    lc[lab] = lc.get(lab, 0.0) + float(c[msk & keep].sum())
        rec = {"halo": hc, "halo_total_common": float(a1[0] * h[keep].sum()), "v4": v4, "v1": v1, "lobes": lc}
        if with_h0:
            rec["h0"] = float(fine.bg_neg2lnl(st, keep).sum()); rec["n_lobe_templates"] = len(lobes)
        out[mname] = rec
    return out


def task(t):
    kind, i, mu, row_key = t
    lib = _W["lib"]; keep = lib["keep"]
    if mu is None:
        counts = lib["counts"]
    else:
        rng = np.random.default_rng([7, i, 20260927])
        counts = np.where(keep[None], rng.poisson(np.clip(mu, 0, None)), lib["counts"]).astype(float)
    return {"kind": kind, "sim": i, "row": list(row_key), "models": analyse(counts, row_key, with_h0=(kind == "data"))}


def ratios(h):
    r = lambda s, n: s / n if abs(n) > 50 else float("nan")
    return (r(h["WS"] + h["ES"], h["WN"] + h["EN"]), r(h["WS"], h["WN"]), r(h["ES"], h["EN"]))


def render(res):
    from scipy.stats import chi2
    lines = [f"# Does the halo's S/N at 20.8 GeV follow the eROSITA lobes? PSF-folded 2° cells, flat mask ({DATE})", "",
             "Halo fitted as four quadrant pieces (W = l<0, native eROSITA; E = l>0, reconstructed). Medians over rows (data: 33-row subset + Totani's baseline; "
             "reference: 20 flat-bubble skies of the true row with a spherical halo, unif lobes true). "
             "Δ(−2lnL) H0 is the 9-bin background-only fit relative to the unif default (negative = better; 92 per extra template is the BIC cost). "
             "Quadrant LR: common amplitude vs 4 quadrants at 20.8 GeV, χ²₃. Lobe S/N and E/W: fitted lobe counts at 20.8 GeV.", ""]
    D = [r for r in res if r["kind"] == "data"]; S = [r for r in res if r["kind"] == "halo_flat"]
    lines += ["## Data", "", "| eROSITA model | lobe templates | Δ(−2lnL) H0 vs unif | halo total | S/N | S/N West | S/N East | quadrant LR p | lobe counts N / S | lobe E/W |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for m in MODEL_ORDER:
        dh0 = [r["models"][m]["h0"] - r["models"]["unif"]["h0"] for r in D]
        H = [r["models"][m]["halo"] for r in D]; rt = np.array([ratios(h) for h in H])
        tot = [sum(h.values()) for h in H]
        lr = [chi2.sf(max(r["models"][m]["v1"] - r["models"][m]["v4"], 0), 3) for r in D]
        lo = [r["models"][m]["lobes"] for r in D]
        lN = np.median([x.get("N", 0) for x in lo]); lS = np.median([x.get("S", 0) for x in lo])
        lEW = np.nanmedian([x.get("E", 0) / x["W"] if abs(x.get("W", 0)) > 1 else np.nan for x in lo])
        lines.append(f"| {m} | {D[0]['models'][m]['n_lobe_templates']} | {np.median(dh0):+.0f} | {np.median(tot):,.0f} | {np.nanmedian(rt[:,0]):.2f} | "
                     f"{np.nanmedian(rt[:,1]):.2f} | {np.nanmedian(rt[:,2]):.2f} | {np.median(lr):.2g} | {lN:,.0f} / {lS:,.0f} | {lEW:.2f} |")
    lines += ["", "## Reference skies (spherical halo injected)", "", "| eROSITA model | halo total | S/N [16–84%] | S/N West [16–84%] | S/N East [16–84%] |", "|---|---|---|---|---|"]
    for m in MODEL_ORDER:
        H = [r["models"][m]["halo"] for r in S]; rt = np.array([ratios(h) for h in H]); tot = [sum(h.values()) for h in H]
        f = lambda c: f"{np.nanmedian(rt[:,c]):.2f} [{np.nanpercentile(rt[:,c],16):.2f}, {np.nanpercentile(rt[:,c],84):.2f}]"
        lines.append(f"| {m} | {np.median(tot):,.0f} | {f(0)} | {f(1)} | {f(2)} |")
    base = [r for r in D if r["row"][0] == "SL_Z6_R30_T150_C2"]
    if base:
        lines += ["", "## Totani's baseline row, quadrant halo counts (WN, WS, EN, ES)", "", "| eROSITA model | WN | WS | EN | ES |", "|---|---|---|---|---|"]
        for m in MODEL_ORDER:
            h = base[0]["models"][m]["halo"]
            lines.append(f"| {m} | {h['WN']:,.0f} | {h['WS']:,.0f} | {h['EN']:,.0f} | {h['ES']:,.0f} |")
    (ROOT / "results" / f"erosita_flex_ns{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); ap.add_argument("--render", action="store_true")
    a = ap.parse_args()
    js = ROOT / "results" / "fits" / f"erosita_flex_ns{fine.tag()}_{DATE}.json"
    if a.render:
        render(json.loads(js.read_text())["results"]); return 0
    build_flex()
    data_rows = [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()] + [("a2020", "-"), ("SL_Z6_R30_T150_C2", "gas1_ics1")]
    tasks = [("data", 0, None, r) for r in data_rows]
    mu = truth_cells("bubble_flat", True)
    tasks += [("halo_flat", i, mu, TRUE_ROW) for i in range(20)]
    res = []
    with Pool(a.workers, initializer=_init) as pool:
        for k, r in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(r); print(f"[{k}/{len(tasks)}] {r['kind']} {r['row']}", flush=True)
    js.write_text(json.dumps({"date": DATE, "results": res}, indent=1))
    render(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
