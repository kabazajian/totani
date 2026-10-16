#!/usr/bin/env python
"""Follow-up to erosita_flex_ns.py: under the flexible eROSITA lobes, does the
halo keep its significance and its peaked spectrum, and where do the lobes
put the counts? PSF-folded 2-deg cells, flat mask.

    FINE_TAG=psf python scripts/erosita_flex_spectra.py --workers 6

Per data row (33-row subset + Totani's baseline) and lobe model (unif,
unif_flex, bright, bright_flex): per bin the N and S halo amplitudes (signed,
joint), the TS of a common halo (vs no halo), the shape test S = r N in every
bin against free N, S (chi2_8), and at 20.8 GeV the fitted counts of every
lobe piece.
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

import erosita_flex_ns as EF  # noqa: E402
import stage_e  # noqa: E402
from ns_halo_spectra import hemi_conversion, joint  # noqa: E402
from totani import fine  # noqa: E402
from totani.likelihood import restrict_stack  # noqa: E402

MODELS = ("unif", "unif_flex", "bright", "bright_flex")
DATE = EF.DATE
E = np.array([1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0])


def analyse(row_key):
    W = EF._W; lib = W["lib"]; keep, J = lib["keep"], lib["J"]; Nm, Em = W["N"], W["E"]
    out = {"row": list(row_key)}
    for m in MODELS:
        lobes = W["models"][m]
        sr, Jr = restrict_stack(EF.stack(lib["counts"], row_key, lobes), keep, J); names = sorted(sr["templates"])
        rec = {"N": [], "S": [], "ts": [], "v_free": []}
        for k in range(9):
            T = np.stack([sr["templates"][nm][k] for nm in names]); n = sr["counts"][k]; h = Jr[k]
            v2, a = joint(T, n, [h * Nm, h * ~Nm]); v1, _ = joint(T, n, [h]); v0, _ = joint(T, n, [h], fixed={0: 0.0})
            rec["N"].append(float(a[0] * (h * Nm)[keep].sum())); rec["S"].append(float(a[1] * (h * ~Nm)[keep].sum()))
            rec["ts"].append(float(v0 - v1)); rec["v_free"].append(v2)
            if k == EF.K20:
                _, bg, _ = EF.joint_full(T, n, [h])
                rec["pieces20"] = {nm: float((amp * sr["templates"][nm][k])[keep].sum()) for nm, amp in zip(names, bg) if nm in lobes}
        def v_ratio(r):
            return sum(joint(np.stack([sr["templates"][nm][k] for nm in names]), sr["counts"][k], [Jr[k] * Nm + r * Jr[k] * ~Nm])[0] for k in range(9))
        rs = np.linspace(0.6, 2.2, 17); vr = [v_ratio(r) for r in rs]
        rec["shape_delta"] = float(min(vr) - sum(rec["v_free"])); rec["best_ratio"] = float(rs[int(np.argmin(vr))])
        out[m] = rec
    return out


def main() -> int:
    from scipy.stats import chi2
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); a = ap.parse_args()
    rows = [(c, v) for v in ("gas1_ics1", "gas5_ics3") for c in stage_e.fit_codes()] + [("a2020", "-"), ("SL_Z6_R30_T150_C2", "gas1_ics1")]
    with Pool(a.workers, initializer=EF._init) as pool:
        res = pool.map(analyse, rows)
    conv = hemi_conversion()
    e2 = lambda c, h: np.asarray(c) * conv[h] * (E * 1e3) ** 2 / 1e3 * 1e7
    (ROOT / "results" / "fits" / f"erosita_flex_spectra{fine.tag()}_{DATE}.json").write_text(json.dumps({"date": DATE, "results": res, "conversion": {k: v.tolist() for k, v in conv.items()}}, indent=1))
    lines = [f"# Halo significance and spectrum under flexible eROSITA lobes, PSF-folded 2° cells, flat mask ({DATE})", "",
             "Medians over 34 data rows. TS20 = common halo vs none at 20.8 GeV. Shape Δ: S = r·N in all 9 bins vs free (χ²₈). "
             "E²dN/dE per hemisphere ×10⁻⁷ GeV cm⁻² s⁻¹, bins 1.5–100 GeV.", "",
             "| lobes | TS20 | halo20 N / S | S/N at 20.8 | shape Δ (p) | best r | N E²dN/dE | S E²dN/dE |", "|---|---|---|---|---|---|---|---|"]
    for m in MODELS:
        R = [r[m] for r in res]
        N = np.median([x["N"] for x in R], 0); S = np.median([x["S"] for x in R], 0)
        sd = np.median([x["shape_delta"] for x in R])
        lines.append(f"| {m} | {np.median([x['ts'][5] for x in R]):.0f} | {N[5]:,.0f} / {S[5]:,.0f} | {np.median([x['S'][5]/x['N'][5] for x in R]):.2f} | "
                     f"{sd:.1f} ({chi2.sf(max(sd,0),8):.2g}) | {np.median([x['best_ratio'] for x in R]):.2f} | "
                     + " ".join(f"{v:+.2f}" for v in e2(N, 'N')) + " | " + " ".join(f"{v:+.2f}" for v in e2(S, 'S')) + " |")
    base = next(r for r in res if r["row"][0] == "SL_Z6_R30_T150_C2")
    for m in ("unif_flex", "bright_flex"):
        lines += ["", f"## Lobe pieces at 20.8 GeV, Totani's baseline row, {m} (fitted counts)", "", "| piece | counts |", "|---|---|"]
        lines += [f"| {k} | {v:,.0f} |" for k, v in sorted(base[m]["pieces20"].items())]
    (ROOT / "results" / f"erosita_flex_spectra{fine.tag()}_{DATE}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
