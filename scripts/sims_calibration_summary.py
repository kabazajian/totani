#!/usr/bin/env python
"""Summary of sims_calibration_2deg.py against the data (oscar_conditional,
2 deg PSF): per configuration, the simulated distributions of the halo at
20.8 GeV, Q^+ (K = 9 and 6) and Q_quad (20.8 GeV; response-valid sum), the
simulated spectra, and the data's position in each null (empirical tail
probability with a 68% Clopper-Pearson range; the independent units are skies).

    python scripts/sims_calibration_summary.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

TAG = "_" + os.environ.get("FINE_TAG", "psf")
DATE = "2026-09-29"
RESP = [3, 4, 5, 6, 7, 8]
E = np.array([1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0])
CONFIGS = ("structured", "flat", "construction")
KINDS = ("null_structured", "null_flat", "halo_structured", "halo_flat")


def emp(obs, null):
    null = np.asarray(null); n = len(null); k = int((null >= obs).sum())
    lo = beta.ppf(0.16, k, n - k + 1) if k > 0 else 0.0; hi = beta.ppf(0.84, k + 1, n - k) if k < n else 1.0
    return f"{k}/{n} ({lo:.2g}–{hi:.2g})"


def main() -> int:
    from halo_spectra_compare import conversion
    conv = conversion(); e2 = lambda c: np.asarray(c) * conv * (E * 1e3) ** 2 / 1e3 * 1e7
    sims = json.loads((ROOT / "results" / "fits" / f"sims_calibration_cells2{TAG}_{DATE}.json").read_text())["results"]
    data = json.loads((ROOT / "results" / "fits" / f"oscar_conditional_cells2{TAG}_{DATE}.json").read_text())["results"]
    D = {r["config"]: r for r in data if tuple(r["row"]) == ("SO_Z6_R30_T150_C5", "gas1_ics1")}
    Dq = None
    qp = ROOT / "results" / "fits" / f"quadrant_calibration_cells2{TAG}_{DATE}.json"
    if qp.exists():
        Dq = {r["config"]: r for r in json.loads(qp.read_text())["results"] if r["kind"] == "data" and tuple(r["row"]) == ("SO_Z6_R30_T150_C5", "gas1_ics1")}
    lines = [f"# Simulated skies at 2° with the PSF, the construction rebuilt per sky: calibration of the reported statistics ({DATE})", "",
             "Truth row SO_Z6_R30_T150_C5 1-gas with the 2-deg companions; 50 skies per truth. Each sky fitted as the data are (fixed configurations; "
             "the construction rebuilt from the sky's own 4.3 GeV pixels on Totani's 10° cells). 'data' = the same row on the real sky. "
             "Tail counts: skies with a statistic ≥ the data's, with the 68% Clopper–Pearson range of the fraction.", ""]
    lines += ["## Distributions (median [16–84%])", "", "| truth | fitted with | halo 20.8 GeV | Q⁺ K=9 | Q⁺ K=6 | Q_quad 20.8 (std) | Q_quad K=6 sum (std) |", "|---|---|---|---|---|---|---|"]
    f = lambda v: f"{np.median(v):,.0f} [{np.percentile(v,16):,.0f}, {np.percentile(v,84):,.0f}]"
    for kind in KINDS:
        S = [r for r in sims if r["kind"] == kind]
        for c in CONFIGS:
            h = [r[c]["counts"][5] for r in S]; q9 = [r[c]["all"]["Q_plus"] for r in S]; q6 = [r[c]["resp"]["Q_plus"] for r in S]
            qq = [r[c]["quad"]["std"][5] for r in S]; qs = [sum(r[c]["quad"]["std"][k] for k in RESP) for r in S]
            lines.append(f"| {kind} ({len(S)}) | {c} | {f(h)} | {f(q9)} | {f(q6)} | {np.median(qq):.1f} [{np.percentile(qq,16):.1f}, {np.percentile(qq,84):.1f}] | {f(qs)} |")
    for c in CONFIGS:
        d = D[c]
        dq = Dq[c] if Dq and c in Dq else None
        lines.append(f"| data | {c} | {d['counts'][5]:,.0f} | {d['all']['Q_plus']:,.0f} | {d['resp']['Q_plus']:,.0f} | "
                     + (f"{dq['std'][5]:.1f} | {sum(dq['std'][k] for k in RESP):,.0f} |" if dq else "— | — |"))
    lines += ["", "## Where the data sit in each simulated scenario (skies ≥ data)", "", "| truth | fitted with | halo 20.8 | Q⁺ K=9 | Q⁺ K=6 |", "|---|---|---|---|---|"]
    for kind in KINDS:
        S = [r for r in sims if r["kind"] == kind]
        for c in CONFIGS:
            lines.append(f"| {kind} | {c} | {emp(D[c]['counts'][5], [r[c]['counts'][5] for r in S])} | {emp(D[c]['all']['Q_plus'], [r[c]['all']['Q_plus'] for r in S])} | "
                         f"{emp(D[c]['resp']['Q_plus'], [r[c]['resp']['Q_plus'] for r in S])} |")
    lines += ["", "## Simulated spectra: median E²dN/dE (×10⁻⁷ GeV cm⁻² s⁻¹) with 16–84% of skies, against the data", "",
              "| truth | fitted with | " + " | ".join(f"{e:g}" for e in E) + " |", "|---|---|" + "---|" * 9]
    for kind in KINDS:
        S = [r for r in sims if r["kind"] == kind]
        for c in CONFIGS:
            M = np.array([e2(r[c]["counts"]) for r in S])
            lines.append(f"| {kind} | {c} | " + " | ".join(f"{np.median(M[:,k]):+.2f} [{np.percentile(M[:,k],16):+.2f}, {np.percentile(M[:,k],84):+.2f}]" for k in range(9)) + " |")
    for c in CONFIGS:
        lines.append(f"| data | {c} | " + " | ".join(f"{x:+.2f}" for x in e2(D[c]["counts"])) + " |")
    out = ROOT / "results" / f"sims_calibration_cells2{TAG}_{DATE}.md"
    out.write_text("\n".join(lines) + "\n"); print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
