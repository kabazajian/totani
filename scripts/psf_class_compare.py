#!/usr/bin/env python
"""CLEAN-table (FINE_TAG=psf) against ULTRACLEAN-table (FINE_TAG=psfuc) results,
for every PSF-folded number the paper quotes (2026-09-30).

    python scripts/psf_class_compare.py

Writes results/psf_class_compare_2026-09-30.md: per analysis, the old and new
values side by side. Missing files (a step still running) are reported as such.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "results" / "fits"
R = ROOT / "results"
OUT = R / "psf_class_compare_2026-09-30.md"


def load(name):
    p = F / name
    return json.loads(p.read_text()) if p.exists() else None


def ratios(h):
    r = lambda s, n: s / n if abs(n) > 50 else np.nan
    return r(h["WS"] + h["ES"], h["WN"] + h["EN"]), r(h["WS"], h["WN"]), r(h["ES"], h["EN"])


def main() -> int:
    L = ["# CLEAN-table vs ULTRACLEAN-table PSF: every PSF-folded headline (2026-09-30)", "",
         "‘old’ = FINE_TAG=psf (gce-fisher CLEAN gtpsf table); ‘new’ = FINE_TAG=psfuc (ULTRACLEAN, this run's livetime cube).", ""]

    # 1. fixed-template conditional table, 2 and 1 deg
    for cell in (2, 1):
        a = load(f"oscar_conditional_cells{cell}_psf_2026-09-29.json"); b = load(f"oscar_conditional_cells{cell}_psfuc_2026-09-29.json")
        L += [f"## Fixed-template table, {cell}° (halo 20.8 GeV [68%], Z₂₀₁₈ K=6)", ""]
        if not (a and b):
            L += ["(not yet available)", ""]; continue
        A = {(tuple(r["row"]), r["config"]): r for r in a["results"]}
        L += ["| row | bubbles | old halo | new halo | Δ% | old Z | new Z |", "|---|---|---|---|---|---|---|"]
        for r in b["results"]:
            o = A[(tuple(r["row"]), r["config"])]
            L.append(f"| {r['row'][0]} {r['row'][1]} | {r['config']} | {o['counts'][5]:,.0f} [{o['lo'][5]:,.0f}, {o['hi'][5]:,.0f}] | "
                     f"{r['counts'][5]:,.0f} [{r['lo'][5]:,.0f}, {r['hi'][5]:,.0f}] | {100*(r['counts'][5]-o['counts'][5])/o['counts'][5]:+.1f} | "
                     f"{o['resp']['Z_2018']:.1f} | {r['resp']['Z_2018']:.1f} |")
        # conditional differences
        L += ["", "| row | flat − structured H1 (all / K=6) old | new | construction − flat H1 (all / K=6) old | new |", "|---|---|---|---|---|"]
        for row in sorted({tuple(r["row"]) for r in b["results"]}):
            def d(src, x, y):
                g = {r["config"]: r for r in src["results"] if tuple(r["row"]) == row}
                f = lambda bins: sum(g[x]["v1"][i] for i in bins) - sum(g[y]["v1"][i] for i in bins)
                return f"{f(range(9)):+.0f} / {f([3,4,5,6,7,8]):+.0f}"
            L.append(f"| {row[0]} {row[1]} | {d(a,'flat','structured')} | {d(b,'flat','structured')} | {d(a,'construction','flat')} | {d(b,'construction','flat')} |")
        L.append("")

    # 2. selections
    L += ["## Selections (construction excluded): selected set; 2012 halo20 median; Z median; ⟨σv⟩ median", ""]
    L += ["| run | split | old | new |", "|---|---|---|---|"]
    for cell in (2, 1):
        for bins in ("all", "resp"):
            a = load(f"stage_select_cells{cell}_psf_{bins}_2026-09-29.json"); b = load(f"stage_select_cells{cell}_psfuc_{bins}_2026-09-29.json")
            for v in ("gas1_ics1", "gas5_ics3"):
                def summ(x):
                    if not x or v not in x["splits"] or "no_residual" not in x["splits"][v]:
                        return "—"
                    d = x["splits"][v]["no_residual"]; s = [y for y in d["scored"] if not y["model"].startswith("G21")]
                    own = sum(1 for o in d["own_paths"] if not o["model"].startswith("G21") and "bubble_structured" in o["chosen"])
                    return (f"{' → '.join(c.replace('bubble_','').replace('erosita_','e:') for c in d['selected'])}; "
                            f"{np.median([y['halo20'] for y in s]):,.0f}; Z {np.median([y['Q']['Z_2018'] for y in s]):.1f}; "
                            f"{np.median([y['sigmav']['sigmav_best'] for y in s])*1e25:+.1f}; own-path structured {own}/{len([o for o in d['own_paths'] if not o['model'].startswith('G21')])}")
                L.append(f"| {cell}°, {bins} | {v} | {summ(a)} | {summ(b)} |")
    L.append("")

    # 3. stage_e_2deg vote
    L += ["## Simulated vote at 2° (flat wins, median over skies; halo under chosen template)", ""]
    for tag in ("psf", "psfuc"):
        p = R / f"stage_e_2deg_{tag}_2026-09-26.md"
        L += [f"**{tag}**", ""] + ([l for l in p.read_text().splitlines() if l.startswith("|")] if p.exists() else ["(not yet available)"]) + [""]

    # 4. N/S spectra, mismatch null
    for name in ("ns_halo_spectra", "ns_mismatch_null"):
        L += [f"## {name}", ""]
        for tag in ("psf", "psfuc"):
            p = R / f"{name}_{tag}_2026-09-27.md"
            if not p.exists():
                L += [f"**{tag}**: (not yet available)", ""]; continue
            lines = p.read_text().splitlines()
            keep = [l for l in lines if l.startswith("Median over rows") or l.startswith("| data →") or l.startswith("| halo_flat") or l.startswith("| null_flat")]
            L += [f"**{tag}**", ""] + keep + [""]

    # 5. flexible lobes
    a = load("erosita_flex_ns_psf_2026-09-27.json"); b = load("erosita_flex_ns_psfuc_2026-09-27.json")
    L += ["## eROSITA lobe models (data medians: halo, S/N, West, East)", ""]
    if a and b:
        L += ["| lobes | old halo | new halo | old S/N W E | new S/N W E |", "|---|---|---|---|---|"]
        for m in ("none", "unif", "bright", "unif_flex", "bright_flex", "unif_mirrorW"):
            def s(x):
                D = [r for r in x["results"] if r["kind"] == "data"]
                h = np.median([sum(r["models"][m]["halo"].values()) for r in D]); rt = np.nanmedian([ratios(r["models"][m]["halo"]) for r in D], 0)
                return f"{h:,.0f}", f"{rt[0]:.2f} {rt[1]:.2f} {rt[2]:.2f}"
            (ha, ra), (hb, rb) = s(a), s(b)
            L.append(f"| {m} | {ha} | {hb} | {ra} | {rb} |")
    else:
        L.append("(not yet available)")
    L.append("")
    for name, date in (("erosita_flex_spectra", "2026-09-27"), ("erosita_flex_truth", "2026-09-27"), ("erosita_flex_score", "2026-09-27"),
                       ("hestia_halo_fits", "2026-09-27"), ("quadrant_calibration_cells2", "2026-09-29"), ("quadrant_calibration_cells1", "2026-09-29"),
                       ("inout_bubbles_cells2", "2026-09-29"), ("inout_bubbles_cells1", "2026-09-29"), ("sims_calibration_cells2", "2026-09-29")):
        L += [f"## {name}", ""]
        for tag in ("psf", "psfuc"):
            p = R / f"{name}_{tag}_{date}.md"
            L += [f"- {tag}: `{p.name}`" + ("" if p.exists() else " (not yet available)")]
        L.append("")
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
