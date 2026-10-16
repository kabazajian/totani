#!/usr/bin/env python
"""Oscar's note of 2026-09-29, items 1 and 2: the decisive configurations at
one binning and response treatment, with the halo spectrum and profile
intervals, the conditional likelihood comparisons, and the template-detection
statistic of Macias et al. (2018) in his convention.

    python scripts/oscar_conditional.py --cell 2                 # unconvolved
    FINE_TAG=psf python scripts/oscar_conditional.py --cell 2    # PSF-folded
    FINE_TAG=psf python scripts/oscar_conditional.py --cell 1

Fixed-template reproduction (no reselection): for each representative row,
three bubble descriptions -- structured map, flat mask, Totani's residual
construction (built on his 10-deg cells from the data at 4.3 GeV, projected;
NOT response-folded, since it is made from counts) -- each with the same
companions (eROSITA uniform lobes N/S, Haslam Loop I, solar IC). Per bin:
background-only fit (H0) and the halo profile fit (H1, signed, refined grid
n_refine = 6, 68% profile interval), TS_e = -2lnL(H0) - (-2lnL(H1)), and
q_e^+ (TS_e if the signed optimum is >= 0, else 0). Summaries over all nine
bins (K = 9) and over the response-valid bins fixed in advance (PSF R68 <=
1/4 of the finest cell used, 1 deg: bins 5.6-130 GeV, K = 6): Q^+, p_K from
the binomial chi-square mixture, Z_2018 and Z_one. Conditional comparisons
between bubble descriptions are differences of -2lnL at THIS binning only.
Also recorded: min expected count in the fitted cells (signed fits must keep
mu > 0) and any materially negative TS (a fit check, not truncated).
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

from halo_spectra_compare import conversion  # noqa: E402
from totani import fine, ranking, stats  # noqa: E402
from totani.likelihood import expected_counts, fit_all_bins, fit_all_bins_profile, restrict_stack  # noqa: E402

DATE = "2026-09-29"
COMP = ("erosita_unif_N", "erosita_unif_S", "loopI_haslam_raw", "sun_moon_ic")
CONFIGS = {"structured": "bubble_structured", "flat": "bubble_flat", "construction": "bubble_totani"}
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("SO_Z6_R30_T150_C5", "gas1_ics1"),
        ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-"), ("G21_SA0_R12", "gas1_ics1")]
ROWS = [r for r in ROWS if ranking.available(*r)]          # the Galp21 row only when its cubes are ingested
RESP_BINS = [3, 4, 5, 6, 7, 8]
E = np.array([1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0])
N_REFINE = 6
_W: dict = {}


def _init(cell):
    _W["cell"] = cell; _W["lib"] = fine.load_library(cell)


def stack_for(counts, row_key, config):
    lib = _W["lib"]; row = fine.load_row(_W["cell"], *row_key)
    st = fine.base_stack(lib, row, counts)
    for n in COMP:
        st["templates"][n] = lib["templates"][n]
    st, _ = fine.add_candidate(st, CONFIGS[config], lib, row, lib["keep"])
    return st


def evaluate(counts, row_key, config, with_intervals=True):
    lib = _W["lib"]; keep, J = lib["keep"], lib["J"]
    st = stack_for(counts, row_key, config)
    sr, Jr = restrict_stack(st, keep, J)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=with_intervals, n_refine=N_REFINE)
    names = pf["component_order"]; pf = pf["per_bin"]
    roi = [float(Jr[k][keep].sum()) for k in range(9)]
    T = np.stack([sr["templates"][n] for n in names])
    out = {"v0": [], "v1": [], "ts": [], "a": [], "lo": [], "hi": [], "counts": [], "min_mu": []}
    for k in range(9):
        a = pf[k]["signal_amplitude"]
        bga = np.array([pf[k]["bg_amplitudes"][n] for n in names])
        mu = expected_counts(np.log10(np.clip(bga, 1e-300, None)), T[:, k]) + a * Jr[k]
        out["v0"].append(bg[k]["neg2lnl"]); out["v1"].append(pf[k]["neg2lnl"]); out["ts"].append(bg[k]["neg2lnl"] - pf[k]["neg2lnl"])
        out["a"].append(a); out["counts"].append(a * roi[k]); out["min_mu"].append(float(mu[keep].min()))
        out["lo"].append(pf[k].get("signal_amplitude_lo", np.nan) * roi[k]); out["hi"].append(pf[k].get("signal_amplitude_hi", np.nan) * roi[k])
    out["all"] = stats.summarize(out["ts"], out["a"])
    out["resp"] = stats.summarize(out["ts"], out["a"], RESP_BINS)
    return out


def task(t):
    row_key, config = t
    return {"row": list(row_key), "config": config, **evaluate(_W["lib"]["counts"], row_key, config)}


def render(res, cell, conv):
    tag = fine.tag()
    e2 = lambda c, k: c * conv[k] * (E[k] * 1e3) ** 2 / 1e3 * 1e7
    lines = [f"# Conditional comparisons at {cell:g}° cells{' (PSF-folded)' if tag else ' (no PSF)'} — Oscar's items 1–2 ({DATE})", "",
             "Fixed-template reproduction; companions eROSITA uniform N/S, Haslam Loop I, solar IC in every configuration. "
             "q⁺ per bin = TS if the signed optimum ≥ 0, else 0. Q⁺ = Σ q⁺; p from the binomial χ² mixture (Macias et al. 2018, Suppl. Eq. 7); "
             "Z₂₀₁₈ = Φ⁻¹(1 − p/2) (the 2018 convention), Z₁ = Φ⁻¹(1 − p) (one-sided). K = 9: all bins; K = 6: response-valid bins "
             "(PSF R68 ≤ ¼ of a 1° cell: 5.6–130 GeV), fixed before any result. Δ columns compare bubble descriptions at this binning only "
             "(negative = the first fits better).", ""]
    rows = list(dict.fromkeys(tuple(r["row"]) for r in res))
    lines += ["## Detection statistic per configuration", "",
              "| row | bubbles | halo 20.8 GeV [68%] | Q⁺ (K=9) | log₁₀p | Z₂₀₁₈ / Z₁ | Q⁺ (K=6) | log₁₀p | Z₂₀₁₈ / Z₁ | bins with a<0 | min μ | min TS |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        for c in CONFIGS:
            x = next(y for y in res if tuple(y["row"]) == r and y["config"] == c)
            A, R = x["all"], x["resp"]
            lines.append(f"| {r[0]} {r[1]} | {c} | {x['counts'][5]:,.0f} [{x['lo'][5]:,.0f}, {x['hi'][5]:,.0f}] | {A['Q_plus']:.1f} | {A['log10_p']:.1f} | "
                         f"{A['Z_2018']:.1f} / {A['Z_one']:.1f} | {R['Q_plus']:.1f} | {R['log10_p']:.1f} | {R['Z_2018']:.1f} / {R['Z_one']:.1f} | "
                         f"{sum(a < 0 for a in x['a'])} | {min(x['min_mu']):.2g} | {min(x['ts']):.2f} |")
    lines += ["", "## Per-bin TS (signed) and q⁺", "", "| row | bubbles | TS by bin (1.5 → 100 GeV) | q⁺ by bin |", "|---|---|---|---|"]
    for r in rows:
        for c in CONFIGS:
            x = next(y for y in res if tuple(y["row"]) == r and y["config"] == c)
            lines.append(f"| {r[0]} {r[1]} | {c} | " + " ".join(f"{t:.1f}" for t in x["ts"]) + " | " + " ".join(f"{q:.1f}" for q in x["all"]["q_plus"]) + " |")
    lines += ["", "## Spectra: E²dN/dE (×10⁻⁷ GeV cm⁻² s⁻¹) with 68% profile intervals", "", "| row | bubbles | " + " | ".join(f"{e:g}" for e in E) + " |",
              "|---|---|" + "---|" * 9]
    for r in rows:
        for c in CONFIGS:
            x = next(y for y in res if tuple(y["row"]) == r and y["config"] == c)
            cells = [f"{e2(x['counts'][k], k):+.2f} [{e2(x['lo'][k], k):+.2f}, {e2(x['hi'][k], k):+.2f}]" for k in range(9)]
            lines.append(f"| {r[0]} {r[1]} | {c} | " + " | ".join(cells) + " |")
    lines += ["", "## Conditional comparisons at this binning (Δ −2lnL)", "",
              "| row | flat − structured, H0 (all / K=6 bins) | flat − structured, H1 | construction − flat, H0 | construction − flat, H1 |", "|---|---|---|---|---|"]
    for r in rows:
        g = {c: next(y for y in res if tuple(y["row"]) == r and y["config"] == c) for c in CONFIGS}
        s = lambda c, key, bins=None: sum(g[c][key][k] for k in (bins or range(9)))
        f = lambda a, b, key: f"{s(a, key) - s(b, key):+.0f} / {s(a, key, RESP_BINS) - s(b, key, RESP_BINS):+.0f}"
        lines.append(f"| {r[0]} {r[1]} | {f('flat', 'structured', 'v0')} | {f('flat', 'structured', 'v1')} | {f('construction', 'flat', 'v0')} | {f('construction', 'flat', 'v1')} |")
    out = ROOT / "results" / f"oscar_conditional_cells{cell:g}{tag}_{DATE}.md"
    out.write_text("\n".join(lines) + "\n"); print("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--cell", type=float, default=2.0); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--render", action="store_true"); a = ap.parse_args()
    js = ROOT / "results" / "fits" / f"oscar_conditional_cells{a.cell:g}{fine.tag()}_{DATE}.json"
    conv = conversion()
    if a.render:
        render(json.loads(js.read_text())["results"], a.cell, conv); return 0
    t0 = time.time()
    tasks = [(r, c) for r in ROWS for c in CONFIGS]
    with Pool(a.workers, initializer=_init, initargs=(a.cell,)) as pool:
        res = []
        for k, x in enumerate(pool.imap_unordered(task, tasks), 1):
            res.append(x); print(f"[{k}/{len(tasks)}] {x['row']} {x['config']} ({time.time()-t0:.0f}s)", flush=True)
    js.write_text(json.dumps({"date": DATE, "cell": a.cell, "psf": fine.tag(), "resp_bins": RESP_BINS, "n_refine": N_REFINE, "results": res}, indent=1, default=float))
    render(res, a.cell, conv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
