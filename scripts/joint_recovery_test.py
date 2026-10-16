#!/usr/bin/env python
"""Oscar's item 8b: joint spectral recovery with an injection drawn from the
EXACT spectral family being fitted (b-bbar, 500 GeV, PPPC4DMID), on the
null structured-bubble skies; the per-bin contributions to the joint TS; and
the joint fit repeated with higher lower-energy thresholds (a sensitivity
check, not a choice of range). The same threshold scan on the DATA for the
representative rows.

    python scripts/joint_recovery_test.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

import stage_e  # noqa: E402
from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import dmspectra, ranking, selection, simsky  # noqa: E402

EDGES = np.asarray(TOTANI_EDGES_GEV[:10], float)
DATE = "2026-09-24"
SV_INJ = 5e-25
N_SKIES = 20
SET = ("bubble_structured", "erosita_unif_NS", "loopI_haslam_raw", "sun_moon_ic")
ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
DATA_ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-")]
THRESH = [0, 1, 2, 3]            # first bin used: 1.16, 1.9, 3.2, 5.4 GeV lower edges


def halo_counts_for_sv(sv, J, edges):
    N = dmspectra.bin_yields("b", 500.0, edges)
    coef = N / (8.0 * np.pi * 500.0 ** 2) / dmspectra.bin_widths_mev(edges)
    return (coef * sv)[:, None, None] * J


def sliced(st, J, k0):
    return ({"counts": st["counts"][k0:], "n_bins": 9 - k0, "templates": {n: t[k0:] for n, t in st["templates"].items()}}, J[k0:], EDGES[k0:])


def bin_contrib(st, J, sv_best):
    """Per-bin TS contribution of the joint best fit: -2lnL(bin | sv=0) - -2lnL(bin | sv_best), backgrounds refit."""
    from totani.cells import halo_cell_mask
    from totani.likelihood import _fit_bg_given_offset, _mu_bg, _poisson_val_grad_mu, restrict_stack
    keep = halo_cell_mask(); sr, Jr = restrict_stack(st, keep, J); names = sorted(sr["templates"])
    N = dmspectra.bin_yields("b", 500.0, EDGES); coef = N / (8.0 * np.pi * 500.0 ** 2) / dmspectra.bin_widths_mev(EDGES)
    out = []
    for k in range(9):
        T = np.stack([sr["templates"][n][k] for n in names]); n = sr["counts"][k]
        sc = np.maximum(np.array([np.abs(T[i]).sum() for i in range(len(names))]), 1e-300); Tb = T / sc[:, None, None]
        seed = np.full(len(names), np.log10(max(n.sum() / (len(names) + 1), 1e-300)))
        vals = []
        for sv in (0.0, sv_best):
            th, _ = _fit_bg_given_offset(Tb, n, coef[k] * sv * Jr[k], seed)
            vals.append(_poisson_val_grad_mu(_mu_bg(th, Tb) + coef[k] * sv * Jr[k], n)[0])
        out.append(float(vals[0] - vals[1]))
    return out


def main() -> int:
    sh = ranking.shared_cells(EDGES); J = sh["J"]; shared = simsky.load_shared()
    truth = dict(np.load(stage_e.truth_path("null_structured"))); rowpix = simsky.load_row(*ROW)
    inj = halo_counts_for_sv(SV_INJ, J, EDGES)
    out = {"date": DATE, "sv_injected": SV_INJ, "injected_counts_per_bin": inj.sum(axis=(1, 2)).tolist(), "skies": [], "data": []}
    lines = [f"# Joint spectral recovery with an exact b b̄ 500 GeV injection at ⟨σv⟩ = {SV_INJ:.0e} ({DATE})", "",
             "Null structured-bubble skies (true row SO_Z6_R30_T150_C5, own structure set), halo injected with the fitted spectral family. "
             "Recovered ⟨σv⟩ (10⁻²⁵) with 68% interval, by lower-energy threshold (first bin used); per-bin TS contributions of the full-range fit.", ""]
    rows = []
    for i in range(N_SKIES):
        rng = np.random.default_rng([stage_e.TRUTH_SEED["null_structured"], i, 20260922])
        cells, pix = simsky.draw(truth, sh["counts"], shared, rng)
        cells = cells + np.where(stage_e.ranking_keep()[None], rng.poisson(inj), 0)      # inject on the analysis cells
        shr = {**sh, "counts": cells, **simsky.residual_cells(pix, rowpix, shared, EDGES)}
        st = selection.base_stack(ROW[0], ROW[1], EDGES, shr)
        for c in SET:
            st, _, _ = selection.add_candidate(st, c, ROW[0], ROW[1], EDGES, shr)
        rec = {"sim": i, "by_threshold": {}}
        for k0 in THRESH:
            s2, J2, E2 = sliced(st, J, k0)
            f = dmspectra.joint_fit(s2, J2, "b", 500, E2, n_scan=31)
            rec["by_threshold"][k0] = {"sv": f["sigmav_best"], "lo": f["sigmav_lo68"], "hi": f["sigmav_hi68"], "ts": f["ts_vs_zero"]}
        rec["bin_ts"] = bin_contrib(st, J, rec["by_threshold"][0]["sv"])
        out["skies"].append(rec); rows.append(rec)
        print(f"[{i+1}/{N_SKIES}] sv by threshold: " + " ".join(f"{rec['by_threshold'][k]['sv']*1e25:+.2f}" for k in THRESH)
              + "  bin TS: " + " ".join(f"{x:5.1f}" for x in rec["bin_ts"]), flush=True)
    lines += ["| first bin (GeV) | recovered ⟨σv⟩ median [16–84% over skies] | median 68% half-width | median TS |", "|---|---|---|---|"]
    lower = [1.16, 1.9, 3.2, 5.4]
    for k0 in THRESH:
        sv = np.array([r["by_threshold"][k0]["sv"] for r in rows]); hw = np.array([(r["by_threshold"][k0]["hi"] - r["by_threshold"][k0]["lo"]) / 2 for r in rows])
        ts = np.array([r["by_threshold"][k0]["ts"] for r in rows])
        lines.append(f"| {lower[k0]} | {np.median(sv)*1e25:+.2f} [{np.percentile(sv,16)*1e25:+.2f}, {np.percentile(sv,84)*1e25:+.2f}] | {np.median(hw)*1e25:.2f} | {np.median(ts):.0f} |")
    bt = np.median([r["bin_ts"] for r in rows], axis=0)
    lines += ["", "Median per-bin TS contribution (1.5 … 100 GeV): " + " ".join(f"{x:.1f}" for x in bt), "",
              "## Data: the same threshold scan", "", "| row | selection | ⟨σv⟩ by first bin 1.16 / 1.9 / 3.2 / 5.4 GeV | per-bin TS of the full-range fit |", "|---|---|---|---|"]
    for row in DATA_ROWS:
        v = row[1] if row[1] != "-" else "gas1_ics1"
        for residual in (True, False):
            cands = tuple(json.loads((ROOT / "results" / "fits" / f"rank_stage_c_2026-09-23_{v}{'' if residual else '_no_residual'}.json").read_text())["consensus_selected"])
            st = selection.base_stack(row[0], row[1], EDGES, sh)
            try:
                for c in cands:
                    st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
            except (FileNotFoundError, KeyError):
                continue
            svs = []
            for k0 in THRESH:
                s2, J2, E2 = sliced(st, J, k0); f = dmspectra.joint_fit(s2, J2, "b", 500, E2, n_scan=31); svs.append(f["sigmav_best"])
            bts = bin_contrib(st, J, svs[0])
            out["data"].append({"model": row[0], "variant": row[1], "residual": residual, "sv_by_threshold": svs, "bin_ts": bts})
            lines.append(f"| {row[0]} {row[1]} | {'with' if residual else 'without'} residual | " + " / ".join(f"{s*1e25:+.1f}" for s in svs) + " | " + " ".join(f"{x:.0f}" for x in bts) + " |")
            print(lines[-1], flush=True)
    (ROOT / "results" / "fits" / f"joint_recovery_test_{DATE}.json").write_text(json.dumps(out, indent=1))
    (ROOT / "results" / f"joint_recovery_test_{DATE}.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
