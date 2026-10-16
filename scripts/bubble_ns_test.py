#!/usr/bin/env python
"""Oscar's item 4: give the Fermi bubbles themselves independent north/south
normalizations, for the flat and the structured morphology.

    python scripts/bubble_ns_test.py

For representative rows: companions = the split's no-residual selection minus
its bubble and eROSITA members (Haslam raw + solar IC), then
  bubble in {flat, structured} x norms in {1, N/S} x eROSITA in {none, uniform N/S, brightness N/S}.
Outputs per configuration: -2lnL(H0) and the gain of the second bubble
normalization, the halo at 20.8 GeV with interval and TS, the full halo
spectrum, and the quadrant gate. Cell space, 10 deg (cells never straddle b=0).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking, selection  # noqa: E402
from totani.cells import cell_centers, halo_cell_mask  # noqa: E402
from totani.likelihood import fit_all_bins, fit_all_bins_profile, quadrant_masks, restrict_stack  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
DATE = "2026-09-24"
K20 = 5
ROWS = [("SL_Z6_R30_T150_C2", "gas1_ics1"), ("SO_Z6_R30_T150_C5", "gas1_ics1"), ("SL_Z10_R30_T150_C5", "gas1_ics1"),
        ("SL_Z6_R30_T150_C2", "gas5_ics3"), ("SS_Z10_R30_T150_C5", "gas5_ics3"), ("a2020", "-"),
        ("G21_SA0_R12+iso", "gas1_ics1"), ("G21_SA100_F98+iso", "gas5_ics3")]
ROWS = [r for r in ROWS if ranking.available(*r)]          # Galp21 rows only when their cubes are ingested
COMPANIONS = ("loopI_haslam_raw", "sun_moon_ic")
EROSITA = {"none": (), "unif": ("erosita_unif_NS",), "bright": ("erosita_bright_NS",)}


def main() -> int:
    sh = ranking.shared_cells(EDGES); J = sh["J"]; keep = halo_cell_mask()
    _, Bc = cell_centers(); north = (Bc > 0)[None]
    for b in ("bubble_flat", "bubble_structured"):
        sh[f"{b}_N"] = sh[b] * north; sh[f"{b}_S"] = sh[b] * ~north
        selection.CANDIDATES[f"{b}_NS"] = (f"{b}_N", f"{b}_S")
    qm = quadrant_masks()
    out = {"date": DATE, "rows": []}
    lines = [f"# Fermi bubbles with independent north/south normalizations ({DATE})", "",
             "Companions: Haslam raw + solar IC (the no-residual selections' non-bubble, non-eROSITA members). "
             "Δ(−2lnL) of the second bubble normalization is background-only (H0), summed over 9 bins; the halo columns are H1.", "",
             "| row | bubble | eROSITA | bubble norms | Δ(−2lnL) 2 vs 1 norms | halo@20.8 [68%] | TS@20.8 | gate | halo spectrum (counts, 9 bins) |",
             "|---|---|---|---|---|---|---|---|---|"]
    for row in ROWS:
        for b in ("bubble_flat", "bubble_structured"):
            for ename, ec in EROSITA.items():
                res = {}
                for norms in ("1", "NS"):
                    cands = ((b,) if norms == "1" else (f"{b}_NS",)) + COMPANIONS + ec
                    st = selection.base_stack(row[0], row[1], EDGES, sh)
                    for c in cands:
                        st, _, _ = selection.add_candidate(st, c, row[0], row[1], EDGES, sh)
                    sr, Jr = restrict_stack(st, keep, J)
                    h0 = float(sum(x["neg2lnl"] for x in fit_all_bins(sr, halo_mask=keep)["per_bin"]))
                    pf = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=True)["per_bin"]
                    roi = [float(Jr[k][keep].sum()) for k in range(9)]
                    halo = [pf[k]["signal_amplitude"] * roi[k] for k in range(9)]
                    ts20 = fit_all_bins(sr, halo_mask=keep)["per_bin"][K20]["neg2lnl"] - pf[K20]["neg2lnl"]
                    quads = {}
                    for qn, m in qm.items():
                        stq, Jq = restrict_stack(st, m, J)
                        quads[qn] = fit_all_bins_profile({"counts": stq["counts"][K20:K20+1], "n_bins": 1, "templates": {n: t[K20:K20+1] for n, t in stq["templates"].items()}},
                                                         Jq[K20:K20+1], halo_mask=m, with_uncertainty=False)["per_bin"][0]["signal_amplitude"] * float(Jq[K20][m].sum())
                    gate = ts20 > 4 and halo[K20] > 0 and all(q > 0 for q in quads.values())
                    res[norms] = {"neg2lnl_H0": h0, "halo": halo, "halo20_lo": pf[K20]["signal_amplitude_lo"] * roi[K20],
                                  "halo20_hi": pf[K20]["signal_amplitude_hi"] * roi[K20], "ts20": ts20, "quadrants": quads, "gate": bool(gate)}
                d = res["1"]["neg2lnl_H0"] - res["NS"]["neg2lnl_H0"]
                rec = {"model": row[0], "variant": row[1], "bubble": b, "erosita": ename, "gain_second_norm": d, **{f"norms_{k}": v for k, v in res.items()}}
                out["rows"].append(rec)
                for norms in ("1", "NS"):
                    r = res[norms]
                    lines.append(f"| {row[0]} {row[1]} | {b[7:]} | {ename} | {norms} | {d:+.0f} | {r['halo'][K20]:,.0f} [{r['halo20_lo']:,.0f}, {r['halo20_hi']:,.0f}] | "
                                 f"{r['ts20']:.1f} | {'Y' if r['gate'] else 'n'} | " + " ".join(f"{x:,.0f}" for x in r["halo"]) + " |")
                print(lines[-2] + "\n" + lines[-1], flush=True)
    (ROOT / "results" / "fits" / f"bubble_ns_test_{DATE}.json").write_text(json.dumps(out, indent=1))
    (ROOT / "results" / f"bubble_ns_test_{DATE}.md").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
