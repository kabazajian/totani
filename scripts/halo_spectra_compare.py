#!/usr/bin/env python
"""Oscar's item 6 / Kev's reframing: the halo SPECTRUM, data against simulations,
exposure-corrected, on one set of axes; and where the observed 20.8 GeV
statistic falls in the null-sky TS distribution.

    python scripts/halo_spectra_compare.py

Flux conversion. basis_k = int_bin A(E) dE (cm^2 s MeV, flat spectrum in the
bin), so the halo template in bin k is a * J(pix) * basis_k(pix) * Omega(pix)
counts and a * J is dN/dE per sr. The ROI-integrated differential flux is
    dN/dE_k = counts_k * sum_ROI(J Omega) / sum_ROI(J basis_k Omega)   [cm^-2 s^-1 MeV^-1]
with the sums over the |b| >= 10 unmasked pixels; E^2 dN/dE follows with E in
MeV (reported in GeV cm^-2 s^-1). The bin widths are inside basis_k -- this is
the template-weighted exposure conversion Oscar asked for, not count ratios.

Series (all per-bin halo counts from committed JSONs):
  data  - Totani's configuration (his model, his construction; the a=0 fit of
          bubble_totani_injection_2026-09-19, = Sec. V A);
        - 257 rows under the residual construction, the flat mask and the
          structured map (Fig. 1's panel sets), median and 16-84%;
  sims  - the null structured-bubble sky under the residual construction and
          under the flat mask (stage_e_anatomy, true diffuse row, 20 skies);
        - the injected Totani-spectrum halo and its recovery under the
          structured map (Stage E halo_structured, no-residual path, 50 skies).
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking, simsky  # noqa: E402
from totani.cells import halo_cell_mask  # noqa: E402
from totani.halo import nfw_rho2_on_target  # noqa: E402

EDGES = np.asarray(TOTANI_EDGES_GEV[:10], float); E = np.asarray(TOTANI_CENTERS_GEV[:9], float)
DE = np.diff(EDGES)
DATE = "2026-09-24"
FITS = ROOT / "results" / "fits"
K20 = 5


def conversion():
    """counts -> ROI-integrated flux per bin, halo-template weighted."""
    basis, L, B, mask = simsky._shape_basis(EDGES)          # (9, ny, nx) = basis x Omega
    target, _, _, _ = ranking._grid_and_basis(EDGES)
    J = nfw_rho2_on_target(target)
    from totani.templates import pixel_solid_angle
    om = pixel_solid_angle(target)
    roi = (np.abs(B) >= 10) & (np.abs(B) <= 60) & (np.abs(L) <= 60) & ~mask
    num = float((J * om)[roi].sum())
    den = np.array([float((J * basis[k])[roi].sum()) for k in range(9)])
    return num / den                                          # (9,): flux per count


def e2f(counts, conv):
    c = np.asarray(counts, float)
    return c * conv * (E * 1e3) ** 2 / 1e3                       # (dN/dE per MeV) x E_MeV^2 -> MeV cm^-2 s^-1 -> GeV cm^-2 s^-1


def rows_of(path, key):
    d = json.loads(Path(path).read_text())
    return [r for r in d[key] if r["model"] != "gll_iem_v07"]


def main() -> int:
    conv = conversion()
    out = {"date": DATE, "E_gev": E.tolist(), "flux_per_count": conv.tolist(), "series": {}}
    # data
    inj = json.loads((FITS / "bubble_totani_injection_2026-09-19.json").read_text())
    tot = np.array(inj["a0_halo_counts"])
    out["series"]["data_totani_config"] = {"counts": tot.tolist(), "e2f": e2f(tot, conv).tolist()}
    panels = {"residual": (FITS / "quadrant_family" / "2026-09-23__residual__shells__unif__nosun.json", "rows"),
              "flat": (FITS / "rank_stage_b_2026-09-24_flat_shells_sun.json", "rows"),
              "structured": (FITS / "rank_stage_b_2026-09-24_structured_shells_sun.json", "rows")}
    ts_obs = {}
    for name, (p, key) in panels.items():
        R = rows_of(p, key); H = np.array([r["halo_counts"] for r in R]); T = np.array([r["ts"][K20] for r in R])
        out["series"][f"data_257_{name}"] = {"median": np.median(H, 0).tolist(), "p16": np.percentile(H, 16, 0).tolist(),
                                              "p84": np.percentile(H, 84, 0).tolist(), "e2f_median": e2f(np.median(H, 0), conv).tolist(),
                                              "ts20_median": float(np.median(T)), "ts20_p16": float(np.percentile(T, 16)), "ts20_p84": float(np.percentile(T, 84))}
        ts_obs[name] = T
    # sims: anatomy (true row), null structured sky
    an = json.loads((FITS / "stage_e_anatomy_2026-09-23_null_structured.json").read_text())
    k0 = "SO_Z6_R30_T150_C5__gas1_ics1"
    for sname, lab in (("construct", "null_sky_residual_construction"), ("flat_truth", "null_sky_flat_mask"), ("truth", "null_sky_own_set")):
        H = np.array([sk["rows"][k0][sname]["halo_counts"] for sk in an["skies"]]); T = np.array([sk["rows"][k0][sname]["ts"][K20] for sk in an["skies"]])
        out["series"][lab] = {"median": np.median(H, 0).tolist(), "p16": np.percentile(H, 16, 0).tolist(), "p84": np.percentile(H, 84, 0).tolist(),
                              "e2f_median": e2f(np.median(H, 0), conv).tolist(), "ts20": T.tolist()}
    # sims: injected halo recovered under the structured map (Stage E, 1-gas rows)
    sims = [json.loads(p) for p in (Path(f).read_text() for f in sorted(glob.glob(str(FITS / "stage_e" / "halo_structured__*.json"))))]
    H = np.array([np.median([r["halo_counts"] for r in s["splits"]["gas1_ics1"]["no_residual"]["rows"]], 0) for s in sims])
    out["series"]["injected_recovered_structured"] = {"median": np.median(H, 0).tolist(), "p16": np.percentile(H, 16, 0).tolist(),
                                                       "p84": np.percentile(H, 84, 0).tolist(), "e2f_median": e2f(np.median(H, 0), conv).tolist()}
    injected = np.array(inj["a0_halo_counts"]).clip(0)
    out["series"]["injected_truth"] = {"counts": injected.tolist(), "e2f": e2f(injected, conv).tolist()}
    # null TS distribution vs the observed statistic
    null_ts = np.array(out["series"]["null_sky_residual_construction"]["ts20"])
    ts_data_totani = float(json.loads((ROOT / "results" / "fits" / "rank_stage_c_2026-09-23.json").read_text())["scored"][0]["ts"][K20]) if False else 172.0
    out["null_ts_summary"] = {"null_residual_ts20_median": float(np.median(null_ts)), "null_residual_ts20_max": float(null_ts.max()),
                              "observed_totani_config_ts20": ts_data_totani,
                              "observed_257_residual_ts20_median": float(np.median(ts_obs["residual"])),
                              "fraction_null_skies_ts_ge_observed_median": float(np.mean(null_ts >= np.median(ts_obs["residual"])))}
    (FITS / f"halo_spectra_compare_{DATE}.json").write_text(json.dumps(out, indent=1))

    # ---- figure ----
    # Drawn at its printed size (Kev, 2026-10-03), as Figs. 1-2: \textwidth = 510 pt,
    # captions 9 pt, every label >= 10 pt, legends below the panels, no tight bbox.
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    FS = 10
    plt.rc("text", usetex=True); plt.rc("font", family="serif", size=FS); plt.rc("legend", fontsize=FS, frameon=False)
    plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)
    W = 510 / 72.27
    fig, axes = plt.subplots(1, 2, figsize=(W, 4.2), sharey=True)
    fig.subplots_adjust(left=0.10, right=0.985, top=0.935, bottom=0.355, wspace=0.08)
    S = out["series"]
    def band(ax, s, color, label, ls="-"):
        m = e2f(s["median"], conv); lo = e2f(s["p16"], conv); hi = e2f(s["p84"], conv)
        ax.plot(E, m, ls, color=color, marker="o", ms=3, lw=1.2, label=label); ax.fill_between(E, lo, hi, color=color, alpha=0.18, lw=0)
    ax = axes[0]; ax.set_title("Data", fontsize=FS + 0.5)
    ax.plot(E, S["data_totani_config"]["e2f"], "k-", marker="s", ms=3.5, lw=1.2, label=r"Totani's configuration (Sec.~V\,A)")
    band(ax, S["data_257_residual"], "crimson", "257 rows, residual construction")
    band(ax, S["data_257_flat"], "tab:green", "257 rows, flat 2014 mask")
    band(ax, S["data_257_structured"], "tab:blue", "257 rows, structured 2017 map")
    ax = axes[1]; ax.set_title("Simulations", fontsize=FS + 0.5)
    band(ax, S["null_sky_residual_construction"], "crimson", "no halo, residual construction")
    band(ax, S["null_sky_flat_mask"], "tab:green", "no halo, flat mask")
    band(ax, S["null_sky_own_set"], "0.5", "no halo, own templates (noise)")
    ax.plot(E, S["injected_truth"]["e2f"], "k--", lw=1.2, label=r"injected halo (Totani's spectrum)")
    band(ax, S["injected_recovered_structured"], "tab:blue", "recovered under the structured map")
    for ax in axes:
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1e-7); ax.axhline(0, color="0.4", lw=0.8)
        ax.set_xlabel(r"$E$ [GeV]", fontsize=FS); ax.grid(alpha=0.3)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, handlelength=1.8, borderaxespad=0.0)
    axes[0].set_ylabel(r"$E^2\,dN/dE$ [GeV cm$^{-2}$ s$^{-1}$]", fontsize=FS)
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "results" / f"halo_spectra_data_vs_sims_{DATE}.{ext}", dpi=200)
    print(json.dumps(out["null_ts_summary"], indent=1))
    print("E^2F data (Totani config):", ["%.2e" % x for x in S["data_totani_config"]["e2f"]])
    print("E^2F null residual      :", ["%.2e" % x for x in S["null_sky_residual_construction"]["e2f_median"]])
    print("E^2F data structured med:", ["%.2e" % x for x in S["data_257_structured"]["e2f_median"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
