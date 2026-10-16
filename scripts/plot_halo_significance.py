#!/usr/bin/env python
"""Significance of the Totani halo vs goodness of fit of the diffuse model,
for three bubble treatments (Kev, 2026-09-19).

    python scripts/plot_halo_significance.py            # BIC-penalized margin
    python scripts/plot_halo_significance.py --no-bic   # raw likelihood margin

Three panels, one per bubble template, everything else identical (the two
Loop I shells + the solar IC halo + isotropic + point sources + the diffuse
model):
  * Totani's 4.3 GeV residual construction (positive + negative templates,
    built per model)  -- from `rank_stage_c_<date>.json` (the consensus set);
  * the Ackermann+ 2014 flat mask        -- `rank_stage_b_<date>_flat_shells_sun.json`;
  * Oscar's 2017 structured map          -- `rank_stage_b_<date>_structured_shells_sun.json`.

x: the model's own margin over the nine bins, relative to the best PHYSICAL
   model in that panel (0 at the right edge, negative to the left) -- the
   "how much worse the interstellar model fits" axis of gce-fisher's
   `fig_evidence_vs_margin`. 🚨 BY DEFAULT IT IS BIC-PENALIZED:

       margin = -(1/2) [ (-2lnL + p ln n) - min over models of the same ]

   with p = (templates per bin) x 9 bins and n = 120 cells x 9 bins = 1080.
   Without the penalty the 13-template rows are compared with the 7-template
   ones on a likelihood that has 54 more free parameters -- the comparison
   this axis exists to avoid. `--no-bic` gives the raw likelihood margin.
y: significance of the Totani halo at 20.8 GeV, sign x sqrt(TS) of the
   NFW-rho^2 amplitude in his headline bin.
Marker: variant. Colour: halo counts at 20.8 GeV. Hollow: the quadrant test
fails at 20.8 GeV (one negative quadrant), i.e. the amplitude is not
morphologically a halo.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-09-24"
N_DATA = 120 * 9

#: The models Totani actually used. His baseline is SLZ6R30T150C2 in his own
#: two-template configuration (sec 2.3); his sec 3.4.4 systematics are the
#: five one-parameter neighbours below, plus the five-ring/three-ICS split of
#: the baseline. Everything else in the figure is grid he never looked at.
TOTANI_BASELINE = ("SL_Z6_R30_T150_C2", "gas1_ics1")
TOTANI_VARIANTS = {
    ("SS_Z6_R30_T150_C2", "gas1_ics1"),        # SNR cosmic-ray source distribution
    ("SL_Z10_R30_T150_C2", "gas1_ics1"),       # z_h = 10 kpc
    ("SL_Z6_R20_T150_C2", "gas1_ics1"),        # R_h = 20 kpc
    ("SL_Z6_R30_T100000_C2", "gas1_ics1"),     # T_S = 1e5 K
    ("SL_Z6_R30_T150_C5", "gas1_ics1"),        # E(B-V) cut 5 mag
    ("SL_Z6_R30_T150_C2", "gas5_ics3"),        # his "gas rings & 3 ICS" split
}
#: Corrected-fitter inputs (2026-09-23/24); DATE names the output figure.
PANELS = [
    # The 09-23 combined consensus took shells+NPS and no eROSITA, so it no longer
    # matches the other panels; the quadrant family scored the residual
    # construction under exactly this panel set (two shells + eROSITA unif N/S).
    (r"residual (4.3\,GeV) bubbles",
     ROOT / "results" / "fits" / "quadrant_family" / "2026-09-23__residual__shells__unif__nosun.json", "rows"),
    (r"flat 2014 mask", ROOT / "results" / "fits" / "rank_stage_b_2026-09-24_flat_shells_sun.json", "rows"),
    (r"structured 2017 map", ROOT / "results" / "fits" / "rank_stage_b_2026-09-24_structured_shells_sun.json", "rows"),
]


def load(path, key, bic: bool):
    d = json.loads(Path(path).read_text())
    out = []
    for r in d[key]:
        ts = r["ts"][5]; amp = r["halo_counts"][5]
        score = r["neg2lnl_total"] + (9 * r["n_free_per_bin"] * np.log(N_DATA) if bic else 0.0)
        out.append({"model": r["model"], "variant": r["variant"], "score": float(score),
                    "neg2lnl": r["neg2lnl_total"], "n_free": r["n_free_per_bin"],
                    "sig": float(np.sign(amp) * np.sqrt(max(ts, 0.0))), "halo": amp,
                    "quad_fail": bool(r["quadrant_any_negative"][5]), "ts_sum": float(sum(r["ts"]))})
    return out, d.get("structures") or d.get("consensus_selected")


#: Galp21 rows (Kev, 2026-09-23: isotropic IC is the fit-quality reference), from
#: the first-look JSONs, which score every Galp21 row under the three panel sets.
GALP21_JSONS = [ROOT / "results" / "fits" / "galp21_first_look_2026-09-23_iso4.json",
                ROOT / "results" / "fits" / "galp21_first_look_last2_2026-09-23.json"]
GALP21_SET = {0: "fig1_residual", 1: "fig1_flat", 2: "fig1_structured"}


def load_galp21(panel: int, bic: bool):
    """Galp21 rows scored under the panel's structure set, as `load` returns them."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from totani.selection import CANDIDATES
    out = []
    for jp in GALP21_JSONS:
        if not jp.exists():
            continue
        for r in json.loads(jp.read_text())["rows"]:
            if not r["model"].endswith("+iso"):
                continue
            st = r["sets"][GALP21_SET[panel]]
            # templates per bin: diffuse + isotropic + point sources + the set's templates
            n_diff = len(r["bare"]["norms"])
            n_set = sum(len(CANDIDATES[c]) if c in CANDIDATES else 1 for c in st["spec"])
            n_free = n_diff + 2 + n_set
            ts, amp = st["ts_20"], st["halo_20"]
            score = st["neg2lnl_total"] + (9 * n_free * np.log(N_DATA) if bic else 0.0)
            out.append({"model": r["model"], "variant": r["variant"], "score": float(score), "sig": float(np.sign(amp) * np.sqrt(max(ts, 0.0))),
                        "halo": amp, "quad_fail": bool(any(v < 0 for v in st["quadrant_counts_20"].values()))})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-bic", action="store_true")
    ap.add_argument("--tag", default="")
    ap.add_argument("--no-title", action="store_true", help="omit the in-figure title (the manuscript's caption carries it)")
    args = ap.parse_args()
    bic = not args.no_bic

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Presentation style, as gce-fisher's figures (Kev, 2026-08-24): real TeX
    # fonts, serif throughout. usetex needs latex + dvipng, both present here.
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath}")      # \tfrac, \text
    # Drawn at its printed size (Kev, 2026-10-02): the manuscript includes it at
    # \textwidth = 510 pt = 7.06 in (revtex4-2 prd), whose captions are 9 pt, so
    # every label here is >= 10 pt and the PDF is saved WITHOUT a tight bbox,
    # which would rescale it.
    FS = 10
    plt.rc("font", family="serif", size=FS)
    plt.rc("legend", fontsize=FS, frameon=False)
    plt.rc("xtick", labelsize=FS); plt.rc("ytick", labelsize=FS)
    from matplotlib.ticker import MaxNLocator

    W = 510 / 72.27
    fig, axes = plt.subplots(1, 3, figsize=(W, 4.1), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.865, top=0.93, bottom=0.36, wspace=0.12)
    vmin, vmax = 0, 9000
    stats = []
    for ax, (title, path, key) in zip(axes, PANELS):
        if not path.exists():
            ax.set_title(f"{title} (missing)"); continue
        rows, structures = load(path, key, bic)
        phys = [r for r in rows if r["model"] != "gll_iem_v07"]
        best = min(r["score"] for r in phys)
        for r in rows:
            r["margin"] = -0.5 * (r["score"] - best)
        n_fav = sum(1 for r in phys if r["sig"] > 2.0 and not r["quad_fail"])
        n_any = sum(1 for r in phys if r["sig"] > 2.0)
        stats.append((title, n_fav, n_any, len(phys)))
        for variant, mk, lab in (("gas1_ics1", "o", r"GALPROP, 1 gas $+$ 1 ICS"),
                                 ("gas5_ics3", "^", r"GALPROP, 5 rings $+$ 3 ICS")):
            for filled in (True, False):
                sub = [r for r in rows if r["variant"] == variant and r["quad_fail"] != filled]
                if not sub:
                    continue
                xs = [r["margin"] for r in sub]; ys = [r["sig"] for r in sub]
                if filled:
                    ax.scatter(xs, ys, c=[r["halo"] for r in sub], cmap="viridis", vmin=vmin, vmax=vmax,
                               marker=mk, s=16, edgecolor="k", linewidth=0.4, label=lab)
                else:
                    ax.scatter(xs, ys, facecolors="none", marker=mk, s=16, edgecolor="k", linewidth=0.6,
                               label=r"hollow: quadrant test fails")
        # ring the models Totani himself used, so his six can be found in the grid
        for r in rows:
            mv = (r["model"], r["variant"])          # not `key`: that names the panel's JSON key
            if mv == TOTANI_BASELINE:
                ax.scatter([r["margin"]], [r["sig"]], marker="o", s=95, facecolors="none",
                           edgecolor="crimson", linewidth=1.5, zorder=5,
                           label=r"Totani's baseline model")
            elif mv in TOTANI_VARIANTS:
                ax.scatter([r["margin"]], [r["sig"]], marker="o", s=70, facecolors="none",
                           edgecolor="crimson", linewidth=0.8, zorder=5,
                           label=r"Totani's \S3.4.4 variants")
        for r in rows:
            if r["model"] == "a2020":
                # hollow when the quadrant test fails, exactly as for the GALPROP rows
                ax.scatter([r["margin"]], [r["sig"]], marker="*", s=80,
                           facecolors="none" if r["quad_fail"] else "k", edgecolor="k",
                           linewidth=0.7, zorder=6,
                           label=r"A2020 ring model" if not r["quad_fail"] else r"hollow: quadrant test fails")
        # Galp21 (isotropic IC): diamonds, on the same margin axis and colour scale
        g21 = load_galp21(list(PANELS).index((title, path, key)), bic)
        for r in g21:
            r["margin"] = -0.5 * (r["score"] - best)
        for filled in (True, False):
            sub = [r for r in g21 if r["quad_fail"] != filled]
            if not sub:
                continue
            xs = [r["margin"] for r in sub]; ys = [r["sig"] for r in sub]
            if filled:
                ax.scatter(xs, ys, c=[r["halo"] for r in sub], cmap="viridis", vmin=vmin, vmax=vmax,
                           marker="D", s=20, edgecolor="darkorange", linewidth=0.8, zorder=4,
                           label=r"Galp21 (v57, isotropic IC)")
            else:
                ax.scatter(xs, ys, facecolors="none", marker="D", s=20, edgecolor="darkorange", linewidth=0.8,
                           zorder=4, label=r"hollow: quadrant test fails")
        xmin = min([r["margin"] for r in phys] + [r["margin"] for r in g21])
        ax.set_xlim(xmin * 1.10, -xmin * 0.08)
        ax.axhline(0, color="0.4", lw=0.9)
        ax.axhline(2, color="0.6", lw=0.7, ls=":")
        ax.set_title(title, fontsize=FS + 0.5)
        # in an empty corner of each panel: upper right, except the flat-mask panel's lower left
        corner = (0.04, 0.10, "left", "bottom") if "flat" in title else (0.97, 0.96, "right", "top")
        ax.text(corner[0], corner[1], rf"{n_fav}/{len(phys)} favor a halo", transform=ax.transAxes, ha=corner[2], va=corner[3],
                fontsize=FS, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
        ax.xaxis.set_major_locator(MaxNLocator(3, prune="upper" if ax is not axes[-1] else None))
        # gll_iem_v07 is off scale on the margin axis; it is reported in the
        # caption and the Markdown, not drawn (Kev, 2026-09-22)
        ref = next((r for r in rows if r["model"] == "gll_iem_v07"), None)
        stats[-1] = stats[-1] + ((ref["sig"], ref["margin"]) if ref is not None else (None, None),)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"significance at 20.8\,GeV, sign$\,\times\sqrt{\mathrm{TS}}$", fontsize=FS)
    handles, labels = [], []
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in labels:
                handles.append(h); labels.append(l)
    # one hollow entry stands for every hollow marker; drawn as a hollow circle
    from matplotlib.lines import Line2D
    keep = [(h, l) for h, l in zip(handles, labels) if not l.startswith("hollow")]
    keep.append((Line2D([], [], marker="o", ls="none", mfc="none", mec="k", ms=4.5, mew=0.6), r"hollow: quadrant test fails"))
    order = [r"GALPROP, 1 gas $+$ 1 ICS", r"GALPROP, 5 rings $+$ 3 ICS", r"A2020 ring model", r"Galp21 (v57, isotropic IC)",
             r"Totani's baseline model", r"Totani's \S3.4.4 variants", r"hollow: quadrant test fails"]
    keep.sort(key=lambda x: order.index(x[1]) if x[1] in order else 99)
    fig.legend([h for h, _ in keep], [l for _, l in keep], loc="upper center", ncol=3, bbox_to_anchor=(0.47, 0.175),
               handletextpad=0.3, columnspacing=1.0, borderaxespad=0.2)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cax = fig.add_axes([0.878, 0.36, 0.015, 0.57])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(r"halo counts at 20.8\,GeV", fontsize=FS); cb.ax.tick_params(labelsize=FS)
    lead = (r"$\left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm model} - \left[\ln\mathcal{L} - \tfrac{1}{2}p\ln n\right]_{\rm best}$"
            if bic else r"$\ln\mathcal{L}_{\rm model} - \ln\mathcal{L}_{\rm best}$")
    fig.text(0.47, 0.25, lead + (r"\quad(BIC-penalized fit quality)" if bic else r"\quad(fit quality)"),
             ha="center", va="center", fontsize=FS)
    rows0, _ = load(PANELS[1][1], PANELS[1][2], bic)
    n_models = len({r["model"] for r in rows0 if r["variant"] != "-"})
    # 🚨 THE COMPANION SET IS READ FROM THE JSONs, NEVER HARD-CODED, AND THE
    # THREE PANELS MUST AGREE. The residual panel takes its structures from
    # the greedy consensus, whose Loop I pick flip-flops between "shells" and
    # "shells + NPS" (they differ by 7 units in 7,000); twice the caption has
    # claimed a companion set the panels did not share. Now it cannot.
    # Stage C records CANDIDATE names ("loopI_shells"), Stage B TEMPLATE names
    # ("loopI_shell1", "loopI_shell2"); expand the former so they compare.
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "src"))
    from totani.selection import CANDIDATES
    comps = []
    for _, path, key in PANELS:
        st = load(path, key, bic)[1] or []
        flat = [t for x in st for t in CANDIDATES.get(x, (x,))]
        comps.append(tuple(sorted(x for x in flat if not x.startswith("bubble"))))
    if len(set(comps)) != 1:
        raise SystemExit("panels do not share companions, so they are not like-for-like:\n"
                         + "\n".join(f"  {p[0]}: {c}" for p, c in zip(PANELS, comps))
                         + "\nRe-score the odd ones out with matching --structures.")
    pretty = {"loopI_shell1": None, "loopI_shell2": r"Loop\,I shells",
              "loopI_haslam_nps": r"Haslam NPS", "loopI_haslam_raw": r"Haslam 408\,MHz",
              "sun_moon_ic": r"solar IC", "erosita_outer_N": None,
              "erosita_outer_S": r"eROSITA outer (N/S)", "loopIV": r"Loop\,IV"}
    companions = " $+$ ".join([pretty.get(c, c) for c in comps[0] if pretty.get(c, c)]
                              + [r"isotropic", r"point sources"])
    if not args.no_title:       # the paper's caption carries this; --no-title for the manuscript
        fig.suptitle(r"Significance of the Totani halo at 20.8\,GeV against the fit quality of the diffuse model, "
                     r"for three \emph{Fermi}-bubble templates" "\n"
                     rf"{n_models} Ackermann et al.\ (2012) GALPROP models $\times$ 2 component splits, plus the A2020 ring model and the six Galp21 (v57) models; "
                     + companions + r" in every panel", fontsize=13, y=1.00)
    tag = args.tag or ("" if bic else "_nobic")
    if args.no_title:
        tag += "_notitle"
    out = ROOT / "results" / f"halo_significance_vs_fit_{DATE}{tag}.png"
    fig.savefig(out, dpi=200)
    fig.savefig(out.with_suffix(".pdf"))            # no tight bbox: the page is exactly \textwidth wide
    print(f"saved {out} (+ .pdf)")
    ref_lines = []
    for t, nf, na, n, (rsig, rmargin) in stats:
        print(f"  {t:42s} {nf} of {n} rows favor a halo ({na} >2 sigma regardless)"
              + (f"; gll_iem_v07 reference {rsig:.1f} sigma at margin {rmargin:,.0f}, off scale" if rsig is not None else ""))
        if rsig is not None:
            ref_lines.append(f"{t.replace(chr(92) + ',', ' ')}: {rsig:.1f}σ at margin {rmargin:,.0f}")

    lines = [f"# Halo significance vs fit quality — {DATE}", "",
             ("BIC-penalized margin = −½[(−2lnL + p ln n) − best], p = templates × 9 bins, n = 1080"
              if bic else "raw likelihood margin = −½[(−2lnL) − best]")
             + "; y = sign × √TS of the NFW-ρ² amplitude at 20.8 GeV; companions: " + companions.replace("\\_", "_")
             + " + isotropic + point sources.", "",
             "gll_iem_v07 (the LAT team's interpolated model; a reference, not a physical model) is off scale on the "
             "margin axis and is not drawn: " + "; ".join(ref_lines) + ".", ""]
    for title, path, key in PANELS:
        if not path.exists():
            continue
        rows, structures = load(path, key, bic)
        phys = [r for r in rows if r["model"] != "gll_iem_v07"]
        best = min(r["score"] for r in phys)
        lines += [f"## {title.replace(chr(92) + ',', ' ')}  (structures: {', '.join(structures)})", "",
                  "| rank | model | variant | n/bin | margin | halo @ 20.8 | √TS (signed) | ΣTS | quadrant |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for i, r in enumerate(sorted(rows, key=lambda r: r["score"]), 1):
            lines.append(f"| {i} | {r['model']} | {r['variant']} | {r['n_free']} | {-0.5*(r['score']-best):+,.0f} | "
                         f"{r['halo']:,.0f} | {r['sig']:+.1f} | {r['ts_sum']:.0f} | "
                         f"{'fail' if r['quad_fail'] else 'pass'} |")
        lines.append("")
    (ROOT / "results" / f"halo_significance_vs_fit_{DATE}{tag}.md").write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
