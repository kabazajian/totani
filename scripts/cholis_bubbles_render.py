"""Render results/cholis_bubbles_2026-10-11.md from the data and simulation JSON of
scripts/cholis_bubbles.py (Cholis, Hu & Zhong bubble model inside our likelihood)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATE = "2026-10-11"
RESP = [3, 4, 5, 6, 7, 8]
E = np.array([1.5, 2.6, 4.3, 7.3, 12.3, 20.8, 35.1, 59.2, 100.0])
LABEL = {"flat": "2014 flat", "structured": "structured", "chz_flatE": "CHZ flat (edges E)",
         "chz_flatE_cocoons": "CHZ flat + cocoons", "chz_full": "CHZ flat + inner + cocoons"}


def _rows(res):
    return list(dict.fromkeys(tuple(r["row"]) for r in res))


def render():
    from halo_spectra_compare import conversion
    conv = conversion()
    e2 = lambda c, k: c * conv[k] * (E[k] * 1e3) ** 2 / 1e3 * 1e7      # counts -> E^2 dN/dE, 1e-7 GeV cm^-2 s^-1
    fits = ROOT / "results" / "fits"
    D = json.loads((fits / f"cholis_bubbles_data_{DATE}.json").read_text())["results"]
    sims_p = fits / f"cholis_bubbles_sims_{DATE}.json"
    S = json.loads(sims_p.read_text()) if sims_p.exists() else None
    L = [f"# The Cholis, Hu & Zhong bubble model in our likelihood — 2° cells, ULTRACLEAN response ({DATE})", "",
         "Their three components (Zenodo 23199004: flat bubbles with edges 'E', Inner Bubbles, Cocoons) each with a free "
         "normalization per bin, replacing our bubble template; companions as in Table I (eROSITA uniform N/S, Haslam Loop I, "
         "solar IC). A(20.8) = NFW counts at 20.8 GeV with the 68% profile interval; Z₂₀₁₈ over 5.6–130 GeV (K = 6). "
         "Δ(−2lnL) over 5.6–130 GeV relative to the 2014 flat template (negative = better), background only (H0) and with NFW (H1). "
         "TS(NFW) = Σ over 5.6–130 GeV of −2lnL(H0) − (−2lnL(H1)).", ""]
    full = [r for r in D if r["kind"] == "data" and r["region"] == "all"]
    L += ["## 1. Full fitted region", "", "| row | bubbles | A(20.8) [68%] | Z₂₀₁₈ | TS(NFW) | Δ H0 vs flat | Δ H1 vs flat |", "|---|---|---|---|---|---|---|"]
    for row in _rows(full):
        g = {r["config"]: r for r in full if tuple(r["row"]) == row}
        s = lambda r, key: sum(r[key][k] for k in RESP)
        for c in LABEL:
            if c not in g:
                continue
            r = g[c]
            L.append(f"| {row[0]} {row[1]} | {LABEL[c]} | {r['counts'][5]:,.0f} [{r['lo'][5]:,.0f}, {r['hi'][5]:,.0f}] | {r['resp']['Z_2018']:.1f} | "
                     f"{s(r, 'ts'):.0f} | {s(r, 'v0') - s(g['flat'], 'v0'):+.0f} | {s(r, 'v1') - s(g['flat'], 'v1'):+.0f} |")
    L += ["", "NFW spectrum with their full model, E²dN/dE (10⁻⁷ GeV cm⁻² s⁻¹) by bin (1.5 → 100 GeV):", "",
          "| row | " + " | ".join(f"{e:g}" for e in E) + " |", "|---|" + "---|" * 9]
    for r in [x for x in full if x["config"] == "chz_full"]:
        L.append(f"| {r['row'][0]} {r['row'][1]} | " + " | ".join(f"{e2(r['counts'][k], k):+.2f}" for k in range(9)) + " |")
    reg = [r for r in D if r["kind"] == "data" and r["region"] != "all"]
    if reg:
        jf = {r["region"]: r["J_frac"] for r in reg}
        L += ["", "## 2. Outside the bubbles", "",
              f"'outside' = cells with no 2014 flat-mask coverage (the region of Table V; {jf.get('outside', float('nan')):.1%} of the folded NFW counts). "
              f"'outside ∪ CHZ' additionally removes every cell their edges-E footprint touches ({jf.get('outside_union', float('nan')):.1%}).", "",
              "| row | bubbles | region | A(20.8) [68%] | Z₂₀₁₈ |", "|---|---|---|---|---|"]
        for row in _rows(reg):
            for c in ("flat", "structured", "chz_full"):
                for g in ("outside", "outside_union"):
                    x = [r for r in reg if tuple(r["row"]) == row and r["config"] == c and r["region"] == g]
                    if x:
                        r = x[0]
                        L.append(f"| {row[0]} {row[1]} | {LABEL[c]} | {'outside ∪ CHZ' if g == 'outside_union' else 'outside'} | "
                                 f"{r['counts'][5]:,.0f} [{r['lo'][5]:,.0f}, {r['hi'][5]:,.0f}] | {r['resp']['Z_2018']:.1f} |")
    Q = [r for r in D if r["kind"] == "dataq"]
    if Q:
        qn = [r["chz_full"] for r in S["results"] if r["kind"] == "qnull"] if S else []
        L += ["", "## 3. Quadrant tests (Q₆ = Σ over 5.6–130 GeV of Δ(−2lnL), common vs four regional amplitudes)", "",
              "| row | bubbles | Q₆ standard | Q₆ rotated | exceedances, CHZ null (std / rot) |", "|---|---|---|---|---|"]
        for r in Q:
            qs, qr = sum(r["std"][k] for k in RESP), sum(r["rot"][k] for k in RESP)
            ex = "—"
            if qn and r["config"] == "chz_full":
                ns = np.array([sum(x["std"][k] for k in RESP) for x in qn]); nr = np.array([sum(x["rot"][k] for k in RESP) for x in qn])
                ex = f"{int((ns >= qs).sum())}/{len(ns)} / {int((nr >= qr).sum())}/{len(nr)}"
            L.append(f"| {r['row'][0]} {r['row'][1]} | {LABEL[r['config']]} | {qs:.1f} | {qr:.1f} | {ex} |")
        if qn:
            ns = np.array([sum(x["std"][k] for k in RESP) for x in qn]); nr = np.array([sum(x["rot"][k] for k in RESP) for x in qn])
            L.append(f"\nNull (CHZ full model, true row + fitted common halo, {len(qn)} skies): standard median {np.median(ns):.1f}, 95% {np.percentile(ns, 95):.1f}; "
                     f"rotated median {np.median(nr):.1f}, 95% {np.percentile(nr, 95):.1f}.")
    if S:
        R = S["results"]; inj = np.array(S["injected_counts"])
        nul = [r for r in R if r["kind"] == "null"]; ij = [r for r in R if r["kind"] == "inj"]
        pct = lambda v: f"{np.median(v):,.0f} [{np.percentile(v, 16):,.0f}, {np.percentile(v, 84):,.0f}]"
        L += ["", f"## 4. Simulations from their model (true row {S['true_row'][0]} {S['true_row'][1]}, {S['skies']} skies each)", "",
              "| skies | fitted with | NFW counts at 20.8 GeV, median [16–84%] | Z₂₀₁₈ median [16–84%] |", "|---|---|---|---|"]
        for c in ("chz_full", "flat"):
            v = [r[c]["counts"][5] for r in nul]; z = [r[c]["Z_2018"] for r in nul]
            L.append(f"| no halo | {LABEL[c]} | {pct(v)} | {np.median(z):.1f} [{np.percentile(z,16):.1f}, {np.percentile(z,84):.1f}] |")
        v = [r["chz_full"]["counts"][5] for r in ij]; z = [r["chz_full"]["Z_2018"] for r in ij]
        L.append(f"| Totani halo injected ({inj[5]:,.0f} at 20.8 GeV) | {LABEL['chz_full']} | {pct(v)} | {np.median(z):.1f} [{np.percentile(z,16):.1f}, {np.percentile(z,84):.1f}] |")
        base = np.array([np.mean([r["chz_full"]["counts"][k] for r in nul]) for k in range(9)])
        rec = np.array([np.median([r["chz_full"]["counts"][k] for r in ij]) for k in range(9)])
        L += ["", "Recovery with their model: (median injected-sky NFW counts − mean no-halo NFW counts) / injected, by bin (1.5 → 100 GeV):", "",
              "| " + " | ".join(f"{e:g}" for e in E) + " |", "|" + "---|" * 9,
              "| " + " | ".join(f"{(rec[k]-base[k])/inj[k]:.2f}" if inj[k] > 0 else "—" for k in range(9)) + " |"]
        L += ["", "Spurious NFW spectrum, E²dN/dE (10⁻⁷ GeV cm⁻² s⁻¹), median over no-halo skies, by bin (1.5 → 100 GeV):", "",
              "| fitted with | " + " | ".join(f"{e:g}" for e in E) + " |", "|---|" + "---|" * 9]
        for c in ("chz_full", "flat"):
            L.append(f"| {LABEL[c]} | " + " | ".join(f"{e2(np.median([r[c]['counts'][k] for r in nul]), k):+.2f}" for k in range(9)) + " |")
        L.append(f"| injected Totani spectrum | " + " | ".join(f"{e2(inj[k], k):+.2f}" for k in range(9)) + " |")
    L += _followups(fits, e2)
    out = ROOT / "results" / f"cholis_bubbles_{DATE}.md"
    out.write_text("\n".join(L) + "\n"); print("\n".join(L))


def _load(fits, name):
    p = fits / f"cholis_bubbles_{name}_{DATE}.json"
    return json.loads(p.read_text())["results"] if p.exists() else None


def _followups(fits, e2):
    """Follow-ups 1-4 (Kev, 2026-10-11): their Loop I, spectral shapes, per-row quadrant nulls, 1 deg."""
    L = []
    D = json.loads((fits / f"cholis_bubbles_data_{DATE}.json").read_text())["results"]
    lp = _load(fits, "loop")
    if lp:
        L += ["", "## 5. Their Loop I shells in place of our Haslam template (2°)", "",
              "| row | A(20.8) [68%], our Loop I | Z₂₀₁₈ | A(20.8) [68%], their Loop I | Z₂₀₁₈ | Z₂₀₁₈ outside ∪ CHZ, their Loop I |", "|---|---|---|---|---|---|"]
        for row in _rows([r for r in lp if r["kind"] == "data"]):
            ours = next(r for r in D if r["kind"] == "data" and r["region"] == "all" and tuple(r["row"]) == row and r["config"] == "chz_full")
            th = next(r for r in lp if r["kind"] == "data" and r["region"] == "all" and tuple(r["row"]) == row)
            tho = next(r for r in lp if r["kind"] == "data" and r["region"] == "outside_union" and tuple(r["row"]) == row)
            L.append(f"| {row[0]} {row[1]} | {ours['counts'][5]:,.0f} [{ours['lo'][5]:,.0f}, {ours['hi'][5]:,.0f}] | {ours['resp']['Z_2018']:.1f} | "
                     f"{th['counts'][5]:,.0f} [{th['lo'][5]:,.0f}, {th['hi'][5]:,.0f}] | {th['resp']['Z_2018']:.1f} | {tho['resp']['Z_2018']:.1f} |")
        L += ["", "| row | Q₆ standard, their Loop I | Q₆ rotated, their Loop I |", "|---|---|---|"]
        for r in [x for x in lp if x["kind"] == "dataq"]:
            L.append(f"| {r['row'][0]} {r['row'][1]} | {sum(r['std'][k] for k in RESP):.1f} | {sum(r['rot'][k] for k in RESP):.1f} |")
    sh = _load(fits, "shape")
    if sh:
        L += ["", "## 6. Annihilation spectra as shapes with their bubbles (Sec. IV B test, 5.6–130 GeV, K = 6, 5 dof)", "",
              "⟨σv⟩ in 10⁻²⁵ cm³ s⁻¹ (our halo normalization); TS₀ = improvement over ⟨σv⟩ = 0; Δ (p) = shape statistic and its χ²₅ p-value.", "",
              "| row | Loop I | b 500 GeV: ⟨σv⟩, TS₀, Δ (p) | W 420 GeV: ⟨σv⟩, TS₀, Δ (p) |", "|---|---|---|---|"]
        for row in _rows(sh):
            for c, lab in (("chz_full", "ours"), ("chz_full_chzloop", "theirs")):
                x = next(y for y in sh if tuple(y["row"]) == row and y["config"] == c)["resp"]["spectra"]
                f = lambda k: f"{x[k]['sigmav_best']*1e25:+.1f}, {x[k]['ts_vs_zero']:.0f}, {x[k]['delta']:.1f} ({x[k]['p_shape']:.2g})"
                L.append(f"| {row[0]} {row[1]} | {lab} | {f('b 500')} | {f('W 420')} |")
    qr = _load(fits, "qnullrows")
    if qr:
        L += ["", "## 7. Quadrant nulls drawn from each row's own fit with their model (200 skies per row)", "",
              "| row | Q₆ standard (data) | exceedances | null median / 95% | Q₆ rotated (data) | exceedances | null median / 95% |", "|---|---|---|---|---|---|---|"]
        for row in _rows(qr):
            N = [r["chz_full"] for r in qr if tuple(r["row"]) == row]
            d = next(r for r in D if r["kind"] == "dataq" and tuple(r["row"]) == row and r["config"] == "chz_full")
            cells = []
            for k in ("std", "rot"):
                obs = sum(d[k][b] for b in RESP); nv = np.array([sum(n[k][b] for b in RESP) for n in N])
                cells += [f"{obs:.1f}", f"{int((nv >= obs).sum())}/{len(nv)}", f"{np.median(nv):.1f} / {np.percentile(nv, 95):.1f}"]
            L.append(f"| {row[0]} {row[1]} | " + " | ".join(cells) + " |")
    p1 = fits / f"cholis_bubbles_data_cells1_{DATE}.json"
    if p1.exists():
        D1 = json.loads(p1.read_text())["results"]
        L += ["", "## 8. 1° cells (outside regions carried from the 2° cells)", "",
              "| row | bubbles | A(20.8) [68%] | Z₂₀₁₈ | Z₂₀₁₈ outside | Z₂₀₁₈ outside ∪ CHZ |", "|---|---|---|---|---|---|"]
        for row in _rows([r for r in D1 if r["kind"] == "data"]):
            for c in ("flat", "structured", "chz_full"):
                g = {r["region"]: r for r in D1 if r["kind"] == "data" and tuple(r["row"]) == row and r["config"] == c}
                if "all" not in g:
                    continue
                a = g["all"]; zo = g.get("outside"); zu = g.get("outside_union")
                L.append(f"| {row[0]} {row[1]} | {LABEL[c]} | {a['counts'][5]:,.0f} [{a['lo'][5]:,.0f}, {a['hi'][5]:,.0f}] | {a['resp']['Z_2018']:.1f} | "
                         f"{zo['resp']['Z_2018']:.1f} | {zu['resp']['Z_2018']:.1f} |" if zo else
                         f"| {row[0]} {row[1]} | {LABEL[c]} | {a['counts'][5]:,.0f} [{a['lo'][5]:,.0f}, {a['hi'][5]:,.0f}] | {a['resp']['Z_2018']:.1f} | — | — |")
        L += ["", "| row | bubbles | Q₆ standard (1°) | Q₆ rotated (1°) |", "|---|---|---|---|"]
        for r in [x for x in D1 if x["kind"] == "dataq"]:
            L.append(f"| {r['row'][0]} {r['row'][1]} | {LABEL[r['config']]} | {sum(r['std'][k] for k in RESP):.1f} | {sum(r['rot'][k] for k in RESP):.1f} |")
    return L
