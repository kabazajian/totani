"""Stage C: greedy forward selection of large structures (plan §5, Kev's
decisions of 2026-09-19).

Every diffuse model starts from the MINIMAL set -- its own diffuse templates,
isotropic, point sources. Candidates are the structure library; each
candidate adds one or more templates with a free normalisation per bin. At
each step the Δ(−2lnL) (background-only, summed over the 9 bins) of every
remaining candidate is computed for every model; the candidate with the
largest MEDIAN Δ across the models is added (the consensus path), the other
members of its morphology GROUP are retired, and the procedure repeats until
the best median Δ falls below the threshold. Each model's OWN greedy path is
recorded too, so disagreement between models about structures is visible.

THE THRESHOLD is BIC-like: a candidate with p new parameters must buy
Δ(−2lnL) > p × ln(n_data), n_data = 120 cells × 9 bins. That is 63 for a
one-template candidate, 126 for a pair. Plan §5 says calibrate this on null
skies; until then it is a stated convention, recorded in the output.

THE TOTANI-RESIDUAL BUBBLE CANDIDATE (Kev: "a legitimate candidate -- it
treats things as he did") is built PER MODEL AND VARIANT from that model's own
4.3 GeV residual (`bubble_totani`), cached, and enters as a pair (positive,
negative) with the negative template's sign chosen per bin by the background
fit (his sec 3.2). It carries the injection recovery curve measured in
`bubble_totani_injection_2026-09-19` (0.41-0.45 at 4.3 GeV, 0.93 at 20.8 GeV,
negative below 3 GeV): a halo fitted beside it is reported with that flag.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import ranking
from .cells import halo_cell_mask
from .likelihood import fit_all_bins, restrict_stack

N_DATA = 120 * 9

#: candidate -> template names it adds. Names resolve through `resolve_template`.
CANDIDATES = {
    "bubble_flat": ("bubble_flat",),
    "bubble_structured": ("bubble_structured",),
    "bubble_totani": ("bubble_totani_pos", "bubble_totani_neg"),     # per model; neg sign free per bin
    "loopI_shells": ("loopI_shell1", "loopI_shell2"),
    "loopI_haslam_nps": ("loopI_haslam_nps",),
    "loopI_haslam_raw": ("loopI_haslam_raw",),
    "loopI_shells+nps": ("loopI_shell1", "loopI_shell2", "loopI_haslam_nps"),
    "erosita_outer": ("erosita_outer",),
    "erosita_outer_NS": ("erosita_outer_N", "erosita_outer_S"),
    # Oscar's brightness lobes and his matching geometry-only controls
    "erosita_bright_NS": ("erosita_bright_N", "erosita_bright_S"),
    "erosita_bright_outer_NS": ("erosita_bright_outer_N", "erosita_bright_outer_S"),
    "erosita_unif_NS": ("erosita_unif_N", "erosita_unif_S"),
    "erosita_unif_outer_NS": ("erosita_unif_outer_N", "erosita_unif_outer_S"),
    "sun_moon_ic": ("sun_moon_ic",),
    "loopII": ("loopII",),
    "loopIII": ("loopIII",),
    "loopIV": ("loopIV",),
}
#: mutually exclusive morphologies: one member per group may be selected
GROUPS = {
    "bubbles": ("bubble_flat", "bubble_structured", "bubble_totani"),
    "loopI": ("loopI_shells", "loopI_haslam_nps", "loopI_haslam_raw", "loopI_shells+nps"),
    "erosita": ("erosita_outer", "erosita_outer_NS", "erosita_bright_NS",
                "erosita_bright_outer_NS", "erosita_unif_NS", "erosita_unif_outer_NS"),
}
SIGN_FREE = {"bubble_totani_neg"}


def bic_threshold(n_params: int) -> float:
    return n_params * float(np.log(N_DATA))


# ---- per-model Totani-residual templates ---------------------------------------

def totani_residual_cells(code: str, variant: str, edges_gev, force: bool = False) -> dict:
    """{'bubble_totani_pos', 'bubble_totani_neg'} cell templates for one model,
    built as in his sec 3.1 with THAT model's diffuse maps, cached."""
    from .bubble_totani import bubble_image, fit_all_cells, pixel_model_for_bin
    from .stack import cell_template_from_shape

    p = ranking.cache_path("totani_residual", f"{code}__{variant}")
    if p.exists() and not force:
        z = np.load(p)
        return {"bubble_totani_pos": z["pos"], "bubble_totani_neg": z["neg"]}
    if code == "gll_iem_v07":
        raise KeyError("no residual construction for gll_iem_v07 (no pixel-level model of it)")
    if code == "a2020":
        # the A2020 ring model has its own pixel-level path in bubble_totani
        pm = pixel_model_for_bin(ranking.RUN, edges_gev, 2, diffuse="a2020", with_gce=True)
    else:
        target, L, B, _ = ranking._grid_and_basis(edges_gev)
        if code.startswith("G21_") or ranking.paperii_tar(code) is not None:
            cubes = ranking.galprop_cubes(code, variant, target, edges_gev)
        else:
            from .galprop import galprop_on_target
            split = variant == "gas5_ics3"
            cubes = galprop_on_target(target, edges_gev, model_set=ranking.code_to_tag(code),
                                      ics_split=split, gas_rings=split)
        pm = pixel_model_for_bin(ranking.RUN, edges_gev, 2, diffuse="galprop_totani", with_gce=True,
                                 galprop_cubes=cubes)
    fit = fit_all_cells(pm)
    img = bubble_image(pm, fit)
    pos = cell_template_from_shape(ranking.RUN, edges_gev, img["positive"])
    neg = cell_template_from_shape(ranking.RUN, edges_gev, img["negative"])
    np.savez_compressed(p, pos=pos, neg=neg, _meta=json.dumps(
        {"code": code, "variant": variant, "amplitudes": fit["amplitudes"],
         "built": time.strftime("%Y-%m-%d %H:%M")}))
    return {"bubble_totani_pos": pos, "bubble_totani_neg": neg}


# ---- evaluation ---------------------------------------------------------------------

def resolve_template(name: str, model: str, variant: str, edges_gev, sh: dict) -> np.ndarray:
    if name in sh:
        return sh[name]
    if name.startswith("bubble_totani_"):
        return totani_residual_cells(model, variant, edges_gev)[name]
    return ranking.structure_cells(name, edges_gev)


def base_stack(model: str, variant: str, edges_gev, sh: dict) -> dict:
    """Diffuse + isotropic + point sources, nothing else."""
    t = {"isotropic": sh["isotropic"], "point_sources": sh["point_sources"]}
    if model == "a2020":
        t.update({k[len("a2020_"):]: v for k, v in sh.items() if k.startswith("a2020_")})
    elif model == "gll_iem_v07":
        t["giem"] = ranking.structure_cells("gll_iem_v07", edges_gev)
    else:
        t.update(ranking.model_cells(model, variant, edges_gev))
    return {"counts": sh["counts"], "n_bins": sh["counts"].shape[0], "templates": t}


def bg_neg2lnl(stack) -> np.ndarray:
    keep = halo_cell_mask()
    sr = restrict_stack(stack, keep)
    return np.array([row["neg2lnl"] for row in fit_all_bins(sr, halo_mask=keep)["per_bin"]])


def with_templates(stack, new: dict) -> dict:
    t = dict(stack["templates"]); t.update(new)
    return {"counts": stack["counts"], "n_bins": stack["n_bins"], "templates": t}


def add_candidate(stack, cand: str, model, variant, edges_gev, sh) -> tuple:
    """The stack with the candidate added. A sign-free template is tried with
    both signs per bin and the better kept; returns (stack, neg2lnl per bin,
    signs)."""
    new = {n: resolve_template(n, model, variant, edges_gev, sh) for n in CANDIDATES[cand]}
    free = [n for n in new if n in SIGN_FREE]
    if not free:
        st = with_templates(stack, new)
        return st, bg_neg2lnl(st), None
    n_free = free[0]
    fixed = {n: v for n, v in new.items() if n != n_free}
    plus = with_templates(stack, {**fixed, n_free: new[n_free]})
    minus = with_templates(stack, {**fixed, n_free: -new[n_free]})
    vp, vm = bg_neg2lnl(plus), bg_neg2lnl(minus)
    signs = np.where(vp <= vm, 1, -1)
    best = np.minimum(vp, vm)
    arr = new[n_free] * signs[:, None, None]
    return with_templates(stack, {**fixed, n_free: arr}), best, signs.tolist()


def greedy_for_model(model, variant, edges_gev, sh, candidates=None, threshold=bic_threshold):
    """One model's own greedy path. Returns the list of steps."""
    cands = list(candidates or CANDIDATES)
    stack = base_stack(model, variant, edges_gev, sh)
    cur = bg_neg2lnl(stack)
    steps, chosen = [], []
    while cands:
        deltas = {}
        for c in cands:
            _, v, _ = add_candidate(stack, c, model, variant, edges_gev, sh)
            deltas[c] = float((cur - v).sum())
        best = max(deltas, key=deltas.get)
        thr = threshold(9 * len(CANDIDATES[best]))
        steps.append({"candidates": deltas, "best": best, "delta": deltas[best], "threshold": thr,
                      "accepted": deltas[best] > thr})
        if deltas[best] <= thr:
            break
        stack, cur, _ = add_candidate(stack, best, model, variant, edges_gev, sh)
        chosen.append(best)
        cands = [c for c in cands if c != best and not any(
            best in g and c in g for g in GROUPS.values())]
    return {"model": model, "variant": variant, "chosen": chosen, "steps": steps,
            "final_neg2lnl": float(cur.sum())}


def fraction_rule_step(phys: list, usable: list, cutoff: float) -> dict:
    """One consensus step by the FRACTION-OF-ROWS rule (Kev, 2026-09-21).

    `phys` is one {candidate: Δ(−2lnL)} dict per physical row. A candidate is
    ELIGIBLE if it clears its own BIC threshold in >= `cutoff` of the rows;
    the one selected is the eligible candidate that is the largest improvement
    in the most rows. With nothing eligible, `best` is the nearest miss and
    `accepted` is False. Shared by `scripts/rank_stage_c.py` (the data) and
    `scripts/stage_e.py` (simulated skies) so both apply the same rule.
    """
    thr = {c: bic_threshold(9 * len(CANDIDATES[c])) for c in usable}
    clears = {c: float(np.mean([r[c] > thr[c] for r in phys])) for c in usable}
    eligible = [c for c in usable if clears[c] >= cutoff]
    if eligible:
        tops = [max(eligible, key=lambda c: r[c]) for r in phys]
        best_frac = {c: float(np.mean([t == c for t in tops])) for c in eligible}
        best, accepted = max(best_frac, key=best_frac.get), True
    else:
        best_frac = {}
        best, accepted = max(clears, key=clears.get), False     # the nearest miss, for the record
    return {"threshold_per_candidate": thr, "clears_fraction": clears, "eligible": eligible,
            "best_fraction_among_eligible": best_frac, "best": best, "accepted": accepted}


def consensus_step(rows, cands):
    """Given per-model candidate Δ dicts for one step, pick by median Δ."""
    med = {c: float(np.median([r[c] for r in rows])) for c in cands}
    lo = {c: float(np.min([r[c] for r in rows])) for c in cands}
    hi = {c: float(np.max([r[c] for r in rows])) for c in cands}
    best = max(med, key=med.get)
    return best, med, lo, hi
