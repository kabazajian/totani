"""Stage D: annihilation spectra from PPPC4DMID and the joint DM fit.

Kev, 2026-09-19: "DM spectra for Stage D: we can get this from PPPC4DMID."
The table is gce-fisher's copy (`data/external/PPPC4DMID/AtProduction_gammas.dat`,
Cirelli et al. 2011, arXiv:1012.4515; electroweak-corrected; sha256 recorded
in `gcefisher.dm_spectrum`), read here with the same conventions:

    columns  mDM [GeV], Log[10,x] with x = E/m, then one column per channel
    value    dN/dLog10(x) per annihilation  ->  dN/dE = value / (E ln 10)

THE JOINT FIT. The per-bin halo amplitude of `likelihood.fit_all_bins_profile`
is a_k such that counts_k = a_k x J_cells_k. 🚨 J_cells is built with
`flat_spectrum_bin_exposure` = int_bin A(E) dE (E in MeV), so a_k is a
DIFFERENTIAL flux per MeV per steradian per unit J, held flat across the
bin -- not the bin-integrated flux. For annihilating DM the bin-averaged
differential flux is

    a_k = <sigma v> / (8 pi m^2)  x  N_k(m, channel) / dE_k[MeV],
    N_k = int_bin dN/dE dE  (photons per annihilation in the bin),

so a spectrum ties the nine a_k to ONE parameter. The first version (morning
of 2026-09-19) omitted the 1/dE_k and returned <sigma v> a factor ~1e6 too
small at 20 GeV; the injected-signal test could not see it because injection
and fit shared the error. `test_scale_matches_totani` now pins the absolute
scale to his configuration. The flat-within-bin approximation is the same
one every energy-independent template in the stack makes. `joint_fit` profiles the
background normalisations at each <sigma v> on a signed grid (as
`fit_bin_profile` does per bin), sums -2lnL over bins, and returns the best
fit, TS relative to <sigma v> = 0, the 68% interval (delta = 1) and the 95%
upper limit (delta = 2.71 from the minimum, positive side). Signed on
purpose: a negative best fit is a statement about the model, and Foster's
600 GeV point is exactly that.

Units: J in GeV^2 cm^-5 sr^-1 (`gcereduce.dm_templates.j_map`), m in GeV,
<sigma v> in cm^3 s^-1, N_k dimensionless: a_k in cm^-2 s^-1 sr^-1 per
(GeV^2 cm^-5 sr^-1), which is what the cell template's counts basis expects.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np

from ._paths import GCE_FISHER_ROOT

#: `AtProduction_gammas.dat` from http://www.marcocirelli.net/PPPC4DMID.html, kept
#: with the other external inputs under the data root; `$PPPC4DMID_TABLE` overrides.
PPPC_PATH = Path(os.environ.get("PPPC4DMID_TABLE",
                                GCE_FISHER_ROOT / "data" / "external" / "PPPC4DMID" / "AtProduction_gammas.dat"))
#: table column -> our channel name
_COLS = {"b": "b", "c": "c", "q": "q", "t": "t", "W": "W", "Z": "Z", "g": "g", "h": "h",
         "\\[Tau]": "tau", "\\[Mu]": "mu", "e": "e", "\\[Gamma]": "gamma"}
CHANNELS = ("b", "tau", "W", "mu", "c", "q", "t", "Z", "g", "h", "e")
_CACHE: dict = {}


def load_pppc(path=PPPC_PATH) -> dict:
    """{'m': masses [GeV], 'log10x': grid, channel: (n_m, n_x) dN/dlog10x}."""
    key = str(path)
    if key in _CACHE:
        return _CACHE[key]
    with open(path) as f:
        header = f.readline().split()
    data = np.loadtxt(path, skiprows=1)
    cols = {name: i for i, name in enumerate(header)}
    masses = np.unique(data[:, cols["mDM"]])
    lx = np.unique(data[:, cols["Log[10,x]"]])
    out = {"m": masses, "log10x": lx}
    for tab, ours in _COLS.items():
        if tab not in cols:
            continue
        col = data[:, cols[tab]]
        grid = np.empty((masses.size, lx.size))
        for i, m in enumerate(masses):
            sel = data[:, cols["mDM"]] == m
            grid[i] = np.interp(lx, data[sel, cols["Log[10,x]"]], col[sel])
        out[ours] = grid
    _CACHE[key] = out
    return out


def dnde(channel: str, m_gev: float, E_gev, path=PPPC_PATH) -> np.ndarray:
    """dN/dE [1/GeV] per annihilation at photon energies E_gev, for DM mass
    m_gev, interpolated in log m and log10 x; zero above E = m."""
    t = load_pppc(path)
    if channel not in t:
        raise KeyError(f"channel {channel!r}; have {[c for c in CHANNELS if c in t]}")
    m = np.asarray(t["m"]); lx = t["log10x"]; g = t[channel]
    if not (m[0] <= m_gev <= m[-1]):
        raise ValueError(f"m = {m_gev} GeV outside the table ({m[0]}-{m[-1]})")
    j = int(np.clip(np.searchsorted(m, m_gev) - 1, 0, m.size - 2))
    f = (np.log(m_gev) - np.log(m[j])) / (np.log(m[j + 1]) - np.log(m[j]))
    row = (1 - f) * g[j] + f * g[j + 1]
    E = np.asarray(E_gev, float)
    x = np.log10(np.clip(E / m_gev, 1e-300, None))
    val = np.interp(x, lx, row, left=0.0, right=0.0)
    val = np.where(E >= m_gev, 0.0, val)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(E > 0, val / (E * np.log(10.0)), 0.0)
    return out


def bin_yields(channel: str, m_gev: float, edges_gev, n_sub: int = 200, path=PPPC_PATH) -> np.ndarray:
    """N_k = int_bin dN/dE dE, photons per annihilation per analysis bin."""
    edges = np.asarray(edges_gev, float)
    out = np.empty(edges.size - 1)
    for k in range(out.size):
        E = np.logspace(np.log10(edges[k]), np.log10(edges[k + 1]), n_sub)
        out[k] = np.trapezoid(dnde(channel, m_gev, E, path), E)
    return out


def bin_widths_mev(edges_gev) -> np.ndarray:
    e = np.asarray(edges_gev, float) * 1e3
    return e[1:] - e[:-1]


def amplitudes(sigmav, m_gev, channel, edges_gev, path=PPPC_PATH) -> np.ndarray:
    """a_k per bin for <sigma v> [cm^3/s]: sigmav/(8 pi m^2) x N_k / dE_k[MeV] --
    the bin-averaged differential flux per MeV per unit J."""
    return (sigmav / (8.0 * np.pi * m_gev ** 2) * bin_yields(channel, m_gev, edges_gev, path=path)
            / bin_widths_mev(edges_gev))


def joint_fit(stack, J, channel, m_gev, edges_gev, halo_mask=None, n_scan: int = 41,
              sv_max=None, path=PPPC_PATH, with_intervals: bool = True) -> dict:
    """Profile -2lnL over a signed <sigma v> grid with all backgrounds free per
    bin. Returns best fit, TS vs zero, 68% interval, 95% UL and the curve.
    `with_intervals=False` skips the interval searches (None in the output):
    the best fit and TS are unchanged, at roughly half the cost."""
    from .cells import halo_cell_mask
    from .likelihood import _fit_bg_given_offset, _mu_bg, _poisson_val_grad_mu, restrict_stack

    keep = halo_cell_mask() if halo_mask is None else halo_mask
    sr, Jr = restrict_stack(stack, keep, J)
    names = sorted(sr["templates"])
    n_bins = sr["n_bins"]
    N = bin_yields(channel, m_gev, edges_gev, path=path)
    if N.sum() <= 0:                                   # e.g. W+W- below the W mass
        return {"channel": channel, "m_gev": m_gev, "sigmav_best": None, "sigmav_lo68": None,
                "sigmav_hi68": None, "sigmav_ul95": None, "ts_vs_zero": None, "neg2lnl_best": None,
                "bin_yields": N.tolist(), "note": "no photon yield in the analysis bins for this channel/mass"}
    per_bin = []
    for k in range(n_bins):
        T_bg = np.stack([sr["templates"][n][k] for n in names]); n_bin = sr["counts"][k]
        scale_bg = np.maximum(np.array([np.abs(T_bg[i]).sum() for i in range(len(names))]), 1e-300)
        Ts = Jr[k]
        seed = np.full(len(names), np.log10(max(n_bin.sum() / (len(names) + 1), 1e-300)))
        per_bin.append((T_bg / scale_bg[:, None, None], n_bin, Ts, seed))
    coef = N / (8.0 * np.pi * m_gev ** 2) / bin_widths_mev(edges_gev)   # a_k = coef_k x sigmav

    cache = {}

    def total(sv):
        key = float(sv)
        if key not in cache:
            v = 0.0
            for k, (Tb, n_bin, Ts, seed) in enumerate(per_bin):
                off = coef[k] * sv * Ts
                theta, _ = _fit_bg_given_offset(Tb, n_bin, off, seed)
                mu = _mu_bg(theta, Tb) + off
                v += _poisson_val_grad_mu(mu, n_bin)[0]
            cache[key] = v
        return cache[key]

    # scale: the <sigma v> at which the halo alone would supply the counts
    # of the bin where it is relatively brightest; scan symmetrically
    if sv_max is None:
        ratios = [per_bin[k][1].sum() / (coef[k] * per_bin[k][2].sum()) for k in range(n_bins) if coef[k] > 0]
        sv_max = 1.5 * float(np.min(ratios))
    grid = np.linspace(-sv_max / 3.0, sv_max, n_scan)
    vals = np.array([total(s) for s in grid])
    i = int(np.argmin(vals))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, n_scan - 1)]
    for _ in range(3):
        sub = np.linspace(lo, hi, 9); sv = np.array([total(s) for s in sub])
        j = int(np.argmin(sv)); lo, hi = sub[max(j - 1, 0)], sub[min(j + 1, 8)]
    best = 0.5 * (lo + hi); vbest = total(best); v0 = total(0.0)

    def crossing(target, direction):
        a, step = best, abs(best) * 0.1 + sv_max / n_scan
        b = a + direction * step
        while total(b) - vbest < target and abs(b) < 100 * sv_max:
            a, b = b, b + direction * step; step *= 1.5
        for _ in range(30):
            c = 0.5 * (a + b)
            if total(c) - vbest < target:
                a = c
            else:
                b = c
        return 0.5 * (a + b)

    if with_intervals:
        lo68, hi68 = float(crossing(1.0, -1)), float(crossing(1.0, +1))
        ul95 = float(crossing(2.71, +1) if best > 0 else crossing(2.71 + (v0 - vbest), +1))
    else:
        lo68 = hi68 = ul95 = None
    return {"channel": channel, "m_gev": m_gev, "sigmav_best": float(best), "sigmav_lo68": lo68,
            "sigmav_hi68": hi68, "sigmav_ul95": ul95, "ts_vs_zero": float(v0 - vbest),
            "neg2lnl_best": float(vbest), "bin_yields": N.tolist(),
            "curve": {"sigmav": grid.tolist(), "neg2lnl": (vals - vbest).tolist()}}
