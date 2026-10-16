"""Totani's diffuse model: GALPROP v54 `SLZ6R30T150C2` from Ackermann et al.
2012 (ApJ 750, 3), read from the pre-computed output the LAT team published.

🚨 WHICH MODEL, AND WHERE IT COMES FROM. Totani (arXiv:2507.07209, sec 2.3)
does NOT use the Ackermann+ 2017 "Sample Model". He adopts the parameter set
SLZ6R30T150C2 from the 128-model grid of Ackermann+ 2012 -- Lorimer pulsar
source distribution, z_h = 6 kpc, R_h = 30 kpc, T_S = 150 K, E(B-V) cut
2 mag -- "the parameter set used for the ICS component of gll_iem_v07". The
2017 paper's five-gas-ring / three-ICS treatment is his sec 3.4.4 VARIANT
("gas rings & 3 ICS"). The pre-computed GALPROP output of every one of the
128 models is public at

    https://galprop.stanford.edu/PaperIISuppMaterial/   (S{S,L,Y,O}_Z*_R*_T*_C*.tar)

and `SL_Z6_R30_T150_C2.tar` (1.36 GB) was downloaded 2026-09-18 into
`external/galprop_paperII/` (gitignored; README there). His five sec-3.4.4
parameter variants are the same archive family.

WHAT IS IN THE FILES (verified 2026-09-18 on `pi0_decay_healpix_*`):
HEALPix NSIDE 128 RING (0.46 deg), 38 energies 50 MeV -> 814 GeV (factor
1.2997), one `Spectra` row per pixel. The values are **dN/dE in
MeV^-1 cm^-2 s^-1 sr^-1**: E^2 dN/dE at the GC peaks at ~700 MeV at
0.19 MeV cm^-2 s^-1 sr^-1 (the pion bump where it should be) and the local
slope above 1 GeV is -2.66. 🚨 `bremss_healpix` is TOTAL bremsstrahlung,
ionized gas INCLUDED: the 51 per-ring files (HIR, H2R, HII x 17 rings) sum
to it exactly, and the 17 HII ring files alone reproduce
`bremss_ionized_healpix` (which is therefore a SUBSET, 5% of the total, not
an addition). The first version of this module (2026-09-18, used in
`totani_diffuse_2026-09-18` run 1) added `bremss_ionized` on top -- a 5%
double count of bremsstrahlung, ~1% of the gas template -- corrected the
same evening and the run redone. `ics_isotropic_healpix` equals the sum of
its three `comp_{1,2,3}` files (CMB / IR / optical 0.70 / 0.17 / 0.13 at
1 GeV).

THE RINGS (his "gas rings & 3 ICS" variant, sec 3.4.4, after Ackermann+ 2017):
gas split into five Galactocentric ranges, 0-1.5, 1.5-3.5, 3.5-8, 8-10 and
10-50 kpc. GALPROP's 17 annuli (Ackermann+ 2012 Table 2, read from
`docs/Ackermann_2012_ApJ_750_3.pdf` p. 88) map onto them exactly:
ring 1 | 2-5 | 6-12 | 13 | 14-17. Each group is pi0 + bremsstrahlung summed
over the H I, H2 and H II phases, ratio fixed to the model.

HOW TOTANI USES IT (sec 2.3): two templates, each with one free normalisation
per energy bin -- `gas` = pi0 + bremsstrahlung with their ratio fixed to the
model, and `ics`. That replaces this project's 10 A2020 gas/dust maps and
4 IC ring groups (14 templates) with 2. His variant splits gas into five
Galactocentric rings and ICS into its three components.

INTEGRATION ONTO THE ANALYSIS BINS: dN/dE is integrated over each bin exactly
(log-log interpolation between the 38 nodes makes every segment a power law,
integrated analytically -- `integrate_bins`) to a per-bin photon flux per
steradian; the template then enters the
stack through `single_energy_bin_exposure` x solid angle, the same route as
the Pohl IC rings (`exposure.py` says why: the cube is already bin-integrated
under its own spectrum). Nearest-pixel sampling from 0.46 deg to 0.125 deg:
a step function at pixel scale, irrelevant on 10-degree cells.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits

from .healpix import ang2pix_ring

#: The delivered archive's naming stem for Totani's baseline set.
TOTANI_SET = "54_Lorimer_z6kpc_R30kpc_Ts150K_EBV2mag"
_ROOT = Path(__file__).resolve().parents[2] / "external" / "galprop_paperII" / "OUTPUT"

#: Totani sec 2.3: the ICS component of gll_iem_v07 uses this set.
GALDEF_NAME = {"SLZ6R30T150C2": TOTANI_SET}

#: Ackermann+ 2012 Table 2: (R_min, R_max) in kpc of GALPROP annuli 1..17.
RING_BOUNDS_KPC = ((0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0),
                   (4.0, 4.5), (4.5, 5.0), (5.0, 5.5), (5.5, 6.5), (6.5, 7.0), (7.0, 8.0),
                   (8.0, 10.0), (10.0, 11.5), (11.5, 16.5), (16.5, 19.0), (19.0, 50.0))
#: Totani sec 3.4.4 / Ackermann+ 2017: five Galactocentric groups -> annulus numbers.
TOTANI_RING_GROUPS = {"R0-1.5": (1,), "R1.5-3.5": (2, 3, 4, 5), "R3.5-8": (6, 7, 8, 9, 10, 11, 12),
                      "R8-10": (13,), "R10-50": (14, 15, 16, 17)}
GAS_PHASES = ("HIR", "H2R", "HII")


def skymap_path(component: str, model_set: str = TOTANI_SET, root=None) -> Path:
    p = (Path(root) if root is not None else _ROOT) / f"{component}_healpix_{model_set}.gz"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} is missing. It is GALPROP output from Ackermann+ 2012, archive "
            f"SL_Z6_R30_T150_C2.tar at https://galprop.stanford.edu/PaperIISuppMaterial/ ; "
            f"see external/galprop_paperII/README.md.")
    return p


def load_skymap(component: str, model_set: str = TOTANI_SET, root=None,
                allow_negative: bool = False) -> tuple:
    """(dN/dE [npix, n_E] in MeV^-1 cm^-2 s^-1 sr^-1, energies [MeV], nside).

    ⚠ THE PER-RING GAS FILES CONTAIN SMALL NEGATIVE VALUES (92 of the 102,
    measured 2026-09-18): the reddening-corrected annular column densities
    can go negative where the E(B-V) residual is negative (Ackermann+ 2012
    sec 3.3 limits, but does not remove, that). Per file the negative sum is
    at most 0.6% of the positive sum (outer H2 rings; 1e-5 to 1e-4 typical),
    no value is non-finite, and the 17 x 3 ring files still sum to the
    non-negative totals exactly. Totals are loaded strictly; rings with
    `allow_negative=True` and are clipped only AFTER grouping (see
    `galprop_on_target`).
    """
    with fits.open(skymap_path(component, model_set, root)) as h:
        hd = h["SKYMAP"].header
        if hd.get("ORDERING") != "RING" or hd.get("PIXTYPE") != "HEALPIX":
            raise ValueError(f"expected HEALPIX RING, got {hd.get('PIXTYPE')}/{hd.get('ORDERING')}")
        nside = int(hd["NSIDE"])
        spec = np.asarray(h["SKYMAP"].data["Spectra"], float)
        E = np.asarray(h["ENERGIES"].data["MeV"], float)
    if spec.shape != (12 * nside * nside, E.size):
        raise ValueError(f"skymap {spec.shape} does not match NSIDE {nside} x {E.size} energies")
    if not np.isfinite(spec).all():
        raise ValueError(f"{component}: non-finite dN/dE")
    if (spec < 0).any() and not allow_negative:
        raise ValueError(f"{component}: negative dN/dE ({int((spec < 0).sum())} entries); "
                         "only the per-ring gas files are expected to have them")
    return spec, E, nside


def integrate_bins(spec, E_mev, edges_gev) -> np.ndarray:
    """`[npix, n_bins]` photon flux per sr in each analysis bin.

    dN/dE is log-log interpolated between the tabulated nodes, so on every
    segment between consecutive nodes (and at the bin edges) it IS a power
    law A E^g, whose integral is analytic: A (E2^(g+1) - E1^(g+1)) / (g+1),
    or A ln(E2/E1) at g = -1. Summing those over the segments inside a bin is
    therefore EXACT for the interpolant -- no quadrature error at all, which
    a trapezoid on a sub-grid did not achieve for steep spectra (2e-4 at
    index -2.7, found by `tests/test_galprop.py` on 2026-09-18).
    """
    edges_mev = np.asarray(edges_gev, float) * 1e3
    E = np.asarray(E_mev, float)
    if edges_mev.min() < E[0] or edges_mev.max() > E[-1]:
        raise ValueError(f"bin edges {edges_mev.min():.3g}-{edges_mev.max():.3g} MeV outside the "
                         f"tabulated {E[0]:.3g}-{E[-1]:.3g} MeV")
    logE = np.log(E)
    with np.errstate(divide="ignore"):
        logS = np.log(spec)                              # -inf where 0

    def interp(x):
        """log dN/dE at log-energy x (vector), and the segment index."""
        j = np.clip(np.searchsorted(logE, x, side="right") - 1, 0, logE.size - 2)
        f = (x - logE[j]) / (logE[j + 1] - logE[j])
        with np.errstate(invalid="ignore"):          # -inf * 0 where dN/dE = 0
            return logS[:, j] * (1 - f) + logS[:, j + 1] * f

    out = np.empty((spec.shape[0], len(edges_mev) - 1))
    for k in range(len(edges_mev) - 1):
        lo, hi = np.log(edges_mev[k]), np.log(edges_mev[k + 1])
        inner = logE[(logE > lo) & (logE < hi)]
        x = np.concatenate([[lo], inner, [hi]])            # segment boundaries
        s = interp(x)                                     # [npix, n_x]
        x1, x2 = x[:-1], x[1:]
        s1, s2 = s[:, :-1], s[:, 1:]
        g = (s2 - s1) / (x2 - x1)                         # power-law index per segment
        A = np.exp(s1)                                    # dN/dE at the segment start
        E1, E2 = np.exp(x1), np.exp(x2)
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            gp1 = g + 1.0
            pl = A * E1 * ((E2 / E1) ** gp1 - 1.0) / gp1  # = A E1^-g (E2^(g+1)-E1^(g+1))/(g+1)
            ln = A * E1 * np.log(E2 / E1)
            seg = np.where(np.abs(gp1) < 1e-9, ln, pl)
        seg = np.where(np.isfinite(seg), seg, 0.0)          # segments with dN/dE = 0
        out[:, k] = seg.sum(axis=1)
    return out


def galprop_on_target(target, edges_gev, model_set: str = TOTANI_SET, root=None,
                      ics_split: bool = False, gas_rings: bool = False) -> dict:
    """Totani's templates on the target grid, (n_bins, ny, nx) in photons
    cm^-2 s^-1 sr^-1 per bin. Default `{'gas', 'ics'}` (his sec 2.3);
    `ics_split` gives `ics_comp1..3`; `gas_rings` gives the five
    `gas_R*` groups of `TOTANI_RING_GROUPS` (his sec 3.4.4).

    Ring groups are summed in SPECTRUM space (17 x 3 raw ring files, pi0 +
    bremsstrahlung) before anything else, so the small negative entries of
    individual ring files cancel against the others in the group; what
    remains negative is clipped to zero and the clipped fraction of the
    group's flux is recorded in `galprop_on_target.last_clipped`. The bin
    integration (log-log, nonlinear) is then done once per group, so the
    five groups do not sum to the single `gas` template to rounding -- to
    about 1e-3, measured in `tests/test_galprop.py`.
    """
    from .templates import grid_lb

    L, B = grid_lb(target)
    pix = None
    out = {}
    galprop_on_target.last_clipped = {}

    def on_target(spec, E, nside):
        nonlocal pix
        if pix is None:
            pix = ang2pix_ring(nside, L, B)
        binned = integrate_bins(spec, E, edges_gev)             # [npix, n_bins]
        return np.transpose(binned[pix], (2, 0, 1))              # (n_bins, ny, nx)

    def cube(component):
        return on_target(*load_skymap(component, model_set, root))

    if gas_rings:
        for tag, rings in TOTANI_RING_GROUPS.items():
            acc = None
            for k in rings:
                for phase in GAS_PHASES:
                    for comp in ("pi0_decay", "bremss"):
                        spec, E, nside = load_skymap(f"{comp}_{phase}_ring_{k}", model_set, root,
                                                     allow_negative=True)
                        acc = spec if acc is None else acc + spec
            neg = acc < 0
            galprop_on_target.last_clipped[tag] = float(-acc[neg].sum() / acc[~neg].sum())
            acc = np.where(neg, 0.0, acc)
            out[f"gas_{tag}"] = on_target(acc, E, nside)
    else:
        out["gas"] = cube("pi0_decay") + cube("bremss")          # bremss INCLUDES ionized
    if ics_split:
        for i in (1, 2, 3):
            out[f"ics_comp{i}"] = cube(f"ics_isotropic_comp_{i}")
    else:
        out["ics"] = cube("ics_isotropic")
    return out
