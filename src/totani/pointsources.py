"""4FGL-DR4 point sources as per-cell expected counts.

🚨 THIS IS DELIBERATELY NOT HOW THE GC ANALYSIS BUILDS POINT SOURCES, AND THE
DIFFERENCE IS THE WHOLE REASON THIS ROI IS AFFORDABLE. The likelihood lives on
10-degree cells (2507.07209 sec 2.2), which are more than ten times the LAT PSF
above 1 GeV (68% containment 0.85 deg at 1 GeV, 0.16 deg at 10 GeV). So what a
cell needs from a source is its FLUX, not its image: catalog spectrum integrated
over the energy bin, times the exposure at its position, deposited into its cell.
A `gtsrcmaps` run over the 1,425 sources in the halo box at 0.125 deg would be
roughly 70 GB of source maps to produce a number that 1,425 scalars already give.

⚠ WHAT THIS APPROXIMATION COSTS, STATED RATHER THAN HIDDEN. Flux from a source
within roughly one PSF width of a cell boundary belongs partly in the
neighbouring cell, and depositing it whole is wrong by that amount. The affected
fraction is ~2 x 0.85 / 100 deg of cell area at 1 GeV and falls with energy.
Totani makes the same approximation and says so ("PSF errors should have little
effect on the results"). `boundary_fraction` below measures it for a given
catalog rather than leaving it as an assertion.

SPECTRAL FORMS -- VERIFIED, NOT REMEMBERED
------------------------------------------
The three 4FGL parameterizations are implemented below, and the SuperExpCutoff
one in particular is easy to get subtly wrong. They are checked against the
catalog's OWN model-derived columns, `Flux1000` (1-100 GeV photon flux) and
`Energy_Flux100` (100 MeV-100 GeV energy flux), which are integrals of these
same models. Measured 2026-09-02 over 400 sources of each type: `Flux1000`
reproduced to a median 5e-5, `Energy_Flux100` to 2e-4, for all three forms.

⚠ **DO NOT VALIDATE AGAINST `Flux_Band`.** Those are the MEASURED band-by-band
fluxes from a separate fit, not integrals of the global model, so they disagree
with it legitimately -- worst for PowerLaw sources, which are the faint ones
whose bands are noisiest. Tried first, 2026-09-02: median ratio 1.24 with a
16-84% spread of 1.00-2.33, which reads as a broken implementation and is not.
That is a published measurement being used as if it were a definition (rule 1b).
"""
from __future__ import annotations

import numpy as np

__all__ = ["SPECTRUM_TYPES", "dnde", "bin_photon_flux", "cell_edges",
           "cell_index", "extended_source_mask", "boundary_fraction", "misplaced_flux_fraction",
           "MEV_PER_ERG"]

#: MeV per erg. Only needed for the Energy_Flux100 cross-check.
MEV_PER_ERG = 624151.0

SPECTRUM_TYPES = ("PowerLaw", "LogParabola", "PLSuperExpCutoff")


def dnde(row, energy_mev):
    """Differential photon flux (ph cm^-2 s^-1 MeV^-1) at `energy_mev`.

    `row` is one record of the `LAT_Point_Source_Catalog` extension.

        PowerLaw          K (E/E0)^-Gamma
        LogParabola       K (E/E0)^(-alpha - beta ln(E/E0))
        PLSuperExpCutoff  K (E/E0)^(-Gamma_S + d/b) exp[(d/b^2)(1 - (E/E0)^b)]

    The last is 4FGL-DR3/DR4's `PLEC4` form, with Gamma_S = PLEC_IndexS,
    d = PLEC_ExpfactorS and b = PLEC_Exp_Index. All three are verified against
    the catalog's own Flux1000 -- see the module docstring.
    """
    E = np.asarray(energy_mev, dtype=float)
    kind = str(row["SpectrumType"]).strip()
    x = E / float(row["Pivot_Energy"])
    if kind == "PowerLaw":
        return float(row["PL_Flux_Density"]) * x ** (-float(row["PL_Index"]))
    if kind == "LogParabola":
        return float(row["LP_Flux_Density"]) * x ** (
            -float(row["LP_Index"]) - float(row["LP_beta"]) * np.log(x))
    if kind == "PLSuperExpCutoff":
        g = float(row["PLEC_IndexS"])
        d = float(row["PLEC_ExpfactorS"])
        b = float(row["PLEC_Exp_Index"])
        if b == 0.0:
            raise ValueError("PLEC_Exp_Index = 0; the PLEC4 form divides by it")
        return (float(row["PLEC_Flux_Density"]) * x ** (-g + d / b)
                * np.exp(d / b**2 * (1.0 - x**b)))
    raise ValueError(
        f"unknown SpectrumType {kind!r}; expected one of {SPECTRUM_TYPES}. "
        "Refusing to fall back to a power law -- a silently substituted "
        "spectrum is exactly the plausible wrong number rule 1 exists for.")


def bin_photon_flux(rows, edges_gev, n_sub: int = 64) -> np.ndarray:
    """Photon flux (ph cm^-2 s^-1) per source per energy bin.

    Integrated on a log-spaced sub-grid inside each bin, which is what the
    Flux1000 cross-check validates. `n_sub = 64` reproduces Flux1000 to 5e-5.
    """
    edges = np.asarray(edges_gev, float) * 1000.0        # -> MeV
    out = np.empty((len(rows), len(edges) - 1))
    for i, row in enumerate(rows):
        for k in range(len(edges) - 1):
            E = np.logspace(np.log10(edges[k]), np.log10(edges[k + 1]), n_sub)
            out[i, k] = np.trapezoid(dnde(row, E), E)
    return out


def cell_edges(half_deg: float = 60.0, cell_deg: float = 10.0) -> np.ndarray:
    """Cell boundaries along one axis, e.g. -60, -50, ..., 60 for the default."""
    n = 2 * half_deg / cell_deg
    if abs(n - round(n)) > 1e-9:
        raise ValueError(f"{2 * half_deg} deg does not divide into "
                         f"{cell_deg} deg cells")
    return np.linspace(-half_deg, half_deg, int(round(n)) + 1)


def cell_index(l_deg, b_deg, half_deg: float = 60.0, cell_deg: float = 10.0):
    """(row, col) cell of each position; -1 for positions outside the ROI.

    Rows index latitude ascending and columns index LONGITUDE ASCENDING -- which
    is the opposite of the display convention the maps use (CDELT1 < 0). Cells
    are a fit-space bookkeeping index, not an image; keeping them in a plain
    ascending frame avoids inheriting the handedness trap, and nothing here is
    ever plotted without going back through the map WCS.
    """
    l = np.asarray(l_deg, float)
    b = np.asarray(b_deg, float)
    l = (l + 180.0) % 360.0 - 180.0
    e = cell_edges(half_deg, cell_deg)
    col = np.digitize(l, e) - 1
    row = np.digitize(b, e) - 1
    n = len(e) - 1
    bad = (col < 0) | (col >= n) | (row < 0) | (row >= n)
    return np.where(bad, -1, row), np.where(bad, -1, col)


def extended_source_mask(l_grid, b_grid, extended, scale: float = 2.0):
    """Boolean mask of pixels to EXCLUDE around catalogued extended sources.

    Totani sec 2.3: "circular regions with a radius twice the semi-major axis of
    the catalog are masked and excluded from the analysis". `extended` is the
    catalog's `ExtendedSources` extension.
    """
    l = np.deg2rad(np.asarray(l_grid, float))
    b = np.deg2rad(np.asarray(b_grid, float))
    out = np.zeros(np.shape(l_grid), dtype=bool)
    for row in extended:
        lc, bc = np.deg2rad(float(row["GLON"])), np.deg2rad(float(row["GLAT"]))
        radius = scale * float(row["Model_SemiMajor"])
        if not np.isfinite(radius) or radius <= 0:
            continue
        cos_sep = (np.sin(b) * np.sin(bc)
                   + np.cos(b) * np.cos(bc) * np.cos(l - lc))
        out |= np.rad2deg(np.arccos(np.clip(cos_sep, -1, 1))) <= radius
    return out


def boundary_fraction(l_deg, b_deg, psf_deg: float = 0.85,
                      half_deg: float = 60.0, cell_deg: float = 10.0) -> float:
    """Fraction of sources within one PSF width of a cell boundary.

    The size of the approximation this module makes, measured on the actual
    source positions rather than argued from cell geometry.
    """
    l = (np.asarray(l_deg, float) + 180.0) % 360.0 - 180.0
    b = np.asarray(b_deg, float)
    inside = (np.abs(l) <= half_deg) & (np.abs(b) <= half_deg)
    if not inside.any():
        return 0.0
    dl = np.abs((l[inside] + half_deg) % cell_deg - cell_deg / 2.0)
    db = np.abs((b[inside] + half_deg) % cell_deg - cell_deg / 2.0)
    near = (dl > cell_deg / 2.0 - psf_deg) | (db > cell_deg / 2.0 - psf_deg)
    return float(near.mean())


def misplaced_flux_fraction(l_deg, b_deg, psf68_deg: float = 0.85,
                            weights=None, half_deg: float = 60.0,
                            cell_deg: float = 10.0) -> float:
    """Flux-weighted fraction of point-source flux deposited in the WRONG cell.

    🚨 THIS, NOT `boundary_fraction`, IS THE SIZE OF THE APPROXIMATION.
    `boundary_fraction` counts sources within a PSF width of a boundary -- 32.1%
    of the halo-box catalog at 1 GeV -- which sounds alarming and is not the
    error. A source 0.8 deg from a boundary has only the far tail of its PSF on
    the other side; only a source sitting ON the boundary misplaces half its
    flux. What matters is the PSF mass that actually crosses, summed over
    sources and weighted by their flux.

    The PSF is approximated as a 2-D Gaussian with the quoted 68% containment,
    so sigma = psf68 / sqrt(2 ln(1/0.32)) ~ psf68 / 1.51. That UNDERSTATES the
    real King-function tails, so this is a lower bound on the misplacement and
    is quoted as such -- `gcefisher`'s own King fit exists precisely because a
    Gaussian misplaces 5.3% against 0.7% in the lowest gcepy bin. Good enough to
    show the effect is small; not good enough to correct with.
    """
    from math import erfc, log, sqrt

    l = (np.asarray(l_deg, float) + 180.0) % 360.0 - 180.0
    b = np.asarray(b_deg, float)
    inside = (np.abs(l) <= half_deg) & (np.abs(b) <= half_deg)
    if not inside.any():
        return 0.0
    w = (np.ones(inside.sum()) if weights is None
         else np.asarray(weights, float)[inside])
    sigma = psf68_deg / sqrt(2.0 * log(1.0 / 0.32))
    # distance to the nearest cell boundary along each axis
    dl = cell_deg / 2.0 - np.abs((l[inside] + half_deg) % cell_deg - cell_deg / 2.0)
    db = cell_deg / 2.0 - np.abs((b[inside] + half_deg) % cell_deg - cell_deg / 2.0)
    tail = lambda d: 0.5 * np.array([erfc(x / (sigma * sqrt(2.0))) for x in d])
    # two independent axes; double counting the corner is negligible and
    # conservative in the direction of over-stating the error.
    frac = tail(dl) + tail(db)
    return float(np.average(np.clip(frac, 0.0, 1.0), weights=w))


def default_catalog_path() -> "Path":
    """`$FERMI_DATA/catalog/gll_psc_v35.fit`, matching the layout `gcereduce`
    already uses for photons and spacecraft data (CONVENTIONS.md)."""
    import os
    from pathlib import Path

    root = Path(os.environ.get("FERMI_DATA", Path.home() / "data" / "fermi"))
    return root / "catalog" / "gll_psc_v35.fit"


def exposure_pixel(l_deg, b_deg, header) -> tuple:
    """(row, col), 0-based, of the exposure/counts-cube pixel containing
    each source position.

    🚨 COMPUTED FROM THE HEADER DIRECTLY, LIKE `templates.grid_lb`, AND NOT
    THROUGH `astropy.wcs`. Same reasoning: these maps are plain plate carree
    with `CRVAL2 = 0`, well inside astropy's valid domain, so the bug that hit
    the delivered A2020 headers does not apply here -- but going through one
    header-native path for every pixel lookup in this project, rather than
    astropy in some places and headers in others, is what keeps that bug from
    coming back in a different module.
    """
    import numpy as np

    l = (np.asarray(l_deg, float) + 180.0) % 360.0 - 180.0
    b = np.asarray(b_deg, float)
    col = np.round((l - float(header["CRVAL1"])) / float(header["CDELT1"])
                   + float(header["CRPIX1"]) - 1.0).astype(int)
    row = np.round((b - float(header["CRVAL2"])) / float(header["CDELT2"])
                   + float(header["CRPIX2"]) - 1.0).astype(int)
    nx, ny = int(header["NAXIS1"]), int(header["NAXIS2"])
    bad = (col < 0) | (col >= nx) | (row < 0) | (row >= ny)
    if np.any(bad):
        raise ValueError(f"{np.sum(bad)} source position(s) fall outside the "
                         f"{nx}x{ny} exposure grid; filter to the ROI first")
    return row, col


def source_expected_counts(rows, l_deg, b_deg, edges_gev, exposure_cube,
                           exposure_energies_mev, header, n_sub: int = 64):
    """Expected counts per source per energy bin: N_ik = int dN/dE * A(E) dE.

    🚨 THE INTEGRATION, NOT A TWO-POINT AVERAGE. `gcereduce.reduce.
    exposure_at_edges` documents exposure varying by a factor of ~1.6 across
    just the widest analysis bins ("biases a bin by tens of percent... worst in
    the top three bins"), and deliberately provides no bin-averaged reducer for
    that reason -- the correct quantity depends on the spectrum, and here the
    spectrum is known (`dnde`). So A(E) is log-log interpolated across ALL 14
    tabulated exposure planes and evaluated on the SAME log-spaced sub-grid
    `bin_photon_flux` uses for the spectrum, and the two are multiplied and
    integrated together -- not exposure-times-flux at the edges only.

    `exposure_cube` is `(n_planes, ny, nx)`, `exposure_energies_mev` its
    `ENERGIES` table (ascending, one per EDGE -- n_bins + 1 values for n_bins
    analysis bins). `header` must be the SAME WCS as the cube (verified
    byte-identical to the counts-cube header for this reduction, 2026-09-03).
    """
    import numpy as np

    edges_mev = np.asarray(edges_gev, float) * 1000.0
    log_eA = np.log(np.asarray(exposure_energies_mev, float))
    if np.any(edges_mev < np.exp(log_eA[0])) or np.any(edges_mev > np.exp(log_eA[-1])):
        raise ValueError(
            "requested energy edges fall outside the tabulated exposure "
            f"range [{np.exp(log_eA[0]):.3g}, {np.exp(log_eA[-1]):.3g}] MeV; "
            "extrapolating exposure is not done here")

    row, col = exposure_pixel(l_deg, b_deg, header)
    out = np.empty((len(rows), len(edges_mev) - 1))
    for i, src_row in enumerate(rows):
        A_planes = exposure_cube[:, row[i], col[i]].astype(float)
        if np.any(A_planes <= 0):
            raise ValueError(
                f"non-positive exposure at source {i} (row={row[i]}, "
                f"col={col[i]}); a source cannot sit where the LAT never "
                "looked, over a 15-year integration")
        log_A = np.log(A_planes)
        for k in range(len(edges_mev) - 1):
            E = np.logspace(np.log10(edges_mev[k]), np.log10(edges_mev[k + 1]), n_sub)
            A = np.exp(np.interp(np.log(E), log_eA, log_A))
            out[i, k] = np.trapezoid(dnde(src_row, E) * A, E)
    return out


def point_source_cell_template(rows, l_deg, b_deg, edges_gev, exposure_cube,
                               exposure_energies_mev, header,
                               half_deg: float = 60.0, cell_deg: float = 10.0,
                               n_sub: int = 64):
    """The SINGLE point-source template, summed over all sources into cells.

    🚨 ONE TEMPLATE, NOT ONE PER SOURCE (Kev, 2026-09-03, confirmed against
    2507.07209 sec 2.3 and sec 3.4.3 -- see CONVENTIONS.md). Every source's
    expected counts are fixed at its catalog spectrum and summed here; the fit
    gives the RESULT one free normalization per bin, exactly like the other
    known components, not 1,425 independent amplitudes.
    """
    import numpy as np

    counts = source_expected_counts(rows, l_deg, b_deg, edges_gev,
                                    exposure_cube, exposure_energies_mev,
                                    header, n_sub=n_sub)
    r, c = cell_index(l_deg, b_deg, half_deg, cell_deg)
    if np.any(r < 0):
        raise ValueError(f"{(r < 0).sum()} source(s) fall outside the "
                         f"{2*half_deg} deg cell grid; filter to the ROI first")
    n = int(round(2 * half_deg / cell_deg))
    flat = r * n + c
    out = np.empty((counts.shape[1], n, n))
    for k in range(counts.shape[1]):
        out[k] = np.bincount(flat, weights=counts[:, k], minlength=n * n).reshape(n, n)
    return out


#: LAT 68% containment (deg) at 1, 10, 100 GeV -- 2507.07209 sec 2.1, quoted.
PSF68_ANCHORS_GEV = (1.0, 10.0, 100.0)
PSF68_ANCHORS_DEG = (0.85, 0.16, 0.10)

#: 95%/68% containment ratio, anchored on a SOURCED number: Ackermann et al.
#: 2014 (ApJ 793,64) sec 3.2 masks point sources "with a radius of 1.48 deg,
#: which corresponds to the 95% containment region of the PSF at 1 GeV",
#: against Totani's 0.85 deg for 68% at the same energy.
#:
#: ⚠ APPLYING THIS ONE RATIO AT ALL ENERGIES IS AN APPROXIMATION -- the true
#: PSF shape is not exactly self-similar in energy, so the 95%/68% ratio drifts.
#: It is anchored EXACTLY at 1 GeV and is approximate elsewhere. That is
#: adequate here because this radius defines a SYSTEMATIC VARIANT (does masking
#: sources instead of modelling them change the answer?), where the conclusion
#: must be robust to the radius anyway -- and the variant is run at a second,
#: deliberately different radius to demonstrate that. It must NOT be reused
#: anywhere the absolute containment matters; compute it from the IRF with
#: gtpsf for that.
PSF95_OVER_PSF68 = 1.48 / 0.85


def psf95_radius_deg(energy_gev) -> float:
    """Approximate 95% containment radius at `energy_gev`. See the note on
    `PSF95_OVER_PSF68` for what "approximate" means and where not to use it."""
    import numpy as np

    lg = np.log10(np.asarray(energy_gev, float))
    p68 = np.interp(lg, np.log10(PSF68_ANCHORS_GEV), PSF68_ANCHORS_DEG)
    return p68 * PSF95_OVER_PSF68


def point_source_mask(l_grid, b_grid, src_l, src_b, radius_deg) -> "np.ndarray":
    """Boolean mask of pixels within `radius_deg` of any source position.

    Totani sec 3.4.3's alternative point-source treatment: mask the source
    regions and exclude them from the analysis, INSTEAD of subtracting a
    catalog-based model. The two must not be combined -- a masked source's
    photons are gone from the data, so its model must go too.
    """
    import numpy as np

    l_grid = np.asarray(l_grid, float)
    b_grid = np.asarray(b_grid, float)
    out = np.zeros(l_grid.shape, dtype=bool)
    src_l = np.asarray(src_l, float)
    src_b = np.asarray(src_b, float)
    cos_r = np.cos(np.deg2rad(radius_deg))

    # 🚨 BOUNDING BOX PER SOURCE, NOT THE WHOLE GRID. The direct version
    # evaluated a full 960x960 angular separation for each of ~2,750 sources --
    # 2.5 BILLION operations per energy bin, ~27 s each, which timed the
    # systematics run out at 10 minutes. A source's mask disc is at most ~1.3
    # deg across here (a 21x21 pixel box), so all but a vanishing fraction of
    # that work was computing large separations in order to discard them.
    # The box is padded by one pixel and widened in longitude by 1/cos(b),
    # since a disc of fixed angular radius spans more longitude at high
    # latitude; the exact spherical test still runs INSIDE the box, so the
    # result is identical, not approximate.
    ny, nx = l_grid.shape
    b_col = b_grid[:, 0]
    l_row = l_grid[0, :]
    db = abs(b_col[1] - b_col[0]) if ny > 1 else 1.0
    dl = abs(l_row[1] - l_row[0]) if nx > 1 else 1.0

    for i in range(len(src_l)):
        rb = radius_deg + db
        j0 = int(np.searchsorted(b_col, src_b[i] - rb) - 1)
        j1 = int(np.searchsorted(b_col, src_b[i] + rb) + 1)
        j0, j1 = max(j0, 0), min(j1, ny)
        if j0 >= j1:
            continue
        cb = max(np.cos(np.deg2rad(min(abs(src_b[i]) + radius_deg, 89.0))), 1e-3)
        rl = radius_deg / cb + dl
        # l_row DECREASES with column index (CDELT1 < 0), so search the
        # reversed axis rather than assuming ascending order.
        lo, hi = src_l[i] - rl, src_l[i] + rl
        inside_l = (l_row >= lo) & (l_row <= hi)
        if not inside_l.any():
            continue
        cols = np.flatnonzero(inside_l)
        i0, i1 = cols[0], cols[-1] + 1

        sub_l = np.deg2rad(l_grid[j0:j1, i0:i1])
        sub_b = np.deg2rad(b_grid[j0:j1, i0:i1])
        cos_sep = (np.sin(sub_b) * np.sin(np.deg2rad(src_b[i]))
                   + np.cos(sub_b) * np.cos(np.deg2rad(src_b[i]))
                   * np.cos(sub_l - np.deg2rad(src_l[i])))
        out[j0:j1, i0:i1] |= cos_sep >= cos_r
    return out
