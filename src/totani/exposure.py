"""Turning a spatial template into an exposure-weighted per-bin cell basis.

🚨 EVERY TEMPLATE HERE CARRIES ITS OWN FREE NORMALIZATION PER BIN (Totani sec
2.3, and CONVENTIONS.md's point-source note for the same reasoning applied to
point sources). That is what makes the unit questions below tractable rather
than blocking:

  * A2020's delivered HI/H2 templates carry `BUNIT = '1e20 cm-2'` -- GAS COLUMN
    DENSITY, not a gamma-ray flux. Standard practice (and A2020's own) is to
    assume emissivity proportional to column density and let the fit's f_l
    absorb the proportionality constant. Checked 2026-09-03: no BUNIT at all on
    the dust or ICS files.
  * The ICS log10-E integration convention is UNCONFIRMED against a BUNIT
    (recorded already in gce-fisher's a2020_stage2_diffuse.py), and gce-fisher's
    own note there says it plainly: "this affects only the in-bin weighting,
    not the overall scale" -- because f_l is free per bin.

⚠ WHAT WOULD NOT BE HARMLESS: a missing factor that varies ACROSS PIXELS
within one bin distorts the fitted SHAPE, not just its scale, and no free f_l
fixes that. 🚨 AND EXACTLY THAT ERROR WAS PRESENT for the first two days of
this module's life. This docstring originally claimed `regrid_car`'s
`solid_angle=True` had the solid angle covered -- but regridding conserves
flux while CONVERTING BETWEEN GRIDS; it never converts an intensity into
counts. The intensity-to-counts step needs the PIXEL SOLID ANGLE (a factor of
cos b, i.e. 2x across |b| <= 60), and no code applied it anywhere until the
2026-09-04 validation step caught it: the fitted isotropic amplitude jumped
five orders of magnitude between adjacent bins, while the point sources --
already integrated fluxes, needing no solid angle -- fit at f = 0.65-0.96.
`templates.pixel_solid_angle` now supplies it, applied in `totani.stack` to
every intensity-like template. The A2020 path never met this bug because
gtsrcmaps applies solid angle internally. What remains below ARE single global
per-bin constants, which f_l absorbs.

TWO EXPOSURE CONVENTIONS, FOR TWO KINDS OF TEMPLATE
----------------------------------------------------
`flat_spectrum_bin_exposure` -- for ENERGY-INDEPENDENT shape maps (gas, dust,
bubble, Loop I, isotropic): the same spatial map is used at every bin, so
turning it into counts needs SOME assumed within-bin spectral shape to average
the exposure over. FLAT (index 0, dN/dE = const) is `a2020_stage2_fit.py`'s own
stated convention for exactly this case ("the diffuse components were declared
with a FLAT PowerLaw... gtsrcmaps supplies the shape per bin and the fit
supplies the scale"), so this reuses it rather than inventing a different one.

`single_energy_bin_exposure` -- for a cube that is ALREADY bin-integrated under
its own physical spectrum (`ics_on_bins` does this internally via GALPROP's
tabulated energy dependence). Re-doing a flat-spectrum sub-integration on top of
an already-integrated physical prediction would conflate two different
spectral assumptions; a single representative-energy sample (the bin's
geometric mean) is the simpler, non-conflicting choice.
"""
from __future__ import annotations

import numpy as np

__all__ = ["flat_spectrum_bin_exposure", "single_energy_bin_exposure"]


def _log_interp_planes(exposure_cube, exposure_energies_mev, E_mev):
    """Exposure map at scalar energy `E_mev`, log-log interpolated between the
    two bracketing tabulated planes. Vectorized over the whole spatial grid --
    the bracket is the SAME for every pixel, since all pixels share one energy
    axis, so this is one array lerp, not a per-pixel search."""
    log_e = np.log(exposure_energies_mev)
    lE = np.log(E_mev)
    j = np.clip(np.searchsorted(log_e, lE) - 1, 0, len(log_e) - 2)
    f = (lE - log_e[j]) / (log_e[j + 1] - log_e[j])
    lo, hi = np.log(exposure_cube[j]), np.log(exposure_cube[j + 1])
    return np.exp(lo + f * (hi - lo))


def flat_spectrum_bin_exposure(exposure_cube, exposure_energies_mev, edges_gev,
                               n_sub: int = 32) -> np.ndarray:
    """`(n_bins, ny, nx)`: `int_bin A(pixel, E) dE` under an assumed flat
    (dN/dE = const) spectrum, so a shape template's counts contribution is
    `shape[pixel] * this[bin, pixel]`.
    """
    edges_mev = np.asarray(edges_gev, float) * 1000.0
    log_e = np.log(np.asarray(exposure_energies_mev, float))
    if edges_mev.min() < np.exp(log_e[0]) or edges_mev.max() > np.exp(log_e[-1]):
        raise ValueError("requested edges fall outside the tabulated exposure "
                         f"range [{np.exp(log_e[0]):.3g}, {np.exp(log_e[-1]):.3g}] MeV")
    ny, nx = exposure_cube.shape[1:]
    out = np.empty((len(edges_mev) - 1, ny, nx))
    for k in range(len(edges_mev) - 1):
        # 🚨 EXACT TRAPEZOIDAL INTEGRATION, NOT `np.gradient`-WEIGHTED
        # ACCUMULATION. An earlier version used `np.gradient(E)` as quadrature
        # weights (mirroring the pattern in `a2020_stage2_diffuse.ics_on_bins`),
        # and it is NOT exact: `np.gradient`'s boundary handling gave a 3.5%
        # error against a CONSTANT integrand, where the true answer is trivial
        # (A * bin width). `np.trapezoid` along the energy axis is exact for
        # that case and matches the convention already used in
        # `bin_photon_flux`/`source_expected_counts` -- one quadrature
        # convention across the whole module, not two.
        E = np.logspace(np.log10(edges_mev[k]), np.log10(edges_mev[k + 1]), n_sub)
        planes = np.stack([_log_interp_planes(exposure_cube, exposure_energies_mev, Ei)
                           for Ei in E])
        out[k] = np.trapezoid(planes, E, axis=0)
    return out


def single_energy_bin_exposure(exposure_cube, exposure_energies_mev,
                               edges_gev) -> np.ndarray:
    """`(n_bins, ny, nx)`: exposure at each bin's geometric-mean energy."""
    edges_mev = np.asarray(edges_gev, float) * 1000.0
    centers = np.sqrt(edges_mev[:-1] * edges_mev[1:])
    log_e = np.log(np.asarray(exposure_energies_mev, float))
    if centers.min() < np.exp(log_e[0]) or centers.max() > np.exp(log_e[-1]):
        raise ValueError("requested bin centers fall outside the tabulated "
                         f"exposure range [{np.exp(log_e[0]):.3g}, "
                         f"{np.exp(log_e[-1]):.3g}] MeV")
    return np.stack([_log_interp_planes(exposure_cube, exposure_energies_mev, E)
                     for E in centers])
