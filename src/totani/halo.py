"""The halo signal template — NFW ρ², Totani's baseline (sec 2.3).

🚨 THIS IS THE MODULE THAT CROSSES THE RULE 6 BOUNDARY. Importing it is safe
(it only builds a spatial map); using its output in a fit against REAL data is
what `scripts/hooks/rule6_halo_guard.py` gates, and it is permitted only
because `validation/RULING_real_data.md` carries Kev's signature (2026-09-02).
Nothing here weakens that boundary -- there is still no route from this module
to a component list matching `gcefisher.metric_a`'s forbidden combination, and
this project's own halo/background split (no bulge templates, no GC ROI) was
never what rule 6a's guard was protecting in the first place.

THE PARAMETERS, FROM THE PAPER (sec 2.3), NOT FROM MEMORY:

    "the density profile of the MW dark matter halo is NFW, and use the
    parameters from the Via Lactea II simulation... rho(r) = rho_s x^-1
    (1+x)^-2, where x = r/rs, rs = 21 kpc, and rho_s = 8.1e6 M_sun kpc^-3. The
    distance from the GC to the solar system is set at 8 kpc, and hence
    rho_sun = 0.42 GeV cm^-3."

    "If dark matter is distributed smoothly, the emissivity obeys eps_gamma
    propto rho^2, and the map template in this case is denoted as NFW-rho^2."

These are Totani's Via Lactea II values, NOT the A2020 defaults (gamma=1.2,
rs=26 kpc, rho_sun=0.28, R_sun=8.18 kpc) that `gcereduce.dm_templates.j_map`
falls back to -- every one of the four physical parameters below is passed
explicitly for exactly that reason, and `j_map` itself refuses an unknown
keyword rather than silently using its default, so a typo here fails loudly.

⚠ gamma = 1 REDUCES THE GENERALIZED NFW TO THE STANDARD ONE. `rho_gnfw`'s
sibling `_rho_s_rs3` is documented "for the gamma=1 NFW", confirming
`profile="gnfw", gamma=1.0` is the object Totani calls NFW, not an
approximation of it -- and `j_map` with this profile computes exactly the
line-of-sight integral of rho^2 that "NFW-rho^2" names.
"""
from __future__ import annotations

import numpy as np

__all__ = ["TOTANI_NFW_PARAMS", "nfw_rho2_on_target"]

#: rs = 21 kpc, rho_s = 8.1e6 Msun/kpc^3 -> rho_sun = 0.42 GeV/cm^3 at
#: R_sun = 8 kpc, all stated directly in sec 2.3 (quoted above). gamma=1 is
#: not separately stated as a number -- it IS the definition of "NFW" as
#: opposed to a generalized profile, which the quoted density formula
#: rho = rho_s x^-1 (1+x)^-2 confirms (x^-1 = x^-gamma at gamma=1).
TOTANI_NFW_PARAMS = dict(profile="gnfw", gamma=1.0, rs=21.0,
                         rho_sun=0.42, R_sun=8.0)


def nfw_rho2_on_target(target, s_max: float = 200.0, n_s: int = 800,
                       row_chunk: int = 80, **overrides) -> np.ndarray:
    """int rho^2 ds over the target grid, Totani's NFW-rho^2 halo template.

    Always non-negative (it is a density-squared integral), so the SPATIAL
    shape carries no sign -- the fit's amplitude is what Totani allows to go
    negative (sec 2.3: "fl < 0 is allowed for the halo components"), searching
    symmetrically for excess or deficit. See `likelihood.fit_bin_with_signal`.

    ⚠ `s_max` IS WIDER THAN `gcereduce`'s A2020 DEFAULT (100 kpc), AND
    DELIBERATELY: Totani's ROI reaches 60 deg from the GC, nearly 3x A2020's
    20 deg half-width. 200 kpc comfortably exceeds the profile's own
    r_vir = 402 kpc / 2; truncation error there is negligible for a
    rho^2-weighted integral this centrally concentrated.

    🚨 `j_map` TAKES 1-D COORDINATE VECTORS, NOT THE FULL 2-D MESHGRID.
    The first version passed `templates.grid_lb`'s (ny, nx) output directly.
    `j_map`/`_los_nodes` builds its OWN 2-D grid internally from two 1-D
    vectors via `[None, :]` / `[:, None]` broadcasting (confirmed against
    `tests/test_dm_profile.py`'s own `j_map(GRID, GRID)` with `GRID` a 1-D
    `linspace`) -- so feeding it an ALREADY-2-D array added a SECOND broadcast
    axis on top of the first, silently turning a (960, 960) grid into a
    (960, 960, 960)-shaped intermediate. Measured: RSS reached 13.85 GB and
    climbing on a single row-chunk, on a MACHINE SHARED WITH ANOTHER SESSION'S
    30+ GB queue -- caught and killed before it repeated the 2026-08-20
    OOM/reboot incident this project's own memory notes exist to prevent, not
    after. Fixed at the source: `templates.axis_coords` returns the two 1-D
    vectors `j_map` actually wants.

    Still chunked over `b_c` (now correctly 1-D) purely as a memory-conscious
    default on a shared machine -- with the shape bug fixed, the full grid at
    n_s=800 is a few GB and would likely run in one call, but there is no
    reason to find that boundary experimentally when chunking costs nothing.
    """
    from gcereduce.dm_templates import j_map

    from .templates import axis_coords

    l_c, b_c = axis_coords(target)
    params = dict(TOTANI_NFW_PARAMS)
    params.update(overrides)
    ny = len(b_c)
    out = np.empty((ny, len(l_c)))
    for r0 in range(0, ny, row_chunk):
        r1 = min(r0 + row_chunk, ny)
        out[r0:r1] = j_map(l_c, b_c[r0:r1], s_max=s_max, n_s=n_s, **params)
    return out
