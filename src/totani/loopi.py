"""Loop I as two uniform-emissivity spherical shells — the only template here
that has no delivered map.

🚨 PARAMETERS, AND WHERE THEY COME FROM (rule 1a: PDF, not a rendered page).
Ackermann et al. 2014, ApJ 793, 64, section 3.1, read from
`gce-fisher/docs/Ackermann_2014_ApJ_793_64.pdf`:

    "The geometric Loop I model assumes synchrotron emission from two shells.
     Each shell is described by five parameters: the center coordinates l, b;
     the distance to the center d; and the inner (r_in) and outer (r_out) radius
     of the shell. The parameters are set to: l1 = 341 deg, b1 = 3 deg,
     d1 = 78 pc, r_in,1 = 62 pc, r_out,1 = 81 pc, l2 = 332 deg, b2 = 37 deg,
     d2 = 95 pc, r_in,2 = 58 pc, r_out,2 = 82 pc."

Ackermann et al. 2017 (ApJ 840, 43) Table 1 names the same object -- "Geometric
template based on radio data (Wolleben 2007)" -- without restating the numbers,
so 2014 is the source. These are Totani's refs [7] and [39], and he notes the
LAT-team values "differ slightly from the values in Wolleben's original paper";
we therefore take the LAT-team numbers, as he did, NOT Wolleben's.

⚠ **THE SUN'S POSITION RELATIVE TO THE TWO SHELLS IS DIFFERENT, AND THE
GEOMETRY CODE MUST HANDLE BOTH.** Shell 1 has d = 78 pc inside r_out = 81 pc, so
the Sun sits INSIDE its outer sphere (and outside the inner one, 62 pc): sight
lines start already within the shell and the near intersection is at s < 0, which
must be clipped to 0. Shell 2 has d = 95 pc beyond r_out = 82 pc, so the Sun is
entirely outside it and both roots are positive. A path-length formula that
assumes one case silently halves the other.

THE MODEL. Uniform emissivity within each shell, so the intensity along a sight
line is proportional to the PATH LENGTH through the shell -- which is why the
two uniform spheres "appear non-uniform with diffuse edges" (2014, section 3.1).
Path length through the shell = (path through the outer sphere) - (path through
the inner sphere), each clipped to s >= 0.

Each shell carries its OWN free normalization in the fit -- 2507.07209 section
2.3: "The normalization factors for each of the two shells are treated as two
independent fit parameters." So the absolute emissivity is never needed and the
returned maps are path lengths in pc.
"""
from __future__ import annotations

import numpy as np

__all__ = ["SHELLS", "sphere_path_length", "shell_map", "loopi_maps"]

#: The two shells, exactly as printed in Ackermann et al. 2014 section 3.1.
#: Longitudes are in the 0-360 convention of the paper (341 deg = -19 deg).
SHELLS = (
    dict(name="loopI_shell1", l_deg=341.0, b_deg=3.0,
         d_pc=78.0, r_in_pc=62.0, r_out_pc=81.0),
    dict(name="loopI_shell2", l_deg=332.0, b_deg=37.0,
         d_pc=95.0, r_in_pc=58.0, r_out_pc=82.0),
)


def sphere_path_length(cos_psi: np.ndarray, d: float, R: float) -> np.ndarray:
    """Path length through a sphere of radius `R` whose centre is at distance
    `d`, along rays from the origin at angle `psi` to the centre.

    The ray is ``s * n_hat`` for ``s >= 0``. Points inside the sphere satisfy
    ``s^2 - 2 s d cos(psi) + d^2 <= R^2``, so

        s± = d cos(psi) ± sqrt(R^2 - d^2 sin^2(psi))

    and the length is ``s+ - max(s-, 0)`` where the discriminant is positive and
    ``s+ > 0``. The ``max(s-, 0)`` is what makes this correct when the origin is
    INSIDE the sphere (shell 1); dropping it loses the near half of every ray.
    """
    cos_psi = np.asarray(cos_psi, dtype=float)
    disc = R**2 - d**2 * (1.0 - cos_psi**2)
    root = np.sqrt(np.clip(disc, 0.0, None))
    s_plus = d * cos_psi + root
    s_minus = d * cos_psi - root
    length = np.clip(s_plus, 0.0, None) - np.clip(s_minus, 0.0, None)
    return np.where(disc > 0.0, length, 0.0)


def shell_map(grid_l_deg, grid_b_deg, l_deg, b_deg, d_pc, r_in_pc, r_out_pc):
    """Path length (pc) through one shell, for each (l, b) on the grid."""
    if not r_in_pc < r_out_pc:
        raise ValueError(f"inner radius {r_in_pc} must be < outer {r_out_pc}")
    lam, beta = np.deg2rad(np.asarray(grid_l_deg, float)), np.deg2rad(np.asarray(grid_b_deg, float))
    lc, bc = np.deg2rad(l_deg), np.deg2rad(b_deg)
    # cos of the angle between each sight line and the shell centre direction.
    cos_psi = (np.cos(beta) * np.cos(bc) * np.cos(lam - lc)
               + np.sin(beta) * np.sin(bc))
    cos_psi = np.clip(cos_psi, -1.0, 1.0)
    return (sphere_path_length(cos_psi, d_pc, r_out_pc)
            - sphere_path_length(cos_psi, d_pc, r_in_pc))


def loopi_maps(grid_l_deg, grid_b_deg) -> dict:
    """Both shells, as `{name: path-length map}`. Two independent templates."""
    return {s["name"]: shell_map(grid_l_deg, grid_b_deg, s["l_deg"], s["b_deg"],
                                 s["d_pc"], s["r_in_pc"], s["r_out_pc"])
            for s in SHELLS}
