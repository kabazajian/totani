#!/usr/bin/env python
"""Item 2 (Kev, 2026-09-27): HESTIA annihilation J-maps over Totani's region,
to test whether a simulated Milky-Way halo is north/south asymmetric enough to
give the data's S/N ~ 1.5 at 20.8 GeV.

    python scripts/hestia_roi60_maps.py            # 6 halos x 36 views -> ~/Data2/hestia_muru/totani_roi60/

DATA AND ATTRIBUTION. The particle data are the HESTIA Milky-Way analogues of
Muru, Silk, Libeskind, Gottloeber & Hoffmann, PRL 135, 161005 (2025), shared
by M. Muru; the kNN densities (k = 32), the disk frames (stellar angular
momentum, within 1 deg of the AHF minor axis) and the loader are those of the
Comment analysis (~/work/gce-hestia-comment/code/stage3_views.py).

ESTIMATOR. J(l,b) = int rho^2 ds, as sum_i m_i rho_i / s_i^2 per solid angle
(s <= 100 kpc, R_sun = 8.18 kpc, observers every 10 deg of azimuth on the
solar circle in the disk plane, b > 0 along the stellar angular momentum).
Over +/-60 deg a raw histogram is dominated by shot noise from particles near
the observer (one particle at 0.5 kpc outshines a 2-deg cell), so:
  * each particle is spread over its own angular kernel, sigma = 0.45 h_i/s_i
    with h_i = (3 k m_i / 4 pi rho_i)^(1/3) its kNN radius (the rms of a
    uniform ball per axis is 0.45 h); particles are grouped in octaves of
    sigma and each group smoothed (FFT, full sky, 0.25 deg CAR,
    sigma_l = sigma / cos b row by row);
  * particles with s < S_MIN = 1 kpc are replaced by the isotropic local
    term <rho^2> S_MIN, <rho^2> = sum m_i rho_i / V over that ball.
Controls: (a) a SPHERICALIZED copy of each halo (every particle rotated to a
random direction about the halo center, r kept), which has the same radial
profile and the same estimator noise but no asymmetry by construction; (b)
the Comment's +/-20 deg histogram maps (stage3_maps), which this estimator
must reproduce after 1-deg smoothing away from the plane.
Output per halo: jmap [36, 480, 480] (float32, 0.25 deg, target orientation:
row = b increasing from -60, column = l DECREASING from +60), the same for
the spherical control, and per-view N/S sums.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.path.expanduser("~/Data2/hestia_muru"))
os.environ.setdefault("MURU_DATA", str(DATA))
sys.path.insert(0, os.path.expanduser("~/work/gce-hestia-comment/code"))
import stage3_views as S3  # noqa: E402

OUT = DATA / "totani_roi60"
R_SUN, S_MAX, S_MIN, K_NN = 8.18, 100.0, 1.0, 32
STEP = 0.25
NL, NB = int(360 / STEP), int(180 / STEP)
AZ = np.arange(0.0, 360.0, 10.0)
SIG_EDGES = np.array([0.0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 1e9])      # deg
B_CTR = -90 + STEP * (np.arange(NB) + 0.5)
COS_B = np.cos(np.radians(B_CTR))
PIX_SR = np.radians(STEP) ** 2 * COS_B[:, None]
FL = np.fft.rfftfreq(NL); FB = np.fft.rfftfreq(2 * NB)


def smooth(img, sig_deg):
    """Separable Gaussian on the full-sky CAR grid: b (zero-padded) then l (periodic, sigma/cos b)."""
    sb = sig_deg / STEP
    pad = np.zeros((2 * NB, NL)); pad[:NB] = img
    x = np.fft.irfft(np.fft.rfft(pad, axis=0) * np.exp(-2 * np.pi ** 2 * sb ** 2 * FB ** 2)[:, None], n=2 * NB, axis=0)[:NB]
    sl = np.minimum(sb / np.maximum(COS_B, 1e-3), 4 * NL)
    return np.fft.irfft(np.fft.rfft(x, axis=1) * np.exp(-2 * np.pi ** 2 * sl[:, None] ** 2 * FL[None] ** 2), n=NL, axis=1)


def view(pos, mass, rho, h, obs, z_hat):
    """Full-sky J map [NB, NL] (Msun^2 kpc^-5 sr^-1), columns l DECREASING from +180."""
    ex = -obs / np.linalg.norm(obs); ey = np.cross(z_hat, ex); ey /= np.linalg.norm(ey)
    rel = pos - obs; s = np.linalg.norm(rel, axis=1)
    near = s < S_MIN; far = (~near) & (s <= S_MAX)
    local = float((mass[near] * rho[near]).sum() / (4 / 3 * np.pi * S_MIN ** 3)) * S_MIN      # isotropic, per sr
    r, sf = rel[far], s[far]
    l = np.degrees(np.arctan2(r @ ey, r @ ex)); b = np.degrees(np.arcsin(np.clip((r @ z_hat) / sf, -1, 1)))
    w = mass[far] * rho[far] / sf ** 2; sig = np.degrees(0.45 * h[far] / sf)
    flat = np.clip(((b + 90) / STEP).astype(int), 0, NB - 1) * NL + np.clip(((180 - l) / STEP).astype(int), 0, NL - 1)
    out = np.full((NB, NL), local)
    g = np.digitize(sig, SIG_EDGES) - 1
    for k in range(len(SIG_EDGES) - 1):
        sel = g == k
        if not sel.any():
            continue
        img = np.bincount(flat[sel], weights=w[sel], minlength=NB * NL).reshape(NB, NL) / PIX_SR
        out += img if k == 0 else smooth(img, min(float(np.median(sig[sel])), 45.0))
    return out


def crop60(m):
    """Full sky -> |l|,|b| <= 60 at 0.25 deg: rows b from -60, columns l from +60 down."""
    r0 = int(30 / STEP); c0 = int(120 / STEP); n = int(120 / STEP)
    return m[r0:r0 + n, c0:c0 + n].astype(np.float32)


def load(name):
    pos, mass, rho, star, old, z_hat, x_hat, fr = S3.load_halo(name)
    dm = rho > 0
    pos, mass, rho = pos[dm], mass[dm], rho[dm]
    h = (3 * K_NN * mass / (4 * np.pi * rho)) ** (1 / 3)
    return pos, mass, rho, h, z_hat / np.linalg.norm(z_hat), x_hat, fr


def gate_vs_stage3(name, jfull):
    """Our map vs the Comment's +/-20 deg histogram map, both smoothed 1 deg, 5 <= |b| <= 20."""
    from scipy.ndimage import gaussian_filter
    ref = np.load(DATA / "stage3_maps" / f"{name}.npz")["jmap"][0].astype(float)     # [ib, il], l INCREASING, 0.1 deg
    ref = gaussian_filter(ref, 10.0)
    bc = -20 + 0.1 * (np.arange(400) + 0.5)
    # ours: sample at the 0.1-deg centers (nearest 0.25 pixel), smoothed 1 deg
    ours = smooth(jfull, 1.0)
    ib = ((bc + 90) / STEP).astype(int); il = ((180 - bc) / STEP).astype(int)        # l = bc on the l axis
    o = ours[np.ix_(ib, il)]
    sel = (np.abs(bc) >= 5)[:, None] & np.ones(400, bool)[None]
    r = o[sel] / ref[sel]
    return {"median_ratio": float(np.median(r)), "p16": float(np.percentile(r, 16)), "p84": float(np.percentile(r, 84))}


def ns_sums(m60):
    bc = -60 + STEP * (np.arange(480) + 0.5); lc = 60 - STEP * (np.arange(480) + 0.5)
    om = np.radians(STEP) ** 2 * np.cos(np.radians(bc))[:, None]
    roi = (np.abs(bc) >= 10)[:, None] & (np.abs(lc) <= 60)[None]
    N = float((m60 * om)[roi & (bc > 0)[:, None]].sum()); S = float((m60 * om)[roi & (bc < 0)[:, None]].sum())
    return N, S


def test(name="G3.2"):
    t0 = time.time(); pos, mass, rho, h, z_hat, x_hat, fr = load(name); print(f"loaded {len(mass)} dm ({time.time()-t0:.0f}s)", flush=True)
    obs = R_SUN * x_hat; t0 = time.time(); full = view(pos, mass, rho, h, obs, z_hat); print(f"view {time.time()-t0:.1f}s", flush=True)
    print("gate", gate_vs_stage3(name, full)); m = crop60(full); print("S/N", ns_sums(m)[1] / ns_sums(m)[0])
    np.save(OUT / f"test_{name}.npy", m)


def main() -> int:
    OUT.mkdir(exist_ok=True)
    if "--test" in sys.argv:
        test(); return 0
    prov = {"script": "totani/scripts/hestia_roi60_maps.py", "date": "2026-09-27", "R_sun": R_SUN, "s_max": S_MAX, "s_min": S_MIN,
            "k_nn": K_NN, "step_deg": STEP, "azimuths": AZ.tolist(),
            "data": "HESTIA MW analogues of Muru et al., PRL 135, 161005 (2025)", "halos": {}}
    rng = np.random.default_rng(20260927)
    for name in S3.FILES:
        t0 = time.time()
        pos, mass, rho, h, z_hat, x_hat, fr = load(name)
        y_hat = np.cross(z_hat, x_hat)
        u = rng.normal(size=pos.shape); u /= np.linalg.norm(u, axis=1)[:, None]
        pos_sph = u * np.linalg.norm(pos, axis=1)[:, None]
        J = np.zeros((len(AZ), 480, 480), np.float32); Jc = np.zeros_like(J); ns = []
        gate = None
        for i, a in enumerate(np.radians(AZ)):
            obs = R_SUN * (np.cos(a) * x_hat + np.sin(a) * y_hat)
            full = view(pos, mass, rho, h, obs, z_hat)
            if i == 0:
                gate = gate_vs_stage3(name, full)
            J[i] = crop60(full); Jc[i] = crop60(view(pos_sph, mass, rho, h, obs, z_hat))
            ns.append({"az": float(AZ[i]), "N": ns_sums(J[i])[0], "S": ns_sums(J[i])[1], "N_sph": ns_sums(Jc[i])[0], "S_sph": ns_sums(Jc[i])[1]})
        np.savez_compressed(OUT / f"{name}.npz", azimuth_deg=AZ, jmap=J, jmap_spherical=Jc)
        rat = np.array([x["S"] / x["N"] for x in ns]); rc = np.array([x["S_sph"] / x["N_sph"] for x in ns])
        prov["halos"][name] = {"n_dm": int(len(mass)), "gate_vs_stage3": gate, "ns": ns, "elapsed_s": round(time.time() - t0, 1)}
        print(f"{name}: gate {gate}; S/N of J over the ROI median {np.median(rat):.3f} [{rat.min():.3f}, {rat.max():.3f}]; "
              f"spherical control {np.median(rc):.3f} [{rc.min():.3f}, {rc.max():.3f}] ({time.time()-t0:.0f}s)", flush=True)
        (OUT / "PROVENANCE.json").write_text(json.dumps(prov, indent=1))
        del pos, mass, rho, h, pos_sph
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
