"""Galp21 (Youyou Li's GALPROP v57, |b| <= 60 delivery of 2026-09-22) onto the
Totani target grid, in Totani's two component splits.

Each model is ONE `.tar.xz` of ~220 FITS mapcubes (167 GiB uncompressed; this
disk has ~40 GiB free), so nothing is extracted: `ingest_model` streams the
archive through `xz -dc` and `tarfile` mode `r|`, and reduces each member as it
passes. The format and its traps were worked out in gce-fisher on the |b| <= 30
delivery (`gcereduce/galp21.py`, `scripts/galp21_tier1.py`); what is reused:

  * native grid: CAR 3600 x 1200 at 0.1 deg, 44 planes log10(E/MeV) = 2.0..6.3,
    float32, differential intensity cm^-2 s^-1 sr^-1 MeV^-1 (no BUNIT; units
    from the setup doc, corroborated by gce-fisher Tier 2 -- and checked here
    against the 2012 grid, `tests/test_galp21.py`);
  * 🚨 the GC sits ON THE WRAP (CRVAL1 = 180 at CRPIX1 = 1800.5, CDELT1 = +0.1:
    longitude INCREASES with column, opposite to fermitools). Asserted per
    member from the header, never assumed;
  * 17 rings on the CO/H2 grid; Totani's five groups (0-1.5, 1.5-3.5, 3.5-8,
    8-10, 10-50 kpc) land on native edges 1.485, 3.500, 7.991, 10.00 -- the
    same ring indices as the 2012 `TOTANI_RING_GROUPS`, map addition only;
  * H2 rings carry small negative pixels (CO noise x X_CO): as the 2012 path,
    rings are summed in SPECTRUM space per group, then clipped, the clipped
    fraction recorded;
  * inverse Compton: the ANISOTROPIC calculation is primary (gce-fisher ruling:
    the better calculation); the isotropic one is carried as an alternative,
    never added to it. ISRF comp_0/1/2 -> `ics_comp1/2/3`.

Sampling is nearest native pixel for every target pixel (0.125 deg). It is
pointwise, so it commutes with the ring sums, the clip and the log-log bin
integration -- sampling each member first is EXACT, and keeps each spectrum
accumulator at (44, 960, 960).
"""
from __future__ import annotations

import io
import json
import re
import subprocess
import tarfile
import time
from pathlib import Path

import numpy as np

from .galprop import TOTANI_RING_GROUPS, integrate_bins

ARCHIVES = Path.home() / "Data2" / "gce-fisher" / "Galp21_60deg"
NATIVE = (44, 1200, 3600)
CRVAL1, CRPIX1, CDELT1 = 180.0, 1800.5, 0.1
CRVAL2, CRPIX2, CDELT2 = 0.0, 600.5, 0.1
N_RINGS, N_ISRF = 17, 3
GROUP_OF_RING = {k: tag for tag, rings in TOTANI_RING_GROUPS.items() for k in rings}

_GAS = re.compile(r"^(pion_decay|bremss)_(HIR|H2R|HII)_mapcube_comp_(\d+)_57_")
_ICS_RING = re.compile(r"^ics_(anisotropic|skymap)_comp_ring_(\d+)_mapcube_comp_(\d)_57_")
_ICS_TOT = re.compile(r"^ics_(anisotropic|skymap)_comp_mapcube_comp_(\d)_57_")


def code(model: str) -> str:
    """The ranking code of a Galp21 model, e.g. 'G21_SA0_R12'."""
    return f"G21_{model}"


def available(model: str) -> bool:
    """Whether the ingested cubes of a Galp21 model (`SA0_R12`, `SA0_R12+iso`,
    `SA0_R12+hp128`, ...) are in the cache. The public release carries no
    GALPROP-derived file (2026-10-15), so every script that lists Galp21 rows
    filters them through this and runs on the 2012 grid alone when they are
    absent."""
    from . import ranking
    base, _, tag = model.partition("+")
    name = base if tag in ("", "iso") else model
    return (ranking.CACHE / f"galp21_cubes__{name}.npz").is_file()


def native_index(L, B) -> tuple:
    """(row, col) of the native pixel whose centre is nearest each (l, b).
    Native centres: l(c) = 180 + (c - 1799.5) * 0.1, b(r) = (r - 599.5) * 0.1."""
    l = np.mod(np.asarray(L, float), 360.0)
    col = np.floor((l - CRVAL1) / CDELT1 + CRPIX1 - 0.5).astype(int) % NATIVE[2]
    row = np.floor((np.asarray(B, float) - CRVAL2) / CDELT2 + CRPIX2 - 0.5).astype(int)
    if row.min() < 0 or row.max() >= NATIVE[1]:
        raise ValueError("target grid extends beyond |b| <= 60, which the delivery does not cover")
    return row, col


def _check_header(h, name):
    hd = h[0].header
    got = (hd["NAXIS1"], hd["NAXIS2"], hd["NAXIS3"], hd["CRVAL1"], hd["CRPIX1"], hd["CDELT1"],
           hd["CRVAL2"], hd["CRPIX2"], hd["CDELT2"], hd["CTYPE1"].strip(), hd["CTYPE2"].strip())
    want = (NATIVE[2], NATIVE[1], NATIVE[0], CRVAL1, CRPIX1, CDELT1, CRVAL2, CRPIX2, CDELT2, "GLON-CAR", "GLAT-CAR")
    if got[:3] != want[:3] or any(abs(a - b) > 1e-9 for a, b in zip(got[3:9], want[3:9])) or got[9:] != want[9:]:
        raise ValueError(f"{name}: header {got} != expected {want}")
    return 10.0 ** (hd["CRVAL3"] + hd["CDELT3"] * (np.arange(NATIVE[0]) + 1 - hd["CRPIX3"]))


def ingest_model(model: str, L, B, edges_gev, out: Path | None, xz_threads: int = 4, log=print,
                 limit: int | None = None) -> dict:
    """Stream one archive; write `out` (npz, float32 cubes (n_bins, ny, nx) in
    photons cm^-2 s^-1 sr^-1 per bin) and `out.with_suffix('.json')` (manifest).
    `limit`: stop after that many used members and return the manifest only
    (a parser dry run; nothing is written)."""
    from astropy.io import fits

    row, col = native_index(L, B)
    shape = (NATIVE[0],) + L.shape
    acc = {f"gas_{g}": np.zeros(shape) for g in TOTANI_RING_GROUPS}
    for flav in ("aniso", "iso"):
        for i in range(N_ISRF):
            acc[f"{flav}_ics_comp{i + 1}"] = np.zeros(shape)
    ics_tot = {}                                           # delivered per-ISRF totals, for the ring-sum check
    manifest, E_ref, seen = [], None, set()
    t0 = time.time()
    proc = subprocess.Popen(["xz", f"-T{xz_threads}", "-dc", str(ARCHIVES / f"{model}.tar.xz")],
                            stdout=subprocess.PIPE, bufsize=1 << 24)
    try:
        with tarfile.open(fileobj=proc.stdout, mode="r|") as tar:
            for m in tar:
                if not m.isfile():
                    continue
                name = Path(m.name).name
                g, r, t = _GAS.match(name), _ICS_RING.match(name), _ICS_TOT.match(name)
                if not (g or r or t):
                    manifest.append({"name": name, "bytes": m.size, "used": False})
                    continue
                raw = tar.extractfile(m).read()
                with fits.open(io.BytesIO(raw)) as h:
                    E = _check_header(h, name)
                    samp = np.asarray(h[0].data[:, row, col], float)     # (44, ny, nx), the only copy kept
                del raw
                if E_ref is None:
                    E_ref = E
                elif not np.allclose(E, E_ref, rtol=1e-9):
                    raise ValueError(f"{name}: energy grid differs")
                if not np.isfinite(samp).all():
                    raise ValueError(f"{name}: non-finite")
                entry = {"name": name, "bytes": m.size, "used": True, "sum_plane20": float(samp[20].sum()),
                         "n_negative": int((samp < 0).sum())}
                if g:
                    ring = int(g.group(3)) + 1
                    if not 1 <= ring <= N_RINGS:
                        raise ValueError(f"{name}: ring {ring}")
                    acc[f"gas_{GROUP_OF_RING[ring]}"] += samp
                    key = ("gas", g.group(1), g.group(2), ring)
                elif r:
                    flav = "aniso" if r.group(1) == "anisotropic" else "iso"
                    ring, isrf = int(r.group(2)), int(r.group(3))
                    if not (1 <= ring <= N_RINGS and 0 <= isrf < N_ISRF):
                        raise ValueError(f"{name}: ring {ring} / ISRF component {isrf}")
                    if (samp < 0).any():
                        raise ValueError(f"{name}: negative inverse Compton")
                    acc[f"{flav}_ics_comp{isrf + 1}"] += samp
                    key = ("ics", flav, ring, isrf)
                else:
                    flav = "aniso" if t.group(1) == "anisotropic" else "iso"
                    ics_tot[f"{flav}_ics_comp{int(t.group(2)) + 1}"] = samp[20].copy()
                    key = ("ics_total", flav, int(t.group(2)))
                if key in seen:
                    raise ValueError(f"{name}: duplicate member {key}")
                seen.add(key)
                manifest.append(entry)
                log(f"  {model}: {len(seen):3d} {name} ({time.time() - t0:.0f}s)")
                if limit is not None and len(seen) >= limit:
                    break
    finally:
        if limit is not None:
            proc.kill()
        proc.stdout.close()
        rc = proc.wait()
    if limit is not None:
        return {"model": model, "manifest": manifest, "E_mev": E_ref.tolist() if E_ref is not None else None}
    if rc != 0:
        raise RuntimeError(f"xz exited {rc} on {model}")

    # completeness: every ring of every gas phase and process, every IC ring x ISRF component
    need = {("gas", p, ph, k) for p in ("pion_decay", "bremss") for ph in ("HIR", "H2R", "HII")
            for k in range(1, N_RINGS + 1)}
    need |= {("ics", f, k, i) for f in ("aniso", "iso") for k in range(1, N_RINGS + 1) for i in range(N_ISRF)}
    missing = sorted(map(str, need - seen))
    if missing:
        raise ValueError(f"{model}: {len(missing)} members missing, e.g. {missing[:5]}")

    checks = {}
    for k, tot in ics_tot.items():                        # per-ring IC sums reproduce the delivered totals
        checks[f"{k}_ring_sum_over_total_plane20"] = float(acc[k][20].sum() / tot.sum())
    # the 1-gas totals are integrated from the summed SPECTRUM, as the 2012 path
    # integrates its pi0 + bremss totals (log-log integration is not additive)
    acc["gas_total"] = sum(acc[f"gas_{g}"] for g in TOTANI_RING_GROUPS)
    for flav in ("aniso", "iso"):
        acc[f"{flav}_ics_total"] = sum(acc[f"{flav}_ics_comp{i + 1}"] for i in range(N_ISRF))
    cubes, clipped = {}, {}
    for k, spec in acc.items():
        if k.startswith("gas_"):
            neg = spec < 0
            clipped[k] = float(-spec[neg].sum() / spec[~neg].sum()) if neg.any() else 0.0
            spec = np.where(neg, 0.0, spec)
        flat = spec.reshape(NATIVE[0], -1).T
        cubes[k] = integrate_bins(flat, E_ref, edges_gev).T.reshape((-1,) + L.shape).astype(np.float32)
    meta = {"model": model, "archive": str(ARCHIVES / f"{model}.tar.xz"), "n_members": len(manifest),
            "n_used": len(seen), "E_mev": [float(E_ref[0]), float(E_ref[-1]), int(E_ref.size)],
            "gas_clipped_fraction": clipped, "checks": checks, "seconds": time.time() - t0,
            "built": time.strftime("%Y-%m-%d %H:%M")}
    np.savez(out, **cubes)
    out.with_suffix(".json").write_text(json.dumps({**meta, "manifest": manifest}, indent=1))
    return meta


def cubes_for(model: str, variant: str, ics: str = "aniso") -> dict:
    """{name: (n_bins, ny, nx)} in the 2012 path's names: `gas`, `ics` (1-gas)
    or `gas_R*`, `ics_comp1..3` (5-ring)."""
    from . import ranking
    # `<model>+iso` takes the ISOTROPIC inverse-Compton alternative (2026-09-23:
    # the anisotropic IC carries most of the Galp21 misfit); any other
    # `<model>+<tag>` reads a filtered copy from scripts/galp21_smoothing_test.py
    base, _, tag = model.partition("+")
    if tag == "iso":
        ics, model = "iso", base
    z = np.load(ranking.CACHE / f"galp21_cubes__{model}.npz")
    if variant == "gas5_ics3":
        out = {f"gas_{g}": z[f"gas_{g}"] for g in TOTANI_RING_GROUPS}
        out.update({f"ics_comp{i}": z[f"{ics}_ics_comp{i}"] for i in (1, 2, 3)})
    elif variant == "gas1_ics1":
        out = {"gas": z["gas_total"], "ics": z[f"{ics}_ics_total"]}
    else:
        raise ValueError(variant)
    return {k: v.astype(float) for k, v in out.items()}
