"""Stage A of the ranking programme: registries and cell-template caches.

`validation/PLAN_model_ranking_2026-09-19.md`, Kev's decisions of 2026-09-19.
Everything the ranking needs is reduced ONCE to (n_bins, 12, 12) cell
templates on the |b| >= 10 analysis grid and cached; Stages B-D are then
cell-space arithmetic, seconds per fit, parallel across cores, with no
memory pressure on the shared machine.

THREE REGISTRIES

* Diffuse models -- `MODEL_CODES`: the 128 Ackermann+ 2012 GALPROP models
  (`S{S,L,Y,O}_Z{4,6,8,10}_R{20,30}_T{150,100000}_C{2,5}`), plus `a2020`
  (this project's baseline ring model) and `gll_iem_v07` (the LAT team's
  interpolated model; Kev: "include as a labelled reference row" -- it is
  not a physical model and already contains patches for the bubbles and
  Loop I, so structure selection on top of it is not meaningful).
  Each Paper II model comes in two VARIANTS: `gas1_ics1` (Totani's sec 2.3:
  one gas, one ICS template) and `gas5_ics3` (his sec 3.4.4: five
  Galactocentric gas groups, three ICS components).

* Structures -- `STRUCTURES`: name -> builder of a (n_bins, 12, 12) cell
  template. Model-independent ones are cached once; the Totani-residual
  bubble templates are built PER DIFFUSE MODEL (his recipe uses the model's
  own 4.3 GeV residual) and cached per model.

* Shared -- counts, isotropic, point sources, from one stack build.

WHERE THE 128 TARBALLS ARE. Six are here (`external/galprop_paperII/`,
extracted); gce-fisher is fetching all 128 to
`/private/nfs/Data2/gce-fisher/galprop_paperII/` (`scripts/paperII_fetch.py`
there). `paperii_tar(code)` looks locally first, then on NFS, and members are
read from the tar in memory -- the same choice as `gcereduce.paperII` -- so
nothing is downloaded or extracted twice. ⚠ The NFS mount hangs `ls` when it
is busy; the lookup is wrapped in a timeout and reports absence rather than
blocking.
"""
from __future__ import annotations

import gzip
import io
import json
import os
import re
import signal
import tarfile
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

from .cells import bin_to_cells, halo_cell_mask
from .exposure import single_energy_bin_exposure
from .galprop import TOTANI_RING_GROUPS, galprop_on_target, integrate_bins
from .stack import _strip_to_2d, build_background_stack, cell_template_from_shape, halo_cell_template
from .templates import grid_lb, pixel_solid_angle

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs" / "fff2de5bc7e0"
CACHE = ROOT / "results" / "fits" / "ranking_cells"
LOCAL_TARS = ROOT / "external" / "galprop_paperII"
NFS_TARS = Path(os.environ.get("PAPERII_NFS", "/private/nfs/Data2/gce-fisher/galprop_paperII"))

SOURCE_CODES = {"SS": "SNR", "SL": "Lorimer", "SY": "Yusifov", "SO": "OBstars"}
_CODE_RE = re.compile(r"^(S[SLYO])_Z(4|6|8|10)_R(20|30)_T(150|100000)_C(2|5)$")
MODEL_CODES = [f"{s}_Z{z}_R{r}_T{t}_C{c}"
               for s in ("SS", "SL", "SY", "SO") for z in (4, 6, 8, 10)
               for r in (20, 30) for t in (150, 100000) for c in (2, 5)]
assert len(MODEL_CODES) == 128
VARIANTS = ("gas1_ics1", "gas5_ics3")
REFERENCE_MODELS = ("a2020", "gll_iem_v07")


def available(code: str, variant: str = None) -> bool:
    """Whether a row's diffuse inputs are present: always for the 2012 grid
    and the reference models (fetched on demand), only with the ingested
    cubes for a Galp21 row (`G21_*`; see `galp21.available`)."""
    if code.startswith("G21_"):
        from . import galp21
        return galp21.available(code[len("G21_"):])
    return True


def code_to_tag(code: str) -> str:
    """`SL_Z6_R30_T150_C2` -> `54_Lorimer_z6kpc_R30kpc_Ts150K_EBV2mag` (the member stem)."""
    m = _CODE_RE.match(code)
    if not m:
        raise ValueError(f"not a Paper II model code: {code!r}")
    s, z, r, t, c = m.groups()
    return f"54_{SOURCE_CODES[s]}_z{z}kpc_R{r}kpc_Ts{t}K_EBV{c}mag"


def _exists_with_timeout(p: Path, seconds: int = 15) -> bool:
    """`p.exists()` that gives up on a hung NFS mount."""
    def handler(signum, frame):
        raise TimeoutError
    old = signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)
    try:
        return p.exists()
    except TimeoutError:
        return False
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def paperii_tar(code: str) -> Path | None:
    for d in (LOCAL_TARS, NFS_TARS):
        p = d / f"{code}.tar"
        if _exists_with_timeout(p):
            return p
    return None


def paperii_available(code: str) -> bool:
    """Extracted locally, or a tarball is reachable."""
    tag = code_to_tag(code)
    if (LOCAL_TARS / "OUTPUT" / f"pi0_decay_healpix_{tag}.gz").exists():
        return True
    return paperii_tar(code) is not None


class TarSkymaps:
    """Read `OUTPUT/<component>_healpix_<tag>.gz` members from a tarball in
    memory, presenting the same (spec, E, nside) tuple as `galprop.load_skymap`."""

    def __init__(self, tar_path: Path, tag: str):
        self.tar = tarfile.open(tar_path)
        self.tag = tag

    def load(self, component: str, allow_negative: bool = False):
        member = f"OUTPUT/{component}_healpix_{self.tag}.gz"
        f = self.tar.extractfile(member)
        if f is None:
            raise FileNotFoundError(member)
        raw = gzip.decompress(f.read())
        with fits.open(io.BytesIO(raw)) as h:
            hd = h["SKYMAP"].header
            nside = int(hd["NSIDE"])
            spec = np.asarray(h["SKYMAP"].data["Spectra"], float)
            E = np.asarray(h["ENERGIES"].data["MeV"], float)
        if not np.isfinite(spec).all():
            raise ValueError(f"{member}: non-finite")
        if (spec < 0).any() and not allow_negative:
            raise ValueError(f"{member}: negative dN/dE")
        return spec, E, nside

    def close(self):
        self.tar.close()


def _galprop_cubes_from_tar(target, edges_gev, tar_path, tag, ics_split, gas_rings):
    """`galprop.galprop_on_target` with the maps read from a tarball."""
    from .healpix import ang2pix_ring
    from .galprop import GAS_PHASES

    L, B = grid_lb(target)
    src = TarSkymaps(tar_path, tag)
    pix = None
    out = {}

    def on_target(spec, E, nside):
        nonlocal pix
        if pix is None:
            pix = ang2pix_ring(nside, L, B)
        return np.transpose(integrate_bins(spec, E, edges_gev)[pix], (2, 0, 1))

    def cube(component):
        return on_target(*src.load(component))

    try:
        if gas_rings:
            for tag_r, rings in TOTANI_RING_GROUPS.items():
                acc = None
                for k in rings:
                    for phase in GAS_PHASES:
                        for comp in ("pi0_decay", "bremss"):
                            spec, E, nside = src.load(f"{comp}_{phase}_ring_{k}", allow_negative=True)
                            acc = spec if acc is None else acc + spec
                acc = np.where(acc < 0, 0.0, acc)
                out[f"gas_{tag_r}"] = on_target(acc, E, nside)
        else:
            out["gas"] = cube("pi0_decay") + cube("bremss")
        if ics_split:
            for i in (1, 2, 3):
                out[f"ics_comp{i}"] = cube(f"ics_isotropic_comp_{i}")
        else:
            out["ics"] = cube("ics_isotropic")
    finally:
        src.close()
    return out


def _grid_and_basis(edges_gev):
    with fits.open(RUN / "ccube.fits") as h:
        target = _strip_to_2d(h[0].header)
    with fits.open(RUN / "expcube.fits") as h:
        exposure_cube = h[0].data.astype(float)
        exposure_mev = h["ENERGIES"].data["Energy"].astype(float)
    L, B = grid_lb(target)
    omega = pixel_solid_angle(target)
    single_sr = single_energy_bin_exposure(exposure_cube, exposure_mev, edges_gev) * omega[None]
    return target, L, B, single_sr


def _ext_mask(L, B):
    from .pointsources import default_catalog_path, extended_source_mask
    with fits.open(default_catalog_path()) as h:
        return extended_source_mask(L, B, h["ExtendedSources"].data)


# ---- caches -------------------------------------------------------------------

def cache_path(kind: str, name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    return CACHE / f"{kind}__{name}.npz"


def model_cells(code: str, variant: str, edges_gev, force: bool = False) -> dict:
    """(n_bins, 12, 12) diffuse cell templates for one model and variant, cached.
    Returns {name: array} with names prefixed `galprop_`."""
    if variant not in VARIANTS:
        raise ValueError(variant)
    p = cache_path("model", f"{code}__{variant}")
    if p.exists() and not force:
        z = np.load(p)
        return {k: z[k] for k in z.files if k != "_meta"}
    target, L, B, single_sr = _grid_and_basis(edges_gev)
    mask = _ext_mask(L, B)
    if code.startswith("G21_"):
        cubes = galprop_cubes(code, variant, target, edges_gev)
        return _save_model_cells(p, code, variant, cubes, single_sr, L, B, mask)
    tag = code_to_tag(code)
    ics_split, gas_rings = (variant == "gas5_ics3"), (variant == "gas5_ics3")
    # Tarball first: it always holds every member (the extracted OUTPUT/ of the
    # five variant models has only the three totals, which broke gas5_ics3 on
    # the first run, 2026-09-19). Extracted files only when no tarball is found.
    tar = paperii_tar(code)
    if tar is not None:
        cubes = _galprop_cubes_from_tar(target, edges_gev, tar, tag, ics_split, gas_rings)
    elif (LOCAL_TARS / "OUTPUT" / f"pi0_decay_healpix_{tag}.gz").exists():
        cubes = galprop_on_target(target, edges_gev, model_set=tag, ics_split=ics_split, gas_rings=gas_rings)
    else:
        raise FileNotFoundError(f"model {code}: no tarball locally or on NFS, and no extracted maps")
    return _save_model_cells(p, code, variant, cubes, single_sr, L, B, mask, tag=tag)


def _save_model_cells(p, code, variant, cubes, single_sr, L, B, mask, tag=None) -> dict:
    out = {}
    for name, cube in cubes.items():
        pix = cube * single_sr
        out[f"galprop_{name}"] = np.stack([bin_to_cells(pix[k], L, B, mask=mask) for k in range(pix.shape[0])])
    np.savez_compressed(p, _meta=json.dumps({"code": code, "variant": variant, "tag": tag,
                                              "built": time.strftime("%Y-%m-%d %H:%M")}), **out)
    return out


def galprop_cubes(code: str, variant: str, target, edges_gev) -> dict:
    """The diffuse cubes (n_bins, ny, nx) of any GALPROP row: a Paper II code
    (from its tarball) or a Galp21 code `G21_<model>` (from the streamed ingest
    cache, `galp21.cubes_for`)."""
    split = variant == "gas5_ics3"
    if code.startswith("G21_"):
        from . import galp21
        return galp21.cubes_for(code[len("G21_"):], variant)
    tar = paperii_tar(code)
    if tar is None:
        raise FileNotFoundError(f"model {code}: no tarball")
    return _galprop_cubes_from_tar(target, edges_gev, tar, code_to_tag(code), split, split)


def shared_cells(edges_gev, force: bool = False) -> dict:
    """Counts and the model-independent components from one stack build:
    bubble_flat, bubble_structured, loopI_shell1/2, isotropic, point_sources,
    and the halo template J. Cached once."""
    p = cache_path("shared", "base")
    if p.exists() and not force:
        z = np.load(p)
        return {k: z[k] for k in z.files if k != "_meta"}
    bg = build_background_stack(RUN, edges_gev, bubble="ackermann")
    from .templates import bubble_macias_on_target
    target, L, B, _ = _grid_and_basis(edges_gev)
    out = {"counts": bg["counts"], "bubble_flat": bg["templates"]["bubble"],
           "loopI_shell1": bg["templates"]["loopI_shell1"], "loopI_shell2": bg["templates"]["loopI_shell2"],
           "isotropic": bg["templates"]["isotropic"], "point_sources": bg["templates"]["point_sources"],
           "bubble_structured": cell_template_from_shape(RUN, edges_gev, bubble_macias_on_target(target)),
           "J": halo_cell_template(RUN, edges_gev)}
    # the A2020 diffuse templates ride along as the `a2020` reference model
    for n in ("HI_ring1", "HI_ring2", "HI_ring3", "HI_ring4", "H2_ring1", "H2_ring2", "H2_ring3", "H2_ring4",
              "dust_pos", "dust_neg", "ics_ring123", "ics_ring4", "ics_ring5", "ics_ring6"):
        out[f"a2020_{n}"] = bg["templates"][n]
    np.savez_compressed(p, _meta=json.dumps({"built": time.strftime("%Y-%m-%d %H:%M")}), **out)
    return out


def structure_cells(name: str, edges_gev, force: bool = False) -> np.ndarray:
    """One model-independent structure template in cell space, cached."""
    from . import structures
    p = cache_path("structure", name)
    if p.exists() and not force:
        return np.load(p)["cells"]
    builder = structures.STRUCTURES[name]
    target, L, B, _ = _grid_and_basis(edges_gev)
    arr = builder(target, edges_gev)
    if arr.ndim == 2:                                   # a spatial shape -> through the flat-spectrum route
        cells = cell_template_from_shape(RUN, edges_gev, arr)
    else:                                               # already a per-bin intensity cube -> single-energy route
        _, _, _, single_sr = _grid_and_basis(edges_gev)
        mask = _ext_mask(L, B)
        pix = arr * single_sr
        cells = np.stack([bin_to_cells(pix[k], L, B, mask=mask) for k in range(pix.shape[0])])
    np.savez_compressed(p, cells=cells, _meta=json.dumps({"name": name, "built": time.strftime("%Y-%m-%d %H:%M")}))
    return cells


# ---- assembling and scoring ---------------------------------------------------

def assemble(model: str, variant: str, structures: tuple, edges_gev) -> dict:
    """A stack dict for `model` (a Paper II code, 'a2020' or 'gll_iem_v07')
    with the named structures plus isotropic and point sources."""
    sh = shared_cells(edges_gev)
    t = {"isotropic": sh["isotropic"], "point_sources": sh["point_sources"]}
    if model == "a2020":
        t.update({k[len("a2020_"):]: v for k, v in sh.items() if k.startswith("a2020_")})
    elif model == "gll_iem_v07":
        t["giem"] = structure_cells("gll_iem_v07", edges_gev)
    else:
        t.update(model_cells(model, variant, edges_gev))
    for s in structures:
        if s in sh:
            t[s] = sh[s]
        else:
            t[s] = structure_cells(s, edges_gev)
    return {"counts": sh["counts"], "n_bins": sh["counts"].shape[0], "templates": t}


def score(stack, J, with_uncertainty: bool = True, with_injection: bool = True) -> dict:
    """The standard scoring (as `scripts/totani_diffuse.py`), in the library.
    `with_injection=False` skips the injection-slope gate (Stage E scores
    simulated skies, where the truth is known and the gate is not needed)."""
    from .likelihood import fit_all_bins, fit_all_bins_profile, injection_slope, quadrant_masks, restrict_stack

    keep = halo_cell_mask()
    n_bins = stack["n_bins"]
    sr, Jr = restrict_stack(stack, keep, J)
    bg = fit_all_bins(sr, halo_mask=keep)["per_bin"]
    fit = fit_all_bins_profile(sr, Jr, halo_mask=keep, with_uncertainty=with_uncertainty)["per_bin"]
    roi_J = [float(Jr[k][keep].sum()) for k in range(n_bins)]
    d = {"neg2lnl_bg_only": [bg[k]["neg2lnl"] for k in range(n_bins)],
         "mean_abs_sigma_halo": [bg[k]["mean_abs_sigma_halo"] for k in range(n_bins)],
         "ts": [bg[k]["neg2lnl"] - fit[k]["neg2lnl"] for k in range(n_bins)],
         "halo_counts": [fit[k]["signal_amplitude"] * roi_J[k] for k in range(n_bins)],
         "n_free_per_bin": len(stack["templates"]), "components": sorted(stack["templates"])}
    if with_uncertainty:
        d["halo_lo"] = [fit[k]["signal_amplitude_lo"] * roi_J[k] for k in range(n_bins)]
        d["halo_hi"] = [fit[k]["signal_amplitude_hi"] * roi_J[k] for k in range(n_bins)]
    if with_injection:
        d["injection_slopes"] = injection_slope(sr, Jr, [0.0, 100e-35, 300e-35], halo_mask=keep)["slopes"]
    quads = {}
    for qn, qm in quadrant_masks().items():
        stq, Jq = restrict_stack(stack, qm, J)
        rq = fit_all_bins_profile(stq, Jq, halo_mask=qm, with_uncertainty=False)["per_bin"]
        quads[qn] = [rq[k]["signal_amplitude"] * float(Jq[k][qm].sum()) for k in range(n_bins)]
    d["quadrant_counts"] = quads
    d["quadrant_sum_over_joint"] = [
        (sum(quads[q][k] for q in quads) / d["halo_counts"][k]) if d["halo_counts"][k] else None
        for k in range(n_bins)]
    d["quadrant_any_negative"] = [bool(any(quads[q][k] < 0 for q in quads)) for k in range(n_bins)]
    north = [quads["l+b+"][k] + quads["l-b+"][k] for k in range(n_bins)]
    south = [quads["l+b-"][k] + quads["l-b-"][k] for k in range(n_bins)]
    d["halo_north"], d["halo_south"] = north, south
    return d
