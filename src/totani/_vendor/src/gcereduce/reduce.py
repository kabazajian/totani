"""Parameterized Fermi-LAT reduction: every data-selection choice is a config field.

Task 6 (Track B). This is the capability gcepy cannot provide -- gcepy ships one
frozen reduction (the Cholis et al. 2112.09706 selection), so no data-selection
axis can be varied from inside it. Everything the Phase-3 blind optimization
wants to scan -- ROI, energy range and binning, pixel size, event class and
type, zenith cut, time span -- is a field of `ReductionConfig` here.

RUNS UNDER THE FERMITOOLS ENVIRONMENT, NOT THE JAX ONE. The two environments
are deliberately separate and talk only through files:

    conda run -n fermi-2.5.3 python -c "from gcefisher.reduce import ..."

(The tools are now launched via `conda run` internally, so callers in the jax
environment need not do this themselves -- see `_run`.)

Nothing in this module imports jax, gcepy, or anything else from `gcefisher`,
so importing it cannot perturb the environment the frozen reference values were
produced in.

Photon data lives OUTSIDE the repository at ``$FERMI_DATA`` (default
``~/data/fermi``) — see `scripts/fetch_weekly.py`. The raw all-sky weekly files
are kept as the generally-usable archive; the 45° GC cone is a derived,
regenerable product.

Caching: every product lands in `runs/<config_hash>/` and is skipped if already
present and readable, so an interrupted pipeline resumes by re-invocation. The
livetime cube is the expensive step and depends only on the time span and the
zenith cut, so it is cached under its own narrower key and shared across every
config that agrees on those two -- which is most of the grid.
"""

from __future__ import annotations

import atexit
import dataclasses
import functools
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

# The Cholis et al. (2112.09706) analysis binning, in GeV.
#
# 🚨 THESE ARE THE EXACT EDGES from the Zenodo release README's explicit
# "#label under 14-ebin setup" table -- NOT the 3-significant-figure values of
# Table III, which is what this project recorded until 2026-07-30 and which are
# WRONG AT A LEVEL THAT MATTERS. The native grid is 38 bins spaced by exactly
# 1.2996566 in energy starting at 0.0438587 GeV; analysis bins 0-10 are native
# 7-17, and bins 11/12/13 group native 18,19,20 / 21,22,23 / 24,25,26.
#
# Using the rounded Table III edges instead misplaces counts by up to 1.35% in
# a single bin (bin 5, whose edge 1323.82 MeV rounds to 1320) while shifting
# the TOTAL by only 0.081% -- i.e. the error is nearly invisible in the totals
# and over 1% per bin. That is the same shape of trap as the two longitude
# flips: a check on the sum cannot see it. See CONVENTIONS.md.
CHOLIS_EDGES_GEV = (
    0.274698, 0.357014, 0.463995, 0.603034, 0.783737, 1.01859, 1.32382,
    1.72051, 2.23607, 2.90612, 3.77696, 4.90875, 10.776, 23.6561, 51.9312,
)

# The rounded Table III values, kept ONLY so the comparison above is
# reproducible and so nobody "restores" them thinking they are the real edges.
CHOLIS_EDGES_GEV_ROUNDED_TABLE3 = (
    0.275, 0.357, 0.464, 0.603, 0.784, 1.02, 1.32, 1.72,
    2.24, 2.91, 3.78, 4.91, 10.8, 23.7, 51.9,
)

#: The lower edge of the Pohl et al. (2022) inverse-Compton grid, in GeV.
#: Read from the file WCS (`CRVAL3 = 2.69897 = log10 500`), and it is a LOWER
#: BIN EDGE rather than a bin centre -- Oscar Macias, asked directly, 2026-08-03.
#: Below it the ring family has no IC information and none will be extrapolated.
ICS_FLOOR_GEV = 0.5

#: AMENDMENT 12 (2026-08-11): the analysis grid, 12 bins from the ICS floor.
#:
#: 🚨 DERIVED, NEVER TYPED. The amendment is "insert an edge at 500 MeV and keep
#: every existing Cholis edge above it", so it is computed from
#: `CHOLIS_EDGES_GEV` rather than transcribed. Hand-typing it to six figures
#: introduced a 1.3e-06 discrepancy that made the Phase 2 external-error gate
#: appear to fail; built this way the eleven unchanged bins reproduce their
#: published sigma BIT-IDENTICALLY (0.00e+00), which is what the gate is for.
AMENDMENT12_EDGES_GEV = (ICS_FLOOR_GEV,) + tuple(
    e for e in CHOLIS_EDGES_GEV if e > ICS_FLOOR_GEV)

#: The Cholis native energy grid, exactly as the Zenodo release README states
#: it: 38 bins spaced by a factor 1.2996566, the first starting at 0.0438587
#: GeV. The published 14-bin analysis grid is native k = 7..18 singly, then
#: three-bin merges at k = 21, 24, 27 -- reproduced from this to 2.7e-06, which
#: is what licenses extending the SAME construction rather than inventing one.
CHOLIS_NATIVE_E0_GEV = 0.0438587
CHOLIS_NATIVE_RATIO = 1.2996566


def cholis_native_edge(k: int) -> float:
    """Native edge `k` of the Cholis grid, in GeV. Derived, never transcribed."""
    return CHOLIS_NATIVE_E0_GEV * CHOLIS_NATIVE_RATIO ** k


#: AMENDMENT 33's floor and ceiling, and they are NOT preferences -- they are the
#: endpoints of the published external-error coverage, below/above which a sigma
#: would have to be invented (rule 1).
#:
#:   floor   0.1 GeV   -- both published tables begin at exactly 100 MeV
#:   ceiling 409.6 GeV -- the BUBBLES table (arXiv:1407.7905) ends there. Usable
#:                        IGRB runs to 579.26 GeV, but 409.6 is the highest
#:                        energy at which EVERY external component has a
#:                        published sigma.
#:
#: 🚨 The 579.26-819.20 GeV IGRB band is an UPPER LIMIT (fIGRB 4.206e-14 against
#: foreground offsets ~2.5e-12, fractional sigma 43.6) and is excluded by
#: `published_bands_mrt`. The paper's TITLE says 820 GeV; the usable measurement
#: stops at 579.26.
A33_FLOOR_GEV = 0.1
A33_CEILING_GEV = 409.6

#: 🚨 THE SECOND FLOOR -- THE INSTRUMENTAL ONE, WHICH AMENDMENT 31's COVERAGE
#: TABLE HAS NO ROW FOR. A31 located the floor of the published INPUTS (100 MeV,
#: where both external-error tables begin). The LAT PSF is not an input whose
#: coverage runs out; it degrades continuously, and measured on `baseline_a33`'s
#: own gtpsf output it passes the registered `r_free = 5 deg` at 274.698 MeV:
#:
#:     E_lo      r68    r90    r95
#:     100.0 MeV 4.93   8.83   10.82   <- also breaches the 10 deg source margin
#:     125.1 MeV 4.21   7.72    9.68
#:     162.6 MeV 3.51   6.66    8.53
#:     211.4 MeV 2.87   5.57    7.24
#:     274.7 MeV 2.37   4.71    6.19   <- clean, and so is everything above
#:
#: ⚠ This value is NOT tuned to the PSF: it is Cholis native edge k = 7, the
#: published table's own lower limit, which is simply where the measurement
#: happens to land. Nothing was fitted to produce it.
#:
#: 🚫 This is NOT an adopted floor. The low edge is 100 MeV by PI ruling (A31,
#: A33). This constant exists so the contingency configuration in
#: `grid_driver.prefreeze` is derived rather than typed, and so the number is
#: findable next to the range it qualifies.
#: Record: `validation/psf_vs_ps_handling_2026-08-18.md`.
A33_PSF_FLOOR_GEV = 0.274698

#: AMENDMENT 33 (2026-08-18): the analysis grid, 21 bins over 100 MeV - 409.6 GeV.
#:
#: 🚨 DERIVED, NEVER TYPED -- the same provision that survives from AMENDMENT 12,
#: and for the same reason: hand-typing A12's edges to six figures introduced a
#: 1.3e-06 discrepancy that made the external-error gate appear to fail.
#:
#: THE CONSTRUCTION, stated so it is reproducible rather than recognisable:
#:
#:   1. keep `CHOLIS_EDGES_GEV` -- the Zenodo release README's own explicit edge
#:      table -- entirely and unchanged;
#:   2. below it, extend with native edges k = 4, 5, 6 (dlog10 = 0.1138, the
#:      single-native spacing the table already uses at its own low end);
#:   3. above it, continue the table's own three-bin merge pattern with native
#:      k = 30, 33 (dlog10 = 0.3415, where the merges already were);
#:   4. insert an edge at the floor and at the ceiling.
#:
#: 🚨 WHICH SOURCE GOVERNS WHERE, AND WHY IT MATTERS. The README states the grid
#: BOTH ways -- as an explicit edge table and as "38 bins spaced by 1.2996566
#: from 0.0438587 GeV" -- and the two disagree at 2.7e-06 relative, one of them
#: being rounded. The TABLE is the published quantity, so it governs everywhere
#: it reaches; the formula is used ONLY to extend past its ends, where there is
#: no tabulated value to prefer. Reconstructing the tabulated region from the
#: formula instead would shift 15 edges by 2.7e-06 -- and a 1.3e-06 shift of
#: exactly this kind is what once made the Phase 2 external-error gate appear to
#: fail. The formula reproducing the table to 2.7e-06 is the WARRANT for the
#: extension, not a substitute for the table.
#:
#: This is AMENDMENT 12's own construction with the floor moved 500 -> 100 MeV
#: (A32 removed the Pohl ICS floor by substituting v57 SA50) and a ceiling added.
#:
#: ⚠ WHAT MOVES AND WHAT DOES NOT. 11 of A12's 12 bins survive BIT-IDENTICALLY
#: (0.603034 GeV upward), so prior results remain comparable bin by bin. The one
#: that moves is A12's bin 0, 0.5-0.603034 GeV: 0.5 was the ICS floor and
#: nothing else, so with the floor gone it becomes the native 0.463995-0.603034.
#: Four bins are added below and three above.
#:
#: ⚠ The top bin (250.26-409.60 GeV) and the bottom bin (0.100-0.125 GeV) are
#: TRUNCATED native bins -- 0.2140 and 0.0974 dex against the regular 0.3415 and
#: 0.1138. That is deliberate and is the same treatment A12 gave its floor bin;
#: the alternative is to discard published coverage.
AMENDMENT33_EDGES_GEV = (
    (A33_FLOOR_GEV,)
    # below the published table: native edges k = 4, 5, 6 above the floor
    + tuple(cholis_native_edge(k) for k in (4, 5, 6)
            if cholis_native_edge(k) > A33_FLOOR_GEV)
    # the published table itself, unchanged and bit-identical
    + tuple(CHOLIS_EDGES_GEV)
    # above it: the merge pattern continued, native k = 30, 33
    + tuple(cholis_native_edge(k) for k in (30, 33)
            if cholis_native_edge(k) < A33_CEILING_GEV)
    + (A33_CEILING_GEV,)
)

# Galactic Centre in J2000, used for the gtselect cone only. The analysis grid
# is defined in GALACTIC coordinates centred on (l, b) = (0, 0).
GC_RA, GC_DEC = 266.404, -28.936


@dataclass(frozen=True)
class ReductionConfig:
    """One point on the Phase-3 grid.

    Defaults reproduce the Cholis/gcepy selection, which is the anchor Task 6
    must hit -- NOT a recommendation for the analysis. The blind optimization
    (Task 8) chooses the analysis config; this default exists so the pipeline
    can be validated against a reduction built independently of ours.
    """

    # --- time -------------------------------------------------------------
    week_min: int = 9
    week_max: int = 670

    # --- energy -----------------------------------------------------------
    # Explicit edges (GeV). A grid point that wants N log-spaced bins between
    # two limits should build the edges and pass them, so that what was
    # actually used is recorded rather than re-derived.
    energy_edges_gev: tuple[float, ...] = CHOLIS_EDGES_GEV

    # --- geometry ---------------------------------------------------------
    roi_half_deg: float = 20.0        # |l| <= this and |b| <= this
    pixel_deg: float = 0.1
    coordsys: str = "GAL"
    proj: str = "CAR"

    # --- event selection --------------------------------------------------
    evclass: int = 256                # P8R3_CLEAN
    evtype: int = 1                   # FRONT
    zmax_deg: float = 100.0
    # Cholis et al. 2112.09706 §2. Note DATA_QUAL==1, not the more common
    # DATA_QUAL>0 -- recorded as the paper states it.
    gtmktime_filter: str = "DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52"
    roicut: str = "no"

    #: Override the IRF that `IRF_BY_EVCLASS` would choose. Needed to reproduce
    #: a published analysis on ITS instrument response: A2020 used
    #: `P8R3_ULTRACLEANVETO_V2` while our default for that event class is V3.
    #: Leave None for our own analysis -- V3 is the current recommendation.
    irf_override: str | None = None

    # --- provenance -------------------------------------------------------
    note: str = ""

    @property
    def npix(self) -> int:
        n = 2 * self.roi_half_deg / self.pixel_deg
        if abs(n - round(n)) > 1e-9:
            raise ValueError(
                f"ROI {2 * self.roi_half_deg} deg does not divide into "
                f"{self.pixel_deg} deg pixels ({n} pixels)"
            )
        return int(round(n))

    @property
    def n_ebins(self) -> int:
        return len(self.energy_edges_gev) - 1

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @property
    def hash(self) -> str:
        """Short, stable hash over everything that changes the products.

        `note` is excluded: it is documentation, and changing it must not
        invalidate a cached reduction.
        """
        d = {k: v for k, v in self.to_dict().items() if k != "note"}
        # 🚨 A FIELD LEFT UNSET MUST NOT PERTURB EXISTING HASHES. `irf_override`
        # was added 2026-08-08 for the A2020 reproduction; including it
        # unconditionally would have re-keyed EVERY config, orphaning the
        # 102-minute baseline reduction and the Cholis anchor that Task 6
        # validates against. Dropping it when None keeps old hashes bit-identical
        # while still recording it whenever it is actually used.
        if d.get("irf_override") is None:
            d.pop("irf_override", None)
        blob = json.dumps(d, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode()).hexdigest()[:12]

    @property
    def ltcube_key(self) -> str:
        """Livetime depends only on the time span and the zenith cut.

        gtltcube is the slow step, so configs that differ only in energy,
        geometry or event class must share one. Keeping this key narrow is
        what makes the Phase-3 grid affordable.
        """
        blob = json.dumps(
            {"week_min": self.week_min, "week_max": self.week_max,
             "zmax_deg": self.zmax_deg, "filter": self.gtmktime_filter},
            sort_keys=True, separators=(",", ":"),
        )
        return hashlib.sha256(blob.encode()).hexdigest()[:12]


#: The conda environment holding fermitools. CONVENTIONS.md: it is DELIBERATELY
#: separate from the jax/gcefisher environment, and the two communicate through
#: files on disk, never in-process.
#:
#: 🚨 `fermi-2.5.3`, NOT `fermi`. CONVENTIONS.md: "Fermitools 2.5.2 is BANNED --
#: it silently mispositions extended sources." The Xerces->RapidXML migration
#: dropped `SpatialFunction::update()`, so a RadialGaussian or RadialDisk source
#: never initialises its centre -- exactly the spatial models 4FGL uses for
#: extended sources. Measured here: on 2.5.2 the map is wrong by 4e8 in
#: normalisation and displaced 28 deg to the ROI corner.
FERMI_ENV = os.environ.get("GCE_FERMI_ENV", "fermi-2.5.3")

#: fermitools versions that must never generate products. See above.
BANNED_FERMITOOLS = ("2.5.2",)
_VERSION_CHECKED: set[str] = set()


def _conda() -> str:
    """Path to the conda executable, or a message explaining what to install.

    🚨 The ST tools must be launched through `conda run -n $GCE_FERMI_ENV`
    (default `fermi-2.5.3`), NOT by calling
    the binary directly: they need the environment's activation hooks, which set
    `PFILES` and `FERMI_DIR`. Calling the binary bare fails with
    `ParFileError: .par file not found`, a message that reads like a broken
    install rather than an unactivated environment (CONVENTIONS.md).
    """
    for c in (os.environ.get("CONDA_EXE"),
              "/opt/homebrew/Caskroom/miniforge/base/bin/conda",
              shutil.which("conda"), shutil.which("mamba")):
        if c and Path(c).exists():
            return c
    raise RuntimeError(
        "conda not found, so the Fermi ST tools cannot be launched. The tools "
        f"live in the {FERMI_ENV!r} environment; create it with "
        "`conda env create -f env/fermitools.yml`. Set CONDA_EXE or "
        "GCE_FERMI_ENV to override.")


def _assert_fermitools_ok() -> None:
    """🚨 Refuse to run against a banned fermitools build. Checked once per process.

    A version guard rather than a comment, because the 2.5.2 failure is SILENT:
    the tools exit 0 and write a map that is displaced 28 degrees. Nothing
    downstream would notice. See CONVENTIONS.md "Fermitools 2.5.2 is BANNED".
    """
    if FERMI_ENV in _VERSION_CHECKED:
        return
    proc = subprocess.run([_conda(), "list", "-n", FERMI_ENV, "^fermitools$"],
                          capture_output=True, text=True)
    line = [l for l in proc.stdout.splitlines()
            if l.strip().startswith("fermitools ")]
    if not line:
        raise RuntimeError(
            f"cannot determine the fermitools version in env {FERMI_ENV!r}; "
            f"refusing to run rather than risk a banned build")
    version = line[0].split()[1]
    if version in BANNED_FERMITOOLS:
        raise RuntimeError(
            f"fermitools {version} in env {FERMI_ENV!r} is BANNED: it silently "
            f"mispositions extended sources (RadialGaussian/RadialDisk never "
            f"initialise their centre; measured 4e8 wrong in normalisation and "
            f"28 deg displaced). Use fermi-2.5.3. See CONVENTIONS.md.")
    _VERSION_CHECKED.add(FERMI_ENV)


@functools.lru_cache(maxsize=1)
def _pfiles_dir() -> str:
    """A parameter-file directory private to THIS process.

    🚨 WITHOUT THIS, TWO REDUCTIONS RUNNING AT ONCE CORRUPT EACH OTHER.
    The Fermi ST tools do not take their parameters on the command line: they
    WRITE them into a `.par` file and then read them back. `PFILES` names where.
    Its default resolves to `.:$HOME/pfiles;<syspfiles>`, whose writable half
    begins with the CURRENT DIRECTORY -- so two processes launched from the repo
    root both write `./gtselect.par`, and whichever reads second gets the other's
    parameters. Nothing errors; the wrong selection is simply applied.

    Per-PROCESS rather than per-reduction is sufficient because reductions are
    sequential within a process; it is concurrent processes that collide.

    Setting `PFILES` to the writable directory ALONE is deliberate: the env's
    `activate_fermitools.sh` appends `<syspfiles>` when it sees a value without a
    ';', giving exactly `<ours>;<syspfiles>`. Passing the read-only half
    ourselves makes it appear twice.
    """
    d = tempfile.mkdtemp(prefix=f"gce_pfiles_{os.getpid()}_")
    atexit.register(lambda: shutil.rmtree(d, ignore_errors=True))
    return d


def tool_env() -> dict:
    """Environment for an ST-tool subprocess, with an isolated `PFILES`.

    Shared by `_run` here and by `gcereduce.srcmaps`, so that a `gtsrcmaps` run
    and a reduction can proceed concurrently.
    """
    env = os.environ.copy()
    env["PFILES"] = _pfiles_dir()
    return env


def conda_prefix(tool: str) -> list[str]:
    """`conda run` prefix for launching `tool` inside the fermi environment."""
    _assert_fermitools_ok()
    return [_conda(), "run", "-n", FERMI_ENV, tool]


def _run(tool: str, **kwargs) -> None:
    """Invoke a Fermi ST tool, failing loudly.

    The tools return 0 on some failures and write the diagnosis to stdout, so
    output is captured and surfaced on any non-zero exit rather than discarded.

    ⚠ Launched via `conda run -n <FERMI_ENV>` (changed 2026-08-08). Previously
    the tool was invoked bare, which requires the caller to already be inside the
    fermi environment -- true when the one existing reduction was made, and false
    for anything driven from the jax environment. The command itself is
    unchanged; only how the process is started differs.
    """
    cmd = conda_prefix(tool) + [f"{k}={v}" for k, v in kwargs.items()]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=tool_env())
    if proc.returncode != 0:
        raise RuntimeError(
            f"{tool} failed (exit {proc.returncode})\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  stdout: {proc.stdout[-2000:]}\n"
            f"  stderr: {proc.stderr[-2000:]}"
        )


def _readable(path: Path) -> bool:
    """Is this cached product COMPLETE, not merely present?

    🚨 THIS USED TO CHECK `hdul[0].header` AND NOTHING ELSE. The primary HDU is
    written first and is essentially empty, so a file truncated anywhere in the
    data -- an interrupted `gtselect`, a killed job, a full disk -- passed, was
    cached as complete, and the pipeline resumed on top of a partial photon set.
    Nothing downstream would notice: the counts would simply be low.

    The fix touches the LAST element of every HDU, which is where truncation
    shows, with warnings escalated to errors because astropy WARNS rather than
    raises on a short data section. Measured cost: 6 ms on a 527 MB FT1, since
    only the final block is read, not the array.

    ⚠ It does not certify physics, only structural completeness. A file that
    finished writing with the wrong selection still passes -- that is what the
    config hash is for.
    """
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        import warnings

        from astropy.io import fits

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with fits.open(path) as hdul:
                for hdu in hdul:
                    if hdu.data is not None and len(hdu.data):
                        _ = hdu.data[-1]
        return True
    except Exception:
        return False


def write_energy_bindef(cfg: ReductionConfig, out: Path) -> Path:
    """Build the gtbindef energy-bin file.

    Required because the Cholis binning is NOT uniform in log E -- the top three
    bins are ~3x wider. `ebinalg=LOG` cannot express that, and silently
    substituting a uniform grid would shift every bin edge.
    """
    if _readable(out):
        return out
    edges_mev = [e * 1000.0 for e in cfg.energy_edges_gev]
    txt = out.with_suffix(".txt")
    txt.write_text(
        "".join(f"{lo:.6f} {hi:.6f}\n" for lo, hi in zip(edges_mev, edges_mev[1:]))
    )
    # NB: the parameter is `binfile`, not `infile` (verified against gtbindef.par;
    # `infile=` is rejected as an invalid parameter name).
    _run("gtbindef", bintype="E", binfile=str(txt), outfile=str(out), energyunits="MeV")
    return out


def select_and_filter(cfg: ReductionConfig, cone_dir: Path, workdir: Path,
                      weekly_dir: Path | None = None) -> Path:
    """gtselect (weeks -> one FT1) then gtmktime (GTI cuts).

    Reads the 45 deg GC cone when the ROI fits inside it, and the raw
    all-sky weeklies when it does not. `weekly_dir` overrides the latter.
    """
    filtered = workdir / "ft1_filtered.fits"
    if _readable(filtered):
        return filtered

    # 🚨 WHICH ARCHIVE SERVES THIS ROI. The cone is a 45 deg cut about the GC,
    # so it can only serve an ROI whose CORNERS fall inside it. An ROI that
    # reaches further is not an error -- it reduces from the raw all-sky
    # weeklies, which is why those are kept (CONVENTIONS.md). Added 2026-09-02
    # for the Totani halo ROI (|l| <= 60, |b| <= 60), whose corners sit at
    # 75.5 deg. STRICTLY ADDITIVE: every ROI that the cone could already serve
    # takes exactly the path it took before, so no cached GC product changes.
    #
    # ⚠ `corner` is the diagonal of the (l, b) COORDINATE square, not the
    # great-circle distance to the corner. For a 60 deg half-width the two are
    # 84.85 and 75.52 deg. The flat value is always the LARGER, so routing on it
    # is conservative -- it never sends an ROI to the cone that the cone cannot
    # cover. It is not the number to quote in prose; see totani's
    # `corner_angles_deg` for the pair computed side by side.
    corner = (2 * cfg.roi_half_deg**2) ** 0.5
    if corner <= CONE_RADIUS_DEG:
        route, source_dir = "cone45", cone_dir
        weeks = _cone_week_files(cfg, cone_dir)
        if not weeks:
            raise FileNotFoundError(
                f"no cone files for weeks {cfg.week_min}-{cfg.week_max} in "
                f"{cone_dir}. Run scripts/fetch_weekly.py first."
            )
    else:
        route = "weekly_allsky"
        source_dir = _default_weekly_dir() if weekly_dir is None else weekly_dir
        weeks = _weekly_week_files(cfg, source_dir)
        if not weeks:
            raise FileNotFoundError(
                f"ROI corners at {corner:.1f} deg exceed the "
                f"{CONE_RADIUS_DEG:.0f} deg ingest cone, so this reduction "
                f"needs the RAW all-sky weeklies -- and none for weeks "
                f"{cfg.week_min}-{cfg.week_max} were found in {source_dir}. "
                f"Set $FERMI_DATA, pass weekly_dir=, or mount the archive."
            )
        # ⚠ A MISSING WEEK IS A DATA GAP, NOT A ROUNDING ERROR. It silently
        # shortens the exposure, and because gtltcube is built from this same
        # FT1 the two stay consistent with each other while both being wrong.
        # Checked only on this NEW route: the cone path keeps its exact former
        # behaviour so that no frozen product is revalidated by this change.
        want = set(range(cfg.week_min, cfg.week_max + 1)) - ABSENT_UPSTREAM_WEEKS
        have = {int(re.search(r"_w(\d+)_", p.name).group(1)) for p in weeks}
        missing = sorted(want - have)
        if missing:
            head = ", ".join(f"w{w:03d}" for w in missing[:10])
            raise FileNotFoundError(
                f"{len(missing)} of {len(want)} weeks are absent from "
                f"{source_dir}: {head}{' ...' if len(missing) > 10 else ''}. "
                f"Refusing to reduce a span with holes in it -- fetch them, or "
                f"narrow week_min/week_max to a contiguous range. (Weeks known "
                f"absent upstream are already excluded: "
                f"{sorted(ABSENT_UPSTREAM_WEEKS)}.)"
            )

    listfile = workdir / "ft1_list.txt"
    listfile.write_text("".join(f"{p}\n" for p in weeks))
    # Provenance: which archive actually fed this run (rule 5).
    (workdir / "ft1_source.json").write_text(json.dumps(
        {"route": route, "source_dir": str(source_dir), "n_weeks": len(weeks),
         "week_min": cfg.week_min, "week_max": cfg.week_max,
         "roi_corner_deg": round(corner, 3)}, indent=2, sort_keys=True))

    selected = workdir / "ft1_selected.fits"
    if not _readable(selected):
        _run(
            "gtselect",
            infile=f"@{listfile}", outfile=str(selected),
            ra=GC_RA, dec=GC_DEC, rad=corner + 1.0,
            tmin="INDEF", tmax="INDEF",
            emin=cfg.energy_edges_gev[0] * 1000.0,
            emax=cfg.energy_edges_gev[-1] * 1000.0,
            zmax=cfg.zmax_deg, evclass=cfg.evclass, evtype=cfg.evtype,
            chatter=2, clobber="yes",
        )

    _run(
        "gtmktime",
        scfile=str(cfg_scfile()), sctable="SC_DATA",
        filter=cfg.gtmktime_filter, roicut=cfg.roicut,
        evfile=str(selected), evtable="EVENTS",
        outfile=str(filtered), chatter=2, clobber="yes",
    )
    return filtered


def cfg_scfile() -> Path:
    import os
    root = Path(os.environ.get("FERMI_DATA", Path.home() / "data" / "fermi"))
    p = root / "spacecraft" / "lat_spacecraft_merged.fits"
    if not p.exists():
        raise FileNotFoundError(
            f"spacecraft file missing: {p} -- required by gtmktime and gtltcube"
        )
    return p


def bin_counts_cube(cfg: ReductionConfig, ft1: Path, workdir: Path) -> Path:
    """gtbin CCUBE onto the analysis grid."""
    ccube = workdir / "ccube.fits"
    if _readable(ccube):
        return ccube
    ebindef = write_energy_bindef(cfg, workdir / "ebounds.fits")
    _run(
        "gtbin",
        algorithm="CCUBE", evfile=str(ft1), outfile=str(ccube),
        scfile="NONE",
        nxpix=cfg.npix, nypix=cfg.npix, binsz=cfg.pixel_deg,
        coordsys=cfg.coordsys, xref=0.0, yref=0.0, axisrot=0.0, proj=cfg.proj,
        ebinalg="FILE", ebinfile=str(ebindef),
        chatter=2, clobber="yes",
    )
    return ccube


# Event class -> IRF name. The G6 answer (task2_gates.md, Oscar) states the
# selection behind gcepy's templates is P8R3_CLEAN_V3; the others follow the
# same naming. Needed by gtexpcube2, which cannot infer the IRF from the file.
IRF_BY_EVCLASS = {
    128: "P8R3_SOURCE_V3",
    256: "P8R3_CLEAN_V3",
    512: "P8R3_ULTRACLEAN_V3",
    1024: "P8R3_ULTRACLEANVETO_V3",
    2048: "P8R3_SOURCEVETO_V3",
}


def livetime_cube(cfg: ReductionConfig, ft1: Path, ltcube_dir: Path) -> Path:
    """gtltcube -- the slow step, cached under `ltcube_key`, not the config hash.

    Livetime depends only on the time span and the zenith cut, so most of the
    Phase-3 grid shares one cube. Over-keying this (e.g. by the full config
    hash) would recompute hours of work for a config that differs only in
    energy binning, and is the difference between an affordable grid and an
    unaffordable one.
    """
    ltcube_dir.mkdir(parents=True, exist_ok=True)
    out = ltcube_dir / f"ltcube_{cfg.ltcube_key}.fits"
    if _readable(out):
        return out
    _run(
        "gtltcube",
        evfile=str(ft1), evtable="EVENTS",
        scfile=str(cfg_scfile()), sctable="SC_DATA",
        outfile=str(out),
        dcostheta=0.025, binsz=1.0, phibins=0,
        zmax=cfg.zmax_deg,
        chatter=2, clobber="yes",
    )
    return out


def exposure_map(cfg: ReductionConfig, ltcube: Path, workdir: Path) -> Path:
    """gtexpcube2 on the analysis grid.

    ⚠ Deliberately the SAME geometry as the counts cube, because the product
    wanted here is exposure per analysis pixel (for turning model fluxes into
    counts, as gcepy's templates already are). That is NOT the same choice a
    standard binned-likelihood fit makes -- there the exposure map is padded
    well beyond the ROI so PSF leakage from outside is handled. Padding belongs
    to template construction (Task 7), where the convolution happens, and the
    Cholis ingest already showed why it matters: convolving on the full 60 deg
    grid and only then cropping to 40 deg closed a 9% total-count deficit that
    cropping first had produced.
    """
    out = workdir / "expcube.fits"
    if _readable(out):
        return out
    irf = cfg.irf_override or IRF_BY_EVCLASS.get(cfg.evclass)
    if irf is None:
        raise ValueError(f"no IRF mapping recorded for evclass={cfg.evclass}")
    _run(
        "gtexpcube2",
        infile=str(ltcube), cmap="none", outfile=str(out),
        irfs=irf, evtype=cfg.evtype,
        nxpix=cfg.npix, nypix=cfg.npix, binsz=cfg.pixel_deg,
        coordsys=cfg.coordsys, xref=0.0, yref=0.0, axisrot=0.0, proj=cfg.proj,
        ebinalg="FILE", ebinfile=str(workdir / "ebounds.fits"),
        chatter=2, clobber="yes",
    )
    return out


def _fermi_data_root() -> Path:
    import os
    return Path(os.environ.get("FERMI_DATA", Path.home() / "data" / "fermi"))


def _default_cone_dir() -> Path:
    return _fermi_data_root() / "cone45"


def _default_weekly_dir() -> Path:
    """The RAW all-sky weeklies -- the generally-usable archive.

    See CONVENTIONS.md "Photon data: outside the repo, and the RAW files are
    kept". The cone is a derived product for GC work; anything whose ROI reaches
    past 45 deg reduces from here instead.
    """
    return _fermi_data_root() / "weekly"


#: 🚨 WEEKS THAT DO NOT EXIST UPSTREAM, and whose absence is therefore NOT a
#: data gap. w512 returns 404 from the FSSC under every processing version
#: because Fermi was in safe hold (the March 2018 solar-array-drive anomaly):
#: w511 stops at 2018-03-16 04:54 and w513 resumes at 2018-04-02 21:17, leaving
#: ~17.7 days with no science data at all. See CONVENTIONS.md "Week 512 does not
#: exist" and `scripts/fetch_weekly.py`.
#:
#: This was prose in two places and a constant in neither, so the completeness
#: check added 2026-09-02 refused a perfectly good 783-week span on its first
#: run. Recording it here makes the rule executable.
ABSENT_UPSTREAM_WEEKS = frozenset({512})


#: Radius of the derived GC cone product, degrees. Named rather than repeated as
#: a literal, because the number appears in the filename, the ROI check and the
#: error message and they must not drift apart.
CONE_RADIUS_DEG = 45.0


def _cone_week_files(cfg: "ReductionConfig", cone_dir: Path) -> list:
    return sorted(
        p for p in cone_dir.glob("w*_cone45.fits")
        if cfg.week_min <= int(p.name[1:4]) <= cfg.week_max
    )


def _weekly_week_files(cfg: "ReductionConfig", weekly_dir: Path) -> list:
    """Raw FSSC weeklies, `lat_photon_weekly_wNNN_pXXX_vNNN.fits`."""
    out = []
    for p in weekly_dir.glob("lat_photon_weekly_w*.fits"):
        m = re.search(r"_w(\d+)_", p.name)
        if m and cfg.week_min <= int(m.group(1)) <= cfg.week_max:
            out.append((int(m.group(1)), p))
    return [p for _, p in sorted(out)]


def reduce(cfg: ReductionConfig, cone_dir: Path | None = None,
           runs_dir: Path = Path("runs"), with_exposure: bool = True,
           weekly_dir: Path | None = None) -> Path:
    """Run the reduction for one config; returns its run directory."""
    cone_dir = _default_cone_dir() if cone_dir is None else cone_dir
    workdir = runs_dir / cfg.hash
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "config.json").write_text(json.dumps(cfg.to_dict(), indent=2, sort_keys=True))

    ft1 = select_and_filter(cfg, cone_dir, workdir, weekly_dir=weekly_dir)
    bin_counts_cube(cfg, ft1, workdir)
    if with_exposure:
        lt = livetime_cube(cfg, ft1, runs_dir / "_ltcubes")
        exposure_map(cfg, lt, workdir)
    return workdir


def counts_cube_array(workdir: Path, gcepy_orientation: bool = True):
    """Load the CCUBE as a (n_ebins, npix, npix) array.

    🚨 THE LONGITUDE AXIS IS FLIPPED RELATIVE TO GCEPY, and this is the second
    time this project has met that bug. gtbin writes a standards-compliant WCS
    with ``CDELT1 = -0.1``: Galactic longitude DECREASES as the column index
    increases (longitude to the left, the astronomical display convention).
    gcepy's arrays run the other way. So ``[:, :, ::-1]`` is required before any
    pixel-by-pixel comparison against gcepy, and ``gcepy_orientation=True``
    (the default) applies it.

    MEASURED, not inferred (9 weeks, 2026-07-30). Against gcepy's 300 brightest
    OFF-PLANE pixels -- off-plane so the longitude-symmetric Galactic plane
    cannot contribute to the discrimination -- the four candidate orientations
    recover this fraction of the expected counts:

        flip lon 0.775  |  as-is 0.261  |  flip lat 0.247  |  flip both 0.240

    and on 0.8 deg-binned maps the correlations are 0.990 / 0.883 / 0.896 /
    0.870. (The winner falls short of 1.0 only because these 9 early-mission
    weeks have very non-uniform exposure; the discrimination is the point.)

    ⚠ A TOTALS CHECK IS BLIND TO THIS -- a flip does not change a sum. That is
    exactly how the Cholis exposure map passed a clean validation while being
    wrong in every pixel (CONVENTIONS.md). Always follow with a morphology test.
    """
    from astropy.io import fits

    with fits.open(workdir / "ccube.fits") as hdul:
        arr = hdul[0].data.copy()
    return arr[:, :, ::-1] if gcepy_orientation else arr


def exposure_at_edges(workdir: Path, gcepy_orientation: bool = True):
    """Return ``(energies_mev, exposure)`` from gtexpcube2.

    🚨 THE EXPOSURE IS EVALUATED AT THE BIN EDGES, NOT AVERAGED WITHIN BINS.
    For 14 analysis bins the cube has **15** planes, one per edge -- confirmed
    from the ENERGIES extension, whose values are exactly our 15 edges. Treating
    plane *i* as "the exposure of bin *i*" silently uses the bin's lower edge.

    Measured on the smoke-test reduction: exposure at 275 MeV is **0.624x** its
    value at 4910 MeV, so that mistake biases a bin by tens of percent, and
    worst in the top three bins -- which are ~3x wider than the rest and are
    exactly where high-mass DM sensitivity lives.

    NO BIN-AVERAGED REDUCER IS PROVIDED ON PURPOSE. The correct quantity is
    ``N = int flux(E) * exposure(E) dE`` over the bin, which depends on the
    assumed spectrum; collapsing exposure to one number per bin first is only
    valid for a flat spectrum. The integration belongs to template construction
    (Task 7), where the spectrum is known.

    The longitude flip of :func:`counts_cube_array` applies here too: the two
    cubes were verified to carry byte-identical WCS keywords, so whatever
    orientation fixes one fixes the other.
    """
    from astropy.io import fits

    with fits.open(workdir / "expcube.fits") as hdul:
        arr = hdul[0].data.copy()
        energies = hdul["ENERGIES"].data["Energy"].copy()
    if gcepy_orientation:
        arr = arr[:, :, ::-1]
    return energies, arr


if __name__ == "__main__":
    cfg = ReductionConfig(note="Cholis/gcepy anchor selection")
    print(f"config hash   : {cfg.hash}")
    print(f"ltcube key    : {cfg.ltcube_key}")
    print(f"grid          : {cfg.n_ebins} x {cfg.npix} x {cfg.npix}")
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        print(f"run directory : {reduce(cfg)}")
