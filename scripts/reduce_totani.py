#!/usr/bin/env python
"""Reduce the Fermi-LAT photons for the Totani halo ROI.

    python scripts/reduce_totani.py --plan          # print, touch nothing
    python scripts/reduce_totani.py --run           # via jobq, never bare

Produces a counts cube and an exposure map on a 960x960 Cartesian Galactic grid
at 0.125 deg over |l| <= 60, |b| <= 60, for weeks 009-791 (2008-08-04 to
2023-08-02) of the UltraClean class at zenith < 100 deg.

🚨 THIS ROI CANNOT USE THE 45 DEG GC CONE. Its corners sit at 75.5 deg. The
reduction reads the RAW all-sky weeklies, which is why they are kept; see
`gcereduce.reduce.select_and_filter`, which routes on the corner distance.

⚠ ALL 13 ENERGY BINS ARE REDUCED, though only bins 0-8 can be FITTED. The
limit is the ring inverse-Compton template (top plane 165.4 GeV), not the data.
Reducing the full range keeps the data product general and puts the restriction
where it belongs -- in the fit -- rather than baking a template limitation into
a photon product. See CONVENTIONS.md "Energy binning".
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import totani  # noqa: F401,E402  -- puts the vendored gcereduce on sys.path

from gcereduce.reduce import (  # noqa: E402
    ABSENT_UPSTREAM_WEEKS,
    CONE_RADIUS_DEG,
    IRF_BY_EVCLASS,
    ReductionConfig,
    _weekly_week_files,
    reduce,
)

ROOT = Path(__file__).resolve().parents[1]

#: 🚨 STATED AS BIN *CENTERS*, NOT EDGES (2507.07209 sec 2.1): "13 energy bins,
#: equally spaced in logarithmic space, with the minimum and maximum energy bin
#: centers at 1.51 and 814 GeV". The ratio follows from those two numbers and
#: the count; the EDGES are then the geometric means of adjacent centers.
#:
#: ⚠ The edge convention is an INFERENCE from "equally spaced in logarithmic
#: space" -- he does not print the edges. It is recorded as UNCONFIRMED in
#: CONVENTIONS.md and must be pinned against the PDF before any flux is quoted.
N_BINS = 13
CENTER_LO_GEV, CENTER_HI_GEV = 1.51, 814.0
_RATIO = (CENTER_HI_GEV / CENTER_LO_GEV) ** (1.0 / (N_BINS - 1))

TOTANI_CENTERS_GEV = tuple(CENTER_LO_GEV * _RATIO**k for k in range(N_BINS))
TOTANI_EDGES_GEV = tuple(
    [c / _RATIO**0.5 for c in TOTANI_CENTERS_GEV]
    + [TOTANI_CENTERS_GEV[-1] * _RATIO**0.5]
)

#: Top energy plane of the ring inverse-Compton cubes, GeV (measured from
#: `Templates_arxiv2003.10416/ics_Ring_1.fits`: 27 planes, 500 MeV x 1.25^k,
#: last = 165436 MeV). Bins whose UPPER EDGE exceeds this have no ICS template.
ICS_CEILING_GEV = 165.436

#: Weeks spanning 2008-08-04 -> 2023-08-02. MET 712,627,205 for the end date
#: falls inside w791 (TSTART 712,112,355, TSTOP 712,715,727), verified
#: 2026-09-02 against the file headers.
WEEK_MIN, WEEK_MAX = 9, 791

#: ⚠ The automounter spelling `/Volumes/auto/Data2` and the local NFS mount
#: `/private/nfs/Data2` are the same share, and which one resolves depends on
#: what is mounted right now. `$FERMI_DATA/weekly` is a symlink to the former,
#: so it DANGLES when only the latter is up -- which reads as "no files found".
#: Pass the path explicitly rather than trusting the default.
DEFAULT_WEEKLY_DIR = Path(os.environ.get(
    "TOTANI_WEEKLY_DIR", "/private/nfs/Data2/gce-fisher/fermi/weekly"))


def usable_bins() -> int:
    """How many bins the ring ICS can actually serve, from the low edge up."""
    return sum(1 for hi in TOTANI_EDGES_GEV[1:] if hi <= ICS_CEILING_GEV)


def totani_config(**overrides) -> ReductionConfig:
    """Totani's selection (sec 2.1) on our reduction machinery.

    `evtype=3` (FRONT+BACK) is UNCONFIRMED -- he states no FRONT/BACK split, and
    both is the natural reading, but it is not written down. Do not quote a
    count from this without pinning it.
    """
    kw = dict(
        week_min=WEEK_MIN, week_max=WEEK_MAX,
        energy_edges_gev=TOTANI_EDGES_GEV,
        roi_half_deg=60.0, pixel_deg=0.125,
        evclass=512,          # UltraClean
        evtype=3,             # FRONT+BACK -- UNCONFIRMED
        zmax_deg=100.0,
        note="Totani 2507.07209 halo ROI; all 13 bins reduced, bins 0-8 fittable",
    )
    kw.update(overrides)
    return ReductionConfig(**kw)


def corner_angles_deg(half_deg: float) -> tuple[float, float]:
    """(flat diagonal, true great-circle distance) from the GC to an ROI corner.

    🚨 THESE ARE DIFFERENT NUMBERS AND ONLY ONE OF THEM IS AN ANGLE ON THE SKY.
    `gcereduce` routes on `sqrt(2) * half`, the diagonal of the (l, b) COORDINATE
    square -- 84.85 deg here. The actual angular separation from (0,0) to
    (60, 60) is `acos(cos b cos l)` = 75.52 deg. The flat value is the larger of
    the two, so using it for the cone test is conservative and safe; but it is
    not the number to quote in prose, and CONVENTIONS.md quotes 75.5.

    Reported as a pair so the two can never be confused again.
    """
    import math
    flat = (2 * half_deg**2) ** 0.5
    r = math.radians(half_deg)
    angular = math.degrees(math.acos(math.cos(r) * math.cos(r)))
    return flat, angular


def describe(cfg: ReductionConfig, weekly_dir: Path) -> dict:
    corner, corner_angular = corner_angles_deg(cfg.roi_half_deg)
    files = _weekly_week_files(cfg, weekly_dir)
    want = set(range(cfg.week_min, cfg.week_max + 1)) - ABSENT_UPSTREAM_WEEKS
    # `_w` also matches inside "weekly" -- use the same regex the upstream
    # helper uses rather than a fragile split.
    have = {int(re.search(r"_w(\d+)_", p.name).group(1)) for p in files}
    return {
        "run_dir": f"runs/{cfg.hash}",
        "config_hash": cfg.hash,
        "ltcube_key": cfg.ltcube_key,
        "grid": f"{cfg.npix}x{cfg.npix} at {cfg.pixel_deg} deg",
        "roi_corner_flat_deg": round(corner, 2),      # what the routing uses
        "roi_corner_angular_deg": round(corner_angular, 2),  # the sky separation
        "archive_route": "weekly_allsky" if corner > CONE_RADIUS_DEG else "cone45",
        "weekly_dir": str(weekly_dir),
        "weeks_found": len(files),
        "weeks_expected": len(want),
        "weeks_missing": sorted(want - have),
        "irf": cfg.irf_override or IRF_BY_EVCLASS[cfg.evclass],
        "n_ebins_reduced": cfg.n_ebins,
        "n_ebins_fittable": usable_bins(),
        "fit_range_gev": [round(TOTANI_EDGES_GEV[0], 4),
                          round(TOTANI_EDGES_GEV[usable_bins()], 3)],
        "ics_ceiling_gev": ICS_CEILING_GEV,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--plan", action="store_true",
                    help="print the configuration and archive state; run nothing")
    ap.add_argument("--run", action="store_true", help="run the reduction")
    ap.add_argument("--weekly-dir", type=Path, default=DEFAULT_WEEKLY_DIR)
    ap.add_argument("--runs-dir", type=Path, default=ROOT / "runs")
    ap.add_argument("--no-exposure", action="store_true",
                    help="skip gtltcube/gtexpcube2 (the multi-hour step)")
    args = ap.parse_args()

    cfg = totani_config()
    info = describe(cfg, args.weekly_dir)
    print(json.dumps(info, indent=2))

    for k, c in enumerate(TOTANI_CENTERS_GEV):
        lo, hi = TOTANI_EDGES_GEV[k], TOTANI_EDGES_GEV[k + 1]
        mark = "fit " if hi <= ICS_CEILING_GEV else "NO-ICS"
        print(f"  bin {k:2d}  {lo:9.3f} - {hi:9.3f} GeV   center {c:8.2f}   {mark}")

    if info["weeks_missing"]:
        print(f"\nBLOCKED: {len(info['weeks_missing'])} weeks absent from "
              f"{args.weekly_dir}: {info['weeks_missing'][:10]}", file=sys.stderr)
        return 1

    if not args.run:
        print("\n--plan only; nothing was run.")
        return 0

    work = reduce(cfg, runs_dir=args.runs_dir,
                  with_exposure=not args.no_exposure,
                  weekly_dir=args.weekly_dir)
    print(f"\nreduction complete: {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
