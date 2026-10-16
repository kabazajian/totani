"""Put the VENDORED gce-fisher code on `sys.path`, and locate gce-fisher's DATA.

**Decision, Kev 2026-09-18: this project no longer imports code from the
`../gce-fisher` checkout.** From 2026-09-02 to 2026-09-18 it did (`gcereduce`
and the A2020 helper scripts were imported from the sibling repo so the
reduction code had one implementation). That coupling was accepted while
gce-fisher was frozen. It stopped being frozen on 2026-09-16 -- 58 modified
files in its working tree while this project is trying to hold its results
steady for a write-up -- and the price of the coupling became "an edit over
there can change our numbers, and nothing here notices until a test runs".

So the code is now VENDORED, byte for byte, in `totani/_vendor/`:

  * `_vendor/src/gcereduce/`  -- `reduce`, `regrid`, `dm_templates`, `forward`,
    `galp21`, taken from gce-fisher commit `62f05c7` (2026-09-16); none of
    them was dirty in that checkout when copied.
  * `_vendor/scripts/`        -- `a2020_templates_official`,
    `a2020_stage2_diffuse`, `a2020_stage2_reduce`: the delivered-file naming,
    the CONFIGS table, the H2 header correction and the log10-E integration
    of the IC cube. Same commit.

The layout mirrors upstream (`src/` and `scripts/` under one root) so the
files import each other unchanged, `diff -r` against a gce-fisher checkout is
meaningful, and `sha256` of every file is pinned in `_vendor/MANIFEST.sha256`
and asserted by `tests/test_imports.py`. Any local edit to a vendored file is
therefore deliberate: edit, re-pin, and say why in `_vendor/UPSTREAM.md`.

⚠ WHAT IS STILL SHARED IS DATA, NOT CODE. The A2020 delivery, the Ackermann
bubble template and the Galp21 maps live under the gce-fisher checkout and are
read from there by path. `GCE_FISHER_ROOT` is that data root, overridable with
`$GCE_FISHER_ROOT`. Data files are not being edited; if they ever are, the
same fix applies (copy and pin). The Python environment is also still
gce-fisher's `.venv` -- an environment coupling, not a code one.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

__all__ = ["GCE_FISHER_ROOT", "VENDOR_ROOT", "gce_fisher_data", "ensure_on_path",
           "ensure_scripts_on_path"]

#: The vendored copy of the gce-fisher code this project uses.
VENDOR_ROOT = Path(__file__).resolve().parent / "_vendor"

#: Sibling of this repo by default: ~/work/gce-fisher next to ~/work/totani.
#: DATA ONLY -- nothing is imported from here any more.
_DEFAULT = Path(__file__).resolve().parents[3] / "gce-fisher"

GCE_FISHER_ROOT = Path(os.environ.get("GCE_FISHER_ROOT", _DEFAULT)).expanduser()


def gce_fisher_data(subpath: str = "data/external") -> Path:
    """A directory under the gce-fisher DATA root, verified to exist.

    Raises with the path it looked at rather than failing later on a missing
    FITS file that does not say where it was expected.
    """
    p = GCE_FISHER_ROOT / subpath
    if not p.is_dir():
        raise FileNotFoundError(
            f"gce-fisher data not found at {p}. This project reads the A2020 "
            f"delivery and the bubble template from there (code is vendored, "
            f"data is not). Set $GCE_FISHER_ROOT to the checkout, or clone it "
            f"beside this repo.")
    return p


def _vendor_dir(name: str, probe: str) -> Path:
    d = VENDOR_ROOT / name
    if not (d / probe).is_file():
        raise FileNotFoundError(
            f"vendored gce-fisher code missing: {d / probe}. The copy lives in "
            f"this repository; a missing file means the checkout is incomplete.")
    return d


def ensure_on_path() -> Path:
    """Idempotently prepend the vendored `src/` (holds `gcereduce`) to `sys.path`."""
    src = _vendor_dir("src", "gcereduce/reduce.py")
    s = str(src)
    if s not in sys.path:
        sys.path.insert(0, s)
    return src


def ensure_scripts_on_path() -> Path:
    """Also expose the vendored A2020 helper scripts for import.

    ⚠ Needed because some knowledge lives in SCRIPTS rather than in the package:
    `a2020_templates_official` holds the delivered-file naming, the CONFIGS
    table and the H2 header correction, and `a2020_stage2_diffuse` holds the
    log10-E integration of the IC cube onto analysis bins. Re-deriving any of
    that here would be re-deriving exactly the things that took a delivery and
    two rounds with Oscar to pin down.
    """
    scripts = _vendor_dir("scripts", "a2020_templates_official.py")
    s = str(scripts)
    if s not in sys.path:
        sys.path.insert(0, s)
    return scripts
