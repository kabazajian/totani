"""The vendored gce-fisher code is what gets imported, and it is what was pinned.

From 2026-09-02 to 2026-09-18 this project imported `gcereduce` from the
sibling `../gce-fisher` checkout, and these tests existed to turn "a refactor
over there breaks this repo silently" into "at test time". The code is now
vendored and pinned (`totani/_vendor/`, `totani._paths`), so the failure modes
these guard are different: the sibling checkout being imported by accident
(it is still on this machine and still has a `gcereduce`), and a vendored file
being edited without the manifest being re-pinned.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def test_gcereduce_is_the_vendored_copy():
    import totani
    import gcereduce.reduce
    import gcereduce.regrid
    import gcereduce.dm_templates

    for mod in (gcereduce.reduce, gcereduce.regrid, gcereduce.dm_templates):
        p = Path(mod.__file__).resolve()
        assert totani.VENDOR_ROOT in p.parents, (
            f"{mod.__name__} imported from {p}, not from {totani.VENDOR_ROOT}")
        assert "gce-fisher" not in p.parts


def test_a2020_helpers_are_the_vendored_copies():
    import totani  # noqa: F401
    import totani.templates  # noqa: F401  (calls ensure_scripts_on_path)
    import a2020_templates_official
    import a2020_stage2_diffuse

    import totani as t
    for mod in (a2020_templates_official, a2020_stage2_diffuse):
        p = Path(mod.__file__).resolve()
        assert t.VENDOR_ROOT in p.parents, f"{mod.__name__} imported from {p}"


def test_vendored_files_match_manifest():
    """Any edit to a vendored file must be deliberate: re-pin and say why in
    `_vendor/UPSTREAM.md`."""
    import totani

    manifest = totani.VENDOR_ROOT / "MANIFEST.sha256"
    pinned = {}
    for line in manifest.read_text().splitlines():
        digest, _, rel = line.partition("  ")
        pinned[rel.strip()] = digest.strip()
    assert len(pinned) == 9, sorted(pinned)

    on_disk = {f"./{p.relative_to(totani.VENDOR_ROOT)}": p
               for p in totani.VENDOR_ROOT.rglob("*.py")}
    assert set(on_disk) == set(pinned), (
        f"vendored files and manifest disagree: "
        f"{set(on_disk) ^ set(pinned)}")
    bad = {rel: p for rel, p in on_disk.items()
           if hashlib.sha256(p.read_bytes()).hexdigest() != pinned[rel]}
    assert not bad, f"vendored files edited without re-pinning: {sorted(bad)}"


def test_missing_data_root_names_the_path(monkeypatch, tmp_path):
    """The error must say where it looked, or the next person re-derives it."""
    import totani._paths as paths

    monkeypatch.setattr(paths, "GCE_FISHER_ROOT", tmp_path / "nowhere")
    with pytest.raises(FileNotFoundError, match="nowhere"):
        paths.gce_fisher_data()


def test_reduction_config_takes_the_totani_selection():
    """Every knob this analysis needs already exists in the vendored code."""
    import totani  # noqa: F401
    from gcereduce.reduce import ReductionConfig

    cfg = ReductionConfig(
        week_min=9, week_max=791, roi_half_deg=60.0, pixel_deg=0.125,
        evclass=512, evtype=3, zmax_deg=100.0,
    )
    # 120 deg / 0.125 deg = 960, exactly -- an inexact division raises.
    assert cfg.npix == 960
    # The livetime cube is keyed on span + zenith + filter, so it is shared
    # across every energy/geometry variant of this analysis.
    assert cfg.ltcube_key == ReductionConfig(
        week_min=9, week_max=791, roi_half_deg=20.0, pixel_deg=0.1,
        evclass=512, evtype=3, zmax_deg=100.0,
    ).ltcube_key


def test_j_map_rejects_misspelled_profile_parameters():
    """Totani's halo parameters are NOT the A2020 defaults, so a typo that
    silently returned the default profile would be a wrong number that looks
    right. The vendored code refuses; assert it, because we rely on it."""
    import totani  # noqa: F401
    import numpy as np
    from gcereduce.dm_templates import j_map

    grid_l, grid_b = np.array([[0.0]]), np.array([[20.0]])
    with pytest.raises(TypeError, match="unknown profile parameter"):
        j_map(grid_l, grid_b, r_sun=8.0)          # the field is R_sun
