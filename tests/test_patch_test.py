"""The Zhang+ patch geometry and the reflection null it depends on.

The whole point of `scripts/patch_test.py` is that reflecting a patch in
l -> -l or b -> -b leaves the angular distance from the Galactic Centre
invariant, so each reflected pair is matched in theta_GC by construction
rather than by a cut. If that were not exactly true the test would be
comparing different annuli and every number in it would be meaningless, so
it is asserted here rather than assumed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from patch_test import FAMILIES, VALIDATION_BOX, _box_mask  # noqa: E402
from totani.pointsources import cell_edges  # noqa: E402
from totani.templates import grid_lb, target_header  # noqa: E402


def _grid():
    return grid_lb(target_header())


def _theta(L, B):
    return np.degrees(np.arccos(np.clip(
        np.cos(np.radians(B)) * np.cos(np.radians(L)), -1, 1)))


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_reflections_preserve_theta_gc(family):
    """Every member of a family samples the SAME angular distances."""
    L, B = _grid()
    th = _theta(L, B)
    ref = None
    for name, box in FAMILIES[family].items():
        vals = np.sort(th[_box_mask(L, B, box)])
        if ref is None:
            ref, ref_name = vals, name
            continue
        assert vals.shape == ref.shape, f"{name} has a different pixel count than {ref_name}"
        assert np.allclose(vals, ref, atol=1e-9), \
            f"{name} does not match {ref_name} in theta_GC"


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_patches_are_disjoint_and_inside_the_roi(family):
    L, B = _grid()
    masks = [_box_mask(L, B, b) for b in FAMILIES[family].values()]
    for i, a in enumerate(masks):
        assert a.any(), "a patch selected no pixels"
        assert np.abs(L[a]).max() <= 60.0 and np.abs(B[a]).max() <= 60.0
        for b in masks[i + 1:]:
            assert not (a & b).any(), "patches overlap"


def test_zhang_boxes_are_the_published_ones():
    """Guard against a silent edit: these are Zhang+ 2408.06312 sec 3.1."""
    assert FAMILIES["R1"]["R1"] == (20.0, 40.0, 45.0, 60.0)
    assert FAMILIES["R3"]["R3"] == (20.0, 40.0, -45.0, -30.0)


def test_validation_box_is_an_exact_union_of_cells():
    """The script's cell-stack cross-check is only valid if this box tiles cells."""
    e = cell_edges(60.0, 10.0)
    for v in VALIDATION_BOX:
        assert np.isclose(v, e).any(), f"{v} is not on the 10-degree cell grid"


def test_zhang_latitudes_do_not_tile_cells():
    """Why the test must sum pixels: R1/R3 latitudes cut cells in half.

    If this ever starts passing the other way, the pixel-sum path could be
    replaced by the cell arrays -- so it is worth knowing.
    """
    e = cell_edges(60.0, 10.0)
    l0, l1, b0, b1 = FAMILIES["R1"]["R1"]
    assert np.isclose(l0, e).any() and np.isclose(l1, e).any()
    assert not np.isclose(b0, e).any(), "R1's lower latitude now lands on a cell edge"
