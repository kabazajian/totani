"""The eROSITA boundary loader and the matched-theta null it depends on.

`scripts/boundary_test.py` compares two regions of genuinely different shape.
The only thing that makes that a null rather than an apples-to-oranges
comparison is the theta_GC matching: reweighted annulus by annulus, any
spherically symmetric intensity must contribute equally to both sides. If that
were not true every number the script produces would be meaningless, so it is
asserted here rather than assumed -- and so is the converse, that a lobe
bounded by the eROSITA edge is NOT cancelled.

Fast: builds the target grid from `target_header()` and samples one all-sky
map. No run directory needed.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.erosita import (  # noqa: E402
    LOBE_MEAN_LON_DEG,
    SPHERICAL_PROFILES,
    boundary_path,
    check_handedness,
    erosita_fraction_on_target,
    matched_theta_weights,
    spherical_contrast,
    theta_gc,
    unmatched_contrast,
)
from totani.templates import grid_lb, pixel_solid_angle, target_header  # noqa: E402

MAP = ROOT / "external/erosita_validation/eROSITA_boundary_Macias.fits"
SUPPORT_CSV = ROOT / "external/erosita_validation/matched_theta_geometric_support.csv"
pytestmark = pytest.mark.skipif(not MAP.exists(), reason="eROSITA boundary not present")


@pytest.fixture(scope="module")
def sky():
    t = target_header()
    L, B = grid_lb(t)
    om = pixel_solid_angle(t)
    # A stand-in for `flat_sr`: exposure x solid angle. The real exposure needs
    # a run directory, and these tests deliberately do not. What matters is
    # that the response is NOT spherically symmetric -- a matching that only
    # cancels under a symmetric response is the bug this file exists to catch
    # -- so the gradient in l and b is the point, not its realism.
    response = om * (1.0 + 0.30 * np.cos(np.radians(B)) + 0.15 * np.sin(np.radians(L)))
    return t, L, B, om, erosita_fraction_on_target(t), theta_gc(L, B), response


@pytest.fixture(scope="module")
def matched(sky):
    _, L, B, om, inside, _, response = sky
    out = {}
    for lobe, sign in (("north", +1), ("south", -1)):
        keep = sign * B >= 10.0
        out[lobe] = matched_theta_weights(L, B, om, inside, keep, response)
    return out


def test_fraction_is_a_coverage_map(sky):
    _, L, _, _, inside, _, _ = sky
    assert inside.shape == L.shape
    assert np.isfinite(inside).all()
    assert inside.min() >= 0.0 and inside.max() <= 1.0
    assert ((inside > 0) & (inside < 1)).any(), "no partially covered pixels: " \
        "the FRACTION extension is not being read, or the grids have shifted"


def test_grid_phase_is_exact_so_sampling_moves_no_flux():
    """The delivered map and the Totani grid share a pixel phase exactly.

    If a future delivery changes CDELT or CRPIX, `erosita_fraction_on_target`
    must raise rather than silently nearest-neighbour a partial-coverage map.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    t = target_header()
    L, B = grid_lb(t)
    with fits.open(boundary_path()) as h:
        x, y = WCS(h[0].header).wcs_world2pix(L.ravel(), B.ravel(), 0)
    assert np.abs(x - np.rint(x)).max() == 0.0
    assert np.abs(y - np.rint(y)).max() == 0.0


def test_handedness_matches_the_publication(sky):
    _, L, B, om, inside, _, _ = sky
    got = check_handedness(L, B, inside, om)
    assert got["north"] < 0 < got["south"], \
        "the lobes must lean opposite ways -- north to negative l"
    for k, v in got.items():
        assert abs(v - LOBE_MEAN_LON_DEG[k]) < 0.5


def test_handedness_rejects_a_longitude_flipped_map(sky):
    """The trap `CONVENTIONS.md` opens with, made to fail loudly."""
    _, L, B, om, inside, _, _ = sky
    with pytest.raises(ValueError, match="handedness"):
        check_handedness(L, B, inside[:, ::-1], om)


@pytest.mark.parametrize("lobe", ["north", "south"])
def test_matched_weights_cancel_every_spherical_profile(matched, sky, lobe):
    """THE NULL. Five unlike radial profiles, all cancelling to 2.2e-5 or better.

    NOT to machine precision, and the difference matters enough to state: the
    weights are normalised per 0.25 deg annulus but each profile is evaluated
    at every pixel's OWN theta, so what survives is the profile's curvature
    across one annulus. Oscar's `matched_theta_geometric_support.csv` reports
    1e-17 because it evaluates profiles at the annulus centres, where the
    cancellation is exact by construction; that number is a property of his
    quadrature, not a tighter version of this one.

    A constant profile still cancels to 6e-16 -- that part is pure area
    matching and is exact. What is left is curvature: worst case 1.6e-5 north
    and 2.1e-5 south, both on the 25 deg Gaussian.

    2e-5 is negligible where it is used. The halo contributes of order 1e3
    counts inside, so a spherical component can leak at most ~0.05 counts into
    a statistic whose error bar is tens to hundreds. `boundary_test.py` records
    that leak in counts per bin rather than leaving it as a ratio.
    """
    _, _, _, _, _, th, response = sky
    w_in, w_out, _ = matched[lobe]
    got = spherical_contrast(w_in, w_out, th, response)
    assert set(got) == set(SPHERICAL_PROFILES)
    for name, v in got.items():
        assert abs(v) < 5e-5, f"{name} does not cancel: contrast {v:.3e}"


def test_matching_on_solid_angle_instead_of_response_does_not_cancel(sky):
    """THE BUG THIS API SHAPE EXISTS TO PREVENT, pinned so it cannot return.

    The first version of `matched_theta_weights` normalised annuli by solid
    angle. Under a non-spherical response that leaves a 14% contrast on a
    CONSTANT profile -- four orders of magnitude worse than the 1e-5 the
    response-weighted version achieves -- because the statistic sums counts,
    and counts carry the exposure. On the real exposure it showed up as the
    fitted NFW halo predicting 34% more counts outside than inside at 20.8 GeV.
    """
    _, L, B, om, inside, th, response = sky
    keep = B >= 10.0
    w_in, w_out, _ = matched_theta_weights(L, B, om, inside, keep, om)
    got = spherical_contrast(w_in, w_out, th, response)
    assert abs(got["constant"]) > 0.05


def test_quarter_degree_annuli_are_the_optimum(sky):
    """Why the default is 0.25 deg and not simply "as fine as possible".

    Coarser annuli lose to the profile's curvature across one annulus; finer
    ones lose because each annulus holds too few pixels for the inside and
    outside theta distributions to match within it. The residual is therefore
    NOT monotone in the step, and 0.25 deg sits at the minimum -- measured
    2026-09-08: 6.4e-4, 1.3e-4, 1.6e-5, 2.2e-5, 3.9e-5 at 1, 0.5, 0.25, 0.125
    and 0.0625 deg. It is also the step Oscar's own quadrature uses.
    """
    _, L, B, om, inside, th, response = sky
    keep = B >= 10.0

    def worst(step):
        w_in, w_out, _ = matched_theta_weights(L, B, om, inside, keep, response,
                                               theta_step_deg=step)
        return max(abs(v) for v in spherical_contrast(w_in, w_out, th, response).values())

    assert worst(0.25) < worst(1.0)
    assert worst(0.25) < worst(0.0625)


@pytest.mark.parametrize("lobe", ["north", "south"])
def test_the_matching_is_what_does_it(sky, lobe):
    """Without matching the same profiles differ by 100%, not by 2e-5.

    Otherwise the null above would pass for a trivial reason -- two regions
    that happened to span the same annuli anyway -- and prove nothing. Taken as
    they are, inside and outside sit at very different mean theta_GC, and a
    falling profile is a factor two brighter on one side than the other.
    """
    _, L, B, om, inside, th, response = sky
    keep = (B >= 10.0) if lobe == "north" else (B <= -10.0)
    got = unmatched_contrast(inside, keep, om, th, response)
    assert max(abs(v) for v in got.values()) > 0.5


@pytest.mark.parametrize("lobe", ["north", "south"])
def test_a_lobe_bounded_by_the_edge_survives_the_matching(matched, sky, lobe):
    """THE POSITIVE CONTROL. A construction that cancels everything is useless.

    Unit intensity inside the eROSITA footprint and zero outside is exactly the
    signal this test exists to find, and the matched difference must recover
    essentially all of it.
    """
    _, _, _, _, inside, _, response = sky
    w_in, w_out, _ = matched[lobe]
    lobe_counts = inside * response
    num = float((w_in * lobe_counts).sum() - (w_out * lobe_counts).sum())
    den = float((w_in * lobe_counts).sum() + (w_out * lobe_counts).sum())
    assert 2 * num / den > 1.9, "the geometric lobe was cancelled away too"


@pytest.mark.parametrize("lobe", ["north", "south"])
def test_both_sides_carry_the_same_total_response(matched, lobe):
    _, _, info = matched[lobe]
    assert info["response_matched_in"] == pytest.approx(
        info["response_matched_out"], rel=1e-12)
    assert info["response_matched_in"] > 0


@pytest.mark.parametrize("lobe", ["north", "south"])
def test_no_support_inside_34_degrees(matched, lobe):
    """The restriction that limits what this test can say, asserted so it
    cannot quietly change: inside ~34.5 deg there is no outside region."""
    _, _, info = matched[lobe]
    assert 33.0 < info["theta_min_deg"] < 36.0
    assert info["inside_area_dropped_sr"] > 0.5 * info["inside_area_sr"] * 0.1


@pytest.mark.skipif(not SUPPORT_CSV.exists(), reason="support CSV not present")
@pytest.mark.parametrize("lobe", ["north", "south"])
def test_support_agrees_with_the_delivered_quadrature(matched, lobe):
    """Our pixel-grid support must agree with Oscar's independent 0.25 deg
    spherical quadrature. Two different constructions, one geometry."""
    rows = [r for r in csv.DictReader(open(SUPPORT_CSV))
            if r["lobe"] == lobe and r["both_sides_supported"] == "1"]
    theirs = min(float(r["theta_center_deg"]) for r in rows)
    _, _, info = matched[lobe]
    assert abs(info["theta_min_deg"] - theirs) < 1.0


def test_matched_weights_reject_a_region_with_no_outside(sky):
    """Restricting to theta < 30 leaves nothing to compare against."""
    _, L, B, om, inside, th, response = sky
    keep = (B >= 10.0) & (th < 30.0)
    with pytest.raises(ValueError, match="no theta annulus"):
        matched_theta_weights(L, B, om, inside, keep, response)


def test_node_weight_choice_moves_significance_but_not_the_null(sky):
    """'harmonic' is the variance-optimal annulus weight, 'min' is Oscar's.
    Both must satisfy the null; they differ only in how they weight annuli."""
    _, L, B, om, inside, th, response = sky
    keep = B >= 10.0
    a = matched_theta_weights(L, B, om, inside, keep, response, node_weight="min")
    b = matched_theta_weights(L, B, om, inside, keep, response, node_weight="harmonic")
    for w_in, w_out, _ in (a, b):
        assert max(abs(v) for v in spherical_contrast(w_in, w_out, th, response).values()) < 5e-5
    assert not np.allclose(a[0], b[0]), "the two weightings should not coincide"


def test_source_hash_guards_against_a_different_raster():
    """A re-delivery built from another figure or another paper version must
    not load silently -- the reconstruction was verified against page 28 of
    our own copy of arXiv:2408.06312v2 byte-for-byte."""
    from astropy.io import fits

    from totani.erosita import SOURCE_SHA256

    with fits.open(boundary_path()) as h:
        assert h[0].header["SRCSHA"] == SOURCE_SHA256
        assert h[0].header["SRCFIG"].startswith("Extended Data Figure 5")


# ---- the generalised weights behind the spatial breakdown (2026-09-18) ------

def test_regions_form_reproduces_the_special_case_exactly(sky):
    """`matched_theta_weights` IS `matched_theta_weights_regions` with
    f_in = inside*keep, f_out = (1-inside)*keep; refactored 2026-09-18 and
    pinned so the 09-08 result cannot drift through the refactor."""
    from totani.erosita import matched_theta_weights_regions

    _, L, B, om, inside, _, response = sky
    keep = B >= 10.0
    a, b, info = matched_theta_weights(L, B, om, inside, keep, response)
    a2, b2, info2 = matched_theta_weights_regions(
        L, B, om, inside * keep, (1.0 - inside) * keep, response)
    assert np.array_equal(a, a2) and np.array_equal(b, b2)
    assert info == info2


@pytest.mark.parametrize("which", ["l_pos", "b_45_60", "far_from_gc"])
def test_a_sub_region_matched_against_the_full_outside_still_cancels(sky, which):
    """Every sub-region in `scripts/boundary_breakdown.py` is matched against
    the SAME full outside. The null must hold for each piece, not only for the
    whole: a sub-region has a narrower theta range and fewer supported annuli,
    so this is a genuinely different check from the whole-inside one."""
    from totani.erosita import matched_theta_weights_regions

    _, L, B, om, inside, th, response = sky
    keep = B >= 10.0
    S = {"l_pos": L > 0,
         "b_45_60": (B >= 45) & (B < 60),
         "far_from_gc": th > 50}[which]
    a, b, info = matched_theta_weights_regions(
        L, B, om, inside * keep * S, (1.0 - inside) * keep, response)
    assert info["supported_nodes"] > 20
    contrast = spherical_contrast(a, b, th, response)
    assert max(abs(v) for v in contrast.values()) < 1e-4, contrast
    naive = unmatched_contrast(inside * S, keep, om, th, response)
    assert max(abs(v) for v in naive.values()) > 0.05, "the matching must be doing something"


def test_a_lobe_confined_to_a_sub_region_shows_up_there_and_not_elsewhere(sky):
    """Positive control for the breakdown: emission planted in the l > 0 half
    of the northern inside is seen by the l > 0 sub-region's weights and NOT
    by the l < 0 sub-region's, both matched to the same outside."""
    from totani.erosita import matched_theta_weights_regions

    _, L, B, om, inside, th, response = sky
    keep = B >= 10.0
    f_out = (1.0 - inside) * keep
    # Planted in FULLY-inside pixels only. A pixel straddling the boundary
    # carries weight on both sides by design, so a lobe drawn with the
    # fractional map would legitimately show a percent-level share on the
    # outside; that is the partial-pixel treatment, not a leak. Fully inside,
    # the outside weight is identically zero and so must the leak be.
    planted = (inside == 1.0) * keep * (L > 0) * response
    seen = {}
    for tag, S in (("l_pos", L > 0), ("l_neg", L < 0)):
        a, b, _ = matched_theta_weights_regions(L, B, om, inside * keep * S, f_out, response)
        seen[tag] = float((a * planted).sum() - (b * planted).sum())
    assert seen["l_pos"] > 0.5 * float(planted.sum())
    assert seen["l_neg"] == 0.0
