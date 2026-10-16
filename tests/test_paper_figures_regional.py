"""Focused checks for the normalization and original quadrant exporter API."""
from pathlib import Path

import numpy as np
import pytest

from totani.paper_figures.regional import (
    export_regional, inspect_regional, partition_masks, recover_profile_intervals,
    solver_to_common, validate_regional_basis,
)


def test_schematic_orientation_comes_from_repository_wcs():
    from totani.paper_figures.regional import _schematics

    schematic = _schematics()
    assert schematic["longitude_edges"][0] == 60
    assert schematic["longitude_edges"][-1] == -60
    assert schematic["latitude_edges"][0] == -60
    assert schematic["latitude_edges"][-1] == 60
    assert schematic["wcs_header"]["CDELT1"] == -1
    # Left is positive longitude E, and top is positive latitude N.
    assert schematic["std"][-1][0] == 2  # EN
    assert schematic["std"][-1][-1] == 0  # WN
    assert schematic["rot"][0][0] == 2  # E tie
    assert schematic["rot"][0][-1] == 3  # W tie
    assert schematic["source_exclusions_shown"] is False


@pytest.mark.parametrize("bins", [[1.] * 8, [1.] * 10, None,
                                  [np.nan] + [1.] * 8,
                                  [1.] * 3 + [np.inf] + [1.] * 5,
                                  [1e308] * 9])
def test_nonfinite_or_wrong_shape_null_cannot_become_zero_exceedances(bins):
    from totani.paper_figures.regional import _calibration

    null = {"kind": "null", "config": "flat", "sim": 0, "std": bins}
    with pytest.raises(ValueError, match="null.*statistic|calibration statistic"):
        _calibration([null], 100., "flat", "std")


@pytest.mark.parametrize("identifier", [None, "0", True, -1])
def test_calibration_requires_independent_integer_sky_identifiers(identifier):
    from totani.paper_figures.regional import _calibration

    null = {"kind": "null", "config": "flat", "sim": identifier, "std": [1.] * 9}
    with pytest.raises(ValueError, match="sky identifiers"):
        _calibration([null], 100., "flat", "std")


def test_calibration_rejects_duplicate_sky_ids_and_nonfinite_observation():
    from totani.paper_figures.regional import _calibration

    null = {"kind": "null", "config": "flat", "sim": 0, "std": [1.] * 9}
    with pytest.raises(ValueError, match="duplicate"):
        _calibration([null, null.copy()], 100., "flat", "std")
    with pytest.raises(ValueError, match="observed"):
        _calibration([null], np.nan, "flat", "std")
    summary = _calibration([null], 6., "flat", "std")
    assert summary["exceedances"] == 1 and summary["definition"] == "null>=observed"


def test_exact_source_partition_orientation_and_diagonal_ties():
    L, B = np.meshgrid([-1.5, -.5, .5, 1.5], [-1.5, -.5, .5, 1.5])
    masks = partition_masks(L, B)
    assert list(masks["std"]) == ["WN", "WS", "EN", "ES"]
    assert list(masks["rot"]) == ["N", "S", "E", "W"]
    assert masks["std"]["WN"][3, 0] and masks["std"]["ES"][0, 3]
    assert masks["rot"]["W"][0, 0] and masks["rot"]["E"][3, 3]
    for pieces in masks.values():
        assert np.all(sum(m.astype(int) for m in pieces.values()) == 1)
    with pytest.raises(ValueError, match="source grid"):
        partition_masks(np.array([0.5]), np.array([0.0]))


def test_unequal_template_weights_with_response_spillover():
    # Each source region contributes to multiple observed cells.  A mistaken
    # conversion based on observed-quadrant photons would fail this fixture.
    response = np.array([[.8, .1, 0, .1], [.2, .7, .1, 0],
                         [0, .2, .8, .1], [0, 0, .1, .8]])
    weights = np.array([2., 5., 11., 17.])
    H = (response * weights).T.reshape(4, 1, 2, 2)
    J = H.sum(0)
    assert validate_regional_basis(J, H)["passed"]
    physical = np.array([-.2, .3, .5, .8])
    solver_counts = physical * H.sum((1, 2, 3))
    np.testing.assert_allclose(solver_to_common(solver_counts, H.sum((1, 2, 3))), physical)
    observed_photons = np.tensordot(physical, H[:, 0], axes=1).ravel()
    assert not np.allclose(observed_photons / weights, physical)
    interval_lo = (physical - .1) * weights
    np.testing.assert_allclose(solver_to_common(interval_lo, weights), physical - .1)


def test_basis_check_includes_observation_cells_outside_the_fitted_cut():
    H = np.ones((4, 1, 2, 2))
    J = H.sum(0)
    H[0, 0, 0, 0] += 1
    with pytest.raises(ValueError, match="sum to the full"):
        validate_regional_basis(J, H)
    with pytest.raises(ValueError, match="positive"):
        solver_to_common([1, 2], [1, 0])


def test_intervals_profile_all_other_amplitudes_in_physical_units():
    scales = np.array([2., 5., 11., 17.])
    signals = [np.array([[scale]]) for scale in scales]
    best = np.array([-.4, .3, .5, .8])
    count_sigma = np.array([.3, .4, .6, .9])
    calls = []

    def original_fixed_api(T, n, sigs, fixed=None):
        # Analytic joint nuisance optimum is represented by this fixture;
        # the adapter must hold only the single signal under examination.
        assert T.shape == (2, 1, 1) and sigs is signals
        assert len(fixed) == 1
        calls.append(dict(fixed))
        i, coefficient = next(iter(fixed.items()))
        delta = ((coefficient - best[i]) * scales[i] / count_sigma[i]) ** 2
        return -10 + delta, best.copy()

    lo, hi, scans = recover_profile_intervals(
        original_fixed_api, np.ones((2, 1, 1)), np.ones((1, 1)), signals, -10., best)
    np.testing.assert_allclose(lo, best - count_sigma / scales, atol=1e-5)
    np.testing.assert_allclose(hi, best + count_sigma / scales, atol=1e-5)
    assert {next(iter(fixed)) for fixed in calls} == {0, 1, 2, 3}
    assert all("other three regional amplitudes" in s["profiled_parameters"] for s in scans)
    assert all(min(s["delta_neg2loglike"]) == 0 for s in scans)
    assert lo[0] < best[0] < hi[0] < 0


def test_actual_repository_joint_solver_with_unequal_weights_and_spillover(monkeypatch):
    """One tiny artificial Poisson bin; this never reads or fits real photons."""
    import sys
    from scipy.ndimage import gaussian_filter

    repo = Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(repo / "scripts"))
    from ns_halo_spectra import joint

    centres = np.array([-2.5, -1.5, -.5, .5, 1.5, 2.5])
    L, B = np.meshgrid(centres, centres)
    J = np.exp(-(L ** 2 + B ** 2) / 8)
    pieces = partition_masks(L, B)["std"]
    weights = [2., 5., 11., 17.]
    signals = [gaussian_filter(J * piece, .6, mode="nearest") * weight
               for piece, weight in zip(pieces.values(), weights)]
    T = np.ones((1, 6, 6))
    true_coefficients = np.array([-2., 3., 2., 4.])
    n = 50. * T[0] + np.tensordot(true_coefficients, np.stack(signals), axes=1)
    assert np.min(n) > 0
    # Gaussian spillover creates nonzero response outside each source region.
    assert all(np.any(signal[~mask] > 0) for signal, mask in zip(signals, pieces.values()))
    value, amplitudes = joint(T, n, signals)
    np.testing.assert_allclose(amplitudes, true_coefficients, rtol=2e-5, atol=2e-5)
    lo, hi, scans = recover_profile_intervals(joint, T, n, signals, value, amplitudes)
    for i in range(4):
        assert lo[i] < amplitudes[i] < hi[i]
        for endpoint in (lo[i], hi[i]):
            constrained_value, _ = joint(T, n, signals, fixed={i: endpoint})
            assert constrained_value - value == pytest.approx(1., abs=5e-4)
    assert len(scans) == 4


def test_default_real_saved_export_is_incomplete_and_does_not_refit(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[1]
    monkeypatch.setattr("totani.paper_figures.regional._recover",
                        lambda *args: pytest.fail("default export must not start fits"))
    config = {"repo_root": repo, "data_root": tmp_path / "missing-owner-data",
              "bundle_dir": tmp_path / "bundle"}
    result = export_regional(config)
    assert result["status"] == "incomplete" and len(result["missing"]) == 6
    assert not config["bundle_dir"].exists() and not config["data_root"].exists()
    expected = {("structured", "std"): 155.9, ("structured", "rot"): 223.1,
                ("flat", "std"): 161.8, ("flat", "rot"): 208.8,
                ("flat_flexlobes", "std"): 15.7, ("flat_flexlobes", "rot"): 174.6}
    for panel in result["panels"]:
        assert round(panel["q_restricted"], 1) == expected[panel["config"], panel["partition"]]
        assert panel["calibration"]["matched_to_data_row"] is False
        assert panel["calibration"]["n_skies"] == 200
        assert panel["coefficients"][5] == [None] * 4
        np.testing.assert_allclose(np.array(panel["neg2loglike_common"]) - np.array(panel["neg2loglike_regions"]), panel["q_per_bin"], atol=1e-8)
    flexstd = next(p for p in result["panels"] if (p["config"], p["partition"]) == ("flat_flexlobes", "std"))
    assert flexstd["calibration"]["exceedances"] == 107
    assert not result["schematics"]["source_exclusions_shown"]
    assert any(not item["available"] for item in inspect_regional(config))


@pytest.mark.parametrize("fine_tag,cell,row", [("psf", 2, None), ("psfuc", 1, None),
                                               ("psfuc", 2, ["SO_Z6_R30_T150_C5", "gas1_ics1"])])
def test_wrong_primary_analysis_is_rejected(tmp_path, fine_tag, cell, row):
    config = {"repo_root": Path(__file__).resolve().parents[1], "bundle_dir": tmp_path,
              "fine_tag": fine_tag, "cell_deg": cell}
    if row is not None:
        config["row"] = row
    with pytest.raises(ValueError, match="primary regional figure"):
        export_regional(config)
