"""Scientific invariants that the portable renderer must enforce."""
from copy import deepcopy
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from totani.paper_figures.bundle import (BundleError, digest, load_bundle,
                                        validate_regional, validate_spatial, write_json)
from totani.paper_figures.fixtures import synthetic_figures


def test_orientation_and_common_coverage_are_checked():
    spatial, arrays, _ = synthetic_figures()
    validate_spatial(spatial, arrays)
    bad = dict(arrays, longitude_edges=arrays["longitude_edges"][::-1])
    with pytest.raises(BundleError, match="WCS axis 1"):
        validate_spatial(spatial, bad)
    bad = dict(arrays, outside_common=arrays["outside_common"][:, ::-1])
    with pytest.raises(BundleError, match="registered|observed sky"):
        validate_spatial(spatial, bad)


def test_finer_region_cannot_be_reselected_independently():
    spatial, arrays, _ = synthetic_figures()
    finer = arrays["outside_common_cells1"].copy()
    finer[0, 0] = ~finer[0, 0]
    with pytest.raises(BundleError, match="contained coarse mask"):
        validate_spatial(spatial, dict(arrays, outside_common_cells1=finer))


def test_signed_input_values_are_preserved():
    spatial, arrays, _ = synthetic_figures()
    checked = validate_spatial(spatial, arrays)
    assert checked["negative_pixels"]["structured"] > 0
    assert arrays["structured"].min() < 0


def test_whole_template_scale_and_provenance_are_required():
    spatial, arrays, _ = synthetic_figures()
    bad = deepcopy(spatial)
    bad["maps"]["erosita"]["display_normalization"]["factor"] *= 0.5
    with pytest.raises(BundleError, match="whole template"):
        validate_spatial(bad, arrays)
    with pytest.raises(BundleError, match="provenance masks"):
        validate_spatial(spatial, dict(arrays, erosita_native=arrays["erosita_reconstructed"]))


def test_profile_intervals_and_common_fit_are_not_replaced():
    _, _, data = synthetic_figures()
    validate_regional(data)
    bad = deepcopy(data)
    bad["panels"][0]["profile_lo"][5] = [None] * 4
    with pytest.raises(BundleError, match="68% profile"):
        validate_regional(bad)
    bad = deepcopy(data)
    bad["panels"][0]["interval_method"] = "hessian"
    with pytest.raises(BundleError, match="profiled"):
        validate_regional(bad)
    bad = deepcopy(data)
    bad["panels"][1]["common"][5] *= 2
    with pytest.raises(BundleError, match="Common fits differ"):
        validate_regional(bad)


def test_q_plus_or_max_partition_cannot_replace_q_part():
    _, _, data = synthetic_figures()
    bad = deepcopy(data)
    bad["panels"][0]["q_restricted"] = bad["panels"][0]["q20"]
    with pytest.raises(BundleError, match="statistic/energy selection"):
        validate_regional(bad)
    bad = deepcopy(data)
    bad["panels"][0]["calibration"]["partition"] = "max"
    with pytest.raises(BundleError, match="fixed partition"):
        validate_regional(bad)


def test_missing_calibration_is_not_zero_and_transfer_is_explicit():
    _, _, data = synthetic_figures()
    bad = deepcopy(data)
    bad["panels"][0]["calibration"]["n_skies"] = None
    with pytest.raises(BundleError, match="Missing calibration"):
        validate_regional(bad)
    bad = deepcopy(data)
    bad["panels"][0]["calibration"]["matched_to_data_row"] = True
    with pytest.raises(BundleError, match="mislabeled"):
        validate_regional(bad)
    bad = deepcopy(data)
    bad["panels"][0]["calibration"]["q_null_restricted"][0] = None
    with pytest.raises(BundleError, match="nonfinite independent"):
        validate_regional(bad)
    bad = deepcopy(data)
    bad["panels"][0]["calibration"]["exceedances"] = 1
    with pytest.raises(BundleError, match="exceedances disagree"):
        validate_regional(bad)


def test_test_data_rejected_even_if_manifest_provenance_is_changed(tmp_path):
    _, _, data = synthetic_figures()
    write_json(tmp_path / "regional.json", data)
    manifest = {"schema_version": 1, "provenance": "test", "figures": {
        "regional": {"metadata": {"path": "regional.json", "sha256": digest(tmp_path / "regional.json")}}}}
    write_json(tmp_path / "manifest.json", manifest)
    with pytest.raises(BundleError, match="TEST DATA"):
        load_bundle(tmp_path, "regional")
    load_bundle(tmp_path, "regional", allow_test=True)
    data["provenance"] = {"source": "unrecognized provenance"}
    write_json(tmp_path / "regional.json", data)
    manifest["figures"]["regional"]["metadata"]["sha256"] = digest(tmp_path / "regional.json")
    write_json(tmp_path / "manifest.json", manifest)
    with pytest.raises(BundleError, match="must match"):
        load_bundle(tmp_path, "regional", allow_test=True)
    data["provenance"] = "test"
    write_json(tmp_path / "regional.json", data)
    manifest["figures"]["regional"]["metadata"]["sha256"] = digest(tmp_path / "regional.json")
    manifest["provenance"] = "scientific"
    write_json(tmp_path / "manifest.json", manifest)
    with pytest.raises(BundleError, match="TEST DATA"):
        load_bundle(tmp_path, "regional")


def test_bundle_integrity_and_path_traversal_are_checked(tmp_path):
    _, _, data = synthetic_figures()
    write_json(tmp_path / "regional.json", data)
    manifest = {"schema_version": 1, "provenance": "test", "figures": {
        "regional": {"metadata": {"path": "regional.json", "sha256": digest(tmp_path / "regional.json")}}}}
    write_json(tmp_path / "manifest.json", manifest)
    (tmp_path / "regional.json").write_text(json.dumps({"modified": True}))
    with pytest.raises(BundleError, match="changed bundle file"):
        load_bundle(tmp_path, "regional", allow_test=True)
    manifest["figures"]["regional"]["metadata"]["path"] = "../regional.json"
    write_json(tmp_path / "manifest.json", manifest)
    with pytest.raises(BundleError, match="inside the bundle"):
        load_bundle(tmp_path, "regional", allow_test=True)


def test_direct_render_calls_also_reject_unwatermarked_test_data(tmp_path):
    from totani.paper_figures.rendering import render_spatial, render_regional
    spatial, arrays, regional = synthetic_figures()
    with pytest.raises(ValueError, match="test|TEST"):
        render_spatial(spatial, arrays, tmp_path)
    with pytest.raises(ValueError, match="test|TEST"):
        render_regional(regional, tmp_path)
    assert not list(tmp_path.iterdir())
