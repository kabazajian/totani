"""Numeric-only exchange format and checks independent of the LAT environment."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from . import SCHEMA_VERSION


class BundleError(ValueError):
    """A bundle is incomplete, inconsistent, or unsuitable for publication."""


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _file(root, entry):
    relative = Path(entry["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise BundleError("Bundle paths must stay inside the bundle directory")
    path = root / relative
    if not path.is_file() or digest(path) != entry["sha256"]:
        raise BundleError(f"Missing or changed bundle file: {relative}")
    return path


def load_bundle(root, figure="all", allow_test=False):
    root = Path(root)
    try:
        manifest = json.loads((root / "manifest.json").read_text())
    except (OSError, ValueError) as exc:
        raise BundleError(f"Cannot read {root / 'manifest.json'}: {exc}") from exc
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise BundleError("Unsupported figure bundle schema")
    provenance = manifest.get("provenance")
    if provenance == "test" and not allow_test:
        raise BundleError("Production mode rejects TEST DATA")
    if provenance not in ("scientific", "test"):
        raise BundleError("Unknown scientific provenance")
    selected = ("spatial", "regional") if figure == "all" else (figure,)
    loaded = {"manifest": manifest}
    for key in selected:
        entry = manifest.get("figures", {}).get(key)
        if not entry:
            raise BundleError(f"No {key} export in this bundle")
        metadata = json.loads(_file(root, entry["metadata"]).read_text())
        if metadata.get("status") != "complete":
            raise BundleError(f"{key} incomplete: " + "; ".join(metadata.get("missing", [])))
        metadata_provenance = metadata.get("provenance", provenance)
        if metadata_provenance == "test" and not allow_test:
            raise BundleError("Production mode rejects TEST DATA in figure metadata")
        if metadata_provenance not in ("scientific", "test") or metadata_provenance != provenance:
            raise BundleError("Figure provenance must match the bundle provenance")
        if key == "spatial":
            with np.load(_file(root, entry["arrays"]), allow_pickle=False) as z:
                arrays = {name: z[name] for name in z.files}
            validate_spatial(metadata, arrays)
            loaded[key] = (metadata, arrays)
        elif key == "regional":
            validate_regional(metadata)
            loaded[key] = metadata
        else:
            raise BundleError(f"Unknown figure {key}")
    return loaded


def _edges(value, name):
    x = np.asarray(value, float)
    if x.ndim != 1 or len(x) < 2 or not np.isfinite(x).all():
        raise BundleError(f"Invalid {name}")
    delta = np.diff(x)
    if not (np.all(delta > 0) or np.all(delta < 0)):
        raise BundleError(f"Nonmonotonic {name}")
    return x


def validate_spatial(metadata, arrays):
    longitude = _edges(arrays["longitude_edges"], "longitude_edges")
    latitude = _edges(arrays["latitude_edges"], "latitude_edges")
    shape = (len(latitude) - 1, len(longitude) - 1)
    map_names = ("nfw", "flat", "structured", "erosita")
    masks = ("observed", "plane_excluded", "source_excluded", "outside_common",
             "flat_support", "erosita_support", "erosita_native", "erosita_reconstructed")
    for name in map_names + masks:
        if name not in arrays or arrays[name].shape != shape:
            raise BundleError(f"Map/mask registration mismatch: {name}")
        if name in masks and arrays[name].dtype.kind != "b":
            raise BundleError(f"Mask must be Boolean: {name}")
    for name in map_names:
        if np.isinf(arrays[name]).any():
            raise BundleError(f"Infinite intensity: {name}")
        description = metadata["maps"][name]
        if description.get("quantity") not in ("input_intensity", "template_weight"):
            raise BundleError(f"Not an input sky-intensity template: {name}")
        norm = description["display_normalization"]
        factor = norm.get("factor")
        if not np.isscalar(factor) or not np.isfinite(factor) or factor <= 0:
            raise BundleError(f"A whole-template scalar display scale is required: {name}")
        if norm.get("rule") != "whole_template_max_abs":
            raise BundleError(f"Unsupported whole-template normalization: {name}")
        if not np.isclose(factor, np.nanmax(np.abs(arrays[name])), rtol=1e-10):
            raise BundleError(f"Display scale does not match the whole template: {name}")
    if metadata["maps"]["flat"]["display_normalization"]["rule"] != metadata["maps"]["structured"]["display_normalization"]["rule"]:
        raise BundleError("Flat and structured templates use different display rules")
    observed = arrays["observed"]
    if np.any(observed & (arrays["plane_excluded"] | arrays["source_excluded"])):
        raise BundleError("Excluded pixels occur inside the observation mask")
    if np.any(arrays["outside_common"] & ~observed):
        raise BundleError("Common outside-bubbles region exceeds observed sky")
    native, reconstructed = arrays["erosita_native"], arrays["erosita_reconstructed"]
    if np.any(native & reconstructed) or not np.array_equal(arrays["erosita_support"], native | reconstructed):
        raise BundleError("eROSITA provenance masks overlap or fail to cover support")
    if not metadata.get("spatial_provenance", {}).get("independent_N_S_normalizations"):
        raise BundleError("Missing inherited eROSITA N/S template-normalization qualification")
    wcs = metadata.get("wcs_header", {})
    for axis, edges in ((1, longitude), (2, latitude)):
        required = (f"CRPIX{axis}", f"CRVAL{axis}", f"CDELT{axis}")
        if not all(k in wcs for k in required):
            raise BundleError("Missing WCS coordinate-orientation evidence")
        expected = float(wcs[f"CRVAL{axis}"]) + (np.arange(len(edges)) + 0.5 - float(wcs[f"CRPIX{axis}"])) * float(wcs[f"CDELT{axis}"])
        expected = (expected + 180) % 360 - 180 if axis == 1 else expected
        if not np.allclose(edges, expected, rtol=0, atol=1e-7):
            raise BundleError(f"Coordinates do not follow WCS axis {axis}")
    common = metadata.get("common_region", {})
    if common.get("weight_definition") != "response_folded_counts":
        raise BundleError("Outside-region weight must use the reported folded-count definition")
    coarse = arrays.get("outside_common_cells2")
    fine = arrays.get("outside_common_cells1")
    if coarse is None or fine is None:
        raise BundleError("Missing coarse/fine common-region coverage evidence")
    if coarse.dtype.kind != "b" or fine.dtype.kind != "b" or not np.array_equal(np.repeat(np.repeat(coarse, 2, 0), 2, 1), fine):
        raise BundleError("Fine common mask is not the exact contained coarse mask")
    if int(coarse.sum()) != common.get("coarse_cells") or int(fine.sum()) != common.get("fine_cells"):
        raise BundleError("Common-region cell counts disagree with arrays")
    # Cells increase in longitude; the native pixel WCS usually decreases.
    lc = (longitude[:-1] + longitude[1:]) / 2
    bc = (latitude[:-1] + latitude[1:]) / 2
    cell_edges = np.linspace(-60, 60, coarse.shape[1] + 1)
    il = np.clip(np.digitize(lc, cell_edges) - 1, 0, coarse.shape[1] - 1)
    ib = np.clip(np.digitize(bc, cell_edges) - 1, 0, coarse.shape[0] - 1)
    registered = coarse[np.ix_(ib, il)] & observed
    if not np.array_equal(arrays["outside_common"], registered):
        raise BundleError("Common region is not registered to the pixel WCS")
    return {"orientation": "WCS verified", "registered_shape": list(shape),
            "common_cells": [int(coarse.sum()), int(fine.sum())],
            "negative_pixels": {name: int((arrays[name] < 0).sum()) for name in map_names}}


def validate_regional(data):
    edges = _edges(data["energy_edges_gev"], "energy_edges_gev")
    n_bins = len(edges) - 1
    k = data["illustrated_bin"]
    restricted = data["restricted_bins"]
    if n_bins != 9 or k != 5 or restricted != [3, 4, 5, 6, 7, 8]:
        raise BundleError("Wrong energy selection for the approved figure")
    center = np.sqrt(edges[k] * edges[k + 1])
    if abs(center - 20.8) > 0.1:
        raise BundleError("Illustrated energy bin is not 20.8 GeV")
    metadata = data["metadata"]
    if metadata.get("row") != ["SL_Z6_R30_T150_C2", "gas1_ics1"] or metadata.get("cell_deg") != 2 or metadata.get("fine_tag") != "psfuc":
        raise BundleError("Wrong diffuse row, response, or spatial binning")
    if not metadata.get("source_partition_before_response") or not metadata.get("basis_sum_validated"):
        raise BundleError("Missing pre-response regional-basis closure evidence")
    expected = {(config, division) for config in ("structured", "flat", "flat_flexlobes") for division in ("std", "rot")}
    panels = data["panels"]
    if {(p["config"], p["partition"]) for p in panels} != expected or len(panels) != 6:
        raise BundleError("Six distinct approved regional panels are required")
    for panel in panels:
        labels = data["partitions"][panel["partition"]]["labels"]
        if panel["labels"] != labels or len(labels) != 4:
            raise BundleError("Regional label order mismatch")
        a = np.asarray(panel["coefficients"], float)
        lo = np.asarray(panel["profile_lo"], float)
        hi = np.asarray(panel["profile_hi"], float)
        if a.shape != (n_bins, 4) or lo.shape != a.shape or hi.shape != a.shape or not np.isfinite(a).all():
            raise BundleError("Invalid physical regional coefficients/profile shape")
        if not np.isfinite(lo[k]).all() or not np.isfinite(hi[k]).all() or np.any(lo[k] > a[k]) or np.any(hi[k] < a[k]):
            raise BundleError("Missing or invalid 68% profile intervals at 20.8 GeV")
        if panel.get("interval_method") != "profile_delta_neg2loglike_1":
            raise BundleError("Only profiled Δ(−2 log L)=1 intervals are accepted")
        v1 = np.asarray(panel["neg2loglike_common"], float)
        v4 = np.asarray(panel["neg2loglike_regions"], float)
        q = np.asarray(panel["q_per_bin"], float)
        common = np.asarray(panel["common"], float)
        if any(x.shape != (n_bins,) or not np.isfinite(x).all() for x in (v1, v4, q, common)):
            raise BundleError("Invalid likelihood or common-coefficient arrays")
        if not np.allclose(v1 - v4, q, rtol=1e-8, atol=1e-5) or np.any(q < -1e-3):
            raise BundleError("Likelihood-statistic sign or nested-fit mismatch")
        if not np.isclose(panel["q20"], q[k], atol=1e-5) or not np.isclose(panel["q_restricted"], q[restricted].sum(), atol=1e-5):
            raise BundleError("Wrong fixed-partition statistic/energy selection")
        calibration = panel.get("calibration")
        if calibration is not None:
            ne, nt = calibration.get("exceedances"), calibration.get("n_skies")
            if not isinstance(ne, int) or not isinstance(nt, int) or nt <= 0 or not 0 <= ne <= nt:
                raise BundleError("Missing calibration must not become zero exceedances")
            if calibration.get("definition") != "null>=observed" or calibration.get("partition", panel["partition"]) != panel["partition"]:
                raise BundleError("Calibration is not for this fixed partition")
            if calibration.get("null_row") != metadata["row"] and calibration.get("matched_to_data_row") is not False:
                raise BundleError("Reference diffuse-row calibration mislabeled as baseline-matched")
            null = np.asarray(calibration.get("q_null_restricted"), float)
            if null.shape != (nt,) or not np.isfinite(null).all():
                raise BundleError("Missing or nonfinite independent calibration statistics")
            if int(np.count_nonzero(null >= panel["q_restricted"])) != ne:
                raise BundleError("Calibration exceedances disagree with the saved null statistics")
    for config in ("structured", "flat", "flat_flexlobes"):
        left, right = [p for p in panels if p["config"] == config]
        if not np.allclose(left["common"], right["common"], rtol=1e-5, atol=0):
            raise BundleError("Common fits differ between identical partition configurations")
    return {"energy_center_gev": float(center), "panel_count": len(panels),
            "normalization": "physical common basis", "intervals": "profile Δ(−2 log L)=1",
            "calibration": "reference row; fixed partitions"}
