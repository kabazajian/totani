"""Small, explicitly synthetic checks of the spatial export adapter."""
from pathlib import Path
import sys

import numpy as np
import pytest
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.paper_figures.spatial import (  # noqa: E402
    _cells_on_pixels, _common_region, export_spatial, inspect_spatial,
)
from totani.templates import target_header, grid_lb  # noqa: E402


def _config(tmp_path):
    return {"repo_root": tmp_path / "repo", "data_root": tmp_path / "data",
            "gce_fisher_root": tmp_path / "gce", "erosita_root": tmp_path / "erosita",
            "fermi_data_root": tmp_path / "fermi", "fine_tag": "psfuc",
            "provenance": "test"}


def _write_image(path, data, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(np.asarray(data), header=header).writeto(path)


def test_inventory_does_not_mistake_lfs_pointer_for_map(tmp_path):
    config = _config(tmp_path)
    path = config["repo_root"] / "external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits"
    path.parent.mkdir(parents=True)
    path.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 100\n")
    inventory = {x["asset"]: x for x in inspect_spatial(config)}
    assert inventory["structured"]["available"] is False
    assert inventory["erosita"]["path"].startswith(str(config["erosita_root"]))
    assert inventory["flat"]["path"].startswith(str(config["gce_fisher_root"]))
    assert inventory["geometry"]["path"].startswith(str(config["data_root"]))


def test_default_export_never_regenerates_nfw_or_invents_source_mask(tmp_path, monkeypatch):
    from totani import halo
    def forbidden(*args, **kwargs):
        raise AssertionError("default export tried to regenerate NFW")
    monkeypatch.setattr(halo, "nfw_rho2_on_target", forbidden)
    meta, arrays = export_spatial(_config(tmp_path))
    assert meta["status"] == "incomplete"
    assert {"nfw", "source_catalog", "geometry", "coarse_unconvolved"} <= set(meta["missing"])
    assert "observed" not in arrays and "source_excluded" not in arrays
    assert arrays["longitude_edges"][0] == 60 and arrays["longitude_edges"][-1] == -60
    assert arrays["latitude_edges"][0] == -60 and arrays["latitude_edges"][-1] == 60
    assert not meta["checks"]["actual_reduction_wcs"]


def test_registration_respects_descending_pixel_longitude():
    L, B = np.meshgrid([45, 15, -15, -45], [-45, -15, 15, 45])
    mask = np.zeros((4, 4), bool)
    mask[0, 0] = True
    registered = _cells_on_pixels(mask, L, B, 30)
    assert registered[0, 3]
    assert registered.sum() == 1
    with pytest.raises(ValueError, match="shape"):
        _cells_on_pixels(np.zeros((3, 3), bool), L, B, 30)


def _common_inputs(tmp_path):
    outside = np.zeros((60, 60), bool)
    # Include an asymmetric pair so accidental longitude reversal is detectable.
    outside[5, 2] = outside[49, 41] = True
    keep = np.ones((60, 60), bool)
    keep[25:35] = False
    flat = np.ones((9, 60, 60))
    flat[:, outside] = 0
    # Exposure/PSF count weights differ strongly from raw intensity/area.
    J2 = np.full((9, 60, 60), 2.0)
    J2[:, outside] = 19.0
    J1 = np.repeat(np.repeat(J2 / 4, 2, axis=1), 2, axis=2)
    paths = {name: tmp_path / f"{name}.npz" for name in
             ("coarse_unconvolved", "coarse_response", "fine_response")}
    np.savez(paths["coarse_unconvolved"], keep=keep, T__bubble_flat=flat)
    np.savez(paths["coarse_response"], keep=keep, J=J2)
    np.savez(paths["fine_response"], keep=np.repeat(np.repeat(keep, 2, 0), 2, 1), J=J1)
    return paths, outside, keep, J2, J1


def test_common_mask_uses_contained_cells_and_folded_count_weight(tmp_path):
    paths, outside, keep, J2, _ = _common_inputs(tmp_path)
    L, B = grid_lb(target_header(npix=120, pixel_deg=1))
    observed = np.abs(B) >= 10
    observed[10, -5] = False  # a masked source pixel, not a dropped entire cell
    arrays = {"observed": observed}
    meta = _common_region(paths, L, B, arrays)
    assert meta["coarse_cells"] == 2 and meta["fine_cells"] == 8
    assert np.array_equal(arrays["outside_common_cells2"], outside)
    assert np.array_equal(arrays["outside_common_cells1"], np.repeat(np.repeat(outside, 2, 0), 2, 1))
    assert np.array_equal(arrays["outside_common"], _cells_on_pixels(outside, L, B, 2) & observed)
    expected = J2[5][outside].sum() / J2[5][keep].sum()
    assert meta["weight_fraction"] == pytest.approx(expected)
    assert meta["weight_fraction"] != pytest.approx(outside.sum() / keep.sum())
    assert meta["weight_definition"] == "response_folded_counts"


def test_common_mask_refuses_finer_cells_outside_declared_common_coverage(tmp_path):
    paths, outside, keep, _, J1 = _common_inputs(tmp_path)
    fine_keep = np.repeat(np.repeat(keep, 2, 0), 2, 1)
    fine_keep[10, 4] = False
    np.savez(paths["fine_response"], keep=fine_keep, J=J1)
    L, B = grid_lb(target_header(npix=120, pixel_deg=1))
    with pytest.raises(ValueError, match="drops pixels"):
        _common_region(paths, L, B, {})


def test_common_mask_refuses_mixed_response_cache(tmp_path):
    paths, _, keep, _, J1 = _common_inputs(tmp_path)
    J1[:, 10:12, 4:6] *= 2
    np.savez(paths["fine_response"], keep=np.repeat(np.repeat(keep, 2, 0), 2, 1), J=J1)
    L, B = grid_lb(target_header(npix=120, pixel_deg=1))
    with pytest.raises(ValueError, match="fractions differ"):
        _common_region(paths, L, B, {})


def test_actual_loaders_preserve_independent_lobe_weights_and_provenance(tmp_path, monkeypatch):
    from totani import halo
    config = _config(tmp_path)
    header = target_header(npix=4, pixel_deg=30)
    L, B = grid_lb(header)
    _write_image(config["data_root"] / "runs/fff2de5bc7e0/ccube.fits", np.ones((4, 4)), header)
    flat = (B > 0).astype(float)
    _write_image(config["gce_fisher_root"] / "data/external/FermiBubblesPublicFiles/bubbleTemp_wcs.fits", flat, header)
    structured = np.arange(16).reshape(4, 4).astype(float)
    _write_image(config["repo_root"] / "external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits", structured, header)
    north = 3.0 * (B > 0)
    south = 11.0 * (B < 0)
    brightness_header = header.copy()
    brightness_header["QUALSTAT"] = "RECONSTRUCTED_APPROXIMATION"
    brightness_header["NATIVE"] = False
    brightness_header["BUNIT"] = "sr^-1"
    hdus = [fits.PrimaryHDU()]
    for name, values in [("T_N", north), ("T_S", south),
                         ("U_N", (B > 0).astype(float)), ("U_S", (B < 0).astype(float))]:
        hdus.append(fits.ImageHDU(values, header=brightness_header, name=name))
    brightness_path = config["erosita_root"] / "erosita_brightness_fullsky_normalized.fits"
    brightness_path.parent.mkdir(parents=True)
    fits.HDUList(hdus).writeto(brightness_path)
    monkeypatch.setattr(halo, "nfw_rho2_on_target", lambda target: np.full((4, 4), 7.0))
    config["allow_nfw_constructor"] = True
    meta, arrays = export_spatial(config)
    assert set(meta["available_maps"]) == {"flat", "structured", "erosita", "nfw"}
    assert np.array_equal(arrays["erosita"], north + south)
    assert meta["maps"]["erosita"]["display_normalization"]["factor"] == 11
    assert meta["maps"]["flat"]["display_normalization"]["rule"] == meta["maps"]["structured"]["display_normalization"]["rule"]
    assert meta["spatial_provenance"]["independent_N_S_normalizations"]
    assert np.array_equal(arrays["erosita_native"], L < 0)
    assert np.array_equal(arrays["erosita_reconstructed"], L >= 0)
    assert meta["provenance"] == "test"
    assert str(tmp_path) not in str(meta)


def test_nfw_count_cache_cannot_be_presented_as_input_intensity(tmp_path):
    config = _config(tmp_path)
    path = tmp_path / "count_map.fits"
    _write_image(path, np.ones((960, 960)), target_header())
    config["nfw_intensity_path"] = path
    meta, arrays = export_spatial(config)
    assert "nfw" in meta["missing"] and "nfw" not in arrays
    assert "intensity units" in meta["errors"]["nfw"]


def test_complete_synthetic_adapter_export_passes_bundle_geometry_checks(tmp_path, monkeypatch):
    from totani import halo
    from totani.paper_figures.bundle import validate_spatial
    config = _config(tmp_path)
    header = target_header(npix=120, pixel_deg=1)
    header["TESTDATA"] = True
    L, B = grid_lb(header)
    _write_image(config["data_root"] / "runs/fff2de5bc7e0/ccube.fits", np.ones(L.shape), header)
    flat = ((np.abs(L) < 20) & (np.abs(B) < 35)).astype(float)
    _write_image(config["gce_fisher_root"] / "data/external/FermiBubblesPublicFiles/bubbleTemp_wcs.fits", flat, header)
    _write_image(config["repo_root"] / "external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits", np.arange(L.size).reshape(L.shape).astype(float), header)
    eh = header.copy()
    eh["QUALSTAT"] = "RECONSTRUCTED_APPROXIMATION"
    eh["NATIVE"] = False
    hdus = [fits.PrimaryHDU()]
    for name, values in [("T_N", 3.0 * (B > 0)), ("T_S", 11.0 * (B < 0)),
                         ("U_N", (B > 0).astype(float)), ("U_S", (B < 0).astype(float))]:
        hdus.append(fits.ImageHDU(values, header=eh, name=name))
    ep = config["erosita_root"] / "erosita_brightness_fullsky_normalized.fits"
    ep.parent.mkdir(parents=True)
    fits.HDUList(hdus).writeto(ep)
    catalog = config["fermi_data_root"] / "catalog/gll_psc_v35.fit"
    catalog.parent.mkdir(parents=True)
    columns = [fits.Column(name=name, format="D", array=[value]) for name, value in
               [("GLON", -55.0), ("GLAT", -49.0), ("Model_SemiMajor", 1.0)]]
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="ExtendedSources")]).writeto(catalog)
    paths, _, _, _, _ = _common_inputs(tmp_path)
    cache = config["data_root"] / "results/fits/ranking_cells"
    for key, directory in [("coarse_unconvolved", "cells2"), ("coarse_response", "cells2_psfuc"),
                           ("fine_response", "cells1_psfuc")]:
        dest = cache / directory / "library.npz"
        dest.parent.mkdir(parents=True)
        dest.write_bytes(paths[key].read_bytes())
    monkeypatch.setattr(halo, "nfw_rho2_on_target", lambda target: np.full(L.shape, 7.0))
    config["allow_nfw_constructor"] = True
    meta, arrays = export_spatial(config)
    assert meta["status"] == "complete", meta
    assert meta["provenance"] == "test"
    checks = validate_spatial(meta, arrays)
    assert checks["orientation"] == "WCS verified"
    assert checks["common_cells"] == [2, 8]
    assert np.any(arrays["source_excluded"])
    assert not np.any(arrays["observed"] & arrays["source_excluded"])
    assert str(tmp_path) not in str(meta)
