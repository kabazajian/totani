"""Export the input intensity maps and the exact fitted footprint.

This adapter never builds exposure/count caches or fits data. Regenerating the
NFW intensity is deliberately opt-in; the adopted constructor is relatively
expensive. Paths below describe existing repository assets, not new data.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from astropy.io import fits


RUN_ID = "fff2de5bc7e0"
STRUCTURED_FILE = "FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits"
BRIGHTNESS_FILE = "erosita_brightness_fullsky_normalized.fits"
FLAT_FILE = "data/external/FermiBubblesPublicFiles/bubbleTemp_wcs.fits"


def _roots(config):
    repo = Path(config["repo_root"]).expanduser().resolve()
    data = Path(config.get("data_root") or repo).expanduser().resolve()
    gce = Path(config.get("gce_fisher_root") or os.environ.get(
        "GCE_FISHER_ROOT", repo.parent / "gce-fisher")).expanduser().resolve()
    return repo, data, gce


def _first_existing(paths):
    return next((p for p in paths if p.is_file()), paths[0])


def _paths(config):
    repo, data, gce = _roots(config)
    run = data / "runs" / config.get("run_id", RUN_ID)
    cache = data / "results" / "fits" / "ranking_cells"
    tag = config.get("fine_tag", "psfuc")
    suffix = f"_{tag}" if tag else ""
    eroot = config.get("erosita_root")
    brightness_candidates = []
    if eroot:
        eroot = Path(eroot).expanduser().resolve()
        brightness_candidates += [eroot / BRIGHTNESS_FILE,
                                  eroot / "eROSITA_fullsky_lat_handoff" / BRIGHTNESS_FILE,
                                  eroot / "external" / "eROSITA_fullsky_lat_handoff" / BRIGHTNESS_FILE]
    brightness_candidates += [data / "external" / "eROSITA_fullsky_lat_handoff" / BRIGHTNESS_FILE,
                              repo / "external" / "eROSITA_fullsky_lat_handoff" / BRIGHTNESS_FILE]
    fermi = Path(config.get("fermi_data_root") or os.environ.get(
        "FERMI_DATA", Path.home() / "data" / "fermi")).expanduser().resolve()
    nfw = config.get("nfw_intensity_path")
    return {
        "geometry": run / "ccube.fits",
        "flat": gce / FLAT_FILE,
        "structured": _first_existing([data / "external" / STRUCTURED_FILE,
                                        repo / "external" / STRUCTURED_FILE]),
        "erosita": _first_existing(brightness_candidates),
        "source_catalog": fermi / "catalog" / "gll_psc_v35.fit",
        "coarse_unconvolved": cache / "cells2" / "library.npz",
        "coarse_response": cache / f"cells2{suffix}" / "library.npz",
        "fine_response": cache / f"cells1{suffix}" / "library.npz",
        "nfw": Path(nfw).expanduser().resolve() if nfw else None,
    }


def _available(path):
    """A Git LFS pointer is not an array, even when named .fits/.npz."""
    if path is None or not path.is_file():
        return False
    with path.open("rb") as stream:
        return not stream.read(64).startswith(b"version https://git-lfs.github.com/spec/")


def inspect_spatial(config):
    paths = _paths(config)
    requirements = {
        "geometry": ("degrees; Galactic CAR", "actual reduction count-cube WCS"),
        "flat": ("dimensionless 0/1 intensity shape", "Ackermann 2014 flat Fermi-bubbles map"),
        "structured": ("arbitrary intensity-template units", "Ackermann 2017 Figure 12; Macias inpainting"),
        "erosita": ("sr^-1 template weights", "normalized eROSITA T_N/T_S and uniform U_N/U_S supports"),
        "source_catalog": ("degrees", "4FGL-DR4 ExtendedSources exclusion geometry"),
        "coarse_unconvolved": ("expected template counts", "cells2 unconvolved flat support for common mask"),
        "coarse_response": ("response-folded template counts", "2-degree keep and J for counts-weight fraction"),
        "fine_response": ("response-folded template counts", "1-degree keep and J to verify common coverage"),
        "nfw": ("GeV^2 cm^-5 sr^-1", "input J intensity FITS, or explicit adopted-constructor opt-in"),
    }
    out = []
    for asset, path in paths.items():
        units, requirement = requirements[asset]
        constructor = asset == "nfw" and path is None and config.get("allow_nfw_constructor", False)
        out.append({"asset": asset, "path": str(path) if path else None,
                    "available": bool(constructor or _available(path)), "units": units,
                    "requirement": requirement,
                    **({"producer": "totani.halo.nfw_rho2_on_target (explicit opt-in)"} if constructor else {})})
    return out


def _source(path, config):
    """Portable provenance; inventory may show local paths, exported data do not."""
    if path is None:
        return "totani.halo.nfw_rho2_on_target"
    roots = [("$data_root", _roots(config)[1]), ("$repo_root", _roots(config)[0]),
             ("$gce_fisher_root", _roots(config)[2])]
    for name in ("erosita_root", "fermi_data_root"):
        if config.get(name):
            roots.append((f"${name}", Path(config[name]).expanduser().resolve()))
    for label, root in roots:
        try:
            return label + "/" + path.relative_to(root).as_posix()
        except ValueError:
            pass
    return "supplied-input/" + path.name


def _cells_on_pixels(mask, longitude, latitude, cell_deg):
    """Register ascending-longitude cells on the original descending-WCS pixels."""
    edges = np.arange(-60.0, 60.0 + cell_deg / 2, cell_deg)
    col = np.digitize(longitude, edges) - 1
    row = np.digitize(latitude, edges) - 1
    if mask.shape != (len(edges) - 1, len(edges) - 1):
        raise ValueError("cell mask shape disagrees with the 120-degree analysis square")
    good = (col >= 0) & (col < mask.shape[1]) & (row >= 0) & (row < mask.shape[0])
    out = np.zeros(longitude.shape, bool)
    out[good] = mask[row[good], col[good]]
    return out


def _header_dict(header):
    names = ("NAXIS", "NAXIS1", "NAXIS2", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2",
             "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "CUNIT1", "CUNIT2")
    return {k: header[k] for k in names if k in header}


def _display(array, stretch):
    factor = float(np.max(np.abs(array)))
    return {"rule": "whole_template_max_abs", "factor": factor if factor > 0 else 1.0,
            "stretch": stretch, "independent_half_normalizations": False}


def _common_region(paths, L, B, arrays):
    with np.load(paths["coarse_unconvolved"], allow_pickle=False) as plain, \
            np.load(paths["coarse_response"], allow_pickle=False) as coarse, \
            np.load(paths["fine_response"], allow_pickle=False) as finer:
        pkeep, ckeep, fkeep = [np.asarray(z["keep"], bool) for z in (plain, coarse, finer)]
        flat = np.asarray(plain["T__bubble_flat"], float)
        if pkeep.shape != (60, 60) or ckeep.shape != (60, 60) or fkeep.shape != (120, 120):
            raise ValueError("common region requires the actual 2-degree and 1-degree libraries")
        if flat.shape != (9, 60, 60) or not np.isfinite(flat).all():
            raise ValueError("unconvolved flat-template cache has unexpected shape or missing values")
        if not np.array_equal(pkeep, ckeep):
            raise ValueError("unconvolved and folded 2-degree observation masks differ")
        out2 = ckeep & (np.abs(flat).sum(axis=0) == 0)
        contained = np.kron(out2, np.ones((2, 2), bool))
        out1 = fkeep & contained
        if not np.array_equal(out1, contained):
            raise ValueError("fine observation mask drops pixels from the declared common coarse region")
        pixel_coarse = _cells_on_pixels(out2, L, B, 2)
        pixel_fine = _cells_on_pixels(out1, L, B, 1)
        if not np.array_equal(pixel_coarse, pixel_fine):
            raise ValueError("coarse and contained finer masks do not cover identical sky pixels")
        J2, J1 = np.asarray(coarse["J"], float), np.asarray(finer["J"], float)
        if J2.shape != (9, 60, 60) or J1.shape != (9, 120, 120):
            raise ValueError("response J arrays have unexpected energy/cell shapes")
        for J in (J2, J1):
            if not np.isfinite(J).all() or np.any(J < 0):
                raise ValueError("response-folded NFW arrays must be finite and nonnegative")
        denominator2, denominator1 = float(J2[5][ckeep].sum()), float(J1[5][fkeep].sum())
        if denominator2 <= 0 or denominator1 <= 0:
            raise ValueError("NFW response template has zero fitted-domain weight")
        fraction2 = float(J2[5][out2].sum() / denominator2)
        fraction1 = float(J1[5][out1].sum() / denominator1)
        if not np.isclose(fraction2, fraction1, rtol=1e-6, atol=1e-8):
            raise ValueError("counts-weighted common-region NFW fractions differ between resolutions")
    if "observed" in arrays:
        arrays["outside_common"] = pixel_coarse & arrays["observed"]
    arrays["outside_common_cells2"] = out2
    arrays["outside_common_cells1"] = out1
    return {"coarse_cells": int(out2.sum()), "fine_cells": int(out1.sum()),
            "coarse_cell_deg": 2, "fine_cell_deg": 1,
            "weight_definition": "response_folded_counts", "weight_fraction": fraction2,
            "fine_weight_fraction": fraction1, "weight_energy_bin": 5,
            "weight_note": "sum J[5] in common region / sum J[5] in fitted domain; already includes exposure, PSF, solid angle and source exclusions; not sky area",
            "coverage_identical": True, "fine_to_coarse_cell_ratio": 4,
            "manuscript_coarse_cells": 2190, "manuscript_fine_cells": 8760,
            "matches_manuscript_cell_counts": bool(out2.sum() == 2190 and out1.sum() == 8760),
            "definition": "zero unconvolved flat-template weight in every energy bin on the 2-degree cells, carried to their contained 1-degree cells"}


def export_spatial(config):
    """Return available real inputs and honest incompleteness without any refits."""
    from totani.templates import (target_header, grid_lb, pixel_solid_angle,
                                 bubble_on_target, bubble_macias_on_target)
    from gcereduce.regrid import car_edges

    paths = _paths(config)
    available = {x["asset"]: x["available"] for x in inspect_spatial(config)}
    meta = {"status": "incomplete", "missing": [], "errors": {}, "maps": {},
            "coordinate_system": "Galactic", "projection": "plate-carree",
            "checks": {}, "provenance": config.get("provenance", "scientific"),
            "spatial_provenance": {}, "observation": {},
            "available_maps": []}
    arrays = {}
    target = target_header()
    if available["geometry"]:
        try:
            from totani.stack import _strip_to_2d
            with fits.open(paths["geometry"], memmap=True) as hdus:
                target = _strip_to_2d(hdus[0].header)
                if hdus[0].header.get("TESTDATA"):
                    meta["provenance"] = "test"
            meta["checks"]["actual_reduction_wcs"] = True
            meta["geometry_source"] = _source(paths["geometry"], config)
        except (OSError, KeyError, ValueError):
            available["geometry"] = False
            target = target_header()
            meta["errors"]["geometry"] = "Actual reduction count-cube WCS could not be read"
    if not available["geometry"]:
        meta["missing"].append("geometry")
        meta["checks"]["actual_reduction_wcs"] = False
        meta["geometry_source"] = "totani.templates.target_header (adopted geometry; actual count-cube WCS unavailable)"
    if target.get("CTYPE1") != "GLON-CAR" or target.get("CTYPE2") != "GLAT-CAR":
        raise ValueError("the figure requires the actual Galactic CAR grid")
    L, B = grid_lb(target)
    arrays["longitude_edges"] = car_edges(target, 1, int(target["NAXIS1"]))
    arrays["latitude_edges"] = car_edges(target, 2, int(target["NAXIS2"]))
    arrays["plane_excluded"] = np.abs(B) < 10.0
    meta["wcs_header"] = _header_dict(target)
    meta["orientation"] = {"longitude_left_deg": float(arrays["longitude_edges"][0]),
                           "longitude_right_deg": float(arrays["longitude_edges"][-1]),
                           "latitude_bottom_deg": float(arrays["latitude_edges"][0]),
                           "latitude_top_deg": float(arrays["latitude_edges"][-1]),
                           "derived_from_wcs": True}
    meta["observation"] = {"plane_exclusion": "|b| < 10 degrees",
                           "source_mask_definition": "catalogued extended-source circles at twice Model_SemiMajor; point sources remain fitted", "source_mask_available": False}

    def missing(asset, error=None):
        if asset not in meta["missing"]:
            meta["missing"].append(asset)
        if error is not None:
            # Full error text can contain private machine paths.
            message = str(error)
            for label, path in paths.items():
                if path:
                    message = message.replace(str(path), _source(path, config))
            message = message.replace(str(Path.home()), "$HOME")
            meta["errors"][asset] = type(error).__name__ + ": " + message

    def add_map(name, values, units, normalization, quantity, source, stretch):
        values = np.asarray(values, float)
        if values.shape != L.shape or not np.isfinite(values).all():
            raise ValueError(f"{name} intensity has invalid shape or missing values")
        arrays[name] = values
        meta["available_maps"].append(name)
        meta["maps"][name] = {"units": units, "original_normalization": normalization,
                               "display_normalization": _display(values, stretch),
                               "quantity": quantity, "source": source,
                               "negative_values_preserved": bool(np.any(values < 0))}

    for name in ("flat", "structured"):
        if not available[name]:
            missing(name)
            continue
        try:
            values = (bubble_on_target(target, root=_roots(config)[2]) if name == "flat"
                      else bubble_macias_on_target(target, path=paths[name]))
            with fits.open(paths[name], memmap=True) as hdus:
                if hdus[0].header.get("TESTDATA"):
                    meta["provenance"] = "test"
            add_map(name, values, "dimensionless 0/1 shape" if name == "flat" else "arbitrary intensity-template units",
                    "delivered 0/1 support" if name == "flat" else "unchanged delivered intensity; no BUNIT in source",
                    "input_intensity", _source(paths[name], config), "linear")
            if name == "flat":
                arrays["flat_support"] = values > 0
        except (OSError, KeyError, ValueError) as error:
            missing(name, error)

    if available["nfw"]:
        try:
            if paths["nfw"] is not None:
                with fits.open(paths["nfw"], memmap=True) as hdus:
                    hdr = hdus[0].header
                    if hdr.get("TESTDATA"):
                        meta["provenance"] = "test"
                    if _header_dict(hdr) != _header_dict(target):
                        raise ValueError("supplied NFW intensity WCS disagrees with reduction")
                    values = np.asarray(hdus[0].data, float)
                    if hdr.get("BUNIT", "").replace(" ", "") not in ("GeV^2cm^-5sr^-1", "GeV2cm-5sr-1"):
                        raise ValueError("supplied NFW FITS must declare its intensity units")
                    expected = {"PROFILE": "gnfw", "GAMMA": 1.0, "RS": 21.0,
                                "RHOSUN": 0.42, "RSUN": 8.0, "SMAX": 200.0, "NS": 800}
                    if any(hdr.get(key) != value for key, value in expected.items()):
                        raise ValueError("supplied NFW intensity lacks matching adopted-constructor provenance (PROFILE, GAMMA, RS, RHOSUN, RSUN, SMAX, NS)")
            else:
                from totani.halo import nfw_rho2_on_target
                values = nfw_rho2_on_target(target)
            from totani.halo import TOTANI_NFW_PARAMS
            add_map("nfw", values, "GeV^2 cm^-5 sr^-1", {"constructor": "totani.halo.nfw_rho2_on_target", "parameters": TOTANI_NFW_PARAMS, "s_max_kpc": 200, "n_s": 800},
                    "input_intensity", _source(paths["nfw"], config), "asinh")
        except (OSError, KeyError, ValueError) as error:
            missing("nfw", error)
    else:
        missing("nfw")
        meta["nfw_owner_instruction"] = "Explicitly allow the adopted NFW constructor, or supply its intensity FITS with matching WCS and BUNIT. Count caches cannot substitute for an intensity."

    if available["erosita"]:
        try:
            from totani.erosita import brightness_on_target, EAST_IS_RECONSTRUCTED
            north = brightness_on_target(target, "whole_N", path=paths["erosita"])
            south = brightness_on_target(target, "whole_S", path=paths["erosita"])
            uniform = (brightness_on_target(target, "unif_whole_N", path=paths["erosita"]) +
                       brightness_on_target(target, "unif_whole_S", path=paths["erosita"]))
            with fits.open(paths["erosita"], memmap=True) as hdus:
                if any(hdu.header.get("TESTDATA") for hdu in hdus):
                    meta["provenance"] = "test"
                header_provenance = {name: {key: hdus[name].header.get(key, hdus[0].header.get(key))
                                           for key in ("QUALSTAT", "NATIVE", "BUNIT")}
                                     for name in ("T_N", "T_S")}
            if not EAST_IS_RECONSTRUCTED or any(h["QUALSTAT"] != "RECONSTRUCTED_APPROXIMATION" for h in header_provenance.values()):
                raise ValueError("eROSITA delivery lacks the analysis reconstructed-approximation metadata")
            arrays["erosita_north"] = north
            arrays["erosita_south"] = south
            arrays["erosita_support"] = uniform > 0
            arrays["erosita_native"] = arrays["erosita_support"] & (L < 0)
            arrays["erosita_reconstructed"] = arrays["erosita_support"] & (L >= 0)
            meta["spatial_provenance"].update({"independent_N_S_normalizations": True,
                                      "erosita": {"native": "Galactic West: longitude < 0; German 0.6-1.0 keV survey",
                                                  "reconstructed": "Galactic East: longitude >= 0; published MPE RGB reconstruction calibrated on West",
                                                  "definition_source": "totani.erosita documented delivery and EAST_IS_RECONSTRUCTED; validated FITS QUALSTAT",
                                                  "fits_metadata": header_provenance}})
            add_map("erosita", north + south, "sr^-1 template weights",
                    "T_N and T_S each have unit solid-angle integral over their own all-sky support; ROI values retained without renormalization",
                    "template_weight", _source(paths["erosita"], config), "asinh")
            meta["maps"]["erosita"]["display_label"] = "independently normalized N/S template weights"
        except (OSError, KeyError, ValueError) as error:
            missing("erosita", error)
    else:
        missing("erosita")

    if available["source_catalog"]:
        from totani.pointsources import extended_source_mask
        try:
            with fits.open(paths["source_catalog"], memmap=True) as hdus:
                arrays["source_excluded"] = extended_source_mask(L, B, hdus["ExtendedSources"].data)
            arrays["observed"] = ~arrays["plane_excluded"] & ~arrays["source_excluded"]
            meta["observation"]["source_mask_available"] = True
            meta["observation"]["source"] = _source(paths["source_catalog"], config)
        except (OSError, KeyError, ValueError) as error:
            missing("source_catalog", error)
    else:
        missing("source_catalog")

    common_assets = ("coarse_unconvolved", "coarse_response", "fine_response")
    if all(available[name] for name in common_assets):
        try:
            meta["common_region"] = _common_region(paths, L, B, arrays)
            meta["checks"]["common_coarse_fine_coverage_identical"] = True
            meta["checks"]["common_counts_weight_matches"] = True
        except (OSError, KeyError, ValueError) as error:
            missing("common_region", error)
    else:
        for name in common_assets:
            if not available[name]:
                missing(name)
    meta["checks"]["solid_angle_measure"] = "exact CAR: abs(delta longitude) * abs(delta sin latitude)"
    meta["checks"]["input_maps_are_count_maps"] = False
    meta["status"] = "complete" if not meta["missing"] else "incomplete"
    return meta, arrays
