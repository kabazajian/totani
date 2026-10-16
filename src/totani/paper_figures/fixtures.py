"""Small, visibly synthetic layout fixtures, never admissible in production."""
from __future__ import annotations

import numpy as np


def synthetic_figures():
    edges_l = np.linspace(60, -60, 61)
    edges_b = np.linspace(-60, 60, 61)
    L, B = np.meshgrid((edges_l[:-1] + edges_l[1:]) / 2, (edges_b[:-1] + edges_b[1:]) / 2)
    plane = np.abs(B) < 10
    source = (L - 43) ** 2 + (B - 35) ** 2 < 7 ** 2
    observed = ~(plane | source)
    flat_support = (np.abs(L) < 20) & (np.abs(B) < 50)
    erosita_support = (np.abs(L) < 30) & (np.abs(B) < 55)
    nfw = 2 + 20 / (1 + (L / 15) ** 2 + (B / 15) ** 2)
    flat = flat_support.astype(float)
    structured = flat * (0.4 + (B + 60) / 120 + 0.1 * np.cos(L / 10))
    structured[(L < -5) & (B < -30) & flat_support] = -0.1
    erosita = erosita_support * (0.2 + (L + 60) / 120) * np.where(B > 0, 1.0, 0.5)
    common = observed & (np.abs(L) > 32)
    # Likelihood cell arrays increase in longitude; pixels decrease.
    coarse = common[:, ::-1]
    arrays = {"longitude_edges": edges_l, "latitude_edges": edges_b,
              "nfw": nfw, "flat": flat, "structured": structured, "erosita": erosita,
              "observed": observed, "plane_excluded": plane, "source_excluded": source,
              "flat_support": flat_support, "erosita_support": erosita_support,
              "erosita_native": erosita_support & (L < 0),
              "erosita_reconstructed": erosita_support & (L >= 0),
              "outside_common": common, "outside_common_cells2": coarse,
              "outside_common_cells1": np.repeat(np.repeat(coarse, 2, 0), 2, 1)}
    maps = {}
    for name in ("nfw", "flat", "structured", "erosita"):
        maps[name] = {"units": "synthetic units", "original_normalization": "TEST DATA",
                      "source": "explicit synthetic fixture",
                      "quantity": "template_weight" if name == "erosita" else "input_intensity",
                      "display_normalization": {"rule": "whole_template_max_abs", "factor": float(np.abs(arrays[name]).max()),
                                                "stretch": "asinh" if name in ("nfw", "erosita") else "linear"}}
    spatial = {"status": "complete", "missing": [], "provenance": "test", "maps": maps,
               "coordinate_system": "Galactic", "projection": "plate-carree",
               "wcs_header": {"CRVAL1": 0, "CRVAL2": 0, "CRPIX1": 30.5, "CRPIX2": 30.5, "CDELT1": -2, "CDELT2": 2},
               "provenance_note": "TEST DATA - NOT SCIENTIFIC RESULTS",
               "provenance_metadata": {"native": "fixture l<0", "reconstructed": "fixture l>=0"},
               "spatial_provenance": {"independent_N_S_normalizations": True},
               "common_region": {"coarse_cells": int(coarse.sum()), "fine_cells": int(coarse.sum()) * 4,
                                 "weight_definition": "response_folded_counts", "weight_fraction": 0.525}}
    ratio = (814 / 1.51) ** (1 / 12)
    energy_edges = [1.51 * ratio ** (i - 0.5) for i in range(10)]
    partitions = {"std": {"labels": ["WN", "WS", "EN", "ES"]}, "rot": {"labels": ["N", "S", "E", "W"]}}
    masks = {"std": [(L < 0) & (B > 0), (L < 0) & (B < 0), (L >= 0) & (B > 0), (L >= 0) & (B < 0)],
             "rot": [B > np.abs(L), B < -np.abs(L), (L >= 0) & (L >= np.abs(B)), (L < 0) & (L <= -np.abs(B))]}
    schematics = {"longitude_edges": edges_l.tolist(), "latitude_edges": edges_b.tolist(), "observed": observed.tolist()}
    for division, values in masks.items():
        schematics[division] = np.argmax(np.stack(values), axis=0).tolist()
    panels = []
    for i, config in enumerate(("structured", "flat", "flat_flexlobes")):
        common_amp = np.linspace(0.5, 1.4, 9) * 1e-25
        for division in ("std", "rot"):
            offsets = np.array([0.4, 1.4, -0.8, -0.2]) if division == "std" else np.array([0.8, 0.6, -1.4, -0.9])
            if config == "flat_flexlobes" and division == "std":
                offsets *= 0.2
            a = common_amp[:, None] + offsets[None, :] * 1e-25
            errors = np.array([0.15, 0.25, 0.4, 0.2])[None, :] * 1e-25
            q = np.arange(1, 10, dtype=float) * (i + 1) * (2 if division == "rot" else 1)
            v1 = np.linspace(1000, 700, 9)
            panels.append({"config": config, "partition": division, "labels": partitions[division]["labels"],
                           "coefficients": a.tolist(), "profile_lo": (a - errors).tolist(), "profile_hi": (a + errors * 1.4).tolist(),
                           "interval_method": "profile_delta_neg2loglike_1", "common": common_amp.tolist(),
                           "neg2loglike_common": v1.tolist(), "neg2loglike_regions": (v1 - q).tolist(),
                           "q_per_bin": q.tolist(), "q20": float(q[5]), "q_restricted": float(q[3:].sum()),
                           "calibration": {"exceedances": 107 if config == "flat_flexlobes" and division == "std" else 0,
                                           "n_skies": 200, "null_row": ["SO_Z6_R30_T150_C5", "gas1_ics1"],
                                           "matched_to_data_row": False, "definition": "null>=observed", "partition": division,
                                           "q_null_restricted": ([float(q[3:].sum() / 2)] * 93 + [float(q[3:].sum() + 1)] * 107)
                                           if config == "flat_flexlobes" and division == "std" else [float(q[3:].sum() / 2)] * 200}})
    regional = {"status": "complete", "missing": [], "provenance": "test", "metadata": {
        "row": ["SL_Z6_R30_T150_C2", "gas1_ics1"], "cell_deg": 2, "fine_tag": "psfuc",
        "coefficient_units": "synthetic units", "source_partition_before_response": True, "basis_sum_validated": True},
        "partitions": partitions, "schematics": schematics, "energy_edges_gev": energy_edges, "illustrated_bin": 5,
        "restricted_bins": [3, 4, 5, 6, 7, 8], "panels": panels}
    return spatial, arrays, regional
