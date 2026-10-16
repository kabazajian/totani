"""Export the existing response-folded quadrant test, without a default refit.

The published calibration exporter discards the four fitted amplitudes.  Its
saved ``common_amp`` is already in the inherited NFW convention; the counts in
``erosita_flex_ns`` use a different (post-response) partition and cannot fill
that gap.  Recovery below uses the original joint solver, including its fixed
signal API, and never generates a calibration ensemble.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

CONFIGS = ("structured", "flat", "flat_flexlobes")
LABELS = {"std": ("WN", "WS", "EN", "ES"), "rot": ("N", "S", "E", "W")}
RESP_BINS = (3, 4, 5, 6, 7, 8)
BASELINE = ("SL_Z6_R30_T150_C2", "gas1_ics1")
NULL_ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
K20 = 5
DEFAULT_CALIBRATION = "results/fits/quadrant_calibration_cells2_psfuc_2026-09-29.json"


def solver_to_common(coefficients, normalizers):
    """For a solver basis H/c, return a=alpha/c, including signed endpoints.

    This is *not* a conversion of photons observed within one quadrant: PSF
    leakage makes that mapping nondiagonal.  ``ns_halo_spectra.joint`` already
    performs this division, so its outputs must not be divided a second time.
    """
    c = np.asarray(normalizers, float)
    if not np.all(np.isfinite(c)) or np.any(c <= 0):
        raise ValueError("regional solver normalizers must be finite and positive")
    return np.asarray(coefficients, float) / c


def source_partition_masks(longitude, latitude):
    """Exact masks from quadrant_calibration.build_quadrants, before response.

    Pixel centres never lie on b=0 in the adopted grid.  Diagonal rotated
    ties belong to E/W.  This function deliberately rejects an unassigned
    equatorial pixel instead of changing the original standard partition.
    """
    L, B = np.broadcast_arrays(np.asarray(longitude), np.asarray(latitude))
    masks = {
        "std": {"WN": (L < 0) & (B > 0), "WS": (L < 0) & (B < 0),
                "EN": (L >= 0) & (B > 0), "ES": (L >= 0) & (B < 0)},
        "rot": {"N": B > np.abs(L), "S": B < -np.abs(L),
                "E": (L >= 0) & (L >= np.abs(B)),
                "W": (L < 0) & (L <= -np.abs(B))},
    }
    for name, pieces in masks.items():
        if not np.all(sum(x.astype(int) for x in pieces.values()) == 1):
            raise ValueError(f"{name} masks do not partition the source grid")
    return masks


# Public short name used by the renderer and independent geometry checks.
partition_masks = source_partition_masks


def _schematics():
    """Compact diagrams of the exact planar source partition, not a mask export.

    The drawing evaluates the original formulas on 1-degree display centres;
    straight boundary lines have the exact source-coordinate equations.  Plane
    shading is geometric only.  The registered source masks remain in Figure A.
    """
    from totani.templates import target_header
    from gcereduce.regrid import car_edges

    header = target_header(npix=120, pixel_deg=1)
    lon_edges = car_edges(header, 1, 120)
    lat_edges = car_edges(header, 2, 120)
    L, B = np.meshgrid((lon_edges[:-1] + lon_edges[1:]) / 2,
                       (lat_edges[:-1] + lat_edges[1:]) / 2)
    masks = partition_masks(L, B)
    out = {"longitude_edges": lon_edges.tolist(), "latitude_edges": lat_edges.tolist(),
           "observed": (np.abs(B) >= 10).tolist(),
           "description": "partition schematic evaluated from exact source formulas; plane exclusion shown; source holes shown in Figure A",
           "wcs_header": {key: header[key] for key in
                          ("NAXIS1", "NAXIS2", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2")},
           "source_exclusions_shown": False, "partition_before_response": True}
    for division in LABELS:
        image = np.zeros(L.shape, int)
        for i, label in enumerate(LABELS[division]):
            image[masks[division][label]] = i
        out[division] = image.tolist()
    return out


def validate_regional_basis(full, regional):
    """Check the *full observation domain*, before applying the fitted-cell cut."""
    J = np.asarray(full, float)
    H = np.asarray(regional, float)
    if H.shape[0] != 4 or H.shape[1:] != J.shape:
        raise ValueError("expected four regional response cubes on the full J grid")
    if not np.all(np.isfinite(J)) or not np.all(np.isfinite(H)):
        raise ValueError("nonfinite NFW basis")
    scale = max(float(np.max(np.abs(J))), np.finfo(float).tiny)
    if not np.allclose(H.sum(0), J, rtol=1e-6, atol=1e-6 * scale):
        raise ValueError("regional response cubes do not sum to the full NFW cube")
    return {"passed": True, "relative_max_absolute_difference":
            float(np.max(np.abs(H.sum(0) - J)) / scale),
            "domain": "all cached observation cells, before fitted-cell mask"}


def _paths(config):
    repo = Path(config["repo_root"])
    data = Path(config.get("data_root", repo))
    paths = config.get("result_paths", {})
    calibration = Path(paths.get("quadrant_calibration", DEFAULT_CALIBRATION))
    if not calibration.is_absolute():
        calibration = repo / calibration
    profiles = Path(config.get("profiles") or paths.get("profiles") or Path(config["bundle_dir"]) / "regional_profiles.json")
    if not profiles.is_absolute():
        profiles = repo / profiles
    cache = data / "results/fits/ranking_cells" / f"cells{float(config.get('cell_deg', 2)):g}_psfuc"
    row = tuple(config.get("row", BASELINE))
    return {"repo": repo, "data": data, "calibration": calibration,
            "profiles": profiles, "cache": cache,
            "library": cache / "library.npz", "quadrants": cache / "halo_quadrants.npz",
            "row": cache / f"row__{row[0]}__{row[1]}.npz",
            "flex": cache / "erosita_flex.npz"}


def _portable(path, config):
    p = Path(path)
    for root, prefix in ((Path(config["repo_root"]), ""),
                         (Path(config.get("data_root", config["repo_root"])), "data_root/")):
        try:
            return prefix + str(p.relative_to(root))
        except ValueError:
            pass
    return "configured_external/" + p.name


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _revision(repo):
    import subprocess
    result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo,
                            capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def inspect_regional(config):
    p = _paths(config)
    specs = (("regional calibration and common coefficients", "calibration", "-2 log likelihood and physical NFW coefficient", "existing final active result"),
             ("saved four-amplitude profiles", "profiles", "physical NFW coefficient", "required for points; recover only with explicit opt-in"),
             ("response-folded NFW and observation mask", "library", "expected counts per unit NFW coefficient", "required for owner recovery/basis validation"),
             ("source-partitioned response cubes", "quadrants", "expected counts per unit NFW coefficient", "required for owner recovery/basis validation"),
             ("baseline single-gas diffuse cubes", "row", "expected counts", "required for owner recovery"),
             ("twelve eROSITA brightness pieces", "flex", "expected counts per template weight", "required for owner recovery"))
    return [{"asset": a, "path": _portable(p[k], config), "available": p[k].is_file(),
             "units": u, "requirement": r} for a, k, u, r in specs]


def _assert_final(config, calibration):
    if float(config.get("cell_deg", 2)) != 2 or config.get("fine_tag", "psfuc") != "psfuc":
        raise ValueError("primary regional figure requires final 2-degree FINE_TAG=psfuc")
    if tuple(config.get("row", BASELINE)) != BASELINE:
        raise ValueError("primary regional figure requires SL_Z6_R30_T150_C2 gas1_ics1")
    if "pixelsrc" in str(calibration).lower() or "_psfuc_" not in calibration.name:
        raise ValueError("reject old pixel-source, CLEAN-proxy, or response-free results")


def _energy_edges(config):
    # Import the actual repository constants, not rounded manuscript energies.
    scripts = str(Path(config["repo_root"]) / "scripts")
    added = scripts not in sys.path
    if added:
        sys.path.insert(0, scripts)
    try:
        from reduce_totani import TOTANI_EDGES_GEV
        return list(TOTANI_EDGES_GEV[:10])
    finally:
        if added:
            sys.path.remove(scripts)


def _calibration(records, obs, cfg, division):
    if not np.isfinite(obs):
        raise ValueError("nonfinite observed calibration statistic")
    skies = [r for r in records if r.get("kind") == "null" and r.get("config") == cfg]
    ids = [r.get("sim") for r in skies]
    if any(not isinstance(i, int) or isinstance(i, bool) or i < 0 for i in ids):
        raise ValueError("independent simulated sky identifiers are missing or invalid")
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate simulated skies in fixed-partition calibration")
    if any(tuple(r.get("row", NULL_ROW)) != NULL_ROW for r in skies):
        raise ValueError("unexpected calibration null diffuse row")
    values = []
    for sky in skies:
        per_bin = np.asarray(sky.get(division), float)
        if per_bin.shape != (9,) or not np.all(np.isfinite(per_bin)):
            raise ValueError("calibration null statistic must have exactly nine finite bins")
        with np.errstate(over="ignore", invalid="ignore"):
            selected = float(per_bin[list(RESP_BINS)].sum())
        if not np.isfinite(selected):
            raise ValueError("nonfinite restricted-range null calibration statistic")
        values.append(selected)
    return {"exceedances": sum(x >= obs for x in values) if skies else None,
            "n_skies": len(skies) if skies else None,
            "null_row": list(NULL_ROW), "matched_to_data_row": False,
            "matching_config_partition_response_energy": True,
            "definition": "null>=observed", "partition": division,
            "selection": "fixed-partition Q(K=6); no maximization over partitions",
            "interpretation": "reference calibration from SO_Z6_R30_T150_C5 gas1_ics1; not a baseline-row null probability",
            "q_null_restricted": values,
            "seed_definition": "numpy.default_rng([13, sim, 20260929, config_index]); original config_index flat=0, structured=1, flat_flexlobes=2",
            "simulated_sky_ids": ids,
            "nuisance_refit": True}


def _basis_from_cache(paths):
    with np.load(paths["library"], allow_pickle=False) as z:
        J, keep = z["J"].astype(float), z["keep"].astype(bool)
    with np.load(paths["quadrants"], allow_pickle=False) as z:
        H = {d: np.stack([z[f"{d}__{label}"] for label in LABELS[d]]) for d in LABELS}
    checks = {d: validate_regional_basis(J, H[d]) for d in LABELS}
    return J, keep, H, checks


@contextmanager
def _analysis_context(config, paths):
    """Temporarily route existing analysis loaders to the configured owner's data."""
    scripts = str(paths["repo"] / "scripts")
    added = scripts not in sys.path
    if added:
        sys.path.insert(0, scripts)
    old_tag = os.environ.get("FINE_TAG")
    os.environ["FINE_TAG"] = "psfuc"
    from totani import ranking, psf
    old_run, old_cache = ranking.RUN, ranking.CACHE
    old_table = psf.PSF_TABLE_ULTRACLEAN
    ranking.RUN = paths["data"] / "runs/fff2de5bc7e0"
    ranking.CACHE = paths["data"] / "results/fits/ranking_cells"
    psf.PSF_TABLE_ULTRACLEAN = ranking.RUN / "psf_P8R3_ULTRACLEAN_V3_evtype3.fits"
    try:
        yield
    finally:
        ranking.RUN, ranking.CACHE = old_run, old_cache
        psf.PSF_TABLE_ULTRACLEAN = old_table
        if old_tag is None:
            os.environ.pop("FINE_TAG", None)
        else:
            os.environ["FINE_TAG"] = old_tag
        if added:
            sys.path.remove(scripts)


def recover_profile_intervals(joint_solver, T, n, signals, best_value, best_amplitudes):
    """Find Delta(-2lnL)=1 crossings, profiling *all* other signals and backgrounds.

    Roots use the original solver's count-space scale for numerical stability.
    Returned endpoints/scans inherit physical coefficients; no curvature or
    Hessian approximation is used.  A missing crossing fails explicitly.
    """
    from scipy.optimize import brentq

    scans, lower, upper = [], [], []
    scales = np.maximum([np.abs(t).sum() for t in signals], 1e-300)
    for i, (a, scale) in enumerate(zip(best_amplitudes, scales)):
        xbest = float(a * scale)
        cache = {xbest: float(best_value)}

        def profile(x):
            x = float(x)
            if x not in cache:
                value, _ = joint_solver(T, n, signals, fixed={i: x / scale})
                if not np.isfinite(value):
                    raise ValueError("nonfinite joint profile evaluation")
                if value < best_value - 0.05:
                    raise ValueError("profile found a better joint optimum; original fit must be investigated")
                cache[x] = float(value)
            return cache[x] - best_value - 1.0

        endpoints = []
        for side in (-1, 1):
            step = max(abs(xbest) * 0.15, np.sqrt(max(float(np.sum(n)), 1.0)) * 0.25, 1.0)
            bound = -1e14 if side < 0 else 1e14
            for _ in range(60):
                x = float(np.clip(xbest + side * step, -1e14, 1e14))
                if profile(x) >= 0:
                    lo, hi = sorted((xbest, x))
                    root = brentq(profile, lo, hi, xtol=1e-5, rtol=1e-10)
                    endpoints.append(root / scale)
                    break
                if x == bound:
                    raise ValueError("profile interval reaches original signed solver bound without a crossing")
                step *= 2
            else:
                raise ValueError("could not bracket 68% profile interval")
        lower.append(float(endpoints[0])); upper.append(float(endpoints[1]))
        scans.append({"amplitude": [x / scale for x in sorted(cache)],
                      "neg2loglike": [cache[x] for x in sorted(cache)],
                      "delta_neg2loglike": [cache[x] - best_value for x in sorted(cache)],
                      "profiled_parameters": "other three regional amplitudes and every foreground amplitude"})
    return lower, upper, scans


def _recover(config, paths, observed):
    required = ("library", "quadrants", "row", "flex")
    absent = [key for key in required if not paths[key].is_file()]
    if absent:
        raise FileNotFoundError("owner recovery requires existing final caches: " + ", ".join(absent))
    J, keep, H, basis_checks = _basis_from_cache(paths)
    with _analysis_context(config, paths):
        import quadrant_calibration as QC
        from ns_halo_spectra import joint
        from totani import fine
        from totani.likelihood import restrict_stack

        # The original stack constructor chooses exactly the original nuisance
        # names/freedoms.  No selection, library rebuilding, or simulation runs.
        QC._init(2.0)
        inputs = {}
        for cfg in CONFIGS:
            st = QC.stack(QC._W["lib"]["counts"], BASELINE, cfg)
            sr, _ = restrict_stack(st, keep, J)
            names = sorted(sr["templates"])
            inputs[cfg] = (sr, names)

        def recover_panel(task):
            cfg, division = task
            sr, names = inputs[cfg]
            saved = observed[cfg]
            coeff, common, vc, vr = [], [], [], []
            lower, upper = [[None] * 4 for _ in range(9)], [[None] * 4 for _ in range(9)]
            scans = {}
            for k in range(9):
                T = np.stack([sr["templates"][name][k] for name in names])
                n = sr["counts"][k]
                v1, a1 = joint(T, n, [J[k] * keep])
                sigs = [cube[k] * keep for cube in H[division]]
                v4, a4 = joint(T, n, sigs)
                if abs(v1 - saved["v_common"][k]) > 0.05 or abs((v1 - v4) - saved[division][k]) > 0.05:
                    raise ValueError(f"recovery differs from final saved likelihood at {cfg}/{division}/bin{k}; investigate inputs")
                if not np.isclose(a1[0], saved["common_amp"][k], rtol=5e-4, atol=1e-42):
                    raise ValueError(f"recovery differs from final saved common coefficient at {cfg}/bin{k}")
                coeff.append([float(x) for x in a4]); common.append(float(a1[0]))
                vc.append(float(v1)); vr.append(float(v4))
                if k == K20:
                    lower[k], upper[k], scans[str(k)] = recover_profile_intervals(joint, T, n, sigs, v4, a4)
            sums = H[division][:, :, keep].sum(-1).T
            return {"config": cfg, "partition": division, "labels": list(LABELS[division]),
                    "coefficients": coeff, "profile_lo": lower, "profile_hi": upper,
                    "common": common, "neg2loglike_common": vc, "neg2loglike_regions": vr,
                    "unit_count_sums": sums.tolist(), "profile_scans": scans,
                    "nuisance_names": names, "basis_sum_check": basis_checks[division],
                    "interval_method": "profile_delta_neg2loglike_1",
                    "interval_description": "joint profile likelihood; other regional and all foreground amplitudes refitted",
                    "solver": "scripts/ns_halo_spectra.py::joint"}

        workers = int(config.get("recovery_workers", 1))
        if workers < 1:
            raise ValueError("recovery_workers must be positive")
        with ThreadPoolExecutor(max_workers=workers) as pool:
            panels = list(pool.map(recover_panel, [(c, d) for c in CONFIGS for d in LABELS]))
    payload = {"schema_version": 1, "provenance": "scientific", "response": "P8R3_ULTRACLEAN_V3 evtype=3, exact-position source placement",
               "cell_deg": 2, "row": list(BASELINE), "source_partition_before_response": True,
               "recovery_revision": _revision(paths["repo"]),
               "active_result_provenance": "active psfuc result rerun after exact-position source rebuild; documented in validation/STATE_2026-09-28.md section 0d and psfuc_chain3_2026-10-01.sh; historical basename retained",
               "calibration_sha256": _sha(paths["calibration"]),
               "source_sha256": {key: _sha(paths[key]) for key in required}, "panels": panels}
    paths["profiles"].parent.mkdir(parents=True, exist_ok=True)
    temporary = paths["profiles"].with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    temporary.replace(paths["profiles"])
    return payload


def _validate_profiles(payload, calibration_path, paths):
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported recovered regional profile schema")
    if payload.get("provenance") != "scientific" or payload.get("row") != list(BASELINE):
        raise ValueError("regional profiles require scientific baseline provenance")
    if payload.get("cell_deg") != 2 or not payload.get("source_partition_before_response"):
        raise ValueError("regional profiles require exact 2-degree pre-response partition")
    if "ULTRACLEAN" not in payload.get("response", "") or "exact-position" not in payload.get("response", ""):
        raise ValueError("regional profiles require final ULTRACLEAN response and exact source placement")
    if payload.get("calibration_sha256") != _sha(calibration_path):
        raise ValueError("regional profiles refer to a different saved calibration result")
    if not all(key in payload.get("source_sha256", {}) for key in ("library", "quadrants", "row", "flex")):
        raise ValueError("regional profiles lack source-cache identity evidence")
    for key, digest in payload.get("source_sha256", {}).items():
        if key in paths and paths[key].is_file() and _sha(paths[key]) != digest:
            raise ValueError(f"regional profile source changed: {key}")
    entries = {(x["config"], x["partition"]): x for x in payload.get("panels", [])}
    if len(entries) != len(payload.get("panels", [])):
        raise ValueError("duplicate recovered regional panels")
    for key, panel in entries.items():
        if key[0] not in CONFIGS or key[1] not in LABELS:
            raise ValueError("unexpected recovered regional configuration")
        if panel.get("labels") != list(LABELS[key[1]]):
            raise ValueError("regional profile labels are out of inherited order")
        if not panel.get("basis_sum_check", {}).get("passed"):
            raise ValueError("regional profiles lack the full-domain regional basis sum check")
        if panel.get("interval_method") != "profile_delta_neg2loglike_1":
            raise ValueError("regional intervals must be joint profile likelihood intervals")
        if panel.get("solver") != "scripts/ns_halo_spectra.py::joint":
            raise ValueError("regional profile recovery must preserve the existing joint solver")
        a = np.asarray(panel["coefficients"], float)
        if a.shape != (9, 4) or not np.all(np.isfinite(a)):
            raise ValueError("nine-bin four-amplitude physical coefficient array required")
        for k in ("profile_lo", "profile_hi"):
            if np.shape(panel[k]) != (9, 4):
                raise ValueError("nine-bin four-amplitude profile endpoint array required")
        lo, hi = np.asarray(panel["profile_lo"][K20], float), np.asarray(panel["profile_hi"][K20], float)
        if not np.all(np.isfinite(lo)) or not np.all(np.isfinite(hi)) or np.any(lo > a[K20]) or np.any(hi < a[K20]):
            raise ValueError("20.8 GeV profile intervals missing or do not contain their estimates")
    return entries


def export_regional(config, recover_profiles=False):
    """Export final saved measurements; refit only on explicit recovery opt-in."""
    paths = _paths(config)
    _assert_final(config, paths["calibration"])
    result = {"status": "incomplete", "missing": [], "metadata": {
        "row": list(BASELINE), "cell_deg": 2, "fine_tag": "psfuc",
        "basis_units": "NFW J = integral rho^2 ds in GeV^2 cm^-5 sr^-1 (repository j_map convention); coefficient multiplies the bin-integrated exposure-folded J",
        "coefficient_units": "ph cm^3 s^-1 GeV^-2, integrated over each energy bin",
        "normalization": "inherited full NFW intensity basis; joint solver already divides independently normalized solver signals by their abs-sum",
        "source_partition_before_response": True,
        "basis_sum_validated": False,
        "masking_order": "source intensity partition; exposure*solid angle; PSF fold; extended-source exclusion when binning; fitted |b|>=10-degree cells",
        "response": "P8R3_ULTRACLEAN_V3 evtype=3; production exact-position point sources after 2026-10-01 rebuild",
        "bounds": {"foreground_solver_count_amplitudes": [0, 1e14], "regional_solver_count_amplitudes": [-1e14, 1e14],
                   "physical_regional_bounds": "[-1e14/c[e,r], +1e14/c[e,r]], c=sum(abs(H[e,r])) over fitted cells"},
        "fit_workload": "opt-in only: 54 four-amplitude bin fits and repeated common fits for six panels; 24 one-dimensional profiles at bin 5; no simulations",
        "cache_validation_method": "reproduce common coefficients and common/regions likelihood differences against active final calibration in every bin; tolerance 0.05 in -2 log L",
        "calibration_warning": "existing independent skies use SO_Z6_R30_T150_C5 gas1_ics1; transparently shown as reference truth-row calibration, not baseline-matched probability",
        "final_provenance_evidence": ["validation/STATE_2026-09-28.md section 0d", "results/fits/psfuc_chain3_2026-10-01.sh"]},
        "partitions": {"std": {"labels": list(LABELS["std"]), "longitude_signs": {"W": "l<0", "E": "l>=0"},
                               "construction": "l<0/l>=0 crossed with b>0/b<0 in source pixels"},
                       "rot": {"labels": list(LABELS["rot"]), "longitude_signs": {"W": "l<0", "E": "l>=0"},
                               "construction": "N:b>|l|; S:b<-|l|; E:l>=0 and l>=|b|; W:l<0 and l<=-|b|; diagonal ties E/W"}},
        "energy_edges_gev": _energy_edges(config), "illustrated_bin": K20,
        "restricted_bins": list(RESP_BINS), "panels": [], "source_paths": [], "diagnostics": {},
        "schematics": _schematics()}
    if not paths["calibration"].is_file():
        result["missing"].append("final active quadrant calibration/common-fit JSON")
        return result
    data = json.loads(paths["calibration"].read_text())
    if data.get("cell") != 2:
        raise ValueError("saved calibration spatial binning differs from requested 2 degrees")
    records = data["results"]
    observed = {r["config"]: r for r in records if r.get("kind") == "data" and tuple(r.get("row", ())) == BASELINE}
    if len(observed) != 3 or set(observed) != set(CONFIGS):
        raise ValueError("saved baseline results do not contain exactly the three approved configurations")
    result["source_paths"].append(_portable(paths["calibration"], config))
    result["metadata"]["calibration_sha256"] = _sha(paths["calibration"])
    payload = _recover(config, paths, observed) if recover_profiles else None
    if payload is None and paths["profiles"].is_file():
        payload = json.loads(paths["profiles"].read_text())
    recovered = _validate_profiles(payload, paths["calibration"], paths) if payload else {}
    if payload:
        result["source_paths"].append(_portable(paths["profiles"], config))
        result["metadata"]["recovery_revision"] = payload.get("recovery_revision")
        result["metadata"]["recovered_source_sha256"] = payload["source_sha256"]
    cache_basis = None
    if paths["library"].is_file() and paths["quadrants"].is_file():
        cache_basis = _basis_from_cache(paths)
        result["metadata"]["basis_sum_checks"] = cache_basis[3]
        result["metadata"]["basis_sum_validated"] = True
        result["metadata"]["max_relative_residual"] = max(x["relative_max_absolute_difference"] for x in cache_basis[3].values())
    elif recovered:
        checks = [p["basis_sum_check"] for p in recovered.values()]
        result["metadata"]["basis_sum_validated"] = all(x["passed"] for x in checks)
        result["metadata"]["max_relative_residual"] = max(x["relative_max_absolute_difference"] for x in checks)
    for cfg in CONFIGS:
        saved = observed[cfg]
        common = np.asarray(saved["common_amp"], float)
        vc = np.asarray(saved["v_common"], float)
        if common.shape != (9,) or vc.shape != (9,) or not np.all(np.isfinite(common)) or not np.all(np.isfinite(vc)):
            raise ValueError("invalid nine-bin saved common fit")
        for division in LABELS:
            q = np.asarray(saved[division], float)
            if q.shape != (9,) or not np.all(np.isfinite(q)) or np.min(q) < -0.001:
                raise ValueError("invalid common-minus-independent likelihood statistic")
            qs = float(q[list(RESP_BINS)].sum())
            entry = {"config": cfg, "partition": division, "labels": list(LABELS[division]),
                     "coefficients": [[None] * 4 for _ in range(9)],
                     "profile_lo": [[None] * 4 for _ in range(9)], "profile_hi": [[None] * 4 for _ in range(9)],
                     "common": common.tolist(), "neg2loglike_common": vc.tolist(),
                     "neg2loglike_regions": (vc - q).tolist(), "q_per_bin": q.tolist(),
                     "q20": float(q[K20]), "q_restricted": qs,
                     "calibration": _calibration(records, qs, cfg, division)}
            if entry["calibration"]["n_skies"] is None:
                entry["calibration"] = None
                entry["calibration_missing"] = "no matching configuration/partition independent simulated skies in saved ensemble"
            if (cfg, division) in recovered:
                panel = recovered[(cfg, division)]
                for key in ("coefficients", "profile_lo", "profile_hi", "profile_scans", "nuisance_names", "basis_sum_check", "interval_method", "solver", "unit_count_sums"):
                    entry[key] = panel[key]
                if not np.allclose(panel["common"], common, rtol=5e-4, atol=1e-42):
                    raise ValueError("saved recovered common coefficients differ from calibrated fit")
                qr = np.asarray(panel["neg2loglike_common"]) - np.asarray(panel["neg2loglike_regions"])
                if np.max(np.abs(qr - q)) > 0.05:
                    raise ValueError("saved recovered Q differs from final calibrated result")
                result["metadata"].setdefault("nuisance_names_by_config", {})[cfg] = panel["nuisance_names"]
            else:
                result["missing"].append(f"{cfg}/{division}: four physical coefficients and 68% joint profiles")
            result["panels"].append(entry)
    result["diagnostics"]["available_other_data_rows"] = [r for r in records if r.get("kind") == "data" and tuple(r.get("row", ())) != BASELINE]
    one_path = Path(config.get("result_paths", {}).get("quadrant_calibration_1deg", "results/fits/quadrant_calibration_cells1_psfuc_2026-09-29.json"))
    if not one_path.is_absolute():
        one_path = paths["repo"] / one_path
    result["diagnostics"]["one_degree"] = {"available": one_path.is_file(), "source": _portable(one_path, config)}
    if one_path.is_file():
        if "pixelsrc" in str(one_path).lower() or "_psfuc_" not in one_path.name:
            raise ValueError("diagnostic 1-degree comparison must use the active final psfuc result")
        other = json.loads(one_path.read_text())
        if other.get("cell") != 1:
            raise ValueError("1-degree diagnostic result has different spatial binning")
        diagnostic_panels = []
        for saved in other["results"]:
            if saved.get("kind") != "data" or tuple(saved.get("row", ())) != BASELINE:
                continue
            for division in LABELS:
                q = np.asarray(saved[division], float)
                qs = float(q[list(RESP_BINS)].sum())
                diagnostic_panels.append({"config": saved["config"], "partition": division,
                                          "common": saved["common_amp"], "q_per_bin": q.tolist(),
                                          "q20": float(q[K20]), "q_restricted": qs,
                                          "calibration": _calibration(other["results"], qs, saved["config"], division)})
        result["diagnostics"]["one_degree"]["panels"] = diagnostic_panels
        result["source_paths"].append(_portable(one_path, config))
    result["diagnostics"]["legacy_observed_quadrant_counts_rejected"] = "erosita_flex_ns partitions folded J in observation cells, unlike this calibrated pre-response test"
    if not result["missing"]:
        result["status"] = "complete"
    return result
