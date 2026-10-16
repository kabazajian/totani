#!/usr/bin/env python3
"""Inspect, export, validate, and render the two approved manuscript figures."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from totani.paper_figures.bundle import BundleError, digest, load_bundle, write_json  # noqa: E402
from totani.paper_figures import SCHEMA_VERSION  # noqa: E402


def configuration(path, args):
    path = Path(path).resolve()
    if path.suffix == ".json":
        cfg = json.loads(path.read_text())
    else:
        import yaml  # optional: default configuration needs no YAML dependency
        cfg = yaml.safe_load(path.read_text())
    cfg["repo_root"] = ROOT
    for name in ("data_root", "gce_fisher_root", "erosita_root", "fermi_data_root", "nfw_intensity_path", "bundle_dir", "output_dir",
                 "test_bundle_dir", "test_output_dir", "profiles"):
        value = getattr(args, name, None) or os.environ.get("TOTANI_" + name.upper()) or cfg.get(name)
        if value:
            value = Path(os.path.expandvars(str(value))).expanduser()
            cfg[name] = value if value.is_absolute() else ROOT / value
    if os.environ.get("GCE_FISHER_ROOT") and not cfg.get("gce_fisher_root"):
        cfg["gce_fisher_root"] = Path(os.environ["GCE_FISHER_ROOT"]).expanduser()
    if getattr(args, "allow_nfw_constructor", False):
        cfg["allow_nfw_constructor"] = True
    if getattr(args, "workers", None) is not None:
        cfg["recovery_workers"] = args.workers
    return cfg


def _revision():
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def export_bundle(cfg, selected, fixture=False, recover=False):
    import numpy as np
    directory = cfg["test_bundle_dir"] if fixture else cfg["bundle_dir"]
    directory.mkdir(parents=True, exist_ok=True)
    if fixture:
        from totani.paper_figures.fixtures import synthetic_figures
        spatial, arrays, regional = synthetic_figures()
    manifest_path = directory / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("provenance") != ("test" if fixture else "scientific"):
            raise BundleError("Do not mix test and scientific bundles")
    else:
        manifest = {"schema_version": SCHEMA_VERSION, "provenance": "test" if fixture else "scientific", "figures": {}}
    manifest["repository"] = {"remote": "https://github.com/kabazajian/totani", "revision": _revision()}
    for figure in selected:
        if figure == "spatial":
            if not fixture:
                from totani.paper_figures.spatial import export_spatial
                spatial, arrays = export_spatial(cfg)
            write_json(directory / "spatial.json", spatial)
            np.savez_compressed(directory / "spatial.npz", **arrays)
            manifest["figures"][figure] = {
                "status": spatial["status"],
                "metadata": {"path": "spatial.json", "sha256": digest(directory / "spatial.json")},
                "arrays": {"path": "spatial.npz", "sha256": digest(directory / "spatial.npz")}}
        else:
            if not fixture:
                from totani.paper_figures.regional import export_regional
                regional = export_regional(cfg, recover_profiles=recover)
            write_json(directory / "regional.json", regional)
            manifest["figures"][figure] = {"status": regional["status"], "metadata": {
                "path": "regional.json", "sha256": digest(directory / "regional.json")}}
        # A completed independent export survives an unavailable later figure.
        write_json(manifest_path, manifest)
    write_json(manifest_path, manifest)
    return directory, manifest


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("inspect", "prepare", "render", "check", "fixture"))
    ap.add_argument("--config", default=str(ROOT / "configs/paper_figures.json"))
    ap.add_argument("--figure", choices=("all", "templates", "spatial", "regional"), default="all")
    ap.add_argument("--data-root")
    ap.add_argument("--gce-fisher-root")
    ap.add_argument("--erosita-root")
    ap.add_argument("--fermi-data-root", help="Root with catalog/gll_psc_v35.fit; default FERMI_DATA")
    ap.add_argument("--bundle-dir")
    ap.add_argument("--output-dir")
    ap.add_argument("--profiles")
    ap.add_argument("--recover-profiles", action="store_true", help="Opt in to the original common/four-amplitude fits and profile intervals; prepare only")
    ap.add_argument("--allow-nfw-constructor", action="store_true", help="Opt in to regenerating the adopted input J map on the owner's machine")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--test-data", action="store_true", help="Use only the segregated, watermarked test bundle")
    args = ap.parse_args(argv)
    if args.recover_profiles and args.command != "prepare":
        ap.error("--recover-profiles is allowed only for prepare")
    if args.test_data and args.command not in ("render", "check"):
        ap.error("--test-data is allowed only for render/check; use fixture to create it")
    cfg = configuration(args.config, args)
    if cfg.get("recovery_workers", 1) < 1:
        ap.error("--workers must be positive")
    figure = "spatial" if args.figure == "templates" else args.figure
    selected = ("spatial", "regional") if figure == "all" else (figure,)
    if args.command == "inspect":
        inventory = []
        if "spatial" in selected:
            from totani.paper_figures.spatial import inspect_spatial
            inventory.extend(inspect_spatial(cfg))
        if "regional" in selected:
            from totani.paper_figures.regional import inspect_regional
            inventory.extend(inspect_regional(cfg))
        print(json.dumps({"revision": _revision(), "assets": inventory}, indent=2, default=str))
        return 0
    if args.command in ("prepare", "fixture"):
        directory, manifest = export_bundle(cfg, selected, fixture=args.command == "fixture", recover=args.recover_profiles)
        print(directory / "manifest.json")
        incomplete = [key for key in selected if manifest["figures"][key]["status"] != "complete"]
        for key in incomplete:
            metadata = json.loads((directory / manifest["figures"][key]["metadata"]["path"]).read_text())
            print(f"{key}: incomplete: " + "; ".join(metadata.get("missing", [])), file=sys.stderr)
        return 2 if incomplete else 0
    directory = cfg["test_bundle_dir"] if args.test_data else cfg["bundle_dir"]
    loaded = load_bundle(directory, figure, allow_test=args.test_data)
    output = cfg["test_output_dir"] if args.test_data else cfg["output_dir"]
    if args.command == "check":
        print(json.dumps({"bundle": str(directory), "provenance": loaded["manifest"]["provenance"],
                          "checked": list(selected), "result": "numeric scientific checks passed" if not args.test_data else "TEST schema checks passed; no scientific validation"}, indent=2))
        return 0
    from totani.paper_figures.rendering import render_spatial, render_regional
    render_manifest_path = Path(output) / "render_manifest.json"
    render_manifest = json.loads(render_manifest_path.read_text()) if render_manifest_path.exists() else {"figures": {}, "bundle_figures": {}}
    if render_manifest.get("provenance", loaded["manifest"]["provenance"]) != loaded["manifest"]["provenance"]:
        raise BundleError("Do not mix test and scientific render directories")
    render_manifest.update(provenance=loaded["manifest"]["provenance"],
                           bundle_manifest_sha256=digest(directory / "manifest.json"),
                           width_inches=cfg["width_inches"])
    for key in selected:
        if key == "spatial":
            paths = render_spatial(*loaded[key], output, test_data=args.test_data, width_inches=cfg["width_inches"])
        else:
            paths = render_regional(loaded[key], output, test_data=args.test_data, width_inches=cfg["width_inches"])
        for path in paths:
            print(path)
        render_manifest["figures"][key] = [{"path": Path(path).name, "sha256": digest(path)} for path in paths]
        render_manifest.setdefault("bundle_figures", {})[key] = loaded["manifest"]["figures"][key]
    write_json(render_manifest_path, render_manifest)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (BundleError, FileNotFoundError, KeyError, ValueError) as exc:
        print(f"paper_figures: {exc}", file=sys.stderr)
        raise SystemExit(2)
