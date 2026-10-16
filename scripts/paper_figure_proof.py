#!/usr/bin/env python3
"""Insert figures in a disposable copy of the supplied current manuscript."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from paper_figures import ROOT, configuration
from totani.paper_figures.bundle import BundleError, digest, load_bundle, write_json


def insertions(source, snippets):
    """Require current editorial anchors; never replace or rewrite original prose."""
    if "\\label{fig:template_overview}" in source or "\\label{fig:regional_nfw_normalizations}" in source:
        raise BundleError("Manuscript already contains one of these figure labels")
    blocks = []
    for name in ("TEMPLATE OVERVIEW", "REGIONAL NORMALIZATIONS"):
        start, stop = "% BEGIN " + name + " INSERTION", "% END " + name + " INSERTION"
        if snippets.count(start) != 1 or snippets.count(stop) != 1:
            raise BundleError("Ambiguous figure insertion snippets")
        blocks.append(snippets.split(start, 1)[1].split(stop, 1)[0])
    anchor = "\\section{Statistical analysis}"
    if source.count(anchor) != 1 or "\\label{sec:halo}" not in source or "\\label{sec:quadrantfinal}" not in source:
        raise BundleError("Use the current copyedited manuscript or paste the snippets at their documented anchors")
    first = "\n% BEGIN GENERATED FIGURE A\n" + blocks[0] + "\n% END GENERATED FIGURE A\n"
    edited = source.replace(anchor, first + anchor, 1)
    matches = list(re.finditer(r"more extreme significance\.\}\s*(?=\\om\{Using twelve brightness-weighted)", edited))
    if len(matches) != 1:
        raise BundleError("Regional insertion anchor differs; use the supplied focused snippet manually")
    second = "\n% BEGIN GENERATED FIGURE B\n" + blocks[1] + "\n% END GENERATED FIGURE B\n"
    position = matches[0].end()
    edited = edited[:position] + second + edited[position:]
    if edited.replace(first, "", 1).replace(second, "", 1) != source:
        raise BundleError("Insertion unexpectedly changed approved manuscript text")
    return edited


def run_checked(command, directory):
    process = subprocess.run(command, cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if process.returncode:
        raise BundleError("Manuscript compiler failed:\n" + process.stdout[-5000:])


def proof_directory(destination, manuscript):
    """Existing generated directories may be reused, but never follow links."""
    proof = Path(destination).expanduser().resolve()
    original = Path(manuscript).resolve().parent
    if proof == original or original.is_relative_to(proof):
        raise BundleError("Proof directory must be separate from the original manuscript directory")
    if proof.exists() and any(path.is_symlink() for path in proof.rglob("*")):
        raise BundleError("Proof destination contains symlinks; use a separate generated directory")
    return proof


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manuscript", required=True)
    ap.add_argument("--proof-dir", default=str(ROOT / "results/paper_figures/proof"))
    ap.add_argument("--config", default=str(ROOT / "configs/paper_figures.json"))
    for flag in ("data-root", "gce-fisher-root", "erosita-root", "fermi-data-root", "bundle-dir", "output-dir", "profiles"):
        ap.add_argument("--" + flag)
    ap.add_argument("--test-data", action="store_true")
    args = ap.parse_args()
    cfg = configuration(args.config, args)
    bundle = cfg["test_bundle_dir"] if args.test_data else cfg["bundle_dir"]
    graphics = cfg["test_output_dir"] if args.test_data else cfg["output_dir"]
    loaded = load_bundle(bundle, "all", allow_test=args.test_data)
    render_manifest = json.loads((graphics / "render_manifest.json").read_text())
    if render_manifest.get("provenance") != loaded["manifest"]["provenance"]:
        raise BundleError("Rendered figure provenance differs from bundle")
    for key in ("spatial", "regional"):
        if render_manifest.get("bundle_figures", {}).get(key) != loaded["manifest"]["figures"][key]:
            raise BundleError("Rerender this bundle before making a manuscript proof")
        for product in render_manifest["figures"][key]:
            if Path(product["path"]).name != product["path"] or digest(graphics / product["path"]) != product["sha256"]:
                raise BundleError("Rendered figure products are missing or changed")
    manuscript = Path(args.manuscript).expanduser().resolve()
    source = manuscript.read_text()
    edited = insertions(source, (ROOT / "docs/paper_figure_insertions.tex").read_text())
    proof = proof_directory(args.proof_dir, manuscript)
    proof.mkdir(parents=True, exist_ok=True)
    # Copy only article assets, leaving arbitrary project data and Git metadata.
    for asset in manuscript.parent.rglob("*"):
        if asset.is_file() and not asset.is_symlink() and asset.suffix.lower() in (".tex", ".bib", ".bbl", ".bst", ".pdf", ".png", ".jpg", ".jpeg"):
            if proof in asset.parents:
                continue
            target = proof / asset.relative_to(manuscript.parent)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(asset, target)
    proof_tex = proof / "figure_proof.tex"
    proof_tex.write_text(edited)
    (proof / "figures").mkdir(exist_ok=True)
    for basename in ("template_overview", "regional_nfw_normalizations"):
        shutil.copy2(graphics / f"{basename}.pdf", proof / "figures" / f"{basename}.pdf")
    if not shutil.which("pdflatex"):
        raise BundleError(f"Figure proof source is ready at {proof_tex}; install/use the owner's existing REVTeX TeX environment to compile")
    command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape", proof_tex.name]
    if manuscript.with_suffix(".bbl").is_file():
        shutil.copy2(manuscript.with_suffix(".bbl"), proof_tex.with_suffix(".bbl"))
    run_checked(command, proof)
    if shutil.which("bibtex") and (proof / "refs.bib").exists():
        run_checked(["bibtex", proof_tex.stem], proof)
    for _ in range(2):
        run_checked(command, proof)
    if manuscript.read_text() != source:
        raise BundleError("Original manuscript changed during proof generation")
    log = proof_tex.with_suffix(".log").read_text(errors="replace")
    warnings = [line for line in log.splitlines() if "Overfull" in line or "undefined" in line]
    write_json(proof / "proof_report.json", {"source_name": manuscript.name, "source_sha256": digest(manuscript),
        "provenance": loaded["manifest"]["provenance"], "approved_text_preserved": True,
        "insertions": ["fig:template_overview", "fig:regional_nfw_normalizations"],
        "width_inches": cfg["width_inches"], "compiler_warnings": warnings,
        "inspection_required": "Inspect the two figure pages at their final size and in grayscale; fixture proofs establish layout only"})
    print(proof_tex.with_suffix(".pdf"))


if __name__ == "__main__":
    try:
        main()
    except (BundleError, OSError, ValueError, KeyError) as exc:
        print(f"paper_figure_proof: {exc}", file=sys.stderr)
        raise SystemExit(2)
