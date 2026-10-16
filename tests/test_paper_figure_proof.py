"""Insertion is bounded and cannot silently alter approved manuscript text."""
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from paper_figure_proof import insertions, proof_directory
from paper_figures import configuration
from totani.paper_figures.bundle import BundleError


def _article():
    return """\\title{Approved title}
\\begin{abstract}Approved abstract.\\end{abstract}
\\label{sec:halo}
Approved NFW equation and numerical result.
\\section{Statistical analysis}
Approved method paragraph.
\\label{sec:quadrantfinal}
\\om{rather than establishing a more extreme significance.}

\\om{Using twelve brightness-weighted eROSITA pieces, the existing result.}
\\bibliography{refs}
"""


def test_insertions_preserve_all_original_content_exactly():
    source = _article()
    updated = insertions(source, (ROOT / "docs/paper_figure_insertions.tex").read_text())
    for label in ("fig:template_overview", "fig:regional_nfw_normalizations"):
        assert updated.count("\\label{" + label + "}") == 1
    import re
    restored = re.sub(r"\n% BEGIN GENERATED FIGURE [AB]\n.*?\n% END GENERATED FIGURE [AB]\n", "", updated, flags=re.S)
    assert restored == source


def test_ambiguous_or_already_modified_manuscript_is_rejected():
    snippets = (ROOT / "docs/paper_figure_insertions.tex").read_text()
    with pytest.raises(BundleError, match="already contains"):
        insertions(_article() + "\\label{fig:template_overview}", snippets)
    with pytest.raises(BundleError, match="current copyedited"):
        insertions(_article().replace("\\section{Statistical analysis}", "\\section{Old methods}"), snippets)
    with pytest.raises(BundleError, match="Regional insertion anchor"):
        insertions(_article().replace("more extreme significance", "some revised significance"), snippets)


def test_proof_directory_cannot_overwrite_original_or_follow_links(tmp_path):
    original = tmp_path / "article"
    original.mkdir()
    manuscript = original / "approved.tex"
    manuscript.write_text(_article())
    with pytest.raises(BundleError, match="separate"):
        proof_directory(original, manuscript)
    with pytest.raises(BundleError, match="separate"):
        proof_directory(tmp_path, manuscript)
    proof = tmp_path / "proof"
    proof.mkdir()
    (proof / "figure_proof.tex").symlink_to(manuscript)
    with pytest.raises(BundleError, match="symlinks"):
        proof_directory(proof, manuscript)
    assert manuscript.read_text() == _article()


def test_owner_paths_are_resolved_independently_of_invoking_directory(tmp_path, monkeypatch):
    from argparse import Namespace
    import json
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({"data_root": ".", "nfw_intensity_path": "${TOTANI_TEST_INPUT}/J.fits",
                                       "fermi_data_root": "owner_fermi"}))
    monkeypatch.setenv("TOTANI_TEST_INPUT", "owner_inputs")
    monkeypatch.chdir(tmp_path)
    config = configuration(config_file, Namespace())
    assert config["data_root"] == ROOT
    assert config["nfw_intensity_path"] == ROOT / "owner_inputs/J.fits"
    assert config["fermi_data_root"] == ROOT / "owner_fermi"
