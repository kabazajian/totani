"""Stage E: the simulated-sky machinery reproduces the data pipeline on the
real counts, and the shared fraction-of-rows rule replays the data run."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_EDGES_GEV  # noqa: E402
from totani import ranking, selection, simsky  # noqa: E402

EDGES = TOTANI_EDGES_GEV[:10]
ROW = ("SO_Z6_R30_T150_C5", "gas1_ics1")
needs_pix = pytest.mark.skipif(not (simsky._row_dir(*ROW) / "done").exists()
                               or not (simsky.PIX / "_shared" / "done").exists(),
                               reason="run scripts/stage_e.py --precompute first")


@pytest.mark.parametrize("tag", ["gas1_ics1", "gas1_ics1_no_residual", "gas5_ics3", "gas5_ics3_no_residual"])
def test_fraction_rule_replays_the_data_run(tag):
    """Every recorded Stage C step is reproduced from its own per-row Δs."""
    p = ROOT / "results" / "fits" / f"rank_stage_c_2026-09-21_{tag}.json"
    d = json.loads(p.read_text())
    for st in d["consensus_trace"]:
        phys = [v for k, v in st["per_model"].items() if not k.startswith("gll_iem_v07")]
        rule = selection.fraction_rule_step(phys, list(st["clears_fraction"]), st["fraction_cutoff"])
        assert rule["best"] == st["best"] and rule["accepted"] == st["accepted"]
        assert rule["eligible"] == st["eligible"]
        assert rule["clears_fraction"] == st["clears_fraction"]


@pytest.mark.skipif(not (ranking.RUN / "ccube.fits").exists(), reason="reduction not present")
def test_shape_to_cells_is_the_stack_route():
    from totani.stack import cell_template_from_shape
    basis, L, B, _ = simsky._shape_basis(EDGES)
    shp = np.random.default_rng(3).random(L.shape)
    ref = cell_template_from_shape(ranking.RUN, EDGES, shp)
    assert np.abs(simsky.shape_to_cells(shp, EDGES) - ref).max() < 1e-12 * np.abs(ref).max()


@needs_pix
def test_fast_residual_rebuild_reproduces_the_cached_templates():
    """On the REAL counts, the precomputed-pixel rebuild is the data run's
    residual construction to rounding."""
    shared = simsky.load_shared()
    r = simsky.residual_cells(shared["counts"], simsky.load_row(*ROW), shared, EDGES)
    z = np.load(ranking.cache_path("totani_residual", f"{ROW[0]}__{ROW[1]}"))
    for k, n in (("pos", "bubble_totani_pos"), ("neg", "bubble_totani_neg")):
        assert np.abs(r[n] - z[k]).max() < 1e-12 * np.abs(z[k]).max()


@needs_pix
def test_truth_pixels_bin_to_the_cell_templates():
    from totani.cells import bin_to_cells
    shared = simsky.load_shared()
    sh = ranking.shared_cells(EDGES)
    for n in ("point_sources", "bubble_structured", "bubble_flat", "erosita_unif_N", "loopI_haslam_raw", "sun_moon_ic"):
        c = bin_to_cells(simsky.structure_pixels(n, EDGES, shared), shared["L"], shared["B"], mask=shared["mask"])
        ref = sh[n][2] if n in sh else ranking.structure_cells(n, EDGES)[2]
        assert np.abs(c - ref).max() < 1e-12 * np.abs(ref).max(), n


@needs_pix
def test_a_simulated_sky_keeps_the_real_disk_and_draws_the_rest():
    from totani.cells import bin_to_cells, halo_cell_mask
    shared = simsky.load_shared()
    sh = ranking.shared_cells(EDGES)
    truth = simsky.build_truth(*ROW, ("bubble_structured",), EDGES)
    cells, pix = simsky.draw(truth, sh["counts"], shared, np.random.default_rng(7))
    keep = halo_cell_mask()
    assert np.array_equal(cells[:, ~keep], sh["counts"][:, ~keep])                 # disk: the real counts
    assert np.array_equal(cells[2], bin_to_cells(pix, shared["L"], shared["B"], mask=shared["mask"]))
    # the analysis cells scatter about the truth by Poisson noise, not by more
    z = (cells[:, keep] - truth["cells"][:, keep]).sum(axis=1) / np.sqrt(truth["cells"][:, keep].sum(axis=1))
    assert np.all(np.abs(z) < 5), z
    # and the truth is a background-only fit: it reproduces the data's totals there
    assert np.allclose(truth["cells"][:, keep].sum(axis=1), sh["counts"][:, keep].sum(axis=1), rtol=1e-4)


def test_fit_ensemble_is_balanced_and_excludes_the_truths():
    import collections
    import stage_e
    codes = stage_e.fit_codes()
    assert len(set(codes)) == 16
    for i, n in enumerate((4, 4, 2, 2, 2)):
        assert set(collections.Counter(c.split("_")[i] for c in codes).values()) == {16 // n}
    assert not {t["model"] for t in stage_e.TRUTHS.values()} & set(codes)
