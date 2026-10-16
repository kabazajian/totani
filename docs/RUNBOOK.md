# Runbook: reproducing Abazajian, Li & Macias (2026)

This file lists, in dependency order, every step from public inputs to each
number, table and figure in the paper (`main_v2.tex`). Steps 1-2 build the
analysis-ready inputs; step 3 runs the fits. **With the Zenodo record of
analysis-ready inputs, a user can start at step 3** and needs neither the LAT
photon files nor Fermitools.

Section 3 maps each paper item to the script and output file that produce
it; every quantitative claim of the paper was traced to those files in an
audit on 2026-10-15, and four of the scripts were rerun in place with
identical results. Tables and figures are named by their LaTeX labels.

## 0. Environment

- Python >= 3.11 with `numpy`, `scipy`, `astropy` (`pip install -e .`);
  `matplotlib` for figures (`pip install -e '.[paper-figures]'`); `pytest`
  for the tests. Developed with the conda environment `gcepy`.
- Step 1 only: Fermitools 2.5.3 (`gtselect`, `gtmktime`, `gtbin`,
  `gtltcube`, `gtexpcube2`, `gtpsf`, and `gtsrcmaps` for the benchmark).
- A TeX installation for the two figures made by `scripts/paper_figures.py`.
- Memory: most steps < 4 GB; library builds and Totani-residual rebuilds up to
  10 GB (run one at a time); simulations scale with `--workers`. Set
  `OMP_NUM_THREADS=1` (and the OpenBLAS/MKL/vecLib equivalents) whenever
  `--workers > 1`.
- `python -m pytest tests/` checks the library (units, grids, statistics,
  vendored-code hashes).

## 1. Inputs

| input | source | read from |
|---|---|---|
| LAT photon and spacecraft files, weeks 009-791, Pass 8 ULTRACLEAN (front+back), zenith < 100 deg | Fermi Science Support Center | `scripts/reduce_totani.py` |
| 4FGL-DR4 catalog `gll_psc_v35.fit`; Galactic model `gll_iem_v07` (reference only); IRFs | FSSC | `~/data/fermi/catalog/`, `external/fssc/` |
| 2012 GALPROP grid (128 models) | galprop.stanford.edu/PaperIISuppMaterial (Ackermann et al. 2012) | `external/galprop_paperII/` |
| 2014 flat Fermi-bubbles template | FSSC public files (Ackermann et al. 2014) | `$GCE_FISHER_ROOT/data/external/FermiBubblesPublicFiles/` |
| PPPC4DMID gamma-ray spectra | Cirelli et al. 2011 | `external/PPPC4DMID/` |
| Haslam 408 MHz map | LAMBDA (Remazeilles et al. 2015) | `external/haslam/` |
| Cholis, Hu & Zhong bubble and Loop I templates | Zenodo 23199004 | `external/cholis2026/` |
| A2020 diffuse templates; solar inverse Compton template | Abazajian et al. 2020; Pohl et al. 2022; from the authors on request | `$GCE_FISHER_ROOT/data/external/Templates_arxiv2003.10416/`, `.../Pohl_etal_2022/` |
| structured Fermi-bubbles map | O. Macias (Fermi-LAT 2017 map, inpainted); from the authors on request | `external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits` |
| eROSITA bubbles brightness templates | O. Macias (eROSITA, Predehl et al. 2020; eastern half reconstructed); from the authors on request | `external/eROSITA_fullsky_lat_handoff/` |
| Galp21 GALPROP v57 outputs | Y. Li (this work; GALPROP, https://galprop.stanford.edu/); from the authors on request, not redistributed | `~/Data2/gce-fisher/Galp21_60deg/` (`totani.galp21.ARCHIVES`) |
| HESTIA particle data | Muru et al. 2025, from its authors on request | `$HESTIA_ROOT` (default `~/Data2/hestia_muru/`) |

`GCE_FISHER_ROOT` (default: a `gce-fisher` directory beside this repository) is
a data root with the layout `data/external/<set>/`; the sets read from it are
`Templates_arxiv2003.10416/`, `Pohl_etal_2022/`, `FermiBubblesPublicFiles/` and
`PPPC4DMID/` (`$PPPC4DMID_TABLE` overrides the last). Every script that lists
Galp21 rows filters them through `totani.ranking.available`, so with the
Galp21 cubes absent the pipeline runs on the 2012 grid alone and the Galp21
rows of the paper's tables are not reproduced.

## 2. Analysis-ready inputs (everything in `runs/` and `results/fits/ranking_cells/`)

Run in this order. Times are on the development machine (Mac Studio).

| step | command | produces | notes |
|---|---|---|---|
| 2a | `python scripts/reduce_totani.py --run` | `runs/fff2de5bc7e0/{ccube,expcube}.fits`, `config.json` | Fermitools; 960 x 960 CAR at 0.125 deg, 13 energy bins (9 fitted, 1.16-130 GeV) |
| 2b | `gtpsf` for P8R3_ULTRACLEAN_V3 evtype 3 (its arguments are recorded in `external/psf/psf_P8R3_ULTRACLEAN_V3_evtype3.json`), then the fit in `src/totani/psf.py`; no driver script, the table is in the Zenodo record | `runs/fff2de5bc7e0/psf_P8R3_ULTRACLEAN_V3_evtype3.{fits,json}` | three-Gaussian fit of the exposure-weighted PSF |
| 2c | `python scripts/rank_stage_a.py` | `results/fits/ranking_cells/{model__*,structure__*,shared__base}.npz` | 10-deg cell caches: 128 models x 2 decompositions, A2020, structures |
| 2d | `python scripts/galp21_ingest.py SA0_F98 SA0_R12 SA50_F98 SA50_R12 SA100_F98 SA100_R12` | `results/fits/ranking_cells/galp21_cubes__*.{npz,json}` | streams the archives; ~12 min and ~5 GB per model, one at a time |
| 2e | `python scripts/rebuild_residual_cache.py` | `results/fits/ranking_cells/totani_residual__*.npz` | Totani's residual bubble templates per configuration; ~75 min, ~8 GB |
| 2f | `python scripts/stage_2deg.py --build`; `python scripts/build_cells.py --cell 1` | `ranking_cells/cells{2,1}/library.npz`, `row__*.npz` | unconvolved 2-deg and 1-deg libraries and rows |
| 2g | `PSF_TAG=psfuc python scripts/build_cells.py --cell 2 --rows existing --psf-only` (and `--cell 1`) | `ranking_cells/cells{2,1}_psfuc/` | PSF-folded libraries |
| 2h | `FINE_TAG=psfuc python scripts/rebuild_pointsources.py --cell 2` (and `--cell 1`) | `cells{2,1}_psfuc/library.npz` | point sources at exact catalog positions (production since 2026-10-01) |
| 2i | built on first use by the step-3 scripts | `cells{2,1}_psfuc/halo_quadrants.npz`, `erosita_flex.npz`, `cholis2026.npz`, `sims_truths/` | quadrant pieces of the NFW template, eROSITA pieces, Cholis et al. components, simulation truths |
| 2j | `python scripts/hestia_roi60_maps.py` | `~/Data2/hestia_muru/totani_roi60/` | HESTIA J-maps (needs the HESTIA data) |
| 2k | `python scripts/stage_e.py --precompute` | `results/fits/stage_e/` truths | 10-deg simulation truths |

## 3. Paper items

Environment for every PSF-folded (2 deg, 1 deg) step: `FINE_TAG=psfuc`. The
10-deg original-binning results use no PSF by design. `--workers 6` is the
setting used; results do not depend on it.

### Main text

| paper item | command | output |
|---|---|---|
| Sec. IV A, `tab:res` | `scripts/oscar_conditional.py --cell 2` and `--cell 1` | `results/oscar_conditional_cells{2,1}_psfuc_2026-09-29.md` |
| Sec. IV A, `tab:sel`; App. F, `tab:twodeg` | `scripts/stage_select.py --cell {2,1} --bins {all,resp}` | `results/stage_select_cells{2,1}_psfuc_{all,resp}_2026-09-29.md` |
| choice of the best configurations | `scripts/rank_stage_d.py` (10 deg) | `results/rank_stage_d_2026-09-23_gas{1_ics1,5_ics3}_nores.md` |
| Sec. IV A reproduction (6412 counts, TS 172) | `scripts/bubble_totani.py` | `results/bubble_totani_2026-09-19.md` |
| Sec. IV A flat/structured at 10 deg (5006, 1575) | `scripts/totani_diffuse.py` | `results/totani_diffuse_2026-09-18.md` |
| Totani's six diffuse models | `scripts/totani_galprop_variants.py` | `results/totani_galprop_variants_2026-09-18.md` |
| Sec. IV B shape tests | `scripts/dm_shape_test.py --cell 2` | `results/dm_shape_test_cells2_psfuc_2026-10-01.md` |
| Secs. IV B, V D; `fig:specsims` | `scripts/sims_calibration_2deg.py --skies 50`; `scripts/sims_calibration_summary.py`; `scripts/oscar_figures.py` | `results/sims_calibration_cells2_psfuc_2026-09-29.md`; `results/oscar_spectra_vs_sims_2026-09-29_psfuc.pdf` |
| Sec. IV C, `tab:outside` | `scripts/inout_bubbles.py --cell 2`; `--cell 1 --common` | `results/inout_bubbles_cells2_psfuc_2026-09-29.md`, `results/inout_bubbles_cells1_psfuc_common_2026-10-01.md` |
| `fig:resid` | `scripts/component_residual_maps.py` | `results/component_residual_maps_cells2_psfuc_2026-10-03.pdf` |
| Sec. V A north/south | `scripts/ns_halo_spectra.py`; `scripts/flat2deg_checks.py --ns` | `results/ns_halo_spectra_psfuc_2026-09-27.md`; `results/flat2deg_checks_psfuc_ns_2026-09-26.md` |
| Sec. V A eROSITA, `fig:erositaflex`, `tab:erositaflex` | `scripts/erosita_flex_ns.py`, `erosita_flex_spectra.py`, `erosita_flex_truth.py`, `erosita_flex_score.py`, `erosita_flex_figure.py` | `results/erosita_flex_{ns,spectra,truth,score}_psfuc_2026-09-27.md`, `results/erosita_flex_figure_2026-09-27_psfuc.pdf` |
| Sec. V A East/West at 10 deg | `scripts/spatial_split_tests.py` | `results/spatial_split_tests_2026-09-24.md` |
| Sec. V B quadrant tests | `scripts/quadrant_calibration.py --cell 2 --skies 200`; `--cell 1 --skies 50` | `results/quadrant_calibration_cells{2,1}_psfuc_2026-09-29.md` |
| `fig:template_overview`, `fig:regional_nfw_normalizations` | `bash scripts/run_paper_figures.sh --data-root . --gce-fisher-root $GCE_FISHER_ROOT --erosita-root external/eROSITA_fullsky_lat_handoff --fermi-data-root ~/data/fermi --allow-nfw-constructor --recover-profiles` | `results/paper_figures/figures/{template_overview,regional_nfw_normalizations}.pdf` |
| Sec. V B HESTIA | `scripts/hestia_halo_fits.py` | `results/hestia_halo_fits_psfuc_2026-09-27.md` |
| Sec. V C, `tab:chz` | `scripts/cholis_bubbles.py --stage data`, `--stage sims --skies 50`, `--stage loop`, `--stage shape`, `--stage qnull_rows --skies 50`, `--stage data --cell 1`; `--render` | `results/cholis_bubbles_2026-10-11.md`; or `results/fits/cholis_chain_2026-10-11.sh` |
| Sec. V D, `tab:injection` | `scripts/injection_refined.py`; `scripts/joint_recovery_test.py` | `results/fits/injection_refined_2026-09-29.json`; `results/joint_recovery_test_2026-09-24.md` |

### Appendices

| paper item | command | output |
|---|---|---|
| App. A point sources, response checks | `scripts/audit_point_sources.py`; `scripts/srcmaps_benchmark.py`; `scripts/spillin_check.py --cell {2,1}` | `results/audit_point_sources_2026-10-01.md`, `results/srcmaps_benchmark_2026-10-01.md`, `results/spillin_check_cells{2,1}_psfuc_2026-10-01.md` |
| App. B Galp21 | `scripts/galp21_first_look.py` -> `scripts/galp21_summary.py`; `scripts/galp21_smoothing_test.py` | `results/galp21_summary_2026-09-23.md`; `results/fits/galp21_smoothing_test_2026-09-23.log` |
| App. C ranking and selection, `fig:ranking` (top row), `tab:selection`, `tab:galp21` | `scripts/rank_stage_b.py`, `rank_stage_c.py`, `rank_stage_d.py`; `scripts/plot_halo_significance.py` | `results/rank_stage_{b,c,d}_2026-09-2*.md`; `results/halo_significance_vs_fit_2026-09-24_notitle.pdf` |
| App. C, `tab:conditional` | `scripts/conditional_table.py` | `results/conditional_table_2026-09-24.md` |
| App. C, `fig:ranking` (nine panels: three bubbles templates at 10, 2 and 1 deg) | top row `scripts/plot_halo_significance.py --no-title`; lower rows `FINE_TAG=psfuc scripts/build_cells.py --cell {2,1} --rows grid` and `--rows galp21`, then `FINE_TAG=psfuc scripts/rank_fine.py --cell {2,1}` and `--summarize` (`results/fits/rank_fine_chain_2026-10-15.sh`, ~6 h per binning); `scripts/plot_rank_all.py --no-title` | `results/fits/rank_fine_cells{2,1}_psfuc_2026-10-15.json`, `results/rank_fine_cells{2,1}_psfuc_2026-10-15.md`, `results/rank_all_figure_2026-10-15_notitle.pdf` |
| App. C accompanying foregrounds, `tab:controls` | `scripts/rank_stage_b.py` (control sets, `audit_fix_chain_2026-10-15.sh`) -> `scripts/controls_table.py`; `scripts/bubble_ns_test.py` | `results/controls_table_2026-10-15.md`; `results/bubble_ns_test_2026-09-24.md` |
| App. C, `tab:family`, quadrant ranges | `scripts/quadrant_family.py`; `scripts/quadrant_range.py`; `scripts/quadrant_lr_{test,null,injected}.py` | `results/quadrant_family_2026-09-23.md`, `results/quadrant_range_2026-09-21_stagec0923.md`, `results/quadrant_lr_test_2026-09-24.md` |
| App. C, `tab:finer` | `scripts/finer_cells_test.py`; `scripts/finer_cells_null.py` | `results/finer_cells_test_2026-09-24.md` |
| App. D, `tab:nullsky` (and the simulated-sky ranking figure, not in the submitted paper) | `scripts/stage_e.py --n-sims 50`; `scripts/stage_e_anatomy.py`; `scripts/sim_ranking_figure.py` | `results/stage_e_2026-09-22.md`, `results/stage_e_anatomy_2026-09-23_null_*.md`, `results/fits/sim_ranking_figure_2026-09-24.pdf` |
| App. D coarse spectra (figure not in the submitted paper) | `scripts/halo_spectra_compare.py` | `results/halo_spectra_data_vs_sims_2026-09-24.pdf` |
| App. D finer-bin simulations | `scripts/stage_e_2deg.py`; `scripts/sims_calibration_2deg.py` | `results/stage_e_2deg_psfuc_2026-09-26.md` |
| App. D reconstruction | `scripts/reconstruction_check.py` | `results/fits/reconstruction_check_2026-09-24.json` |
| App. E eROSITA and HESTIA | `scripts/erosita_flex_*.py`; `scripts/ns_mismatch_null.py`; `scripts/hestia_halo_fits.py` | `results/erosita_flex_*_psfuc_2026-09-27.md`, `results/ns_mismatch_null_psfuc_2026-09-27.md`, `results/hestia_halo_fits_psfuc_2026-09-27.md` |
| App. F full-range fits | `scripts/stage_select.py --cell 2 --bins all` | `results/stage_select_cells2_psfuc_all_2026-09-29.md` |
| Numbers no other script writes (catalog share, eROSITA East share, CHZ cells, PSF containment, cell redistribution, pion spectrum, Galp21 gas ratio and IC/gas swaps, unit injection, spectral slopes, eROSITA E/W ratios, HESTIA shapes) | `scripts/paper_numbers.py` (~2 min, < 4 GB) | `results/paper_numbers_2026-10-15.md`, `results/fits/paper_numbers_2026-10-15.json` |

The production PSF-folded runs were made by `results/fits/psfuc_chain3_2026-10-01.sh`,
which runs the main-text steps above in order.
The audit reruns of 2026-10-15 (the representative five-ring row of `tab:outside`
at both binnings, and rows 1-3 of `tab:controls` with the corrected fitter) were made
by `results/fits/audit_fix_chain_2026-10-15.sh`.

## 4. Where to start with the archived inputs

The Zenodo record of analysis-ready inputs holds three archives, laid out
relative to the repository root:

| archive | contents | from step |
|---|---|---|
| `totani_runs.tar` | `runs/fff2de5bc7e0/`: counts cube, exposure, energy bounds, PSF table, selection configuration | 2a, 2b |
| `totani_ranking_cells.tar` | `results/fits/ranking_cells/`: the shared and structure 10-deg caches, Totani's residual constructions, the 2-deg and 1-deg libraries (PSF-folded and unconvolved) with their eROSITA, Cholis et al. and quadrant pieces, the simulation truths | 2c, 2e, 2f-2i, 2k |
| `totani_stage_e.tar` | `results/fits/stage_e/`: 10-deg simulation truths and skies | 2k |

Unpack all three into the repository root. The record contains no
GALPROP-derived file: the per-model caches (`model__*.npz` at 10 deg and
`row__*.npz` at 2 and 1 deg) are rebuilt from the public 2012 grid by steps
2c and 2f-2h (about 5 GB of downloads from galprop.stanford.edu and some
hours of CPU), after which step 3 runs for every 2012-grid row. The Galp21
rows need the Galp21 cubes (section 1) and step 2d.
