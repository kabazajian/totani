#!/bin/bash
# ULTRACLEAN chain, rerun 2026-10-01 after the point-source template of the
# PSF-folded libraries was rebuilt with exact source positions
# (fine.point_sources_cells, exact=True; scripts/rebuild_pointsources.py; the
# previous libraries are kept as library_pixelsrc.npz). Same steps as
# psfuc_chain_2026-09-30.sh minus the library builds, plus the 2026-10-01
# additions (common-mask outside-footprint fit, spill-in check, DM shape test).
cd /Users/aba/work/totani
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
L=results/fits/psfuc_logs3; mkdir -p $L
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
# keep every pre-rebuild output for the before/after comparison
B=results/pixelsrc_2026-10-01; mkdir -p $B/fits
cp results/*_psfuc_*.md results/*_psfuc_*.json results/*_psfuc*.pdf results/*_psfuc*.png $B/ 2>/dev/null
cp results/fits/*_psfuc_*.json results/fits/*_psfuc*.json $B/fits/ 2>/dev/null
export FINE_TAG=psfuc
step cond2      $PY scripts/oscar_conditional.py --cell 2 --workers 6
step cond1      $PY scripts/oscar_conditional.py --cell 1 --workers 6
step stage_e    $PY scripts/stage_e_2deg.py --workers 6
step flat2      $PY scripts/flat2deg_checks.py --workers 6
step flat2ns    $PY scripts/flat2deg_checks.py --ns --workers 6
step nsspec     $PY scripts/ns_halo_spectra.py --workers 6
step nsnull     $PY scripts/ns_mismatch_null.py --workers 6
step eflex      $PY scripts/erosita_flex_ns.py --workers 6
step eflexspec  $PY scripts/erosita_flex_spectra.py --workers 6
step eflextruth $PY scripts/erosita_flex_truth.py --workers 6
step eflexfig   $PY scripts/erosita_flex_figure.py
step hestia     $PY scripts/hestia_halo_fits.py --workers 6
step quad2      $PY scripts/quadrant_calibration.py --cell 2 --skies 200 --workers 6
# the cell-level representative models are fitted with the library's point-source template: rebuild them
[ -d results/fits/ranking_cells/cells2_psfuc/sims_truths ] && mv results/fits/ranking_cells/cells2_psfuc/sims_truths results/fits/ranking_cells/cells2_psfuc/sims_truths_pixelsrc
step sims2      $PY scripts/sims_calibration_2deg.py --skies 50 --workers 6
step simsum     $PY scripts/sims_calibration_summary.py
step inout2     $PY scripts/inout_bubbles.py --cell 2 --workers 6
step inout1c    $PY scripts/inout_bubbles.py --cell 1 --common --workers 6
step resid2     $PY scripts/residual_maps.py --cell 2
step figs       $PY scripts/oscar_figures.py
step quad1      $PY scripts/quadrant_calibration.py --cell 1 --skies 50 --workers 6
step inout1     $PY scripts/inout_bubbles.py --cell 1 --workers 6
step resid1     $PY scripts/residual_maps.py --cell 1
step eflexscore $PY scripts/erosita_flex_score.py --workers 6
step sel2resp   $PY scripts/stage_select.py --cell 2 --bins resp --workers 6
step sel2all    $PY scripts/stage_select.py --cell 2 --bins all --workers 6
step sel1resp   $PY scripts/stage_select.py --cell 1 --bins resp --workers 6
step sel1all    $PY scripts/stage_select.py --cell 1 --bins all --workers 6
step spill2     $PY scripts/spillin_check.py --cell 2
step spill1     $PY scripts/spillin_check.py --cell 1
step shape2     $PY scripts/dm_shape_test.py --cell 2 --workers 6
step compare    $PY scripts/psf_class_compare.py
echo "chain done $(date)"
