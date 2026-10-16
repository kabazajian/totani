#!/bin/zsh
# ULTRACLEAN PSF rerun (Kev, 2026-09-30): rebuild the PSF-folded libraries with the data's own PSF
# class and livetime, then rerun every PSF-dependent analysis under FINE_TAG=psfuc.
cd /Users/aba/work/totani
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=2
L=results/fits/psfuc_logs; mkdir -p $L
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
export PSF_TAG=psfuc
step build2 $PY scripts/build_cells.py --cell 2 --rows existing --psf-only
step build1 $PY scripts/build_cells.py --cell 1 --rows existing --psf-only
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
step sims2      $PY scripts/sims_calibration_2deg.py --skies 50 --workers 6
step simsum     $PY scripts/sims_calibration_summary.py
step inout2     $PY scripts/inout_bubbles.py --cell 2 --workers 6
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
echo "chain done $(date)"
