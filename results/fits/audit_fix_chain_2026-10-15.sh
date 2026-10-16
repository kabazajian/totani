#!/bin/bash
# Audit follow-ups (Kev, 2026-10-15): the outside-bubbles fits with the representative
# five-ring configuration added, and the three tab:controls rows rerun with the
# post-2026-09-22 fitter.
cd /Users/aba/work/totani
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=results/fits/audit_fix_logs_2026-10-15; mkdir -p $L
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
FINE_TAG=psfuc step inout2  $PY scripts/inout_bubbles.py --cell 2 --workers 6
FINE_TAG=psfuc step inout1c $PY scripts/inout_bubbles.py --cell 1 --common --workers 6
S="bubble_structured,loopI_shell1,loopI_shell2"
step ctrl_solaric   $PY scripts/rank_stage_b.py --workers 6 --structures "$S,sun_moon_ic" --tag structured_shells_SOLARIC_ctrl_rerun2026-10-15
step ctrl_geomns    $PY scripts/rank_stage_b.py --workers 6 --structures "$S,erosita_outer_N,erosita_outer_S" --tag structured_shells_GEOMNS_ctrl_rerun2026-10-15
step ctrl_unifouter $PY scripts/rank_stage_b.py --workers 6 --structures "$S,erosita_unif_outer_N,erosita_unif_outer_S" --tag structured_shells_UNIFOUTER_ctrl_rerun2026-10-15
echo "done $(date)"
