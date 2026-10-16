#!/bin/bash
# Reproducibility spot check (2026-10-15): rerun representative scripts in place
# from the current code and caches; `git diff` then shows any change from the
# committed outputs.
cd /Users/aba/work/totani
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=results/fits/repro_logs_2026-10-15; mkdir -p $L
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
FINE_TAG=psfuc step cond2  $PY scripts/oscar_conditional.py --cell 2 --workers 6
FINE_TAG=psfuc step inout2 $PY scripts/inout_bubbles.py --cell 2 --workers 6
FINE_TAG=psfuc step shape2 $PY scripts/dm_shape_test.py --cell 2 --workers 6
step condtab $PY scripts/conditional_table.py
echo "done $(date)"
