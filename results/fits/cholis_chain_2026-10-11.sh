#!/bin/bash
# Kev 2026-10-11: follow-ups 1-4 to the Cholis, Hu & Zhong comparison.
set -u
cd /Users/aba/work/totani
export FINE_TAG=psfuc OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
P=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
for args in "--stage loop" "--stage shape" "--stage qnull_rows --skies 50" "--stage data --cell 1"; do
  echo "=== $args $(date)"
  /usr/bin/time -l $P scripts/cholis_bubbles.py $args --workers 6 2>&1 | grep -vE "RuntimeWarning|w = np.exp" | grep -E "^\[|built|Error|Traceback|real|maximum resident|File|line" | tail -8
done
echo "=== done $(date)"
