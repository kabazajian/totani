#!/bin/bash
# The finer-binning ranking figure (Kev, 2026-10-15): score the cached rows under the
# three bubbles templates, build the rest of the 2012 grid at this cell size
# (unconvolved, then PSF-folded), score the new rows, summarize.
#   results/fits/rank_fine_chain_2026-10-15.sh 2      # or 1
cd /Users/aba/work/totani
C=$1
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=results/fits/rank_fine_logs_2026-10-15; mkdir -p $L
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
FINE_TAG=psfuc step score${C}_existing $PY scripts/rank_fine.py --cell $C --workers 6
PSF_TAG=psfuc OMP_NUM_THREADS=2 step build${C}_grid $PY scripts/build_cells.py --cell $C --rows grid
FINE_TAG=psfuc step score${C}_grid $PY scripts/rank_fine.py --cell $C --workers 6
FINE_TAG=psfuc step summarize${C} $PY scripts/rank_fine.py --cell $C --summarize
echo "done $C $(date)"
