#!/bin/bash
# After the 1-deg grid chain: the 24 Galp21 rows at 1 deg, scored and summarized (Kev, 2026-10-15).
cd /Users/aba/work/totani
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=results/fits/rank_fine_logs_2026-10-15
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
until grep -q "^done 1" results/fits/rank_fine_chain1_2026-10-15.log; do sleep 60; done
PSF_TAG=psfuc OMP_NUM_THREADS=2 step build1_galp21 $PY scripts/build_cells.py --cell 1 --rows galp21
FINE_TAG=psfuc step score1_galp21 $PY scripts/rank_fine.py --cell 1 --workers 6
FINE_TAG=psfuc step summarize1b $PY scripts/rank_fine.py --cell 1 --summarize
echo "done 1b $(date)"
