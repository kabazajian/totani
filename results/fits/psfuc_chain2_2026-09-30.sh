#!/bin/zsh
# Rerun the four simulation scripts whose seeds used Python's per-process hash() (now zlib.crc32).
cd /Users/aba/work/totani
PY=/opt/homebrew/Caskroom/miniforge/base/envs/gcepy/bin/python
export OMP_NUM_THREADS=2 FINE_TAG=psfuc
L=results/fits/psfuc_logs
until grep -q "chain done" results/fits/psfuc_chain_2026-09-30.log; do sleep 60; done
step() { local name=$1; shift; echo "start $name $(date)"; "$@" > $L/$name.log 2>&1; echo "end $name exit $? $(date)"; }
step stage_e_seed $PY scripts/stage_e_2deg.py --workers 6
step flat2_seed   $PY scripts/flat2deg_checks.py --workers 6
step flat2ns_seed $PY scripts/flat2deg_checks.py --ns --workers 6
step nsspec_seed  $PY scripts/ns_halo_spectra.py --workers 6
step nsnull_seed  $PY scripts/ns_mismatch_null.py --workers 6
echo "chain2 done $(date)"
