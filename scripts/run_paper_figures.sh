#!/usr/bin/env bash
# Use Kev's existing analysis Python. This wrapper never creates a new ensemble.
set -uo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_python="${PYTHON:-python3}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# A writable cache also makes this usable in the restricted Codex environment.
task_cache="$(mktemp -d "${TMPDIR:-/tmp}/totani-figures.XXXXXXXX")"
trap 'rm -rf "$task_cache"' EXIT
export MPLCONFIGDIR="$task_cache/matplotlib" XDG_CACHE_HOME="$task_cache"
task_prepare=()
task_common=()
task_figure="all"
task_manuscript=""
task_proof_output="$task_root/results/paper_figures/proof"
while (($#)); do
  case "$1" in
    --recover-profiles|--allow-nfw-constructor)
      task_prepare+=("$1"); shift ;;
    --workers)
      task_prepare+=("$1" "${2:?--workers needs a value}"); shift 2 ;;
    --figure)
      task_figure="${2:?--figure needs a value}"; shift 2 ;;
    --manuscript)
      task_manuscript="${2:?--manuscript needs a .tex path}"; shift 2 ;;
    --proof-output)
      task_proof_output="${2:?--proof-output needs a directory}"; shift 2 ;;
    --config|--data-root|--gce-fisher-root|--erosita-root|--fermi-data-root|--bundle-dir|--output-dir|--profiles)
      task_common+=("$1" "${2:?option needs a value}"); shift 2 ;;
    --help|-h)
      "$task_python" "$task_root/scripts/paper_figures.py" --help
      printf '%s\n' 'Wrapper extras: --manuscript CURRENT.tex --proof-output DIR; all inputs are local.'
      exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
done
case "$task_figure" in
  all) task_figures=(spatial regional) ;;
  templates|spatial) task_figures=(spatial) ;;
  regional) task_figures=(regional) ;;
  *) printf 'Unknown figure: %s\n' "$task_figure" >&2; exit 2 ;;
esac
task_status=0
for task_current in "${task_figures[@]}"; do
  "$task_python" "$task_root/scripts/paper_figures.py" inspect --figure "$task_current" "${task_common[@]}" || task_status=2
  if "$task_python" "$task_root/scripts/paper_figures.py" prepare --figure "$task_current" "${task_common[@]}" "${task_prepare[@]}"; then
    "$task_python" "$task_root/scripts/paper_figures.py" check --figure "$task_current" "${task_common[@]}" || { task_status=2; continue; }
    "$task_python" "$task_root/scripts/paper_figures.py" render --figure "$task_current" "${task_common[@]}" || task_status=2
  else
    task_status=2
  fi
done
if [[ -n "$task_manuscript" && "$task_status" == 0 ]]; then
  "$task_python" "$task_root/scripts/paper_figure_proof.py" --manuscript "$task_manuscript" --proof-dir "$task_proof_output" "${task_common[@]}" || task_status=2
fi
exit "$task_status"
