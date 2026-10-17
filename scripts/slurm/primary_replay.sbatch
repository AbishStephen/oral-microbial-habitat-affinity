#!/usr/bin/env bash
#SBATCH --job-name=HOMD13_replay
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
repo=${1:?repository checkout required}
run_dir=${2:?separate run directory required}
python_bin=${3:-python3}
cd -- "$repo"
export PYTHONPATH="$repo/src"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
"$python_bin" -m homd_signatures.reproduce --output "$run_dir"
