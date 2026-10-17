#!/bin/bash
#SBATCH --job-name=HOMD13_NULL
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/null_%A_%a.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/null_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python3 __PROJECT_ROOT__/study13/scripts/60_expanded13_permutation_null.py --chunk "$SLURM_ARRAY_TASK_ID" --per-chunk 5 --seed 20261011
