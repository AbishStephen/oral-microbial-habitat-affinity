#!/bin/bash
#SBATCH -J HOMD_PERM
#SBATCH -p computeshort
#SBATCH -t 00:20:00
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --array=0-199%12
#SBATCH -o __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/perm_%A_%a.out
#SBATCH -e __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/perm_%A_%a.err
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 __PROJECT_ROOT__/study13/scripts/106_nested_permutation.py "$SLURM_ARRAY_TASK_ID"
