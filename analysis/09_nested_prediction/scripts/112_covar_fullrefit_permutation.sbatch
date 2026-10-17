#!/bin/bash
#SBATCH -J HOMD_COVNULL
#SBATCH -p computeshort
#SBATCH -t 00:35:00
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --array=0-99%12
#SBATCH -o __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/covperm_%A_%a.out
#SBATCH -e __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/covperm_%A_%a.err
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 __PROJECT_ROOT__/study13/scripts/112_covar_fullrefit_permutation.py "$SLURM_ARRAY_TASK_ID"
