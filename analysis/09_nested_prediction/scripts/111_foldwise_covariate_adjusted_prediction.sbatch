#!/bin/bash
#SBATCH -J HOMD_COVPRED
#SBATCH -p computeshort
#SBATCH -t 00:45:00
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -o __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/adjusted_111_%j.out
#SBATCH -e __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/adjusted_111_%j.err
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 __PROJECT_ROOT__/study13/scripts/111_foldwise_covariate_adjusted_prediction.py
