#!/bin/bash
#SBATCH --job-name=HOMD13_GBAL_META
#SBATCH --partition=computeshort
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --output=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/crossgenus_96_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/crossgenus_96_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=6 OMP_NUM_THREADS=6 MKL_NUM_THREADS=6
python3 __PROJECT_ROOT__/study13/scripts/96_crossgenus_balanced_permutation_scan.py 4096
