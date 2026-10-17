#!/bin/bash
#SBATCH --job-name=HOMD14_LINK
#SBATCH --partition=computeshort
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=9G
#SBATCH --output=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/linkage14_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/linkage14_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python3 __PROJECT_ROOT__/study13/scripts/100_crossgenus_candidate_physical_linkage.py
