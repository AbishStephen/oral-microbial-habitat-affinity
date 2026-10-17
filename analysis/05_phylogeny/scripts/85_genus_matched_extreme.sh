#!/bin/bash
#SBATCH --job-name=HOMD13_EXTREMES
#SBATCH --partition=computeshort
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/extreme13_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/extreme13_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
python3 __PROJECT_ROOT__/study13/scripts/85_genus_matched_extreme.py
