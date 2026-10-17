#!/bin/bash
#SBATCH --job-name=HOMD13_WG
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --output=__PROJECT_ROOT__/study13/logs/within_genus_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/within_genus_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
python3 __PROJECT_ROOT__/study13/scripts/61_expanded13_within_genus.py
