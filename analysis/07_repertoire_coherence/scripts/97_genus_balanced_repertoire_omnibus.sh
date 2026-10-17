#!/bin/bash
#SBATCH --job-name=HOMD13_REPERTOIRE
#SBATCH --partition=computeshort
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/repertoire_97_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/repertoire_97_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8
python3 __PROJECT_ROOT__/study13/scripts/97_genus_balanced_repertoire_omnibus.py
