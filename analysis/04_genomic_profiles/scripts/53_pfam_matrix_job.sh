#!/bin/bash
#SBATCH --job-name=HOMD13_PFMAT
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/pfmat_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/pfmat_%j.err
set -euo pipefail
python3 __PROJECT_ROOT__/study13/scripts/53_build_Pfam_matrices.py
