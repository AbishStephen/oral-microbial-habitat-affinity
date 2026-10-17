#!/bin/bash
#SBATCH --job-name=HOMD_FUNC_PROGRAM
#SBATCH --partition=computeshort
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/programs_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/10_additional_analyses_20261012/programs_%j.err
set -euo pipefail
python3 __PROJECT_ROOT__/study13/scripts/79_hypothesis_guided_programs.py
