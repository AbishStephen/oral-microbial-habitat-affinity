#!/bin/bash
#SBATCH --job-name=HOMD13_X34
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/x34_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/x34_%j.err
set -euo pipefail
python3 __PROJECT_ROOT__/study13/scripts/54_build_X3_X4.py
