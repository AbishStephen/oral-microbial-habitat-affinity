#!/bin/bash
#SBATCH --job-name=HOMD13_X5
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --output=__PROJECT_ROOT__/study13/logs/mmseqs_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/mmseqs_%j.err
set -euo pipefail
bash __PROJECT_ROOT__/study13/scripts/49_mmseqs_linclust.sh
