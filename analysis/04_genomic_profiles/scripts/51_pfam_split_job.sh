#!/bin/bash
#SBATCH --job-name=HOMD13_PFSPLIT
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/pfsplit_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/pfsplit_%j.err
set -euo pipefail
rm -rf __WORK_ROOT__/HOMD_convergence_expanded13_work/pfam/chunks __WORK_ROOT__/HOMD_convergence_expanded13_work/pfam/results
python3 __PROJECT_ROOT__/study13/scripts/51_split_pfam_fasta.py
