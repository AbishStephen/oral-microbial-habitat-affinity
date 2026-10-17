#!/bin/bash
#SBATCH --job-name=SnoaL_long_allPfam
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --output=__PROJECT_ROOT__/study13/logs/SnoaL_long_allPfam_%j.log
set -euo pipefail
S=__WORK_ROOT__/SnoaL2_subfamilies_20261016
H=hmmscan
"$H" --cpu 4 --domE 0.01 --noali --domtblout "$S/domains/long5_allPfam.domtblout" -o "$S/domains/long5_allPfam.report" __WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm "$S/SnoaL2_long_cluster_reps.faa"
echo DONE $(date)
