#!/bin/bash
#SBATCH --job-name=SnoaL_Actino17
#SBATCH --partition=compute
#SBATCH --time=00:40:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --output=__PROJECT_ROOT__/study13/logs/SnoaL_Actino17_%j.log
set -euo pipefail
S=__WORK_ROOT__/SnoaL2_subfamilies_20261016
H=hmmscan
"$H" --cpu 4 --domE 0.01 --noali --domtblout "$S/domains/Actinomyces_326aa_all17_allPfam.domtblout" -o "$S/domains/Actinomyces_326aa_all17_allPfam.report" __WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm "$S/Actinomyces_326aa_all17.faa"
echo DONE $(date)
