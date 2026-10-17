#!/bin/bash
#SBATCH --job-name=HxlR_HOMD631
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=8
#SBATCH --output=__PROJECT_ROOT__/study13/logs/HxlR_HOMD631_%j.log
set -euo pipefail
T=__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014
F="$T/HOMD_631_HMTrepresentatives_2025_PROKKA.faa"
H=hmmsearch
O=__PROJECT_ROOT__/study13/14_HxlR_fullHOMD_controls_20261014
for P in PF01638 PF12680; do
    echo "START $P $(date)"
    "$H" --cpu 8 --cut_ga --noali --domtblout "$O/HOMD631_${P}_domtblout.tsv" -o "$O/HOMD631_${P}_summary.txt" "$T/$P.hmm" "$F"
    echo "DONE $P $(date)"
done
