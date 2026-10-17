#!/bin/bash
#SBATCH --job-name=SnoaL2_Clusters
#SBATCH --partition=compute
#SBATCH --time=00:25:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=2
#SBATCH --output=__PROJECT_ROOT__/study13/logs/SnoaL2_Clusters_%j.log
set -euo pipefail
M=mmseqs
S=__WORK_ROOT__/SnoaL2_subfamilies_20261016
mkdir -p "$S/tmp"
"$M" easy-cluster "$S/SnoaL2_all524_proteins.faa" "$S/SnoaL2_c50_cov80" "$S/tmp" --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads 2 >/dev/null 2>"$S/mmseqs.stderr.log"
echo DONE $(date)
