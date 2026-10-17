#!/bin/bash
#SBATCH --job-name=HOMD13_DLCHECK
#SBATCH --partition=compute
#SBATCH --time=00:15:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/download_check_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/download_check_%j.err
set -euo pipefail
E=__PROJECT_ROOT__/study13
S=__WORK_ROOT__/HOMD_convergence_expanded13_work
: > "$E/01_genome_manifests/MISSING_DOWNLOADS.txt"
while read -r g; do
  [ -s "$S/raw_genomes/$g.fna" ] || echo "$g" >> "$E/01_genome_manifests/MISSING_DOWNLOADS.txt"
done < "$E/01_genome_manifests/NEW_GCA_LIST.txt"
n=$(wc -l < "$E/01_genome_manifests/MISSING_DOWNLOADS.txt")
echo "missing_downloads=$n"
if [ "$n" -gt 0 ]; then
  cat "$E/01_genome_manifests/MISSING_DOWNLOADS.txt"
  exit 1
fi
