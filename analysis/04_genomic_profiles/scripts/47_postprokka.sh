#!/bin/bash
#SBATCH --job-name=HOMD13_CAT
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/catalogue_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/catalogue_%j.err
set -euo pipefail
E=__PROJECT_ROOT__/study13
PY=python3
mkdir -p "$E/work/proteins"
for sub in faa tsv gff; do
  n=$(find "$E/02_prokka/$sub" -maxdepth 1 \( -type f -o -type l \) | wc -l)
  echo "$sub=$n"
  [ "$n" -eq 1121 ]
done
"$PY" "$E/scripts/47_build_protein_catalogue.py"
cat $(awk 'NR>1{print $1}' "$E/01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv" | sed "s#^#$E/02_prokka/faa/#; s#$#.faa#") > "$E/work/proteins/HOMD_1121_all_proteins.faa"
echo proteins=$(grep -c '^>' "$E/work/proteins/HOMD_1121_all_proteins.faa")
"$PY" "$E/scripts/48_build_function_matrices.py"
