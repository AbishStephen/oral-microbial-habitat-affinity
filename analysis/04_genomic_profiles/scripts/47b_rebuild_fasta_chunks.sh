#!/bin/bash
#SBATCH --job-name=HOMD13_FASTA
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/fasta_rebuild_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/fasta_rebuild_%j.err
set -euo pipefail
E=__PROJECT_ROOT__/study13
OUT=$E/work/proteins/HOMD_1121_all_proteins.faa
: > "$OUT"
tail -n +2 "$E/01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv" | cut -f1 | while read -r g; do
  f="$E/02_prokka/faa/$g.faa"
  [ -s "$f" ] || { echo "Missing FAA: $g" >&2; exit 1; }
  cat "$f" >> "$OUT"
done
n=$(grep -c '^>' "$OUT")
echo "proteins=$n"
[ "$n" -eq 2386061 ]
ls -lh "$OUT"
rm -rf __WORK_ROOT__/HOMD_convergence_expanded13_work/pfam/chunks __WORK_ROOT__/HOMD_convergence_expanded13_work/pfam/results
python3 "$E/scripts/51_split_pfam_fasta.py"
echo chunks=$(find __WORK_ROOT__/HOMD_convergence_expanded13_work/pfam/chunks -name 'chunk_*.faa' | wc -l)
