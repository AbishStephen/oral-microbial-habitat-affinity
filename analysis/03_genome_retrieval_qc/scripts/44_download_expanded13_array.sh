#!/bin/bash
#SBATCH --job-name=HOMD13_DL
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/download_%A_%a.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/download_%A_%a.err
set -euo pipefail
E=__PROJECT_ROOT__/study13
S=__WORK_ROOT__/HOMD_convergence_expanded13_work
DATASETS=datasets
id=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")
list=$S/download_lists/chunk_$id
zip=$S/download_chunks/chunk_$id.zip
out=$S/download_chunks/chunk_$id
mkdir -p "$S/download_chunks" "$S/raw_genomes"
rm -rf "$out" "$zip"
"$DATASETS" download genome accession --inputfile "$list" --include genome --filename "$zip" --no-progressbar
mkdir -p "$out"
unzip -q "$zip" -d "$out"
missing=0
while read -r g; do
  [ -n "$g" ] || continue
  f=$(find "$out/ncbi_dataset/data/$g" -maxdepth 1 -type f -name "*.fna" | head -1 || true)
  if [ -z "$f" ]; then
    echo "MISSING_FASTA $g" >&2
    missing=$((missing+1))
  else
    cp "$f" "$S/raw_genomes/$g.fna"
  fi
done < "$list"
echo "chunk=$id requested=$(wc -l < "$list") missing=$missing"
exit 0
