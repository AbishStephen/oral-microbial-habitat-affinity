#!/bin/bash
#SBATCH --job-name=HxlR_RefTree
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=8
#SBATCH --output=__PROJECT_ROOT__/study13/logs/HxlR_RefTree_%j.log
set -euo pipefail
S=__WORK_ROOT__/HOMD_HxlR_reference_20261014
O=__PROJECT_ROOT__/study13/16_HxlR_reference_phylogeny_20261014
MAFFT=mafft
RAXML=raxml-ng
echo "START 80 reference backbone" $(date)
"$MAFFT" --auto --thread 8 "$S/HxlR_X5_cluster_reps_plus_refs.faa" > "$S/HxlR_X5_reps_plus_refs.aligned.faa" 2> "$S/HxlR_X5_alignment.log"
echo "ALIGN80" $(date)
"$RAXML" --all --msa "$S/HxlR_X5_reps_plus_refs.aligned.faa" --model LG+G --bs-trees 100 --threads auto{8} --seed 20261014 --prefix "$S/HxlR_refX5" --force > "$S/HxlR_refX5_raxml.stdout.log" 2>&1
echo "TREE80" $(date)
cp "$S/HxlR_X5_reps_plus_refs.aligned.faa" "$O/HxlR_X5_reference_alignment.faa"
for ext in raxml.bestTree raxml.support raxml.log raxml.raxml.log raxml.bootstraps raxml.raxml.bestTree raxml.raxml.support raxml.raxml.log; do
  test ! -f "$S/HxlR_refX5.$ext" || cp "$S/HxlR_refX5.$ext" "$O/HxlR_refX5.$ext"
done
echo "DONE" $(date)
