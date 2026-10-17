#!/bin/bash
set -euo pipefail

BASE=__PROJECT_ROOT__/study13
SCR=__WORK_ROOT__/HOMD_convergence_expanded13_work

MM=mmseqs
FASTA=$BASE/work/proteins/HOMD_1121_all_proteins.faa

OUT=$SCR/mmseqs
TMP=$SCR/tmp/mmseqs_x5

mkdir -p "$OUT" "$TMP"

echo "=== MMseqs2 ==="
"$MM" version

echo "=== INPUT ==="
grep -c '^>' "$FASTA"
du -sh "$FASTA"

echo "=== LINCLUST START ==="
date

"$MM" easy-linclust \
  "$FASTA" \
  "$OUT/HOMD13_x5" \
  "$TMP" \
  --min-seq-id 0.50 \
  -c 0.80 \
  --cov-mode 0 \
  --threads "${SLURM_CPUS_PER_TASK:-12}"

echo "=== LINCLUST COMPLETE ==="
date

echo
echo "=== OUTPUT ==="
ls -lh "$OUT"/HOMD13_x5*

echo
echo "Cluster memberships:"
wc -l "$OUT/HOMD13_x5_cluster.tsv"

echo
echo "Unique clusters:"
cut -f1 "$OUT/HOMD13_x5_cluster.tsv" | sort -u | wc -l

# copy essential outputs to permanent storage
mkdir -p "$BASE/03_orthogroups/mmseqs_x5"

cp "$OUT/HOMD13_x5_cluster.tsv" \
   "$BASE/03_orthogroups/mmseqs_x5/"

cp "$OUT/HOMD13_x5_rep_seq.fasta" \
   "$BASE/03_orthogroups/mmseqs_x5/"

echo
echo "=== PERMANENT COPY ==="
ls -lh "$BASE/03_orthogroups/mmseqs_x5/"
