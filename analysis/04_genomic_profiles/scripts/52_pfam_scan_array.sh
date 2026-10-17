#!/bin/bash
set -euo pipefail

BASE=__PROJECT_ROOT__/study13
SCR=__WORK_ROOT__/HOMD_convergence_expanded13_work

PFAM=__WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm
CHUNK=$SCR/pfam/chunks/chunk_$(printf "%02d" $SLURM_ARRAY_TASK_ID).faa
OUT=$SCR/pfam/results

mkdir -p "$OUT"

echo "Task: $SLURM_ARRAY_TASK_ID"
echo "Input: $CHUNK"
grep -c '^>' "$CHUNK"
date

hmmscan \
  --cpu "${SLURM_CPUS_PER_TASK:-8}" \
  --cut_ga \
  --noali \
  --domtblout "$OUT/chunk_$(printf "%02d" $SLURM_ARRAY_TASK_ID).domtbl" \
  "$PFAM" \
  "$CHUNK" \
  > /dev/null

echo "Hits:"
grep -vc '^#' "$OUT/chunk_$(printf "%02d" $SLURM_ARRAY_TASK_ID).domtbl"

date
