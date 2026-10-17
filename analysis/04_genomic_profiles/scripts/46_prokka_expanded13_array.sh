#!/bin/bash
#SBATCH --job-name=HOMD13_PROKKA
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=__PROJECT_ROOT__/study13/logs/prokka_%A_%a.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/prokka_%A_%a.err
set -euo pipefail
E=__PROJECT_ROOT__/study13
S=__WORK_ROOT__/HOMD_convergence_expanded13_work
CONT=__WORK_ROOT__/sintermedius_comparative_2026/containers/prokka_1.15.6.sif
g=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$E/01_genome_manifests/NEW_GCA_LIST.txt")
[ -n "$g" ]
fna="$S/raw_genomes/$g.fna"
[ -s "$fna" ]
tmp="$S/prokka_tmp/$g"
rm -rf "$tmp"; mkdir -p "$tmp"
tag=$(echo "$g" | tr -cd '[:alnum:]' | tail -c 18)
apptainer exec "$CONT" prokka --outdir "$tmp" --prefix "$g" --locustag "$tag" --kingdom Bacteria --cpus 2 --force --compliant --mincontiglen 200 "$fna"
cp "$tmp/$g.faa" "$E/02_prokka/faa/$g.faa"
cp "$tmp/$g.tsv" "$E/02_prokka/tsv/$g.tsv"
cp "$tmp/$g.gff" "$E/02_prokka/gff/$g.gff"
rm -rf "$tmp"
echo "DONE $g"
