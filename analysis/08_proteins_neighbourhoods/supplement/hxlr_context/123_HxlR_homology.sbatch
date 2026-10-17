#!/bin/bash
#SBATCH --job-name=HxlR_Homology
#SBATCH --partition=compute
#SBATCH --time=01:15:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=__PROJECT_ROOT__/study13/logs/HxlR_Homology_%j.log
set -euo pipefail
B=__PROJECT_ROOT__/study13
O=$B/13_HxlR_broad_and_cassette_20261014
M=mmseqs
T=__WORK_ROOT__/HOMD_HxlR_20261014_tmp
mkdir -p "$T" "$O"
echo START $(date)
"$M" easy-search "$O/shared_1789nt_reference.fna" "$O/P_Prev_all243_contigs.fna" "$O/P_Prev_all243_nucleotide_local_homology.tsv" "$T/dna" --search-type 3 -s 6.5 --min-seq-id 0.4 -c 0.12 --cov-mode 0 --max-seqs 30000 --format-output 'query,target,fident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits' --threads 4
echo NUC_DONE $(date)
"$M" easy-search "$O/P_Prev_HxlR_all166_proteins.faa" "$O/P_Prev_HxlR_all166_proteins.faa" "$O/P_Prev_HxlR_protein_all_vs_all.tsv" "$T/protein" --search-type 1 -s 7.5 -c 0.5 --cov-mode 0 --max-seqs 1000 --format-output 'query,target,fident,alnlen,qlen,tlen,qcov,tcov,evalue,bits' --threads 4
echo ALL_DONE $(date)
