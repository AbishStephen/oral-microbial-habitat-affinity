#!/bin/bash
#SBATCH --job-name=HOMD13_TREE
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=__PROJECT_ROOT__/study13/logs/phylogeny_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/phylogeny_%j.err
set -euo pipefail
BASE=__PROJECT_ROOT__/study13
OUT=$BASE/07_phylogeny
MAFFT=mafft
RAXML=raxml-ng
python3 $BASE/scripts/55_build_phylogeny.py
rm -rf $OUT/aligned
mkdir -p $OUT/aligned
while read g; do $MAFFT --auto --thread 1 $OUT/$g.faa > $OUT/aligned/$g.aln.faa 2>>$OUT/mafft.log; done < $OUT/marker_genes.txt
python3 - <<'PY2'
from pathlib import Path
B=Path("__PROJECT_ROOT__/study13/07_phylogeny")
genes=[x.strip() for x in open(B/"marker_genes.txt") if x.strip()]
seq={}; order=[]
for i,g in enumerate(genes):
    d={}; k=None; buf=[]
    for line in open(B/"aligned"/f"{g}.aln.faa"):
        if line.startswith(">"):
            if k:d[k]="".join(buf)
            k=line[1:].strip();buf=[]
        else:buf.append(line.strip())
    if k:d[k]="".join(buf)
    if i==0:order=sorted(d)
    assert set(d)==set(order),(g,len(d),len(order))
    for x in order:seq[x]=seq.get(x,"")+d[x]
with open(B/"HOMD47_ribosomal_concat.faa","w") as f:
    for x in order:f.write(f">{x}\n{seq[x]}\n")
print("Taxa",len(order),"Alignment aa",len(next(iter(seq.values()))))
assert len(order)==47
PY2
cd $OUT
rm -f HOMD47.raxml.*
$RAXML --search --msa HOMD47_ribosomal_concat.faa --model LG+G --threads 8 --seed 20261011 --prefix HOMD47
cp HOMD47.raxml.bestTree $BASE/handoff/Z_HOMD47_ribosomal_tree.newick
echo TREE
cat HOMD47.raxml.bestTree
