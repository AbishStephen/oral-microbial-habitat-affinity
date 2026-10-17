#!/usr/bin/env python3
from pathlib import Path
import csv
B=Path('__PROJECT_ROOT__/study13')
O=B/'12_HxlR_SnoaL_regulatory_genomics_20261013'
S=Path('__WORK_ROOT__/HOMD_convergence_expanded13_work/raw_genomes')
rows={r['GCA']:r for r in csv.DictReader(open(O/'linked_locus_audit.tsv'),delimiter='\t')}
def fasta(f):
 d={};n=None;buf=[]
 for l in open(f):
  if l.startswith('>'):
   if n:d[n]=''.join(buf)
   n=l[1:].split()[0];buf=[]
  elif n:buf.append(l.strip().upper())
 if n:d[n]=''.join(buf)
 return d
def comp(x):return x.translate(str.maketrans('ACGT','TGCA'))[::-1]
def get(g):
 r=rows[g];seqs=fasta(S/f'{g}.fna' if (S/f'{g}.fna').exists() else O/'genomes'/f'{g}.fna')
 return seqs[r['accession_contig']],int(r['start_1based']),int(r['end_1based'])
def paired(a,b):
 aseq,astart,aend=get(a);bseq,bstart,bend=get(b)
 seed=aseq[astart-1:aend]
 bidx=bseq.find(seed);direction='+'
 if bidx<0:
  bseq=comp(bseq);bidx=bseq.find(seed);direction='-'
 if bidx<0:return dict(A=a,B=b,matched_intergenic=False)
 ai=astart-1; bi=bidx
 if bseq.find(seed,bidx+1)>=0:print('WARN multiple seed',b)
 l=0
 while ai-l-1>=0 and bi-l-1>=0 and aseq[ai-l-1]==bseq[bi-l-1]:l+=1
 r=0
 while ai+len(seed)+r<len(aseq) and bi+len(seed)+r<len(bseq) and aseq[ai+len(seed)+r]==bseq[bi+len(seed)+r]:r+=1
 return dict(A=a,B=b,matched_intergenic=True,relative_orientation=direction,identical_nt_span=len(seed)+l+r,
             nt_identical_to_left=l,nt_identical_to_right=r,intergenic_bp=len(seed),
             complete_5kb_left=l>=5000,complete_5kb_right=r>=5000)
pairs=[
('GCA_002892595.1','GCA_001548195.1'),
('GCA_002892595.1','GCA_002753835.1'),
('GCA_000271945.1','GCA_000261025.1'),
('GCA_001263815.1','GCA_001548195.1')
]
result=[paired(a,b) for a,b in pairs]
import json
(O/'cross_genus_identical_tracts.json').write_text(json.dumps(result,indent=2))
for z in result:print(json.dumps(z))
