#!/usr/bin/env python3
"""Extract oriented intergenic regions of divergent PF01638-like regulator–redox-enzyme cassettes.
Status: candidates for shared regulation, NOT demonstrated binding sites/promoters.
Read archived raw genome FASTA and match Prokka GFF contigs by unique length only
when original IDs cannot be matched. Require CDS coordinates and expected gap.
"""
from pathlib import Path
from collections import defaultdict,Counter
import csv,json
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'15_HxlR_regulator_target_ligand_20261014';O.mkdir(exist_ok=True)
R=Path('__WORK_ROOT__/HOMD_convergence_expanded13_work/raw_genomes')
n=pd.read_csv(B/'13_HxlR_broad_and_cassette_20261014/HxlR_allgenome_flanking_3CDS.tsv',sep='\t')
modules={'QorR_like_quinone_reductase':('GCA0028475551_01508','Quinone oxidoreductase 2'),
'CoA_disulfide_reductase':('GCA9380307751_00100','Coenzyme A disulfide reductase'),
'NADH_oxidase':('GCA0378938651_00779','NADH oxidase'),
'MhqN_nitroreductase':('GCA_000339975.1_00196','Putative NAD(P)H nitroreductase MhqN')}
def fasta(path):
 seq={};name=None;parts=[]
 for line in open(path):
  if line.startswith('>'):
   if name:seq[name]=''.join(parts).upper()
   name=line[1:].split()[0];parts=[]
  else:parts.append(line.strip())
 if name:seq[name]=''.join(parts).upper()
 return seq
def rev(x):return x.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def annot(gca):
 d={};lengths={}
 for line in open(B/'02_prokka/gff'/f'{gca}.gff'):
  if line.startswith('##sequence-region'):
   v=line.split()
   if len(v)>=4:lengths[v[1]]=int(v[3])
   continue
  if line.startswith('#'):continue
  cols=line.rstrip('\n').split('\t')
  if len(cols)!=9 or cols[2]!='CDS':continue
  attrs=dict(c.split('=',1) for c in cols[8].split(';') if '=' in c)
  if attrs.get('ID'):d[attrs['ID']]=dict(contig=cols[0],start=int(cols[3]),end=int(cols[4]),strand=cols[6])
 return d,lengths
stats=[];fail=[]
for key,(cluster,prod) in modules.items():
 cand=n[(n.HxlR_X5cluster==cluster)&(n.neighbor_product==prod)&(n.flank_offset_CDS.abs()==1)&(n.intergenic_distance_bp>=40)&(n.intergenic_distance_bp<=500)].copy()
 cand=cand.sort_values(['HMT','GCA']).drop_duplicates(['HMT','GCA'])
 selected=[];records=[]
 for rr in cand.itertuples():
  g=rr.GCA
  if not (R/f'{g}.fna').exists():continue
  d,lens=annot(g)
  a=d.get(rr.HxlR_protein);b=d.get(rr.neighbor_protein)
  if not a or not b or a['contig']!=b['contig']:continue
  left,right=sorted([a,b],key=lambda x:x['start'])
  gap=right['start']-left['end']-1
  if not (40<=gap<=500):continue
  if not (left['strand']=='-' and right['strand']=='+'):continue
  seqs=fasta(R/f'{g}.fna')
  seq=seqs.get(left['contig']) or seqs.get(left['contig'].split('|')[-1])
  how='direct'
  if seq is None:
   choices=[s for s in seqs.values() if len(s)==lens.get(left['contig'],-1)]
   if len(choices)!=1:fail.append((g,'ambiguous_contig_mapping'));continue
   seq=choices[0];how='unique_contig_length'
  nuc=seq[left['end']:right['start']-1]
  if len(nuc)!=gap:fail.append((g,'mismatch_gap'));continue
  oriented=nuc
  # If HxlR is on right, reverse-complement to put regulator consistently on left.
  if right['start']==a['start']:oriented=rev(nuc)
  else:oriented=nuc
  if 'N' in oriented:continue
  selected.append({'key':key,'GCA':g,'HMT':rr.HMT,'Genus':rr.Genus,'protein':rr.HxlR_protein,'target':rr.neighbor_protein,
   'gap_bp':gap,'GFF_contig':left['contig'],'GFF_coordinate_Lend':left['end'],'GFF_coordinate_Rstart':right['start'],'contig_mapping':how,'sequence':oriented})
  if len(selected)%100==0:print('DONE',key,len(selected),flush=True)
 # preserve one representative sequence per HMT; include unique alternative gaps as research duplicates separate
 kept=[];seen=set()
 for row in selected:
  if row['HMT'] in seen:continue
  kept.append(row);seen.add(row['HMT'])
 with open(O/f'{key}_intergenic_1perHMT.fna','w') as fout:
  for row in kept:fout.write('>'+row['GCA']+'|'+row['HMT']+'|'+key+'\n'+row['sequence']+'\n')
 records=kept
 with open(O/f'{key}_intergenic_1perHMT.tsv','w',newline='') as f:
  if kept:
   w=csv.DictWriter(f,fieldnames=list(kept[0]),delimiter='\t');w.writeheader();w.writerows(kept)
 stats.append(dict(module=key,n_genomes_with_candidate=len(cand),n_HMT_candidate=cand.HMT.nunique(),n_genome_DNA_valid=len(selected),n_HMT_recovered=len(kept),
  n_unique_intergenic_sequences=len(set(z['sequence'] for z in kept))))
print('COUNTS',json.dumps(stats,indent=2))
print('EXTRACTION_ERRORS',len(fail),fail[:20])
(O/'intergenic_candidate_modules_summary.json').write_text(json.dumps({'modules':stats,'errors':fail,'warning':'Promoter motifs not inferred; only candidate regulatory regions in divergent loci. Rely on unique contig mapping check and validate before experimental use.'},indent=2))
