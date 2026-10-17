#!/usr/bin/env python3
"""Physical linkage and protein architecture of cross-genus candidate Pfams in 1121 HOMD genomes.
Pairwise closest same-contig distinct CDS, 3 kb/10 kb; counts of independent genus/HMT occurrences.
"""
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter
import gzip,csv,json,os
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012';O.mkdir(exist_ok=True)
pfams=['PF01638','PF04072','PF13577','PF13676','PF22166','PF01745','PF01337','PF13413','PF18075','PF14470','PF12680','PF10604','PF07366']
pfset=set(pfams)
manifest=pd.read_csv(B/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')
lookup={r.GCA:(r.HMT,r.Genus,r.Species) for r in manifest.itertuples()}
P=defaultdict(set);names={};seen_pergenome=defaultdict(set)
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  pf=r['Pfam']
  if pf in pfset and r['GCA'] in lookup:
   P[r['GCA']].add((r['protein_id'],pf))
   names[pf]=r['Pfam_name']
   seen_pergenome[r['GCA']].add(pf)
print('HMM_LOADED','genes',sum(len(x) for x in P.values()),'genomes',len(P),'Pfam_names',names,flush=True)
# coordinates from all Prokka GFF using prefixed composite headers (the contig includes accession prefix)
loc={};all_genes=defaultdict(list);meta={}
target_set={gene for pairs in P.values() for gene,pf in pairs}
for i,(gca,ss) in enumerate(P.items()):
 gf=B/'02_prokka/gff'/f'{gca}.gff'
 if not gf.exists():continue
 with gf.open(errors='replace') as f:
  for line in f:
   if line.startswith('#'):continue
   z=line.rstrip().split('\t')
   if len(z)!=9 or z[2]!='CDS':continue
   attr=z[8]; pre='ID='
   k=attr.find(pre)
   if k<0:continue
   eid=attr[k+len(pre):].split(';')[0]
   st=int(z[3]);en=int(z[4]);contig=z[0];strand=z[6]
   all_genes[gca].append((contig,st,en,strand,eid))
   if eid in target_set:
    product=''
    for t in attr.split(';'):
     if t.startswith('product='):product=t[8:];break
    loc[eid]=(contig,st,en,strand,gca,product)
 if i%250==0:print('GFFS_PROCESSED',i,flush=True)
print('COORD_MATCH',len(loc),'of',len(target_set),'total_GFF_genes',sum(len(v) for v in all_genes.values()),flush=True)
# candidate genes positions, proteins carrying domains, other domains in proteins
gene_domains=defaultdict(set)
for gca,rows in P.items():
 for pr,pf in rows:gene_domains[pr].add(pf)
gene_tsv=[]
for prot,doms in gene_domains.items():
 if prot not in loc:continue
 contig,st,en,strand,gca,prod=loc[prot];h,g,s=lookup[gca]
 for pf in doms:
  gene_tsv.append(dict(GCA=gca,HMT=h,Genus=g,Species=s,protein_id=prot,contig=contig,start=st,end=en,strand=strand,Pfam=pf,Pfam_name=names[pf],Prokka_product=prod,domains_on_protein=';'.join(sorted(doms))))
pd.DataFrame(gene_tsv).to_csv(O/'crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t',index=False)
# all 13 PFAM joint pairs
pair_rows=[];near_rows=[]
for a,b in combinations(pfams,2):
 cg=[gca for gca in P if a in seen_pergenome[gca] and b in seen_pergenome[gca]]
 same=0;close3=0;close10=0;samecontig=0;indep3=0
 hitgenus3=set();hitgenus10=set();hitgca3=set();hithmt3=set();close3genes=set();examples=[]
 for gca in cg:
  aset=[(pr,loc[pr]) for pr,dom in P[gca] if dom==a and pr in loc]
  bset=[(pr,loc[pr]) for pr,dom in P[gca] if dom==b and pr in loc]
  if not aset or not bset:continue
  sameprotein=any(x==y for x,_ in aset for y,_ in bset)
  if sameprotein:same+=1
  nearest=None
  for x,Lx in aset:
   for y,Ly in bset:
    if x==y or Lx[0]!=Ly[0]:continue
    dist=max(0,max(Lx[1],Ly[1])-min(Lx[2],Ly[2])-1)
    if nearest is None or dist<nearest[0]:nearest=(dist,x,y,Lx,Ly)
  if nearest:
   samecontig+=1;delta,pr1,pr2,Lx,Ly=nearest
   if delta<=10000:
    close10+=1;hitgenus10.add(lookup[gca][1])
   if delta<=3000:
    close3+=1;indep3+=1;hithmt3.add(lookup[gca][0]);hitgenus3.add(lookup[gca][1]);hitgca3.add(gca)
    near_rows.append(dict(Pfam1=a,Pfam2=b,GCA=gca,HMT=lookup[gca][0],Genus=lookup[gca][1],protein1=pr1,protein2=pr2,distance_bp=delta,
    strand1=Lx[3],strand2=Ly[3],same_strand=Lx[3]==Ly[3],product1=Lx[5],product2=Ly[5]))
    if len(examples)<4:examples.append(f'{lookup[gca][1]}:{lookup[gca][0]}:{delta}bp')
 # suppress tiny cooccur pairs later
 pair_rows.append(dict(Pfam1=a,Pfam2=b,name1=names[a],name2=names[b],
  n_genomes_both=len(cg),n_genomes_same_protein=same,n_genomes_distinct_genes_same_contig=samecontig,
  n_genomes_distinct_genes_nearest_le3kb=close3,n_genomes_distinct_genes_nearest_le10kb=close10,
  n_genera_distinct_genes_le3kb=len(hitgenus3),n_HMT_distinct_genes_le3kb=len(hithmt3),
  n_genera_distinct_genes_le10kb=len(hitgenus10),
  linked_genus_list=';'.join(sorted(hitgenus3)),
  examples=';'.join(examples)))
pfpairs=pd.DataFrame(pair_rows)
pfpairs.sort_values(['n_genera_distinct_genes_le3kb','n_genomes_distinct_genes_nearest_le3kb','n_genomes_same_protein'],ascending=False).to_csv(O/'crossgenus_all_candidate_Pfam_pairwise_physical_linkage.tsv',sep='\t',index=False)
pd.DataFrame(near_rows).to_csv(O/'crossgenus_candidate_distinct_genes_within3kb_examples.tsv',sep='\t',index=False)
print('COORDINATES',len(gene_tsv),'LINK_PAIRS_ANY_3KB',(pfpairs.n_genomes_distinct_genes_nearest_le3kb>0).sum(),flush=True)
print('LINK_PAIRS_TOP\n',pfpairs[['Pfam1','Pfam2','n_genomes_both','n_genomes_same_protein','n_genomes_distinct_genes_nearest_le3kb','n_genera_distinct_genes_le3kb','n_HMT_distinct_genes_le3kb','linked_genus_list']].head(35).to_string(index=False),flush=True)
for pf in pfams:
 vv=[x for x in gene_tsv if x['Pfam']==pf]
 products=Counter(x['Prokka_product'] for x in vv)
 grp=Counter(x['Genus'] for x in vv)
 print('ANNOT',pf,names[pf],'n_proteins',len(set(x['protein_id'] for x in vv)),'n_genomes',len(set(x['GCA'] for x in vv)),
       'n_genera',len(grp),'PRODUCTS',products.most_common(5),flush=True)
