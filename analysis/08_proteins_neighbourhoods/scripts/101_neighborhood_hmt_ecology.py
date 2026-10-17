#!/usr/bin/env python3
"""Within-HMT gene-neighborhood carriage, ecological association, HxlR-SnoaL architecture overlap."""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import numpy as np,pandas as pd,csv,gzip
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
pairs=pd.read_csv(O/'crossgenus_candidate_distinct_genes_within3kb_examples.tsv',sep='\t')
genomes=pd.read_csv(B/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')
hmt=pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t')
G=genomes.groupby(['HMT','Genus','Ecology']).GCA.nunique().rename('genome_total').reset_index()
focus=[('PF01638','PF12680'),('PF01638','PF07366'),('PF12680','PF10604'),('PF01638','PF04072'),('PF01638','PF22166')]
out=[]
for aa,bb in focus:
 sub=pairs[(pairs.Pfam1==aa)&(pairs.Pfam2==bb)]
 D=sub.groupby('HMT').GCA.nunique().rename('genomes_linked')
 t=G.join(D,on='HMT');t['genomes_linked']=t.genomes_linked.fillna(0).astype(int)
 t['fraction_linked']=t.genomes_linked/t.genome_total
 t['feature_pair']=aa+'|'+bb
 out.append(t)
 print('\nPAIR',aa,bb,'nGENOMES',sub.GCA.nunique(),'nHMT',sub.HMT.nunique(),'nGENUS',sub.Genus.nunique(),flush=True)
 print(t[t.genomes_linked>0][['HMT','Genus','Ecology','genomes_linked','genome_total','fraction_linked']].to_string(index=False),flush=True)
 if len(sub):
  print('DISTANCES','<=500bp',int((sub.distance_bp<=500).sum()),'<=1kb',int((sub.distance_bp<=1000).sum()),'median',sub.distance_bp.median(),'same orientation',float(sub.same_strand.mean()),flush=True)
r=pd.concat(out);r.to_csv(O/'crossgenus_candidate_neighborhood_HMT_ecology_prevalence.tsv',sep='\t',index=False)
# Exact genus-balanced ecological association for linked HxlR + SnoaL_2; all genus HMT included.
target=r[r.feature_pair=='PF01638|PF12680'].copy()
bygroup=list(target.groupby('Genus'))
diff=[]
for genus,d in bygroup:
 p=d[d.Ecology=='plaque'].fraction_linked; m=d[d.Ecology=='mucosal'].fraction_linked
 diff.append(dict(genus=genus,plaque_mean=p.mean(),mucosal_mean=m.mean(),difference=p.mean()-m.mean()))
out=pd.DataFrame(diff)
out.to_csv(O/'HxlR_SnoaL2_linkage_bygenus_effect.tsv',sep='\t',index=False)
print('\nHxlR-SnoaL2 GENUS ECOLOGY EFFECT\n',out.to_string(index=False),flush=True)
print('MEAN ALL GENERA',out.difference.mean(),'positive',sum(out.difference>0),'negative',sum(out.difference<0),'equal',sum(out.difference==0),flush=True)
# estimated exact permutation P via within genus shuffling, 40000 perms
rg=np.random.default_rng(20261012);n=50000
Y=target.fraction_linked.to_numpy()
groups=[target[target.Genus==g] for g in sorted(target.Genus.unique())]
null=np.zeros(n)
for group in groups:
 a=group.fraction_linked.to_numpy()
 npos=int(sum(group.Ecology=='plaque'));nn=len(a)
 pp=np.zeros(n)
 for k in range(n):
  perm=rg.permutation(nn)
  pp[k]=a[perm[:npos]].mean()-a[perm[npos:]].mean()
 null+=pp/len(groups)
obs=float(out.difference.mean())
print('HXL_SNOA_L2_GENUS_BALANCED_PERM',n,'obs',obs,'P_ONE_DIRECTION',(1+int((null>=obs-1e-12).sum()))/(n+1),'P_TWO',(1+int((abs(null)>=abs(obs)-1e-12).sum()))/(n+1),flush=True)
# inspect PF12680/PF07366 protein-level domain alignment coordinates, not protein distance
found=defaultdict(dict)
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for d in csv.DictReader(f,delimiter='\t'):
  if d['Pfam'] in {'PF12680','PF07366'}:
   pr=d['protein_id'];pf=d['Pfam'];lo=int(d['ali_from']);hi=int(d['ali_to'])
   if pf not in found[pr] or (hi-lo)>(found[pr][pf][1]-found[pr][pf][0]):found[pr][pf]=(lo,hi)
ov=[]
for pr,v in found.items():
 if len(v)==2:
  (a,b),(c,d)=v['PF12680'],v['PF07366']
  over=max(0,min(b,d)-max(a,c)+1);union=max(b,d)-min(a,c)+1
  ov.append((pr,a,b,c,d,over,union,over/union))
if ov:
 t=pd.DataFrame(ov,columns=['protein_id','snoal2_start','snoal2_end','snoal_start','snoal_end','overlap_aa','spanning_region_length_aa','frac_overlap_union'])
 t.to_csv(O/'PF12680_PF07366_sameprotein_domain_overlap.tsv',sep='\t',index=False)
 print('SAMEPROTEIN_SnoaL2_SnoaL',len(t),'median overlapAA',t.overlap_aa.median(),'fraction >50% overlapping region',float((t.frac_overlap_union>.5).mean()),flush=True)
