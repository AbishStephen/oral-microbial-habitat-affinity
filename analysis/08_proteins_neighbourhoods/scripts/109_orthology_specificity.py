#!/usr/bin/env python3
"""Disaggregate candidate Pfam families into MMseqs-X5 clusters and map locus synteny.
Descriptive HMT/genus effects; do not infer HGT without gene tree or nucleotide evidence.
"""
from pathlib import Path
from collections import defaultdict
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13');O=B/'11_manuscript_completion_20261013'
dom=pd.read_csv(O/'candidate_Pfam_fullprotein_orthogroup_architecture.tsv',sep='\t')
dic=pd.read_csv(B/'06_hmt_matrices/X5_cluster_dictionary.tsv',sep='\t')
dic=dic.set_index('representative_protein')
dt=dom.merge(dic[['cluster_id']],left_on='MMseqs_X5_representative',right_index=True,how='left')
needed=set(dt.cluster_id.dropna()); p=pd.read_csv(B/'06_hmt_matrices/X5_HMT_prevalence_ALL.tsv',sep='\t',usecols=lambda c: c in {'HMT','Genus','Ecology'} or c in needed)
print('MATRIX',p.shape,'COLS',p.columns[:8].tolist(),flush=True)
m=pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t')[['HMT','Genus','Ecology','Ecology_tier']]
q=p
s=dt.groupby(['Pfam','cluster_id']).agg(n_proteins=('protein_id','nunique'),n_genera=('Genus','nunique'),n_HMT=('HMT','nunique')).reset_index()
rows=[]
for r in s.itertuples():
 if r.n_genera<2 or r.cluster_id not in q.columns:continue
 d=[]
 for g,grp in q.groupby('Genus'):
  pp=grp[grp.Ecology=='plaque'][r.cluster_id];mm=grp[grp.Ecology=='mucosal'][r.cluster_id]
  if len(pp) and len(mm):
   d.append((g,float(pp.mean()-mm.mean())))
 pll=[g for g,x in d if x>.05]; muc=[g for g,x in d if x<-.05]
 rows.append(dict(Pfam=r.Pfam,MMseqs_cluster=r.cluster_id,cluster_representative=dic.index[dic.cluster_id==r.cluster_id][0],
  n_genomes=int(dt[(dt.Pfam==r.Pfam)&(dt.cluster_id==r.cluster_id)].GCA.nunique()),
  n_HMT=r.n_HMT,n_genera=r.n_genera,mean_equal_genus_effect=np.mean([x for g,x in d]),
  plaque_genus_support=len(pll),mucosal_genus_support=len(muc),
  plaque_genus_names=';'.join(pll),mucosal_genus_names=';'.join(muc)))
R=pd.DataFrame(rows).sort_values(['Pfam','n_genera','n_HMT'],ascending=[True,False,False])
R.to_csv(O/'candidate_Pfam_crossgenus_orthogroup_effects.tsv',sep='\t',index=False)
print('X5_AT_LEAST_2_GENERA',len(R),'TOP\n',R.sort_values(['n_genera','plaque_genus_support'],ascending=False).head(20).to_string(index=False),flush=True)
near=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_candidate_distinct_genes_within3kb_examples.tsv',sep='\t')
v=near[(near.Pfam1=='PF01638')&(near.Pfam2=='PF12680')].copy()
lookup=dt.drop_duplicates(['protein_id','GCA']).set_index(['GCA','protein_id']).MMseqs_X5_representative
v['HxlR_X5_representative']=[lookup.get((x.GCA,x.protein1),'-') for x in v.itertuples()]
v['SnoaL2_X5_representative']=[lookup.get((x.GCA,x.protein2),'-') for x in v.itertuples()]
v.to_csv(O/'HxlR_SnoaL2_linked_gene_orthogroup_pairs.tsv',sep='\t',index=False)
pair=(v.groupby(['HxlR_X5_representative','SnoaL2_X5_representative'])
      .agg(n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),n_genera=('Genus','nunique'),genera=('Genus',lambda x:';'.join(sorted(set(x)))))
      .reset_index().sort_values('n_genomes',ascending=False))
pair.to_csv(O/'HxlR_SnoaL2_linked_orthogroup_pair_summary.tsv',sep='\t',index=False)
print('LOCUS_ORTHOGROUP_PAIR_COUNT',len(pair),'MULTIGENUS',int(sum(pair.n_genera>1)))
print(pair.head(15).to_string(index=False))
