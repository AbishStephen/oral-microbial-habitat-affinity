#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
import pandas as pd
O=Path('__PROJECT_ROOT__/study13/10_additional_analyses_20261012')
near=pd.read_csv(O/'crossgenus_candidate_distinct_genes_within3kb_examples.tsv',sep='\t')
a=pd.read_csv(O/'crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t').drop_duplicates(['GCA','protein_id']).set_index(['GCA','protein_id'])
out=[]
for rr in near.itertuples():
 u=a.loc[(rr.GCA,rr.protein1)];v=a.loc[(rr.GCA,rr.protein2)]
 if u.start<v.start: left,right=u,v;ln='A'
 else:left,right=v,u;ln='B'
 orientation='divergent' if left.strand=='-' and right.strand=='+' else ('convergent' if left.strand=='+' and right.strand=='-' else 'same-strand')
 out.append(dict(Pfam1=rr.Pfam1,Pfam2=rr.Pfam2,GCA=rr.GCA,HMT=rr.HMT,Genus=rr.Genus,
  orientation=orientation,order=ln,distance_bp=rr.distance_bp,
  product1=rr.product1,product2=rr.product2))
s=pd.DataFrame(out)
s.to_csv(O/'crossgenus_candidate_neighborhood_orientation.tsv',sep='\t',index=False)
for aa,bb in [('PF01638','PF12680'),('PF01638','PF22166'),('PF12680','PF10604')]:
 v=s[(s.Pfam1==aa)&(s.Pfam2==bb)]
 print('\n',aa,bb,'N',len(v),'ORIENTATION',v.orientation.value_counts().to_dict(),'GENERA',v.Genus.value_counts().to_dict())
 print(v.groupby(['Genus','orientation']).agg(n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),median_bp=('distance_bp','median')).to_string())
