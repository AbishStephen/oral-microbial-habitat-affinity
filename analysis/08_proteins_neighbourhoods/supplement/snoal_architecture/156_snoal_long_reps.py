#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'18_SnoaL2_function_gene_integrity_20261016'
S=Path('__WORK_ROOT__/SnoaL2_subfamilies_20261016')
aa=pd.read_csv(O/'SnoaL2_original524_protein_sizes.tsv',sep='\t')
c=pd.read_csv(B/'17_SnoaL2_independent_targets_20261016/SnoaL2_all524_protein_clusters.tsv',sep='\t')
aa=aa.merge(c[['protein_id','cluster']],on='protein_id',how='left',validate='one_to_one')
seqs={};n=None;v=[]
for line in open(S/'SnoaL2_all524_proteins.faa'):
 if line.startswith('>'):
  if n:seqs[n]=''.join(v)
  n=line[1:].strip();v=[]
 else:v.append(line.strip())
if n:seqs[n]=''.join(v)
long=aa[aa.protein_length_aa>=200]
print('LONG>=200',len(long),'unique_seq',long.protein_id.map(lambda x:seqs[x]).nunique(),'HMT',long.HMT.nunique(),'clusters',long.cluster.nunique(),flush=True)
print('LENGTH_COUNTS',long.protein_length_aa.value_counts().head(25).to_dict(),flush=True)
reps=long.drop_duplicates('cluster').sort_values('protein_length_aa',ascending=False)
reps.to_csv(O/'SnoaL2_long_cluster_representatives.tsv',sep='\t',index=False)
with open(S/'SnoaL2_long_cluster_reps.faa','w') as f:
 for r in reps.itertuples():f.write('>'+r.protein_id+'\n'+seqs[r.protein_id]+'\n')
print('TOP_CLUSTERS\n'+long.groupby('cluster').agg(n=('protein_id','size'),median_length=('protein_length_aa','median'),max_length=('protein_length_aa','max'),genera=('Genus','nunique'),n_HMT=('HMT','nunique')).sort_values('n',ascending=False).head(20).to_string(),flush=True)
print('REPS\n'+reps[['protein_id','Genus','HMT','protein_length_aa','cluster','linked_HxlR_3kb']].to_string(index=False),flush=True)
