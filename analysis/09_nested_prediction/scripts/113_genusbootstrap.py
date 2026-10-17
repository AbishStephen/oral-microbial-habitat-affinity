#!/usr/bin/env python3
"""Genus-cluster bootstrap CIs of macro pairwise AUC and paired model differences.
Not an independent validation CI; 13 genus blocks selected for ecological extremes.
"""
from pathlib import Path
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13');O=B/'11_manuscript_completion_20261013'
r={}
for name,path in {
 'original_pooled_rank':B/'08_validation/SPARSE_DEDUP_outer_predictions.tsv',
 'original_genusbalanced':O/'nested_genusbalanced_original_HMT_predictions.tsv',
 'original_genusbalanced_covadj':O/'nested_genusbalanced_covariate_residualized_original.tsv',
 'sourcebalanced_pooled_rank':B/'08_validation/GROUPBAL_SPARSE_DEDUP_outer_predictions.tsv',
 'sourcebalanced_genusbalanced':O/'nested_genusbalanced_sourcebalanced_HMT_predictions.tsv',
 'sourcebalanced_genusbalanced_covadj':O/'nested_genusbalanced_covariate_residualized_sourcebalanced.tsv'
}.items():
 d=pd.read_csv(path,sep='\t'); arr=[]
 for g,z in d.groupby('Genus'):
  plus=z[z.Ecology=='plaque'].Y_pred.values;minus=z[z.Ecology=='mucosal'].Y_pred.values
  auc=np.mean([1. if a>b else .5 if a==b else 0. for a in plus for b in minus])
  arr.append((g,auc))
 r[name]=dict(arr)
assert all(set(x)==set(next(iter(r.values()))) for x in r.values())
g=sorted(next(iter(r.values())))
A=pd.DataFrame({k:[v[x] for x in g] for k,v in r.items()},index=g)
A.index.name='Genus';A.to_csv(O/'all_predictive_models_per_genus_AUC.tsv',sep='\t')
rng=np.random.default_rng(20261013)
ndraw=50000;ind=rng.integers(len(g),size=(ndraw,len(g)))
d=[]
for key in A.columns:
 v=A[key].values;vv=v[ind].mean(axis=1)
 d.append(dict(Model=key,macro_auc=v.mean(),bootstrap_genus_lo95=np.quantile(vv,.025),bootstrap_genus_hi95=np.quantile(vv,.975)))
for a,b in [('original_pooled_rank','original_genusbalanced'),('original_genusbalanced','original_genusbalanced_covadj'),('sourcebalanced_pooled_rank','sourcebalanced_genusbalanced'),('sourcebalanced_genusbalanced','sourcebalanced_genusbalanced_covadj')]:
 z=(A[a].values-A[b].values)[ind].mean(axis=1)
 d.append(dict(Model='Difference '+a+' minus '+b,macro_auc=np.mean(A[a]-A[b]),bootstrap_genus_lo95=np.quantile(z,.025),bootstrap_genus_hi95=np.quantile(z,.975)))
S=pd.DataFrame(d)
S.to_csv(O/'predictive_macroAUC_genus_bootstrap95.tsv',sep='\t',index=False)
print(S.to_string(index=False),flush=True)
