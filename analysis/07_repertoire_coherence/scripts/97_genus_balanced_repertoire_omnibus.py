#!/usr/bin/env python3
"""Genus-balanced omnibus repertoire convergence. Null ecological-label shuffles within each genus."""
from pathlib import Path
import csv
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13');O=B/'10_additional_analyses_20261012'
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');h=next(rd);rows=list(rd)
features=np.array(h[5:],dtype=object);X=np.array([[float(v) for v in row[5:]] for row in rows],dtype=np.float32);ids=[row[0] for row in rows]
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT').loc[ids]
genus=d.Genus.to_numpy();group_names=sorted(set(genus));idx=[np.flatnonzero(genus==g) for g in group_names]
pidx=[np.flatnonzero((genus==g)&(d.Ecology.to_numpy()=='plaque')) for g in group_names]
nidx=[np.flatnonzero((genus==g)&(d.Ecology.to_numpy()=='mucosal')) for g in group_names]
C=d[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float);C=(C-C.mean(0))/(C.std(0)+1e-12)
Z=np.column_stack([np.ones(len(ids)),pd.get_dummies(d.Genus,drop_first=True,dtype=float),C]);q=np.linalg.qr(Z)[0]
W=np.zeros((len(idx),len(ids)),np.float32)
for j,(a,b) in enumerate(zip(pidx,nidx)):W[j,a]=1./len(a);W[j,b]=-1./len(b)
layers=['PFAM','COG','EC','KMOD','KPATH','X5']
rng=np.random.default_rng(20261013)
nsim=5000
# Permutations shared across layers so BH and max layer correction meaningful.
perm=np.zeros((nsim,len(idx),len(ids)),dtype=np.float32)
for s in range(nsim):
 for j,g in enumerate(idx):
  pp=len(pidx[j]);nn=len(nidx[j]);v=rng.permutation(g)
  perm[s,j,v[:pp]]=1./pp;perm[s,j,v[pp:]]=-1./nn
# Score is off-diagonal directional similarity of mean contrasts, normalized by total effect sum squares.
def calc(DD):
  sumv=DD.sum(axis=-2)
  a=np.sum(sumv*sumv,axis=-1)
  b=np.sum(DD*DD,axis=(-2,-1))
  return (a-b)/(len(idx)-1)/np.maximum(b,1e-15)
result=[]
for layer in layers:
 sel=np.flatnonzero([k.startswith(layer+':') for k in features])
 xx=X[:,sel].astype(np.float64)
 if layer in ('KMOD','KPATH','EC'):pass
 # remove non-variable categories
 v=xx.std(0);xx=xx[:,v>1e-8]
 for transform in ['unscaled','unit_feature_variance']:
  x=xx.copy()
  if transform=='unit_feature_variance':x/=np.maximum(x.std(0),0.1)
  # adjust genome size, GC, CDS, genus without using ecological Y
  x=x-q@(q.T@x)
  obs=calc((W@x)[None,:,:])[0]
  ps=np.empty(nsim)
  for start in range(0,nsim,100):
   k=min(100,nsim-start)
   D=(perm[start:start+k].reshape(k*len(idx),len(ids))@x).reshape(k,len(idx),x.shape[1])
   ps[start:start+k]=calc(D)
  p=(1+np.sum(ps>=obs-1e-12))/(nsim+1)
  result.append(dict(layer=layer,feature_scale=transform,n_features=x.shape[1],n_genera=len(idx),global_direction_similarity=float(obs),perm_p=float(p),
    null_mean=float(ps.mean()),null_95pct=float(np.quantile(ps,.95))))
  print(layer,transform,'features',x.shape[1],'score',round(float(obs),5),'P',round(float(p),5),flush=True)
r=pd.DataFrame(result)
# BH across all 12 layer/scale tests, only descriptive due dependence
def bh(v):
 n=len(v);o=np.argsort(v);q=np.zeros(n);cur=1.
 for j in range(n-1,-1,-1):
  i=o[j];cur=min(cur,v[i]*n/(j+1));q[i]=cur
 return q
r['BH_q_12_views']=bh(r.perm_p.to_numpy())
r.to_csv(O/'genus_balanced_repertoire_omnibus_six_layers.tsv',sep='\t',index=False)
print('FINISHED\n',r.to_string(index=False))
