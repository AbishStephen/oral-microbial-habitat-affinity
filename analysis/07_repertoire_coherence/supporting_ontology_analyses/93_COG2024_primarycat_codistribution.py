#!/usr/bin/env python3
"""Conditional co-distribution of COG2024 functional classes; genus-block permutations."""
from pathlib import Path
from scipy.stats import pearsonr
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
S=pd.read_csv(O/'COG2024_primarycat_scores_HMT47.tsv',sep='\t').set_index('HMT')
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT').loc[S.index]
labs=list(S.columns);raw=S.to_numpy(float); n,m=raw.shape
catname={}
for line in (O/'ontology/cog-24.fun.tab').read_text().splitlines():
 z=line.split('\t')
 if len(z)>=4 and len(z[0])==1:catname[z[0]]=z[3]
g=d.Genus.to_numpy();groups=[np.flatnonzero(g==x) for x in sorted(set(g))]
C=d[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float)
C=(C-C.mean(0))/(C.std(0)+1e-12)
G=pd.get_dummies(d.Genus,drop_first=True,dtype=float).to_numpy()
Z=np.column_stack([np.ones(n),G,C]);q=np.linalg.qr(Z)[0]
project=lambda x: x-q@(q.T@x)
A=project(raw);A=A/(np.linalg.norm(A,axis=0)+1e-15)
YY=project(d.Y.to_numpy(float))
Yscore=YY/(np.linalg.norm(YY)+1e-15)
Aa=project(np.column_stack([raw,d.Y.to_numpy(float)]))
# residualize taxonomy, genome size, GC, CDS, plus ecology Y
qq=np.linalg.qr(np.column_stack([Z,d.Y.to_numpy(float)]))[0]
Bmat=raw-qq@(qq.T@raw);Bmat=Bmat/(np.linalg.norm(Bmat,axis=0)+1e-15)
pairs=[(i,j) for i in range(m) for j in range(i+1,m)]
obs=A.T@A
ob2=Bmat.T@Bmat
rawr=np.corrcoef(raw,rowvar=False)
# Randomize one category within genus and condition back on covariates (marginal test).
rng=np.random.default_rng(20261012);iters=5000
pcounts=np.zeros((m,m),dtype=int);pycounts=np.zeros((m,m),dtype=int)
for j in range(m):
 for rep in range(iters):
  x=A[:,j].copy()
  for ix in groups:x[ix]=rng.permutation(x[ix])
  x=project(x);x/=np.linalg.norm(x)+1e-15
  r=x@A
  pcounts[:j,j]+=np.abs(r[:j])>=np.abs(obs[:j,j])-1e-12
  # second residual: shuffled profile projected again on Y, preserving phenotype-conditioned null
  xy=x-qq@(qq.T@x);xy/=np.linalg.norm(xy)+1e-15
  ry=xy@Bmat
  pycounts[:j,j]+=np.abs(ry[:j])>=np.abs(ob2[:j,j])-1e-12
print('NULLS_DONE',iters,'per category',flush=True)
rows=[]
for i,j in pairs:
 x=d.drop(columns=[],errors='ignore')
 loo=[]
 for gg in groups:
  keep=np.setdiff1d(np.arange(n),gg)
  if len(keep)<5:continue
  vv=np.corrcoef(A[keep,i],A[keep,j])[0,1]
  loo.append(vv)
 rows.append(dict(class1=labs[i],class2=labs[j],name1=catname[labs[i]],name2=catname[labs[j]],
  raw_correlation=float(rawr[i,j]),adjusted_correlation=float(obs[i,j]),
  adjusted_ecology_conditioned_correlation=float(ob2[i,j]),
  p_blockperm=(pcounts[i,j]+1)/(iters+1),
  p_conditioned_blockperm=(pycounts[i,j]+1)/(iters+1),
  leave_one_genus_out_r_min=float(min(loo)),leave_one_genus_out_r_max=float(max(loo)),
  same_sign_13=int(sum((x>0)==(obs[i,j]>0) for x in loo))))
df=pd.DataFrame(rows)
def bh(v):
 z=np.asarray(v);n=len(z);order=np.argsort(z);out=np.empty(n)
 cur=1.
 for k in range(n-1,-1,-1):
  index=order[k];cur=min(cur,z[index]*n/(k+1));out[index]=cur
 return out
df['q_blockperm_BH']=bh(df.p_blockperm);df['q_conditioned_BH']=bh(df.p_conditioned_blockperm)
df=df.sort_values('p_conditioned_blockperm')
df.to_csv(O/'COG2024_primarycat_codistribution_genus_blockperm.tsv',sep='\t',index=False)
print('EDGES_TESTED',len(df),'PERSIST_ECOLOGY_Q05',sum(df.q_conditioned_BH<.05),'PERSIST_ECOLOGY_Q10',sum(df.q_conditioned_BH<.1),'NORMAL_Q05',sum(df.q_blockperm_BH<.05),flush=True)
print('STRONGEST ECOLOGY CONDITIONED\n',df[['class1','class2','name1','name2','adjusted_correlation','adjusted_ecology_conditioned_correlation','p_conditioned_blockperm','q_conditioned_BH','same_sign_13']].head(25).to_string(index=False),flush=True)
