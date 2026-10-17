#!/usr/bin/env python3
"""Full nested covariate-residualized model refitted under within-genus Y and label permutation.
"""
import sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
src=Path('__PROJECT_ROOT__/study13/scripts/111_foldwise_covariate_adjusted_prediction.py')
ss=src.read_text().splitlines();cut=next(i for i,z in enumerate(ss) if z.startswith("for label,y in "));exec(compile("\n".join(ss[:cut])+"\n",str(src),"exec"))
seed=int(sys.argv[1]);rng=np.random.default_rng(20261013+seed);reindex=np.arange(len(hmt))
for g in gs:
 idx=np.where(gen==g)[0];reindex[idx]=rng.permutation(idx)
eco=eco[reindex]
out=[]
for label,y in [('original',yp[reindex]),('sourcebalanced',y_alt[reindex])]:
 pred=np.full(len(y),np.nan)
 for og in gs:
  tr=np.where(gen!=og)[0];te=np.where(gen==og)[0]
  inn={}
  for ig in np.unique(gen[tr]):
   it=tr[gen[tr]==ig];ix=tr[gen[tr]!=ig];zz=adjusted_X(ix);inn[ig]=(ix,it,zz,rank_adjusted(zz,ix))
  best=None
  for k in Ks:
   for a in alphas:
    pp=np.full(len(y),np.nan)
    for ig,(ix,it,zz,rr) in inn.items():pp[it]=pred_adjusted(zz,ix,it,rr,k,a,y)
    score=auc_mean(tr,pp[tr]);rho=float(spearmanr(y[tr],pp[tr]).statistic)
    if not np.isfinite(rho):rho=-1
    loss=float(np.mean(abs(y[tr]-pp[tr])));key=(score,rho,-loss,-k,-a)
    if best is None or key>best[0]:best=(key,k,a)
  k,a=best[1:];zz=adjusted_X(tr);rr=rank_adjusted(zz,tr);pred[te]=pred_adjusted(zz,tr,te,rr,k,a,y)
 out.append(dict(perm_id=seed,label=label,macro_auc=auc_mean(np.arange(len(y)),pred),rho=float(spearmanr(y,pred).statistic)))
pd.DataFrame(out).to_csv(O/f'covadjusted_nested_perm_{seed:04d}.tsv',sep='\t',index=False)
print(out,flush=True)
