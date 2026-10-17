#!/usr/bin/env python3
"""Full nested genus-balanced refitting under within-genus phenotype permutation.
Each Slurm array task is one permuted ecological dataset, both primary outcomes.
"""
from pathlib import Path
import sys,os,time
import numpy as np,pandas as pd
src=Path('__PROJECT_ROOT__/study13/scripts/104_nested_genusbalanced.py')
exec(compile(src.read_text().split("for label,y in [('original',yp),('sourcebalanced',y_alt)]:")[0],str(src),'exec'))
seed=int(sys.argv[1]); rg=np.random.default_rng(20261013+seed); ix_all=np.arange(len(hmt))
for g in gs:
 ii=np.where(gen==g)[0];ix_all[ii]=rg.permutation(ii)
# ecological preference and abundance scores shuffled JOINTLY among HMTs within genus
eco=eco[ix_all]; targets=[('original',yp[ix_all]),('sourcebalanced',y_alt[ix_all])]
out=[]
for label,y in targets:
 preds=np.full(len(y),np.nan)
 for outer in gs:
  tr=np.where(gen!=outer)[0];te=np.where(gen==outer)[0]
  inner={}
  for inn in np.unique(gen[tr]):
   it=tr[gen[tr]==inn];ix=tr[gen[tr]!=inn];inner[inn]=(ix,it,rank(ix))
  best=None
  for k in Ks:
   for a in alphas:
    ixpred=np.full(len(y),np.nan)
    for inn,(ix,it,rr) in inner.items():ixpred[it]=fitpredict(ix,it,rr,k,a,y)
    score=auc_mean(tr,ixpred[tr])
    rho=float(spearmanr(y[tr],ixpred[tr]).statistic)
    if not np.isfinite(rho):rho=-1
    err=float(np.mean(abs(y[tr]-ixpred[tr])))
    key=(score,rho,-err,-k,-a)
    if best is None or key>best[0]:best=(key,k,a)
  k,a=best[1:]
  rr=rank(tr)
  preds[te]=fitpredict(tr,te,rr,k,a,y)
 out.append(dict(perm_id=seed,label=label,macro_auc=auc_mean(np.arange(len(y)),preds),pooled_rho=float(spearmanr(y,preds).statistic)))
pd.DataFrame(out).to_csv(O/f'nested_perm_{seed:04d}.tsv',sep='\t',index=False)
print(seed,out,flush=True)
