#!/usr/bin/env python3
"""Nested leave-one-genus-out genus-balanced phenotype prediction after removing genome
size, GC, CDS associations from features INSIDE EVERY inner/outer training fold.
"""
from pathlib import Path
import csv,numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
src=Path('__PROJECT_ROOT__/study13/scripts/104_nested_genusbalanced.py')
exec(compile(src.read_text().split("for label,y in [('original',yp),('sourcebalanced',y_alt)]:")[0],str(src),'exec'))
m=pd.read_csv(V/'confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT').loc[hmt]
C=m[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float)
def adjusted_X(train):
 mu=C[train].mean(axis=0);sd=C[train].std(axis=0)
 sd[sd<1e-9]=1
 Z=np.column_stack([np.ones(len(C)),(C-mu)/sd])
 bet=np.linalg.pinv(Z[train])@X[train]
 return (X-Z@bet).astype(np.float32)
def rank_adjusted(XX,train):
 xx=XX[train];gg=gen[train];ee=eco[train]
 ds=[]
 for g in np.unique(gg):
  zz=xx[gg==g];lab=ee[gg==g];ds.append(zz[lab=='plaque'].mean(axis=0)-zz[lab=='mucosal'].mean(axis=0))
 d=np.stack(ds);v=d.mean(axis=0);pos=(d>.05).sum(0);neg=(d<-.05).sum(0)
 maj=np.maximum(pos,neg);opp=np.minimum(pos,neg)
 sc=abs(v)/(0.15+d.std(axis=0))*np.maximum(0,(maj-opp)/len(ds))**2*(.2+np.sqrt(maj))
 sc[maj<2]=0
 sc[np.ptp(xx,axis=0)<1e-7]=0
 return np.lexsort((np.arange(len(sc)),-sc))
def pred_adjusted(XX,tr,te,sel,k,a,y):
 feat=sel[:k]
 model=make_pipeline(StandardScaler(),Ridge(alpha=a))
 model.fit(XX[tr][:,feat],y[tr])
 return model.predict(XX[te][:,feat])
for label,y in [('original',yp),('sourcebalanced',y_alt)]:
 outer_pred=np.full(len(y),np.nan);params=[]
 for og in gs:
  tr=np.where(gen!=og)[0];te=np.where(gen==og)[0]
  inn={}
  for ig in np.unique(gen[tr]):
   it=tr[gen[tr]==ig];ix=tr[gen[tr]!=ig]
   zz=adjusted_X(ix);inn[ig]=(ix,it,zz,rank_adjusted(zz,ix))
  best=None
  for k in Ks:
   for a in alphas:
    p=np.full(len(y),np.nan)
    for ig,(ix,it,zz,r) in inn.items():p[it]=pred_adjusted(zz,ix,it,r,k,a,y)
    auc=auc_mean(tr,p[tr]);rho=float(spearmanr(y[tr],p[tr]).statistic)
    if not np.isfinite(rho):rho=-1
    loss=float(np.mean(abs(y[tr]-p[tr])))
    key=(auc,rho,-loss,-k,-a)
    if best is None or key>best[0]:best=(key,k,a)
  k,a=best[1:];zz=adjusted_X(tr);sel=rank_adjusted(zz,tr);outer_pred[te]=pred_adjusted(zz,tr,te,sel,k,a,y)
  params.append(dict(Genus=og,k=k,alpha=a,inner_auc=best[0][0],test_auc=auc_mean(te,outer_pred[te])))
  print('FOLD',label,og,k,a,params[-1]['test_auc'],flush=True)
 pd.DataFrame(dict(HMT=hmt,Genus=gen,Ecology=eco,Y=y,Y_pred=outer_pred)).to_csv(O/f'nested_genusbalanced_covariate_residualized_{label}.tsv',sep='\t',index=False)
 pd.DataFrame(params).to_csv(O/f'nested_genusbalanced_covariate_residualized_{label}_folds.tsv',sep='\t',index=False)
 print('COMPLETE',label,'AUC',auc_mean(np.arange(len(y)),outer_pred),'rho',spearmanr(y,outer_pred).statistic,flush=True)
