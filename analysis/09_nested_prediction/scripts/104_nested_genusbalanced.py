#!/usr/bin/env python3
"""Leakage-free nested leave-one-GENUS-out genome-function transfer, genus balanced feature rank.
Train folds select and tune using 12 training genera only. No global feature pre-selection.
Evaluate separately original and source-family-balanced continuous Y.
"""
from pathlib import Path
import csv, numpy as np, pandas as pd, argparse, time
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr
B=Path('__PROJECT_ROOT__/study13')
V=B/'08_validation';O=B/'11_manuscript_completion_20261013';O.mkdir(exist_ok=True)
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 r=csv.reader(f,delimiter='\t');hd=next(r);rows=list(r)
names=np.asarray(hd[5:]); hmt=np.array([x[0] for x in rows]);gen=np.array([x[1] for x in rows])
tier=np.array([x['Ecology_tier'] for x in csv.DictReader(open(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv'),delimiter='\t')]) # manifest same row order? verify next
eco=np.array([x[3] for x in rows]);yp=np.asarray([float(x[4]) for x in rows]);X=np.array([[float(z) for z in r[5:]] for r in rows],dtype=np.float32)
src={r['HMT']:float(r['groupbalanced_Y_1em06']) for r in csv.DictReader(open(V/'confound_audit/ecology_source_family_pseudocount_sensitivity.tsv'),delimiter='\t')}
y_alt=np.array([src[h] for h in hmt])
meta={r['HMT']:r for r in csv.DictReader(open(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv'),delimiter='\t')}
tier=np.array([meta[h]['Ecology_tier'] for h in hmt])
gs=np.unique(gen)
def rank(train_idx):
 xx=X[train_idx];gg=gen[train_idx];yy=eco[train_idx]
 ds=[]
 for g in np.unique(gg):
  v=xx[gg==g];w=yy[gg==g]
  assert np.any(w=='plaque') and np.any(w=='mucosal')
  ds.append(v[w=='plaque'].mean(axis=0)-v[w=='mucosal'].mean(axis=0))
 d=np.stack(ds,axis=0)
 m=d.mean(axis=0);pos=np.sum(d>.05,axis=0);neg=np.sum(d<-.05,axis=0)
 maj=np.maximum(pos,neg); opp=np.minimum(pos,neg)
 # geometric mean combines large effect and directional consistency; feature support across genera
 score=(np.abs(m)/(0.15+d.std(axis=0)))*np.maximum(0,(maj-opp)/len(ds))**2*(0.2+np.sqrt(maj))
 score[maj<2]=0
 score[np.ptp(xx,axis=0)<1e-7]=0
 return np.lexsort((np.arange(len(score)),-score))
def auc_mean(idx,preds):
 gg=gen[idx];ep=eco[idx]; pp=preds
 per=[]
 for g in np.unique(gg):
  pr=pp[(gg==g)&(ep=='plaque')]; mu=pp[(gg==g)&(ep=='mucosal')]
  if len(pr) and len(mu):per.append(float(np.mean([(1. if a>b else (.5 if a==b else 0.)) for a in pr for b in mu])))
 return float(np.mean(per)) if per else np.nan
def fitpredict(train,test,ranked,k,alpha,y):
 sel=ranked[:min(k,len(ranked))]
 pipe=make_pipeline(StandardScaler(),Ridge(alpha=alpha))
 pipe.fit(X[train][:,sel],y[train])
 return pipe.predict(X[test][:,sel])
Ks=(5,15,40,100);alphas=(10.,100.,1000.)
for label,y in [('original',yp),('sourcebalanced',y_alt)]:
 print('BEGIN',label,flush=True)
 start=time.time()
 preds=np.full(len(y),np.nan);fold_rows=[];chosen=[];stability={}
 for outer in gs:
  tr=np.where(gen!=outer)[0];te=np.where(gen==outer)[0]
  inner={}
  # feature selection happens separately for every inner training split
  for inn in np.unique(gen[tr]):
   it=tr[gen[tr]==inn];ix=tr[gen[tr]!=inn]
   inner[inn]=(ix,it,rank(ix))
  best=None
  for k in Ks:
   for a in alphas:
    ix_pred=np.full(len(y),np.nan)
    for inn,(ix,it,rr) in inner.items():ix_pred[it]=fitpredict(ix,it,rr,k,a,y)
    score=auc_mean(tr,ix_pred[tr])
    rho=float(spearmanr(y[tr],ix_pred[tr]).statistic)
    if not np.isfinite(rho):rho=-1
    err=float(np.mean(abs(y[tr]-ix_pred[tr])))
    key=(score,rho,-err,-k,-a)
    if best is None or key>best[0]:best=(key,k,a)
  _,k,a=best
  rr=rank(tr);sel=rr[:k]
  preds[te]=fitpredict(tr,te,rr,k,a,y)
  for x in sel:stability[names[x]]=stability.get(names[x],0)+1
  fold_rows.append(dict(heldout_genus=outer,n_HMT=len(te),selected_k=k,alpha=a,genus_auc=auc_mean(te,preds[te]),n_correct_HMT=sum((preds[te]>0)==(eco[te]=='plaque')),inner_auc=best[0][0]))
  print('FOLD',label,outer,'k',k,'alpha',a,'inner',round(best[0][0],3),'test_auc',round(fold_rows[-1]['genus_auc'],3),'elapsed',round(time.time()-start,1),flush=True)
  for j in sel: chosen.append(dict(heldout_genus=outer,feature=names[j],rank_in_selected=list(sel).index(j)+1))
 assert np.isfinite(preds).all()
 rr=pd.DataFrame(dict(HMT=hmt,Genus=gen,Ecology=eco,Ecology_tier=tier,Y=y,Y_pred=preds))
 rr.to_csv(O/f'nested_genusbalanced_{label}_HMT_predictions.tsv',sep='\t',index=False)
 pd.DataFrame(fold_rows).to_csv(O/f'nested_genusbalanced_{label}_fold_results.tsv',sep='\t',index=False)
 pd.DataFrame(chosen).to_csv(O/f'nested_genusbalanced_{label}_selected_features_by_fold.tsv',sep='\t',index=False)
 pd.DataFrame([dict(feature=f,n_folds=n) for f,n in sorted(stability.items(),key=lambda z:-z[1])]).to_csv(O/f'nested_genusbalanced_{label}_feature_selection_stability.tsv',sep='\t',index=False)
 print('FINISHED',label,'macro_auc',auc_mean(np.arange(len(y)),preds),'allrho',spearmanr(y,preds).statistic,'high_tier_macro_auc',auc_mean(np.where(tier=='High')[0],preds[tier=='High']),'time',round(time.time()-start,1),flush=True)
