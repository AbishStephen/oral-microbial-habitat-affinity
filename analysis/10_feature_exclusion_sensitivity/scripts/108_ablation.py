#!/usr/bin/env python3
"""Full-nested domain/locus ablation; evaluates prediction reliance on linked annotated features.
Each scenario excludes annotated family globally (as a pre-specified sensitivity).
"""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
src=Path('__PROJECT_ROOT__/study13/scripts/104_nested_genusbalanced.py')
exec(compile(src.read_text().split("for label,y in [('original',yp),('sourcebalanced',y_alt)]:")[0],str(src),'exec'))
C=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_feature_integrated_raw_adjusted_annotated.tsv',sep='\t')
cand=set(C[C.support_selection==True].feature)
cases={
 'remove_LCM':{'PFAM:PF04072'},
 'remove_HxlR':{'PFAM:PF01638'},
 'remove_SnoaL2_SnoaL':{'PFAM:PF07366','PFAM:PF12680','PFAM:PF13577','PFAM:PF10604'},
 'remove_HxlR_SnoaL_locus':{'PFAM:PF01638','PFAM:PF07366','PFAM:PF12680','PFAM:PF13577','PFAM:PF10604'},
 'remove_LCM_plus_HxlR_SnoaL':{'PFAM:PF04072','PFAM:PF01638','PFAM:PF07366','PFAM:PF12680','PFAM:PF13577','PFAM:PF10604'},
 'remove_all_14_candidates':cand | {'PFAM:PF07366','PFAM:PF12680','PFAM:PF10604'}
}
base_rank=rank
rows_out=[];outerpred=[]
for case,ban in cases.items():
 ban_idx=np.isin(names,list(ban))
 assert ban_idx.any()
 def rank(train_idx):
  ranked=base_rank(train_idx)
  return ranked[~ban_idx[ranked]]
 for label,y in [('original',yp),('sourcebalanced',y_alt)]:
  preds=np.full(len(y),np.nan)
  for outer in gs:
   tr=np.where(gen!=outer)[0];te=np.where(gen==outer)[0];inner={}
   for inn in np.unique(gen[tr]):
    it=tr[gen[tr]==inn];ix=tr[gen[tr]!=inn];inner[inn]=(ix,it,rank(ix))
   best=None
   for k in Ks:
    for a in alphas:
     ip=np.full(len(y),np.nan)
     for inn,(ix,it,rr) in inner.items():ip[it]=fitpredict(ix,it,rr,k,a,y)
     score=auc_mean(tr,ip[tr]);rho=float(spearmanr(y[tr],ip[tr]).statistic)
     if not np.isfinite(rho):rho=-1
     err=float(np.mean(abs(y[tr]-ip[tr])));key=(score,rho,-err,-k,-a)
     if best is None or key>best[0]:best=(key,k,a)
   k,a=best[1:];rr=rank(tr);preds[te]=fitpredict(tr,te,rr,k,a,y)
  rows_out.append(dict(ablated=case,label=label,n_excluded=len(ban),macro_within_genus_auc=auc_mean(np.arange(len(y)),preds),pooled_spearman=float(spearmanr(y,preds).statistic)))
  outerpred.extend(dict(HMT=hmt[i],Genus=gen[i],Ecology=eco[i],label=label,ablated=case,Y=y[i],Y_pred=preds[i]) for i in range(len(y)))
  print('DONE',case,label,rows_out[-1],flush=True)
pd.DataFrame(rows_out).to_csv(O/'nested_genusbalanced_locus_ablation_summary.tsv',sep='\t',index=False)
pd.DataFrame(outerpred).to_csv(O/'nested_genusbalanced_locus_ablation_predictions.tsv',sep='\t',index=False)
