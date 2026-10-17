#!/usr/bin/env python3
"""All 14 candidate feature site-definition robustness; descriptive, within-genus permutations.
Site abundances are correlated with original label: NOT independent validation.
"""
from pathlib import Path
import csv,numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'11_manuscript_completion_20261013';O.mkdir(exist_ok=True)
A=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_feature_integrated_raw_adjusted_annotated.tsv',sep='\t')
C=A[A.support_selection==True].copy();assert len(C)==14
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 r=csv.reader(f,delimiter='\t');header=next(r);rows=list(r)
names=header[5:];idx={name:i+5 for i,name in enumerate(names)}
hmt=[x[0] for x in rows]
mat=pd.DataFrame([{c:x[idx[c]] for c in C.feature} for x in rows],index=hmt,dtype=float)
site=pd.read_csv(B/'10_additional_analyses_20261012/site_gradient_HMT47_scores.tsv',sep='\t').set_index('HMT').loc[hmt]
assert (site.index==mat.index).all()
gen=site.Genus.values
groups=[np.where(gen==g)[0] for g in sorted(set(gen))]
X=mat.values
Xc=X.copy()
for ids in groups:Xc[ids]-=X[ids].mean(axis=0)
nperm=10000
rg=np.random.default_rng(20261013)
site_names=['SUPP_vs_TDO','SUBP_vs_TDO','SUPP_vs_BMU','SUBP_vs_SUPP','TDO_vs_BMU']
direction=C.set_index('feature').direction
raw=[]
for st in site_names:
 y=site[st].values.astype(float);yc=y.copy()
 for ids in groups:yc[ids]-=y[ids].mean()
 norm=np.sqrt(np.sum(Xc**2,axis=0)*np.sum(yc**2))
 obs=Xc.T@yc/np.maximum(norm,1e-10)
 null=np.zeros((nperm,len(C)),dtype=float)
 for k in range(nperm):
  pp=yc.copy()
  for ids in groups:pp[ids]=yc[ids][rg.permutation(len(ids))]
  null[k]=Xc.T@pp/np.maximum(norm,1e-10)
 maxAbs=np.max(abs(null),axis=1)
 for i,f in enumerate(C.feature):
  raw.append(dict(feature=f,annotation=C.iloc[i].annotation,expected_ecology=direction[f],site_contrast=st,
  genus_centered_pearson_r=obs[i],positive_direction_agreement=(obs[i]>0)==(direction[f]=='plaque'),
  within_genus_perm_two_sided=(1+np.count_nonzero(abs(null[:,i])>=abs(obs[i])-1e-12))/(nperm+1),
  within_14_candidate_maxT_P=(1+np.count_nonzero(maxAbs>=abs(obs[i])-1e-12))/(nperm+1),
  n_HMT=len(hmt),n_genera=len(groups)))
T=pd.DataFrame(raw)
T.to_csv(O/'all14_site_gradient_withingenus_perm.tsv',sep='\t',index=False)
summary=[]
for feature,z in T.groupby('feature'):
 q=z.set_index('site_contrast')
 expected=direction[feature]
 main=q.loc[['SUPP_vs_TDO','SUBP_vs_TDO','SUPP_vs_BMU']]
 summary.append(dict(feature=feature,annotation=z.iloc[0].annotation,direction=expected,n_expected_direction=int(main.positive_direction_agreement.sum()),
  n_main_contrasts=3,mean_abs_correlation=float(main.genus_centered_pearson_r.abs().mean()),
  n_nominal_p_lt_005=int((main.within_genus_perm_two_sided<.05).sum()),
  n_maxT14_p_lt_005=int((main.within_14_candidate_maxT_P<.05).sum()),
  site_gradient_effects='; '.join(f'{nm}:r={q.loc[nm,"genus_centered_pearson_r"]:+.2f}' for nm in ['SUPP_vs_TDO','SUBP_vs_TDO','SUPP_vs_BMU'])))
S=pd.DataFrame(summary).sort_values(['n_expected_direction','mean_abs_correlation'],ascending=False)
S.to_csv(O/'all14_site_gradient_robustness_summary.tsv',sep='\t',index=False)
print('DONE 14x5 contrast site sensitivity and restricted-candidate maxT; no cohort independence',flush=True)
print(S.to_string(index=False),flush=True)
