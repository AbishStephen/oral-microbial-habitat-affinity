#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import spearmanr
B=Path('__PROJECT_ROOT__/study13');O=B/'10_additional_analyses_20261012'
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t')
rows=[]
for g,x in d.groupby('Genus'):
 p=x[x.Y>0];m=x[x.Y<0]
 n=len(p)*len(m)
 auc=float(sum(a>b for a in p.Y_pred for b in m.Y_pred)+.5*sum(a==b for a in p.Y_pred for b in m.Y_pred))/n
 sz=float(sum(a>b for a in p.size_predict for b in m.size_predict)+.5*sum(a==b for a in p.size_predict for b in m.size_predict))/n
 rho=float(spearmanr(x.Y,x.Y_pred).statistic) if len(x)>2 else np.nan
 rows.append(dict(Genus=g,n_HMT=len(x),n_plaque=len(p),n_mucosal=len(m),n_pairs=n,LOGO_pairwise_auc=auc,genomesize_pairwise_auc=sz,LOGO_rho_within_genus=rho))
z=pd.DataFrame(rows).sort_values('Genus')
z.to_csv(O/'existing_nestedLOGO_per_genus_performance.tsv',sep='\t',index=False)
print(z.to_string(index=False))
print('MACRO',z.LOGO_pairwise_auc.mean(),'SIZE',z.genomesize_pairwise_auc.mean(),'WEIGHT',np.average(z.LOGO_pairwise_auc,weights=z.n_pairs),'TOTAL PAIRS',z.n_pairs.sum())
