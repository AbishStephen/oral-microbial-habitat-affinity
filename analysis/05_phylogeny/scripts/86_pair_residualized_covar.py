#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
from itertools import product
B=Path('__PROJECT_ROOT__/study13');O=B/'10_additional_analyses_20261012'
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT')
g=d.Genus.to_numpy();G=pd.get_dummies(g,drop_first=True,dtype=float).to_numpy()
C=d[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(dtype=float)
C=(C-C.mean(0))/(C.std(0)+1e-12)
Z=np.column_stack([np.ones(len(d)),G,C])
Q,_,_=np.linalg.svd(Z,full_matrices=False)
Zrank=np.linalg.matrix_rank(Z);Q=Q[:,:Zrank]
candidate=['PF04072','PF12680','PF10604','PF08020','PF13199','PF22815']
X=d[candidate].to_numpy(float);XR=X-Q@(Q.T@X)
res=pd.DataFrame(XR,index=d.index,columns=candidate)
pairfiles=['sister8_pairwise_HMT_genome_domain_contrasts.tsv','extreme13_genus_matched_HMT_pairs.tsv']
out=[]
for fname in pairfiles:
 x=pd.read_csv(O/fname,sep='\t')
 for feat in candidate:
  vals=(res.loc[x.plaque_HMT,feat].to_numpy()-res.loc[x.mucosal_HMT,feat].to_numpy())
  raw=(d.loc[x.plaque_HMT,feat].to_numpy()-d.loc[x.mucosal_HMT,feat].to_numpy())
  sgns=np.array(list(product((-1,1),repeat=len(vals))))
  z=float(vals.sum()/np.sqrt((vals*vals).sum())) if np.linalg.norm(vals)>1e-12 else 0.
  zs=(sgns@vals)/(np.linalg.norm(vals)+1e-12)
  p=float(np.mean(abs(zs)>=abs(z)-1e-12))
  out.append(dict(panel=fname.replace('_HMT_genome_domain_contrasts.tsv','').replace('_genus_matched_HMT_pairs.tsv',''),Pfam=feat,n_pairs=len(vals),raw_plaque_higher=int(sum(raw>1e-9)),raw_mucosa_higher=int(sum(raw< -1e-9)),adjusted_plaque_higher=int(sum(vals>1e-9)),adjusted_mucosa_higher=int(sum(vals< -1e-9)),adjusted_mean_delta=float(vals.mean()),exact_adjusted_exploratory_p=p))
pd.DataFrame(out).to_csv(O/'sister_extreme_covariate_residualized_candidate_sensitivity.tsv',sep='\t',index=False)
print('RESIDUALIZED',Zrank,'parameters\n',pd.DataFrame(out).to_string(index=False))
