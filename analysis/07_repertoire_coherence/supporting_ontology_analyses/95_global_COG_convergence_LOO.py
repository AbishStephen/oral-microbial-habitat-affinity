#!/usr/bin/env python3
from pathlib import Path
from itertools import product
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
P=pd.read_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t')
res=[]
for label,f in [('COG_all','COG2024_category_adjusted_scores_HMT47.tsv'),('COG_primary','COG2024_primarycat_adjusted_scores_HMT47.tsv'),('GO_BP','PFAM_GO_BP_adjusted_category_scores_HMT47.tsv')]:
 d=pd.read_csv(O/f,sep='\t').set_index('HMT')
 dd=d.loc[P.plaque_HMT].to_numpy(float)-d.loc[P.mucosal_HMT].to_numpy(float)
 for scale in ['unscaled','feature_standardized']:
  X=dd if scale=='unscaled' else dd/np.maximum(d.to_numpy(float).std(0),1e-5)
  for omit in ['NONE']+P.Genus.tolist():
   keep=np.array([True]*len(X)) if omit=='NONE' else (P.Genus!=omit).to_numpy()
   Y=X[keep];Y=Y/np.maximum(np.linalg.norm(Y,axis=1,keepdims=True),1e-15)
   ix,jx=np.triu_indices(len(Y),1);cos=(Y@Y.T)[ix,jx]
   signs=np.array(list(product([-1.,1.],repeat=len(Y))))
   null=((signs[:,ix]*signs[:,jx])@cos)/len(cos)
   obs=float(cos.mean())
   res.append(dict(ontology=label,scale=scale,excluded_genus=omit,n_pairs=len(Y),cosine=obs,one_sided_p=float(np.mean(null>=obs-1e-12))))
out=pd.DataFrame(res);out.to_csv(O/'global_COG_sister_convergence_leave_one_genus_out.tsv',sep='\t',index=False)
print(out[(out.ontology=='COG_primary')&(out.scale=='feature_standardized')].to_string(index=False))
