#!/usr/bin/env python3
from pathlib import Path
from itertools import product
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
P=pd.read_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t')
E=pd.read_csv(O/'extreme13_genus_matched_HMT_pairs.tsv',sep='\t')
tabs={'COG2024':'COG2024_category_adjusted_scores_HMT47.tsv','COG2024_primary':'COG2024_primarycat_adjusted_scores_HMT47.tsv','PFAM_GO_BP':'PFAM_GO_BP_adjusted_category_scores_HMT47.tsv'}
rows=[]
for label,fname in tabs.items():
 D0=pd.read_csv(O/fname,sep='\t').set_index('HMT');a=D0.to_numpy(float)
 for pn,pp in [('sister8',P),('extreme13',E),('extreme7strong',E[E.both_strong_2p5])]:
  D=D0.loc[pp.plaque_HMT].to_numpy(float)-D0.loc[pp.mucosal_HMT].to_numpy(float)
  for variant in ['unscaled','feature_standardized']:
   d=D.copy()
   if variant=='feature_standardized':d=d/np.maximum(a.std(axis=0),.00001)
   good=np.linalg.norm(d,axis=1)>1e-12
   d=d[good]/np.linalg.norm(d[good],axis=1)[:,None]
   n=len(d);sim=(d@d.T)[np.triu_indices(n,1)]
   ii,jj=np.triu_indices(n,1)
   signs=np.array(list(product([-1.,1.],repeat=n)),dtype=float)
   null=(signs[:,ii]*signs[:,jj])@sim/len(sim)
   ob=float(sim.mean())
   rows.append(dict(ontology=label,panel=pn,scaling=variant,n_pairs=n,n_categories=len(D0.columns),
    mean_crossgenus_cosine=ob,median_cosine=float(np.median(sim)),
    cosine_positive=int(sum(sim>0)),cosine_negative=int(sum(sim<0)),
    exact_one_sided_convergence_p=float(np.mean(null>=ob-1e-12)),
    exact_two_sided_p=float(np.mean(abs(null)>=abs(ob)-1e-12)),
    null_95pct=float(np.quantile(null,.95))))
r=pd.DataFrame(rows);r.to_csv(O/'global_crossgenus_ontology_convergence_exact.tsv',sep='\t',index=False)
print(r.to_string(index=False))
