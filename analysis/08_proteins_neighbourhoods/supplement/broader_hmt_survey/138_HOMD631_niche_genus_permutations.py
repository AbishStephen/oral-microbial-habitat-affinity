#!/usr/bin/env python3
"""Genus-balanced test of expanded HOMD HxlR/SnoaL domain and locus prevalence
using conservatively site-labelled HOMD taxa.
Not an independent clinical cohort; derived from published related 16S studies.
One reference genome per HMT, not an all-strain prevalence estimate.
"""
from pathlib import Path
import pandas as pd,numpy as np,json
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
x=pd.read_csv(O/'HOMD631_HxlR_SnoaL2_4source_ecology_annotations.tsv',sep='\t')
q=pd.read_csv(O/'HOMD631_PHYSICAL_HxlR_SnoaL2_loci_within3kb.tsv',sep='\t')
x['physically_linked']=x.GCA.isin(q.GCA)
x['tight_divergent']=x.GCA.isin(q.loc[q.tight_divergent,'GCA'])
x['HxlR_present']=x.PF01638>0
x['SnoaL2_present']=x.PF12680>0
x['both_present']=x.both
orig=pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t')
x['new_HMT']=~x.HMT_link.isin(set(orig.HMT))
oral=x[(x.hmt_primary_body_site_w_abundance=='Oral')&(x.ecology_4source_conservative.isin(['plaque','mucosal']))].copy()
all_stats=[]
rng=np.random.default_rng(20261014)
features=['HxlR_present','SnoaL2_present','both_present','physically_linked','tight_divergent']
for label,sub in [('all_conservatively_labeled_oral',oral),('new_HMT_only',oral[oral.new_HMT])]:
 tab=sub.groupby('ecology_4source_conservative')[features].agg(['sum','count'])
 print(label,'n=',len(sub),'site counts=',sub.ecology_4source_conservative.value_counts().to_dict(),flush=True)
 print('PREVALENCE\n'+sub.groupby('ecology_4source_conservative')[features].agg(['sum','mean']).to_string(),flush=True)
 groups=[]
 for genus,df in sub.groupby('main_genus'):
  if set(df.ecology_4source_conservative)=={'plaque','mucosal'}:
   groups.append((genus,df))
 print('BOTH_SITE_GENERA',len(groups),[g for g,d in groups],flush=True)
 for feature in features:
  sign=[]
  for g,d in groups:
   y=(d.ecology_4source_conservative=='plaque').to_numpy()
   f=d[feature].to_numpy(dtype=float)
   obs=f[y].mean()-f[~y].mean()
   sign.append((g,obs,len(d),int(y.sum()),int((~y).sum())))
  if not sign:continue
  stat=float(np.mean([z[1] for z in sign]))
  sims=np.zeros(10000)
  for g,d in groups:
   f=d[feature].to_numpy(dtype=float); k=sum(d.ecology_4source_conservative=='plaque')
   for j in range(len(sims)):
    p=rng.permutation(len(d));sims[j]+=(f[p[:k]].mean()-f[p[k:]].mean())/len(groups)
  p2=float((sum(abs(sims)>=abs(stat)-1e-12)+1)/(len(sims)+1))
  all_stats.append({'set':label,'feature':feature,'n_HMT':len(sub),'n_genus_with_both_site_classes':len(groups),
  'genus_balanced_mean_plaque_minus_mucosa_presence':stat,'two_sided_within_genus_perm_p':p2,
  'plaque_positive':int(sub[sub.ecology_4source_conservative=='plaque'][feature].sum()),
  'plaque_total':int(sum(sub.ecology_4source_conservative=='plaque')),
  'mucosal_positive':int(sub[sub.ecology_4source_conservative=='mucosal'][feature].sum()),
  'mucosal_total':int(sum(sub.ecology_4source_conservative=='mucosal'))})
 print('GENUSPAIR_SITE_DIFFERENCES',feature,sign,flush=True)
r=pd.DataFrame(all_stats)
r.to_csv(O/'HOMD631_HxlR_SnoaL2_expanded_ecological_genus_control.tsv',sep='\t',index=False)
print('FINAL\n'+r.to_string(index=False),flush=True)
