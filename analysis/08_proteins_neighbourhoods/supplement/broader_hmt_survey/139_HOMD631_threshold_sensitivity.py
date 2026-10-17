#!/usr/bin/env python3
"""Sensitivity of HxlR/SnoaL niche patterns to conservative ecology-label thresholds.
One representative genome per HMT; site abundance studies are correlated;
all analyses exploratory and not phylogenetically independent.
"""
from pathlib import Path
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
t=pd.read_csv(O/'HOMD631_HxlR_SnoaL2_4source_ecology_annotations.tsv',sep='\t')
q=pd.read_csv(O/'HOMD631_PHYSICAL_HxlR_SnoaL2_loci_within3kb.tsv',sep='\t')
d=set(pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t').HMT)
t=t[t.hmt_primary_body_site_w_abundance=='Oral'].copy()
t['new']=~t.HMT_link.isin(d)
t['linked']=t.GCA.isin(q.GCA)
sources=['erenv1v3','erenv3v5','hmprefseqv1v3','hmprefseqv3v5']
A=t[[x+'_score' for x in sources]].to_numpy(float).copy()
U=t[[x+'_usable' for x in sources]].fillna(False).to_numpy(bool)
A[~U]=np.nan
features={'HxlR':(t.PF01638>0).to_numpy(int),'SnoaL_2':(t.PF12680>0).to_numpy(int),'physical_locus':t.linked.to_numpy(int)}
results=[]
rng=np.random.default_rng(20261014)
for threshold in [0.5,1,2]:
 for min_support in [2,3]:
  pos=(A>=threshold).sum(axis=1);neg=(A<=-threshold).sum(axis=1)
  label=np.where((pos>=min_support)&(neg==0),'plaque',np.where((neg>=min_support)&(pos==0),'mucosal','undetermined'))
  for only_new in [False,True]:
   keep=(label!='undetermined')&((t.new.to_numpy(bool)) if only_new else True)
   taxa=t.loc[keep].copy();taxa['eco']=label[keep]
   if len(taxa)<10:continue
   grp=[z for g,z in taxa.groupby('main_genus') if z.eco.nunique()==2]
   if len(grp)<2:continue
   for feat,x in features.items():
    o=x[keep];taxa['present']=o
    groups=[z for g,z in taxa.groupby('main_genus') if z.eco.nunique()==2]
    obs=np.mean([z[z.eco=='plaque'].present.mean()-z[z.eco=='mucosal'].present.mean() for z in groups])
    null=np.zeros(2000)
    for z in groups:
     values=z.present.to_numpy(int);n=(z.eco=='plaque').sum()
     for ii in range(2000):
      idx=rng.permutation(len(values))
      null[ii]+=(values[idx[:n]].mean()-values[idx[n:]].mean())/len(groups)
    p=(sum(abs(null)>=abs(obs)-1e-10)+1)/2001
    results.append(dict(threshold_log2=threshold,min_source_support=min_support,new_HMT_only=only_new,n_labeled_HMT=len(taxa),n_genus=len(groups),feature=feat,genus_diff=float(obs),perm_p_2sided=float(p),
      n_plaque=int(sum(taxa.eco=='plaque')),n_mucosal=int(sum(taxa.eco=='mucosal'))))
r=pd.DataFrame(results);r.to_csv(O/'HOMD631_HxlR_SnoaL2_ecology_threshold_sensitivity.tsv',sep='\t',index=False)
print(r.to_string(index=False))
