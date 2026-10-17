#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'17_SnoaL2_independent_targets_20261016'
a=pd.read_csv(O/'SnoaL2_47HMT_presence_by_architecture.tsv',sep='\t')
a['either']= (a.PF12680_only>0)|(a['PF07366+PF12680']>0)|(a['PF12680+PF13577']>0)
cols=['PF12680_only','PF07366+PF12680','linked_HxlR_3kb','either']
rng=np.random.default_rng(20261016);out=[]
for feature in cols:
 sub=a.copy();sub['positive']=sub[feature].astype(bool)
 groups=[x for name,x in sub.groupby('Genus') if x.Ecology.nunique()==2]
 obs=np.mean([float(x.loc[x.Ecology=='plaque','positive'].mean()-x.loc[x.Ecology=='mucosal','positive'].mean()) for x in groups])
 z=np.zeros(10000)
 for t in groups:
  n=int((t.Ecology=='plaque').sum());v=t.positive.to_numpy(float)
  for k in range(10000):
   perm=rng.permutation(len(v));z[k]+=(v[perm[:n]].mean()-v[perm[n:]].mean())/len(groups)
 p=(int(sum(abs(z)>=abs(obs)-1e-12))+1)/10001
 out.append(dict(feature=feature,plaque_HMT=int(sum(sub.loc[sub.Ecology=='plaque','positive'])),plaque_total=int(sum(sub.Ecology=='plaque')),
 mucosal_HMT=int(sum(sub.loc[sub.Ecology=='mucosal','positive'])),mucosal_total=int(sum(sub.Ecology=='mucosal')),n_informative_genera=len(groups),genus_balanced_difference=obs,permutation_p=p))
res=pd.DataFrame(out);res.to_csv(O/'SnoaL2_47HMT_genus_block_null.tsv',sep='\t',index=False)
print(res.to_string(index=False))
