#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'17_SnoaL2_independent_targets_20261016'
H=B/'14_HxlR_fullHOMD_controls_20261014'
d=pd.read_csv(H/'HOMD631_HxlR_SnoaL2_independent_HMM_protein_hits.tsv',sep='\t')
p=O/'HOMD631_PF07366_domtblout.tsv'
out=[]
for line in open(p):
 if line.startswith('#'):continue
 z=line.split()
 if len(z)>=19:out.append({'GCA':z[0].split('|')[0],'protein_ID':z[0],'Pfam':'PF07366','domain_iEvalue':float(z[12])})
p73=pd.DataFrame(out)
s=d[d.Pfam=='PF12680'].drop_duplicates(['GCA','protein_ID'])
pf73=set(p73.protein_ID)
dual=s[s.protein_ID.isin(pf73)].copy()
single=s[~s.protein_ID.isin(pf73)].copy()
c=pd.read_csv(H/'HOMD631_HxlR_SnoaL2_4source_ecology_annotations.tsv',sep='\t')
c['dual_PF07366_PF12680']=c.GCA.isin(dual.GCA)
c['single_PF12680_only']=c.GCA.isin(single.GCA)
c['PF07366_present']=c.GCA.isin(p73.GCA)
c.to_csv(O/'HOMD631_SnoaL2_architecture_and_ecology.tsv',sep='\t',index=False)
orals=c[c.hmt_primary_body_site_w_abundance=='Oral']
e=orals[orals.ecology_4source_conservative.isin(['mucosal','plaque'])].copy()
orig=set(pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t').HMT)
e['new_HMT']=~e.HMT_link.isin(orig)
rng=np.random.default_rng(20261016);records=[]
for view,sub in [('all_conservative_oral',e),('new_HMT_only',e[e.new_HMT])]:
 groups=[x for _,x in sub.groupby('main_genus') if x.ecology_4source_conservative.nunique()==2]
 for feature in ['dual_PF07366_PF12680','single_PF12680_only','PF07366_present']:
  obs=float(np.mean([x.loc[x.ecology_4source_conservative=='plaque',feature].mean()-x.loc[x.ecology_4source_conservative=='mucosal',feature].mean() for x in groups]))
  null=np.zeros(5000)
  for x in groups:
   v=x[feature].to_numpy(float);n=int((x.ecology_4source_conservative=='plaque').sum())
   for j in range(len(null)):
    idx=rng.permutation(len(v));null[j]+=(v[idx[:n]].mean()-v[idx[n:]].mean())/len(groups)
  p=(np.sum(abs(null)>=abs(obs)-1e-12)+1)/(len(null)+1)
  records.append({'view':view,'feature':feature,'HMT_count':len(sub),'informative_genera':len(groups),
    'plaque_positive':int(sub[sub.ecology_4source_conservative=='plaque'][feature].sum()),'plaque_total':int(sum(sub.ecology_4source_conservative=='plaque')),
    'mucosal_positive':int(sub[sub.ecology_4source_conservative=='mucosal'][feature].sum()),'mucosal_total':int(sum(sub.ecology_4source_conservative=='mucosal')),
    'genus_balanced_difference':obs,'two_sided_perm_p':float(p)})
res=pd.DataFrame(records)
res.to_csv(O/'HOMD631_SnoaL2_architecture_genus_permutation.tsv',sep='\t',index=False)
print('TOTAL_PF07366',p73.protein_ID.nunique(),'DUAL',len(dual),'SINGLE',len(single),flush=True)
print('GENUS_CONTROL\n'+res.to_string(index=False),flush=True)
