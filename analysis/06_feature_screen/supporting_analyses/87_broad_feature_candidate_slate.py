#!/usr/bin/env python3
"""Feature candidate slate across 8 sister and 13 extreme paired contrasts; descriptive, no multiple-testing claims."""
from pathlib import Path
import gzip,csv
import numpy as np
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
pairs=pd.read_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t')
ex=pd.read_csv(O/'extreme13_genus_matched_HMT_pairs.tsv',sep='\t')
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');h=next(rd);rr=list(rd)
feat=np.array(h[5:],object);ids=[z[0] for z in rr];ii={z:i for i,z in enumerate(ids)}
mat=np.array([[float(v) for v in z[5:]] for z in rr],dtype=float)
D1=mat[[ii[z] for z in pairs.plaque_HMT]]-mat[[ii[z] for z in pairs.mucosal_HMT]]
D2=mat[[ii[z] for z in ex.plaque_HMT]]-mat[[ii[z] for z in ex.mucosal_HMT]]
mt=pd.read_csv(B/'08_validation/confound_audit/genome_covariate_adjusted_maxT_all_features.tsv',sep='\t').set_index('feature')
pfam={}
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  if row['Pfam'] not in pfam:pfam[row['Pfam']]=row['Pfam_name']
names={}
for k in feat:
 if k.startswith('PFAM:'):names[k]=pfam.get(k[5:],'')
for file,pref in [('kegg_pathway_names.tsv','KPATH:'),('kegg_module_names.tsv','KMOD:')]:
 for line in (B/'05_function/mappings'/file).read_text().splitlines():
  z=line.split('\t',1)
  if len(z)>1:names[pref+z[0]]=z[1]
others=pd.read_csv(B/'08_validation/within_genus_expanded13/all_features_expanded13_within_genus.tsv',sep='\t',usecols=['feature','annotation']).drop_duplicates('feature')
names.update({k:v for k,v in zip(others.feature,others.annotation) if isinstance(v,str) and v})
mu1=D1.mean(0);mu2=D2.mean(0)
p1=(D1>1e-8).sum(0);n1=(D1< -1e-8).sum(0)
p2=(D2>1e-8).sum(0);n2=(D2< -1e-8).sum(0)
scores=(np.maximum(p1,n1)/8+.4*np.maximum(p2,n2)/13)*np.sqrt(np.abs(mu1)*np.abs(mu2))
out=pd.DataFrame(dict(feature=feat,layer=[v.split(':',1)[0] for v in feat],name=[names.get(v,'') for v in feat],sister_plaque=p1,sister_mucosal=n1,sister_nonzero=p1+n1,sister_mean=mu1,extreme_plaque=p2,extreme_mucosal=n2,extreme_nonzero=p2+n2,extreme_mean=mu2,rank_score=scores))
out=out.merge(mt,how='left',left_on='feature',right_index=True)
out['abs_adjr']=out.genus_genomesize_GC_CDS_partial_r.abs()
out=out.sort_values('rank_score',ascending=False)
out.to_csv(O/'broad_feature_candidate_slate_all13917.tsv',sep='\t',index=False)
for key,z in [('PFAM',out[(out.layer=='PFAM') & (out.sister_nonzero>=4) & (out.extreme_nonzero>=6)].head(110)),('COG',out[(out.layer=='COG') & (out.sister_nonzero>=4)].head(45)),('KMOD',out[(out.layer=='KMOD') & (out.sister_nonzero>=3)].head(40)),('KPATH',out[(out.layer=='KPATH') & (out.sister_nonzero>=3)].head(35)),('EC',out[(out.layer=='EC') & (out.sister_nonzero>=3)].head(40))]:
 print('\n',key,'N',len(z),flush=True)
 print(z[['feature','name','sister_plaque','sister_mucosal','sister_mean','extreme_plaque','extreme_mucosal','extreme_mean','genus_genomesize_GC_CDS_partial_r','maxT_FWER_p']].to_string(index=False,max_rows=120),flush=True)
print('FINISHED',len(out),flush=True)
