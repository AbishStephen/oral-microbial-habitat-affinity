#!/usr/bin/env python3
from pathlib import Path
import csv,numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
from sklearn.metrics import roc_auc_score
B=Path('__PROJECT_ROOT__');E=B/'study13'
O=E/'08_validation/confound_audit';O.mkdir(exist_ok=True,parents=True)
d=pd.read_csv(O/'hmt_covariates_and_predictions.tsv',sep='\t')
files={
'Dewhirst':'abundance_Dewhirst.tsv',
'Eren_V1V3':'abundance_Eren_V1V3.tsv',
'Eren_V3V5':'abundance_Eren_V3V5.tsv',
'HMP_MetaPhlan':'abundance_HMP_MetaPhlan.tsv',
'HMP_RefSeq_V1V3':'abundance_HMP_RefSeq_V1V3.tsv',
'HMP_RefSeq_V3V5':'abundance_HMP_RefSeq_V3V5.tsv'}
groups={'Dewhirst':'Dewhirst','Eren_V1V3':'Eren','Eren_V3V5':'Eren',
'HMP_MetaPhlan':'HMP','HMP_RefSeq_V1V3':'HMP','HMP_RefSeq_V3V5':'HMP'}
wanted=set(d.HMT)
allr=[]
for source,fn in files.items():
 with open(B/'00_metadata'/fn,encoding='utf-8-sig') as f:
  for line in f:
   if line.startswith('TAX\tHMT\t'):break
  records=csv.DictReader(f,fieldnames=line.strip().split('\t'),delimiter='\t')
  for r in records:
   try:h=f'HMT-{int(float(r["HMT"])):03d}'
   except:continue
   if h not in wanted or (h=='HMT-575' and source=='HMP_MetaPhlan'):continue
   vals={}
   for key in ['SUBP','SUPP','BMU','TDO']:
    try:vals[key]=float(r.get(key+'_mean','0') or 0)
    except:vals[key]=0
   pl=vals['SUBP']+vals['SUPP'];mu=vals['BMU']+vals['TDO']
   allr.append(dict(HMT=h,source=source,source_group=groups[source],plaque=pl,mucosa=mu,focal=pl+mu))
long=pd.DataFrame(allr)
bad=pd.read_csv(B/"08_validation/ecology_source_sensitivity/source_ambiguity_exclusions.tsv",sep="\t")
ban=set(zip(bad.source,bad.HMT))
long=long.loc[[ (s,h) not in ban for s,h in zip(long.source,long.HMT) ]].copy()
assert long.HMT.nunique()==47
for eps in [1e-8,1e-6,1e-4,1e-2,1e-1]:
 nm='Y_'+format(eps,'.0e').replace('-','m')
 long[nm]=np.log2((long.plaque+eps)/(long.mucosa+eps))
# source group median => 3 cohort families equally weighted
res=[]
for h,z in long.groupby('HMT'):
 row=dict(HMT=h,n_sources=len(z),n_source_groups=z.source_group.nunique())
 for eps in [1e-8,1e-6,1e-4,1e-2,1e-1]:
  nm='Y_'+format(eps,'.0e').replace('-','m')
  row['six_'+nm]=float(z[nm].median())
  row['groupbalanced_'+nm]=float(z.groupby('source_group')[nm].median().median())
 row['group_source_labels']=';'.join(sorted(z.source_group.unique()))
 row['nonzero_plaque_count']=int((z.plaque>0).sum());row['nonzero_mucosa_count']=int((z.mucosa>0).sum())
 res.append(row)
alt=pd.DataFrame(res).merge(d[['HMT','Genus','Ecology','Y','Y_pred','PF04072','PF12680','PF10604','PF13199','PF22815']],on='HMT',validate='one_to_one')
baseline=alt['six_Y_1em06']
print('CHECK original vs recomputed',float(np.max(np.abs(alt.Y-baseline))))
summary=[]
for name in [c for c in alt.columns if c.startswith(('six_Y','groupbalanced_Y'))]:
 x=alt[name].to_numpy()
 correct=((x>0)==(alt.Y>0)).sum()
 # within-genus pair directional concordance for existing predictor (frozen model only)
 pear=pearsonr(x,alt.Y)[0]
 rank=spearmanr(x,alt.Y).statistic
 centered=x-alt.groupby('Genus')[name].transform('mean').to_numpy()
 m='groupbalanced' if name.startswith('groupbalanced') else 'six'
 rs=[name,int(correct),pear,rank,spearmanr(x,alt.Y_pred).statistic,float(pearsonr(centered,alt['PF04072']-alt.groupby('Genus')['PF04072'].transform('mean')).statistic),float(pearsonr(centered,alt['PF12680']-alt.groupby('Genus')['PF12680'].transform('mean')).statistic)]
 summary.append(rs)
alt.to_csv(O/'ecology_source_family_pseudocount_sensitivity.tsv',sep='\t',index=False)
pd.DataFrame(summary,columns=['variant','n_matching_original_sign','pearson_origY','spearman_origY','spearman_frozen_pred','genus_resid_LCM_pearson','genus_resid_SnoaL_pearson']).to_csv(O/'ecology_source_family_sensitivity_summary.tsv',sep='\t',index=False)
print('N_SOURCES_COUNTS',alt.n_sources.value_counts().to_dict(),'N_SOURCE_GROUPS',alt.n_source_groups.value_counts().to_dict())
print('ALT_SUMMARY')
print(pd.read_csv(O/'ecology_source_family_sensitivity_summary.tsv',sep='\t').to_string(index=False))
print('BIGGEST SHIFTS groupbalanced epsilon1e-6')
k='groupbalanced_Y_1em06'
print(alt.assign(delta=alt[k]-alt.Y).sort_values('delta',key=abs,ascending=False)[['HMT','Genus','Ecology','Y',k,'delta']].head(15).to_string(index=False))
