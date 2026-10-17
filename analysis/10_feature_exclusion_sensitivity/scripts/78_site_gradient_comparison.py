#!/usr/bin/env python3
"""Site-specific reanalysis of six correlated HOMD summaries; exploratory contrasts only."""
from pathlib import Path
import csv,numpy as np,pandas as pd
from scipy.stats import pearsonr,spearmanr
B=Path('__PROJECT_ROOT__');E=B/'study13'
O=E/'10_additional_analyses_20261012';O.mkdir(exist_ok=True)
ref=pd.read_csv(E/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t')
sources={'Dewhirst':'Dewhirst','Eren_V1V3':'Eren','Eren_V3V5':'Eren','HMP_MetaPhlan':'HMP','HMP_RefSeq_V1V3':'HMP','HMP_RefSeq_V3V5':'HMP'}
ban=pd.read_csv(B/'08_validation/ecology_source_sensitivity/source_ambiguity_exclusions.tsv',sep='\t')
bad=set(zip(ban.source,ban.HMT))
rr=[]
for source,grp in sources.items():
 f=B/'00_metadata'/('abundance_'+source+'.tsv')
 with open(f,encoding='utf-8-sig') as fi:
  for ln in fi:
   if ln.startswith('TAX\tHMT\t'):break
  rd=csv.DictReader(fi,fieldnames=ln.strip().split('\t'),delimiter='\t')
  for row in rd:
   try:h=f'HMT-{int(float(row["HMT"])):03d}'
   except:continue
   if h not in set(ref.HMT) or (source,h) in bad or (h=='HMT-575' and source=='HMP_MetaPhlan'):continue
   z={'HMT':h,'source':source,'group':grp}
   for si in ['SUBP','SUPP','BMU','TDO']:
    try:z[si]=float(row[si+'_mean'] or 0)
    except:z[si]=0.
   rr.append(z)
x=pd.DataFrame(rr)
pairs={'SUPP_vs_TDO':('SUPP','TDO'),'SUBP_vs_TDO':('SUBP','TDO'),'SUPP_vs_BMU':('SUPP','BMU'),'SUBP_vs_SUPP':('SUBP','SUPP'),'TDO_vs_BMU':('TDO','BMU')}
for nm,(a,b) in pairs.items():
 x[nm]=np.log2((x[a]+1e-6)/(x[b]+1e-6))
out=[]
for h,z in x.groupby('HMT'):
 rec=dict(HMT=h,n_sources=z.source.nunique(),n_source_groups=z.group.nunique())
 for nm in pairs:
  rec[nm]=float(z.groupby('group')[nm].median().median())
 out.append(rec)
d=pd.DataFrame(out).merge(ref[['HMT','Genus','Y','Y_pred','PF04072','PF12680','PF10604']],on='HMT')
rows=[]
for nm in pairs:
 for feat in ['Y_pred','PF04072','PF12680','PF10604']:
  a=d[nm]-d.groupby('Genus')[nm].transform('mean')
  b=d[feat]-d.groupby('Genus')[feat].transform('mean')
  rows.append(dict(site_contrast=nm,feature=feat,genus_centered_pearson=pearsonr(a,b).statistic,genus_centered_spearman=spearmanr(a,b).statistic,pooled_spearman=spearmanr(d[nm],d[feat]).statistic,n_HMT=len(d),n_zero_at_both_sites=int(sum((x[pairs[nm][0]]==0)&(x[pairs[nm][1]]==0)))))
d.to_csv(O/'site_gradient_HMT47_scores.tsv',sep='\t',index=False)
pd.DataFrame(rows).to_csv(O/'site_gradient_genus_centered_correlations.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False))
