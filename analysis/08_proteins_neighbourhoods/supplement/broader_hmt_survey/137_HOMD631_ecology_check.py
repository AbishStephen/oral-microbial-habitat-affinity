#!/usr/bin/env python3
"""External ecological source concordance audit (descriptive), using
four site-abundance summaries for HMTs beyond the 47 primary panel.
This is not an independently recruited clinical validation set.
Plaque = mean of SUPP/SUBP; mucosa = mean of TDO/BMU; source cutoff
and score pseudocounts are disclosed. Retain unknown/unassigned taxa.
"""
from pathlib import Path
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
c=pd.read_csv(O/'HOMD631_HxlR_SnoaL2_oneHMT_oneGenome_census.tsv',sep='\t')
sources=['erenv1v3','erenv3v5','hmprefseqv1v3','hmprefseqv3v5']
avail=[]
for src in sources:
 f=O/f'HOMD_sites_{src}.tsv'
 t=pd.read_csv(f,sep='\t',skiprows=1,low_memory=False)
 n=pd.to_numeric(t.HMT,errors='coerce')
 t=t[n.notna()&(n>0)].copy();t['HMT_link']=['HMT-'+str(int(v)).zfill(3) for v in t.HMT]
 t=t.drop_duplicates('HMT_link')
 for x in ['SUPP_mean','SUBP_mean','TDO_mean','BMU_mean']:
  if x not in t:raise ValueError((src,x,list(t.columns)))
  t[x]=pd.to_numeric(t[x],errors='coerce').fillna(0)
 pl=(t.SUPP_mean+t.SUBP_mean)/2
 mu=(t.TDO_mean+t.BMU_mean)/2
 p=0.01
 t[src+'_score']=np.log2((pl+p)/(mu+p))
 t[src+'_usable']=(np.maximum(pl,mu)>=0.01)
 avail.append(t[['HMT_link',src+'_score',src+'_usable']])
 print(src,'hmt',len(t),'usable',int(t[src+'_usable'].sum()),flush=True)
ec=c.copy()
for t in avail:ec=ec.merge(t,on='HMT_link',how='left',validate='one_to_one')
score=ec[[src+'_score' for src in sources]].to_numpy(dtype=float)
usable=ec[[src+'_usable' for src in sources]].fillna(False).to_numpy(dtype=bool)
score[~usable]=np.nan
ec['n_usable_sources']=np.isfinite(score).sum(axis=1)
ec['median_plaque_mucosa_log2']=np.nanmedian(score,axis=1)
ec['n_plaque_support']=(score>=1).sum(axis=1)
ec['n_mucosa_support']=(score<=-1).sum(axis=1)
ec['n_near_equal']=((score<1)&(score>-1)).sum(axis=1)
ec['ecology_4source_conservative']='undetermined'
ec.loc[(ec.n_usable_sources>=3)&(ec.n_plaque_support>=3)&(ec.n_mucosa_support==0),'ecology_4source_conservative']='plaque'
ec.loc[(ec.n_usable_sources>=3)&(ec.n_mucosa_support>=3)&(ec.n_plaque_support==0),'ecology_4source_conservative']='mucosal'
ec.to_csv(O/'HOMD631_HxlR_SnoaL2_4source_ecology_annotations.tsv',sep='\t',index=False)
oral=ec[ec.hmt_primary_body_site_w_abundance=='Oral']
groups=oral.groupby('ecology_4source_conservative').agg(n_HMT=('GCA','size'),HxlR_positive=('PF01638',lambda s:int((s>0).sum())),SnoaL2_positive=('PF12680',lambda s:int((s>0).sum())),both_domains=('both','sum')).reset_index()
groups['HxlR_pct']=100*groups.HxlR_positive/groups.n_HMT
groups['SnoaL2_pct']=100*groups.SnoaL2_positive/groups.n_HMT
groups.to_csv(O/'HOMD631_HxlR_SnoaL2_conservative_ecology_by_HMT.tsv',sep='\t',index=False)
# Primary-label consistency only assesses agreement between abundance summaries, not independent corroboration.
orig=pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t')
d=ec.merge(orig[['HMT','Ecology']],left_on='HMT_link',right_on='HMT',how='inner')
conf=d[d.ecology_4source_conservative!='undetermined']
print('ORAL\n'+groups.to_string(index=False),flush=True)
print('ORIGINAL_47_COMPARISON',len(d),'strongly_classified',len(conf),'agreement',sum(conf.Ecology==conf.ecology_4source_conservative),flush=True)
print('MUCOSAL_HXLR_AND_SNOAL',oral[(oral.ecology_4source_conservative=='mucosal')&(oral.both)][['HMT_link','main_genus','main_species','PF01638','PF12680','median_plaque_mucosa_log2']].head(45).to_string(index=False),flush=True)
