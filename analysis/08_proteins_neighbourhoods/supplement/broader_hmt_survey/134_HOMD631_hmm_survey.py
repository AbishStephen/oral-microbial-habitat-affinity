#!/usr/bin/env python3
"""HOMD 631 HMT one-reference-per-taxon independent Pfam HMM screen summary.
Report broad oral site versus nonoral, never call "oral" = plaque/mucosa.
"""
from pathlib import Path
import csv,json,pandas as pd,collections
B=Path('__PROJECT_ROOT__/study13');O=B/'14_HxlR_fullHOMD_controls_20261014'
panel=pd.read_csv(O/'HOMD_2025_FTP_taxon_representative_panel.tsv',sep='\t')
def parse(pfam):
 t=[]
 for line in open(O/f'HOMD631_{pfam}_domtblout.tsv'):
  if line.startswith('#'):continue
  z=line.split()
  if len(z)<19:continue
  sid=z[0];gca=sid.split('|')[0]
  t.append({'GCA':gca,'protein_ID':sid,'Pfam':pfam,'full_Evalue':float(z[6]),'domain_iEvalue':float(z[12]),'domain_score':float(z[13]),'HMM_start':int(z[15]),'HMM_end':int(z[16]),'protein_start':int(z[17]),'protein_end':int(z[18])})
 return pd.DataFrame(t)
h=parse('PF01638');s=parse('PF12680');both=pd.concat([h,s],ignore_index=True)
both.to_csv(O/'HOMD631_HxlR_SnoaL2_independent_HMM_protein_hits.tsv',sep='\t',index=False)
counts=both.groupby(['GCA','Pfam']).protein_ID.nunique().unstack(fill_value=0)
res=panel.merge(counts,on='GCA',how='left')
for k in ('PF01638','PF12680'):res[k]=res[k].fillna(0).astype(int)
res['both']= (res.PF01638>0)&(res.PF12680>0)
res['oral']=res.hmt_primary_body_site_w_abundance=='Oral'
res.to_csv(O/'HOMD631_HxlR_SnoaL2_oneHMT_oneGenome_census.tsv',sep='\t',index=False)
summary=res.groupby('hmt_primary_body_site_w_abundance').agg(HMTs=('GCA','size'),HxlR=('PF01638',lambda x:int((x>0).sum())),SnoaL2=('PF12680',lambda x:int((x>0).sum())),both=('both','sum')).reset_index()
summary['HxlR_pct']=100*summary.HxlR/summary.HMTs
summary['SnoaL2_pct']=100*summary.SnoaL2/summary.HMTs
summary.to_csv(O/'HOMD631_domain_prevalence_primary_body_site.tsv',sep='\t',index=False)
g=res.groupby(['main_genus','hmt_primary_body_site_w_abundance']).agg(HMTs=('GCA','size'),HxlR=('PF01638',lambda x:int((x>0).sum())),SnoaL2=('PF12680',lambda x:int((x>0).sum())),both=('both','sum')).reset_index()
g.to_csv(O/'HOMD631_domain_prevalence_genus.tsv',sep='\t',index=False)
oral=res[res.oral]
print('HMMHITS',len(h),'HxlR proteins',h.protein_ID.nunique(),'SNO',len(s),'protein',s.protein_ID.nunique(),flush=True)
print('TOTAL HMT',len(res),'ORAL HMT',len(oral),'HX ORAL',int((oral.PF01638>0).sum()),'SNO ORAL',int((oral.PF12680>0).sum()),'BOTH ORAL',int(oral.both.sum()),flush=True)
print('SITE TABLE\n'+summary.to_string(index=False),flush=True)
print('ORAL GENERA both PFAM\n'+g[(g.hmt_primary_body_site_w_abundance=='Oral')&(g.both>0)].sort_values('both',ascending=False).to_string(index=False),flush=True)
print('ORAL non-primary HMT HxlR',len(oral[(oral.PF01638>0) & (~oral.HMT_link.isin(pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t').HMT))]),flush=True)
