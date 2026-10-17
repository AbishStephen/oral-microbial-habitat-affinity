#!/usr/bin/env python3
"""Cross-lineage high-priority nominations and nonredundant profile clusters."""
from pathlib import Path
import gzip,csv
import numpy as np,pandas as pd
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
ra=pd.read_csv(O/'crossgenus_allfeatures_raw_permutation_results.tsv',sep='\t')
ad=pd.read_csv(O/'crossgenus_allfeatures_genome_adjusted_permutation_results.tsv',sep='\t')
choose=['feature','effect_genus_mean','n_plaque_genus_material','n_mucosal_genus_material','n_plaque_genus_any_delta','n_mucosal_genus_any_delta','n_phyla_plaque','n_phyla_mucosa','phylum_plaque','phylum_mucosa','perm_p_consistency','maxT_all_feature_FWER_p_consistency','perm_p_t','maxT_all_feature_FWER_p_t']
rr=ra[choose].rename(columns={k:'raw_'+k for k in choose if k!='feature'})
aa=ad[choose].rename(columns={k:'adj_'+k for k in choose if k!='feature'})
z=rr.merge(aa,on='feature',validate='one_to_one')
z['direction']=np.where((z.raw_effect_genus_mean>0)&(z.adj_effect_genus_mean>0),'plaque',np.where((z.raw_effect_genus_mean<0)&(z.adj_effect_genus_mean<0),'mucosal','discordant'))
z['genus_support']=np.where(z.direction=='plaque',np.minimum(z.raw_n_plaque_genus_material,z.adj_n_plaque_genus_material),np.minimum(z.raw_n_mucosal_genus_material,z.adj_n_mucosal_genus_material))
z['phyla_support']=np.where(z.direction=='plaque',np.minimum(z.raw_n_phyla_plaque,z.adj_n_phyla_plaque),np.minimum(z.raw_n_phyla_mucosa,z.adj_n_phyla_mucosa))
z['support_selection']=((z.direction!='discordant')&(z.genus_support>=5)&(z.phyla_support>=3)&(z.raw_perm_p_consistency<.05)&(z.adj_perm_p_consistency<.05))
pfname={}
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  pfname.setdefault(row['Pfam'],row['Pfam_name'])
cogname={}
for l in (O/'ontology/cog-24.def.tab').read_text().splitlines():
 e=l.split('\t')
 if len(e)>2:cogname[e[0]]=e[2]
kpath={}
for fname,pre in [('kegg_pathway_names.tsv','KPATH:'),('kegg_module_names.tsv','KMOD:')]:
 for l in (B/'05_function/mappings'/fname).read_text().splitlines():
  e=l.split('\t',1)
  if len(e)==2:kpath[pre+e[0]]=e[1]
def annotate(f):
 if f.startswith('PFAM:'):return pfname.get(f.split(':')[1],'PFAM unknown')
 if f.startswith('COG:'):return cogname.get(f.split(':')[1],'COG unknown')
 if f.startswith('KMOD:') or f.startswith('KPATH:'):return kpath.get(f,'')
 if f.startswith('X5:'):return 'MMseqs2 orthogroup (requires sequence-based annotation)'
 return ''
z.insert(1,'annotation',z.feature.map(annotate))
z['layer']=z.feature.str.split(':').str[0]
z=z.sort_values(['support_selection','genus_support','phyla_support','raw_perm_p_consistency','adj_perm_p_consistency'],ascending=[False,False,False,True,True])
z.to_csv(O/'crossgenus_feature_integrated_raw_adjusted_annotated.tsv',sep='\t',index=False)
good=z[z.support_selection].copy()
print('CONCORDANT NOMINATIONS',len(good),'LAYERS',good.layer.value_counts().to_dict())
print(good[['feature','annotation','direction','genus_support','phyla_support','raw_effect_genus_mean','raw_n_plaque_genus_material','raw_n_mucosal_genus_material','raw_perm_p_consistency','adj_perm_p_consistency']].head(80).to_string(index=False),flush=True)
nom=good.head(160).feature.tolist()+['PFAM:PF04072','PFAM:PF12680','PFAM:PF10604','PFAM:PF07366','PFAM:PF22696']
nom=list(dict.fromkeys(nom))
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');h=next(rd);dat=list(rd)
colidx=[h.index(f) for f in nom];X=np.array([[float(row[j]) for j in colidx] for row in dat])
corr=np.corrcoef(X,rowvar=False);corr=np.nan_to_num(corr)
distance=np.clip(1-np.abs(corr),0,2)
distance=(distance+distance.T)/2;np.fill_diagonal(distance,0)
clusters=fcluster(linkage(squareform(distance),method='complete'),t=.10,criterion='distance')
m=pd.DataFrame({'feature':nom,'cluster':clusters,'annotation':[annotate(x) for x in nom]})
stats=z.set_index('feature')
m['direction']=[stats.loc[f,'direction'] for f in nom]
m['genus_support']=[int(stats.loc[f,'genus_support']) for f in nom]
m['max_fwer_raw']=[float(stats.loc[f,'raw_maxT_all_feature_FWER_p_consistency']) for f in nom]
m['max_fwer_adj']=[float(stats.loc[f,'adj_maxT_all_feature_FWER_p_consistency']) for f in nom]
m.sort_values(['cluster','genus_support'],ascending=[True,False]).to_csv(O/'crossgenus_candidate_correlation_modules.tsv',sep='\t',index=False)
print('CLUSTERS',len(set(clusters)),'FOR',len(nom),'FEATURES',flush=True)
print('LARGEST',m.groupby('cluster').size().sort_values(ascending=False).head(16).to_string(),flush=True)
