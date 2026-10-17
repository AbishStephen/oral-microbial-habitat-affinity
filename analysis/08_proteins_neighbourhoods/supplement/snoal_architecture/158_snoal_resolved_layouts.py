#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, collections
B=Path('__PROJECT_ROOT__/study13')
O=B/'18_SnoaL2_function_gene_integrity_20261016'
S=Path('__WORK_ROOT__/SnoaL2_subfamilies_20261016')
a=pd.read_csv(O/'SnoaL2_original524_protein_sizes.tsv',sep='\t')
clusters=pd.read_csv(B/'17_SnoaL2_independent_targets_20261016/SnoaL2_all524_protein_clusters.tsv',sep='\t')
a=a.merge(clusters[['protein_id','cluster']],on='protein_id',validate='one_to_one')
hits=[]
for pf in ['PF12680','PF07366','PF13577']:
 for line in open(S/'domains'/f'{pf}.long524.domtblout'):
  if line.startswith('#'):continue
  z=line.split()
  if len(z)<22:continue
  hits.append({'protein_id':z[0],'pfam':pf,'full_evalue':float(z[6]),'i_evalue':float(z[12]),'hmm_from':int(z[15]),'hmm_to':int(z[16]),'env_start':int(z[19]),'env_end':int(z[20]),'ali_start':int(z[17]),'ali_end':int(z[18])})
d=pd.DataFrame(hits)
d=d[d.i_evalue<=1e-3]
d.to_csv(O/'SnoaL2_524_PF07366_PF12680_PF13577_domain_segments.tsv',sep='\t',index=False)
rows=[]
for r in a.itertuples():
 sub=d[d.protein_id==r.protein_id].sort_values('ali_start')
 # group distinct regions across PF models if at least 50% of smaller span overlaps
 seg=[]
 for x in sub.itertuples():
  hitspec={'start':x.ali_start,'end':x.ali_end,'pfams':{x.pfam},'evalue':x.i_evalue}
  overlapgroups=[j for j,y in enumerate(seg) if max(0,min(y['end'],hitspec['end'])-max(y['start'],hitspec['start'])+1)>=0.5*min(y['end']-y['start']+1,hitspec['end']-hitspec['start']+1)]
  if overlapgroups:
   k=overlapgroups[0];y=seg[k];y['start']=min(y['start'],hitspec['start']);y['end']=max(y['end'],hitspec['end']);y['pfams'].add(x.pfam);y['evalue']=min(y['evalue'],x.i_evalue)
  else:seg.append(hitspec)
 seg.sort(key=lambda x:x['start'])
 row={'protein_id':r.protein_id,'Genus':r.Genus,'HMT':r.HMT,'cluster':r.cluster,'length_aa':r.protein_length_aa,
      'linked_HxlR':r.linked_HxlR_3kb,'n_independent_snoal_regions':len(seg),
      'n_original_hmm_hits':len(sub),
      'domain_layout':';'.join(f"{s['start']}-{s['end']}:{'+'.join(sorted(s['pfams']))}" for s in seg),
      'domains_cover_fraction':round(sum(x['end']-x['start']+1 for x in seg)/r.protein_length_aa,3)}
 rows.append(row)
out=pd.DataFrame(rows)
out.to_csv(O/'SnoaL2_524_resolved_independent_domain_layouts.tsv',sep='\t',index=False)
long=out[out.length_aa>=200]
print('ALL',out.n_independent_snoal_regions.value_counts().to_dict())
print('LONG',len(long),'REGIONS',long.n_independent_snoal_regions.value_counts().to_dict())
print('LONG_CLUSTERS',long.groupby('cluster').agg(n=('protein_id','size'),median_length=('length_aa','median'),median_domains=('n_independent_snoal_regions','median'),linked_HxlR=('linked_HxlR','sum')).to_string())
for cl,z in long.groupby('cluster'):
 print('GROUP',cl)
 print(z[['protein_id','Genus','length_aa','n_independent_snoal_regions','domain_layout']].head(4).to_string(index=False))
