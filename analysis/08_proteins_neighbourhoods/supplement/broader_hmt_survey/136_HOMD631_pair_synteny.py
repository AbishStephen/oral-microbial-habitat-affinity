#!/usr/bin/env python3
"""Physical HxlR/SnoaL_2 linkage and +/-3 neighboring gene architecture in
HOMD 631 taxon-representative genomes with both Pfam HMM hits.
Species/HMT-level replication; note primary site Oral != plaque or mucosa.
"""
from pathlib import Path
import csv,collections,json
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
G=Path('__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014/GFF_both')
dom=pd.read_csv(O/'HOMD631_HxlR_SnoaL2_independent_HMM_protein_hits.tsv',sep='\t')
census=pd.read_csv(O/'HOMD631_HxlR_SnoaL2_oneHMT_oneGenome_census.tsv',sep='\t').set_index('GCA')
hits={g:{'PF01638':set(),'PF12680':set()} for g in census.index}
for d in dom.itertuples():
 if d.GCA in hits:hits[d.GCA][d.Pfam].add(d.protein_ID.split('|',1)[1])
rows=[];neighborhood=[];missing=[];hit_missed=[]
for n,(g,hit) in enumerate(hits.items(),1):
 if not hit['PF01638'] or not hit['PF12680']:continue
 p=G/f'{g}.gff'
 if not p.exists():missing.append(g);continue
 records=collections.defaultdict(list);loc={}
 with p.open() as f:
  for line in f:
   if line.startswith('##FASTA'):break
   if line.startswith('#'):continue
   z=line.rstrip('\n').split('\t')
   if len(z)!=9 or z[2]!='CDS':continue
   attrs=dict(a.split('=',1) for a in z[8].split(';') if '=' in a)
   pid=attrs.get('ID',attrs.get('locus_tag',''))
   if not pid:continue
   v={'id':pid,'contig':z[0],'start':int(z[3]),'end':int(z[4]),'strand':z[6],'product':attrs.get('product','').replace('%20',' ')}
   records[z[0]].append(v);loc[pid]=v
 for ct in records:records[ct].sort(key=lambda r:r['start'])
 pos={p['id']:(ct,i) for ct,genes in records.items() for i,p in enumerate(genes)}
 hr=[loc[p] for p in hit['PF01638'] if p in loc]
 sn=[loc[p] for p in hit['PF12680'] if p in loc]
 bad=len(hit['PF01638'])-len(hr)+len(hit['PF12680'])-len(sn)
 if bad:hit_missed.append({'GCA':g,'hits_not_in_GFF':bad})
 for a in hr:
  for b in sn:
   if a['contig']!=b['contig'] or a['id']==b['id']:continue
   left,right=sorted([a,b],key=lambda d:d['start'])
   d=max(0,right['start']-left['end']-1)
   if d>3000:continue
   orient=left['strand']+right['strand']
   met=census.loc[g]
   pair={'GCA':g,'HMT':met.HMT_link,'Genus':met.main_genus,'Species':met.main_species,
    'Primary_site':met.hmt_primary_body_site_w_abundance,'HxlR_protein':a['id'],'SnoaL2_protein':b['id'],
    'genomic_contig':a['contig'],'HxlR_start':a['start'],'SnoaL2_start':b['start'],
    'distance_bp':int(d),'left_gene_domain':'HxlR' if left['id']==a['id'] else 'SnoaL_2',
    'orientation':orient,'tight_divergent':orient=='-+' and d<=500}
   rows.append(pair);ct=a['contig'];genes=records[ct]
   s=min(pos[a['id']][1],pos[b['id']][1]);e=max(pos[a['id']][1],pos[b['id']][1])
   for side,kvals in [('left',range(max(0,s-3),s)),('right',range(e+1,min(len(genes),e+4)))]:
    for k in kvals:
     gene=genes[k]
     dd=(left['start']-gene['end']-1) if side=='left' else (gene['start']-right['end']-1)
     neighbor={**{i:pair[i] for i in ['GCA','HMT','Genus','Species','Primary_site']},
       'HxlR_protein':a['id'],'SnoaL2_protein':b['id'],'side':side,'offset_CDS':k-s if side=='left' else k-e,
       'distance_bp_to_core':int(dd),'gene_id':gene['id'],'product':gene['product'],'strand':gene['strand']}
     neighborhood.append(neighbor)
 if n%40==0:print('PROGRESS',n,'/',len(hits),flush=True)
r=pd.DataFrame(rows);nn=pd.DataFrame(neighborhood)
r.to_csv(O/'HOMD631_PHYSICAL_HxlR_SnoaL2_loci_within3kb.tsv',sep='\t',index=False)
nn.to_csv(O/'HOMD631_PHYSICAL_HxlR_SnoaL2_flanks3CDS.tsv',sep='\t',index=False)
x=census.reset_index()
x['physical_pair']=x.GCA.isin(r.GCA) if len(r) else False
x['tight_divergent_pair']=x.GCA.isin(r[r.tight_divergent].GCA) if len(r) else False
group=x.groupby('hmt_primary_body_site_w_abundance').agg(n_HMT=('GCA','size'),n_cooccur=('both','sum'),n_physical_3kb=('physical_pair','sum'),n_tight_divergent=('tight_divergent_pair','sum')).reset_index()
group.to_csv(O/'HOMD631_PHYSICAL_HxlR_SnoaL2_site_summary.tsv',sep='\t',index=False)
rby=r.groupby(['Primary_site','Genus']).agg(HMTs_with_pair=('HMT','nunique'),physical_gene_pairs=('GCA','size'),tight_divergent=('tight_divergent','sum')).reset_index() if len(r) else pd.DataFrame()
rby.to_csv(O/'HOMD631_PHYSICAL_HxlR_SnoaL2_genus_summary.tsv',sep='\t',index=False)
(O/'HOMD631_PHYSICAL_audit_summary.json').write_text(json.dumps({'representative_HMT':len(x),'both_domains':int(x.both.sum()),'physical_3kb_HMT':int(x.physical_pair.sum()),'tight_divergent_HMT':int(x.tight_divergent_pair.sum()),'pair_instances':len(r),'unretrieved_GFF':len(missing),'unmatched_proteins':hit_missed},indent=2))
print('TOTAL',len(x),'BOTH',int(x.both.sum()),'PAIR',len(r),'PAIRHMT',int(x.physical_pair.sum()),'TIGHT',int(x.tight_divergent_pair.sum()),'MISSING',missing,flush=True)
print('SITE\n'+group.to_string(index=False),flush=True)
print('ORAL PAIRS\n'+rby[rby.Primary_site=='Oral'].sort_values('HMTs_with_pair',ascending=False).to_string(index=False),flush=True)
