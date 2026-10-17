#!/usr/bin/env python3
"""Symmetric physical-neighborhood audit around both flanks of every linked HxlR/SnoaL2 pair.
Preserve CDS orientation and actual gene order; compare by HMT/ecology.
"""
from pathlib import Path
import pandas as pd,collections,csv,json
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014';O.mkdir(exist_ok=True)
P=B/'13_HxlR_broad_and_cassette_20261014'
l=pd.read_csv(P/'HxlR_SnoaL2_all_gene_pairs_3kb.tsv',sep='\t')
meta=pd.read_csv(B/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')
l=l.merge(meta[['GCA','Ecology','Ecology_tier']],on='GCA',how='left')
def categorize(v):
 x=v.lower()
 if any(y in x for y in ['transposase','integrase','recombinase','insertion sequence','resolvase']):return 'mobile'
 if any(y in x for y in ['reductase','oxidase','dehydrogenase','nitroreductase','oxidoreductase']):return 'redox'
 if any(y in x for y in ['transcription','regulator','repressor']):return 'regulator'
 if 'hypothetical' in x or 'uncharacterized' in x:return 'unknown'
 if any(y in x for y in ['transporter','permease','pump','export','import']):return 'transporter'
 return 'other'
flank=[]
for g,sub in l.groupby('GCA'):
 contig=collections.defaultdict(list)
 for line in open(B/'02_prokka/gff'/f'{g}.gff'):
  if line.startswith('#'):continue
  v=line.rstrip('\n').split('\t')
  if len(v)<9 or v[2]!='CDS':continue
  a=dict(t.split('=',1) for t in v[8].split(';') if '=' in t)
  contig[v[0]].append({'gene':a.get('ID',''),'start':int(v[3]),'end':int(v[4]),'strand':v[6],'product':a.get('product','').replace('%20',' ')})
 for c in contig:contig[c].sort(key=lambda z:z['start'])
 for k,row in sub.iterrows():
  genes=contig[row.contig];ind={gene['gene']:i for i,gene in enumerate(genes)}
  if row.HxlR_protein not in ind or row.SnoaL2_protein not in ind:continue
  i=ind[row.HxlR_protein];j=ind[row.SnoaL2_protein];left=min(i,j);right=max(i,j)
  for side,choices in [('left',[left-3,left-2,left-1]),('right',[right+1,right+2,right+3])]:
   for z in choices:
    if z<0 or z>=len(genes):continue
    d=genes[z];offset=z-left if side=='left' else z-right
    flank.append({'GCA':g,'HMT':row.HMT,'Genus':row.Genus,'Species':row.Species,'Ecology':row.Ecology,
      'pair_ordinal':int(k),'side':side,'position_from_core':abs(offset),'gene':d['gene'],'strand':d['strand'],
      'product':d['product'],'function_group':categorize(d['product']),
      'distance_to_core_bp':(genes[left]['start']-d['end']-1) if side=='left' else (d['start']-genes[right]['end']-1)})
o=pd.DataFrame(flank)
o.to_csv(O/'HxlR_SnoaL_all70_pair_flanking_six_CDS.tsv',sep='\t',index=False)
p=o.groupby(['Genus','HMT','Ecology','side','position_from_core','function_group']).agg(n_gene_instances=('GCA','size'),n_genomes=('GCA','nunique')).reset_index()
p.to_csv(O/'HxlR_SnoaL_flank_functions_by_HMT.tsv',sep='\t',index=False)
hm=l.groupby(['Genus','Species','HMT','Ecology']).agg(genomes_with_pair=('GCA','nunique'),pair_instances=('GCA','size'),
  tight_divergent_pairs=('tight_divergent','sum'),mean_gap_bp=('distance_bp','mean')).reset_index()
hm.to_csv(O/'HxlR_SnoaL_HMT_ecology_summary.tsv',sep='\t',index=False)
mobile=o[o.function_group=='mobile'].groupby(['GCA','pair_ordinal']).agg(n_mobile_within_3_CDS=('gene','size')).reset_index()
pairmob=l.merge(mobile,left_on=['GCA',l.index],right_on=['GCA','pair_ordinal'],how='left') if False else None
print('HMT SUMMARY\n'+hm.to_string(index=False),flush=True)
for (genus,HMT,ecology),sub in o[o.position_from_core<=2].groupby(['Genus','HMT','Ecology']):
 if genus in ('Porphyromonas','Prevotella','Streptococcus','Actinomyces'):
  print('CONTEXT',genus,HMT,ecology,'n_genomes',sub.GCA.nunique(),flush=True)
  for (side,pos),s in sub.groupby(['side','position_from_core']):
   print(' ',side,pos,[(a,int(b)) for a,b in s['product'].value_counts().head(4).items()],flush=True)
print('MOBILE',o[o.function_group=='mobile'].groupby(['Genus','HMT','Ecology','side']).GCA.nunique().to_string(),flush=True)
