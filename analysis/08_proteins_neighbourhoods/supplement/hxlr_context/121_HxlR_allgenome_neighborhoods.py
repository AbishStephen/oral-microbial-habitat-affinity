#!/usr/bin/env python3
"""HxlR all-genome neighbor audit: +/- 3 CDS around 1524 PF01638 proteins."""
from pathlib import Path
from collections import defaultdict,Counter
import csv, pandas as pd,json,sys,time
B=Path('__PROJECT_ROOT__/study13'); O=B/'13_HxlR_broad_and_cassette_20261014'
D=B/'10_additional_analyses_20261012'
h=pd.read_csv(O/'HxlR_all_proteins_annotated.tsv',sep='\t')
domain=pd.read_csv(D/'crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t')
pfams=domain.groupby(['GCA','protein_id']).Pfam.apply(lambda x:';'.join(sorted(set(x)))).to_dict()
targets={g:subset for g,subset in h.groupby('GCA')}
out=[];miss=[]
def classify(prod):
 x=prod.lower()
 for cat,keys in [('mobile_DNA',['integrase','recombinase','transposase','insertion sequence','resolvase']),('transport',['transport','permease','efflux','porin','channel','pump']),('DNA_regulation',['transcription','regulator','dna-binding','repressor','activator','sigma factor']),('redox_metabolism',['oxidoreductase','dehydrogenase','reductase','oxidase','flavoprotein','redox']),('glycan_cell_envelope',['glycosyl','polysaccharide','lipopolysaccharide','glycan','peptidoglycan','glucosidase','glycosidase']),('defense_toxin',['antitoxin','toxin','nuclease','restriction','crispr','phage']),('central_metabolism',['kinase','transferase','synthase','isomerase','hydrolase','lyase'])]:
  if any(y in x for y in keys):return cat
 if 'hypothetical' in x or 'uncharacterized' in x:return 'uncharacterized'
 return 'other_annotated'
for i,(g,subset) in enumerate(targets.items(),1):
 p=B/'02_prokka/gff'/f'{g}.gff'
 names=set(subset.protein_id)
 coords={};bycont=defaultdict(list)
 with open(p) as f:
  for line in f:
   if line.startswith('#'):continue
   z=line.rstrip('\n').split('\t')
   if len(z)<9 or z[2]!='CDS':continue
   a=dict(c.split('=',1) for c in z[8].split(';') if '=' in c)
   ident=a.get('ID','')
   if not ident:continue
   v={'gene_id':ident,'GCA':g,'contig':z[0],'start':int(z[3]),'end':int(z[4]),'strand':z[6],'product':a.get('product','').replace('%20',' ')}
   coords[ident]=v;bycont[z[0]].append(v)
 for cc,genes in bycont.items():
  genes.sort(key=lambda d:d['start'])
  for j,z in enumerate(genes):z['gene_index']=j
 for hr in subset.itertuples():
  target=coords.get(hr.protein_id)
  if not target:
   miss.append((g,hr.protein_id));continue
  gg=bycont[target['contig']]; j=target['gene_index']
  for neighbor_offset in (-3,-2,-1,1,2,3):
   ix=j+neighbor_offset
   if ix<0 or ix>=len(gg):continue
   z=gg[ix]
   distance=target['start']-z['end']-1 if neighbor_offset<0 else z['start']-target['end']-1
   orientation='divergent' if (neighbor_offset<0 and z['strand']=='-' and target['strand']=='+') or (neighbor_offset>0 and target['strand']=='-' and z['strand']=='+') else ('convergent' if (neighbor_offset<0 and z['strand']=='+' and target['strand']=='-') or (neighbor_offset>0 and target['strand']=='+' and z['strand']=='-') else 'codirectional')
   out.append({'GCA':g,'Genus':hr.Genus,'Species':hr.Species,'HMT':hr.HMT,'HxlR_protein':hr.protein_id,'HxlR_X5cluster':hr.MMseqs_X5_representative,'HxlR_strand':target['strand'],
         'flank_offset_CDS':neighbor_offset,'neighbor_protein':z['gene_id'],'neighbor_strand':z['strand'],
         'intergenic_distance_bp':distance,'relative_orientation':orientation,
         'neighbor_product':z['product'],'neighbor_product_group':classify(z['product']),
         'neighbor_Pfam_in_candidate_catalog':pfams.get((g,z['gene_id']),'')})
 if i%100==0:print('PROCESSED',i,'/',len(targets),flush=True)
result=pd.DataFrame(out)
result.to_csv(O/'HxlR_allgenome_flanking_3CDS.tsv',sep='\t',index=False)
immediate=result[result.flank_offset_CDS.isin([-1,1])].copy()
immediate['neighbor_group']=immediate.neighbor_product_group
sumtab=immediate.groupby(['Genus','neighbor_group']).agg(n_HxlR_flanks=('HxlR_protein','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique')).reset_index()
sumtab.to_csv(O/'HxlR_immediate_neighbor_functional_patterns_by_genus.tsv',sep='\t',index=False)
group=immediate.groupby(['Genus','neighbor_product']).agg(n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),n_HxlR_proteins=('HxlR_protein','nunique')).reset_index()
group.to_csv(O/'HxlR_immediate_neighbor_products_by_genus.tsv',sep='\t',index=False)
summary={'HxlR_proteins_total':int(len(h)),'HxlR_proteins_mapped':int(result.HxlR_protein.nunique()),'GFF_targets_missing':len(miss),'neighbor_rows':int(len(out))}
(O/'HxlR_neighborhood_summary.json').write_text(json.dumps(summary,indent=2))
print('SUMMARY',summary,'IMMEDIATE_CLASSES',immediate.neighbor_group.value_counts().head(12).to_dict(),flush=True)
