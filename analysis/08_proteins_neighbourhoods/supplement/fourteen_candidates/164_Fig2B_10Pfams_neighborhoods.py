#!/usr/bin/env python3
"""Neighborhood context for all ten Figure 2B Pfam candidates.
GFF CDS-based nearest genes; store detailed rows in scratch, compact summaries under /data.
No target regulation asserted solely from gene adjacency.
"""
from pathlib import Path
import pandas as pd,collections,urllib.parse,math
B=Path('__PROJECT_ROOT__/study13')
S=Path('__WORK_ROOT__/Figure2B_14factor_20261017');S.mkdir(exist_ok=True)
O=B/'19_Figure2B_14factor_mechanisms_20261017';O.mkdir(exist_ok=True)
coord=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t')
pfs=['PF01638','PF04072','PF13577','PF13676','PF22166','PF01745','PF01337','PF13413','PF18075','PF14470']
coord=coord[coord.Pfam.isin(pfs)].drop_duplicates(['GCA','protein_id','Pfam'])
cols=['GCA','HMT','Genus','Species','protein_id','Pfam','Pfam_name','contig','start','end','strand']
coord=coord[cols]
out=[]
missing=collections.Counter()
for j,(g,z) in enumerate(coord.groupby('GCA'),1):
 fn=B/'02_prokka/gff'/f'{g}.gff'
 if not fn.exists():missing['missing_gff']+=len(z);continue
 gff=collections.defaultdict(list)
 for line in fn.open():
  if line.startswith('##FASTA'):break
  if line.startswith('#'):continue
  x=line.rstrip().split('\t')
  if len(x)<9 or x[2]!='CDS':continue
  name=x[0]
  try:start,end=int(x[3]),int(x[4])
  except:continue
  attr={}
  for a in x[8].split(';'):
   if '=' in a:
    k,v=a.split('=',1)
    if k in ('ID','product','gene','Name'):attr[k]=urllib.parse.unquote(v)
  gff[name].append(dict(id=attr.get('ID',''),product=attr.get('product',''),gene=attr.get('gene',''),start=start,end=end,strand=x[6]))
 index={}
 for ct,genes in gff.items():
  genes.sort(key=lambda v:v['start'])
  for idx,gene in enumerate(genes):index[gene['id']]=(ct,idx)
 for r in z.itertuples():
  if r.protein_id not in index:
   missing['no_protein_in_gff']+=1;continue
  ct,i=index[r.protein_id];genes=gff[ct];source=genes[i]
  for offset in (-3,-2,-1,1,2,3):
   k=i+offset
   if k<0 or k>=len(genes):continue
   target=genes[k]
   left,right=(target,source) if target['start']<source['start'] else (source,target)
   gap=right['start']-left['end']-1
   # divergence means transcription directions away from centre: left '-' and right '+'
   divergent=(left['strand']=='-' and right['strand']=='+')
   convergent=(left['strand']=='+' and right['strand']=='-')
   out.append(dict(GCA=g,HMT=r.HMT,Genus=r.Genus,protein_id=r.protein_id,Pfam=r.Pfam,
        neighbor_id=target['id'],neighbor_product=target['product'],neighbor_gene=target['gene'],
        CDS_offset=offset,intergenic_bp=gap,divergent=divergent,convergent=convergent,
        same_strand=target['strand']==source['strand'],target_aa_length=(target['end']-target['start']+1)/3,
        neighbor_coord=f"{ct}:{target['start']}-{target['end']}:{target['strand']}"))
 if j%150==0:print('SCANNED GENOMES',j,'/',coord.GCA.nunique(),flush=True)
detail=pd.DataFrame(out)
detail.to_csv(S/'Fig2B_10Pfams_all_3CDS_context.tsv',sep='\t',index=False)
adj=detail[detail.CDS_offset.abs()==1].copy()
# classification by product - categories for hypothesis generation
def cat(v):
 v=str(v).lower()
 if any(s in v for s in ['hypothetical','uncharacterized']):return 'hypothetical'
 if any(s in v for s in ['reductase','oxidase','dehydrogenase','oxygenase','peroxidase','quinone']):return 'redox enzyme'
 if any(s in v for s in ['transcription','regulator','repressor','activator','hth','sigma']):return 'regulator'
 if any(s in v for s in ['transport','permease','abc','efflux','channel']):return 'transporter'
 if any(s in v for s in ['transferase','kinase','methylase']):return 'transferase/kinase'
 if any(s in v for s in ['peptidase','protease','hydrolase','nuclease']):return 'hydrolase'
 return 'other annotated'
adj['neighbor_category']=adj.neighbor_product.map(cat)
adj.to_csv(O/'Fig2B_10Pfams_immediately_adjacent_CDS.tsv',sep='\t',index=False)
summary=adj.groupby(['Pfam','neighbor_category']).agg(n_gene_pairs=('neighbor_id','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),n_genera=('Genus','nunique')).reset_index()
summary.to_csv(O/'Fig2B_10Pfams_neighborhood_category_counts.tsv',sep='\t',index=False)
sno=adj[adj.Pfam=='PF13577'].copy()
sno.to_csv(O/'SnoaL4_PF13577_immediately_adjacent_genes.tsv',sep='\t',index=False)
print('ROWS',len(detail),'ADJ',len(adj),'SNOAL4',len(sno),'MISSING',dict(missing),flush=True)
print('SNOAL4_NEIGHBORS\n',sno.groupby('neighbor_product').agg(n=('neighbor_id','size'),HMT=('HMT','nunique'),genera=('Genus','nunique')).sort_values('n',ascending=False).head(20).to_string(),flush=True)
