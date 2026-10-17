#!/usr/bin/env python3
"""Genus and genome matched random-neighborhood negative-control test.
Preserves per-genome HxlR count and protein length distribution.
Tests redox and mobility enrichment of immediate +/-1 annotated neighbors.
Not a transcription assay or ecological enrichment test.
"""
import csv,json,collections,random
from pathlib import Path
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
H=pd.read_csv(B/'13_HxlR_broad_and_cassette_20261014/HxlR_all_proteins_annotated.tsv',sep='\t')
def classify(prod):
 x=prod.lower()
 return {'redox':any(y in x for y in ['reductase','oxidase','dehydrogenase','oxidoreductase','flavoprotein','nitroreductase']),
         'mobile':any(y in x for y in ['integrase','recombinase','transposase','insertion sequence','resolvase']),
         'unknown':any(y in x for y in ['hypothetical','uncharacterized'])}
def getcds(gca):
 co=collections.defaultdict(list)
 for line in open(B/'02_prokka/gff'/f'{gca}.gff'):
  if line.startswith('#'):continue
  p=line.rstrip('\n').split('\t')
  if len(p)<9 or p[2]!='CDS':continue
  attrs=dict(t.split('=',1) for t in p[8].split(';') if '=' in t)
  i=attrs.get('ID','')
  if i:co[p[0]].append({'id':i,'len':int(p[4])-int(p[3])+1,'prod':attrs.get('product','')})
 for ct in co:co[ct].sort(key=lambda x:x['id']) # GFF CDS typically sorted by coordinate; GFF IDs may track coordinate; fix below
 return co
# Explicitly parse genomic coordinates not identifier order
def parse(gca):
 co=collections.defaultdict(list)
 for line in open(B/'02_prokka/gff'/f'{gca}.gff'):
  if line.startswith('#'):continue
  p=line.rstrip('\n').split('\t')
  if len(p)<9 or p[2]!='CDS':continue
  attrs=dict(t.split('=',1) for t in p[8].split(';') if '=' in t)
  i=attrs.get('ID','')
  if i:co[p[0]].append((int(p[3]),int(p[4]),i,attrs.get('product','')))
 targets=[]
 forcont=[]
 for ct,genes in co.items():
  genes.sort(key=lambda x:x[0]);ids={p[2]:k for k,p in enumerate(genes)}
  forcont.append((ct,genes,ids))
 return forcont
lookup={g:set(s.protein_id) for g,s in H.groupby('GCA')}
groups=collections.defaultdict(list)
skip=0
for k,(g,ids) in enumerate(lookup.items(),1):
 geno=H.loc[H.GCA==g,'Genus'].iloc[0]
 for ct,genes,index in parse(g):
  blocked={j for h in ids if h in index for j in (index[h]-2,index[h]-1,index[h],index[h]+1,index[h]+2)}
  positions=[j for j in range(1,len(genes)-1) if j not in blocked]
  if not positions:continue
  contexts=np.array([[int(classify(genes[j-1][3])[z])+int(classify(genes[j+1][3])[z]) for z in ('redox','mobile','unknown')] for j in range(1,len(genes)-1)])
  # j positions constrained, map index j-1
  for protein in ids:
   if protein not in index:continue
   j=index[protein]
   if j<=0 or j>=len(genes)-1:continue
   leng=genes[j][1]-genes[j][0]+1
   match=[z for z in positions if 0.6*leng <=genes[z][1]-genes[z][0]+1<=1.4*leng]
   if len(match)<5:match=positions
   if not match:continue
   observed=np.array([int(classify(genes[j-1][3])[z])+int(classify(genes[j+1][3])[z]) for z in ('redox','mobile','unknown')])
   cands=np.array([[int(classify(genes[z-1][3])[q])+int(classify(genes[z+1][3])[q]) for q in ('redox','mobile','unknown')] for z in match],dtype=np.int8)
   groups[geno].append((observed,cands))
 if k%100==0:print('PARSED',k,len(lookup),flush=True)
rng=np.random.default_rng(20261014)
records=[]
for genus,loc in sorted(groups.items()):
 observed=np.sum([v[0] for v in loc],axis=0)
 simulated=np.zeros((500,3))
 for actual,cands in loc:
  ix=rng.integers(0,len(cands),size=500)
  simulated+=cands[ix]
 for j,label in enumerate(['redox','mobile','unknown']):
  x=float(observed[j]);null=simulated[:,j];p=(np.sum(null>=x)+1)/501
  records.append({'Genus':genus,'type':label,'n_target_HxlR_proteins_matched':len(loc),'n_flank_sites':2*len(loc),
    'observed_adjacent_genes':int(x),'obs_fraction':x/(2*len(loc)),'null_mean':float(np.mean(null)),
    'null_fraction':float(np.mean(null))/(2*len(loc)),'enrichment_ratio':float(x/max(np.mean(null),.001)),
    'one_sided_p_overrepresentation':float(p),'null_95pct_upper_count':float(np.quantile(null,.95))})
r=pd.DataFrame(records);r.to_csv(O/'HxlR_neighbor_genome_length_matched_null500.tsv',sep='\t',index=False)
print(r[r.type=='redox'].to_string(index=False),flush=True)
print('NOTE: null controls address generic annotation-neighborhood enrichment, not biochemical function or regulatory binding.',flush=True)
