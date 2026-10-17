#!/usr/bin/env python3
"""Recompute all physical PF01638 / PF12680 gene-pair loci, allowing >1 per genome.
Fixes original per-genome nearest-pair summarization and reports copy-number heterogeneity."""
from pathlib import Path
import pandas as pd, json, collections
B=Path('__PROJECT_ROOT__/study13'); O=B/'13_HxlR_broad_and_cassette_20261014'
d=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t')
h=d[d.Pfam=='PF01638'].drop_duplicates(['GCA','protein_id'])
s=d[d.Pfam=='PF12680'].drop_duplicates(['GCA','protein_id'])
out=[]
for g,he in h.groupby('GCA'):
 ss=s[s.GCA==g]
 for a in he.itertuples():
  for b in ss.itertuples():
   if a.protein_id==b.protein_id or a.contig!=b.contig:continue
   start=max(a.start,b.start);end=min(a.end,b.end)
   dist=0 if start<=end else start-end-1
   if dist>3000:continue
   left,right=(a,b) if a.start<=b.start else (b,a)
   orient=left.strand+right.strand
   out.append({'GCA':g,'Genus':a.Genus,'HMT':a.HMT,'Species':a.Species,
      'HxlR_protein':a.protein_id,'SnoaL2_protein':b.protein_id,'contig':a.contig,
      'distance_bp':dist,'orientation':orient,'divergent':orient=='-+',
      'tight_divergent':orient=='-+' and dist<=500,'HxlR_start':a.start,'HxlR_end':a.end,
      'SnoaL_start':b.start,'SnoaL_end':b.end})
out=pd.DataFrame(out).drop_duplicates(['GCA','HxlR_protein','SnoaL2_protein'])
out.to_csv(O/'HxlR_SnoaL2_all_gene_pairs_3kb.tsv',sep='\t',index=False)
summary=out.groupby(['Genus','Species','HMT']).agg(n_genomes=('GCA','nunique'),n_distinct_gene_pairs=('HxlR_protein','size'),n_tight_divergent=('tight_divergent','sum'),max_pairs_per_genome=('GCA',lambda x:int(x.value_counts().max()))).reset_index()
summary.to_csv(O/'HxlR_SnoaL2_allpairs_genus_HMT_summary.tsv',sep='\t',index=False)
print('PAIRS',len(out),'GENOMES',out.GCA.nunique(),'TIGHT_PAIRS',int(out.tight_divergent.sum()),flush=True)
print(summary.to_string(index=False),flush=True)
print('MULTI',out.GCA.value_counts()[out.GCA.value_counts()>1].to_dict(),flush=True)
