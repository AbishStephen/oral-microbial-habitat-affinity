#!/usr/bin/env python3
"""Interpret MMseqs2 DNA and protein homology with strain/HMT provenance and copy counts."""
from pathlib import Path
import pandas as pd, json, itertools
B=Path('__PROJECT_ROOT__/study13')
O=B/'13_HxlR_broad_and_cassette_20261014'
h=pd.read_csv(O/'HxlR_all_proteins_annotated.tsv',sep='\t')
h=h[h.Genus.isin(['Porphyromonas','Prevotella'])]
p=h.set_index('protein_id')
cols=['query','target','fident','alnlen','qlen','tlen','qcov','tcov','evalue','bits']
aln=pd.read_csv(O/'P_Prev_HxlR_protein_all_vs_all.tsv',sep='\t',names=cols)
aln=aln[(aln.qcov>=0.8)&(aln.tcov>=0.8)&(aln.evalue<1e-5)]
aln=aln[aln['query']!=aln['target']].copy()
for typ,col in [('query','query'),('target','target')]:
 for f in ['GCA','Genus','Species','HMT','MMseqs_X5_representative']:
  aln[f'{typ}_{f}']=aln[col].map(p[f])
aln['same_GCA']=aln.query_GCA==aln.target_GCA
aln['same_species']=aln.query_Species==aln.target_Species
aln.to_csv(O/'P_Prev_HxlR_all_protein_homology_annotated.tsv',sep='\t',index=False)
# For each protein choose strongest hit in each different species group; summarize directional median and best identity.
best=aln[(aln.query_Species!=aln.target_Species)].sort_values(['query','target_Genus','target_Species','fident'],ascending=[True,True,True,False]).drop_duplicates(['query','target_Genus','target_Species'])
pairs=best.groupby(['query_Genus','query_Species','target_Genus','target_Species']).agg(n_query_proteins_with_homolog=('query','nunique'),max_identity=('fident','max'),median_best_identity=('fident','median'),mean_best_identity=('fident','mean'),n_identical_aa=('fident',lambda z:int((z==1).sum()))).reset_index()
pairs.to_csv(O/'P_Prev_HxlR_cross_species_best_identity.tsv',sep='\t',index=False)
spec=h.groupby(['Genus','Species','HMT']).agg(n_HxlR_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_clusters=('MMseqs_X5_representative','nunique')).reset_index()
spec.to_csv(O/'P_Prev_HxlR_protein_subfamily_species_summary.tsv',sep='\t',index=False)
cluster=h.groupby(['MMseqs_X5_representative','Genus','Species']).agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique')).reset_index()
cluster.to_csv(O/'P_Prev_HxlR_cluster_lineage_map.tsv',sep='\t',index=False)
shared=cluster.groupby('MMseqs_X5_representative').Genus.nunique()
for c in shared[shared>=2].index:
 print('SHARED_CLUSTER',c,cluster[cluster.MMseqs_X5_representative==c].to_dict('records'),flush=True)
print('HXLRSPECIES\n'+spec.to_string(index=False),flush=True)
print('BEST_CROSS_GENUS\n'+pairs[(pairs.query_Genus=='Porphyromonas')&(pairs.target_Genus=='Prevotella')].sort_values('max_identity',ascending=False).to_string(index=False),flush=True)
nuccols=['query','target','fident','alnlen','qlen','tlen','qstart','qend','tstart','tend','evalue','bits']
nuc=pd.read_csv(O/'P_Prev_all243_nucleotide_local_homology.tsv',sep='\t',names=nuccols)
nuc['GCA']=nuc.target.str.split('__',regex=False).str[0]
nuc=nuc.merge(pd.read_csv(B/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')[['GCA','Genus','Species','HMT']],on='GCA',how='left')
nuc['query_cov']=nuc.alnlen/nuc.qlen
nuc.to_csv(O/'P_Prev_cassette_nucleotide_alignments_annotated.tsv',sep='\t',index=False)
nucstrong=nuc[(nuc.alnlen>=1000)&(nuc.fident>=0.9)]
print('NUCLEOTIDE_STRONG\n'+nucstrong[['GCA','Genus','Species','fident','alnlen','qstart','qend','tstart','tend']].to_string(index=False),flush=True)
print('NUC_COPY_COUNTS',nucstrong.GCA.value_counts().to_dict(),flush=True)
