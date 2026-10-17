#!/usr/bin/env python3
"""Domain-wide HxlR motif and subfamily/neighborhood context audit (1524 proteins)."""
from pathlib import Path
import pandas as pd,json,collections
B=Path('__PROJECT_ROOT__/study13')
O=B/'13_HxlR_broad_and_cassette_20261014'
h=pd.read_csv(O/'HxlR_all_proteins_annotated.tsv',sep='\t')
n=pd.read_csv(O/'HxlR_allgenome_flanking_3CDS.tsv',sep='\t')
n=n[n.flank_offset_CDS.isin([-1,1])].copy()
lookup={}
for g,grp in h.groupby('GCA'):
 expected=set(grp.protein_id)
 name=None;seq=[]
 for line in open(B/'02_prokka/faa'/f'{g}.faa'):
  if line.startswith('>'):
   if name in expected:lookup[name]=''.join(seq)
   name=line[1:].split()[0];seq=[]
  else:seq.append(line.strip())
 if name in expected:lookup[name]=''.join(seq)
print('RECOVERED',len(lookup),len(h),flush=True)
h['length_aa']=[len(lookup.get(i,'')) for i in h.protein_id]
h['Bx_C11_K13']=[bool(len(lookup.get(i,''))>=13 and lookup[i][10]=='C' and lookup[i][12]=='K') for i in h.protein_id]
h['any_CxK_first35']=[any(s[j]=='C' and s[j+2]=='K' for j in range(min(33,len(s)-2))) for i in h.protein_id for s in [lookup.get(i,'')]]
h[['GCA','Genus','Species','HMT','protein_id','MMseqs_X5_representative','length_aa','Bx_C11_K13','any_CxK_first35']].to_csv(O/'HxlR_all1524_protein_motif_audit.tsv',sep='\t',index=False)
agg=h.groupby(['Genus','MMseqs_X5_representative']).agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),n_formaldehyde_like_C11K13=('Bx_C11_K13','sum'),n_any_CxK_first35=('any_CxK_first35','sum'),median_length=('length_aa','median')).reset_index()
agg.to_csv(O/'HxlR_subfamily_CxK_by_genus.tsv',sep='\t',index=False)
gene=(h.groupby(['Genus','Species']).agg(n_proteins=('protein_id','size'),n_CxK=('any_CxK_first35','sum')).reset_index())
gene.to_csv(O/'HxlR_CxK_species_summary.tsv',sep='\t',index=False)
kind=n.groupby(['HxlR_X5cluster','Genus','neighbor_product_group']).agg(n_neighbour_genes=('neighbor_protein','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique')).reset_index()
kind.to_csv(O/'HxlR_sequence_cluster_by_neighbor_group.tsv',sep='\t',index=False)
top=h.groupby('MMseqs_X5_representative').agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_genera=('Genus','nunique'),n_HMT=('HMT','nunique'),n_CxK=('any_CxK_first35','sum')).sort_values('n_genomes',ascending=False)
top.to_csv(O/'HxlR_top_clusters_all_genomes.tsv',sep='\t')
report={'n_HxlR':len(h),'n_protein_sequences':len(lookup),'C11K13_count':int(h.Bx_C11_K13.sum()),'any_CxK_first35_count':int(h.any_CxK_first35.sum()),'family_clusters':h.MMseqs_X5_representative.nunique()}
(O/'HxlR_CxK_motif_summary.json').write_text(json.dumps(report,indent=2))
print('RESULT',report,flush=True)
print('CXK_BY_GENUS',gene.to_string(index=False),flush=True)
