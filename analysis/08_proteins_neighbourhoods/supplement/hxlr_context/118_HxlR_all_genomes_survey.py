#!/usr/bin/env python3
"""HxlR global survey: all 1121 nested HOMD genomes; 13 genus / 47 HMT universe.
Does not treat reference genomes as independent habitat replicates."""
from pathlib import Path
import pandas as pd, collections, json
B=Path('__PROJECT_ROOT__/study13')
D=B/'10_additional_analyses_20261012'
OUT=B/'13_HxlR_broad_and_cassette_20261014'; OUT.mkdir(exist_ok=True)
meta=pd.read_csv(B/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')
domains=pd.read_csv(D/'crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t')
prot=pd.read_csv(B/'11_manuscript_completion_20261013/candidate_Pfam_fullprotein_orthogroup_architecture.tsv',sep='\t')
h=domains[domains.Pfam=='PF01638'].copy()
hp=prot[prot.Pfam=='PF01638'].copy()
print('META',len(meta),meta.Genus.nunique(),meta.HMT.nunique(),'HxlR',len(h),'domain proteins',len(hp),flush=True)
assert h.protein_id.nunique()==len(h)==len(hp)
h=h.merge(hp[['protein_id','MMseqs_X5_representative','all_Pfam_domains','n_domain_families']],on='protein_id',how='left',validate='one_to_one')
h=h.merge(meta[['GCA','Ecology','Ecology_tier']],on='GCA',how='left',validate='many_to_one')
print('genus annotations',h.Genus.value_counts().to_dict())
h.to_csv(OUT/'HxlR_all_proteins_annotated.tsv',sep='\t',index=False)
counts=h.groupby('GCA').size().rename('HxlR_copy_count')
clusters=h.groupby('GCA').MMseqs_X5_representative.nunique().rename('HxlR_cluster_count')
m=meta.merge(counts,on='GCA',how='left').merge(clusters,on='GCA',how='left')
m[['HxlR_copy_count','HxlR_cluster_count']]=m[['HxlR_copy_count','HxlR_cluster_count']].fillna(0).astype(int)
m.to_csv(OUT/'HxlR_all1121_genome_copy_number.tsv',sep='\t',index=False)
tbl=m.groupby(['Genus','HMT','Species','Ecology']).agg(genomes=('GCA','size'),HxlR_positive_genomes=('HxlR_copy_count',lambda x:int((x>0).sum())),HxlR_mean_copies=('HxlR_copy_count','mean'),HxlR_max_copies=('HxlR_copy_count','max')).reset_index()
unique=h.groupby(['Genus','HMT']).MMseqs_X5_representative.nunique().rename('distinct_X5_clusters').reset_index()
tbl=tbl.merge(unique,on=['Genus','HMT'],how='left').fillna({'distinct_X5_clusters':0})
tbl['positive_fraction']=tbl.HxlR_positive_genomes/tbl.genomes
tbl.to_csv(OUT/'HxlR_prevalence_47HMT_13genus.tsv',sep='\t',index=False)
gen=m.groupby('Genus').agg(genomes=('GCA','size'),HxlR_positive_genomes=('HxlR_copy_count',lambda z:int((z>0).sum())),copies=('HxlR_copy_count','sum'),mean_copies=('HxlR_copy_count','mean'),max_copies=('HxlR_copy_count','max')).reset_index()
gen=gen.merge(h.groupby('Genus').MMseqs_X5_representative.nunique().rename('distinct_X5_protein_clusters'),on='Genus',how='left').fillna(0)
gen['positive_fraction']=gen.HxlR_positive_genomes/gen.genomes
gen.to_csv(OUT/'HxlR_genus_survey.tsv',sep='\t',index=False)
sub=h.groupby(['Genus','Species','HMT','MMseqs_X5_representative']).agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_arch=('all_Pfam_domains','nunique')).reset_index()
sub.to_csv(OUT/'HxlR_subfamily_distribution_by_species.tsv',sep='\t',index=False)
# HxlR target cassette annotations in all 67 genomes; partition by species name correctly.
old=pd.read_csv(B/'12_HxlR_SnoaL_regulatory_genomics_20261013/linked_locus_audit.tsv',sep='\t')
old=old.merge(meta[['GCA','Species','Ecology']],on='GCA',how='left')
old.to_csv(OUT/'HxlR_SnoaL_67_genome_locus_metadata.tsv',sep='\t',index=False)
p2=tbl[tbl.Genus.isin(['Porphyromonas','Prevotella'])]
print('SPECIES_MATRIX\n'+p2[['Genus','Species','HMT','genomes','HxlR_positive_genomes','HxlR_mean_copies','distinct_X5_clusters']].to_string(index=False),flush=True)
print('GENUS\n'+gen.to_string(index=False),flush=True)
print('CASSETTE\n'+old.groupby(['Genus','Species','HMT']).size().to_string(),flush=True)
summary={'universe_genomes':int(len(meta)),'HxlR_proteins':int(len(h)),'HxlR_positive_genomes':int((m.HxlR_copy_count>0).sum()),'HxlR_positive_genera':int((gen.HxlR_positive_genomes>0).sum()),'protein_clusters':int(h.MMseqs_X5_representative.nunique()),'scoped_genomes_porphyromonas':int(sum(m.Genus=='Porphyromonas')),'scoped_genomes_prevotella':int(sum(m.Genus=='Prevotella')),'scope':'Curated HOMD 13-genus, 47 HMT, 1121 genomes; not all public oral genomes'}
(OUT/'HxlR_global_survey_summary.json').write_text(json.dumps(summary,indent=2))
print('DONE',summary,flush=True)
