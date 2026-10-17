#!/usr/bin/env python3
"""Evidence-led candidates for PF01638-like regulator targets.
Use adjacent genes in 13-genus collection; collapse to distinct HMTs within
sequence cluster; orientation-aware promoter hypothesis; preserve caveat.
"""
from pathlib import Path
import pandas as pd,numpy as np
B=Path('__PROJECT_ROOT__/study13')
O=B/'15_HxlR_regulator_target_ligand_20261014';O.mkdir(exist_ok=True)
h=pd.read_csv(B/'13_HxlR_broad_and_cassette_20261014/HxlR_allgenome_flanking_3CDS.tsv',sep='\t')
h=h[h.flank_offset_CDS.abs()==1].copy()
h['product_norm']=h.neighbor_product.fillna('').str.replace('%2C',',',regex=False)
h['candidate_divergent_promoter']=(
 ((h.flank_offset_CDS==1)&(h.HxlR_strand=='-')&(h.neighbor_strand=='+'))|
 ((h.flank_offset_CDS==-1)&(h.HxlR_strand=='+')&(h.neighbor_strand=='-'))
)
h['candidate_direct_promoter']=h['candidate_divergent_promoter'] & h.intergenic_distance_bp.between(0,500)
h['redox']=h.neighbor_product_group=='redox_metabolism'
h['is_unknown']=h.neighbor_product_group=='uncharacterized'
# one adjacent neighbour class observed per HMT/cluster/genus, instead of treating duplicate genomes as independent
summary=(h.groupby(['HxlR_X5cluster','Genus','product_norm','neighbor_product_group'])
   .agg(n_gene_instances=('HxlR_protein','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),
       n_close_divergent_instances=('candidate_direct_promoter','sum')).reset_index()
   .sort_values(['n_HMT','n_genomes'],ascending=[False,False]))
summary.to_csv(O/'HxlR_X5cluster_adjacent_target_family_evidence.tsv',sep='\t',index=False)
g=(h.groupby(['HxlR_X5cluster','Genus','neighbor_product_group'])
   .agg(n_instances=('GCA','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),
        n_close_divergent=('candidate_direct_promoter','sum')).reset_index()
   .sort_values(['n_HMT','n_genomes'],ascending=[False,False]))
g.to_csv(O/'HxlR_X5cluster_neighbor_categories.tsv',sep='\t',index=False)
print('SUMMARY_SELECTED',flush=True)
for product in ['Quinone oxidoreductase 2','Coenzyme A disulfide reductase','NADH oxidase',
                'Putative NAD(P)H nitroreductase MhqN','NADPH dehydrogenase',
                'Hydroxyacylglutathione hydrolase GloB','hypothetical protein']:
 z=summary[summary.product_norm==product].head(8)
 print('TARGET',product,flush=True)
 print(z[['HxlR_X5cluster','Genus','n_genomes','n_HMT','n_close_divergent_instances']].to_string(index=False),flush=True)
print('TOTAL',len(summary),'target combinations',flush=True)
