#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np,json
B=Path('__PROJECT_ROOT__/study13')
O=B/'17_SnoaL2_independent_targets_20261016'
n=pd.read_csv(B/'13_HxlR_broad_and_cassette_20261014/HxlR_allgenome_flanking_3CDS.tsv',sep='\t')
n=n[n.flank_offset_CDS.abs()==1]; n=n.merge(pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='	')[['HMT','Ecology']].drop_duplicates(),on='HMT',how='left')
# These are anchor specific sequence subfamilies, not all enzymes of the same type
modules={
'Actinomyces quinone oxidoreductase':('GCA0028475551_01508','Quinone oxidoreductase 2'),
'Fusobacterium CoA disulfide reductase':('GCA9380307751_00100','Coenzyme A disulfide reductase'),
'Fusobacterium NADH oxidase':('GCA0378938651_00779','NADH oxidase'),
'Streptococcus MhqN-like nitroreductase':('GCA_000339975.1_00196','Putative NAD(P)H nitroreductase MhqN')}
hh=pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t')
res=[];detail=[]
for mod,(c,product) in modules.items():
 z=n[(n.HxlR_X5cluster==c)&(n.neighbor_product==product)].copy()
 # 0-500bp divergent orientations point outward
 z['divergent']=((z.flank_offset_CDS==-1)&(z.HxlR_strand=='+')&(z.neighbor_strand=='-')) |((z.flank_offset_CDS==1)&(z.HxlR_strand=='-')&(z.neighbor_strand=='+'))
 z['close_divergent']=z.divergent & z.intergenic_distance_bp.between(0,500)
 z['module']=mod;detail.append(z)
 res.append({'module':mod,'n_HMT':z.HMT.nunique(),'n_genomes':z.GCA.nunique(),'n_HxlR_proteins':z.HxlR_protein.nunique(),'n_close_divergent':int(z.close_divergent.sum()),'median_intergenic_gap_bp':float(z.intergenic_distance_bp.median()) if len(z) else None,
 'plaque_HMT':z[z.Ecology=='plaque'].HMT.nunique(),'mucosal_HMT':z[z.Ecology=='mucosal'].HMT.nunique()})
result=pd.DataFrame(res)
result.to_csv(O/'HxlR_four_target_modules_cross_HMT_context.tsv',sep='\t',index=False)
pd.concat(detail).to_csv(O/'HxlR_four_target_modules_evidence_per_genome.tsv',sep='\t',index=False)
print(result.to_string(index=False))
