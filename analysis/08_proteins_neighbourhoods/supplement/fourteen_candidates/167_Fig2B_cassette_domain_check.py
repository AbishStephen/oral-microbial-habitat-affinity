#!/usr/bin/env python3
import pandas as pd,gzip,collections
from pathlib import Path
B=Path('__PROJECT_ROOT__/study13')
O=B/'19_Figure2B_14factor_mechanisms_20261017'
d=pd.read_csv(O/'Fig2B_10Pfams_immediately_adjacent_CDS.tsv',sep='\t')
candidate={
'Barstar_pair_RNase': d[(d.Pfam=='PF01337')&d.neighbor_product.str.contains('ribonuclease',case=False,na=False)].copy(),
'SnoaL4_pair_MDR':d[(d.Pfam=='PF13577')&d.neighbor_product.str.contains('multidrug resistance protein 3',case=False,na=False)].copy(),
'WHD_ArnR1_pair_redox':d[(d.Pfam=='PF22166')&d.neighbor_product.str.contains('quinone oxidoreductase|ubiquinone biosynthesis',case=False,na=False)].copy(),
}
ids=set()
for t in candidate.values():ids.update(t.neighbor_id);ids.update(t.protein_id)
pfams=collections.defaultdict(list)
src=B/'04_domains/Pfam_domain_hits.tsv.gz'
with gzip.open(src,'rt') as f:
 cols=next(f).rstrip().split('\t')
 for line in f:
  z=line.rstrip().split('\t')
  if z[0] not in ids:continue
  pfams[z[0]].append(z[3]+':'+z[4])
rows=[]
for kind,t in candidate.items():
 for r in t.itertuples():
  rows.append(dict(system=kind,GCA=r.GCA,HMT=r.HMT,Genus=r.Genus,
   focal_protein=r.protein_id,neighbor_id=r.neighbor_id,
   source_domain=r.Pfam,neighbor_product=r.neighbor_product,
   gap_bp=r.intergenic_bp,same_strand=r.same_strand,divergent=r.divergent,
   focal_protein_Pfam=';'.join(sorted(set(pfams[r.protein_id]))),
   neighbor_Pfam=';'.join(sorted(set(pfams[r.neighbor_id])))))
df=pd.DataFrame(rows)
df.to_csv(O/'Figure2B_three_cassettes_Pfam_crossvalidation.tsv',sep='\t',index=False)
for name,t in df.groupby('system'):
 print('CASSETTE',name,'genomes',t.GCA.nunique(),'HMT',t.HMT.nunique(),'genera',t.Genus.nunique(),flush=True)
 print('NEIGHBOR_PFAM',t.neighbor_Pfam.value_counts().head(8).to_dict(),flush=True)
 print('CONSERVED_STRAND',t.same_strand.value_counts().to_dict(),'INTERGENIC_MEDIAN',t.gap_bp.median(),flush=True)
