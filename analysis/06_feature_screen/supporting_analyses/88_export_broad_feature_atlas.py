#!/usr/bin/env python3
"""Export precise HOMD47 paired feature matrices to compact, transferable JSON."""
from pathlib import Path
import csv,json,zlib,base64
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13');O=B/'10_additional_analyses_20261012'
groups={
'Candidate / poorly characterized':[
('PFAM:PF04072','LCM domain'),('PFAM:PF12680','SnoaL_2 domain'),('PFAM:PF10604','Polyketide_cyc2 domain'),('PFAM:PF08020','DUF1706'),('PFAM:PF13635','DUF4143'),('PFAM:PF30184','DUF8612')],
'Regulation & sensing':[
('PFAM:PF01638','HxlR-like regulator'),('PFAM:PF02311','AraC ligand-binding'),('PFAM:PF09339','IclR-like HTH'),('PFAM:PF13443','HTH_26'),('PFAM:PF13188','PAS_8 sensor'),('PFAM:PF00376','MerR family')],
'Envelope / interface':[
('PFAM:PF09374','PG-binding_3'),('PFAM:PF04226','Transgly_assoc'),('PFAM:PF18075','FtsX extracellular'),('PFAM:PF08139','LPAM_1'),('PFAM:PF02129','Peptidase_S15')],
'Metabolism & redox':[
('PFAM:PF00903','Glyoxalase domain'),('PFAM:PF12681','Glyoxalase_2'),('PFAM:PF12682','Flavodoxin_4'),('PFAM:PF03358','FMN reductase'),('PFAM:PF02737','3HCDH N-terminal'),('PFAM:PF04191','PEMT domain'),('PFAM:PF13386','DsbD_2'),('KMOD:M00045','Histidine degradation'),('KMOD:M00530','Nitrate → ammonia'),('KMOD:M00023','Tryptophan biosynthesis'),('KPATH:map00300','Lysine biosynthesis'),('KPATH:map00052','Galactose metabolism')],
'Mobile elements & defence':[
('PFAM:PF07275','ArdA antirestriction'),('PFAM:PF02534','T4SS-DNA transfer'),('PFAM:PF19044','TraG P-loop'),('PFAM:PF03432','Relaxase'),('PFAM:PF10145','Phage minor-tail'),('PFAM:PF09827','CRISPR Cas2'),('PFAM:PF13333','Integrase rve_2'),('PFAM:PF07751','Abi_2'),('PFAM:PF06414','Zeta toxin')],
'DNA repair & maintenance':[
('PFAM:PF08713','DNA alkylation repair'),('PFAM:PF03852','Vsr mismatch repair'),('PFAM:PF02565','RecO C-terminal')]}
features=[{'id':key,'short':label,'group':grp} for grp,entries in groups.items() for key,label in entries]
sister=pd.read_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t').fillna('')
extreme=pd.read_csv(O/'extreme13_genus_matched_HMT_pairs.tsv',sep='\t').fillna('')
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');h=next(rd);data=list(rd)
needed=[x['id'] for x in features]
missing=set(needed)-set(h)
assert not missing,missing
cols=[h.index(key) for key in needed]
values={row[0]:[float(row[i]) for i in cols] for row in data}
def build(P,typ):
 row=[]
 for _,p in P.iterrows():
  a,b=p.plaque_HMT,p.mucosal_HMT
  v=np.array(values[a])-np.array(values[b])
  obs=dict(type=typ,genus=p.Genus,p=a,m=b,delta_Y=float(p.delta_Y if 'delta_Y' in p else p.Y_separation),
          bootstrap=int(p.support) if 'support' in p else None,high=bool(p.both_strong_2p5) if 'both_strong_2p5' in p else None,delta=[round(float(z),5) for z in v])
  row.append(obs)
 return row
s=build(sister,'sister');e=build(extreme,'extreme')
s_stats=pd.read_csv(O/'sister8_all13917_exact_signflip_maxT.tsv',sep='\t').set_index('feature')
adj=pd.read_csv(B/'08_validation/confound_audit/genome_covariate_adjusted_maxT_all_features.tsv',sep='\t').set_index('feature')
for f in features:
 k=f['id'];z=s_stats.loc[k];v=adj.loc[k]
 f['sister_plaque']=int(z.plaque_higher);f['sister_mucosal']=int(z.mucosal_higher);f['sister_nonzero']=int(z.n_informative)
 f['sister_p']=round(float(z.exact_signflip_p_two_sided),6)
 f['sister_fwer']=round(float(z.exact_maxT_FWER_p),6)
 f['adjusted_r']=round(float(v.genus_genomesize_GC_CDS_partial_r),4)
 f['adjusted_fwer']=round(float(v.maxT_FWER_p),4)
obj=dict(features=features,sister=s,extreme=e,
         note='Each value is genomic annotation prevalence in plaque HMT minus mucosal HMT; KMOD/KPATH scores not directly identical to protein-family prevalence; pathway scores are annotation-derived, not measured metabolic flux.')
jsonout=O/'feature_atlas_41feature_pair_data.json'
jsonout.write_text(json.dumps(obj,separators=(',',':')))
print('NFEATURES',len(features),'N_PAIRS',len(s),len(e),'JSON_BYTES',jsonout.stat().st_size)
print('ENCODED_BEGIN')
print(base64.b64encode(zlib.compress(jsonout.read_bytes(),level=9)).decode('ascii'))
print('ENCODED_END')
