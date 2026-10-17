#!/usr/bin/env python3
"""Figure 2B fourteen-feature mechanistic evidence audit. No recasting post hoc
features as replicated ecological associations. Incorporates HMT-level breadth,
1121-genome protein homology, and immediate gene neighborhoods."""
from pathlib import Path
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'19_Figure2B_14factor_mechanisms_20261017';O.mkdir(exist_ok=True)
f=pd.read_csv(B/'11_manuscript_completion_20261013/Figure2/Figure2_B_14focal_statistics.tsv',sep='\t')
sub=pd.read_csv(B/'11_manuscript_completion_20261013/candidate_Pfam_protein_subfamily_summary.tsv',sep='\t')
mapping={
'PF01638':('chemical-responsive HTH transcription factor family','redox/electrophile/carbonylic sensing varies by regulator subfamily','medium','already characterized HxlR oral subfamilies; chemical identity unverified'),
'PF04072':('LCM-family SAM-dependent methyltransferase','methylation of protein or other substrate; assess conserved active-site context','medium','test sequence subfamilies and neighbouring genes for substrate specificity'),
'PF13577':('SnoaL_4, SnoaL/NTF2-like fold','small molecule-related fold with unknown biochemical activity in oral proteins','high','not interchangeable with PF12680; examine efflux-associated streptococcal loci and holo-fold'),
'PF13676':('TIR_2 domain proteins','putative bacterial NADase/defense signalling; many paralogues','high','test defense-system architectures, NADase active sites and genetic linkage'),
'PF22166':('ArnR1-like winged-helix','transcription-related DNA binding, unknown target specificity','high','inspect redox, methylation and metabolic gene contexts by independent subfamily'),
'PF01745':('IPT isopentenyltransferase','often tRNA isopentenylation; broader pathway interpretation requires orthology','medium','test miaA/miaB tRNA-modification cassette and split by homologues'),
'PF01337':('Barstar RNase inhibitor family','possible protection against cognate secreted/intracellular RNases','high','confirm RNase-Barstar conserved adjacency and shared HMT-independent origins'),
'PF13413':('HTH_25 helix-turn-helix regulator domain','broad DNA-binding family with diverse regulons','medium','cluster paralogues and split into specific operon-linked families'),
'PF18075':('FtsX extracellular domain','FtsEX cell-envelope/division machine; peptidoglycan hydrolase activation','medium','FtsE and EnvC adjacency validate annotation; likely housekeeping confound'),
'PF14470':('bPH_3 bacterial PH-like domain','protein interaction / membrane-associated fold conjectural','medium-high','test SHOCT pairing and mobile/phage island genomic context'),
'COG:COG1216':('GT2 glycosyltransferase family','glycan synthesis, substrate and polymer specificity unresolved','medium','split GT2 by catalytic/domain architecture and polysaccharide locus'),
'RHEA:12012':('thiamine kinase','thiamine salvage to thiamine monophosphate','medium','check thiamine kinase orthology, gene names and availability of alternate pathways'),
'RHEA:10532':('tryptophan synthase','TrpA/TrpB-associated tryptophan biosynthesis','medium','check trpA/trpB gene-level completeness and isolate auxotrophy; no enzymatic activity inferred'),
'M00545':('trans-cinnamate degradation module','multi-step aromatic degradation to acetyl-CoA','high audit priority','SOURCE SCORE IS PARTIAL EC COVERAGE, NO HMT COMPLETES EVEN THE OBSERVED 4 EC STEPS; replace pathway interpretation')}
row=[]
for r in f.itertuples():
 id=r.feature.replace('PFAM:','').replace('KMOD:','').replace('RHEA:RHEA:','RHEA:')
 ann= mapping.get(id,('unresolved','','',''))
 g=dict(feature=r.feature,feature_id=id,Figure2B_original_label=r.annotation,
        ecology_direction=r.direction,genus_support=r.genus_support,phyla_support=r.phyla_support,
        raw_genus_mean_effect=r.raw_effect_genus_mean,adjusted_genus_mean_effect=r.adj_effect_genus_mean,
        original_raw_p=r.raw_perm_p_t,original_adjusted_p=r.adj_perm_p_t,
        original_raw_maxT_FWER=r.raw_maxT_all_feature_FWER_p_t,
        original_adjusted_maxT_FWER=r.adj_maxT_all_feature_FWER_p_t,
        likely_function=ann[0],mechanistic_hypothesis=ann[1],work_priority=ann[2],required_validation=ann[3],
        method_caution='selected after screening 13,917 features, none survives FWER')
 if id.startswith('PF') and len(sub[sub.Pfam==id]):
  s=sub[sub.Pfam==id].iloc[0]
  g.update(original_n_proteins=s.n_proteins,original_n_genomes=s.n_genomes,
      original_n_HMT=s.n_HMT,original_n_genera=s.n_genera,
      original_MMseqs50_clusters=s.n_distinct_MMseqs_X5_clusters,
      most_broad_homology_genera=s.max_genera_in_a_single_cluster)
 row.append(g)
x=pd.DataFrame(row)
x.to_csv(O/'Figure2B_all14_mechanistic_evidence_audit.tsv',sep='\t',index=False)
assert len(x)==14,(len(x),x.feature)
assert all(x.original_raw_maxT_FWER>=0.05) and all(x.original_adjusted_maxT_FWER>=0.05)
print(x[['feature_id','ecology_direction','original_n_proteins','original_MMseqs50_clusters','original_adjusted_p','original_adjusted_maxT_FWER','likely_function']].to_string(index=False))
