#!/usr/bin/env python3
"""Independent oral-HMT challenge for ten Figure2B Pfams using 631 representative
HOMD assemblies, HMMER --cut_ga, 4-source ecological labels. Genus-balanced
within-genus label permutations; exclude 47 primary HMTs for holdout.
All results exploratory: earlier 13,917-wide FWER tests were negative.
"""
from pathlib import Path
import numpy as np,pandas as pd,collections
B=Path('__PROJECT_ROOT__/study13')
O=B/'19_Figure2B_14factor_mechanisms_20261017'
inputfile=O/'HOMD631_Fig2B_12Pfams.domtblout'
ref=B/'14_HxlR_fullHOMD_controls_20261014/HOMD631_HxlR_SnoaL2_4source_ecology_annotations.tsv'
refs=pd.read_csv(ref,sep='\t')
assert refs.GCA.nunique()==631
original=set(pd.read_csv(B/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t').HMT)
name={'PF01638':'HxlR','PF04072':'LCM','PF13577':'SnoaL_4','PF13676':'TIR_2','PF22166':'WHD_ArnR1',
'PF01745':'IPT','PF01337':'Barstar','PF13413':'HTH_25','PF18075':'FtsX_ECD',
'PF14470':'bPH_3','PF12680':'SnoaL_2','PF07366':'SnoaL'}
rows=[]
for line in inputfile.open():
 if line.startswith('#'):continue
 z=line.split()
 if len(z)<22:continue
 acc=z[1].split('.')[0];protein=z[3]
 if acc not in name:continue
 # - GA threshold determined by hmmscan --cut_ga, domain evalue for ranking only
 rows.append(dict(GCA=protein.split('|')[0],protein_id=protein,pfam=acc,domain_i_evalue=float(z[12]),
   aln_start=int(z[17]),aln_end=int(z[18]),sequence_aa=int(z[5])))
hits=pd.DataFrame(rows)
hits.to_csv(O/'HOMD631_Fig2B_12Pfams_hits_all.tsv',sep='\t',index=False)
assert set(hits.GCA).issubset(set(refs.GCA))
counts=hits.groupby(['GCA','pfam']).protein_id.nunique().unstack(fill_value=0)
x=refs[['GCA','HMT_link','main_genus','main_species','oral','n_usable_sources','quality_score','ecology_4source_conservative']].merge(counts,left_on='GCA',right_index=True,how='left')
for pf in name:x[pf]=x.get(pf,pd.Series(0,index=x.index)).fillna(0).astype(int)
x['HMT_in_original_47']=x.HMT_link.isin(original)
x.to_csv(O/'HOMD631_Fig2B_12Pfams_presence_and_ecology.tsv',sep='\t',index=False)
small=x[(x.oral==True)&x.ecology_4source_conservative.isin(['plaque','mucosal'])].copy()
assert len(small)==155,len(small)
rng=np.random.default_rng(20261017)
out=[]
N=10000
for case,z in [('all_155_oral_HMTs',small),('new_HMTs_only',small[~small.HMT_in_original_47])]:
 gs=[d for gen,d in z.groupby('main_genus') if d.ecology_4source_conservative.nunique()==2]
 for pf,desc in name.items():
  gene_groups=[];rawdiff=[]
  for t in gs:
   vals=(t[pf]>0).to_numpy(float)
   n=int((t.ecology_4source_conservative=='plaque').sum())
   obs=float((t.loc[t.ecology_4source_conservative=='plaque',pf]>0).mean()-(t.loc[t.ecology_4source_conservative=='mucosal',pf]>0).mean())
   rawdiff.append(obs);gene_groups.append((vals,n))
  observed=float(np.mean(rawdiff)) if rawdiff else float('nan')
  null=np.zeros(N)
  for vals,n in gene_groups:
   for i in range(N):
    idx=rng.permutation(len(vals));null[i]+=(vals[idx[:n]].mean()-vals[idx[n:]].mean())/len(gene_groups)
  p=float((int((np.abs(null)>=abs(observed)-1e-12).sum())+1)/(N+1)) if gene_groups else float('nan')
  plaque=z[z.ecology_4source_conservative=='plaque'];mucosal=z[z.ecology_4source_conservative=='mucosal']
  out.append(dict(case=case,feature=pf,annotation=desc,n_HMT=len(z),n_mixed_habitat_genera=len(gs),
  positive_plaque=int((plaque[pf]>0).sum()),n_plaque=len(plaque),
  positive_mucosal=int((mucosal[pf]>0).sum()),n_mucosal=len(mucosal),
  raw_prevalence_difference=float((plaque[pf]>0).mean()-(mucosal[pf]>0).mean()),
  genus_balanced_difference=observed,genus_stratified_permutation_p=p))
res=pd.DataFrame(out)
for case,z in res.groupby('case'):
 p=z.genus_stratified_permutation_p.to_numpy();idx=np.argsort(p);q=np.empty(len(p));q[idx]=np.minimum.accumulate((p[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];res.loc[z.index,'BH_q_within_12_posthoc_Pfams']=np.minimum(q,1)
res.to_csv(O/'HOMD631_Fig2B_12Pfams_genus_balanced_ecology.tsv',sep='\t',index=False)
print('GENE_HITS',len(hits),'GENOMES',hits.GCA.nunique(),'N_ORAL_LABELED',len(small),'ORIGINAL_HMTs',small.HMT_in_original_47.sum(),flush=True)
print('RESULTS\n'+res.sort_values(['case','genus_stratified_permutation_p']).to_string(index=False),flush=True)
# CF: family co-hit same protein vs same genome and contrast with HxlR. PF13577 vs PF12680.
sets={pf:set(hits.loc[hits.pfam==pf,'protein_id']) for pf in name}
pf4=sets['PF13577'];pf2=sets['PF12680']
within=len(pf4&pf2)
by_genome=hits.groupby('pfam').GCA.apply(set).to_dict()
print('SNOAL4_SNOAL2',len(pf4),len(pf2),'same_protein',within,
 'same_genome',len(by_genome['PF13577']&by_genome['PF12680']),flush=True)
