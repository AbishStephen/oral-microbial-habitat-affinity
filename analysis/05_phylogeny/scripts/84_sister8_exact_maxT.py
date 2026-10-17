#!/usr/bin/env python3
"""HOMD47 eight habitat-discordant terminal sister contrasts, exact signflip all-feature maxT."""
from pathlib import Path
from scipy.stats import binomtest
import numpy as np,pandas as pd,csv
from itertools import product
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012';O.mkdir(exist_ok=True,parents=True)
P=pd.read_csv(O/'HOMD47_terminal_sister_pairs_bootstrap200.tsv',sep='\t')
P=P[P.discordant==True].copy()
assert len(P)==8
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 reader=csv.reader(f,delimiter='\t');hdr=next(reader);dat=list(reader)
ids=np.array([r[0] for r in dat]);names=np.array(hdr[5:])
X=np.array([[float(z) for z in r[5:]] for r in dat],dtype=np.float64)
print('MATRIX',X.shape,flush=True)
meta=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT')
idx={h:i for i,h in enumerate(ids)}
P['plaque_HMT']=np.where(P.Ecology1=='plaque',P.HMT1,P.HMT2)
P['mucosal_HMT']=np.where(P.Ecology1=='mucosal',P.HMT1,P.HMT2)
P['Genus']=P.genus_1
P['delta_Y']=P.apply(lambda z:meta.loc[z.plaque_HMT,'Y']-meta.loc[z.mucosal_HMT,'Y'],axis=1)
for v in ['GenomeSize_Mb','GC_percent','CDS','n_genomes','Y_pred','PF04072','PF12680','PF10604','PF13199','PF22815']:
 P['delta_'+v]=P.apply(lambda z:float(meta.loc[z.plaque_HMT,v])-float(meta.loc[z.mucosal_HMT,v]),axis=1)
pidx=[idx[z] for z in P.plaque_HMT];midx=[idx[z] for z in P.mucosal_HMT]
Delta=X[pidx]-X[midx];support=P.support.to_numpy()
assert all(P.plaque_HMT.map(meta.Ecology)=='plaque')
P.to_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t',index=False)
# Candidate features, plus possible other labels.
targets=['PFAM:PF04072','PFAM:PF12680','PFAM:PF10604','PFAM:PF08020','PFAM:PF14534','PFAM:PF13199','PFAM:PF22815','KMOD:M00117','RHEA:RHEA:27770']
def summarize(F,kind,signs,dd,names,threshold_mask=None):
 # exact randomization: under pair-exchangeability, flip ecology labels within each phylogenetic cherry.
 # score per feature is signed sum divided by observed L2 norm, avoiding a single high-magnitude gene dominating as strongly as a mean.
 norms=np.sqrt(np.sum(dd**2,axis=0))
 obs=np.divide(dd.sum(axis=0),norms,out=np.zeros(len(names)),where=norms>1e-12)
 counter=np.zeros(len(names),dtype=np.int32);max_stats=[]
 for sign in signs:
  scores=np.divide(np.asarray(sign)@dd,norms,out=np.zeros(len(names)),where=norms>1e-12)
  ab=np.abs(scores)
  counter+=(ab>=np.abs(obs)-1e-12)
  max_stats.append(ab.max())
 max_stats=np.array(max_stats)
 raw=counter/len(signs)
 adj=np.array([np.mean(max_stats>=v-1e-12) for v in np.abs(obs)])
 output=[]
 for j,f in enumerate(names):
  delta=dd[:,j];v=delta[np.abs(delta)>1e-12]
  output.append(dict(feature=f,test=kind,n_pairs=len(dd),n_informative=len(v),plaque_higher=int((delta>1e-12).sum()),mucosal_higher=int((delta< -1e-12).sum()),zero_pairs=int((np.abs(delta)<=1e-12).sum()),mean_prevalence_difference=float(delta.mean()),median_prevalence_difference=float(np.median(delta)),median_nonzero_difference=float(np.median(v)) if len(v) else 0.,direction_consistency=max(sum(v>0),sum(v<0))/len(v) if len(v) else np.nan,exact_signflip_p_two_sided=float(raw[j]),exact_maxT_FWER_p=float(adj[j]),signed_L2_stat=float(obs[j])))
 return pd.DataFrame(output).sort_values(['exact_signflip_p_two_sided','exact_maxT_FWER_p','signed_L2_stat'],ascending=[True,True,False])
sgns=np.array(list(product((-1.,1.),repeat=len(Delta))))
out=summarize(Delta,'sister8',sgns,Delta,names)
out.to_csv(O/'sister8_all13917_exact_signflip_maxT.tsv',sep='\t',index=False)
sel=out[out.feature.isin(targets)].copy()
sel.to_csv(O/'sister8_candidate_features.tsv',sep='\t',index=False)
print('SISTER8 CANDIDATES\n',sel.to_string(index=False),flush=True)
print('SISTER8 TOP20\n',out.head(20).to_string(index=False),flush=True)
print('SISTER8 GENUS COUNTS\n',P[['Genus','plaque_HMT','mucosal_HMT','support','delta_Y','delta_GenomeSize_Mb','delta_CDS','delta_PF04072','delta_PF12680','delta_PF10604']].to_string(index=False),flush=True)
print('ALL FEATURE FWER<0.05',sum(out.exact_maxT_FWER_p<.05),flush=True)
# Exclude unsupported 59% Gemella: exact 2^7; with p floor 2/128
keep=P.support>=70
DK=Delta[keep]
sk=np.array(list(product((-1.,1.),repeat=len(DK))))
sm=summarize(DK,'sister7_bootstrap70',sk,DK,names)
sm[sm.feature.isin(targets)].to_csv(O/'sister7_bootstrap70_candidate_sensitivity.tsv',sep='\t',index=False)
print('SISTER7 CANDIDATES\n',sm[sm.feature.isin(targets)].to_string(index=False),flush=True)
