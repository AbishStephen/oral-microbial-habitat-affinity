#!/usr/bin/env python3
"""Cross-lineage 13-genus feature scan: genus-balanced effects, site-within-genus maxT permutations.
Panels limited to 47 HMTs. Genomes per HMT are NOT independent ecology replicates.
"""
from pathlib import Path
import csv,sys
import numpy as np,pandas as pd
from scipy.stats import t as student
B=Path('__PROJECT_ROOT__/study13');O=B/'10_additional_analyses_20261012'
N_PERM=int(sys.argv[1]) if len(sys.argv)>1 else 4096
SEED=20261012;batch=32
PHYLUM={'Actinomyces':'Actinomycetota','Rothia':'Actinomycetota','Campylobacter':'Campylobacterota',
'Fusobacterium':'Fusobacteriota','Leptotrichia':'Fusobacteriota',
'Gemella':'Bacillota','Selenomonas':'Bacillota','Streptococcus':'Bacillota','Veillonella':'Bacillota',
'Neisseria':'Pseudomonadota','Porphyromonas':'Bacteroidota','Prevotella':'Bacteroidota','Segatella':'Bacteroidota'}
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');hdr=next(rd);rows=list(rd)
names=np.array(hdr[5:],dtype=object);ids=np.array([r[0] for r in rows])
X=np.array([[float(v) for v in r[5:]] for r in rows],dtype=np.float32)
D=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT').loc[ids]
assert len(D)==47 and len(PHYLUM)==13 and all(D.Genus.isin(PHYLUM))
group_names=sorted(set(D.Genus));groups=[np.flatnonzero(D.Genus.to_numpy()==g) for g in group_names]
n=len(groups);num_feat=X.shape[1]
pidx=[np.flatnonzero((D.Genus.to_numpy()==g)&(D.Ecology.to_numpy()=='plaque')) for g in group_names]
midx=[np.flatnonzero((D.Genus.to_numpy()==g)&(D.Ecology.to_numpy()=='mucosal')) for g in group_names]
assert all(len(p)>0 and len(m)>0 for p,m in zip(pidx,midx))
W=np.zeros((n,len(D)),np.float32)
for i,(p,m) in enumerate(zip(pidx,midx)):W[i,p]=1/len(p);W[i,m]=-1/len(m)
phylum=np.array([PHYLUM[g] for g in group_names])
G=pd.get_dummies(D.Genus,drop_first=True,dtype=float).to_numpy()
C=D[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float)
C=(C-C.mean(0))/(C.std(0)+1e-12)
Z=np.column_stack([np.ones(len(D)),G,C]); Q=np.linalg.qr(Z)[0]
Xadj=(X-Q@(Q.T@X)).astype(np.float32)
# Within-genus prevalence contrasts; adjust entire HMT feature matrix for genomic covariates.
variants={'raw':X, 'genome_adjusted':Xadj}
f_sd=np.std(X.astype(float),axis=0)
effect_floor=np.maximum(0.05,0.25*f_sd).astype(np.float32)
def stat(Diff,delta):
 # Diff shape [batch, 13, 13917]
 pos=Diff>delta;neg=Diff< -delta
 score=(pos.sum(axis=1)-neg.sum(axis=1)).astype(np.float32)
 dmu=Diff.mean(axis=1);sd=Diff.std(axis=1,ddof=1)
 t=dmu/(sd/np.sqrt(n)+1e-7)
 # Cap t for negligible feature variation across genera
 return np.abs(score),np.abs(t)
observed={}
for label,XX in variants.items():
 E=W@XX
 score,t=stat(E[None,:,:],effect_floor[None,None,:])
 observed[label]={'E':E,'score':score[0],'t':t[0]}
 npl=(E>effect_floor[None,:]).sum(0);nmu=(E< -effect_floor[None,:]).sum(0)
 nzi=(np.abs(E)>effect_floor[None,:]).sum(0)
 npl0=(E>1e-7).sum(0);nmu0=(E< -1e-7).sum(0)
 phyla_pos=[];phyla_neg=[]
 for j in range(num_feat):
  phyla_pos.append(';'.join(sorted(set(phylum[E[:,j]>effect_floor[j]]))))
  phyla_neg.append(';'.join(sorted(set(phylum[E[:,j]<-effect_floor[j]]))))
 mean=E.mean(0);se=E.std(0,ddof=1)/np.sqrt(n)
 tmp=pd.DataFrame({'feature':names,'effect_genus_mean':mean,'effect_between_genus_sd':E.std(0,ddof=1),
 'effect_95CI_lo_descriptive':mean-student.ppf(.975,n-1)*se,
 'effect_95CI_hi_descriptive':mean+student.ppf(.975,n-1)*se,
 'n_plaque_genus_material':npl,'n_mucosal_genus_material':nmu,'n_informative_genus':nzi,
 'n_plaque_genus_any_delta':npl0,'n_mucosal_genus_any_delta':nmu0,
 'phylum_plaque':phyla_pos,'phylum_mucosa':phyla_neg,
 'n_phyla_plaque':[len(v.split(';')) if v else 0 for v in phyla_pos],
 'n_phyla_mucosa':[len(v.split(';')) if v else 0 for v in phyla_neg],
 'abs_material_directional_consistency':np.abs(npl-nmu),
 'equal_genus_studentized_abs_t':t[0],
 'meaningful_delta_floor':effect_floor,
 'genome_adjustment':label})
 for gi,g in enumerate(group_names):tmp['diff_'+g]=E[gi,:]
 tmp.to_csv(O/('crossgenus_allfeatures_'+label+'_effects.tsv'),sep='\t',index=False)
 observed[label]['result']=tmp
 print('OBS',label,'FEATURES',len(tmp),'>=6 CONSISTENT',int((tmp.abs_material_directional_consistency>=6).sum()),flush=True)
# Permute ecological labels within each genus, preserving nPlaque/nMucosa. Each genus is weighted once.
rng=np.random.default_rng(SEED)
exceed={name:{k:np.zeros(num_feat,dtype=np.int32) for k in ('consistency','t')} for name in variants}
maxstats={name:{k:[] for k in ('consistency','t')} for name in variants}
for start in range(0,N_PERM,batch):
 b=min(batch,N_PERM-start);WW=np.zeros((b,n,len(ids)),dtype=np.float32)
 for i,indices in enumerate(groups):
  p=len(pidx[i]);m=len(midx[i])
  for j in range(b):
   order=rng.permutation(indices)
   WW[j,i,order[:p]]=1./p
   WW[j,i,order[p:]]=-1./m
 mat=WW.reshape(b*n,len(ids))
 for lab,XX in variants.items():
  diff=(mat@XX).reshape(b,n,num_feat)
  cs,tt=stat(diff,effect_floor[None,None,:])
  exceed[lab]['consistency']+=np.count_nonzero(cs>=observed[lab]['score'][None,:]-1e-8,axis=0)
  exceed[lab]['t']+=np.count_nonzero(tt>=observed[lab]['t'][None,:]-1e-8,axis=0)
  maxstats[lab]['consistency'].extend(cs.max(axis=1))
  maxstats[lab]['t'].extend(tt.max(axis=1))
 if (start+batch)%512==0 or start+b==N_PERM:print('PERMS',start+b,'/',N_PERM,flush=True)
for lab in variants:
 r=observed[lab]['result']
 for k,ob in [('consistency',observed[lab]['score']),('t',observed[lab]['t'])]:
  mp=np.array(maxstats[lab][k]);point=(exceed[lab][k]+1)/(N_PERM+1)
  # +1 permits the observed statistic in the permutation universe
  fwer=(np.count_nonzero(mp[:,None]>=ob[None,:]-1e-8,axis=0)+1)/(N_PERM+1)
  r['perm_p_'+k]=point;r['maxT_all_feature_FWER_p_'+k]=fwer
 r['passes_crossphylum_3_and_genus_5']=((r.n_phyla_plaque>=3)&(r.n_plaque_genus_material>=5)) | ((r.n_phyla_mucosa>=3)&(r.n_mucosal_genus_material>=5))
 r.sort_values(['perm_p_consistency','perm_p_t'],inplace=True)
 r.to_csv(O/('crossgenus_allfeatures_'+lab+'_permutation_results.tsv'),sep='\t',index=False)
 print('RESULT',lab,'FWER_SIGN',int((r.maxT_all_feature_FWER_p_consistency<.05).sum()),'FWER_T',int((r.maxT_all_feature_FWER_p_t<.05).sum()))
 print('TOP\n',r[['feature','effect_genus_mean','n_plaque_genus_material','n_mucosal_genus_material','n_phyla_plaque','n_phyla_mucosa','perm_p_consistency','maxT_all_feature_FWER_p_consistency','perm_p_t','maxT_all_feature_FWER_p_t']].head(15).to_string(index=False),flush=True)
