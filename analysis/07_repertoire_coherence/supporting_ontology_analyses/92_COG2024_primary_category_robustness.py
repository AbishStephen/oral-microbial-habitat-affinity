#!/usr/bin/env python3
"""Predefined 2024 COG ontology groups: exact paired sign-flip class-level enrichment.

Population: 47 phenotype-selected HMTs. Distinct analyses:
1) Eight disjoint opposite-habitat tree cherries.
2) One extreme plaque and mucosal HMT per genus (13 genus-matched contrasts).
All inference conditional on the selected panel and COG2024 annotation labels.
"""
from pathlib import Path
from itertools import product
import csv, numpy as np, pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012';R=O/'ontology'
meta=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t')
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');head=next(rd);rows=list(rd)
H=[r[0] for r in rows];assert set(H)==set(meta.HMT)
X=np.array([[float(x) for x in z[5:]] for z in rows])
cols=np.array(head[5:]);cidx=np.flatnonzero([f.startswith('COG:') for f in cols])
F=X[:,cidx];fn=cols[cidx]
cats={};groups={};cogdesc={}
for line in (R/'cog-24.fun.tab').read_text().splitlines():
 a=line.split('\t')
 if len(a)>=4 and len(a[0])==1 and a[0].isalpha():cats[a[0]]=a[3];groups[a[0]]=a[1]
for line in (R/'cog-24.def.tab').read_text().splitlines():
 a=line.split('\t')
 if len(a)>=3:cogdesc[a[0]]=(a[1],a[2])
n=len(cats);category_names=list(cats)
M=np.zeros((len(fn),n))
for i,key in enumerate(fn):
 t=cogdesc.get(key.split(':')[1])
 if t is None:continue
 valid=[k for k in t[0][:1] if k in cats]
 for k in valid:M[i,category_names.index(k)]=1./len(valid)
nfeat=(M>0).sum(axis=0)
# Fixed ontology weights; each category is the mean annotation prevalence across its genes.
S=F@M / np.maximum(M.sum(axis=0),1)
d=meta.set_index('HMT').loc[H]
gg=pd.get_dummies(d.Genus,drop_first=True,dtype=float).to_numpy()
cc=d[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float)
cc=(cc-cc.mean(0))/(cc.std(0)+1e-12)
Z=np.column_stack([np.ones(len(d)),gg,cc]);q,_=np.linalg.qr(Z)
Sadj=S-q@(q.T@S)
keep=(nfeat>=10)&(S.std(axis=0)>1e-7)&(Sadj.std(axis=0)>1e-7)
S=S[:,keep];Sadj=Sadj[:,keep]
k=np.array(category_names)[keep];nfeat=nfeat[keep]
name={j:cats[j] for j in k}
C=pd.DataFrame(S,columns=k);C.insert(0,'HMT',H)
A=pd.DataFrame(Sadj,columns=k);A.insert(0,'HMT',H)
C.to_csv(O/'COG2024_primarycat_scores_HMT47.tsv',sep='\t',index=False)
A.to_csv(O/'COG2024_primarycat_adjusted_scores_HMT47.tsv',sep='\t',index=False)
pairs=pd.read_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t')
extreme=pd.read_csv(O/'extreme13_genus_matched_HMT_pairs.tsv',sep='\t')
idx={h:i for i,h in enumerate(H)}
def exact(D):
 norms=np.linalg.norm(D,axis=0)
 obs=np.divide(D.sum(0),norms,out=np.zeros(D.shape[1]),where=norms>1e-15)
 signs=np.array(list(product([-1.,1.],repeat=D.shape[0])))
 allstats=signs@D/np.maximum(norms,1e-15)
 pv=(np.abs(allstats)>=np.abs(obs)[None,:]-1e-11).mean(0)
 mx=np.abs(allstats).max(1)
 qv=(mx[:,None]>=np.abs(obs)[None,:]-1e-11).mean(0)
 return pv,qv,obs
records=[]
for label,pp in [('sister8',pairs),('extreme13',extreme),('extreme_strong_2p5',extreme[extreme.both_strong_2p5])]:
 ai=[idx[v] for v in pp.plaque_HMT];bi=[idx[v] for v in pp.mucosal_HMT]
 for adjustment,T in [('unadjusted',S),('genus_size_GC_CDS_residualized',Sadj)]:
  Delta=T[ai]-T[bi]
  pv,fw,stat=exact(Delta)
  for j,letter in enumerate(k):
   records.append(dict(panel=label,adjustment=adjustment,COG_category=letter,
     ontology_label=name[letter],n_annotated_COGs=int(nfeat[j]),n_pairs=len(pp),
     plaque_higher=int(np.sum(Delta[:,j]>1e-9)),
     mucosa_higher=int(np.sum(Delta[:,j]< -1e-9)),
     mean_category_delta=float(Delta[:,j].mean()),
     median_category_delta=float(np.median(Delta[:,j])),
     normalized_stat=float(stat[j]),exact_two_sided_p=float(pv[j]),
     exact_maxT_category_FWER_p=float(fw[j])))
r=pd.DataFrame(records).sort_values(['panel','adjustment','exact_two_sided_p'])
r.to_csv(O/'COG2024_ontology_set_enrichment_exact_signflip.tsv',sep='\t',index=False)
for panel in ['sister8','extreme13','extreme_strong_2p5']:
 a=r[(r.panel==panel)&(r.adjustment=='genus_size_GC_CDS_residualized')]
 print('\n',panel,'ADJUSTED_RESULTS\n',a[['COG_category','ontology_label','n_annotated_COGs','plaque_higher','mucosa_higher','mean_category_delta','exact_two_sided_p','exact_maxT_category_FWER_p']].to_string(index=False),flush=True)
print('CATEGORIES',len(k),'GROUP SIZES',dict(zip(k,nfeat)),'MATRIX',S.shape,flush=True)
