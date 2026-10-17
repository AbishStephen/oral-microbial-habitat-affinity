#!/usr/bin/env python3
"""13 independent genus-matched extreme ecological contrasts. No pilot comparison."""
from pathlib import Path
import pandas as pd,numpy as np,csv
from itertools import product
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012';O.mkdir(exist_ok=True)
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t')
rows=[]
for genus,g in d.groupby('Genus'):
 a=g[g.Y>0].sort_values('Y',ascending=False).iloc[0]
 b=g[g.Y<0].sort_values('Y').iloc[0]
 z={'Genus':genus,'plaque_HMT':a.HMT,'mucosal_HMT':b.HMT,'plaque_Y':a.Y,'mucosal_Y':b.Y,'Y_separation':a.Y-b.Y,
    'both_strong_1p5':(a.Y>=1.5 and b.Y<=-1.5),
    'both_strong_2p5':(a.Y>=2.5 and b.Y<=-2.5),
    'both_strong_4p0':(a.Y>=4.0 and b.Y<=-4.0)}
 for feat in ['GenomeSize_Mb','GC_percent','CDS','PF04072','PF12680','PF10604','PF13199','PF22815']:
  z['delta_'+feat]=float(a[feat]-b[feat])
 rows.append(z)
P=pd.DataFrame(rows)
P.to_csv(O/'extreme13_genus_matched_HMT_pairs.tsv',sep='\t',index=False)
print('STRONG_PAIRS',len(P),int(P.both_strong_1p5.sum()),int(P.both_strong_2p5.sum()),int(P.both_strong_4p0.sum()),flush=True)
print(P[['Genus','plaque_HMT','mucosal_HMT','plaque_Y','mucosal_Y','both_strong_2p5','delta_PF04072','delta_PF12680','delta_PF10604']].to_string(index=False),flush=True)
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 reader=csv.reader(f,delimiter='\t');hdr=next(reader);data=list(reader)
names=np.array(hdr[5:]);ids=np.array([r[0] for r in data]);X=np.array([[float(v) for v in r[5:]] for r in data])
idx={h:i for i,h in enumerate(ids)}
D=X[[idx[h] for h in P.plaque_HMT]]-X[[idx[h] for h in P.mucosal_HMT]]
targets=['PFAM:PF04072','PFAM:PF12680','PFAM:PF10604','PFAM:PF08020','PFAM:PF14534','PFAM:PF13199','PFAM:PF22815','KMOD:M00117','RHEA:RHEA:27770']
target_idx=np.flatnonzero(np.isin(names,targets))
def test_delta(DD,subset,norm_corr=True):
 n=DD.shape[0]; norms=np.sqrt((DD*DD).sum(axis=0))
 obs=np.divide(DD.sum(axis=0),norms,out=np.zeros(DD.shape[1]),where=norms>1e-12)
 count=np.zeros(len(names),dtype=np.int32)
 max_obs=np.zeros(2**n,dtype=float)
 # 2^13 sign-flips exactly, stable process memory <10 MiB
 for i,sgn in enumerate(product((-1.,1.),repeat=n)):
  val=np.divide(np.asarray(sgn)@DD,norms,out=np.zeros(len(names)),where=norms>1e-12)
  a=np.abs(val);count+=(a>=np.abs(obs)-1e-12);max_obs[i]=a.max()
 out=[]
 for j in range(len(names)):
  raw=DD[:,j];nonzero=raw[np.abs(raw)>1e-12]
  if j in target_idx: out.append(dict(panel=subset,feature=names[j],n_pairs=n,n_informative=len(nonzero),plaque_higher=int((raw>1e-12).sum()),mucosal_higher=int((raw< -1e-12).sum()),zero_pairs=int((np.abs(raw)<=1e-12).sum()),mean_delta=float(raw.mean()),exact_unadjusted_two_sided=count[j]/(2**n),exact_allfeature_maxT_p=np.mean(max_obs>=abs(obs[j])-1e-12)))
 outdf=pd.DataFrame(out).sort_values('exact_unadjusted_two_sided')
 return outdf,int(np.count_nonzero(np.array([np.mean(max_obs>=val-1e-12) for val in abs(obs)])<0.05))
allreport=[]
for attr in [None,'both_strong_1p5','both_strong_2p5','both_strong_4p0']:
 mask=np.ones(len(P),dtype=bool) if attr is None else P[attr].to_numpy(dtype=bool)
 n=mask.sum()
 if n<3:continue
 if n<=13:
  df,nums=test_delta(D[mask],attr or 'one_extreme_pair_per_13genera')
  allreport.append(df)
  print('TEST',attr,'N_PAIRS',n,'GENOME_WIDE_FWER_HITS',nums,'CANDIDATES\n',df.to_string(index=False),flush=True)
pd.concat(allreport).to_csv(O/'extreme13_strongpref_threshold_candidate_exact.tsv',sep='\t',index=False)
