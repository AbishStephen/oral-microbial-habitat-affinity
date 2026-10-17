#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import csv
import numpy as np,pandas as pd
from scipy.stats import pearsonr
B=Path('__PROJECT_ROOT__/study13')
O=B/'08_validation/confound_audit';O.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(O/'hmt_covariates_and_predictions.tsv',sep='\t')
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rr=csv.reader(f,delimiter='\t');hd=next(rr);rows=list(rr)
ids=np.array([z[0] for z in rows]);g=np.array([z[1] for z in rows]);y=np.array([float(z[4]) for z in rows])
dd=d.set_index('HMT').loc[ids]
X=np.asarray([[float(t) for t in z[5:]] for z in rows],dtype=np.float64)
feature=np.array(hd[5:])
# Residualize within each genus, then jointly adjust genome-size, GC, CDS.
def residual(z,groups):
 a=z.copy()
 for v in np.unique(groups):
  sel=groups==v;a[sel]-=a[sel].mean(axis=0)
 return a
Y0=residual(y,g)
XX=residual(X,g)
C=residual(dd[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(dtype=float),g)
assert np.isfinite(C).all()
C=(C/(np.std(C,axis=0,ddof=1)+1e-9))
U,ss,vt=np.linalg.svd(C,full_matrices=False)
rank=int(np.sum(ss>ss.max()*1e-8));Q=U[:,:rank]
def project(z):return z-Q@(Q.T@z)
Y=project(Y0);F=project(XX)
norm=np.sqrt((F*F).sum(0));valid=norm>1e-12
raw=(Y@F)/np.sqrt((Y@Y)*norm**2+1e-30);raw[~valid]=0
groups=[np.where(g==v)[0] for v in np.unique(g)]
rng=np.random.default_rng(20261012)
nperm=5000;counter=np.zeros(len(feature),dtype=np.int32);maxvals=[]
# Freedman-Lane: permute residuals from Y ~ genus + size+GC+CDS, reproject nuisance
for first in range(0,nperm,50):
 k=min(50,nperm-first);Yp=np.tile(Y,(k,1))
 for j in range(k):
  for ids0 in groups:Yp[j,ids0]=rng.permutation(Yp[j,ids0])
 Yp=Yp-(Yp@Q)@Q.T
 co=Yp@F
 den=np.sqrt((Yp*Yp).sum(1))[:,None]*norm[None,:]
 r=np.divide(co,den,out=np.zeros_like(co),where=den>1e-15)
 abs_r=np.abs(r)
 counter+=(abs_r>=np.abs(raw)[None,:]-1e-14).sum(0)
 maxvals.extend(abs_r.max(axis=1).tolist())
mx=np.array(maxvals)
adj=np.array([(np.count_nonzero(mx>=v-1e-14)+1)/(nperm+1) for v in np.abs(raw)])
point=(counter+1)/(nperm+1)
cands=['PFAM:PF04072','PFAM:PF12680','PFAM:PF10604','PFAM:PF08020','PFAM:PF14534','PFAM:PF13199','PFAM:PF22815','KMOD:M00117','RHEA:RHEA:27770']
out=pd.DataFrame({'feature':feature,'genus_genomesize_GC_CDS_partial_r':raw,'pointwise_perm_p':point,'maxT_FWER_p':adj})
out.sort_values('maxT_FWER_p').to_csv(O/'genome_covariate_adjusted_maxT_all_features.tsv',sep='\t',index=False)
sel=out.set_index('feature').reindex(cands)
sel.to_csv(O/'genome_covariate_adjusted_candidates.tsv',sep='\t')
# genus jackknife for candidates, with covariates re-fit on reduced samples
lrows=[]
for fn in cands:
 if fn not in feature:continue
 j=int(np.where(feature==fn)[0][0])
 for omit in np.unique(g):
  t=g!=omit
  gg=g[t];yy=residual(y[t],gg);ff=residual(X[t,j],gg);cc=residual(dd[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(dtype=float)[t],gg)
  cc=cc/(np.std(cc,axis=0,ddof=1)+1e-9)
  uu,ss,v=np.linalg.svd(cc,full_matrices=False);q=uu[:,:int(np.sum(ss>ss.max()*1e-8))]
  yr=yy-q@(q.T@yy);fr=ff-q@(q.T@ff)
  v=float(pearsonr(yr,fr).statistic) if np.std(fr)>1e-12 else np.nan
  lrows.append([fn,omit,v])
pd.DataFrame(lrows,columns=['feature','omitted_genus','partial_r']).to_csv(O/'feature_leave_one_genus_out_adjusted.tsv',sep='\t',index=False)
print('NUISANCE_RANK',rank,'residual_df',len(y)-len(np.unique(g))-rank-1,'maxT_null_95%',np.quantile(mx,.95))
print('CANDIDATES:\n',sel.to_string())
print('LEAVE-ONE-GENUS:\n',pd.DataFrame(lrows,columns=['f','omit','r']).groupby('f').r.agg(['min','max']).to_string())
