#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
import numpy as np
from scipy.stats import rankdata,spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_absolute_error,r2_score,roc_auc_score,balanced_accuracy_score
B=Path("__PROJECT_ROOT__/study13")
H=B/"handoff"; O=B/"08_validation"/"permutation_null_expanded13"; O.mkdir(parents=True,exist_ok=True)
with open(H/"X_primary_transferable_DEDUP.tsv") as f:
    r=csv.reader(f,delimiter="\t");hd=next(r);rows=list(r)
hmts=np.array([x[0] for x in rows]); genera=np.array([x[1] for x in rows])
y0=np.array([float(x[4]) for x in rows]); X=np.array([[float(v) for v in x[5:]] for x in rows],float)
Ks=[10,25,50,100,250,500,1000];alphas=[.1,1,10,100,1000]
def crossmask(A,G):
    var=np.var(A,axis=0)>1e-12; keep=np.zeros(A.shape[1],bool)
    for j in np.where(var)[0]:
        if sum(np.any(np.abs(A[G==g,j])>1e-12) for g in np.unique(G))>=2:keep[j]=True
    return keep
def prep(A):
    R=rankdata(A,axis=0);R=R-R.mean(axis=0);return R,np.sqrt(np.sum(R*R,axis=0))
cache=[]
for og in np.unique(genera):
    ote=genera==og;otr=~ote;Xi=X[otr];Gi=genera[otr];inn=[]
    for ig in np.unique(Gi):
        ite=Gi==ig;itr=~ite;idx=np.where(crossmask(Xi[itr],Gi[itr]))[0];R,dx=prep(Xi[itr][:,idx]);inn.append((ite,itr,idx,R,dx))
    idx=np.where(crossmask(Xi,Gi))[0];R,dx=prep(Xi[:,idx]);cache.append((og,ote,otr,inn,idx,R,dx))
def sc(R,dx,y):
    yr=rankdata(y);yr=yr-yr.mean();dy=np.sqrt(np.sum(yr*yr));den=dx*dy;num=np.sum(R*yr[:,None],axis=0)
    return np.divide(np.abs(num),den,out=np.zeros_like(num,dtype=float),where=den>0)
def rho(a,b):
    z=spearmanr(a,b).statistic;return -9. if not np.isfinite(z) else float(z)
def one(y):
    pred=np.zeros(len(y))
    for og,ote,otr,inn,oidx,oR,odx in cache:
        Xi=X[otr];Y=y[otr];orders=[]
        for ite,itr,idx,R,dx in inn:
            order=idx[np.argsort(-sc(R,dx,Y[itr]))];orders.append((ite,itr,order))
        best=None
        for k in Ks:
            for a in alphas:
                pr=np.zeros(len(Y))
                for ite,itr,order in orders:
                    sel=order[:min(k,len(order))];m=make_pipeline(StandardScaler(),Ridge(alpha=a));m.fit(Xi[itr][:,sel],Y[itr]);pr[ite]=m.predict(Xi[ite][:,sel])
                key=(rho(Y,pr),-mean_absolute_error(Y,pr))
                if best is None or key>best[0]:best=(key,k,a)
        k,a=best[1],best[2];order=oidx[np.argsort(-sc(oR,odx,Y))];sel=order[:min(k,len(order))]
        m=make_pipeline(StandardScaler(),Ridge(alpha=a));m.fit(Xi[:,sel],Y);pred[ote]=m.predict(X[ote][:,sel])
    b=(y>0).astype(int)
    return dict(rho=rho(y,pred),MAE=float(mean_absolute_error(y,pred)),R2=float(r2_score(y,pred)),AUC=float(roc_auc_score(b,pred)),BAcc=float(balanced_accuracy_score(b,pred>0)))
def perm(seed):
    rng=np.random.default_rng(seed);y=y0.copy()
    for g in np.unique(genera):
        ix=np.where(genera==g)[0];y[ix]=rng.permutation(y[ix])
    return y
ap=argparse.ArgumentParser();ap.add_argument("--chunk",type=int,required=True);ap.add_argument("--per-chunk",type=int,default=5);ap.add_argument("--seed",type=int,default=20261011);a=ap.parse_args()
rec=[]
for q in range(a.chunk*a.per_chunk,(a.chunk+1)*a.per_chunk):
    m=one(perm(a.seed+1000+q));m.update(permutation=q,seed=a.seed+1000+q);rec.append(m);print(q,m,flush=True)
with open(O/f"chunk_{a.chunk:03d}.tsv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=["permutation","seed","rho","MAE","R2","AUC","BAcc"],delimiter="\t");w.writeheader();w.writerows(rec)
