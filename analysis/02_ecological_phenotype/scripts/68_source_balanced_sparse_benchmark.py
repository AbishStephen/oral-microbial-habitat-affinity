import csv
from pathlib import Path
import numpy as np
from scipy.stats import rankdata,spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_absolute_error,r2_score

B=Path('__PROJECT_ROOT__/study13')
H=B/'handoff'; V=B/'08_validation'
with open(H/'X_primary_transferable_DEDUP.tsv') as f:
    r=csv.reader(f,delimiter='\t'); hd=next(r); rows=list(r)
names=np.array(hd[5:],dtype=object)
hmts=np.array([x[0] for x in rows]); genera=np.array([x[1] for x in rows]); ecology=np.array([x[3] for x in rows])
y=np.array([float(x[4]) for x in rows])
sourceY={r["HMT"]:float(r["groupbalanced_Y_1em06"]) for r in csv.DictReader(open(V/"confound_audit/ecology_source_family_pseudocount_sensitivity.tsv"),delimiter="\t")}
y=np.array([sourceY[h] for h in hmts])
print("Using equally weighted three provenance groups, groups B=HMP/Eren/Dewhirst")
X=np.array([[float(v) for v in x[5:]] for x in rows],float)
Ks=[10,25,50,100,250,500,1000]
alphas=[0.1,1,10,100,1000]

def crossmask(X,G):
    var=np.var(X,axis=0)>1e-12
    idx=np.where(var)[0]; keep=np.zeros(X.shape[1],bool)
    for j in idx:
        ng=sum(np.any(np.abs(X[G==g,j])>1e-12) for g in np.unique(G))
        if ng>=2: keep[j]=True
    return keep

def abs_spearman(X,y):
    yr=rankdata(y); yr=yr-yr.mean(); yd=np.sqrt(np.sum(yr*yr))
    out=np.zeros(X.shape[1])
    for j in range(X.shape[1]):
        xr=rankdata(X[:,j]); xr=xr-xr.mean(); den=np.sqrt(np.sum(xr*xr))*yd
        out[j]=0 if den==0 else abs(np.sum(xr*yr)/den)
    return out

pred=np.zeros(len(y)); folds=[]; stability={}
for og in np.unique(genera):
    ote=genera==og; otr=~ote
    Xin=X[otr]; Yin=y[otr]; Gin=genera[otr]
    # Precompute inner rankings once per held-out inner genus
    inn=[]
    for ig in np.unique(Gin):
        ite=Gin==ig; itr=~ite
        m=crossmask(Xin[itr],Gin[itr]); idx=np.where(m)[0]
        sc=abs_spearman(Xin[itr][:,idx],Yin[itr]); order=idx[np.argsort(-sc)]
        inn.append((ite,itr,order))
    best=None
    for k in Ks:
        for a in alphas:
            pr=np.zeros(len(Yin))
            for ite,itr,order in inn:
                sel=order[:min(k,len(order))]
                if len(sel)==0: pr[ite]=np.mean(Yin[itr]); continue
                mod=make_pipeline(StandardScaler(),Ridge(alpha=a))
                mod.fit(Xin[itr][:,sel],Yin[itr]); pr[ite]=mod.predict(Xin[ite][:,sel])
            rho=spearmanr(Yin,pr).statistic
            if np.isnan(rho): rho=-9
            mae=mean_absolute_error(Yin,pr)
            key=(rho,-mae)
            if best is None or key>best[0]: best=(key,k,a)
    k,a=best[1],best[2]
    m=crossmask(Xin,Gin); idx=np.where(m)[0]; sc=abs_spearman(Xin[:,idx],Yin); order=idx[np.argsort(-sc)]
    sel=order[:min(k,len(order))]
    mod=make_pipeline(StandardScaler(),Ridge(alpha=a)); mod.fit(Xin[:,sel],Yin); pred[ote]=mod.predict(X[ote][:,sel])
    folds.append([og,int(np.sum(ote)),k,a,len(sel)])
    for j in sel: stability[names[j]]=stability.get(names[j],0)+1

rho=spearmanr(y,pred).statistic; mae=mean_absolute_error(y,pred); r2=r2_score(y,pred)
with open(V/'GROUPBAL_SPARSE_DEDUP_outer_predictions.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['HMT','Genus','Ecology','Y','Y_pred'])
    for i in range(len(y)): w.writerow([hmts[i],genera[i],ecology[i],y[i],pred[i]])
with open(V/'GROUPBAL_SPARSE_DEDUP_fold_parameters.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['heldout_genus','n_test','selected_k','ridge_alpha','n_features']); w.writerows(folds)
with open(V/'GROUPBAL_SPARSE_DEDUP_selection_stability.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['feature','outer_folds_selected'])
    for n,c in sorted(stability.items(),key=lambda z:(-z[1],z[0])): w.writerow([n,c])
print('SPARSE_DEDUP rho',rho,'MAE',mae,'R2',r2)
for x in folds: print(*x)
print('selected >=4 folds',sum(c>=4 for c in stability.values()))
print('selected all outer folds',sum(c==len(np.unique(genera)) for c in stability.values()))
