import csv, math
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_absolute_error,r2_score,roc_auc_score,balanced_accuracy_score

B=Path("__PROJECT_ROOT__/study13")
O=B/"06_hmt_matrices"; V=B/"08_validation"; H=B/"handoff"
V.mkdir(exist_ok=True); H.mkdir(exist_ok=True)
layers=[
 ("COG",O/"COG_HMT_prevalence.tsv"),
 ("EC",O/"EC_HMT_prevalence.tsv"),
 ("PFAM",O/"X2_Pfam_HMT_prevalence_TRANSFERABLE.tsv"),
 ("X5",O/"X5_HMT_prevalence_TRANSFERABLE.tsv"),
 ("RHEA",O/"X3_Rhea_HMT_prevalence_TRANSFERABLE.tsv"),
 ("KMOD",O/"X4_KEGG_module_HMT_ECcoverage_TRANSFERABLE.tsv"),
 ("KPATH",O/"X4_KEGG_pathway_HMT_ECcoverage_TRANSFERABLE.tsv"),
]
def load(path):
    with open(path) as f:
        r=csv.reader(f,delimiter="\t"); hd=next(r); rr=list(r)
    meta=[x[:6] for x in rr]
    X=np.array([[float(v) for v in x[6:]] for x in rr],float)
    return meta,hd[6:],X
base_meta=None; mats={}; names={}
for tag,p in layers:
    m,n,x=load(p)
    if base_meta is None: base_meta=m
    else: assert [z[0] for z in m]==[z[0] for z in base_meta]
    mats[tag]=x; names[tag]=n
hmts=np.array([m[0] for m in base_meta]); genera=np.array([m[1] for m in base_meta])
species=np.array([m[2] for m in base_meta]); ecology=np.array([m[3] for m in base_meta])
y=np.array([float(m[4]) for m in base_meta]); ybin=(ecology=="plaque").astype(int)
# global cross-genus filter (outcome-independent); PFAM/X5/RHEA/K* already filtered but recheck all
filtered={}; fnames={}
for tag in mats:
    X=mats[tag]; keep=[]
    for j in range(X.shape[1]):
        ng=sum(np.any(np.abs(X[genera==g,j])>1e-12) for g in np.unique(genera))
        if ng>=2 and np.var(X[:,j])>1e-12: keep.append(j)
    filtered[tag]=X[:,keep]; fnames[tag]=[f"{tag}:{names[tag][j]}" for j in keep]
    print(tag,"raw",X.shape[1],"transferable+variable",len(keep))
Xall=np.hstack([filtered[t] for t,_ in layers])
Nall=sum((fnames[t] for t,_ in layers),[])
# write integrated primary handoff matrix
with open(H/"X_primary_transferable.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["HMT","Genus","Species","Ecology","Y"]+Nall)
    for i in range(len(hmts)):
        w.writerow([hmts[i],genera[i],species[i],ecology[i],y[i]]+Xall[i].tolist())
# feature dictionary
with open(H/"feature_dictionary_primary.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["feature","layer","source_feature","n_HMTs_nonzero","n_genera_nonzero","variance"])
    k=0
    for tag,_ in layers:
        for j,nm in enumerate(fnames[tag]):
            col=filtered[tag][:,j]
            nh=int(np.sum(np.abs(col)>1e-12))
            ng=sum(np.any(np.abs(col[genera==g])>1e-12) for g in np.unique(genera))
            w.writerow([nm,tag,nm.split(":",1)[1],nh,ng,float(np.var(col))]); k+=1
assert k==Xall.shape[1]
def train_filter(X,tr_groups):
    mask=np.var(X,axis=0)>1e-12
    for j in np.where(mask)[0]:
        if sum(np.any(np.abs(X[tr_groups==g,j])>1e-12) for g in np.unique(tr_groups))<2: mask[j]=False
    return mask
alphas=[0.01,0.1,1,10,100,1000]
Cs=[0.001,0.01,0.1,1,10,100]
def inner_reg(X,Y,G):
    best=(None,-1e99)
    for a in alphas:
        pr=np.zeros(len(Y))
        for g in np.unique(G):
            te=G==g; tr=~te; mask=train_filter(X[tr],G[tr])
            if not np.any(mask): pr[te]=np.mean(Y[tr]); continue
            model=make_pipeline(StandardScaler(),Ridge(alpha=a))
            model.fit(X[tr][:,mask],Y[tr]); pr[te]=model.predict(X[te][:,mask])
        rho=spearmanr(Y,pr).statistic
        if np.isnan(rho): rho=-1e99
        if rho>best[1]: best=(a,rho)
    return best[0]
def inner_clf(X,Y,G):
    best=(None,-1e99)
    for c in Cs:
        pr=np.zeros(len(Y))
        for g in np.unique(G):
            te=G==g; tr=~te; mask=train_filter(X[tr],G[tr])
            if not np.any(mask): pr[te]=np.mean(Y[tr]); continue
            model=make_pipeline(StandardScaler(),LogisticRegression(C=c,solver="liblinear",class_weight="balanced",max_iter=10000))
            model.fit(X[tr][:,mask],Y[tr]); pr[te]=model.predict_proba(X[te][:,mask])[:,1]
        auc=roc_auc_score(Y,pr)
        if auc>best[1]: best=(c,auc)
    return best[0]
def benchmark(tag,X):
    preg=np.zeros(len(y)); pcl=np.zeros(len(y)); fold=[]
    for g in np.unique(genera):
        te=genera==g; tr=~te
        mask=train_filter(X[tr],genera[tr])
        a=inner_reg(X[tr][:,mask],y[tr],genera[tr]) if np.any(mask) else None
        c=inner_clf(X[tr][:,mask],ybin[tr],genera[tr]) if np.any(mask) else None
        if np.any(mask):
            mr=make_pipeline(StandardScaler(),Ridge(alpha=a)); mr.fit(X[tr][:,mask],y[tr]); preg[te]=mr.predict(X[te][:,mask])
            mc=make_pipeline(StandardScaler(),LogisticRegression(C=c,solver="liblinear",class_weight="balanced",max_iter=10000)); mc.fit(X[tr][:,mask],ybin[tr]); pcl[te]=mc.predict_proba(X[te][:,mask])[:,1]
        else:
            preg[te]=np.mean(y[tr]); pcl[te]=np.mean(ybin[tr])
        fold.append([tag,g,int(np.sum(te)),int(np.sum(mask)),a,c])
    rho=spearmanr(y,preg).statistic; mae=mean_absolute_error(y,preg); r2=r2_score(y,preg)
    auc=roc_auc_score(ybin,pcl); bacc=balanced_accuracy_score(ybin,(pcl>=.5).astype(int))
    with open(V/f"{tag}_outer_predictions.tsv","w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(["HMT","Genus","Ecology","Y","Y_pred","P_plaque"])
        for i in range(len(y)): w.writerow([hmts[i],genera[i],ecology[i],y[i],preg[i],pcl[i]])
    return [tag,X.shape[1],rho,mae,r2,auc,bacc],fold
results=[]; folds=[]
for tag,_ in layers:
    z,fd=benchmark(tag,filtered[tag]); results.append(z); folds+=fd
z,fd=benchmark("INTEGRATED",Xall); results.append(z); folds+=fd
with open(V/"benchmark_summary.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["layer","n_features","LOGO_Spearman_Y","LOGO_MAE_Y","LOGO_R2_Y","LOGO_AUC_ecology","LOGO_balanced_accuracy"]); w.writerows(results)
with open(V/"benchmark_fold_parameters.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["layer","heldout_genus","n_test","n_features_train","ridge_alpha","logistic_C"]); w.writerows(folds)
# exploratory cross-genus feature ranking
rank=[]
for j,nm in enumerate(Nall):
    x=Xall[:,j]; rho=spearmanr(x,y).statistic
    if np.isnan(rho): rho=0.0
    overall=np.mean(x[ybin==1])-np.mean(x[ybin==0]); s=1 if overall>=0 else -1
    diffs=[]
    for g in np.unique(genera):
        q=genera==g
        diffs.append(np.mean(x[q & (ybin==1)])-np.mean(x[q & (ybin==0)]))
    conc=sum((d*s)>0 for d in diffs)/len(diffs)
    rank.append((nm,nm.split(":",1)[0],rho,overall,conc,sum(abs(d)>1e-12 for d in diffs)))
rank.sort(key=lambda z:(abs(z[2])*z[4],abs(z[2])),reverse=True)
with open(V/"exploratory_crossgenus_feature_ranking.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["feature","layer","Spearman_with_Y","plaque_minus_mucosal_mean","genus_direction_concordance","genera_nonzero_difference"]); w.writerows(rank)
# Y/Z handoff
with open(H/"Y_ecology.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["HMT","Y","Ecology"]); [w.writerow([hmts[i],y[i],ecology[i]]) for i in range(len(y))]
with open(H/"Z_lineage.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["HMT","Genus","Species"]); [w.writerow([hmts[i],genera[i],species[i]]) for i in range(len(y))]
print("Integrated features:",Xall.shape[1])
for r in results: print(*r)
print("Top 10 exploratory:")
for r in rank[:10]: print(*r)
