#!/usr/bin/env python3
from pathlib import Path
import numpy as np,csv
from scipy.stats import t as tdist,binomtest
B=Path("__PROJECT_ROOT__/study13")
H=B/"handoff";O=B/"08_validation"/"within_genus_expanded13";O.mkdir(parents=True,exist_ok=True)
with open(H/"X_primary_transferable_DEDUP.tsv") as f:
 r=csv.reader(f,delimiter="\t");hd=next(r);rows=list(r)
features=np.array(hd[5:],dtype=object); X=np.array([[float(v) for v in z[5:]] for z in rows],float)
hmts=np.array([z[0] for z in rows]); G=np.array([z[1] for z in rows]); eco=np.array([z[3] for z in rows]); y=np.array([float(z[4]) for z in rows])
assert len(np.unique(G))==13
Xr=X.copy();yr=y.copy()
for g in np.unique(G):
 z=G==g;Xr[z]-=X[z].mean(0);yr[z]-=y[z].mean()
ssx=np.sum(Xr*Xr,axis=0);ssy=np.sum(yr*yr);cov=Xr.T@yr
r=np.divide(cov,np.sqrt(ssx*ssy),out=np.zeros(X.shape[1]),where=ssx>1e-15)
df=len(y)-len(np.unique(G))-1;t=r*np.sqrt(df/np.maximum(1-r*r,1e-15));p=2*tdist.sf(np.abs(t),df)
# direction recurrence
agree=[];inform=[];signp=[]
for j in range(X.shape[1]):
 overall=np.sign(r[j]);a=0;n=0
 for g in np.unique(G):
  z=G==g;d=X[z&(eco=="plaque"),j].mean()-X[z&(eco=="mucosal"),j].mean()
  if abs(d)>1e-12:n+=1;a+=int(d*overall>0)
 agree.append(a);inform.append(n);signp.append(binomtest(a,n,.5,alternative="greater").pvalue if n else np.nan)
# maxT genus-preserving perm
rng=np.random.default_rng(20261011);nperm=5000;obs=np.abs(r);count=np.zeros(len(features),int);maxabs=[]
denx=np.sqrt(ssx)
for st in range(0,nperm,100):
 m=min(100,nperm-st);Y=np.empty((m,len(y)))
 for b in range(m):
  yp=y.copy()
  for g in np.unique(G):
   ix=np.where(G==g)[0];yp[ix]=rng.permutation(yp[ix])
  for g in np.unique(G):
   z=G==g;yp[z]-=yp[z].mean()
  Y[b]=yp
 num=Y@Xr;den=np.sqrt(np.sum(Y*Y,axis=1))[:,None]*denx[None,:]
 RR=np.divide(num,den,out=np.zeros_like(num),where=den>1e-15);A=np.abs(RR);count+=(A>=obs[None,:]-1e-15).sum(0);maxabs.extend(A.max(1).tolist())
maxabs=np.array(maxabs)
point=(count+1)/(nperm+1);fwer=np.array([(np.sum(maxabs>=v-1e-15)+1)/(nperm+1) for v in obs])
cands={
"PFAM:PF12680":"SnoaL_2","PFAM:PF04072":"LCM","RHEA:RHEA:27770":"UbiG/CoQ3 EC 2.1.1.222",
"KMOD:M00117":"Ubiquinone biosynthesis","PFAM:PF02661":"Fic","PFAM:PF13199":"GH66",
"PFAM:PF22815":"CatAgl_D1","PFAM:PF14534":"DUF4440","PFAM:PF20730":"YetF_N",
"RHEA:RHEA:13085":"EC 2.3.1.89","PFAM:PF08020":"DUF1706","PFAM:PF10604":"Polyketide_cyc2"}
rec=[]
for j,f in enumerate(features):
 rec.append([f,cands.get(f,""),r[j],t[j],p[j],point[j],fwer[j],agree[j],inform[j],signp[j]])
rec.sort(key=lambda z:(z[6],z[5],-abs(z[2])))
with open(O/"all_features_expanded13_within_genus.tsv","w",newline="") as f:
 w=csv.writer(f,delimiter="\t");w.writerow(["feature","annotation","within_genus_r","t","p","perm_p_pointwise","perm_p_maxT_FWER","genera_agree","genera_informative","sign_p_one_sided"]);w.writerows(rec)
with open(O/"original_candidates_expanded13.tsv","w",newline="") as f:
 w=csv.writer(f,delimiter="\t");w.writerow(["feature","annotation","within_genus_r","p","perm_p_pointwise","perm_p_maxT_FWER","genera_agree","genera_informative","sign_p_one_sided"])
 for z in rec:
  if z[0] in cands:w.writerow([z[0],z[1],z[2],z[4],z[5],z[6],z[7],z[8],z[9]])
with open(O/"summary.tsv","w",newline="") as f:
 w=csv.writer(f,delimiter="\t");w.writerow(["n_HMTs","n_genera","n_features","n_permutations","maxT_q95","maxT_q99","FWER_lt_0.05"]);w.writerow([len(y),len(np.unique(G)),len(features),nperm,np.quantile(maxabs,.95),np.quantile(maxabs,.99),sum(fwer<.05)])
print("HMTs",len(y),"genera",len(np.unique(G)),"features",len(features),"FWER<.05",sum(fwer<.05))
print("Top 25")
for z in rec[:25]:print(z[:10])
print("Candidates")
for z in rec:
 if z[0] in cands:print(z[:10])
