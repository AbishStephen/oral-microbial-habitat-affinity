#!/usr/bin/env python3
"""Expanded HMT47 phylogenetic GLS sensitivity; exploratory, not genomic GWAS."""
from pathlib import Path
import csv,re
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t as student_t, chi2
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'; O.mkdir(exist_ok=True)
class Node:
 def __init__(self,name='',length=0.0,children=None):
  self.name=name;self.length=length;self.children=children or [];self.parent=None
  for c in self.children:c.parent=self
def parse(s):
 s=s.strip().rstrip(';');i=0
 def rec():
  nonlocal i
  children=[]
  if s[i]=='(':
   i+=1
   while True:
    children.append(rec())
    if s[i]==',':i+=1;continue
    if s[i]==')':i+=1;break
  label=''
  while i<len(s) and s[i] not in ',():':label+=s[i];i+=1
  length=0.
  if i<len(s) and s[i]==':':
   i+=1;num=''
   while i<len(s) and s[i] not in ',()':num+=s[i];i+=1
   length=float(num) if num else 0.
  return Node(label,length,children)
 return rec()
tree=parse((B/'handoff/Z_HOMD47_ribosomal_tree.newick').read_text())
tips={}
def walk(node,depth):
 node.depth=depth
 if not node.children:tips[node.name.split('|')[0]]=node
 for c in node.children:walk(c,depth+max(c.length,0))
walk(tree,0.)
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 r=csv.reader(f,delimiter='\t');hd=next(r);rows=list(r)
h=np.array([x[0] for x in rows]);g=np.array([x[1] for x in rows]);y=np.array([float(x[4]) for x in rows])
assert set(h)==set(tips), (set(h)-set(tips),set(tips)-set(h))
fkeys=hd[5:]
feature={f:np.array([float(x[5+j]) for x in rows]) for j,f in enumerate(fkeys) if f in ['PFAM:PF04072','PFAM:PF12680','PFAM:PF10604','PFAM:PF08020','PFAM:PF14534','PFAM:PF13199','PFAM:PF22815','KMOD:M00117','RHEA:RHEA:27770']}
cov=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT').loc[h]
pnames=['GenomeSize_Mb','GC_percent','CDS']
meta=cov[pnames].to_numpy(float)
meta=(meta-meta.mean(axis=0))/(meta.std(axis=0)+1e-9)
base=pd.get_dummies(g,drop_first=True,dtype=float).to_numpy()
X0=np.column_stack([np.ones(len(y)),base,meta])
def anc(node):
 out=[]
 while node is not None:out.append(node);node=node.parent
 return out
ancs={t:anc(tips[t]) for t in h}
C=np.empty((len(h),len(h)))
for i,a in enumerate(h):
 aa=set(ancs[a])
 for j,b in enumerate(h):
  m=next(n for n in ancs[b] if n in aa)
  C[i,j]=m.depth
# Positive definite Brownian covariance from arbitrary rooted ML topology.
C+=np.eye(len(y))*1e-7
def fit(X,lam):
 v=lam*C+(1-lam)*np.diag(C.diagonal())
 v+=np.eye(len(y))*1e-7
 l=np.linalg.cholesky(v)
 xx=np.linalg.solve(l,X);yy=np.linalg.solve(l,y)
 bet=np.linalg.lstsq(xx,yy,rcond=None)[0];r=yy-xx@bet
 rss=r@r; n=len(y);p=np.linalg.matrix_rank(xx)
 nll=.5*(2*np.log(np.diag(l)).sum()+n*np.log(max(rss/n,1e-15)))
 se=np.sqrt(np.diag(np.linalg.pinv(xx.T@xx)) * rss / max(n-p,1))
 return nll,bet,se,rss,p
opt=minimize_scalar(lambda lam:fit(X0,lam)[0],bounds=(0,1),method='bounded',options={'xatol':1e-6})
lam=float(opt.x);nll0,bet0,se0,rss0,p0=fit(X0,lam)
out=[]
for k,x in feature.items():
 X=np.column_stack([X0,x])
 nll,b,se,rss,p=fit(X,lam)
 t=b[-1]/se[-1] if se[-1]>0 else np.nan
 out.append(dict(feature=k,beta=b[-1],se=se[-1],t=t,p_nominal=2*student_t.sf(abs(t),len(y)-p),partial_R2=max(0,1-rss/rss0),pgls_lambda_fixed=lam,df=len(y)-p))
pd.DataFrame(out).sort_values('p_nominal').to_csv(O/'PGLS47_genus_genome_covariates_candidates.tsv',sep='\t',index=False)
pd.DataFrame([dict(n_HMT=len(h),n_genera=len(set(g)),lambda_null=lam,nll_lambda0=fit(X0,0)[0],nll_lambda1=fit(X0,1)[0],nll_lambdaMLE=nll0,tree_root='arbitrary; sensitivity only')]).to_csv(O/'PGLS47_diagnostics.tsv',sep='\t',index=False)
print('47 HMT PGLS null lambda',lam,'df',len(y)-p0,'NLL',nll0)
print(pd.DataFrame(out).sort_values('p_nominal').to_string(index=False))
