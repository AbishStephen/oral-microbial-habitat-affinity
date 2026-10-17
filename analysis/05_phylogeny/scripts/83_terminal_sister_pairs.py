#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
s=(B/'07_phylogeny/HOMD47_BS200_SUP.raxml.support').read_text().strip().rstrip(';')
class N:
 def __init__(self,kids,name,length):
  self.kids=kids;self.name=name;self.length=length
def parse(s):
 i=0
 def node():
  nonlocal i
  kids=[]
  if s[i]=='(':
   i+=1
   while True:
    kids.append(node())
    if s[i]==',':i+=1;continue
    assert s[i]==')',s[i:i+20]
    i+=1;break
  name=''
  while i<len(s) and s[i] not in ':,()':
   name+=s[i];i+=1
  length=0.
  if i<len(s) and s[i]==':':
   i+=1;start=i
   while i<len(s) and s[i] not in ',()':i+=1
   length=float(s[start:i])
  return N(kids,name,length)
 return node()
root=parse(s)
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT')
pairs=[]
def walk(n):
 if len(n.kids)==2 and all(not k.kids for k in n.kids):
  a,b=(k.name.split('|')[0] for k in n.kids)
  if a in d.index and b in d.index:
   x,y=d.loc[a],d.loc[b]
   pairs.append(dict(HMT1=a,HMT2=b,genus_1=x.Genus,genus_2=y.Genus,
     Y1=x.Y,Y2=y.Y,Ecology1=x.Ecology,Ecology2=y.Ecology,
     support=float(n.name) if n.name.replace('.','',1).isdigit() else np.nan,
     branch1=n.kids[0].length,branch2=n.kids[1].length,
     discordant=x.Ecology!=y.Ecology))
 for k in n.kids:walk(k)
walk(root)
allpair=pd.DataFrame(pairs)
print('TOTAL_TERMINAL_CHERRIES',len(allpair),'DISCORDANT',sum(allpair.discordant))
print(allpair[['HMT1','HMT2','genus_1','genus_2','Ecology1','Ecology2','Y1','Y2','support','discordant']].to_string(index=False))
out=B/'10_additional_analyses_20261012';out.mkdir(exist_ok=True)
allpair.to_csv(out/'HOMD47_terminal_sister_pairs_bootstrap200.tsv',sep='\t',index=False)
