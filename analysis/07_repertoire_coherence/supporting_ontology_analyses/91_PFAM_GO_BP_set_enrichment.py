#!/usr/bin/env python3
"""GO biological-process enrichment from Pfam2GO using GO-basic ancestry; exact genus-matched sign flips."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import csv, re
import numpy as np,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012';G=O/'ontology'
go_bp=[
 'GO:0005975', # carbohydrate metabolic
 'GO:0006520', # cellular amino acid metabolic
 'GO:0006091', # precursor metabolites and energy
 'GO:0009058', # biosynthetic process
 'GO:0009056', # catabolic process
 'GO:0055085', # transmembrane transport
 'GO:0006810', # transport
 'GO:0007155', # cell adhesion
 'GO:0071554', # cell wall organization
 'GO:0006950', # response to stress
 'GO:0006281', # DNA repair
 'GO:0006310', # DNA recombination
 'GO:0042710', # biofilm formation
 'GO:0007165', # signal transduction
 'GO:0006355', # regulation of DNA-templated transcription
 'GO:0006629', # lipid metabolic
 'GO:0006508', # proteolysis
 'GO:0051607', # defense response to virus
 'GO:0006260', # DNA replication
 'GO:0006811', # ion transport
 'GO:0055114', # oxidation-reduction
 'GO:0006468', # protein phosphorylation
 'GO:0045333', # cellular respiration
 'GO:0015980', # energy derivation by oxidation
]
parents={};names={};alt={}
term=None
def commit(t):
 if t and 'id' in t:
  parents[t['id']]=t['p'];names[t['id']]=t.get('name','')
  for a in t.get('alt',[]):alt[a]=t['id']
with (G/'go-basic.obo').open() as f:
 for line in f:
  line=line.strip()
  if line=='[Term]':commit(term);term={'p':[],'alt':[]};continue
  if line.startswith('[') and term:commit(term);term=None;continue
  if not term:continue
  if line.startswith('id: '):term['id']=line[4:]
  if line.startswith('name: '):term['name']=line[6:]
  if line.startswith('alt_id: '):term['alt'].append(line[8:])
  if line.startswith('is_a: '):term['p'].append(line[6:].split()[0])
  if line.startswith('relationship: part_of '):term['p'].append(line[22:].split()[0])
commit(term)
@lru_cache(None)
def ancestors(id):
 id=alt.get(id,id)
 result={id}
 for x in parents.get(id,[]):
  if x!=id:result|=ancestors(x)
 return frozenset(result)
pfam_go={}
with (G/'pfam2go').open() as f:
 for line in f:
  if not line.startswith('Pfam:PF'):continue
  mm=re.search(r'Pfam:(PF\d+).*?; (GO:\d+)',line)
  if mm:pfam_go.setdefault(mm.group(1),set()).add(mm.group(2))
with open(B/'handoff/X_primary_transferable_DEDUP.tsv') as f:
 rd=csv.reader(f,delimiter='\t');hd=next(rd);r=list(rd)
h=[z[0] for z in r]
ii={z:i for i,z in enumerate(h)}
j=np.flatnonzero([v.startswith('PFAM:') for v in hd[5:]])
feature=np.array(hd[5:],object)[j]
matrix=np.array([[float(z[5+k]) for k in j] for z in r])
maprows=[];M=np.zeros((len(feature),len(go_bp)))
for i,v in enumerate(feature):
 base=v.split(':',1)[1]
 A=set()
 for go in pfam_go.get(base,[]):A.update(ancestors(go))
 for k,goterm in enumerate(go_bp):
  if goterm in A:M[i,k]=1
  if goterm in A:maprows.append((v,goterm))
count=M.sum(0)
S=matrix@M/np.maximum(count,1)
d=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t').set_index('HMT').loc[h]
g=pd.get_dummies(d.Genus,drop_first=True,dtype=float).to_numpy()
cov=d[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float);cov=(cov-cov.mean(0))/(cov.std(0)+1e-12)
Z=np.column_stack([np.ones(len(d)),g,cov]);q=np.linalg.qr(Z)[0]
Sadj=S-q@(q.T@S)
keep=(count>=8)&(Sadj.std(0)>1e-6)
go_bp=np.array(go_bp)[keep].tolist();count=count[keep];S=S[:,keep];Sadj=Sadj[:,keep]
name=[names.get(x,'Unknown') for x in go_bp]
pd.DataFrame(S,columns=go_bp,index=h).rename_axis('HMT').to_csv(O/'PFAM_GO_BP_category_scores_HMT47.tsv',sep='\t')
pd.DataFrame(Sadj,columns=go_bp,index=h).rename_axis('HMT').to_csv(O/'PFAM_GO_BP_adjusted_category_scores_HMT47.tsv',sep='\t')
sister=pd.read_csv(O/'sister8_pairwise_HMT_genome_domain_contrasts.tsv',sep='\t')
extreme=pd.read_csv(O/'extreme13_genus_matched_HMT_pairs.tsv',sep='\t')
out=[]
for panel,P in [('sister8',sister),('extreme13',extreme),('extreme_strong_2p5',extreme[extreme.both_strong_2p5])]:
 for condition,T in [('raw',S),('genus_size_GC_CDS_adjusted',Sadj)]:
  D=T[[ii[z] for z in P.plaque_HMT]]-T[[ii[z] for z in P.mucosal_HMT]]
  norm=np.linalg.norm(D,axis=0);stat=D.sum(0)/np.maximum(norm,1e-15)
  sign=np.array(list(product([-1.,1.],repeat=len(D))))
  pp=(sign@D)/np.maximum(norm,1e-15)
  pv=(np.abs(pp)>=np.abs(stat)[None,:]-1e-12).mean(0)
  mx=np.abs(pp).max(axis=1)
  fw=(mx[:,None]>=np.abs(stat)[None,:]-1e-12).mean(0)
  for i,v in enumerate(go_bp):
   out.append(dict(panel=panel,condition=condition,GO=v,process=name[i],n_Pfam_families=int(count[i]),
    n_pairs=len(P),plaque_higher=int(sum(D[:,i]>1e-8)),mucosal_higher=int(sum(D[:,i]<-1e-8)),
    mean_prevalence_delta=float(D[:,i].mean()),exact_two_sided_p=float(pv[i]),exact_maxT_FWER_p=float(fw[i])))
r=pd.DataFrame(out)
r.to_csv(O/'PFAM_GO_BP_exact_set_enrichment.tsv',sep='\t',index=False)
for label in ['sister8','extreme13']:
 a=r[(r.panel==label)&(r.condition=='genus_size_GC_CDS_adjusted')].sort_values('exact_two_sided_p')
 print('\n',label,'GO_BP RESULTS\n',a.to_string(index=False),flush=True)
print('GO_TERMS_TESTED',len(go_bp),'PFAM_MATRIX',matrix.shape,'ANNOTATED_PFAM',len(set(x for x,_ in maprows)),flush=True)
