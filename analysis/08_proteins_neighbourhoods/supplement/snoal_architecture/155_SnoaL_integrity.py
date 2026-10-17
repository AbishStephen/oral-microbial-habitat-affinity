from pathlib import Path
import pandas as pd,numpy as np,collections
B=Path('__PROJECT_ROOT__/study13')
O=B/'18_SnoaL2_function_gene_integrity_20261016';O.mkdir(exist_ok=True)
meta=pd.read_csv(B/'17_SnoaL2_independent_targets_20261016/SnoaL2_all524_domain_architecture_linkage.tsv',sep='\t')
seq={};n=None;buf=[]
for l in open('__WORK_ROOT__/SnoaL2_subfamilies_20261016/SnoaL2_all524_proteins.faa'):
 if l.startswith('>'):
  if n:seq[n]=''.join(buf)
  n=l[1:].strip();buf=[]
 else:buf.append(l.strip())
if n:seq[n]=''.join(buf)
meta['protein_length_aa']=meta.protein_id.map(lambda x:len(seq[x]))
meta.to_csv(O/'SnoaL2_original524_protein_sizes.tsv',sep='\t',index=False)
print('PROTEIN_LENGTH_PER_ARCHITECTURE',flush=True)
print(meta.groupby('architecture').protein_length_aa.describe(percentiles=[.1,.25,.5,.75,.9]).round(1).to_string(),flush=True)
print('LONG_PROTEINS',meta.sort_values('protein_length_aa',ascending=False)[['protein_id','Genus','architecture','protein_length_aa','linked_HxlR_3kb']].head(20).to_string(index=False),flush=True)
def hits(p):
 d={}
 for l in open(p):
  if l.startswith('#'):continue
  x=l.split()
  if len(x)<23:continue
  key=x[0]
  d.setdefault(key,[]).append((int(x[17]),int(x[18]),float(x[12])))
 return d
D=B/'14_HxlR_fullHOMD_controls_20261014'
z1=hits(D/'HOMD631_PF12680_domtblout.tsv');z2=hits(B/'17_SnoaL2_independent_targets_20261016/HOMD631_PF07366_domtblout.tsv')
out=[]
for key in set(z1)&set(z2):
 for s1,e1,p1 in z1[key]:
  for s2,e2,p2 in z2[key]:
   overlap=max(0,min(e1,e2)-max(s1,s2)+1)
   span=min(e1-s1+1,e2-s2+1)
   out.append(dict(protein_id=key,pf12680_start=s1,pf12680_end=e1,pf07366_start=s2,pf07366_end=e2,
      overlap_aa=overlap,overlap_fraction_smaller=overlap/span,
      nonoverlapping=overlap==0))
df=pd.DataFrame(out)
df.to_csv(O/'HOMD631_PF07366_PF12680_domain_overlaps.tsv',sep='\t',index=False)
print('OVERLAP',len(df),'PROTEINS',df.protein_id.nunique(),'SAME_REGION_FRAC_GE_0_5',int((df.overlap_fraction_smaller>=.5).sum()),'NONOVERLAPPING',int(df.nonoverlapping.sum()),flush=True)
print('EXAMPLES',df.head(8).to_string(index=False),flush=True)
