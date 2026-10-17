#!/usr/bin/env python3
from pathlib import Path
import gzip,csv
from collections import defaultdict,Counter
import re,pandas as pd
B=Path('__PROJECT_ROOT__/study13');O=B/'10_additional_analyses_20261012'
np=(O/'crossgenus_candidate_neighborhood_orientation.tsv')
p=pd.read_csv(np,sep='\t')
d=p[(p.Pfam1=='PF01638') & (p.Pfam2=='PF12680')]
nr=pd.read_csv(O/'crossgenus_candidate_distinct_genes_within3kb_examples.tsv',sep='\t')
h=nr[(nr.Pfam1=='PF01638')&(nr.Pfam2=='PF12680')]
protein_gca=defaultdict(set)
for r in h.itertuples():protein_gca[r.GCA].add(r.protein1)
def iterseq(path):
 name=None;seq=[]
 for line in open(path):
  line=line.strip()
  if not line:continue
  if line.startswith('>'):
   if name:yield name,''.join(seq)
   name=line[1:].split()[0];seq=[]
  else:seq.append(line)
 if name:yield name,''.join(seq)
data=[]
for r in h.itertuples():
 prot=r.protein1
 path=B/'02_prokka/faa'/f'{r.GCA}.faa'
 if not path.exists():continue
 if prot not in protein_gca[r.GCA]:continue
for gca, proteins in protein_gca.items():
 path=B/'02_prokka/faa'/f'{gca}.faa'
 if not path.exists():continue
 for p_id,seq in iterseq(path):
  if p_id not in proteins:continue
  hits=[m.start()+1 for m in re.finditer('C.K',seq[:36])]
  data.append(dict(GCA=gca,protein_id=p_id,Genus=str(h[h.GCA==gca].iloc[0].Genus),len_aa=len(seq),Nterm_26=seq[:26],CxK_positions=';'.join(map(str,hits)),first_35_has_CxK=bool(hits),first20_cys=seq[:20].count('C'),first20_lys=seq[:20].count('K')))
out=pd.DataFrame(data)
out.to_csv(O/'HxlR_SnoaL2_linked_Nterminal_CxK_motif_audit.tsv',sep='\t',index=False)
print('HXL_CXK_CANDIDATES',len(out),'WITH CxK near Nterm',int(out.first_35_has_CxK.sum()))
print('BY GENUS',out.groupby('Genus')[['len_aa','first_35_has_CxK']].agg(['median','mean']).round(3).to_string())
print('EXAMPLES',out.drop_duplicates('Genus')[['Genus','len_aa','Nterm_26','CxK_positions']].to_string(index=False))
