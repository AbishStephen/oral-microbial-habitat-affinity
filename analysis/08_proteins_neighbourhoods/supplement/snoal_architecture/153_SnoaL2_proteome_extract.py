#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'17_SnoaL2_independent_targets_20261016'; S=Path('__WORK_ROOT__/SnoaL2_subfamilies_20261016');S.mkdir(exist_ok=True)
a=pd.read_csv(O/'SnoaL2_all524_domain_architecture_linkage.tsv',sep='\t')
found=set()
with open(S/'SnoaL2_all524_proteins.faa','w') as w:
 for g,grp in a.groupby('GCA'):
  ids=set(grp.protein_id)
  n=None;buf=[]
  for line in open(B/'02_prokka/faa'/f'{g}.faa'):
   if line.startswith('>'):
    if n in ids:w.write('>'+n+'\n'+''.join(buf)+'\n');found.add(n)
    n=line[1:].split()[0];buf=[]
   else:buf.append(line.strip())
  if n in ids:w.write('>'+n+'\n'+''.join(buf)+'\n');found.add(n)
assert len(found)==524,(len(found),524)
print('FASTA',len(found),S/'SnoaL2_all524_proteins.faa')
