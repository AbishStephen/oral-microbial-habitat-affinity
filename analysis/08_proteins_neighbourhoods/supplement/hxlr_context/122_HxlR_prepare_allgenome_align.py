#!/usr/bin/env python3
"""Prepare complete 243-genome nucleotide alignment and PF01638 protein query.
Run on the study workflow compute node, not login."""
from pathlib import Path
import csv,pandas as pd,json, collections
B=Path('__PROJECT_ROOT__/study13')
O=B/'13_HxlR_broad_and_cassette_20261014'
S=Path('__WORK_ROOT__/HOMD_convergence_expanded13_work/raw_genomes')
OLD=B/'12_HxlR_SnoaL_regulatory_genomics_20261013/genomes'
A=pd.read_csv(O/'Porphyromonas_Prevotella_all243_cassette_nucleotide_screen.tsv',sep='\t')
h=pd.read_csv(O/'HxlR_all_proteins_annotated.tsv',sep='\t')
subset=h[h.Genus.isin(['Porphyromonas','Prevotella'])]
dirs=(S,OLD,O/'genome_fasta')
outnuc=O/'P_Prev_all243_contigs.fna'
maprows=[]
with open(outnuc,'w') as fout:
 for j,z in enumerate(A.itertuples(index=False),1):
  g=z.GCA
  path=next((d/f'{g}.fna' for d in dirs if (d/f'{g}.fna').exists()),None)
  if path is None:raise RuntimeError('missing FASTA '+g)
  with open(path) as f:
   for line in f:
    if line.startswith('>'):
     raw=line[1:].split()[0]
     newid=g+'__'+raw.replace('|','_').replace(' ','_')
     maprows.append({'contig_combined':newid,'GCA':g,'original_contig':raw,'Species':z.Species,'Genus':z.Genus})
     fout.write('>'+newid+'\n')
    else:fout.write(line)
  if j%50==0:print('COMBINED_GENOMES',j,flush=True)
pd.DataFrame(maprows).to_csv(O/'P_Prev_all243_contig_dictionary.tsv',sep='\t',index=False)
prows=[]
for g,sub in subset.groupby('GCA'):
 ids=set(sub.protein_id)
 fn=B/'02_prokka/faa'/f'{g}.faa'
 with open(fn) as f:
  name=None; seq=[]
  for line in f:
   if line.startswith('>'):
    if name in ids:prows.append((name,''.join(seq)))
    name=line[1:].split()[0];seq=[]
   else:seq.append(line.strip())
  if name in ids:prows.append((name,''.join(seq)))
with open(O/'P_Prev_HxlR_all166_proteins.faa','w') as f:
 for pid,seq in prows:f.write('>'+pid+'\n'+seq+'\n')
pd.DataFrame({'protein_id':[z[0] for z in prows],'length':[len(z[1]) for z in prows]}).to_csv(O/'P_Prev_HxlR_fasta_audit.tsv',sep='\t',index=False)
assert len(prows)==len(subset),('not all proteins',len(prows),len(subset))
print('RESULT',len(A),'genomes',len(maprows),'contigs',len(prows),'HxlR proteins','nt_file_bytes',outnuc.stat().st_size,flush=True)
