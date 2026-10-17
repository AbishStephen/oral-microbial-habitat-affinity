#!/usr/bin/env python3
import csv,gzip,os
from pathlib import Path
from collections import defaultdict,Counter
B=Path('__PROJECT_ROOT__/study13')
O=B/'09_evolution';O.mkdir(parents=True,exist_ok=True)
targets=['PF04072','PF12680','PF10604','PF08020','PF14534','PF13199','PF22815','PF04239','PF20730','PF04191','PF13333']
with open(B/'02_prokka/master_protein_catalogue.tsv') as f:
 prot={x['protein_id']:(x['GCA'],x['HMT'],x['Genus'],x['Species'],x['product']) for x in csv.DictReader(f,delimiter='\t')}
best={};n_all=Counter();unique=defaultdict(set);gcas=defaultdict(set)
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for z in csv.DictReader(f,delimiter='\t'):
  fam=z['Pfam']
  if fam not in targets:continue
  pid=z['protein_id']
  if pid not in prot:continue
  gca,hmt,genus,species,product=prot[pid]
  key=fam,hmt
  n_all[fam]+=1;unique[fam].add(pid);gcas[fam].add(gca)
  try:score=float(z['domain_score'])
  except:continue
  if key not in best or score>best[key][0]:
   best[key]=(score,pid,int(z['ali_from']),int(z['ali_to']),product,gca,genus,species)
# load peptide FASTA for only target proteins (best hits) from directory
wanted=defaultdict(set)
for (fam,hmt),z in best.items(): wanted[z[5]].add(z[1])
seq={}
def fasta(path,want):
 d={};key=None;buf=[]
 with open(path) as f:
  for line in f:
   if line.startswith('>'):
    if key in want:d[key]=''.join(buf)
    key=line[1:].split()[0];buf=[]
   else:buf.append(line.strip())
  if key in want:d[key]=''.join(buf)
 return d
for gca,want in wanted.items():
 path=B/'02_prokka/faa'/f'{gca}.faa'
 if not path.exists():path=B/'02_prokka/faa'/f'{gca}.faa' # most older in original; fallback
 if not path.exists():
  base=Path('__PROJECT_ROOT__')
  path=base/'02_prokka/faa'/f'{gca}.faa'
 if not path.exists():
  print('MISSING FASTA',gca);continue
 seq.update(fasta(path,want))
rec=[]
for (fam,hmt),v in sorted(best.items()):
 score,pid,lo,hi,product,gca,genus,species=v
 if pid not in seq:continue
 full=seq[pid];dom=full[lo-1:hi]
 if len(dom)<28:continue
 name=f'{hmt}|{genus}_{species}|{pid}'.replace(' ','_')
 rec.append((fam,hmt,genus,species,gca,pid,score,lo,hi,len(full),len(dom),product))
 with open(O/f'{fam}.domains.faa','a') as out:out.write(f'>{name}\n{dom}\n')
with open(O/'domain_phylogeny_catalogue.tsv','w',newline='') as out:
 w=csv.writer(out,delimiter='\t');w.writerow(['Pfam','HMT','Genus','Species','GCA','protein_id','domain_score','ali_start','ali_end','protein_length','domain_length','product']);w.writerows(rec)
with open(O/'domain_phylogeny_summary.tsv','w',newline='') as out:
 w=csv.writer(out,delimiter='\t');w.writerow(['Pfam','total_domain_hits','total_proteins','total_genomes','HMT_representatives_in_domain_tree'])
 for fam in targets:w.writerow([fam,n_all[fam],len(unique[fam]),len(gcas[fam]),sum(z[0]==fam for z in rec)])
print('TARGET PROTEINS',len(rec));print(open(O/'domain_phylogeny_summary.tsv').read())
