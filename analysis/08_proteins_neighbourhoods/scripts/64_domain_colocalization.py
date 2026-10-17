#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import gzip,csv, itertools
B=Path('__PROJECT_ROOT__');E=B/'study13'
O=E/'09_evolution'/'colocalization';O.mkdir(parents=True,exist_ok=True)
families=['PF04072','PF12680','PF10604','PF08020','PF14534','PF13199','PF22815']
hits=defaultdict(lambda:defaultdict(set))
with gzip.open(E/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['Pfam'] in families: hits[r['GCA']][r['Pfam']].add(r['protein_id'])
pairlist=list(itertools.combinations(families,2))
genomes={k:{p for s in dd.values() for p in s} for k,dd in hits.items()}
count=defaultdict(int);examples=defaultdict(list)
for gca,ff in hits.items():
 if len(ff)<2:continue
 gf=E/'02_prokka/gff'/f'{gca}.gff'
 if not gf.exists():gf=B/'02_prokka/gff'/f'{gca}.gff'
 if not gf.exists():continue
 idx={};byctg=defaultdict(list)
 with open(gf,errors='replace') as f:
  for line in f:
   if line.startswith('#'):continue
   p=line.rstrip().split('\t')
   if len(p)<9 or p[2]!='CDS':continue
   attr={z.split('=',1)[0]:z.split('=',1)[1] for z in p[8].split(';') if '=' in z}
   pid=attr.get('ID') or attr.get('locus_tag')
   if pid in genomes[gca]:byctg[p[0]].append((int(p[3]),int(p[4]),pid))
 for ctg,genes in byctg.items():
  for j,(start,end,pid) in enumerate(sorted(genes)):idx[pid]=(ctg,j,start,end)
 # NOTE: byctg currently only target genes; true +/-10 rank requires all CDS. Correct below.
 for x,y in pairlist:
  if x not in ff or y not in ff:continue
  count[(x,y,'co_genomes')]+=1
  if ff[x]&ff[y]:
   count[(x,y,'same_protein')]+=1
  poss=[]
  for a in ff[x]:
   for b in ff[y]:
    if a not in idx or b not in idx:continue
    ca,ja,sa,ea=idx[a];cb,jb,sb,eb=idx[b]
    if ca!=cb:continue
    d=max(0,max(sa,sb)-min(ea,eb))
    poss.append((d,a,b))
  if poss:
   dis,a,b=min(poss)
   count[(x,y,'same_contig')]+=1
   if dis<=10000:count[(x,y,'within_10kb')]+=1
   if dis<=3000:count[(x,y,'within_3kb')]+=1
   if len(examples[(x,y)])<12: examples[(x,y)].append((gca,dis,a,b))
# count all component frequencies genome-wide
for x,y in pairlist:
 for ff in hits.values():
  if x in ff:count[(x,y,'x_genomes')]+=1
  if y in ff:count[(x,y,'y_genomes')]+=1
with open(O/'pairwise_colocalization.tsv','w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['Pfam1','Pfam2','x_genomes','y_genomes','co_genomes','same_protein','same_contig','within_3kb','within_10kb'])
 for x,y in pairlist:w.writerow([x,y]+[count[(x,y,s)] for s in ['x_genomes','y_genomes','co_genomes','same_protein','same_contig','within_3kb','within_10kb']])
with open(O/'example_neighbours.tsv','w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['Pfam1','Pfam2','GCA','distance_bp','protein1','protein2'])
 for (x,y),ex in examples.items():
  for g,d,a,b in ex:w.writerow([x,y,g,d,a,b])
print('TOP3 COLLOC:')
with open(O/'pairwise_colocalization.tsv') as f:
 for line in f:
  if 'PF04072' in line or ('PF12680' in line and 'PF10604' in line):print(line.rstrip())
print('COMPLETE',len(hits),'genomes with targeted hits')
