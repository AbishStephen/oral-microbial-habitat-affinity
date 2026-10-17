#!/usr/bin/env python3
import csv,collections
from pathlib import Path
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012'
m=list(csv.DictReader(open(B/'08_validation/confound_audit/genome_metadata_provenance.tsv'),delimiter='\t'))
by={x['GCA']:(x['HMT'],x['Genus']) for x in m}
r=list(csv.DictReader(open(B/'09_evolution/colocalization/example_neighbours.tsv'),delimiter='\t'))
out=[];s=collections.defaultdict(list)
for x in r:
 if x['Pfam1']=='PF12680' and x['Pfam2']=='PF10604':
  H,G=by[x['GCA']];s[G].append(x)
for G,v in sorted(s.items(),key=lambda z:-len(z[1])):
 out.append(dict(Genus=G,n_genomes=len(v),n_HMT=len(set(by[z['GCA']][0] for z in v)),median_separation_bp=sorted(int(z['distance_bp']) for z in v)[len(v)//2]))
print('NEIGHBOR_GENE_ALLOCATION')
for x in out:print(x)
with open(O/'PF12680_PF10604_linked_genus_summary.tsv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
