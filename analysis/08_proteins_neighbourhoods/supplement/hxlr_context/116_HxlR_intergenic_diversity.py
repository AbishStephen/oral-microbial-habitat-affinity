#!/usr/bin/env python3
"""Descriptive regulatory-region diversity and HxlR protein evidence; not a promoter caller."""
from pathlib import Path
import collections, csv, json,statistics,itertools,random
B=Path('__PROJECT_ROOT__/study13')
O=B/'12_HxlR_SnoaL_regulatory_genomics_20261013'
def readfa(path):
 a={};name=None;buf=[]
 for line in open(path):
  if line.startswith('>'):
   if name:a[name]=''.join(buf).upper()
   name=line[1:].rstrip();buf=[]
  else:buf.append(line.strip())
 if name:a[name]=''.join(buf).upper()
 return a
def comp(seq):return seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
def table(name,records):
 if not records:return
 with open(O/name,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
meta=list(csv.DictReader(open(O/'linked_locus_audit.tsv'),delimiter='\t'))
seqs=readfa(O/'HxlR_SnoaL2_intergenic_all_available.fna')
assert len(meta)==67 and len(seqs)==67
by_g=collections.defaultdict(list)
for title,seq in seqs.items():
 g,genus,hmt=title.split('|')[:3]
 by_g[genus].append((g,seq,hmt))
groups=[];reps={}
for genus,arr in sorted(by_g.items()):
 coll=collections.Counter(s for g,s,hmt in arr)
 for seq,n in coll.most_common():
  if seq not in reps:reps[seq]=set()
  reps[seq].add(genus)
 rows={'Genus':genus,'n_genomes':len(arr),'n_distinct_intergenic_sequences':len(coll),
       'most_common_sequence_count':coll.most_common(1)[0][1], 'gap_lengths':','.join(str(z) for z in sorted({len(s) for _,s,_ in arr})),
       'n_HMTs':len(set(hmt for _,_,hmt in arr))}
 groups.append(rows)
table('intergenic_diversity_by_genus.tsv',groups)
# collapse identical DNA within each genus, require support in >= 2 distinct genus families
rep_per_genus={g:list(set(seq for gc,seq,h in arr)) for g,arr in by_g.items() if g!='Neisseria'}
motifs=[]
for size in [8,9,10,11,12]:
 count=collections.defaultdict(set)
 for genus,uniq in rep_per_genus.items():
  km=set()
  for seq in uniq:
   for i in range(len(seq)-size+1):
    word=seq[i:i+size]
    if 'N' not in word:km.add(min(word,comp(word)))
  for x in km:count[x].add(genus)
 for seq,gs in count.items():
  if len(gs)>=3:motifs.append({'k':size,'canonical_sequence':seq,'n_genera':len(gs),'genera':';'.join(sorted(gs))})
motifs=sorted(motifs,key=lambda z:(-z['n_genera'],-z['k'],z['canonical_sequence']))
table('exploratory_shared_8to12mers.tsv',motifs)
# Store experimentally undefined; no claimed statistical enrichment
proteins=[]
for genus,arr in by_g.items():
 for g,_,_ in arr:
  record=next(x for x in meta if x['GCA']==g)
  raw=readfa(B/'02_prokka/faa'/f'{g}.faa')
  def seq_for(pid):
   matches=[s for hd,s in raw.items() if hd.split()[0]==pid]
   return matches[0] if matches else ''
  hx=seq_for(record['HxlR_protein'])
  sn=seq_for(record['SnoaL_protein'])
  proteins.append({'GCA':g,'Genus':genus,'HMT':record['HMT'],'HxlR_id':record['HxlR_protein'],
      'HxlR_length':len(hx),'HxlR_CxK_first35':any(hx[i]=='C' and hx[i+2]=='K' for i in range(min(33,len(hx)-2))) if hx else None,
      'SnoaL2_id':record['SnoaL_protein'],'SnoaL2_length':len(sn),
      'HxlR_seq':hx,'SnoaL2_seq':sn})
table('linked_protein_lengths_and_motifs.tsv',[{k:v for k,v in row.items() if not k.endswith('_seq')} for row in proteins])
for gene,key in [('HxlR','HxlR_seq'),('SnoaL2','SnoaL2_seq')]:
 with open(O/f'linked_{gene}_proteins.faa','w') as f:
  for row in proteins:
   seq=row[key]
   if seq:f.write('>'+row['GCA']+'|'+row['Genus']+'|'+row['HMT']+'\n'+seq+'\n')
summ={'n_genomes':len(meta),'distinct_intergenic_sequences':len(reps),'protein_sequences_HxlR':sum(bool(x['HxlR_seq']) for x in proteins),
 'protein_sequences_SnoaL2':sum(bool(x['SnoaL2_seq']) for x in proteins),
 'n_shared_identical_sequences_across_genera':sum(len(gs)>1 for gs in reps.values()),
 'shared_short_8to12mer_candidates_genera3plus':len(motifs),
 'caution':'Short k-mers are not validated promoter sequences, motifs, or HxlR DNA-binding sites; genus-aware null tests and independent negative controls required.'}
(O/'sequence_diversity_summary.json').write_text(json.dumps(summ,indent=2))
print(json.dumps(summ,indent=2));print('GENERA',groups);print('TOP_MOTIFS',motifs[:12])
