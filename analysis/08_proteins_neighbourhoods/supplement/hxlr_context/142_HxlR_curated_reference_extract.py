#!/usr/bin/env python3
"""Build small, provenance-controlled experimentally anchored HxlR/DUF24 reference dataset;
reuse existing 1,524 oral proteins. Bulk FASTA and trees to scratch only."""
from pathlib import Path
import pandas as pd,urllib.request,json,csv,collections,time,hashlib
B=Path('__PROJECT_ROOT__/study13')
D=B/'13_HxlR_broad_and_cassette_20261014'
O=B/'16_HxlR_reference_phylogeny_20261014';O.mkdir(parents=True,exist_ok=True)
S=Path('__WORK_ROOT__/HOMD_HxlR_reference_20261014');S.mkdir(parents=True,exist_ok=True)
refs=[
 ('Bsub_HxlR','P42406','Bacillus subtilis 168','formaldehyde','hxlAB','experimental','15978081;26973631'),
 ('Cglu_QorR','Q8NQQ5','Corynebacterium glutamicum ATCC13032','diamide,H2O2,cumene hydroperoxide','qor2','candidate_crosswalk','19403527'),
 ('Bsub_YodB','O34844','Bacillus subtilis 168','quinone,diamide,catechol','azoR1,yodC,mhq','experimental',''),
 ('Bsub_HypR_YybR','P37486','Bacillus subtilis 168','diamide,hypochlorite','hypO','experimental',''),
 ('Saureus_QsrR','A0A0H3JWI0','Staphylococcus aureus Mu50','quinone,menadione','catE2,azoR1,yodC','experimental','23479646'),
 ('Bsub_YtcD','O34533','Bacillus subtilis 168','unknown','unknown','predicted_regulog',''),
 ('Bsub_YdeP','P96673','Bacillus subtilis 168','unknown','unknown','predicted_regulog',''),
]
def fetcha(acc):
 url=f'https://rest.uniprot.org/uniprotkb/{acc}.fasta'
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'HOMD-HxlR-phylogeny/1.0'}),timeout=30) as h:raw=h.read().decode()
 seq=''.join(z.strip() for z in raw.splitlines() if not z.startswith('>'))
 if len(seq)<70 or len(seq)>400 or any(x not in 'ACDEFGHIKLMNPQRSTVWYUXO' for x in seq):raise ValueError(acc+' invalid sequence')
 return seq,raw.splitlines()[0]
records=[];refseq={}
for short,acc,org,eff,target,status,pmid in refs:
 seq,desc=fetcha(acc);refseq[short]=seq
 records.append(dict(label=short,accession=acc,organism=org,effectors_from_literature=eff,
  known_or_inferred_targets=target,evidence_class=status,pubmed=pmid,length=len(seq),
  first_35aa=seq[:35],cysteines=','.join(str(i+1) for i,a in enumerate(seq) if a=='C'),
  accession_description=desc))
 print('REF',short,acc,len(seq),[(i+1) for i,x in enumerate(seq) if x=='C'],flush=True)
pd.DataFrame(records).to_csv(O/'reference_regulators_curated_provenance.tsv',sep='\t',index=False)
def readfa(path,lookup):
 ids={x for x in lookup}
 seq={};name=None, 
 name=None;chunks=[]
 for line in open(path):
  if line.startswith('>'):
   if name in ids:seq[name]=''.join(chunks)
   name=line[1:].split()[0];chunks=[]
  else:chunks.append(line.strip())
 if name in ids:seq[name]=''.join(chunks)
 return seq
df=pd.read_csv(D/'HxlR_all_proteins_annotated.tsv',sep='\t')
oral={};missing=[]
for g,subset in df.groupby('GCA'):
 pidlist=set(subset.protein_id)
 seq=readfa(B/'02_prokka/faa'/f'{g}.faa',pidlist)
 for pid in pidlist:
  if pid in seq:oral[pid]=seq[pid]
  else:missing.append((g,pid))
assert not missing,(len(missing),missing[:3])
assert len(oral)==1524,len(oral)
label_map={short:'REF__'+short for short in refseq}
with open(S/'HxlR_all_oral_plus_verified_refs.faa','w') as f:
 for short,seq in refseq.items():f.write('>'+label_map[short]+'\n'+seq+'\n')
 for pid,seq in oral.items():f.write('>'+pid+'\n'+seq+'\n')
# use MMseqs X5 representatives: one per established family cluster
clusters=df.groupby('MMseqs_X5_representative').protein_id.first().to_dict()
rep_ids={cluster: (cluster if cluster in oral else pid) for cluster,pid in clusters.items()}
with open(S/'HxlR_X5_cluster_reps_plus_refs.faa','w') as f:
 for short,seq in refseq.items():f.write('>'+label_map[short]+'\n'+seq+'\n')
 for cluster,pid in rep_ids.items():f.write('>'+pid+'\n'+oral[pid]+'\n')
assert len(rep_ids)==73,len(rep_ids)
pd.DataFrame([dict(MMseqs_X5_representative=c,protein_id=p) for c,p in rep_ids.items()]).to_csv(O/'oral_73_X5_family_representative_selection.tsv',sep='\t',index=False)
data={'oral_proteins':len(oral),'oral_X5_clusters':len(rep_ids),'reference_proteins':len(refseq),
 'reference_labels':list(refseq),'scratch_dir':str(S),'compact_dir':str(O),
 'notes':'Q8NQQ5 (Cgl1371/cg1552) is provisionally crosswalked; verify Cys17 and genomic coordinates against Ehira 2009 before classifying as experimentally anchored.'}
(O/'run_manifest.json').write_text(json.dumps(data,indent=2))
print('FINISH',json.dumps(data),flush=True)
