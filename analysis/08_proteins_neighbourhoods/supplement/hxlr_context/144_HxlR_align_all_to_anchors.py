#!/usr/bin/env python3
"""Quantify nearest experimentally anchored sequence relatives for all 1524 oral HxlR
proteins using local BLOSUM62 alignments; no sensory functional transfer.
This is descriptive complement to the 73-cluster phylogenetic backbone.
"""
from pathlib import Path
import pandas as pd,collections,json
from Bio.Align import PairwiseAligner,substitution_matrices
B=Path('__PROJECT_ROOT__/study13')
S=Path('__WORK_ROOT__/HOMD_HxlR_reference_20261014')
O=B/'16_HxlR_reference_phylogeny_20261014'
h=pd.read_csv(B/'13_HxlR_broad_and_cassette_20261014/HxlR_all_proteins_annotated.tsv',sep='\t')
ref=pd.read_csv(O/'reference_regulators_curated_provenance.tsv',sep='\t')
def readfa(file):
 d={};name=None;seq=[]
 for line in open(file):
  if line.startswith('>'):
   if name:d[name]=''.join(seq)
   name=line[1:].split()[0];seq=[]
  else:seq.append(line.strip())
 if name:d[name]=''.join(seq)
 return d
seqs=readfa(S/'HxlR_all_oral_plus_verified_refs.faa')
refs={str(r.label):seqs['REF__'+str(r.label)] for r in ref.itertuples()}
aligner=PairwiseAligner()
aligner.mode='local'
aligner.substitution_matrix=substitution_matrices.load('BLOSUM62')
aligner.open_gap_score=-10
aligner.extend_gap_score=-0.5
records=[]
for i,entry in enumerate(h.itertuples(),1):
 p=entry.protein_id;query=seqs[p]
 for label,target in refs.items():
  aln=aligner.align(query,target)[0]
  iq=[];it=[]
  for (qa,qb),(ta,tb) in zip(aln.aligned[0],aln.aligned[1]):
   iq.append((qa,qb));it.append((ta,tb))
  matches=aligned=0
  for (qa,qb),(ta,tb) in zip(iq,it):
   fragQ=query[qa:qb];fragT=target[ta:tb]
   matches+=sum(a==b for a,b in zip(fragQ,fragT))
   aligned+=qb-qa
  qcov=aligned/len(query);tcov=aligned/len(target)
  ident=matches/aligned if aligned else 0
  records.append({'protein_id':p,'GCA':entry.GCA,'Genus':entry.Genus,'Species':entry.Species,'HMT':entry.HMT,
       'Ecology':entry.Ecology,'X5_cluster':entry.MMseqs_X5_representative,'reference':label,
       'n_aligned_positions':aligned,'local_identity':round(ident,4),'query_coverage':round(qcov,4),
       'reference_coverage':round(tcov,4),'blosum_score':round(aln.score,2),
       'strong_coverage':qcov>=0.7 and tcov>=0.7})
 if i%300==0:print('PAIRS',i,'/',len(h),flush=True)
tbl=pd.DataFrame(records)
tbl.to_csv(O/'HxlR_1524_protein_reference_all_pairwise_scores.tsv',sep='\t',index=False)
valid=tbl[tbl.strong_coverage].sort_values(['protein_id','local_identity','blosum_score'],ascending=[True,False,False]).drop_duplicates('protein_id')
closest=h[['protein_id','GCA','Genus','Species','HMT','Ecology','MMseqs_X5_representative']].merge(valid[['protein_id','reference','local_identity','query_coverage','reference_coverage','blosum_score']],on='protein_id',how='left')
closest['alignment_confidence']='uncertain_no_strong_coverage'
closest.loc[(closest.local_identity>=0.3)&closest.local_identity.notna(),'alignment_confidence']='family_similarity_ge30pct'
closest.loc[(closest.local_identity>=0.45)&closest.local_identity.notna(),'alignment_confidence']='close_family_similarity_ge45pct'
closest.to_csv(O/'HxlR_1524_oral_closest_reference_cautious.tsv',sep='\t',index=False)
group=closest.groupby(['reference','alignment_confidence'],dropna=False).agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),n_genera=('Genus','nunique')).reset_index()
group.to_csv(O/'HxlR_nearest_reference_family_census.tsv',sep='\t',index=False)
clusters=closest.groupby(['MMseqs_X5_representative','reference'],dropna=False).agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),median_identity=('local_identity','median')).reset_index()
clusters.to_csv(O/'HxlR_X5_closest_anchored_reference.tsv',sep='\t',index=False)
print('FINISHED',len(h),'REFERENCE_GROUPS\n'+group.to_string(index=False),flush=True)
