#!/usr/bin/env python3
"""Domain-to-protein and MMseqs-X5-family audit for 14 cross-genus candidates;
identifies whether recurrent Pfams correspond to homologous protein clusters.
"""
from pathlib import Path
import csv,gzip,collections,pandas as pd
B=Path('__PROJECT_ROOT__/study13'); O=B/'11_manuscript_completion_20261013'
A=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_feature_integrated_raw_adjusted_annotated.tsv',sep='\t')
C=A[A.support_selection==True].copy();targets={x.split(':')[1] for x in C.feature if x.startswith('PFAM:')}
targets|={'PF12680','PF07366','PF10604'}
protein_hits=collections.defaultdict(set); protinfo={}; protein_hits_count=0
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['Pfam'] in targets:
   p=r['protein_id'];protein_hits[p].add(r['Pfam']);protinfo[p]=(r['GCA'],r['HMT']);protein_hits_count+=1
print('PASS1',len(protein_hits),'protein ids',protein_hits_count,'hits',flush=True)
protein_arch=collections.defaultdict(set)
with gzip.open(B/'04_domains/Pfam_domain_hits.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['protein_id'] in protinfo:protein_arch[r['protein_id']].add(r['Pfam'])
print('PASS2',flush=True)
p_to_rep={}
with open(B/'03_orthogroups/mmseqs_x5/HOMD13_x5_cluster.tsv') as f:
 for line in f:
  z=line.rstrip().split('\t')
  if len(z)>1 and z[1] in protinfo:p_to_rep[z[1]]=z[0]
print('MMSEQS_MATCH',len(p_to_rep),len(protinfo),flush=True)
gcagenus={}
for r in csv.DictReader(open(B/'08_validation/confound_audit/genome_metadata_provenance.tsv'),delimiter='\t'):
 gcagenus[r['GCA']]=r['Genus']
co=collections.defaultdict(lambda:collections.Counter())
by=collections.defaultdict(list)
for pr,hits in protein_hits.items():
 gca,hmt=protinfo[pr]; genus=gcagenus.get(gca,'UNKNOWN');cluster=p_to_rep.get(pr,'UNMAPPED')
 for pf in hits:
  rec=dict(Pfam=pf,protein_id=pr,GCA=gca,HMT=hmt,Genus=genus,MMseqs_X5_representative=cluster,all_Pfam_domains=';'.join(sorted(protein_arch.get(pr,set()))),n_domain_families=len(protein_arch.get(pr,set())))
  by[pf].append(rec)
flat=[r for v in by.values() for r in v]
pd.DataFrame(flat).to_csv(O/'candidate_Pfam_fullprotein_orthogroup_architecture.tsv',sep='\t',index=False)
summary=[];clrows=[]
for pf,rr in sorted(by.items()):
 t=pd.DataFrame(rr)
 h=collections.defaultdict(list)
 for r in rr:h[r['MMseqs_X5_representative']].append(r)
 no_unmapped={k:v for k,v in h.items() if k!='UNMAPPED'}
 support=sorted([(k,len(v),len({z['Genus'] for z in v}),len({z['HMT'] for z in v})) for k,v in no_unmapped.items()],key=lambda v:(-v[2],-v[3],-v[1]))
 for rep,n,ng,nh in support:
  clrows.append(dict(Pfam=pf,cluster=rep,n_proteins=n,n_genera=ng,n_HMT=nh,Genus_list=';'.join(sorted({z['Genus'] for z in h[rep]}))))
 multi=[(x,y) for x,y in t.groupby('all_Pfam_domains').size().items() if x!=pf]
 summary.append(dict(Pfam=pf,n_proteins=len(t),n_genomes=t.GCA.nunique(),n_HMT=t.HMT.nunique(),n_genera=t.Genus.nunique(),
  n_distinct_MMseqs_X5_clusters=len(no_unmapped),max_genera_in_a_single_cluster=support[0][2] if support else 0,
  n_clusters_reaching_ge3_genera=sum(s[2]>=3 for s in support),
  fraction_proteins_single_Pfam=(t.n_domain_families==1).mean(),n_distinct_domain_architectures=t.all_Pfam_domains.nunique(),
  top_cluster_id=support[0][0] if support else 'NA',top_cluster_genus_count=support[0][2] if support else 0,
  top_architecture=';'.join(t.all_Pfam_domains.value_counts().index[:4])))
S=pd.DataFrame(summary).sort_values('Pfam')
S.to_csv(O/'candidate_Pfam_protein_subfamily_summary.tsv',sep='\t',index=False)
pd.DataFrame(clrows).to_csv(O/'candidate_Pfam_MMseqs_family_crossgenus.tsv',sep='\t',index=False)
print(S.to_string(index=False),flush=True)
