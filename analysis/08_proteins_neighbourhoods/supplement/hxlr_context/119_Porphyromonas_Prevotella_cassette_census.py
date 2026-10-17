#!/usr/bin/env python3
"""Provenance-aware Porphyromonas / Prevotella HxlR/SnoaL locus homology census.
Reference shared tract ~1789 nt GCA_002892595.1; screen DNA in original
assemblies and include all species/HMTs in 13-genus manifest.
Exact-seed matches are discovery screens, not stand-in for full ANI/phylogeny."""
from pathlib import Path
import csv,json,collections
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'13_HxlR_broad_and_cassette_20261014';O.mkdir(exist_ok=True)
OLD=B/'12_HxlR_SnoaL_regulatory_genomics_20261013'
S=Path('__WORK_ROOT__/HOMD_convergence_expanded13_work/raw_genomes')
MM=B/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv'
df=pd.read_csv(MM,sep='\t');t=df[df.Genus.isin(['Porphyromonas','Prevotella'])].copy()
hx=pd.read_csv(O/'HxlR_all_proteins_annotated.tsv',sep='\t')
d=pd.read_csv(B/'10_additional_analyses_20261012/crossgenus_candidate_Pfam_genomic_coordinates.tsv',sep='\t')
lk=pd.read_csv(OLD/'linked_locus_audit.tsv',sep='\t')
def readfa(p):
 out={};n=None;buf=[]
 for z in open(p):
  if z.startswith('>'):
   if n:out[n]=''.join(buf).upper()
   n=z[1:].split()[0];buf=[]
  else:buf.append(z.strip())
 if n:out[n]=''.join(buf).upper()
 return out
def rev(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
r=lk[lk.GCA=='GCA_002892595.1'].iloc[0]
refg=r.GCA;refseqs=readfa(S/f'{refg}.fna')
refcontig=refseqs[r.accession_contig]
interstart=int(r.start_1based)-1;interend=int(r.end_1based)
reference=refcontig[interstart-849:interend+808]
assert len(reference)==1789,len(reference)
# fragment until full exact tracts are characterized. 45nt anchors distributed every ~150bp.
anchors=[(i,reference[i:i+45]) for i in range(0,len(reference)-44,150)]
revref=rev(reference)
result=[]; missing=[]
for ii,z in enumerate(t.itertuples(index=False),1):
 g=z.GCA;paths=[S/f'{g}.fna',OLD/'genomes'/f'{g}.fna',O/'genome_fasta'/f'{g}.fna']
 f=next((q for q in paths if q.exists()),None)
 if f is None:
  missing.append({'GCA':g,'Genus':z.Genus,'Species':z.Species,'HMT':z.HMT,'reason':'FNA not archived'})
  result.append(dict(GCA=g,Genus=z.Genus,Species=z.Species,HMT=z.HMT,has_genome_fna=False,
                     full_1789nt_exact=False,anchor45_nt_hits=0,region_best_contig='',region_direction=''))
  continue
 contigs=readfa(f)
 hits=0; full=False;top='';direc='';best=0
 for contig,seq in contigs.items():
  if reference in seq or revref in seq:
   full=True;top=contig;direc='+' if reference in seq else '-';hits=len(anchors);break
  for strand,q in [('+',seq),('-',rev(seq))]:
   h=sum(fragment in q for offset,fragment in anchors)
   if h>best:best=h;top=contig;direc=strand
 if not full:hits=best
 result.append(dict(GCA=g,Genus=z.Genus,Species=z.Species,HMT=z.HMT,has_genome_fna=True,full_1789nt_exact=full,
   anchor45_nt_hits=hits,region_best_contig=top,region_direction=direc))
 if ii%30==0:print('COMPARED',ii,'/',len(t),flush=True)
survey=pd.DataFrame(result)
hxcounts=hx.groupby('GCA').size().rename('n_HxlR_PF01638')
sncount=d[d.Pfam=='PF12680'].groupby('GCA').protein_id.nunique().rename('n_SnoaL2_PF12680')
survey=survey.merge(hxcounts,on='GCA',how='left').merge(sncount,on='GCA',how='left')
survey[['n_HxlR_PF01638','n_SnoaL2_PF12680']]=survey[['n_HxlR_PF01638','n_SnoaL2_PF12680']].fillna(0).astype(int)
survey['HxlR_SnoaL_linked']=survey.GCA.isin(lk.GCA)
survey.to_csv(O/'Porphyromonas_Prevotella_all243_cassette_nucleotide_screen.tsv',sep='\t',index=False)
summary=survey.groupby(['Genus','HMT','Species']).agg(n_genomes=('GCA','size'),n_with_nucleotide=('has_genome_fna','sum'),
 exact_1789nt=('full_1789nt_exact','sum'),with_any_anchor45=('anchor45_nt_hits',lambda x:int((x>0).sum())),
 positive_HxlR=('n_HxlR_PF01638',lambda x:int((x>0).sum())),
 positive_SnoaL2=('n_SnoaL2_PF12680',lambda x:int((x>0).sum())),
 linked_HxlR_SnoaL2=('HxlR_SnoaL_linked','sum')).reset_index()
summary.to_csv(O/'Porphyromonas_Prevotella_nucleotide_homology_by_species.tsv',sep='\t',index=False)
with open(O/'shared_1789nt_reference.fna','w') as f:f.write(f'>{refg}|HxlR-SnoaL2|shared_segment\n'+reference+'\n')
(O/'shared_1789nt_reference_metadata.json').write_text(json.dumps({'reference_gca':refg,'contig':r.accession_contig,'segment_1based_start':interstart-849+1,'segment_1based_end':interend+808,'length_bp':1789,'anchors_45nt':len(anchors),'reference_scope':'Published shared exact tract; check all homologous/divergent segments separately'},indent=2))
print('SUMMARY\n'+summary.to_string(index=False),flush=True)
print('TOTAL',len(survey),'available_fna',int(survey.has_genome_fna.sum()),'exact',int(survey.full_1789nt_exact.sum()),'anchor_positive',int((survey.anchor45_nt_hits>0).sum()),flush=True)
