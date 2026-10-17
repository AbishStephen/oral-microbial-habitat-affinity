#!/usr/bin/env python3
"""Recover Prokka-renamed contigs by exact contig length and validate against translated loci."""
from pathlib import Path
import csv,json,collections,statistics
B=Path('__PROJECT_ROOT__/study13')
O=B/'12_HxlR_SnoaL_regulatory_genomics_20261013'
S=Path('__WORK_ROOT__/HOMD_convergence_expanded13_work/raw_genomes')
rows={r['GCA']:r for r in csv.DictReader(open(B/'10_additional_analyses_20261012/crossgenus_candidate_distinct_genes_within3kb_examples.tsv'),delimiter='\t') if r['Pfam1']=='PF01638' and r['Pfam2']=='PF12680'}
missing=list(csv.DictReader(open(O/'missing_genomes.tsv'),delimiter='\t'))
meta=list(csv.DictReader(open(O/'linked_locus_audit.tsv'),delimiter='\t'))
flanks=list(csv.DictReader(open(O/'flanking_genes_2kb.tsv'),delimiter='\t'))
def fasta(f):
 d={};n=None;buf=[]
 for line in open(f):
  if line.startswith('>'):
   if n:d[n]=''.join(buf).upper()
   n=line[1:].split()[0];buf=[]
  elif n:buf.append(line.strip())
 if n:d[n]=''.join(buf).upper()
 return d
def gff(gca):
 d={};lengths={}
 for line in open(B/'02_prokka/gff'/f'{gca}.gff'):
  if line.startswith('##sequence-region'):
   z=line.split();lengths[z[1]]=int(z[3]);continue
  if line.startswith('#'):continue
  z=line.rstrip('\n').split('\t')
  if len(z)<9 or z[2]!='CDS':continue
  attrs=dict(x.split('=',1) for x in z[8].split(';') if '=' in x)
  pid=attrs.get('ID')
  if pid:d[pid]={'contig':z[0],'start':int(z[3]),'end':int(z[4]),'strand':z[6],'product':attrs.get('product','')}
 return d,lengths
def tsv(path,rows):
 if not rows:return
 with open(path,'w',newline='') as h:
  writer=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
recover=[];unresolved=[]
for r in missing:
 g=r['GCA'];pair=rows[g];d,lens=gff(g)
 a=d[pair['protein1']];b=d[pair['protein2']]
 if a['contig']!=b['contig']:
  unresolved.append({'GCA':g,'reason':'different annotated contigs'});continue
 seq=fasta(S/f'{g}.fna' if (S/f'{g}.fna').exists() else O/'genomes'/f'{g}.fna')
 contig=a['contig']; n=lens[contig]
 options=[(name,s) for name,s in seq.items() if len(s)==n]
 if len(options)!=1:
  unresolved.append({'GCA':g,'reason':'contig length ambiguous or absent','n_matches':len(options)});continue
 originalname,nt=options[0]
 left=max(a['start'],b['start']) if False else min(a['end'],b['end'])+1
 right=max(a['start'],b['start'])-1
 gap=right-left+1
 if gap!=int(pair['distance_bp']):
  unresolved.append({'GCA':g,'reason':f'gap mismatch {gap} vs {pair["distance_bp"]}'});continue
 strandleft= a['strand'] if a['start']<b['start'] else b['strand']
 strandright=b['strand'] if a['start']<b['start'] else a['strand']
 dna=nt[left-1:right]
 if len(dna)!=gap:
  unresolved.append({'GCA':g,'reason':'sequence extraction length mismatch'});continue
 x={'GCA':g,'Genus':pair['Genus'],'HMT':pair['HMT'],'gap_bp':gap,
 'orientation_left_to_right':strandleft+strandright,'divergent':str(strandleft=='-' and strandright=='+'),
 'accession_contig':originalname,'start_1based':left,'end_1based':right,
 'sequence_access':'contig_length_mapped','GC_intergenic':round((dna.count('C')+dna.count('G'))/len(dna),4),
 'HxlR_protein':pair['protein1'],'SnoaL_protein':pair['protein2']}
 meta.append(x)
 recover.append((f'{g}|{pair["Genus"]}|{pair["HMT"]}|gap={gap}|orientation={strandleft+strandright}',dna))
 for pid,gene in d.items():
  if gene['contig']==contig and gene['end']>=min(a['start'],b['start'])-2000 and gene['start']<=max(a['end'],b['end'])+2000:
   tp='HxlR' if pid==pair['protein1'] else 'SnoaL2' if pid==pair['protein2'] else 'flanking'
   flanks.append({'GCA':g,'Genus':pair['Genus'],'HMT':pair['HMT'],'gene':pid,'gene_type':tp,'contig':originalname,
       'start':gene['start'],'end':gene['end'],'strand':gene['strand'],'product':gene['product']})
 print('RECOVERED',g,originalname,gap,flush=True)
with open(O/'HxlR_SnoaL2_intergenic_all_available.fna','a') as h:
 for name,seq in recover:h.write('>'+name+'\n'+seq+'\n')
tsv(O/'linked_locus_audit.tsv',meta);tsv(O/'flanking_genes_2kb.tsv',flanks);tsv(O/'unresolved_contig_mapping.tsv',unresolved)
summary={'total':67,'recovered':len(meta),'new_recovered':len(recover),'unresolved':len(unresolved),
         'tight_divergent_le500':sum(str(r['divergent']).lower()=='true' and int(r['gap_bp'])<=500 for r in meta),
         'by_genus':dict(collections.Counter(r['Genus'] for r in meta)),
         'gaps':{'median':statistics.median(int(r['gap_bp']) for r in meta)}}
(O/'audit_summary.json').write_text(json.dumps(summary,indent=2))
print('DONE',json.dumps(summary),flush=True)
