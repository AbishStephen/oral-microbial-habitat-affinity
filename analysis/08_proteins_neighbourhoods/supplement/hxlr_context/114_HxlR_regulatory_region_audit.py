#!/usr/bin/env python3
"""HxlR/SnoaL2 regulatory-region audit, genome-sequence provenance, no unvalidated promoter claims."""
from pathlib import Path
import csv,collections,json,gzip,re,urllib.request,io,zipfile,time,statistics
B=Path('__PROJECT_ROOT__/study13')
O=B/'12_HxlR_SnoaL_regulatory_genomics_20261013';O.mkdir(exist_ok=True)
S=Path('__WORK_ROOT__/HOMD_convergence_expanded13_work/raw_genomes')
E=B/'10_additional_analyses_20261012'
rows=list(csv.DictReader(open(E/'crossgenus_candidate_distinct_genes_within3kb_examples.tsv'),delimiter='\t'))
rows=[r for r in rows if r['Pfam1']=='PF01638' and r['Pfam2']=='PF12680']
# Assertions: one gene pair per genome expected, flag duplicate explicitly.
cnt=collections.Counter(r['GCA'] for r in rows)
print('ROWS',len(rows),'DISTINCT_GENOMES',len(cnt),'DUPLICATE_GENOMES',[(x,y) for x,y in cnt.items() if y>1],flush=True)
def gff(gca):
 p=B/'02_prokka/gff'/f'{gca}.gff'
 d={}
 for line in open(p):
  if line.startswith('#'):continue
  z=line.rstrip('\n').split('\t')
  if len(z)<9 or z[2]!='CDS':continue
  attr=dict(t.split('=',1) for t in z[8].split(';') if '=' in t)
  pr=attr.get('ID','')
  if pr:d[pr]={'contig':z[0].split('|',1)[-1],'start':int(z[3]),'end':int(z[4]),'strand':z[6],'product':attr.get('product','').replace('%20',' ')}
 return d
def fasta(f):
 seq={};name=None;buf=[]
 with open(f) as h:
  for line in h:
   if line.startswith('>'):
    if name:seq[name]=''.join(buf).upper()
    name=line[1:].strip().split()[0];buf=[]
   elif name:buf.append(line.strip())
  if name:seq[name]=''.join(buf).upper()
 return seq
def revcomp(s):
 return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fetch(gca):
 out=O/'genomes'/f'{gca}.fna';out.parent.mkdir(exist_ok=True)
 if out.exists() and out.stat().st_size>100000:return out,'existing_target'
 if (S/f'{gca}.fna').exists():return S/f'{gca}.fna','scratch_archive'
 url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{gca}/download?include_annotation_type=GENOME_FASTA'
 try:
  with urllib.request.urlopen(url,timeout=32) as resp:body=resp.read()
  z=zipfile.ZipFile(io.BytesIO(body));found=[x for x in z.namelist() if x.endswith(('.fna','.fasta'))]
  if not found:raise ValueError('no FASTA in NCBI package')
  seq=z.read(found[0])
  if len(seq)<100000 or not seq.startswith(b'>'):raise ValueError('invalid FASTA')
  out.write_bytes(seq);return out,'downloaded_NCBI'
 except Exception as e:return None,str(e)[:140]
meta=[];fasta_out=[];flanks=[];missing=[];gene_proteins={}
for i,r in enumerate(rows):
 g=r['GCA'];G=gff(g);a=G.get(r['protein1']);b=G.get(r['protein2'])
 if not a or not b or a['contig']!=b['contig']:
  missing.append({'GCA':g,'reason':'GFF_gene_or_contig_mismatch'});continue
 if a['end']<b['start']:
  left,right=a,b;inter_left=a['end']+1;inter_right=b['start']-1
 elif b['end']<a['start']:
  left,right=b,a;inter_left=b['end']+1;inter_right=a['start']-1
 else:
  missing.append({'GCA':g,'reason':'overlapping_CDS'});continue
 gap=inter_right-inter_left+1
 orientation=left['strand']+' '+right['strand']
 f,origin=fetch(g)
 if not f:
  missing.append({'GCA':g,'reason':origin});print('MISSING',g,origin,flush=True);continue
 try:
  seqs=fasta(f);nt=seqs.get(a['contig'])
  if not nt:raise ValueError('contig not in fasta: '+a['contig'])
  dna=nt[inter_left-1:inter_right]
  if len(dna)!=gap:raise ValueError('gap-length mismatch')
 except Exception as e:
  missing.append({'GCA':g,'reason':str(e)});continue
 # The intergenic sequence runs from lower to higher coordinate; shared promoters are distinct hypotheses
 conserved=(orientation=='- +')
 name=f'{g}|{r["Genus"]}|{r["HMT"]}|gap={gap}|orientation={orientation.replace(" ","")}'
 fasta_out.append((name,dna))
 row={'GCA':g,'Genus':r['Genus'],'HMT':r['HMT'],'gap_bp':gap,'orientation_left_to_right':orientation.replace(' ',''),'divergent':conserved,'accession_contig':a['contig'],'start_1based':inter_left,'end_1based':inter_right,'sequence_access':origin,'GC_intergenic':round((dna.count('C')+dna.count('G'))/len(dna),4) if dna else None,'HxlR_protein':r['protein1'],'SnoaL_protein':r['protein2']}
 meta.append(row)
 leftstart=min(a['start'],b['start']);rightend=max(a['end'],b['end'])
 nearby=sorted([(pid,d) for pid,d in G.items() if d['contig']==a['contig'] and d['end']>=leftstart-2000 and d['start']<=rightend+2000],key=lambda x:x[1]['start'])
 for pid,d in nearby:
  typ='HxlR' if pid==r['protein1'] else 'SnoaL2' if pid==r['protein2'] else 'flanking'
  flanks.append({**{'GCA':g,'Genus':r['Genus'],'HMT':r['HMT'],'gene':pid,'gene_type':typ},**d})
 if (i+1)%10==0:print('PROCESSED',i+1,'/',len(rows),flush=True)
with open(O/'HxlR_SnoaL2_intergenic_all_available.fna','w') as h:
 for n,s in fasta_out:h.write('>'+n+'\n'+ '\n'.join(s[k:k+80] for k in range(0,len(s),80))+'\n')
def table(path,data):
 if not data:return
 with open(path,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
table(O/'linked_locus_audit.tsv',meta)
table(O/'flanking_genes_2kb.tsv',flanks)
table(O/'missing_genomes.tsv',missing)
s=collections.Counter((r['Genus'],r['divergent']) for r in meta)
summ={'expected_linked_rows':len(rows),'sequence_recovered':len(meta),'genus_orientation_count':{str(k):v for k,v in s.items()},'missed':len(missing),'tight_divergent_le500':sum(r['divergent'] and r['gap_bp']<=500 for r in meta),'gap_distribution':dict(min=min((r['gap_bp'] for r in meta),default=0),median=statistics.median((r['gap_bp'] for r in meta)) if meta else 0,max=max((r['gap_bp'] for r in meta),default=0))}
(O/'audit_summary.json').write_text(json.dumps(summ,indent=2))
print('FINISHED',json.dumps(summ),flush=True)
