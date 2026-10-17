#!/usr/bin/env python3
"""Download 631 one-per-HMT annotated proteomes from publicly released HOMD
2025 PROKKA archive; keep versioned, bounded to max six HTTP connections.
Download stage only, no heavy computation on the study workflow login node."""
from pathlib import Path
import csv,concurrent.futures,urllib.request,time,json,gzip
from collections import Counter
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
S=Path('__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014')
S.mkdir(exist_ok=True)
t=list(csv.DictReader(open(O/'HOMD_2025_FTP_taxon_representative_panel.tsv'),delimiter='\t'))
def fetch(row):
 g=row['GCA'];p=S/(g+'.faa.gz')
 if p.exists() and p.stat().st_size>2500:return dict(GCA=g,HMT=row['HMT_link'],status='cached',bytes=p.stat().st_size)
 url=f'https://www.homd.org/ftp/genomes/PROKKA/current/faa/{g}.faa'
 err=''
 for attempt in range(3):
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'eHOMD-HxlR-audit/1.0'}),timeout=65) as f:content=f.read()
   if len(content)<10000 or not content.startswith(b'>'):raise ValueError('invalid fasta')
   with gzip.open(p,'wb',compresslevel=3) as out:out.write(content)
   return dict(GCA=g,HMT=row['HMT_link'],status='downloaded',bytes=p.stat().st_size)
  except Exception as e:
   err=str(e)[:140];time.sleep((attempt+1)*1.5)
 return dict(GCA=g,HMT=row['HMT_link'],status='failed',bytes=0,error=err)
result=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
 for i,r in enumerate(ex.map(fetch,t),1):
  result.append(r)
  if i%50==0:
   print('PROGRESS',i,'/',len(t),Counter(v['status'] for v in result),flush=True)
with open(O/'HOMD_631_protein_fasta_download_audit.tsv','w') as f:
 keys=['GCA','HMT','status','bytes','error']
 w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(result)
print('FINISHED',Counter(v['status'] for v in result),flush=True)
