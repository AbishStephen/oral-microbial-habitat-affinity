#!/usr/bin/env python3
"""Download versioned HOMD PROKKA GFF for representative HMT genomes
carrying both PF01638 and PF12680; bounded network I/O only."""
from pathlib import Path
import csv,urllib.request,concurrent.futures,time
from collections import Counter
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
S=Path('__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014/GFF_both')
S.mkdir(exist_ok=True)
rows=list(csv.DictReader(open(O/'HOMD631_HxlR_SnoaL2_oneHMT_oneGenome_census.tsv'),delimiter='\t'))
rows=[r for r in rows if int(r['PF01638'])>0 and int(r['PF12680'])>0]
print('TO_FETCH',len(rows),flush=True)
def worker(r):
 g=r['GCA'];p=S/f'{g}.gff'
 if p.exists() and p.stat().st_size>10000:return dict(GCA=g,site=r['hmt_primary_body_site_w_abundance'],status='cached',bytes=p.stat().st_size)
 for attempt in range(3):
  try:
   url=f'https://www.homd.org/ftp/genomes/PROKKA/current/gff/{g}.gff'
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'HxlR-eHOMD-audit/1.0'}),timeout=70) as h:content=h.read()
   if len(content)<10000 or not (content.startswith(b'##gff') or b'|##gff-version' in content[:80]):raise ValueError('unexpected GFF bytes')
   p.write_bytes(content)
   return dict(GCA=g,site=r['hmt_primary_body_site_w_abundance'],status='OK',bytes=len(content))
  except Exception as e:
   if attempt==2:return dict(GCA=g,site=r['hmt_primary_body_site_w_abundance'],status='FAILED',bytes=0,error=str(e)[:110])
   time.sleep(1.5*(attempt+1))
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as e:
 for i,result in enumerate(e.map(worker,rows),1):
  results.append(result)
  if i%30==0:print('GFF',i,'/',len(rows),Counter(r['status'] for r in results),flush=True)
with open(O/'HOMD631_both_domain_GFF_download_audit.tsv','w') as h:
 writer=csv.DictWriter(h,fieldnames=['GCA','site','status','bytes','error'],delimiter='\t');writer.writeheader();writer.writerows(results)
print('COMPLETE',Counter(r['status'] for r in results),flush=True)
