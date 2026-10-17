#!/usr/bin/env python3
"""NCBI genome FASTA retrieval for remaining P/Pre assemblies only.
Uses bounded, throttled concurrent requests, with per-accession provenance logs."""
from pathlib import Path
import pandas as pd, urllib.request,zipfile,io,time,concurrent.futures,json,threading
B=Path('__PROJECT_ROOT__/study13')
O=B/'13_HxlR_broad_and_cassette_20261014'; W=O/'genome_fasta';W.mkdir(exist_ok=True)
tab=pd.read_csv(O/'Porphyromonas_Prevotella_all243_cassette_nucleotide_screen.tsv',sep='\t')
targets=tab.loc[~tab.has_genome_fna,'GCA'].tolist()
def fetch(g):
 out=W/f'{g}.fna'
 if out.exists() and out.stat().st_size>100000:return dict(GCA=g,status='already_done',bytes=out.stat().st_size)
 url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{g}/download?include_annotation_type=GENOME_FASTA'
 for attempt in range(1,4):
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'HOMD-genomics-/1.0'}),timeout=55) as resp:content=resp.read()
   with zipfile.ZipFile(io.BytesIO(content)) as z:
    candidates=[k for k in z.namelist() if k.endswith(('.fna','.fasta'))]
    if not candidates:raise RuntimeError('no genomic FASTA in datasets zip')
    body=z.read(candidates[0])
   if not body.startswith(b'>') or len(body)<100000:raise RuntimeError('bad or short FASTA')
   out.write_bytes(body)
   return dict(GCA=g,status='OK',bytes=len(body))
  except Exception as e:
   if attempt==3:return dict(GCA=g,status='FAILED',reason=str(e)[:150])
   time.sleep(1.5*attempt)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 rec=[]
 for i,r in enumerate(pool.map(fetch,targets),1):
  rec.append(r)
  if i%10==0:print('NCBI',i,'/',len(targets),'OK',sum(x['status']!='FAILED' for x in rec),flush=True)
pd.DataFrame(rec).to_csv(O/'NCBI_nucleotide_retrieval_audit.tsv',sep='\t',index=False)
print('DONE',len(rec),'OK',sum(z['status']!='FAILED' for z in rec),'FAIL',sum(z['status']=='FAILED' for z in rec),flush=True)
