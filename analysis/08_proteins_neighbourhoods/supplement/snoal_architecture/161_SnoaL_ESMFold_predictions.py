#!/usr/bin/env python3
"""Fold five nonredundant public oral SnoaL-family representatives via ESM Atlas
public ESMFold API. Sequential requests to respect service rate limits.
No new model weights or genome datasets. Saves PDBs in lab project directory."""
from pathlib import Path
import urllib.request,urllib.error,time,json,hashlib
B=Path('__PROJECT_ROOT__/study13')
S=Path('__WORK_ROOT__/SnoaL2_subfamilies_20261016')
O=B/'18_SnoaL2_function_gene_integrity_20261016'
names=['GCA_001956355.1_01713','GCA_000191105.1_00970','GCA0028848951_00983','GCA_000176755.2_01530','GCA_900573985.1_01698']
records=[]
for name in names:
 seq=(S/f'{name}.sequence.txt').read_text().strip()
 target=O/f'{name}_ESMFold_v1.pdb'
 if target.exists() and target.stat().st_size>5000 and 'ATOM ' in target.read_text()[:1500]:
  records.append(dict(protein_id=name,length=len(seq),status='model_preexisting',file=str(target)))
  continue
 print('FOLD_START',name,len(seq),flush=True)
 req=urllib.request.Request('https://api.esmatlas.com/foldSequence/v1/pdb/',data=seq.encode(),headers={'Content-Type':'application/x-www-form-urlencoded','User-Agent':'OralProtein_StructureAnalysis/1.0'},method='POST')
 try:
  start=time.monotonic()
  with urllib.request.urlopen(req,timeout=105) as r:
   blob=r.read()
   status=r.status
  elapsed=round(time.monotonic()-start,1)
  if status==200 and blob.count(b'ATOM  ')>len(seq):
   target.write_bytes(blob)
   print('FOLD_OK',name,len(blob),'seconds',elapsed,flush=True)
   records.append(dict(protein_id=name,length=len(seq),status='predicted',http=status,seconds=elapsed,file=str(target),sha256=hashlib.sha256(blob).hexdigest()))
  else:
   print('FOLD_INVALID',name,status,blob[:130],flush=True)
   records.append(dict(protein_id=name,length=len(seq),status='invalid_response',http=status,response_head=blob[:150].decode(errors='replace')))
 except Exception as e:
  print('FOLD_ERROR',name,str(e)[:200],flush=True)
  records.append(dict(protein_id=name,length=len(seq),status='error',error=str(e)[:250]))
 (O/'ESMFold_v1_prediction_manifest.json').write_text(json.dumps(records,indent=2))
 time.sleep(4)
print('FINISHED',len(records),flush=True)
