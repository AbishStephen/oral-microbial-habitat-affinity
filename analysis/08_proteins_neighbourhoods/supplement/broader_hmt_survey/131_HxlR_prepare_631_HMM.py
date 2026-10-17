#!/usr/bin/env python3
"""Build per-HMT protein challenge panel and isolate Pfam-A HMM for PF01638+PF12680."""
from pathlib import Path
import csv,gzip,json
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
S=Path('__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014')
m=list(csv.DictReader(open(O/'HOMD_2025_FTP_taxon_representative_panel.tsv'),delimiter='\t'))
faa=S/'HOMD_631_HMTrepresentatives_2025_PROKKA.faa'
ct=0
with open(faa,'w') as o:
 for i,x in enumerate(m,1):
  g=x['GCA'];p=S/f'{g}.faa.gz'
  if not p.exists():raise FileNotFoundError(p)
  with gzip.open(p,'rt') as inp:
   for line in inp:
    if line.startswith('>'):
     protein=line[1:].strip().split()[0]
     o.write('>'+g+'|'+protein+'\n');ct+=1
    else:o.write(line)
  if i%100==0:print('COMBINED',i,'HMT','proteins',ct,flush=True)
print('DONEFAA',faa.stat().st_size,'proteins',ct,flush=True)
source=Path('__WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm')
models={'PF01638':None,'PF12680':None}
buffer=[];acc=''
with source.open() as f:
 for line in f:
  if line.startswith('HMMER3/f'):buffer=[line];acc='';continue
  if not buffer:continue
  buffer.append(line)
  if line.startswith('ACC'):acc=line.split()[1].split('.')[0]
  if line.startswith('//'):
   if acc in models and models[acc] is None:
    dest=S/f'{acc}.hmm';dest.write_text(''.join(buffer));models[acc]=str(dest);print('MODEL',acc,'bytes',dest.stat().st_size,flush=True)
   buffer=[]
   if all(models.values()):break
assert all(models.values()),models
(O/'HOMD_631_representative_proteome_scan_manifest.json').write_text(json.dumps({'n_genomes':len(m),'n_proteins':ct,'ftp_prokka_release':'2025-07-28','families':list(models),'hmm_source':str(source)},indent=2))
