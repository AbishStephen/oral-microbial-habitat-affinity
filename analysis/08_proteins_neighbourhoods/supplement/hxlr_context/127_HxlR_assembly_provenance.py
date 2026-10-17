#!/usr/bin/env python3
from pathlib import Path
import subprocess,json,pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'13_HxlR_broad_and_cassette_20261014'
a=pd.read_csv(O/'Porphyromonas_Prevotella_all243_cassette_nucleotide_screen.tsv',sep='\t')
gca=a[a.anchor45_nt_hits>0].GCA.tolist()
exe='datasets'
p=subprocess.run([exe,'summary','genome','accession',*gca,'--as-json-lines'],capture_output=True,text=True,timeout=90)
if p.returncode:raise RuntimeError(p.stderr[:300])
results=[]
for line in p.stdout.splitlines():
 if not line.strip() or not line.startswith('{'):continue
 r=json.loads(line)
 ai=r.get('assembly_info') or {}
 org=r.get('organism') or {}
 bi=ai.get('biosample') or {}
 results.append({'GCA':r.get('accession'),'organism':org.get('organism_name'),
  'strain':(org.get('infraspecific_names') or {}).get('strain'),'taxid':org.get('tax_id'),
  'BioSample':bi.get('accession'),'BioProject':ai.get('bioproject_accession'),
  'assembly_name':ai.get('assembly_name'),'assembly_method':ai.get('assembly_method'),
  'sample_owner':(bi.get('owner') or {}).get('name')})
rr=pd.DataFrame(results).merge(a[['GCA','full_1789nt_exact','anchor45_nt_hits']],on='GCA',how='left')
rr.to_csv(O/'cassette_NCBI_assembly_provenance.tsv',sep='\t',index=False)
print(rr.to_string(index=False),flush=True)
print('distinct strains by species')
print(rr.groupby('organism').strain.nunique().to_string())
