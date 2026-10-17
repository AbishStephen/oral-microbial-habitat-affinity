#!/usr/bin/env python3
"""Archive HxlR-family comparative-genomics TF-binding sites from RegPrecise.
Read-only HTML transfer, compact 11 regulog datasets, label as INFERRED unless
independently experimentally verified; do not combine different motif families."""
from pathlib import Path
import urllib.request,re,html,time,csv,collections,json
B=Path('__PROJECT_ROOT__/study13')
O=B/'16_HxlR_reference_phylogeny_20261014'
S=Path('__WORK_ROOT__/HOMD_HxlR_reference_20261014')
url='https://regprecise.lbl.gov/collection_tffam.jsp?tffamily_id=30'
def fetch(url):
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'OralHxlR-citation-audit/1.0'}),timeout=20) as f:return f.read().decode()
coll=fetch(url)
ids=sorted(set(re.findall(r'regul og_id=([0-9]+)'.replace(' ',''),coll)))
print('REGULOGS',ids,flush=True)
def clean(x):
 x=re.sub(r'<[^>]+>',' ',x)
 return html.unescape(' '.join(x.split()))
out=[];meta=[]
for k in ids:
 page=fetch('https://regprecise.lbl.gov/sites.jsp?regulog_id='+k)
 (S/f'regprecise_sites_{k}.html').write_text(page)
 title=clean(re.search(r'<h1[^>]*>(.*?)</h1>',page,re.S).group(1)) if re.search(r'<h1[^>]*>(.*?)</h1>',page,re.S) else ''
 # infer readable regulator name from header h2/h3 also
 text=clean(page)
 m=re.search(r'Profile of regulator\s+(.+?)\s+in\s+(.+?)\s+Properties',text)
 name=m.group(1).strip() if m else ('unresolved_'+k)
 lineage=m.group(2).strip() if m else 'unknown'
 motifs=re.findall(r'<td\s+class=["\']sequence["\']>([ACGT]{14,35})</td>',page,re.I)
 # extract sequence-bearing rows, regulator/site position target fields in table
 rows=re.findall(r'<tr\b[^>]*>(.*?)</tr>',page,re.S|re.I)
 added=0
 for row in rows:
  matched=re.search(r'<td\s+class=["\']sequence["\']>([ACGT]{14,35})</td>',row,re.I)
  if not matched:continue
  tds=[clean(t) for t in re.findall(r'<td\b[^>]*>(.*?)</td>',row,re.S|re.I)]
  seq=matched.group(1).upper()
  out.append({'regulog_id':k,'regulator':name,'lineage':lineage,'sequence':seq,
    'full_table_cells':' | '.join(tds),'source_url':'https://regprecise.lbl.gov/sites.jsp?regulog_id='+k})
  added+=1
 meta.append(dict(regulog_id=k,name=name,lineage=lineage,n_sites=added,n_unique_motifs=len(set(motifs))))
 print('REGULOG',k,name,lineage,'sites',added,flush=True)
df=__import__('pandas').DataFrame(out)
df.to_csv(O/'RegPrecise_HxlR11_regulog_annotated_binding_sites.tsv',sep='\t',index=False)
__import__('pandas').DataFrame(meta).to_csv(O/'RegPrecise_HxlR11_regulog_metadata.tsv',sep='\t',index=False)
(O/'RegPrecise_motifs_provenance.json').write_text(json.dumps({'source':url,'regulog_count':len(ids),'binding_sites_extracted':len(out),
 'caution':'Sites curated computationally in RegPrecise, not all direct biochemical binding evidence; do not pool motifs across regulogs or infer chemical sensor specificity from binding motif alone.'},indent=2))
