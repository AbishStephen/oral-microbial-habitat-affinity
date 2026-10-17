#!/usr/bin/env python3
"""Hypothesis-guided functional programs from Prokka product text; audit candidate scope explicitly."""
from pathlib import Path
import csv,re,collections,numpy as np,pandas as pd
from scipy.stats import pearsonr
B=Path('__PROJECT_ROOT__/study13')
O=B/'10_additional_analyses_20261012';O.mkdir(exist_ok=True)
mm=list(csv.DictReader(open(B/'08_validation/confound_audit/genome_metadata_provenance.tsv'),delimiter='\t'))
patterns={
 'glucan_synthesis':r'glucosyltransferase|glucan synthase|dextransucrase',
 'adhesion_surface':r'adhesin|fimbrial|fimbria|sortase|collagen[- ]binding|fibronectin[- ]binding|surface[- ]associated',
 'oxygen_catalase':r'catalase(?![- ]peroxidase)|hydroperoxidase',
 'cytochrome_bd':r'cytochrome [bd]{1,2}(?:-type|[- ]ubiquinol|[- ]oxidase)|cytochrome bd|cytochrome d oxidase',
 'nitrate_nitrite':r'nitrate reductase|nitrite reductase|nitric oxide reductase',
 'mucin_glycosidase':r'sialidase|neuraminidase|fucosidase|hexosaminidase',
 'IgA_protease':r'IgA[ -]protease|immunoglobulin A[ -]protease'}
patterns={k:re.compile(v,re.I) for k,v in patterns.items()}
genome_to_groups=collections.defaultdict(set);num_missing=0;proof=collections.defaultdict(list)
for i,rec in enumerate(mm):
 gca=rec['GCA']
 f=B/'02_prokka/tsv'/f'{gca}.tsv'
 if not f.exists():num_missing+=1;continue
 with f.open() as fh:
  for r in csv.DictReader(fh,delimiter='\t'):
   if r['ftype']!='CDS':continue
   prod=(r.get('product') or '').strip()
   for k,pat in patterns.items():
    if pat.search(prod):
     genome_to_groups[gca].add(k)
     if len(proof[k])<30:proof[k].append((gca,rec['HMT'],r.get('gene',''),prod))
print('MISSING_PROKKA_TSV',num_missing,'GENOMES',len(mm),'PATTERN_HITS', {k:sum(k in v for v in genome_to_groups.values()) for k in patterns},flush=True)
hms=pd.read_csv(B/'08_validation/confound_audit/hmt_covariates_and_predictions.tsv',sep='\t')
data=[]
for h,group in pd.DataFrame(mm).groupby('HMT'):
 row={'HMT':h,'n_genomes':len(group)}
 for k in patterns:
  row[k]=np.mean([k in genome_to_groups.get(gca,set()) for gca in group.GCA])
 data.append(row)
d=hms.merge(pd.DataFrame(data),on='HMT',validate='one_to_one')
g=d.Genus.to_numpy(); y=d.Y.to_numpy(float)
groups=[np.flatnonzero(g==v) for v in np.unique(g)]
def cen(x):
 z=np.array(x,float).copy()
 for idx in groups:z[idx]-=z[idx].mean(axis=0)
 return z
c=cen(d[['GenomeSize_Mb','GC_percent','CDS']].to_numpy(float))
u,ss,v=np.linalg.svd(c,full_matrices=False);rk=sum(ss>ss.max()*1e-8);q=u[:,:rk]
proj=lambda x: x-q@(q.T@x)
yy=proj(cen(y))
out=[];mat=[]
for k in patterns:
 z=d[k].to_numpy(float)
 zz=proj(cen(z))
 r=float(pearsonr(yy,zz).statistic) if np.std(zz)>1e-12 else 0.
 out.append(dict(program=k,n_genomes_with_hit=sum(k in v for v in genome_to_groups.values()),n_HMT_nonzero=int(sum(z>0)),n_genus_nonzero=int(sum(d.groupby('Genus')[k].max()>0)),partial_r=r))
 mat.append(zz)
F=np.column_stack(mat);F/=np.linalg.norm(F,axis=0,keepdims=True)+1e-20
rng=np.random.default_rng(20261012);maxval=[]
raw=np.array([z['partial_r'] for z in out])
cnt=np.zeros(len(out),int)
for it in range(5000):
 yp=yy.copy()
 for ix in groups:yp[ix]=rng.permutation(yp[ix])
 yp=proj(yp);yp/=np.linalg.norm(yp)+1e-20
 co=np.abs(yp@F)
 cnt+=(co>=abs(raw)-1e-12)
 maxval.append(max(co))
mx=np.array(maxval)
for j,o in enumerate(out):
 o['p_pointwise']=(cnt[j]+1)/5001
 o['p_familywise_7programs']=(1+sum(mx>=abs(raw[j])-1e-12))/5001
pd.DataFrame(out).sort_values('p_pointwise').to_csv(O/'hypothesis_guided_genus_size_GC_CDS.tsv',sep='\t',index=False)
d[['HMT','Genus','Y']+list(patterns)].to_csv(O/'hypothesis_guided_HMT_program_prevalence.tsv',sep='\t',index=False)
with (O/'hypothesis_guided_example_product_matches.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['program','GCA','HMT','gene','annotated_product'])
 for k,v in proof.items():
  for x in v:w.writerow([k,*x])
print(pd.DataFrame(out).sort_values('p_pointwise').to_string(index=False))
