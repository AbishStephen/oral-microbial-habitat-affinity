#!/usr/bin/env python3
"""Rebuild the final 13-genus primary panel from frozen metadata and documented exclusions.
The cohort uses deterministic selection/capping and records the five unavailable
assemblies excluded from the final 1,121-genome analyses.
"""
import argparse, csv, json, math, statistics
from collections import defaultdict
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument('--metadata',type=Path,required=True)
ap.add_argument('--exclusions',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
cli=ap.parse_args(); M=cli.metadata; O=cli.output; O.mkdir(parents=True,exist_ok=True)
F=["abundance_Eren_V1V3.tsv","abundance_Eren_V3V5.tsv","abundance_Dewhirst.tsv",
   "abundance_HMP_MetaPhlan.tsv","abundance_HMP_RefSeq_V1V3.tsv","abundance_HMP_RefSeq_V3V5.tsv"]

def rd(p,h):
    z=open(p,encoding="utf-8-sig",errors="replace").readlines()
    i=next(i for i,x in enumerate(z) if x.startswith(h))
    return list(csv.DictReader(z[i:],delimiter="\t"))
def fl(x):
    try:return float(x)
    except:return 0.0
def hid(x):
    try:return "HMT-%03d"%int(float(x))
    except:return None

T={}
for r in rd(M/"HOMD_taxon_table.tsv","HMT-ID\t"):
    h=(r.get("HMT-ID") or "").strip()
    if h.startswith("HMT-"): T[h]=((r.get("Genus") or "").strip(),(r.get("Species") or "").strip())

R=[]
for fn in F:
    src=fn[10:-4]
    for r in rd(M/fn,"TAX\tHMT\t"):
        h=hid(r.get("HMT"))
        if not h or h=="HMT-000" or h not in T: continue
        v=tuple(fl(r.get(x+"_mean")) for x in ("SUBP","SUPP","BMU","TDO"))
        R.append((src,h,T[h][0],v))

D=defaultdict(list)
for src,h,g,v in R:
    if any(v): D[(src,g,v)].append(h)
bad={(s,h) for (s,g,v),hs in D.items() if len(set(hs))>1 for h in hs}
bad.add(("HMP_MetaPhlan","HMT-575"))

E=defaultdict(list)
for s,h,g,v in R:
    if (s,h) in bad: continue
    p=v[0]+v[1]; m=v[2]+v[3]; tot=p+m
    E[h].append((math.log2((p+1e-6)/(m+1e-6)),tot))

C={}
for h,a in E.items():
    ys=[x[0] for x in a]; med=statistics.median(ys)
    d=1 if med>0 else -1 if med<0 else 0
    agree=sum(y*d>0 for y,t in a)/len(a) if d else 0
    adequate=sum(t>=0.1 for y,t in a)
    tier=None
    if len(a)>=5 and adequate>=4 and agree>=0.8 and abs(med)>=1.5:
        tier="High"
    elif len(a)>=4 and adequate>=2 and agree>=0.8 and abs(med)>=1.0:
        tier="Moderate"
    if tier:
        C[h]=(T[h][0],T[h][1],"plaque" if med>0 else "mucosal",med,tier,len(a),adequate,agree)

G=rd(M/"HOMD_genomic_metadata.tsv","GENOME_ID\t")
Q=[]
for r in G:
    h=(r.get("HMT_ID") or "").strip(); g=(r.get("GENOME_ID") or "").strip()
    try:c=float(r["checkM2_completeness"]); x=float(r["checkM2_contamination"])
    except: continue
    if h in C and g and c>=95 and x<=5: Q.append((h,g,c,x,r))

bh=defaultdict(list)
for x in Q: bh[x[0]].append(x)
H={h for h,v in bh.items() if len({x[1] for x in v})>=5}
gd=defaultdict(set)
for h in H: gd[C[h][0]].add(C[h][2])
pg={g for g,d in gd.items() if {"plaque","mucosal"}<=d}
H=sorted(h for h in H if C[h][0] in pg)

# Apply the specified genome cap and accession deduplication in metadata order.
cap=[]; n=defaultdict(int)
for x in Q:
    if x[0] in H and n[x[0]]<50:
        cap.append(x); n[x[0]]+=1
seen=set(); uq=[]
for x in cap:
    if x[1] not in seen:
        seen.add(x[1]); uq.append(x)
uh=defaultdict(list)
for x in uq: uh[x[0]].append(x)
FH=sorted(h for h in H if len(uh[h])>=5)
FIN=[x for x in uq if x[0] in FH]

candidate_count=len(FIN)
with cli.exclusions.open() as f:
    excluded={r['GCA']:r['reason'] for r in csv.DictReader(f,delimiter='\t')}
matched_exclusions={x[1] for x in FIN}&set(excluded)
if matched_exclusions != set(excluded):
    raise ValueError('Exclusion accessions do not all occur in the candidate panel')
FIN=[x for x in FIN if x[1] not in excluded]
uh=defaultdict(list)
for x in FIN: uh[x[0]].append(x)
FH=sorted(uh)
assert len(pg)==13 and len(FH)==47 and len(FIN)==1121
assert all(len(uh[h])>=5 for h in FH)
with (O/'PRIMARY_HMT_MANIFEST.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t')
    w.writerow(['HMT','Genus','Species','Ecology','Median_Y','Ecology_tier','n_sources','n_adequate_sources','direction_agreement','n_genomes'])
    for h in FH:
        g,sp,e,y,t,ns,na,ag=C[h]
        w.writerow([h,g,sp,e,y,t,ns,na,ag,len(uh[h])])
with (O/'PRIMARY_GENOME_MANIFEST.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t')
    w.writerow(['GCA','HMT','Genus','Species','Ecology','Median_Y','Ecology_tier','CheckM2_complete','CheckM2_contam'])
    for h,g,c,x,r in FIN:
        genus,sp,e,y,t,*_=C[h]
        w.writerow([g,h,genus,sp,e,y,t,c,x])
summary=dict(candidate_genomes=candidate_count,excluded_accessions=sorted(excluded),primary_genomes=len(FIN),primary_HMTs=len(FH),primary_genera=len(pg),selection_history='Five documented unavailable genomes removed after original cap; no replacement sampling')
(O/'manifest_rebuild_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
