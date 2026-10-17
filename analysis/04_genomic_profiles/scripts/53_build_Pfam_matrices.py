import csv,gzip
from pathlib import Path
from collections import defaultdict,Counter
B=Path("__PROJECT_ROOT__/study13")
S=Path("__WORK_ROOT__/HOMD_convergence_expanded13_work")
R=S/"pfam/results"; O=B/"06_hmt_matrices"; D=B/"04_domains"
O.mkdir(exist_ok=True); D.mkdir(exist_ok=True)
files=sorted(R.glob("chunk_*.domtbl")); assert len(files)==48
with open(B/"01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv") as f:
    rows=list(csv.DictReader(f,delimiter="\t"))
meta={r["GCA"]:r for r in rows}; genomes=sorted(meta)
byh=defaultdict(list)
for g,r in meta.items(): byh[r["HMT"]].append(g)
hmts=sorted(byh)

# Authoritative protein -> genome mapping from the master catalogue.
protein_to_gca={}
with open(B/"02_prokka/master_protein_catalogue.tsv") as f:
    rr=csv.DictReader(f,delimiter="\t")
    for z in rr:
        protein_to_gca[z["protein_id"]]=z["GCA"]

gc=defaultdict(Counter); fhg=defaultdict(lambda:defaultdict(set))
names={}; hits=Counter(); proteins=set(); total=unknown=0
with gzip.open(D/"Pfam_domain_hits.tsv.gz","wt") as out:
    w=csv.writer(out,delimiter="\t")
    w.writerow(["protein_id","GCA","HMT","Pfam","Pfam_name","i_evalue","domain_score","hmm_from","hmm_to","ali_from","ali_to","env_from","env_to"])
    for fn in files:
        with open(fn,errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"): continue
                p=line.rstrip().split(None,22)
                if len(p)<22: continue
                pfname,acc,protein=p[0],p[1],p[3]
                pf=acc.split(".")[0] if acc!="-" else pfname
                gca=protein_to_gca.get(protein)
                if gca not in meta: unknown+=1; continue
                h=meta[gca]["HMT"]; names[pf]=pfname
                gc[gca][pf]+=1; fhg[pf][h].add(gca); hits[pf]+=1
                proteins.add(protein); total+=1
                w.writerow([protein,gca,h,pf,pfname,p[12],p[13],p[15],p[16],p[17],p[18],p[19],p[20]])
assert unknown==0, unknown
features=sorted(hits)
transfer=[]
for pf in features:
    genera={meta[byh[h][0]]["Genus"] for h in fhg[pf]}
    if len(genera)>=2: transfer.append(pf)
with open(O/"X2_Pfam_dictionary.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["Pfam","Pfam_name","domain_hits","n_genomes","n_HMTs","n_genera","plaque_HMTs","mucosal_HMTs","transferable_ge2genera"])
    for pf in features:
        hs=set(fhg[pf]); gs={meta[byh[h][0]]["Genus"] for h in hs}
        gen=set().union(*(fhg[pf][h] for h in hs))
        pl=sum(meta[byh[h][0]]["Ecology"]=="plaque" for h in hs)
        mu=sum(meta[byh[h][0]]["Ecology"]=="mucosal" for h in hs)
        w.writerow([pf,names[pf],hits[pf],len(gen),len(hs),len(gs),pl,mu,"YES" if len(gs)>=2 else "NO"])
def write(path,feats,mode):
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["HMT","Genus","Species","Ecology","Median_Y","n_genomes"]+feats)
        for h in hmts:
            gs=byh[h]; m=meta[gs[0]]; n=len(gs)
            if mode=="prev": vals=[len(fhg[pf].get(h,set()))/n for pf in feats]
            else: vals=[sum(gc[g].get(pf,0) for g in gs)/n for pf in feats]
            w.writerow([h,m["Genus"],m["Species"],m["Ecology"],m["Median_Y"],n]+vals)
with open(O/"X2_Pfam_genome_copy_number.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["GCA","HMT","Genus","Species","Ecology","Median_Y"]+features)
    for g in genomes:
        m=meta[g]; w.writerow([g,m["HMT"],m["Genus"],m["Species"],m["Ecology"],m["Median_Y"]]+[gc[g].get(pf,0) for pf in features])
write(O/"X2_Pfam_HMT_prevalence_ALL.tsv",features,"prev")
write(O/"X2_Pfam_HMT_mean_copy_ALL.tsv",features,"copy")
write(O/"X2_Pfam_HMT_prevalence_TRANSFERABLE.tsv",transfer,"prev")
write(O/"X2_Pfam_HMT_mean_copy_TRANSFERABLE.tsv",transfer,"copy")
print("Pfam result files:",len(files))
print("Total significant domain hits:",total)
print("Proteins with >=1 Pfam hit:",len(proteins))
print("Pfam families:",len(features))
print("HMTs:",len(hmts))
print("Pfams in >=2 genera:",len(transfer))
print("Single-genus Pfams:",len(features)-len(transfer))
print("Unknown protein mappings:",unknown)
