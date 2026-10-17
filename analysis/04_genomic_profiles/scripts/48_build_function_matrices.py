import csv
from pathlib import Path
from collections import defaultdict

B=Path("__PROJECT_ROOT__/study13")
CAT=B/"02_prokka/master_protein_catalogue.tsv"
MAN=B/"01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv"
OUT=B/"06_hmt_matrices"
OUT.mkdir(exist_ok=True)

# genome metadata
with open(MAN) as f:
    rows=list(csv.DictReader(f,delimiter="\t"))

meta={r["GCA"]:r for r in rows}
genomes=sorted(meta)

# feature counts per genome
cog=defaultdict(lambda:defaultdict(int))
ec=defaultdict(lambda:defaultdict(int))

with open(CAT) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        g=r["GCA"]

        c=(r.get("COG") or "").strip()
        if c:
            for x in c.replace(";"," ").replace(","," ").split():
                if x.startswith("COG"):
                    cog[g][x]+=1

        e=(r.get("EC_number") or "").strip()
        if e:
            for x in e.replace(";"," ").replace(","," ").split():
                if x:
                    ec[g][x]+=1

def build(name,data):

    features=sorted({x for g in data for x in data[g]})

    # genome x feature copy-number matrix
    fn=OUT/f"{name}_genome_copy_number.tsv"
    with open(fn,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["GCA","HMT","Genus","Species","Ecology","Median_Y"]+features)

        for g in genomes:
            m=meta[g]
            w.writerow([
                g,m["HMT"],m["Genus"],m["Species"],
                m["Ecology"],m["Median_Y"]
            ]+[data[g].get(x,0) for x in features])

    # group genomes by HMT
    byh=defaultdict(list)
    for g in genomes:
        byh[meta[g]["HMT"]].append(g)

    # HMT prevalence
    fn=OUT/f"{name}_HMT_prevalence.tsv"
    with open(fn,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["HMT","Genus","Species","Ecology","Median_Y","n_genomes"]+features)

        for h in sorted(byh):
            gs=byh[h]
            m=meta[gs[0]]

            vals=[
                sum(data[g].get(x,0)>0 for g in gs)/len(gs)
                for x in features
            ]

            w.writerow([
                h,m["Genus"],m["Species"],m["Ecology"],
                m["Median_Y"],len(gs)
            ]+vals)

    # HMT mean copy number
    fn=OUT/f"{name}_HMT_mean_copy.tsv"
    with open(fn,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["HMT","Genus","Species","Ecology","Median_Y","n_genomes"]+features)

        for h in sorted(byh):
            gs=byh[h]
            m=meta[gs[0]]

            vals=[
                sum(data[g].get(x,0) for g in gs)/len(gs)
                for x in features
            ]

            w.writerow([
                h,m["Genus"],m["Species"],m["Ecology"],
                m["Median_Y"],len(gs)
            ]+vals)

    print(name)
    print("  genomes :",len(genomes))
    print("  HMTs    :",len(byh))
    print("  features:",len(features))

build("COG",cog)
build("EC",ec)
