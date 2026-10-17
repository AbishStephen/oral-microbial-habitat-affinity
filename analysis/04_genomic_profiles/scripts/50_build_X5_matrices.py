import csv
from pathlib import Path
from collections import defaultdict, Counter

B=Path("__PROJECT_ROOT__/study13")
CL=B/"03_orthogroups/mmseqs_x5/HOMD13_x5_cluster.tsv"
MAN=B/"01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv"
OUT=B/"06_hmt_matrices"
OUT.mkdir(exist_ok=True)

# ----------------------------------------------------------
# Genome/HMT metadata
# ----------------------------------------------------------
with open(MAN) as f:
    rows=list(csv.DictReader(f,delimiter="\t"))

meta={r["GCA"]:r for r in rows}

by_hmt=defaultdict(list)
for g,r in meta.items():
    by_hmt[r["HMT"]].append(g)

hmts=sorted(by_hmt)

# ----------------------------------------------------------
# Cluster membership
# rep protein = cluster identity
# ----------------------------------------------------------
cluster_counts=defaultdict(Counter)      # cluster -> HMT -> protein copies
cluster_genomes=defaultdict(set)         # cluster -> genomes containing it
cluster_hmt_genomes=defaultdict(lambda:defaultdict(set))

# Do not infer genome accession from Prokka locus tags: new annotations may
# use sanitized locustags. Use the authoritative master catalogue mapping.
protein_to_gca={}
with open(B/"02_prokka/master_protein_catalogue.tsv") as f:
    rr=csv.DictReader(f,delimiter="\t")
    for z in rr:
        protein_to_gca[z["protein_id"]]=z["GCA"]

n_members=0
unknown=0

with open(CL) as f:
    for line in f:
        rep,protein=line.rstrip("\n").split("\t")[:2]
        gca=protein_to_gca.get(protein)

        if gca not in meta:
            unknown+=1
            continue

        h=meta[gca]["HMT"]

        cluster_counts[rep][h]+=1
        cluster_genomes[rep].add(gca)
        cluster_hmt_genomes[rep][h].add(gca)
        n_members+=1

assert unknown==0, f"{unknown} proteins could not be mapped to genomes"

reps=sorted(cluster_counts)

# Stable short IDs
cid={rep:f"X5_{i:06d}" for i,rep in enumerate(reps,1)}

# ----------------------------------------------------------
# Cluster feature dictionary
# ----------------------------------------------------------
transferable=[]

with open(OUT/"X5_cluster_dictionary.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow([
        "cluster_id","representative_protein","n_proteins",
        "n_genomes","n_HMTs","n_genera",
        "plaque_HMTs","mucosal_HMTs","transferable_ge2genera"
    ])

    for rep in reps:
        hs=set(cluster_counts[rep])
        gs={meta[g]["Genus"] for g in cluster_genomes[rep]}

        plaque=sum(
            meta[by_hmt[h][0]]["Ecology"]=="plaque"
            for h in hs
        )
        mucosal=sum(
            meta[by_hmt[h][0]]["Ecology"]=="mucosal"
            for h in hs
        )

        ok=len(gs)>=2
        if ok:
            transferable.append(rep)

        w.writerow([
            cid[rep],rep,
            sum(cluster_counts[rep].values()),
            len(cluster_genomes[rep]),
            len(hs),len(gs),
            plaque,mucosal,
            "YES" if ok else "NO"
        ])

# ----------------------------------------------------------
# HMT matrices
# ----------------------------------------------------------
def write_matrix(path, features, mode):

    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t")

        w.writerow([
            "HMT","Genus","Species","Ecology",
            "Median_Y","n_genomes"
        ]+[cid[x] for x in features])

        for h in hmts:
            gs=by_hmt[h]
            m=meta[gs[0]]
            n=len(gs)

            vals=[]

            for rep in features:
                if mode=="prevalence":
                    v=len(cluster_hmt_genomes[rep].get(h,set()))/n
                else:
                    v=cluster_counts[rep].get(h,0)/n
                vals.append(v)

            w.writerow([
                h,m["Genus"],m["Species"],m["Ecology"],
                m["Median_Y"],n
            ]+vals)

write_matrix(
    OUT/"X5_HMT_prevalence_ALL.tsv",
    reps,"prevalence"
)

write_matrix(
    OUT/"X5_HMT_mean_copy_ALL.tsv",
    reps,"copy"
)

write_matrix(
    OUT/"X5_HMT_prevalence_TRANSFERABLE.tsv",
    transferable,"prevalence"
)

write_matrix(
    OUT/"X5_HMT_mean_copy_TRANSFERABLE.tsv",
    transferable,"copy"
)

print("Proteins mapped:",n_members)
print("Clusters:",len(reps))
print("HMTs:",len(hmts))
print("Clusters present in >=2 genera:",len(transferable))
print("Single-genus clusters:",len(reps)-len(transferable))
