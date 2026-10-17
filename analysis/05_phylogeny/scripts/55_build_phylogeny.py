import csv
from pathlib import Path
from collections import defaultdict,Counter
B=Path("__PROJECT_ROOT__/study13")
MAN=B/"01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv"; P=B/"02_prokka"; OUT=B/"07_phylogeny"
OUT.mkdir(exist_ok=True)
with open(MAN) as f: rows=list(csv.DictReader(f,delimiter="\t"))
byh=defaultdict(list)
for r in rows: byh[r["HMT"]].append(r)
# choose best CheckM2 representative: highest completeness, then lowest contamination
rep={}
for h,rs in byh.items():
    rep[h]=sorted(rs,key=lambda r:(-float(r["CheckM2_complete"]),float(r["CheckM2_contam"]),r["GCA"]))[0]
cand=set(("rplA rplB rplC rplD rplE rplF rplI rplJ rplK rplL rplM rplN rplO rplP rplQ rplR rplS rplT rplU rplV rplW rplX rplY "
          "rpmA rpmB rpmC rpmD rpmE rpmF rpmG rpmH rpmI rpmJ "
          "rpsA rpsB rpsC rpsD rpsE rpsF rpsG rpsH rpsI rpsJ rpsK rpsL rpsM rpsN rpsO rpsP rpsQ rpsR rpsS rpsT rpsU "
          "tuf fusA infA infB infC tsf").split())
calls={}; counts=defaultdict(Counter)
for h,r in rep.items():
    g=r["GCA"]; d=defaultdict(list)
    with open(P/"tsv"/f"{g}.tsv") as f:
        for x in csv.DictReader(f,delimiter="\t"):
            if x.get("ftype")=="CDS" and (x.get("gene") or "") in cand:
                d[x["gene"]].append(x["locus_tag"])
    calls[h]=d
    for gene,loci in d.items(): counts[gene][h]=len(loci)
genes=sorted(g for g in cand if all(counts[g][h]==1 for h in rep))
assert len(genes)>=10, f"Too few universal markers: {genes}"
# load fasta sequences per representative genome
seqs={}
for h,r in rep.items():
    g=r["GCA"]; dd={}; key=None; buf=[]
    with open(P/"faa"/f"{g}.faa") as f:
        for line in f:
            if line.startswith(">"):
                if key: dd[key]="".join(buf)
                key=line[1:].split()[0]; buf=[]
            else: buf.append(line.strip())
        if key: dd[key]="".join(buf)
    seqs[h]=dd
for gene in genes:
    with open(OUT/f"{gene}.faa","w") as f:
        for h in sorted(rep):
            r=rep[h]; locus=calls[h][gene][0]
            label=f"{h}|{r['Genus']}_{r['Species']}|{r['GCA']}".replace(" ","_")
            f.write(f">{label}\n{seqs[h][locus]}\n")
with open(OUT/"representatives.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["HMT","GCA","Genus","Species","Ecology","Median_Y","CheckM2_complete","CheckM2_contam"])
    for h in sorted(rep):
        r=rep[h]; w.writerow([h,r["GCA"],r["Genus"],r["Species"],r["Ecology"],r["Median_Y"],r["CheckM2_complete"],r["CheckM2_contam"]])
with open(OUT/"marker_genes.txt","w") as f: f.write("\n".join(genes)+"\n")
print("HMT representatives:",len(rep))
print("Universal single-copy markers:",len(genes))
print("Markers:",",".join(genes))
