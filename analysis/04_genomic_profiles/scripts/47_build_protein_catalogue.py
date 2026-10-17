import csv, os
from pathlib import Path

B=Path("__PROJECT_ROOT__/study13")
P=B/"02_prokka"

with open(B/"01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv") as f:
    M={r["GCA"]:r for r in csv.DictReader(f,delimiter="\t")}

out=open(P/"master_protein_catalogue.tsv","w",newline="")
qc=open(P/"protein_qc.tsv","w",newline="")

w=csv.writer(out,delimiter="\t")
w.writerow(["protein_id","GCA","HMT","Genus","Species","Ecology","Median_Y",
            "gene","EC_number","COG","product"])

q=csv.writer(qc,delimiter="\t")
q.writerow(["GCA","HMT","faa_proteins","tsv_CDS","match"])

tot=gene_n=cog_n=ec_n=0
mismatches=0

for gca,m in sorted(M.items()):
    tsv=P/"tsv"/f"{gca}.tsv"
    faa=P/"faa"/f"{gca}.faa"

    cds=0
    with open(tsv,errors="replace") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if r.get("ftype")!="CDS":
                continue
            cds+=1
            gene=(r.get("gene") or "").strip()
            ec=(r.get("EC_number") or "").strip()
            cog=(r.get("COG") or "").strip()
            product=(r.get("product") or "").strip()

            w.writerow([
                r["locus_tag"],gca,m["HMT"],m["Genus"],m["Species"],
                m["Ecology"],m["Median_Y"],gene,ec,cog,product
            ])

            tot+=1
            gene_n+=bool(gene)
            cog_n+=bool(cog)
            ec_n+=bool(ec)

    faa_n=0
    with open(faa,errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                faa_n+=1

    ok=(faa_n==cds)
    if not ok:mismatches+=1
    q.writerow([gca,m["HMT"],faa_n,cds,"YES" if ok else "NO"])

out.close(); qc.close()

print("Genomes:",len(M))
print("Total CDS/proteins:",tot)
print("With gene symbol:",gene_n)
print("With COG:",cog_n)
print("With EC:",ec_n)
print("FAA/TSV genome mismatches:",mismatches)
