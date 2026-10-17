import csv,re
from pathlib import Path
from collections import defaultdict,Counter
B=Path("__PROJECT_ROOT__/study13")
O=B/"06_hmt_matrices"; M=B/"05_function/mappings"
ECM=O/"EC_genome_copy_number.tsv"
with open(ECM) as f:
    r=csv.reader(f,delimiter="\t"); hdr=next(r); rows=list(r)
meta_cols=6; ecs=hdr[meta_cols:]; complete=[e for e in ecs if re.fullmatch(r"\d+\.\d+\.\d+\.\d+",e)]
idx={e:i+meta_cols for i,e in enumerate(ecs)}
genome_meta={row[0]:row[:meta_cols] for row in rows}
byh=defaultdict(list)
for row in rows: byh[row[1]].append(row)
# mappings
rhea=defaultdict(set)
with open(M/"rhea2ec.tsv") as f:
    rr=csv.DictReader(f,delimiter="\t")
    for x in rr:
        ec=x["ID"].strip()
        if ec in idx: rhea[ec].add("RHEA:"+x["MASTER_ID"].strip())
def links(path,prefix1,prefix2):
    d=defaultdict(set)
    with open(path) as f:
        for line in f:
            a,b=line.rstrip().split("\t")[:2]
            a=a.removeprefix(prefix1); b=b.removeprefix(prefix2)
            if a in idx: d[a].add(b)
    return d
kp=links(M/"kegg_ec_pathway.tsv","ec:","path:")
km=links(M/"kegg_ec_module.tsv","ec:","md:")
def names(path):
    d={}
    with open(path) as f:
        for line in f:
            a,b=line.rstrip().split("\t",1); d[a]=b
    return d
pnames=names(M/"kegg_pathway_names.tsv"); mnames=names(M/"kegg_module_names.tsv")
# Rhea genome presence
rxn_to_ec=defaultdict(set)
for ec,rs in rhea.items():
    for z in rs: rxn_to_ec[z].add(ec)
rxns=sorted(rxn_to_ec)
g_rxn={}
for row in rows:
    s=set()
    for ec in complete:
        if float(row[idx[ec]])>0:
            s.update(rhea.get(ec,()))
    g_rxn[row[0]]=s
def transfer_features(features, present_func):
    keep=[]
    for ft in features:
        genera=set()
        for h,hrs in byh.items():
            if any(present_func(r,ft) for r in hrs): genera.add(hrs[0][2])
        if len(genera)>=2: keep.append(ft)
    return keep
trx=transfer_features(rxns,lambda r,x:x in g_rxn[r[0]])
def write_prev(path,features,present_func):
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow(["HMT","Genus","Species","Ecology","Median_Y","n_genomes"]+features)
        for h in sorted(byh):
            rs=byh[h]; m=rs[0]
            vals=[sum(present_func(r,x) for r in rs)/len(rs) for x in features]
            w.writerow([h,m[2],m[3],m[4],m[5],len(rs)]+vals)
write_prev(O/"X3_Rhea_HMT_prevalence_ALL.tsv",rxns,lambda r,x:x in g_rxn[r[0]])
write_prev(O/"X3_Rhea_HMT_prevalence_TRANSFERABLE.tsv",trx,lambda r,x:x in g_rxn[r[0]])
with open(O/"X3_Rhea_dictionary.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["reaction","supporting_ECs","n_HMTs","n_genera","transferable_ge2genera"])
    for z in rxns:
        hs=[h for h,rs in byh.items() if any(z in g_rxn[r[0]] for r in rs)]
        gs={byh[h][0][2] for h in hs}
        w.writerow([z,";".join(sorted(rxn_to_ec[z])),len(hs),len(gs),"YES" if len(gs)>=2 else "NO"])
# pathway/module EC coverage
def build_cov(tag,mapping,name_map):
    feature_ec=defaultdict(set)
    for ec,fs in mapping.items():
        for ft in fs: feature_ec[ft].add(ec)
    # keep features with >=2 linked ECs represented in our EC universe
    features=sorted(ft for ft,s in feature_ec.items() if len(s)>=2)
    g_cov=defaultdict(dict)
    for row in rows:
        g=row[0]
        for ft in features:
            denom=feature_ec[ft]
            obs=sum(float(row[idx[e]])>0 for e in denom)
            g_cov[g][ft]=obs/len(denom)
    transferable=[]
    for ft in features:
        genera=set()
        for h,rs in byh.items():
            if any(g_cov[r[0]][ft]>0 for r in rs): genera.add(rs[0][2])
        if len(genera)>=2: transferable.append(ft)
    def write(path,fts):
        with open(path,"w",newline="") as f:
            w=csv.writer(f,delimiter="\t"); w.writerow(["HMT","Genus","Species","Ecology","Median_Y","n_genomes"]+fts)
            for h in sorted(byh):
                rs=byh[h]; m=rs[0]
                vals=[sum(g_cov[r[0]][ft] for r in rs)/len(rs) for ft in fts]
                w.writerow([h,m[2],m[3],m[4],m[5],len(rs)]+vals)
    write(O/f"X4_KEGG_{tag}_HMT_ECcoverage_ALL.tsv",features)
    write(O/f"X4_KEGG_{tag}_HMT_ECcoverage_TRANSFERABLE.tsv",transferable)
    with open(O/f"X4_KEGG_{tag}_dictionary.tsv","w",newline="") as f:
        w=csv.writer(f,delimiter="\t"); w.writerow([tag,"name","linked_ECs_in_dataset","transferable_ge2genera"])
        for ft in features: w.writerow([ft,name_map.get(ft,""),";".join(sorted(feature_ec[ft])),"YES" if ft in transferable else "NO"])
    print(tag,"features",len(features),"transferable",len(transferable))
    return len(features),len(transferable)
print("Complete EC features:",len(complete),"of",len(ecs))
print("Rhea reactions:",len(rxns),"transferable",len(trx))
build_cov("pathway",kp,pnames)
build_cov("module",km,mnames)
