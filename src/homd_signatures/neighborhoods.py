"""Select and extract coordinate-supported examples for all ten focal Pfams.

Selection is descriptive: prioritize the pre-existing annotated partner hypotheses,
then maximize genus breadth and choose one representative per HMT. Genome counts
are not treated as ecological replicates. Inputs and all choices are exported.
"""
from pathlib import Path
from collections import defaultdict, Counter
import argparse, csv, gzip, json, urllib.parse, hashlib

PFAMS = ['PF01638','PF04072','PF13577','PF13676','PF22166',
         'PF13413','PF18075','PF01745','PF01337','PF14470']
NAMES = dict(zip(PFAMS,['HxlR','LCM','SnoaL_4','TIR_2','WHD_ArnR1',
                       'HTH_25','FtsX_ECD','IPT / MiaA','Barstar','bPH_3']))
PREFERRED = {
 'PF01638': ['quinone','disulfide reductase','nadh oxidase','nitroreductase'],
 'PF04072': ['malate dehydrogenase'],
 'PF13577': ['multidrug resistance','mdtd'],
 'PF13676': ['defense','restriction','nadase'],
 'PF22166': ['quinone','methyltransferase','epimerase'],
 'PF13413': ['rodz','cell shape'],
 'PF18075': ['ftse'],
 'PF01745': ['miab'],
 'PF01337': ['ribonuclease'],
 'PF14470': ['phage','sh oct','shoct'],
}

def read_tsv(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p,rows):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if not rows: raise ValueError(f'Empty export: {p}')
    cols=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)

def preference(pf,neighbors):
    return sum(any(k in (r['neighbor_product']+' '+r.get('neighbor_gene','')).lower()
                   for k in PREFERRED[pf]) for r in neighbors)

def choose_loci(coords,adjacent,manifest):
    meta={r['GCA']:r for r in manifest}
    neighbors=defaultdict(list)
    for r in adjacent:
        neighbors[(r['Pfam'],r['GCA'],r['protein_id'])].append(r)
    groups=defaultdict(list);seen=set()
    for r in coords:
        key=(r['Pfam'],r['GCA'],r['protein_id'])
        if r['Pfam'] not in PFAMS or key in seen:continue
        seen.add(key)
        z=dict(r);m=meta[z['GCA']]
        z.update(Species=m['Species'],Ecology=m['Ecology'])
        z['partner_score']=preference(z['Pfam'],neighbors[key])
        z['annotated_neighbors']=sum('hypothetical' not in v['neighbor_product'].lower() for v in neighbors[key])
        groups[r['Pfam']].append(z)
    selected=[]
    for pf in PFAMS:
        items=groups[pf]
        if not items:raise ValueError(f'No coordinate evidence for {pf}')
        support=Counter(r['HMT'] for r in items)
        ranked=sorted(items,key=lambda r:(-r['partner_score'],-support[r['HMT']],-r['annotated_neighbors'],r['GCA'],r['protein_id']))
        picked=[];used_hmt=set();used_genus=set()
        # Two known streptococcal SnoaL4/MFS examples plus a distinct genus.
        if pf=='PF13577':
            for r in ranked:
                if r['partner_score'] and r['HMT'] not in used_hmt:
                    picked.append(r);used_hmt.add(r['HMT']);used_genus.add(r['Genus'])
                    if len(picked)==2:break
        while len(picked)<3:
            pool=[r for r in ranked if r['HMT'] not in used_hmt]
            diverse=[r for r in pool if r['Genus'] not in used_genus]
            if not pool:break
            r=(diverse or pool)[0];picked.append(r);used_hmt.add(r['HMT']);used_genus.add(r['Genus'])
        for i,r in enumerate(picked):
            r=dict(r);r['locus_id']=pf+'__'+r['protein_id'];r['display_order']=i
            r['n_domain_genomes_in_HMT']=len({x['GCA'] for x in items if x['HMT']==r['HMT']})
            r['selection_rule']='annotated partner priority; genus breadth; one representative per HMT; deterministic accession tie break'
            selected.append(r)
    return selected

def parse_gff(p,gca):
    rows=[]
    with Path(p).open() as f:
        for line in f:
            if line.startswith('##FASTA'):break
            if line.startswith('#'):continue
            x=line.rstrip().split('\t')
            if len(x)<9 or x[2]!='CDS':continue
            a={k:urllib.parse.unquote(v) for k,v in (s.split('=',1) for s in x[8].split(';') if '=' in s)}
            rows.append(dict(GCA=gca,contig=x[0],start=int(x[3]),end=int(x[4]),strand=x[6],
                             protein_id=a.get('ID',''),gene=a.get('gene',''),product=a.get('product','')))
    return rows

def extract(selected,all_genes):
    by_contig=defaultdict(list)
    for r in all_genes:
        z=dict(r);z['start']=int(z['start']);z['end']=int(z['end']);by_contig[z['GCA'],z['contig']].append(z)
    genes=[];loci=[]
    for s in selected:
        track=sorted(by_contig[s['GCA'],s['contig']],key=lambda r:(r['start'],r['end'],r['protein_id']))
        ii=[i for i,r in enumerate(track) if r['protein_id']==s['protein_id']]
        if len(ii)!=1:raise ValueError(('Nonunique or missing focal gene',s,ii))
        i=ii[0];focal=track[i]
        assert all(str(focal[k])==str(s[k]) for k in ['start','end','strand']),('GFF/coordinate mismatch',s)
        z=dict(s);z['reverse_complement_display']=focal['strand']=='-'
        z['n_left_CDS_available']=min(3,i);z['n_right_CDS_available']=min(3,len(track)-i-1)
        loci.append(z)
        for j in range(max(0,i-3),min(len(track),i+4)):
            r=dict(track[j]);r.update(locus_id=s['locus_id'],Pfam=s['Pfam'],HMT=s['HMT'],Species=s['Species'],Ecology=s['Ecology'],CDS_offset=j-i,focal=j==i)
            # GFF uses 1-based inclusive coordinates; display uses interbase lengths.
            if focal['strand']=='+':
                r.update(display_start=r['start']-focal['start'],display_end=r['end']-focal['start']+1,display_strand=r['strand'])
            else:
                r.update(display_start=focal['end']-r['end'],display_end=focal['end']-r['start']+1,display_strand='+' if r['strand']=='-' else '-')
            assert r['display_end']-r['display_start']==r['end']-r['start']+1
            genes.append(r)
    return loci,genes

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--coordinates',type=Path,required=True)
    ap.add_argument('--adjacent',type=Path,required=True);ap.add_argument('--manifest',type=Path,required=True)
    group=ap.add_mutually_exclusive_group(required=True);group.add_argument('--gff-dir',type=Path);group.add_argument('--cds-table',type=Path)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    selected=choose_loci(read_tsv(a.coordinates),read_tsv(a.adjacent),read_tsv(a.manifest))
    all_genes=[]
    if a.gff_dir:
        for g in sorted({r['GCA'] for r in selected}):all_genes.extend(parse_gff(a.gff_dir/(g+'.gff'),g))
        with gzip.open(a.output/'selected_genome_CDS.tsv.gz','wt') as f:
            w=csv.DictWriter(f,fieldnames=list(all_genes[0]),delimiter='\t');w.writeheader();w.writerows(all_genes)
    else:all_genes=read_tsv(a.cds_table)
    loci,genes=extract(selected,all_genes)
    write_tsv(a.output/'selected_loci.tsv',loci);write_tsv(a.output/'locus_genes.tsv',genes)
    summary=dict(candidate_Pfams=len(PFAMS),loci=len(loci),HMTs=len({r['HMT'] for r in loci}),genomes=len({r['GCA'] for r in loci}),gene_records=len(genes),interpretation='Illustrative genomic arrangements; no causal or genome-wide significant feature claim')
    (a.output/'extraction_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))

if __name__=='__main__':main()
