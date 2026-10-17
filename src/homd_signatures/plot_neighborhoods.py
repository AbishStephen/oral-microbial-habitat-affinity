"""Render ten candidate-domain neighbourhood panels from measured GFF coordinates."""
from .figures import validate_export
from pathlib import Path
from collections import defaultdict
import argparse, json, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from .neighborhoods import read_tsv, PFAMS, NAMES, PREFERRED

COL={'focal':'#224961','partner':'#D7A747','annotated':'#AFC5CE','unknown':'#DEE3E5',
     'plaque':'#CD725C','mucosal':'#278F91','text':'#203D4C','muted':'#627983'}

def arrow(ax,start,end,y,strand,color):
    lo,hi=-3500,5500
    s,e=max(lo,start),min(hi,end)
    if e<=s:return
    h=.13;head=min(140,(e-s)*.35)
    if strand=='+':xy=[(s,y-h/2),(e-head,y-h/2),(e,y),(e-head,y+h/2),(s,y+h/2)]
    else:xy=[(e,y-h/2),(s+head,y-h/2),(s,y),(s+head,y+h/2),(e,y+h/2)]
    ax.add_patch(Polygon(xy,closed=True,facecolor=color,edgecolor='white',lw=.4,zorder=3))
    if start<lo:ax.plot([lo,lo],[y-h/2,y+h/2],color=COL['muted'],lw=1)
    if end>hi:ax.plot([hi,hi],[y-h/2,y+h/2],color=COL['muted'],lw=1)

def label_for(r,focal,pf):
    if focal:return NAMES[pf].split(' /')[0]
    gene=re.sub(r'_\d+$','',r['gene'])
    if gene:return gene
    t=r['product'].lower()
    for sub,label in [('ribonuclease','RNase'),('quinone','quinone red.'),('disulfide reductase','disulfide red.'),('nadh oxidase','NADH ox.'),('methyltransferase','methyltransf.'),('malate dehydrogenase','malate DH'),('multidrug','MFS transporter')]:
        if sub in t:return label
    return ''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    loci=read_tsv(a.input/'selected_loci.tsv');genes=read_tsv(a.input/'locus_genes.tsv')
    bypf=defaultdict(list);bylocus=defaultdict(list)
    for r in loci:bypf[r['Pfam']].append(r)
    for r in genes:bylocus[r['locus_id']].append(r)
    assert set(bypf)==set(PFAMS)
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','font.size':8,'axes.unicode_minus':False})
    fig,axes=plt.subplots(5,2,figsize=(10,13.3));fig.subplots_adjust(left=.055,right=.975,bottom=.105,top=.9,hspace=.43,wspace=.12)
    fig.text(.055,.969,'SUPPLEMENTARY GENOMIC CONTEXTS',fontsize=9,fontweight='bold',color=COL['mucosal'])
    fig.text(.055,.946,'Gene arrangements around ten candidate protein domains',fontsize=17,fontweight='bold',color=COL['text'])
    fig.text(.055,.926,'30 representative loci from 16 HMTs; real CDS coordinates and strand orientations',fontsize=9,color=COL['muted'])
    order=PFAMS
    for panel,(pf,ax) in enumerate(zip(order,axes.flat)):
        orientation='plaque-oriented' if pf in PFAMS[:7] else 'mucosa-oriented'
        ax.set_title(f'{chr(65+panel)}   {NAMES[pf]}  |  {pf}',loc='left',fontsize=10,fontweight='bold',color=COL['text'],pad=15)
        ax.text(0,1.035,orientation+' in the exploratory screen',transform=ax.transAxes,fontsize=6.5,color=COL['muted'])
        rows=sorted(bypf[pf],key=lambda r:int(r['display_order']))
        for j,s in enumerate(rows):
            y=2.7-1.15*j
            taxon=s['Genus']+' '+s['Species']
            if len(taxon)>38:taxon=s['Genus'][0]+'. '+s['Species']
            ax.text(-3500,y+.58,taxon+'  '+s['HMT'],fontsize=7.2,fontstyle='italic',color=COL['text'],va='bottom')
            ax.text(5500,y+.58,'P' if s['Ecology']=='plaque' else 'M',fontsize=7.5,fontweight='bold',color=COL[s['Ecology']],ha='right',va='bottom')
            ax.plot([-3500,5500],[y,y],lw=.5,color='#CCD7DA',zorder=0)
            immediate=[]
            for r in bylocus[s['locus_id']]:
                focal=r['focal']=='True';off=int(r['CDS_offset']);t=(r['product']+' '+r['gene']).lower()
                pref=(not focal and abs(off)==1 and any(k in t for k in PREFERRED[pf]))
                color=COL['focal'] if focal else COL['partner'] if pref else COL['unknown'] if 'hypothetical' in t else COL['annotated']
                st,en=int(r['display_start']),int(r['display_end'])
                arrow(ax,st,en,y,r['display_strand'],color)
                if focal or pref or abs(off)==1:
                    label=label_for(r,focal,pf)
                    if label and -3400<(st+en)/2<5400:immediate.append(((st+en)/2,label,focal,pref))
            # Focal labels above; neighbour symbols below avoid adjacent-gene collisions.
            for x,label,focal,pref in sorted(immediate):
                ax.text(x,y+.12 if focal else y-.12,label,ha='center',va='bottom' if focal else 'top',fontsize=6.2,color=COL['text'],fontweight='bold' if focal else 'normal')
        ax.set_xlim(-3500,5500);ax.set_ylim(-.03,3.49);ax.set_yticks([])
        ax.set_xticks([-2000,0,2000,4000],['-2','0','2','4']);ax.tick_params(axis='x',length=2,labelsize=6.5,colors=COL['muted'],pad=2)
        ax.set_xlabel('Position relative to focal CDS start (kb)',fontsize=6.5,color=COL['muted'],labelpad=3)
        ax.spines[['left','right','top']].set_visible(False);ax.spines['bottom'].set_color('#CAD6DA');ax.spines['bottom'].set_linewidth(.6)
    handles=[Patch(facecolor=COL[k],label=l) for k,l in [('focal','Domain-bearing CDS'),('partner','Prioritized annotated neighbour'),('annotated','Other annotated CDS'),('unknown','Uncharacterized CDS')]]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.055,.028),ncol=2,frameon=False,fontsize=8,columnspacing=2.3)
    fig.text(.57,.049,'P / M: plaque / mucosal HMT affinity',fontsize=7,color=COL['muted'])
    fig.text(.57,.035,'Focal CDS oriented rightward; end bars indicate clipping.',fontsize=7,color=COL['muted'])
    fig.text(.055,.013,'Illustrative arrangements from an exploratory shortlist; adjacency does not establish shared regulation, activity or ecological causation.',fontsize=7,color=COL['muted'])
    stem=a.output/'Supplementary_candidate_domain_neighborhoods'
    exports={}
    for ext in ['pdf','svg','png']:
        export=stem.with_suffix('.'+ext)
        fig.savefig(export,dpi=240 if ext=='png' else 160,facecolor='white')
        exports[export.name]=validate_export(export)
    plt.close(fig)
    caption='''Supplementary figure | Genomic contexts of ten exploratory Pfam-domain candidates. Three representative loci are shown per domain, selected using the archived immediate-neighbour annotations, genus breadth and one assembly per HMT. Two SnoaL_4/MFS examples illustrate the previously observed streptococcal arrangement; the third shows a different genus. Arrows use actual 1-based inclusive GFF coordinates converted to length-preserving interbase display coordinates. Reverse-strand focal genes and their neighbourhoods are reverse-complemented so focal transcription points rightward. The view spans -3.5 to +5.5 kb relative to the displayed focal CDS start; end bars mark clipped genes. Up to three CDS on each side were extracted; contig ends can reduce the available context. Dark blue marks the domain-bearing gene, gold marks an immediately adjacent gene meeting a predefined annotation-based partner priority, blue-grey marks other annotated genes, and pale grey marks uncharacterized genes. P/M denotes the assigned HMT-level habitat affinity, not an experimentally measured phenotype of the displayed isolate. Gene symbols and product names are predicted annotations. Representative selection is descriptive and the displayed arrangements are not a prevalence estimate. None of the original 14 candidates passed all-profile family-wise correction. The four non-Pfam annotations in the shortlist are not included as if each represented one domain-bearing locus. Full accessions, contigs, original and display coordinates, products and selection rules are supplied in the figure source tables.'''
    (a.output/'Supplementary_candidate_domain_neighborhoods_caption.md').write_text(caption+'\n')
    import zipfile
    with zipfile.ZipFile(a.output/'verified_supplementary_exports.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name,data in exports.items():z.writestr(name,data)
        z.writestr('Supplementary_candidate_domain_neighborhoods_caption.md',caption+'\n')
    print(json.dumps(dict(figure=str(stem),panels=10,loci=len(loci),gene_records=len(genes))))

if __name__=='__main__':main()
