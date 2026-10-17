#!/usr/bin/env python3
"""Figure 4 | Coordinate-accurate synteny and alternative locus architectures.
Rebuild using bundled GFF-derived examples, the 67-linkage table and Figure 3
sequence-cluster audit. Representative tracks have true CDS lengths, order,
relative genomic distances and transcriptional strands; HxlR + / SnoaL_2 -
examples are reverse-complemented for comparable display.
"""
from pathlib import Path
import json, math, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D

P=Path(__file__).parent
rows=json.loads((P/'Figure4_real_GFF_loci.json').read_text())
allpairs=json.loads((P/'Figure4_all67_linkages.json').read_text())
types=json.loads((P/'Figure4_8_cluster_pair_types.json').read_text());assert len(types)==8 and sum(t['n_genomes'] for t in types)==67
assert len(allpairs)==67
assert len(rows)==8

NAVY='#16364c'; MUTED='#596f7c'; LINE='#d4e0e4'; RED='#d46b54'; GOLD='#dda53e';
TEAL='#178f95'; PURPLE='#8e76b6'; BLUELIGHT='#c1cbd0'; BLUEDARK='#739cab'; BG='#f7f9fa'
phy_col={'Actinomyces':'#d8a23c','Porphyromonas':'#be7774','Prevotella':'#be7774',
         'Streptococcus':'#9280bd','Neisseria':'#549baa','Fusobacterium':'#ac628b'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(17.1,11.3),dpi=210,facecolor='white')
fig.text(.041,.968,'EVOLUTIONARY GENOMIC ARCHITECTURES',color=TEAL,weight='bold',fontsize=9.7)
fig.text(.041,.940,'Recurrent synteny, different protein families',color=NAVY,weight='bold',fontsize=18.7)
fig.text(.041,.918,'HxlR–SnoaL_2 neighbourhoods recur across genera despite distinct sequence-cluster combinations',color=MUTED,fontsize=10.0)
fig.add_artist(Line2D([.041,.969],[.904,.904],transform=fig.transFigure,lw=0.85,color=LINE))

def panel_title(ax,l,t,sub):
    ax.text(0,1.050,l,transform=ax.transAxes,fontsize=18,weight='bold',color=NAVY,ha='left',va='bottom')
    ax.text(.055,1.061,t,transform=ax.transAxes,fontsize=11.2,weight='bold',color=NAVY,ha='left',va='bottom')
    ax.text(.055,1.019,sub,transform=ax.transAxes,fontsize=7.8,color=MUTED,ha='left',va='bottom')

def arrow_patch(start,end,y,strand,col,h=.17,opacity=1):
    if end-start<=1:return
    # polygon with real extent in genomic bp, direction indicated by arrowhead
    gap=end-start;head=min(130,gap*.34);lower=y-h/2;upper=y+h/2
    if strand=='+':
        pts=[(start,lower),(end-head,lower),(end,y),(end-head,upper),(start,upper)]
    else:
        pts=[(end,lower),(start+head,lower),(start,y),(start+head,upper),(end,upper)]
    return Polygon(pts,closed=True,facecolor=col,edgecolor='white',linewidth=.7,alpha=opacity,zorder=3)

def product_short(s):
    s=s.replace('%2C',',')
    if 'transposase' in s.lower():return 'transposase'
    if 'recombinase' in s.lower():return 'XerC'
    if 'repressor' in s.lower():return 'repressor'
    if 'nitroreductase' in s.lower():return 'nitroreductase'
    if 'efflux' in s.lower():return 'efflux'
    return ''

def normalize(row):
    # gene centers + orientation to reflect the coordinate-reversal choice
    center=(row['anchor_start']+row['anchor_end'])/2.
    orientation=-1 if row['type']=='HxlR/SnoaL2' and row['strand1']=='+' and row['strand2']=='-' else +1
    converted=[]
    for g in row['genes']:
        a,b=g['start']-center,g['end']-center
        if orientation<0:
            a,b=-b,-a
        strand=g['strand'] if orientation>0 else ('+' if g['strand']=='-' else '-')
        role='gene1' if g['gene']==row['protein1'] else ('gene2' if g['gene']==row['protein2'] else 'other')
        converted.append(dict(g,offset_start=a,offset_end=b,display_strand=strand,role=role))
    converted.sort(key=lambda d:d['offset_start'])
    return converted,orientation

# A: 6 true coordinate GFF neighbourhoods
axA=fig.add_axes([.055,.467,.618,.395]);axA.set_facecolor('white')
panel_title(axA,'A','GFF-derived HxlR–SnoaL_2 genomic neighbourhoods','Actual CDS lengths, strand and spacing; the five genera are represented by six assemblies')
core=rows[:6]
for k,row in enumerate(core):
    y=5.4-k
    axA.add_patch(Rectangle((-3760,y-.35),7520,.71,facecolor=BG,edgecolor=LINE,lw=.65,zorder=0))
    genes,ori=normalize(row)
    for g in genes:
        if g['offset_end']<-3500 or g['offset_start']>3500:continue
        color=RED if g['role']=='gene1' else GOLD if g['role']=='gene2' else (
            PURPLE if ('recombinase' in g['product'].lower() or 'transposase' in g['product'].lower()) else BLUELIGHT)
        patch=arrow_patch(g['offset_start'],g['offset_end'],y,g['display_strand'],color,h=.21)
        if patch:axA.add_patch(patch)
        mid=(g['offset_start']+g['offset_end'])/2
        if g['role']=='gene1':axA.text(mid,y+.20,'HxlR',fontsize=8.2,color=NAVY,ha='center',va='bottom',weight='bold')
        elif g['role']=='gene2':axA.text(mid,y-.19,'SnoaL_2',fontsize=8.05,color=NAVY,ha='center',va='top',weight='bold')
        elif (g['offset_end']-g['offset_start'])>420 and product_short(g['product']) and abs(mid)<2900:
            axA.text(mid,y-.18,product_short(g['product']),fontsize=6.7,color=PURPLE if color==PURPLE else MUTED,ha='center',va='top')
    axA.text(-5040,y+.06,row['Genus'],fontsize=9.9,color=phy_col[row['Genus']],weight='bold',ha='left',va='center')
    axA.text(-5040,y-.19,f"{row['GCA']}  ·  {row['HMT']}",fontsize=7.4,color=MUTED,ha='left',va='center')
    axA.text(3890,y+.05,f"{row['distance_bp']:,} bp",fontsize=9.4,color=TEAL if row['distance_bp']<=500 else PURPLE,weight='bold',ha='left',va='center')
    axA.text(3890,y-.19,'divergent' if row['distance_bp']<=500 else 'same-strand',fontsize=7.6,color=MUTED,ha='left',va='center')
axA.set_xlim(-5200,5000);axA.set_ylim(-.05,6)
axA.set_xticks([-3000,-2000,-1000,0,1000,2000,3000]);axA.set_xticklabels(['−3','−2','−1','0','1','2','3'])
axA.set_xlabel('Genomic distance from central locus midpoint (kb)',color=MUTED,fontsize=8.5,labelpad=7)
axA.set_yticks([]);axA.tick_params(axis='x',labelsize=7.9,length=2.5,color=LINE,labelcolor=MUTED)
axA.spines[['top','right','left']].set_visible(False);axA.spines['bottom'].set_color(LINE)
axA.axvline(0,color='#dae6e9',lw=.9,zorder=-1)
fig.text(.055,.412,'■',fontsize=9,color=RED);fig.text(.067,.412,'HxlR',fontsize=8.1,color=MUTED)
fig.text(.107,.412,'■',fontsize=9,color=GOLD);fig.text(.119,.412,'SnoaL_2',fontsize=8.1,color=MUTED)
fig.text(.184,.412,'■',fontsize=9,color=BLUELIGHT);fig.text(.196,.412,'other CDS',fontsize=8.1,color=MUTED)
fig.text(.269,.412,'■',fontsize=9,color=PURPLE);fig.text(.281,.412,'recombinase / transposase',fontsize=8.1,color=MUTED)
fig.text(.055,.392,'Flanking labels are putative Prokka product annotations, not confirmed functions or demonstrated horizontal-transfer events.',fontsize=7.5,color=MUTED)

# B: full 67-distance distribution
axB=fig.add_axes([.737,.596,.231,.266]);panel_title(axB,'B','Genomic spacing across all 67 linked genomes','Divergent ≤500 bp in four genera; one loose, same-strand Neisseria example')
order=['Actinomyces','Porphyromonas','Prevotella','Streptococcus','Neisseria']
for i,g in enumerate(order):
    s=[p for p in allpairs if p['genus']==g]
    yy=len(order)-1-i
    for j,p in enumerate(s):
        jitter=((j*7)%13-6)*.014
        axB.scatter(p['bp'],yy+jitter,s=37 if p['bp']>500 else 21,color=phy_col[g],alpha=.74,
                    edgecolor='white',linewidth=.35,zorder=4)
    axB.text(3000,yy,f'n={len(s)}',fontsize=7.4,color=MUTED,ha='right',va='center')
axB.axvline(500,color=PURPLE,linewidth=1,linestyle='--',zorder=1)
axB.axvspan(35,500,color='#ecf5f5',alpha=.88,zorder=0)
axB.text(480,4.65,'500 bp',ha='right',fontsize=7.7,color=PURPLE)
axB.set_xscale('log');axB.set_xlim(38,3200);axB.set_ylim(-.6,4.9)
axB.set_yticks(range(4,-1,-1));axB.set_yticklabels(order,fontsize=8.1,color=NAVY)
axB.set_xticks([50,100,200,500,1000,2000]);axB.set_xticklabels(['50','100','200','500','1k','2k'])
axB.set_xlabel('Intergenic separation (bp; log scale)',color=MUTED,fontsize=8.0,labelpad=5)
axB.spines[['top','right','left']].set_visible(False);axB.spines['bottom'].set_color(LINE)
axB.tick_params(axis='both',length=0,labelsize=8)
fig.text(.741,.550,'66 / 67',fontsize=19,color=RED,weight='bold')
fig.text(.850,.555,'tightly divergent (≤500 bp)',fontsize=8.4,color=MUTED)

# C: two other loci, exact GFF coords
axC=fig.add_axes([.055,.155,.615,.175]);panel_title(axC,'C','Different locus arrangements in individual genera','Alternative linked domain pairs show same-strand organization, with intervening genes')
for i,row in enumerate(rows[6:]):
    y=1.6-i
    gene,ori=normalize(row)
    center=(row['anchor_start']+row['anchor_end'])/2.
    axC.add_patch(Rectangle((-2950,y-.34),5900,.68,facecolor=BG,edgecolor=LINE,lw=.7))
    for g in gene:
        if g['offset_end']<-2900 or g['offset_start']>2900:continue
        if g['role']=='gene1':
            color=RED if row['type'].startswith('PF01638') else GOLD
        elif g['role']=='gene2':
            color=PURPLE if row['type'].endswith('PF22166') else BLUEDARK
        else:
            color=BLUELIGHT
        # Keep flanking gene tails within the common genomic window.
        lo=max(-2890,g['offset_start']);hi=min(2890,g['offset_end'])
        patch=arrow_patch(lo,hi,y,g['display_strand'],color,h=.23)
        if patch:axC.add_patch(patch)
        mid=(g['offset_start']+g['offset_end'])/2
        if g['role']=='gene1':
            label='HxlR' if row['type'].startswith('PF01638') else 'SnoaL_2'
            axC.text(mid,y+.23,label,fontsize=8.2,weight='bold',ha='center',color=NAVY)
        elif g['role']=='gene2':
            label='WHD_ArnR1' if row['type'].endswith('PF22166') else 'Polyketide_cyc2'
            axC.text(mid,y+.23,label,fontsize=8.1,weight='bold',ha='center',color=NAVY)
    axC.text(-4200,y+.065,row['Genus'],fontsize=9.9,color=phy_col[row['Genus']],weight='bold')
    axC.text(-4200,y-.18,f'{row["GCA"]} · {row["HMT"]}',fontsize=7.1,color=MUTED)
    axC.text(3200,y,f'{row["distance_bp"]} bp  ·  same-strand',color=TEAL,fontsize=8.3,ha='left',va='center')
axC.set_xlim(-4250,4800);axC.set_ylim(.06,2.20);axC.set_yticks([])
axC.set_xticks([-2000,-1000,0,1000,2000]);axC.set_xticklabels(['−2','−1','0','1','2'])
axC.spines[['top','right','left']].set_visible(False);axC.spines['bottom'].set_color(LINE)
axC.tick_params(axis='x',length=2.5,labelsize=7.5,colors=MUTED)
axC.set_xlabel('Genomic distance (kb; actual relative coordinates)',fontsize=8.1,color=MUTED)

# D: 8 protein cluster pair identities and ecology caveat
axD=fig.add_axes([.739,.170,.228,.22]);panel_title(axD,'D','Eight alternative protein-cluster pairs','Most exact HxlR/SnoaL_2 pairs do not extend across genera')
order_type=list(enumerate(types))
for i,(idx,r) in enumerate(order_type):
    y=7-i
    axD.barh(y,r['n_genomes'],height=.62,color=TEAL if r['n_genera']>1 else BLUELIGHT)
    axD.text(-1,y,f'{idx+1}',fontsize=8,color=MUTED,ha='center',va='center')
    axD.text(r['n_genomes']+.6,y,str(r['n_genomes']),fontsize=8.4,color=NAVY,ha='left',va='center',weight='bold')
axD.set_yticks([]);axD.set_xticks([0,10,20]);axD.tick_params(axis='x',colors=MUTED,labelsize=7.5,length=2)
axD.set_xlim(0,32);axD.set_ylim(-.3,7.7)
axD.spines[['top','right','left']].set_visible(False);axD.spines['bottom'].set_color(LINE)
axD.set_xlabel('Number of genomes (exact cluster pair)',fontsize=8,color=MUTED)
fig.text(.742,.125,'Only 1/8 pair types spans two genera: Porphyromonas + Prevotella.',fontsize=7.6,color=NAVY)
fig.text(.742,.105,'Linked locus versus plaque preference: P = 0.217 (two-sided).',fontsize=7.7,color=MUTED)
fig.text(.055,.065,'INTERPRETATION  |  Conserved genomic orientation ≠ identical protein orthologues, shared regulation or demonstrated ecological causation.',fontsize=8.65,color=NAVY,weight='medium')
fig.add_artist(Line2D([.055,.967],[.085,.085],transform=fig.transFigure,lw=.8,color=LINE))

caption=("Figure 4 | Conservation of gene orientation across heterogeneous protein families. "
"(A) True gene coordinates, strand orientation and CDS lengths from Prokka GFF files for six representative genomic neighbourhoods containing linked HxlR (PF01638) and SnoaL_2 (PF12680) domains. Tracks are aligned by the midpoint of the linked pair; the Streptococcus HMT-411 track is reverse-complemented for visualization. Coloured arrows encode annotated domain-bearing genes, while grey and purple arrows encode surrounding CDSs and putatively mobile genetic elements. The single Neisseria locus has a 2,124-bp gap and is same-strand, unlike the tighter divergent motifs. (B) Intergenic separation of all 67 linked HxlR–SnoaL_2 genome records from five genera, 66 of which were divergent within 500 bp; the remaining Neisseria example was loosely linked. (C) Actual representative genomic contexts of two distinct locus architectures restricted to Fusobacterium (HxlR–WHD_ArnR1, 38 genomes in the complete scan) and Actinomyces (SnoaL_2–Polyketide_cyc2, 35 genomes). (D) The 67 linked HxlR–SnoaL_2 genomes contain eight MMseqs2 sequence-cluster pair types; just one pair type occurred in more than one genus. This is consistent with recurrent organization despite protein sequence diversity, not necessarily independent evolution, shared regulation, or convergence. Genus-balanced linkage enrichment with plaque specialization was not statistically significant (two-sided P≈0.217). All observed CDS spans and linkage distances are from the accompanying GFF-derived data, while product annotation and functional inference remain putative.")
(P/'Figure4_caption.txt').write_text(caption+'\n')
# data export
pd.DataFrame([{k:v for k,v in r.items() if k!='genes'} for r in rows]).to_csv(P/'Figure4_A_C_reference_loci.tsv',sep='\t',index=False)
pd.DataFrame(allpairs).to_csv(P/'Figure4_B_all67_distances.tsv',sep='\t',index=False)
pd.DataFrame(types).to_csv(P/'Figure4_D_8cluster_pair_types.tsv',sep='\t',index=False)
for ex in ['png','pdf','svg']:
    f=P/('HOMD_Figure4_GFF_genomic_architectures.'+ex)
    fig.savefig(f,dpi=240 if ex=='png' else 210,facecolor='white',bbox_inches='tight',pad_inches=.20)
    print(f.name,f.stat().st_size)
plt.close(fig)
