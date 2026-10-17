#!/usr/bin/env python3
"""Figure 3, HOMD 13-genus paper. Runs from Figure3_source_aggregate.json.

Data provenance: the study workflow study13:
11_manuscript_completion_20261013/candidate_Pfam_protein_subfamily_summary.tsv
11_manuscript_completion_20261013/candidate_Pfam_MMseqs_family_crossgenus.tsv
11_manuscript_completion_20261013/candidate_Pfam_fullprotein_orthogroup_architecture.tsv
10_additional_analyses_20261012/genus_balanced_repertoire_omnibus_six_layers.tsv
Protein sequence clusters: MMseqs2 easy-linclust 50% id, 80% cov, cov-mode 0.

CAUTION: X5 sequence clusters are not experimentally validated orthologues;
Pfam-derived functions remain putative, domain architectures in this graphic
are co-annotations, not genomic or N-to-C domain coordinate diagrams.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
import pandas as pd

OUT=Path(__file__).resolve().parent
J=json.loads((OUT/'Figure3_source_aggregate.json').read_text())
D=pd.DataFrame(J['domains']).set_index('Pfam')
OM=pd.DataFrame(J['omnibus'])

# palette coherent with other panels in manuscript
INK='#17354B'; NAVY='#12334B'; SUB='#5B7280'; GREY='#D3E0E3'; GRID='#E4ECEE'; LIGHT='#F6F9FA'
CORAL='#D6755D'; CORAL_D='#A9523D'; TEAL='#138D95'; TEAL_D='#09747F'; GOLD='#DEAF68'
PHY={'Fusobacterium':'Fusobacteriota','Leptotrichia':'Fusobacteriota',
     'Streptococcus':'Bacillota','Gemella':'Bacillota','Veillonella':'Bacillota','Selenomonas':'Bacillota',
     'Actinomyces':'Actinomycetota','Rothia':'Actinomycetota',
     'Campylobacter':'Campylobacterota','Neisseria':'Pseudomonadota',
     'Porphyromonas':'Bacteroidota','Prevotella':'Bacteroidota','Segatella':'Bacteroidota'}
PHYC={'Fusobacteriota':'#AE7196','Bacillota':'#927EC5','Actinomycetota':'#D69F4D',
     'Campylobacterota':'#5C9EBB','Pseudomonadota':'#63A591','Bacteroidota':'#BF737A'}
GENERA=['Fusobacterium','Leptotrichia','Streptococcus','Gemella','Veillonella','Selenomonas',
        'Actinomyces','Rothia','Campylobacter','Neisseria','Porphyromonas','Prevotella','Segatella']
assert set(GENERA)==set(J['genus_names'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.unicode_minus':False,
                     'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(19.2,13.0),dpi=200,facecolor='white')
fig.text(.045,.967,'ORAL HABITAT SPECIALIZATION',color=TEAL,weight='bold',fontsize=10)
fig.text(.045,.933,'Recurrent domains, divergent protein contexts',color=NAVY,weight='bold',fontsize=22)
fig.text(.045,.906,'Cross-lineage habitat associations are embedded in heterogeneous sequence families and protein architectures',
         color=SUB,fontsize=10.5)
fig.add_artist(Line2D([.045,.960],[.893,.893],transform=fig.transFigure,color=GREY,lw=1.2))

def title(x,y,letter,headline,description):
    fig.text(x,y,letter,fontsize=20,weight='bold',color=NAVY)
    fig.text(x+.037,y+.003,headline,fontsize=13.0,weight='bold',color=NAVY)
    fig.text(x+.037,y-.020,description,fontsize=8.8,color=SUB)

def vline(x,y1,y2):fig.add_artist(Line2D([x,x],[y1,y2],color=GREY,lw=1.05,transform=fig.transFigure))
def hline(y):fig.add_artist(Line2D([.045,.961],[y,y],color=GREY,lw=.95,transform=fig.transFigure))

# PANEL A — quantify fragmentation across sequence clusters
# all 13 tested PFAM annotations, including the linked SnoaL and polyketide-context domains
lab={
 'PF01638':'HxlR regulator','PF04072':'LCM','PF13676':'TIR_2','PF12680':'SnoaL_2',
 'PF01745':'IPT / MiaA','PF18075':'FtsX ECD','PF13413':'HTH_25','PF22166':'WHD_ArnR1',
 'PF07366':'SnoaL','PF13577':'SnoaL_4','PF10604':'Polyketide_cyc2','PF01337':'Barstar','PF14470':'bPH_3'}
order=['PF01638','PF04072','PF13676','PF12680','PF01745','PF18075','PF13413','PF22166',
       'PF07366','PF13577','PF10604','PF01337','PF14470']
title(.045,.849,'A','A domain does not define a single protein family',
      'Across genera, sequence clusters span fewer lineages than their parent domain annotations')
axA=fig.add_axes([.184,.512,.245,.274])
axA.set_xlim(0,13.7);axA.set_ylim(-.55,12.55);axA.invert_yaxis()
axA.axvspan(0,13.05,color='#F4F8F8',zorder=0)
for y in range(13):axA.hlines(y,0,13,color='white',lw=1.8,zorder=1)
for y,pf in enumerate(order):
    z=D.loc[pf]
    nt=int(z.n_genera); nn=int(z.max_genera_in_a_single_cluster)
    axA.plot([nn,nt],[y,y],color='#AAC2CA',lw=3.0,zorder=2,solid_capstyle='round')
    axA.scatter([nt],[y],s=80,facecolors=TEAL,edgecolors='white',linewidths=.8,zorder=4)
    axA.scatter([nn],[y],s=92,facecolors=CORAL,edgecolors='white',linewidths=.8,zorder=5,marker='D')
    axA.text(-.33,y,lab[pf],fontsize=8.7,color=NAVY,ha='right',va='center',clip_on=False)
axA.set_xticks([0,3,6,9,13]);axA.set_xlabel('Number of genera',color=SUB,fontsize=9)
axA.tick_params(axis='both',colors=SUB,length=0,labelsize=8.8)
axA.set_yticks([])
for sp in axA.spines.values():sp.set_visible(False)
fig.text(.197,.798,'◆ Broadest shared sequence cluster',color=CORAL_D,fontsize=8.4)
fig.text(.378,.798,'● Domain observed',color=TEAL_D,fontsize=8.4)

# right aligned X5 count encoded as lollipop lengths, not arbitrary ordinal rank
axAc=fig.add_axes([.459,.512,.072,.274],sharey=axA)
axAc.set_xlim(0,84);axAc.set_yticks([])
for y,pf in enumerate(order):
    n=int(D.loc[pf,'n_distinct_MMseqs_X5_clusters'])
    axAc.barh(y,n,height=.43,edgecolor='none',color='#A2B7C7')
    axAc.text(n+1.0,y,str(n),ha='left',va='center',color=NAVY,fontsize=8.2,fontweight='bold')
axAc.set_xticks([0,40,80]);axAc.tick_params(length=0,labelsize=8.0,colors=SUB)
for sp in axAc.spines.values():sp.set_visible(False)
axAc.set_xlabel('Sequence clusters',color=SUB,fontsize=8.3)
fig.text(.462,.798,'MMseqs X5',fontsize=8.1,color=SUB)

# PANEL B — linked constellations: each row a real X5 sequence cluster with coloured genus marks
vline(.552,.499,.871)
title(.571,.849,'B','Protein family mosaics across genera',
      'Selected sequence clusters; coloured circles indicate genera carrying each cluster')
axB=fig.add_axes([.625,.520,.327,.293])
axB.set_xlim(-.5,12.5);axB.set_ylim(-.7,14.35);axB.invert_yaxis()
axB.set_xticks(range(13));axB.set_xticklabels([x[:6]+'.' if len(x)>8 else x for x in GENERA],
                                         rotation=42,ha='left',fontsize=7.0,color=SUB)
axB.xaxis.set_ticks_position('top');axB.tick_params(axis='x',length=0,pad=4)
axB.set_yticks([])
for sp in axB.spines.values():sp.set_visible(False)
for x,g in enumerate(GENERA):
    if x in (2,6,8,9,10):axB.axvline(x-.5,color='white',lw=3,zorder=2)
    axB.axvspan(x-.44,x+.44,color=to_rgba(PHYC[PHY[g]],.095),zorder=-5)
rows=[]
for pf,num in [('PF01638',6),('PF04072',4),('PF13676',3)]:
    rows.extend([(pf,i,entry) for i,entry in enumerate(J['network_top_clusters'][pf]['clusters'][:num])])
# room for boundaries and intentional breaks
positions=[0,1,2,3,4,5,6.6,7.6,8.6,9.6,11.2,12.2,13.2]
for (pf,i,e),y in zip(rows,positions):
    gs=e['genera'];xx=np.array([GENERA.index(g) for g in gs]);pcol=CORAL if pf=='PF01638' else TEAL if pf=='PF04072' else '#8F7FBE'
    axB.hlines(y,np.min(xx),np.max(xx),color=to_rgba(pcol,.38),lw=2.4,zorder=0)
    for x in xx:
        axB.scatter([x],[y],s=70,facecolors=PHYC[PHY[GENERA[x]]],edgecolors='white',linewidths=.8,zorder=4)
    axB.text(-.54,y,'',fontsize=8)
    axB.text(12.75,y,f"{e['n_proteins']} proteins",fontsize=7.2,color=SUB,ha='left',va='center',clip_on=False)
    # actual sequence family accession stored in TSV, simplified index display in figure
    axB.text(-1.1,y,f'{i+1:02d}',ha='right',va='center',fontsize=8.2,color=pcol,clip_on=False)
for y0,label,color in [(-.45,'HxlR  ·  73 clusters',CORAL_D),(6.25,'LCM  ·  15 clusters',TEAL_D),(10.87,'TIR_2  ·  69 clusters','#7962AE')]:
    axB.text(-2.15,y0,label,color=color,fontsize=9.4,ha='left',weight='bold',clip_on=False)
axB.text(-2.15,14.9,'Rows: top clusters ranked by genus breadth and protein count; omitted clusters not shown.',
    color=SUB,fontsize=8.1,ha='left',clip_on=False)

hline(.481)
# PANEL C — exact co-occurrence proportions (Pfam co-annotations, not sequence or order)
title(.045,.439,'C','Protein-domain combinations reveal distinct organizational patterns',
      'Annotation co-occurrence for each protein; not a positional map of N–C terminal domains')
axC=fig.add_axes([.19,.166,.354,.215]);axC.set_xlim(0,100);axC.set_ylim(-.55,4.48);axC.invert_yaxis()
cohort=[('PF04072','LCM','PF04072 alone',CORAL),
        ('PF01638','HxlR','PF01638 alone',CORAL),
        ('PF13676','TIR_2','PF13676 alone','#8F7FBE'),
        ('PF01745','IPT / MiaA','PF01715 + PF01745',TEAL),
        ('PF18075','FtsX ECD','PF02687 + PF18075',TEAL)]
readout=[]
for y,(pf,name,pattern,color) in enumerate(cohort):
    arr=J['architectures'][pf];n=int(D.loc[pf,'n_proteins'])
    primary=next(a['count'] for a in arr if a['architecture']==pattern.replace(' + ',';')) if False else None
    lookup={'PF04072':'PF04072','PF01638':'PF01638','PF13676':'PF13676',
            'PF01745':'PF01715;PF01745','PF18075':'PF02687;PF18075'}
    primary=next(a['count'] for a in arr if a['architecture']==lookup[pf]);other=n-primary
    p=primary/n*100
    axC.barh(y,p,height=.43,color=color,edgecolor='none',zorder=3)
    axC.barh(y,100-p,height=.43,left=p,color='#DDE7E9',edgecolor='none')
    axC.text(-3,y,name,color=NAVY,fontsize=10,weight='bold',va='center',ha='right',clip_on=False)
    axC.text(1.0,y-.33,pattern,color=SUB,fontsize=7.4,va='center',ha='left',clip_on=False)
    axC.text(102,y,f'{p:.1f}%  ({primary:,}/{n:,})',color=NAVY,fontsize=9.3,va='center',ha='left',clip_on=False)
    readout.append(dict(Pfam=pf,main_coannotation=pattern,n_primary=primary,n_proteins=n,percent=p,
                        other_architectures=other,n_distinct_architectures=int(D.loc[pf,'n_distinct_domain_architectures'])))
axC.set_xticks([0,25,50,75,100]);axC.tick_params(axis='x',length=0,labelsize=8,color=SUB)
axC.grid(axis='x',color=GRID,lw=.8,zorder=0)
axC.set_yticks([])
for sp in axC.spines.values():sp.set_visible(False)
axC.set_xlabel('Percentage of proteins with indicated domain co-annotation',color=SUB,fontsize=9)
fig.text(.19,.397,'Dominant combination → colour    ·    Other co-annotations → grey',color=SUB,fontsize=8.35)
fig.text(.086,.106,'LCM and HxlR are often standalone Pfam annotations; TIR_2 is more heterogeneous.',color=SUB,fontsize=8.3)

# PANEL D — absence of whole-repertoire directional signature
vline(.552,.087,.465)
title(.573,.439,'D','No shared genome-wide functional programme',
      'Directional coherence tested across six annotation layers, equal genus weights')
axD=fig.add_axes([.695,.165,.253,.213])
ordr=['PFAM','COG','EC','KMOD','KPATH','X5']
for y,layer in enumerate(ordr):
    r=OM.loc[OM.layer==layer].iloc[0];p=float(r.perm_p)
    axD.barh(y,p,height=.46,color='#AABEC8',zorder=2)
    axD.scatter(p,y,s=48,color=TEAL,edgecolors='white',lw=.7,zorder=4)
    axD.text(-.035,y,{'PFAM':'Pfam','COG':'COG','EC':'EC','KMOD':'KEGG modules','KPATH':'KEGG pathways','X5':'Sequence clusters'}[layer],
             fontsize=9.25,ha='right',va='center',color=NAVY,clip_on=False)
    axD.text(1.03,y,f'p = {p:.3f}',fontsize=8.3,ha='left',va='center',color=SUB,clip_on=False)
axD.set_xlim(0,1);axD.set_ylim(-.65,5.65);axD.invert_yaxis()
axD.axvline(.05,color=CORAL,ls='--',lw=1.2,zorder=4)
axD.set_xticks([0,.25,.5,.75,1]);axD.tick_params(axis='x',labelsize=8,colors=SUB,length=0)
axD.set_xlabel('Within-genus permutation P-value',fontsize=9,color=SUB)
axD.set_yticks([])
axD.set_axisbelow(True);axD.grid(axis='x',color=GRID,lw=.85)
for sp in axD.spines.values():sp.set_visible(False)
fig.text(.695,.392,'Dashed line: P = 0.05',color=SUB,fontsize=8.4)
fig.text(.579,.108,'Omnibus nonsignificance does not exclude a sparse subset of recurrent functions.',color=SUB,fontsize=8.35)

# Legend and explicit caveat
fig.add_artist(Line2D([.045,.960],[.080,.080],transform=fig.transFigure,color=GREY,lw=.8))
fig.text(.046,.064,'DESCRIPTIVE GENOMIC EVIDENCE',weight='bold',color=NAVY,fontsize=8.4)
fig.text(.046,.044,'MMseqs2 clusters: ≥50% sequence identity and ≥80% coverage. Counts are among 1,121 genome annotations; genomes are nested within 47 HMTs.',
         color=SUB,fontsize=8.25)
fig.text(.046,.026,'Protein sequence clusters are not established orthologues, catalytic functions or independent evolutionary gains. No individual candidate survived genome-wide maxT correction.',
         color=SUB,fontsize=8.25)

prefix=OUT/'HOMD_Figure3_domain_to_protein_context'
for ext in ['svg','pdf','png']:
    fig.savefig(str(prefix)+'.'+ext,dpi=220 if ext=='png' else 180,
                bbox_inches=None,facecolor='white')
plt.close(fig)

# Flat source tables for all chart panels, preserve original precision
pd.DataFrame([dict(Pfam=pf,**{k:v for k,v in D.loc[pf].to_dict().items()}) for pf in order]).to_csv(OUT/'Figure3_A_all_13_domain_cluster_counts.tsv',sep='\t',index=False)
rr=[]
for pf in ['PF01638','PF04072','PF13676','PF22166']:
    for i,q in enumerate(J['network_top_clusters'][pf]['clusters']):
        rr.append(dict(Pfam=pf,rank=i+1,cluster_id=q['id'],n_proteins=q['n_proteins'],n_genera=q['n_genera'],genera=';'.join(q['genera'])))
pd.DataFrame(rr).to_csv(OUT/'Figure3_B_selected_cluster_genus_memberships.tsv',sep='\t',index=False)
pd.DataFrame(readout).to_csv(OUT/'Figure3_C_domain_architecture_proportions.tsv',sep='\t',index=False)
OM.to_csv(OUT/'Figure3_D_omnibus_repertoire_tests.tsv',sep='\t',index=False)
caption='''Figure 3 | Recurrent domain annotations encompass heterogeneous protein sequence clusters and domain combinations. (a) For 13 Pfam families studied in the candidate/locus analysis, a domain can occur across a larger number of oral bacterial genera (teal circles) than any single sequence cluster (coral diamonds). Grey bars show the number of sequence clusters within that Pfam family. (b) Selected high-breadth MMseqs2 clusters bearing HxlR, LCM or TIR_2 annotations, arranged by the taxa in which the corresponding protein sequences occurred. Dot colours correspond to bacterial phyla and line segments connect the observed genus membership of each cluster; intervening taxa are not implied to contain the cluster. The plots display only the highest-ranked clusters for readability, not all clusters. Clusters derive from MMseqs2 easy-linclust at 50% sequence identity and 80% coverage and therefore do not necessarily represent distinct evolutionary orthologues or biochemical specificities. (c) Dominant combinations of Pfam annotations among proteins carrying five representative domains. Percentages are calculated over the protein counts in the source annotation data. Architectures are represented as co-annotations and do not establish N-to-C domain order or transcriptional units. (d) Genus-balanced omnibus tests of directional functional-repertoire coherence across six annotation layers: Pfam (4,850 features), COG (1,575), EC (246), KEGG modules (251), KEGG pathways (126) and 50% identity protein clusters (5,745). All displayed permutation P values are nonsignificant. These results do not exclude sparse repeated signals within particular feature families, and the current 14-feature shortlist did not pass the genome-wide maxT multiple-testing correction. Genomic observations were aggregated within HMT, with no individual reference assembly treated as an independent ecological replicate.'''
(OUT/'Figure3_caption.txt').write_text(caption+'\n')
print('DONE',[(x.name,x.stat().st_size) for x in OUT.glob('HOMD_Figure3*')])
