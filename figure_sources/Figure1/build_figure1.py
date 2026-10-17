#!/usr/bin/env python3
"""Figure 1, HOMD oral habitat preference across 47 taxa and 13 genera.

Inputs are the original the study workflow files transferred with lossless zlib encoding:
- HOMD47_ribosomal_tree.newick
- HOMD47_HMT_manifest.tsv
- No simulated HMT or trait values.
"""
from pathlib import Path
from collections import Counter, defaultdict
import re, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
from matplotlib import colors
from Bio import Phylo

BASE=Path(__file__).resolve().parent
M=pd.read_csv(BASE/'HOMD47_HMT_manifest.tsv', sep='\t')
M['Ecology']=M['Ecology'].map({'p':'plaque','m':'mucosal'})
M['Ecology_tier']=M['Ecology_tier'].map({'H':'High','M':'Moderate'})
tree=Phylo.read(str(BASE/'HOMD47_ribosomal_tree.newick'),'newick')
leaf=list(tree.get_terminals())
leaf_hmt=[z.name.split('|')[0] for z in leaf]
assert len(leaf)==len(M)==47 and len(set(leaf_hmt))==47
assert set(leaf_hmt)==set(M.HMT)
D=M.set_index('HMT').loc[leaf_hmt].reset_index()
assert D.n_genomes.sum()==1121 and D.Genus.nunique()==13
assert sum(D.Ecology=='plaque')==24 and sum(D.Ecology=='mucosal')==23
PHY={
    'Fusobacterium':'Fusobacteriota','Leptotrichia':'Fusobacteriota',
    'Streptococcus':'Bacillota','Gemella':'Bacillota','Veillonella':'Bacillota','Selenomonas':'Bacillota',
    'Actinomyces':'Actinomycetota','Rothia':'Actinomycetota',
    'Campylobacter':'Campylobacterota','Neisseria':'Pseudomonadota',
    'Porphyromonas':'Bacteroidota','Prevotella':'Bacteroidota','Segatella':'Bacteroidota'}
D['Phylum']=D.Genus.map(PHY)
assert D.Phylum.nunique()==6
NAVY='#17334A'; DARK='#1D3040'; SUB='#627887'; FAINT='#E3E8EA'; PALE='#F8FAFA'; PAPER='#FFFFFF'
PLA='#DF795D'; MUC='#168A86'; PLAL='#F8DCD4'; MUCL='#D6EEEA'
PHYC={
 'Fusobacteriota':'#A56688', 'Bacillota':'#9D7BBD','Actinomycetota':'#D2A054',
 'Campylobacterota':'#6697AE','Pseudomonadota':'#5C9C86','Bacteroidota':'#B87069'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                     'ps.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False})

# Canvas in high-res editorial format. Rows y=46,...,0.0 (tree native order).
fig=plt.figure(figsize=(18.2,13.3),dpi=175,facecolor=PAPER)
ax=fig.add_axes([.025,.065,.95,.790]);ax.set_xlim(0,21.8);ax.set_ylim(-1.8,50.5);ax.axis('off')
fig.text(.044,.970,'ORAL HABITAT SPECIALIZATION',size=10.4,color=MUC,weight='bold',va='center')
fig.text(.044,.942,'A phylogenetic atlas of plaque and mucosal specialists',size=21.5,
         color=NAVY,weight='bold',va='center')
fig.text(.044,.918,'Contrasting plaque and mucosal preferences sampled within every genus across six phyla',
         size=11.3,color=SUB,va='center')
# Grand summary drawn in header far right
header=[('13','GENERA'),('47','ORAL HMTs'),('1,121','GENOMES'),('6','PHYLA')]
for i,(n,label) in enumerate(header):
 x=.575+i*.101
 fig.text(x,.959,n,color=NAVY,size=15,weight='bold',ha='center',va='center')
 fig.text(x,.939,label,color=SUB,size=7.2,weight='bold',ha='center',va='center')
fig.add_artist(Line2D([.044,.958],[.902,.902],color='#D7E0E2',lw=.9,transform=fig.transFigure))

# Panel A dimensions
x_tree0=.68;x_tree1=5.48
x_taxon=5.90;x_id=10.05
x_plot0=10.74;x_plot1=15.46;x_zero=12.98
x_genomes=15.95
# Panel B x 17..21.7
bx0=17.05;bxbar=19.23;bw=.246;bxmax=21.38

# Column captions
ax.text(.32,49.3,'A',size=21,color=NAVY,weight='bold',ha='center')
ax.text(.74,49.4,'18-marker phylogeny',size=10.0,color=NAVY,weight='bold')
ax.text(x_taxon,49.4,'Oral taxon',size=10.0,color=NAVY,weight='bold')
ax.text(x_id,49.4,'HMT',size=9.3,color=NAVY,weight='bold')
ax.text(13.10,49.45,'Ecological affinity',size=10,color=NAVY,weight='bold',ha='center')
ax.text(x_genomes+.13,49.4,'Genomes',size=8.7,color=NAVY,weight='bold',ha='center')
ax.text(bx0,49.3,'B',size=21,color=NAVY,weight='bold')
ax.text(bx0+.52,49.45,'Matched genera',size=10,color=NAVY,weight='bold')
ax.text(bx0+.52,47.8,'Plaque / mucosal HMTs',size=8.1,color=SUB)
# separators
ax.plot([16.59,16.59],[-.45,47.5],color='#D8E0E2',lw=1.05)
ax.plot([.35,16.25],[47.78,47.78],color='#CCD8DC',lw=1.0)
ax.plot([bx0,21.6],[47.0,47.0],color='#CCD8DC',lw=1.0)

# Y locations, tree root at depth0, tips by tree's given terminal order
ypos={c:len(leaf)-i-1 for i,c in enumerate(leaf)}
y_node={}
for c in tree.find_clades(order='postorder'):
 y_node[c]=ypos[c] if c in ypos else np.mean([y_node[d] for d in c.clades])
depths=tree.depths(); scale=(x_tree1-x_tree0)/max(depths.values())
for c in tree.find_clades(order='preorder'):
 if not c.clades:continue
 cx=x_tree0+depths[c]*scale
 yy=[y_node[d] for d in c.clades]
 ax.plot([cx,cx],[min(yy),max(yy)],color='#879BA7',lw=.9,zorder=2,solid_capstyle='round')
 for ch in c.clades:
  xx=x_tree0+depths[ch]*scale
  ax.plot([cx,xx],[y_node[ch],y_node[ch]],color='#879BA7',lw=.88,zorder=2,solid_capstyle='round')
# phylum patches, contiguous blocks. genus separators subtle.
gen_idx=defaultdict(list);phy_idx=defaultdict(list)
for i,r in D.iterrows():
 y=46-i;gen_idx[r.Genus].append(y);phy_idx[r.Phylum].append(y)
# phylum background bands and strips
for p,yy in phy_idx.items():
 top=max(yy);bottom=min(yy)
 ax.add_patch(Rectangle((.16,bottom-.49),16.12,top-bottom+.98,
   facecolor=colors.to_rgba(PHYC[p],.041),edgecolor='none',zorder=-1))
 ax.add_patch(Rectangle((.19,bottom-.49),.105,top-bottom+.98,
   facecolor=PHYC[p],edgecolor='none',zorder=1))
 # thin label (only once per phylum, near center rotated)
 ax.text(.075,np.mean(yy),('Campyl.' if p=='Campylobacterota' else ('Pseudom.' if p=='Pseudomonadota' else p)),color=PHYC[p],rotation=90,ha='center',va='center',size=7.0 if len(yy)>3 else 6.1,weight='bold')
for g,yy in gen_idx.items():
 top=max(yy);bottom=min(yy)
 # genus distinction line
 if bottom>0:
  ax.plot([.32,16.30],[bottom-.52,bottom-.52],color='#E1E6E7',lw=.65,zorder=0)
 # mark genus at phylogram rightmost via small slash? no labels repeated in taxon col

# Ecological affinity axis [-10,12], 0 in middle-ish
minv,maxv=-10.,12.
zero=x_plot0+(0-minv)/(maxv-minv)*(x_plot1-x_plot0)
def xx(value):return x_plot0+(value-minv)/(maxv-minv)*(x_plot1-x_plot0)
for t in [-8,-4,0,4,8,12]:
 xt=xx(t)
 ax.plot([xt,xt],[47.5,-.5],color='#D3DEE1' if t==0 else '#EAF0F1',lw=1.15 if t==0 else .55,zorder=0)
 ax.text(xt,48.22,str(t),size=7.6,ha='center',color=SUB if t else NAVY,weight='bold' if t==0 else 'normal')
ax.text(xx(-5.8),47.75,'Mucosa',size=7.9,color=MUC,ha='center',weight='bold')
ax.text(xx(6),47.75,'Plaque',size=7.9,color=PLA,ha='center',weight='bold')
# data rows
for i,r in D.iterrows():
 y=46-i
 col=PLA if r.Ecology=='plaque' else MUC
 # dotted lead line to taxa labels, so branch lengths remain uncluttered
 endx=x_tree0+depths[leaf[i]]*scale
 ax.plot([endx,x_taxon-.13],[y,y],color='#C5D0D4',lw=.75,linestyle=(0,(1.4,2.6)),zorder=1)
 # ecological class itself in phylogenetic tip marker
 ax.scatter([endx],[y],s=24,c=[col],edgecolors='white',linewidths=.7,zorder=7)
 # genus+species. printed full for legibility
 name=f'{r.Genus} {r.Species}'
 if r.Genus=='Segatella':name=f'{r.Genus} {r.Species}'
 ax.text(x_taxon,y,name,ha='left',va='center',size=8.0,
         color=DARK,fontstyle='italic',zorder=5)
 ax.text(x_id,y,r.HMT.replace('HMT-',''),ha='left',va='center',size=7.5,color=SUB,zorder=5)
 # age of ecology estimation is not a new inference; high confidence solid circle, moderate hollow center
 xv=xx(r.Median_Y)
 ax.plot([zero,xv],[y,y],color=col,lw=2.4,alpha=.75,solid_capstyle='round',zorder=3)
 if r.Ecology_tier=='High':
  ax.scatter([xv],[y],s=44,facecolors=col,edgecolors='white',linewidths=.85,zorder=9)
 else:
  ax.scatter([xv],[y],s=48,facecolors='white',edgecolors=col,linewidths=1.65,zorder=9)
 ax.text(x_genomes+.10,y,str(r.n_genomes),ha='center',va='center',size=7.3,color=SUB)

# phylogenetic scale bar: branch 0.1 substitution/site
bar=.10*scale
ax.plot([x_tree0+.15,x_tree0+.15+bar],[48.3,48.3],color='#849AA6',lw=1.55)
ax.plot([x_tree0+.15,x_tree0+.15],[48.18,48.42],color='#849AA6',lw=.95)
ax.plot([x_tree0+.15+bar,x_tree0+.15+bar],[48.18,48.42],color='#849AA6',lw=.95)
ax.text(x_tree0+.15+bar/2,48.6,'0.1 subs. / site',size=6.8,ha='center',color=SUB)

# B panel, grouped as tree order, if any genus order interleaving use first appearance
genus_order=list(dict.fromkeys(D.Genus))
assert len(genus_order)==13
# rows spaced approx 3.2 high = 43, full extent
ys=np.linspace(44.8,6.0,len(genus_order))
for z,(g,y) in enumerate(zip(genus_order,ys)):
 rr=D[D.Genus==g]; np_=sum(rr.Ecology=='plaque');nm_=sum(rr.Ecology=='mucosal');ph=PHY[g]
 if z%2==0:
  ax.add_patch(Rectangle((bx0-.03,y-1.22),4.5,2.43,color='#F7F9F9',lw=0,zorder=0))
 ax.scatter([bx0+.13],[y],c=[PHYC[ph]],s=32,zorder=4)
 ax.text(bx0+.37,y,g,size=8.9,va='center',color=NAVY,weight='medium')
 # two-tone segmented count bar, each segment separate square
 for u in range(np_):
  ax.add_patch(FancyBboxPatch((bxbar+u*bw,y-.37),bw-.034,.72,
              boxstyle='round,pad=0,rounding_size=.06',facecolor=PLA,edgecolor='none'))
 for u in range(nm_):
  ax.add_patch(FancyBboxPatch((bxbar+(np_+u)*bw,y-.37),bw-.034,.72,
              boxstyle='round,pad=0,rounding_size=.06',facecolor=MUC,edgecolor='none'))
 ax.text(bxbar+7*bw+.12,y,f'{np_} : {nm_}',size=7.7,ha='left',va='center',color=SUB)
ax.text(bx0+.39,3.6,'Every genus contains both ecological classes',size=8.4,weight='bold',color=NAVY)
ax.text(bx0+.39,2.2,'Counts are HMTs, not genomes',size=8.2,color=SUB)

# bottom key in fig coordinates
fig.add_artist(Line2D([.044,.956],[.058,.058],color='#D7E0E2',lw=.8,transform=fig.transFigure))
fig.text(.045,.037,'Ecological phenotype:',size=8,color=SUB,weight='bold',va='center')
fig.text(.155,.037,'●',size=13,color=PLA,va='center')
fig.text(.169,.037,'Plaque (n=24)',size=8,color=DARK,va='center')
fig.text(.242,.037,'●',size=13,color=MUC,va='center')
fig.text(.256,.037,'Mucosal (n=23)',size=8,color=DARK,va='center')
fig.text(.350,.037,'●',size=12,color=SUB,va='center')
fig.text(.365,.037,'High confidence',size=8,color=DARK,va='center')
fig.text(.447,.037,'○',size=14,color=SUB,va='center')
fig.text(.463,.037,'Moderate confidence',size=8,color=DARK,va='center')
fig.text(.585,.037,'Score = median log2[(SUBP + SUPP + ε)/(BMU + TDO + ε)]',size=8,color=SUB,va='center')
fig.text(.045,.015,'Phylogram from 18 conserved marker proteins; display-rooted for illustration, with no ancestral-state reconstruction. Each HMT is one ecological unit.',
         size=7.9,color=SUB,va='center')
fig.text(.957,.015,'HOMD  ·  13-genus comparative cohort',size=7.8,color=NAVY,ha='right',va='center',weight='bold')
for ext in ['pdf','svg','png']:
 fn=BASE/('HOMD_Figure1_13genus_phylogenetic_ecological_atlas.'+ext)
 fig.savefig(fn,facecolor='white',dpi=220 if ext=='png' else 175,bbox_inches=None)
 print('WRITTEN',fn,fn.stat().st_size)
plt.close(fig)
D.to_csv(BASE/'HOMD_Figure1_plotted_HMT_order.tsv',index=False,sep='\t')
# Structured genus counts, family labels, manuscript caption
counts=D.groupby(['Genus','Phylum','Ecology'],sort=False).size().unstack(fill_value=0).reset_index()
counts['n_HMT']=counts.plaque+counts.mucosal
counts.to_csv(BASE/'HOMD_Figure1_13genus_counts.tsv',sep='\t',index=False)
caption='''Fig. 1 | Oral habitat specialization is represented across 13 bacterial genera spanning six phyla. (a) Branch-length phylogram inferred from 18 conserved marker proteins for 47 Human Microbial Taxa (HMTs), ordered as in the input tree and coloured by their plaque-associated (coral) or mucosa-associated (teal) ecological classification. The horizontal dot-bar shows the HMT-specific median log2 plaque-to-mucosal site-affinity score derived from publicly available HOMD abundance summaries (SUBP and SUPP versus BMU and TDO; a small pseudocount was applied). Solid affinity symbols designate high-confidence classifications (n=32); open symbols designate moderate-confidence classifications (n=15). Numbers at right indicate the reference genomes supporting each HMT (1,121 in total), not independently sampled ecological observations. (b) Numbers of HMTs assigned to the plaque and mucosal classes in each of the 13 genera; each genus includes representatives of both classes. Phylum colours are used only as taxonomic guides. The tree is display-rooted and no ancestral habitat states or independent evolutionary transitions are inferred. Habitat categories are intentionally enriched for strong specialists rather than being an unbiased census of oral bacteria.'''
(BASE/'HOMD_Figure1_draft_caption.txt').write_text(caption+'\n')
