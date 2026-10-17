#!/usr/bin/env python3
"""Figure 2: fully data-backed HOMD 13-genus functional atlas.

Inputs: Figure2_plotting_data.json, exported from the original
crossgenus_allfeatures_raw_permutation_results.tsv and
crossgenus_allfeatures_genome_adjusted_permutation_results.tsv,
crossgenus_feature_integrated_raw_adjusted_annotated.tsv, and
feature_atlas_41feature_pair_data.json.
Outputs: SVG, PDF, PNG, panel-level CSV/TSV and figure caption.
"""
from pathlib import Path
from collections import Counter
import json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch
from matplotlib.lines import Line2D

P=Path(__file__).resolve().parent
D=json.loads((P/"Figure2_plotting_data.json").read_text())
G=list(D["genera"])
atlas=D["atlas"]
focal=D["focal"]
assert len(G)==13 and len(atlas)==41 and len(focal)==14 and D["n_features"]==13917
assert D["n_raw_nominal"]==119 and D["n_both_nominal"]==23
names=["Fusobacterium","Leptotrichia","Streptococcus","Gemella","Veillonella","Selenomonas",
       "Actinomyces","Rothia","Campylobacter","Neisseria","Porphyromonas","Prevotella","Segatella"]
assert set(names)==set(G)
phylum={
"Fusobacterium":"Fusobacteriota","Leptotrichia":"Fusobacteriota",
"Streptococcus":"Bacillota","Gemella":"Bacillota","Veillonella":"Bacillota","Selenomonas":"Bacillota",
"Actinomyces":"Actinomycetota","Rothia":"Actinomycetota","Campylobacter":"Campylobacterota",
"Neisseria":"Pseudomonadota","Porphyromonas":"Bacteroidota","Prevotella":"Bacteroidota","Segatella":"Bacteroidota"}
colorsphy={"Fusobacteriota":"#A66387","Bacillota":"#9C81BA","Actinomycetota":"#D69F48",
"Campylobacterota":"#5498B8","Pseudomonadota":"#579E87","Bacteroidota":"#BF7373"}
nHMT={"Fusobacterium":5,"Leptotrichia":2,"Streptococcus":7,"Gemella":3,"Veillonella":5,
"Selenomonas":2,"Actinomyces":4,"Rothia":2,"Campylobacter":2,"Neisseria":4,
"Porphyromonas":2,"Prevotella":7,"Segatella":2}
assert sum(nHMT.values())==47

NAVY="#16344B"; TEXT="#23394C"; MID="#617784"; LINE="#D7E1E4"; BG="#F6F9F9"
PLA="#D86E58"; MUC="#168D92"; PALE_P="#FCEBE6"; PALE_M="#E5F2F2"
CMAP=mcolors.LinearSegmentedColormap.from_list("habitat",[MUC,"#B7DCDD","#FBFCFC","#F2B9A8","#B72E40"],N=256)
NORM=mcolors.TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.7,"pdf.fonttype":42,
 "ps.fonttype":42,"svg.fonttype":"none","axes.unicode_minus":False})
fig=plt.figure(figsize=(19.2,13.8),dpi=175,facecolor="white")
fig.text(.042,.971,"ORAL HABITAT SPECIALIZATION",color=MUC,weight="bold",fontsize=10.2)
fig.text(.042,.941,"A cross-lineage atlas of habitat-associated genomic functions",color=NAVY,
         weight="bold",fontsize=19.1)
fig.text(.042,.916,"Genus-balanced genomic contrasts across 47 oral taxa, 13 genera and six phyla",
         color=MID,fontsize=10.4)
fig.add_artist(Line2D([.042,.964],[.904,.904],transform=fig.transFigure,lw=.9,color=LINE))
def tag(ax,letter,title,subtitle=""):
 ax.text(0,1.035,letter,color=NAVY,fontsize=19,fontweight="bold",va="bottom",transform=ax.transAxes)
 ax.text(.048,1.051,title,color=NAVY,fontsize=11.4,fontweight="bold",va="bottom",transform=ax.transAxes)
 if subtitle:ax.text(.048,1.013,subtitle,color=MID,fontsize=8.15,va="bottom",transform=ax.transAxes)
def panel_line(x0,y0,x1,y1):
 fig.add_artist(Line2D([x0,x1],[y0,y1],transform=fig.transFigure,color=LINE,lw=.8))
def box(ax,x,y,w,h,main,sub,edge="#D6E2E4",fill="#F7FAFB",maincol=NAVY,size=9):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.009,rounding_size=.025",
         fc=fill,ec=edge,lw=.78,transform=ax.transAxes))
 ax.text(x+w/2,y+h*.64,main,fontsize=size,fontweight="bold",color=maincol,
          ha="center",va="center",transform=ax.transAxes)
 ax.text(x+w/2,y+h*.28,sub,fontsize=7.0,color=MID,ha="center",va="center",transform=ax.transAxes)
def arrow(ax,x1,y1,x2,y2):
 ax.annotate("",xy=(x2,y2),xytext=(x1,y1),xycoords=ax.transAxes,textcoords=ax.transAxes,
             arrowprops=dict(arrowstyle="->",color="#8CADB8",lw=1.4,shrinkA=2,shrinkB=2))

# A. Computational pipeline, faithfully reflects six HOMD sources and actual filtering
axA=fig.add_axes([.042,.665,.616,.206]);axA.axis("off")
tag(axA,"A","Computational design and selection","Ecological evidence and genome annotations converge at the HMT level")
box(axA,.022,.594,.282,.338,"6 eHOMD abundance summaries","SUBP/SUPP versus BMU/TDO")
box(axA,.357,.594,.282,.338,"1,121 annotated genomes","Pfam, COG, EC, RHEA, KEGG, MMseqs")
box(axA,.692,.594,.286,.338,"47 HMTs · 13 genera","13,917 nonredundant prevalence profiles",size=8.4)
arrow(axA,.642,.755,.682,.755);arrow(axA,.163,.591,.163,.464)
box(axA,.022,.135,.282,.323,"Equal-weight genus contrasts","Plaque minus mucosa; one effect / genus")
box(axA,.357,.135,.282,.323,"4,096 within-genus permutations","Repeated after size, GC and CDS adjustment")
box(axA,.692,.135,.286,.323,"119  →  23  →  14","Raw nominal  →  dual nominal  →  breadth",edge="#E7C2B9",fill="#FCF0EB",maincol="#A44938",size=11.0)
arrow(axA,.690,.601,.610,.466)
arrow(axA,.307,.300,.347,.300);arrow(axA,.644,.300,.683,.300)
axA.text(.02,.018,"Separately: 41 post hoc illustrative features grouped into six functional themes; not a relaxed-P discovery set.",
         transform=axA.transAxes,color=MID,fontsize=7.3,va="bottom")
panel_line(.671,.59,.671,.895)

# B. 14-focal horizontal directional support graphic — not a second heatmap
axB=fig.add_axes([.690,.604,.273,.265]);axB.set_xlim(-7.5,19.9);axB.set_ylim(-.95,15.25)
axB.axis("off")
tag(axB,"B","Lineage recurrence in the 14 focal candidates","Descriptive supporting genera and phyla, not genome-wide significance")
def short_id(k):
 if k.startswith("PFAM:"):return k.replace("PFAM:","")
 return k.replace("RHEA:RHEA:","RHEA:").replace("KMOD:","")
def short_description(k,a):
 lab=str(a)
 if k=="PFAM:PF01638":return "HxlR"
 if k=="PFAM:PF04072":return "LCM"
 if k=="COG:COG1216":return "GT2"
 if k=="KMOD:M00545":return "cinnamate pathway"
 if k.startswith("RHEA"):return "reaction"
 return lab if len(lab)<=12 else lab[:10]+"…"
F=sorted(focal,key=lambda z:(0 if z["direction"]=="plaque" else 1,-int(z["genus_support"])))
assert [z["direction"] for z in F].count("plaque")==8
for i,r in enumerate(F):
 y=14.1-i
 color=PLA if r["direction"]=="plaque" else MUC
 axB.scatter([-7.0],[y],color=color,s=16,zorder=5)
 lab=short_id(r["feature"])+" · "+short_description(r["feature"],r["annotation"])
 axB.text(-6.4,y,lab,color=TEXT,fontsize=7.1,va="center",ha="left")
 axB.barh(y,13,.53,color="#EDF1F2",left=0,edgecolor="none")
 axB.barh(y,int(r["genus_support"]),.53,color=color,left=0,edgecolor="none")
 axB.text(14.0,y,f'{int(r["genus_support"])}/13',color=NAVY,ha="right",va="center",fontsize=7.6)
 axB.text(17.8,y,f'{int(r["phyla_support"])}/6',color=MID,ha="right",va="center",fontsize=7.6)
 if i==7:axB.plot([-7.4,19.0],[y-.57,y-.57],lw=.8,color=LINE)
axB.text(.3,15.03,"SUPPORTING GENERA",fontsize=7.1,color=MID)
axB.text(17.8,15.03,"PHYLA",fontsize=7.1,color=MID,ha="right")
axB.text(-7.2,15.03,"ANNOTATION",fontsize=7.1,color=MID)
axB.text(-7.3,-.75,"Coral: plaque direction  •  Teal: mucosal direction",color=MID,fontsize=7.0)

# C. 41 curated features across 13 genera; effects calculated from all HMTs in each genus
col=[G.index(k) for k in names];Z=np.asarray(D["raw41"],dtype=float)[:,col]
C=fig.add_axes([.328,.107,.397,.438])
C.imshow(Z,cmap=CMAP,norm=NORM,origin="upper",interpolation="none",aspect="auto")
C.set_xlim(-.5,12.5);C.set_ylim(40.5,-.5);C.set_yticks([])
C.set_xticks(range(13));C.set_xticklabels([k if len(k)<13 else k[:11]+'.' for k in names],
 rotation=42,ha="right",rotation_mode="anchor",fontsize=7.65,color=TEXT)
C.xaxis.set_ticks_position('top');C.tick_params(axis='x',which='both',length=0,pad=14)
for s in C.spines.values():s.set_visible(False)
# fine separation of related phyla
for v in [1.5,5.5,7.5,8.5,9.5]:C.axvline(v,color="white",lw=2.1)
for yy in range(1,41):C.axhline(yy-.5,color="white",lw=.38,alpha=.62)
physeq=np.asarray([phylum[g] for g in names])
# phylum bars placed directly above heatmap without obscuring text
for ph in list(dict.fromkeys(physeq)):
 inds=np.where(physeq==ph)[0];i0=int(inds[0]);i1=int(inds[-1])
 C.plot([i0-.46,i1+.46],[-1.9,-1.9],lw=5.8,color=colorsphy[ph],
        solid_capstyle="butt",clip_on=False)
 C.text((i0+i1)/2,-2.35,{"Campylobacterota":"Campyl.","Pseudomonadota":"Pseudom."}.get(ph,ph),
        color=colorsphy[ph],fontsize=7.3,weight="bold",ha="center",va="bottom",clip_on=False)
# 41 feature names in aligned label axis
L=fig.add_axes([.050,.107,.194,.438],sharey=C)
L.set_xlim(0,1);L.set_ylim(40.5,-.5);L.axis("off")
theme_colors={
"Candidate / poorly characterized":"#467FAD","Regulation & sensing":"#34A3A1",
"Envelope / interface":"#689CB6","Metabolism & redox":"#63A28B",
"Mobile elements & defence":"#9A79B8","DNA repair & maintenance":"#C49B43"}
grp=Counter(x["group"] for x in atlas);assert sum(grp.values())==41
row_i=0
for group,k in grp.items():
 first=row_i;last=row_i+k-1;color=theme_colors[group]
 L.add_patch(Rectangle((.000,first-.45),.015,k-.1,fc=color,ec="none"))
 short=group.replace("Candidate / poorly characterized","Candidate / uncharacterized").replace("DNA repair & maintenance","DNA repair")
 L.text(.030,(first+last)/2,short+"\n("+str(k)+")",fontsize=7.05,color=color,ha="left",va="center",weight="bold",linespacing=1.3)
 row_i+=k
# actual ID/description printed over heatmap left edge in dedicated overlay column
I=fig.add_axes([.166,.107,.065,.438],sharey=C)
I.set_xlim(0,1);I.set_ylim(40.5,-.5);I.axis("off")
focids={x["feature"] for x in focal}
for i,item in enumerate(atlas):
 fid=short_id(item["id"])
 I.text(.01,i,fid,fontsize=6.7,color=NAVY,va="center",
        weight="bold" if item["id"] in focids else "normal")
# Function descriptions, directly adjacent to feature identifiers
J=fig.add_axes([.233,.107,.092,.438],sharey=C)
J.set_xlim(0,1);J.set_ylim(40.5,-.5);J.axis("off")
for i,item in enumerate(atlas):
 lab=str(item["short"]).replace(" domain", "").replace(" family", "")
 if len(lab)>20:lab=lab[:19]+"…"
 J.text(.01,i,lab,fontsize=6.35,color=TEXT,va="center")
# thematic boundaries with aligned horizontal rules
cc=0
for grpname,k in grp.items():
 cc+=k
 if cc<41:
  C.axhline(cc-.5,color="white",lw=2.4)
  L.plot([0,.99],[cc-.5,cc-.5],color=LINE,lw=.8)
# C panel heading outside top axis, editorial
fig.text(.043,.625,"C",color=NAVY,fontsize=19,weight="bold")
fig.text(.073,.629,"Genus-specific functional landscape",color=NAVY,fontsize=11.9,weight="bold")
fig.text(.073,.607,"41 curated features across six themes · original plaque-minus-mucosa contrasts",
         color=MID,fontsize=8.5)
fig.text(.050,.557,"THEME",color=MID,size=7.1,weight="bold")
fig.text(.165,.557,"ID",color=MID,size=7.1,weight="bold")
fig.text(.234,.557,"PUTATIVE FUNCTION",color=MID,size=7.1,weight="bold")
# legends and method caveat
fig.text(.047,.077,"Rows are post hoc illustrative annotations; Pfam/COG and derived biochemical-profile scores should not be treated as identical enzyme activities.",
         color=MID,fontsize=7.05)
fig.text(.047,.058,"No individual profile passes all-feature maxT FWER. Differences show equal genus effects, not counts of independent genomes.",
         color=NAVY,fontsize=7.45)
cb_ax=fig.add_axes([.285,.037,.215,.013])
cb=plt.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cb_ax,orientation="horizontal")
cb.ax.tick_params(labelsize=7.5,length=2,pad=2);cb.outline.set_visible(False)
fig.text(.507,.042,"mucosa ←  mean prevalence/score difference  → plaque",color=MID,fontsize=7.55)

# D. Four representative candidate distributions across genera, no heatmap
axD=fig.add_axes([.770,.142,.190,.397])
axD.text(.00,1.044,"D",transform=axD.transAxes,color=NAVY,fontsize=19,weight="bold")
axD.text(.10,1.052,"Four genomic effects across genera",transform=axD.transAxes,color=NAVY,fontsize=9.45,weight="bold")
axD.text(.10,1.014,"One dot per genus; zero denotes no difference",transform=axD.transAxes,color=MID,fontsize=7.45)
axD.set_xlim(-1.05,1.05);axD.set_ylim(-.55,4.25)
axD.spines[['top','right','left']].set_visible(False)
axD.spines['bottom'].set_color("#CAD7D9")
axD.set_xticks([-1,-.5,0,.5,1]);axD.tick_params(axis="x",labelsize=8,colors=MID)
axD.set_yticks([]);axD.axvline(0,color="#9DB2BA",lw=.9,ls="--",zorder=0)
key=["PFAM:PF01638","PFAM:PF04072","PFAM:PF13577","COG:COG1216"]
fidx={r["feature"]:i for i,r in enumerate(focal)}
genus_cols=np.asarray([colorsphy[phylum[g]] for g in names])
rng=np.random.default_rng(82)
for i,f in enumerate(key):
 y=3.65-i*1.06
 rr=np.asarray(D["raw14"][fidx[f]],dtype=float)[col]
 # deterministic, barely perturbed points for overlapping values
 jitter=np.random.default_rng(i+12).uniform(-.12,.12,13)
 axD.hlines(y,-1,1,color="#E2E9E9",lw=.75,zorder=0)
 axD.scatter(rr,y+jitter,s=31,c=genus_cols,ec="white",lw=.60,zorder=3)
 mu=rr.mean()
 axD.scatter([mu],[y],marker="D",s=66,color=PLA if mu>0 else MUC,ec=NAVY,lw=.7,zorder=5)
 nm=next(x for x in focal if x["feature"]==f)
 title={"PFAM:PF01638":"HxlR","PFAM:PF04072":"LCM","PFAM:PF13577":"SnoaL_4","COG:COG1216":"GT2 glycosyltransferase"}[f]
 axD.text(-1.01,y+.33,title,ha="left",va="bottom",fontsize=8.85,color=NAVY,weight="bold")
 axD.text(1.02,y+.33,f'{int(nm["genus_support"])} genera · {int(nm["phyla_support"])} phyla',
     ha="right",va="bottom",fontsize=7.15,color=MID)
axD.set_xlabel("Plaque – mucosal difference within genus",fontsize=8.0,color=MID,labelpad=6)
fig.text(.771,.075,"●  Genus effect (coloured by phylum)",color=MID,fontsize=8)
fig.text(.771,.059,"◆  Equal-genus mean",color=MID,fontsize=8)
panel_line(.740,.08,.740,.570)

caption="""Fig. 2 | Exploratory functional-genomic recurrence across 13 oral bacterial genera.
(a) Six eHOMD abundance summaries were used to assign plaque versus mucosal ecological preference to 47 HMTs. Genomic features were annotated in 1,121 reference assemblies and aggregated as HMT-level prevalence profiles. Each genus contributed one plaque-minus-mucosal effect per genomic profile. The 13,917 deduplicated profiles were assessed using 4,096 within-genus ecological-label permutations, before and after residualizing annotation prevalence against genus, genome size, genomic GC content and predicted CDS count. Nominal consistency P<0.05 occurred for 119 features in the raw scan; 23 had nominal support in both raw and genome-adjusted scans, and 14 met an additional support-breadth screen (at least five genera and three phyla) with the same directional interpretation. The 41-feature thematic atlas was separately selected post hoc to illustrate biological diversity rather than by an unadjusted statistical cutoff. (b) Descriptive support breadth across 13 genera and six phyla for the 14 exploratory shortlist annotations. The shown counts quantify directional genomic-feature differences, not per-feature genome-wide significance. (c) Signed within-genus plaque-minus-mucosa differences for 41 post hoc features across 13 genera, grouped into six annotation themes. The values are differences in reference-genome feature prevalence for Pfam profiles and differences in annotation-derived module/pathway scores for other layers; they do not demonstrate enzyme activity. The atlas was not selected under a prespecified statistical threshold. (d) Distributions of raw genus-specific effects for four illustrative annotations; each dot corresponds to one genus and the diamond denotes the equal-genus mean. None of the 13,917 individual profiles survived the all-feature maxT family-wise error correction. Genomes were nested within HMTs rather than treated as independent ecological replicates. Functional interpretations remain putative."""
(P/"Figure2_caption.txt").write_text(caption+"\n")
# aligned source data exports
columns=["feature","annotation","theme"]+names
data_rows=[]
for i,item in enumerate(atlas):
 data_rows.append([item["id"],item["short"],item["group"]]+[float(z) for z in Z[i]])
pd.DataFrame(data_rows,columns=columns).to_csv(P/"Figure2_C_41atlas_genus_effects.tsv",sep="\t",index=False)
fq=pd.DataFrame(focal);fq.to_csv(P/"Figure2_B_14focal_statistics.tsv",sep="\t",index=False)
counts=pd.DataFrame([dict(stage="genomic profiles",n=13917),dict(stage="raw consistency p<0.05",n=119),
 dict(stage="raw and genome-adjusted nominal",n=23),dict(stage="nominal plus 5 genera 3 phyla",n=14),
 dict(stage="curated illustrative atlas, independent branch",n=41),
 dict(stage="all-feature FWER at p<0.05",n=0)])
counts.to_csv(P/"Figure2_A_pipeline_counts.tsv",sep="\t",index=False)
for ext in ["pdf","svg","png"]:
 fn=P/("HOMD_Figure2_functional_genomic_atlas."+ext)
 fig.savefig(fn,facecolor="white",dpi=205 if ext=="png" else 175)
 print("CREATED",fn,fn.stat().st_size,flush=True)
plt.close(fig)
