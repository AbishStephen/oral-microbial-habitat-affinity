from pathlib import Path
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
R=Path(__file__).resolve().parent
(R/'figures').mkdir(exist_ok=True)
D=json.loads((R/'Figure5_model_predictions_and_sister_pairs.json').read_text())
P=pd.DataFrame(D['pred']);A=pd.DataFrame(D['auc'])
M=pd.read_csv(R/'model_performance_and_uncertainty.tsv',sep='\t')
P.to_csv(R/'Figure5_HMT_predictions.tsv',sep='\t',index=False)
A.to_csv(R/'Figure5_genus_level_auc.tsv',sep='\t',index=False)
NAV='#15384D';TEXT='#284457';GRAY='#68808B';GRID='#D9E3E6';PCLA='#D7785B';MUCA='#278D92';BGC='#F8FBFB'
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'font.size':8.5})
fig=plt.figure(figsize=(17,11.5),dpi=175,facecolor='white')
fig.text(.04,.965,'CROSS-LINEAGE TRANSFER',color=MUCA,fontsize=10,weight='bold')
fig.text(.04,.938,'Genomic habitat prediction generalizes imperfectly to unseen genera',color=NAV,fontsize=18,weight='bold')
fig.text(.04,.914,'Nested genus-held-out evaluation with source-weighting, compositional adjustment and full-refitting permutation nulls',color=GRAY,fontsize=9)
fig.add_artist(Line2D([.04,.963],[.9,.9],color=GRID,lw=.9,transform=fig.transFigure))

def title(ax,letter,label,sub=''):
    ax.text(0,1.07,letter,fontweight='bold',fontsize=17,color=NAV,va='bottom',transform=ax.transAxes)
    ax.text(.061,1.09,label,fontweight='bold',fontsize=10.5,color=NAV,transform=ax.transAxes,va='bottom')
    if sub: ax.text(.061,1.045,sub,color=GRAY,fontsize=7.9,transform=ax.transAxes,va='bottom')

# A: 1/4 wide small design diagram: no overclaim on 13,917 feature selected
ax=fig.add_axes([.047,.663,.906,.195]);ax.axis('off');title(ax,'A','Genus-held-out prediction design','Every held-out genus is excluded from feature ranking, hyperparameter selection and model fitting')
ax.set_xlim(0,1);ax.set_ylim(0,1)
steps=[(.02,.29,'47 oral HMTs','13 genera / 13,917 profiles'),(.275,.29,'Outer leave-one-genus-out','12 training genera → 1 unseen genus'),(.533,.29,'Inner genus-based tuning','Feature ranking; 5–100 features; ridge'),(.791,.29,'Held-out evaluation','Genus-level pairwise AUC')]
for i,(x,y,l,s) in enumerate(steps):
    w=.19
    fc=BGC if i<3 else '#EEF6F5'
    edge=GRID if i<3 else '#AED2CE'
    ax.add_patch(FancyBboxPatch((x,y),w,.47,boxstyle='round,pad=.016,rounding_size=.023',fc=fc,ec=edge,lw=.95))
    ax.text(x+w/2,y+.30,l,color=NAV,fontsize=8.75,weight='bold',ha='center',va='center')
    ax.text(x+w/2,y+.13,s,color=GRAY,fontsize=7.2,ha='center',va='center')
    if i<3: ax.add_patch(FancyArrowPatch((x+w+.008,y+.236),(steps[i+1][0]-.012,y+.236),arrowstyle='-|>',mutation_scale=14,lw=1.25,color='#8EADB6'))
ax.text(.02,.13,'Feature selection repeated within training folds · 200 full-refitting within-genus permutations per outcome · adjusted models: 100 permutations',fontsize=8,color=GRAY)

# B: Prediction per HMT (major panel): x gen labels across 13 matched genera, y signed predicted affinity
order=['Fusobacterium','Leptotrichia','Streptococcus','Gemella','Veillonella','Selenomonas','Actinomyces','Rothia','Campylobacter','Neisseria','Porphyromonas','Prevotella','Segatella']
assert len(P)==47 and set(order)==set(A.Genus)
ax=fig.add_axes([.055,.335,.90,.257]);title(ax,'B','Predictions for all 47 held-out taxa','Source-family-balanced outcome; points are individual HMTs, not genomes')
ax.axhline(0,color=NAV,lw=.8,alpha=.6,ls='--')
for i,g in enumerate(order):
    z=P[P.Genus==g]
    ax.axvspan(i-.45,i+.45,color='#F5F8F8' if i%2==0 else 'white',zorder=-3)
    for j,eco in enumerate(['mucosal','plaque']):
        zz=z[z.Ecology==eco].sort_values('HMT')
        # deterministic x offsets avoiding visual overlap
        offs=np.linspace(-.19,.19,len(zz)) if len(zz)>1 else np.array([0.])
        ax.scatter(i+offs,zz.Y_pred,s=30,c=MUCA if eco=='mucosal' else PCLA,
                   edgecolor='white',linewidth=.65,zorder=5,alpha=.95)
    if len(z)>2:
        lo=z[z.Ecology=='mucosal'].Y_pred.mean();hi=z[z.Ecology=='plaque'].Y_pred.mean()
        ax.plot([i,i],[lo,hi],color='#A1B5BA',lw=2,zorder=3)
ax.set_xlim(-.55,len(order)-.45);ax.set_xticks(range(len(order)))
ax.set_xticklabels(order,rotation=28,ha='right',fontsize=8.1,color=TEXT)
ax.set_ylabel('Out-of-genus predicted ecological affinity (a.u.)',color=TEXT)
ax.tick_params(axis='y',colors=GRAY,labelsize=7.8)
ax.spines[['top','right']].set_visible(False)
ax.spines[['bottom','left']].set_color(GRID)
ax.legend([Line2D([0],[0],marker='o',linestyle='none',color=PCLA,markersize=6),Line2D([0],[0],marker='o',linestyle='none',color=MUCA,markersize=6)],['Plaque-associated HMT','Mucosa-associated HMT'],loc='upper right',frameon=False,fontsize=8,ncol=2)
# AUC per genus below: no ax reference crossing
axC=fig.add_axes([.055,.087,.565,.154]);title(axC,'C','Transferability is heterogeneous across genera','One macro-AUC component per genus; chance = 0.5')
Z=A.set_index('Genus').loc[order]
for i,g in enumerate(order):
    y=12-i
    x0=Z.loc[g,'sourcebalanced_genusbalanced']
    x1=Z.loc[g,'sourcebalanced_genusbalanced_covadj']
    axC.hlines(y,min(x0,x1),max(x0,x1),color='#B1C7CC',lw=2.1)
    axC.scatter([x0],[y],c=PCLA,s=38,zorder=4)
    axC.scatter([x1],[y],c=MUCA,s=38,zorder=4,marker='D')
axC.axvline(.5,color=GRAY,lw=.8,ls='--')
axC.set_ylim(-.7,12.7);axC.set_xlim(-.1,1.06)
axC.set_yticks(np.arange(13));axC.set_yticklabels(order[::-1],fontsize=6.8,color=TEXT)
axC.set_xticks([0,.25,.5,.75,1]);axC.tick_params(axis='x',labelsize=8,colors=GRAY)
axC.spines[['top','right']].set_visible(False)
axC.spines[['left','bottom']].set_color(GRID)
axC.set_xlabel('Within-genus paired AUC',fontsize=8,color=TEXT)
axC.legend([Line2D([0],[0],marker='o',linestyle='none',color=PCLA,markersize=6),Line2D([0],[0],marker='D',linestyle='none',color=MUCA,markersize=6)],['Unadjusted','Genome-covariate residualized'],loc='lower right',frameon=False,fontsize=7.1,ncol=2)

# D: benchmark point estimates with 95% genus bootstrap CIs, P values from true perm nulls
axD=fig.add_axes([.704,.087,.25,.154]);title(axD,'D','Prediction and covariate sensitivity','Macro AUC and genus-bootstrap 95% intervals')
scores=M[M.model.str.startswith('Genus-balanced')]
variants=[('original','Genus-balanced nested'),('original','Genus-balanced nested, genome covariate-residualized'),('source-family-balanced','Genus-balanced nested'),('source-family-balanced','Genus-balanced nested, genome covariate-residualized')]
labels=['Original','Original, adjusted','Source balanced','Source balanced, adjusted']
for i,((eco,model),label) in enumerate(zip(variants,labels)):
    r=scores[(scores.ecology_outcome==eco)&(scores.model==model)]
    assert len(r)==1, (eco,model)
    r=r.iloc[0];y=3-i
    c=MUCA if 'adjusted' in label else PCLA
    axD.hlines(y,float(r.genus_bootstrap_95_low),float(r.genus_bootstrap_95_high),color=c,lw=2.5,zorder=2)
    axD.scatter(float(r.macro_within_genus_AUC),y,c=c,s=46,zorder=4,edgecolor='white',linewidth=.6)
    p=r.permutation_empirical_p
    axD.text(1.035,y,f'P={p:.3f}',va='center',fontsize=7,color=TEXT)
axD.axvline(.5,color=GRAY,lw=.8,ls='--')
axD.set_xlim(.34,1.20);axD.set_ylim(-.6,3.75)
axD.set_yticks([3,2,1,0]);axD.set_yticklabels(labels,fontsize=7.4,color=TEXT)
axD.set_xticks([.5,.7,.9,1]);axD.tick_params(axis='x',colors=GRAY,labelsize=7)
axD.spines[['top','right','left']].set_visible(False);axD.spines['bottom'].set_color(GRID)
axD.set_xlabel('Macro within-genus AUC',fontsize=7.5,color=TEXT)
fig.text(.045,.032,'Raw modelling: P=0.00995 (200 full refits per outcome). Covariate-adjusted: P=0.099 original; P=0.059 source balanced (100 full refits).',fontsize=8,color=NAV)
for ext in ['png','pdf','svg']:
    path=R/'figures'/f'HOMD_Figure5_prediction_only.{ext}'
    fig.savefig(path,dpi=240 if ext=='png' else 175,bbox_inches='tight',facecolor='white')
plt.close(fig)
print('Created Figure 5 prediction only; files in',R/'figures')
