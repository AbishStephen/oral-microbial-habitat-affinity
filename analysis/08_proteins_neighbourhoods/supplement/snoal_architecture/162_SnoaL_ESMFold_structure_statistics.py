#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from mpl_toolkits.mplot3d import Axes3D
B=Path('__PROJECT_ROOT__/study13')
O=B/'18_SnoaL2_function_gene_integrity_20261016'
models=[
('GCA_001956355.1_01713','Actinomyces | 326 aa','Rossmann-like (9–95)',(9,95),'SnoaL_2 (215–312)',(215,312)),
('GCA_000191105.1_00970','Streptococcus | 254 aa','SnoaL-like I (8–115)',(8,115),'SnoaL-like II (143–240)',(143,240)),
('GCA0028848951_00983','Fusobacterium | 256 aa','SnoaL-like I (11–105)',(11,105),'SnoaL-like II (142–243)',(142,243)),
('GCA_000176755.2_01530','Neisseria | 226 aa','N flank (1–57)',(1,57),'Central SnoaL-like (58–178)',(58,178)),
('GCA_900573985.1_01698','Campylobacter | 209 aa','N flank (1–37)',(1,37),'Central SnoaL-like (38–151)',(38,151)),
]
parser=PDBParser(QUIET=True)
rows=[]
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
fig=plt.figure(figsize=(15,12),facecolor='white')
for i,(name,tit,label1,ab,label2,cd) in enumerate(models):
 path=O/f'{name}_ESMFold_v1.pdb'
 r=list(parser.get_structure(name,str(path)).get_residues())
 res=[x for x in r if 'CA' in x]
 idx=np.array([x.id[1] for x in res])
 xyz=np.array([x['CA'].coord for x in res])
 conf=np.array([x['CA'].bfactor for x in res])
 scale=100 if np.nanmedian(conf)<2 else 1
 conf=conf*scale
 first=(idx>=ab[0])&(idx<=ab[1])
 second=(idx>=cd[0])&(idx<=cd[1])
 x,y=xyz[first],xyz[second]
 dist=np.sqrt(((x[:,None,:]-y[None,:,:])**2).sum(axis=2))
 contacts8=int(np.sum(dist<8));contacts12=int(np.sum(dist<12))
 pair8=int(np.sum(np.min(dist,axis=1)<8)) 
 min_dist=float(dist.min())
 centers=float(np.linalg.norm(x.mean(axis=0)-y.mean(axis=0)))
 t={'protein_id':name,'organism_label':tit,'length':len(r),
 'region1':label1,'region2':label2,
 'mean_pLDDT_full':float(conf.mean()),'mean_pLDDT_region1':float(conf[first].mean()),
 'mean_pLDDT_region2':float(conf[second].mean()),
 'interregion_min_CA_A':min_dist,
 'interregion_CA_pairs_lt8A':contacts8,
 'interregion_CA_pairs_lt12A':contacts12,
 'region1_residues_contacting_region2_lt8A':pair8,
 'interregion_centroid_distance_A':centers,
 'model_method':'ESMFold v1 (single-sequence prediction)',
 'no_PAE_note':'relative region orientations cannot be verified from pLDDT alone'}
 rows.append(t)
 # show trace with domain coloring in five panels
 ax=fig.add_subplot(2,3,i+1,projection='3d')
 masks=[(idx<ab[0],'#8897A8'),(first,'#B87651'),
        ((idx>ab[1])&(idx<cd[0]),'#9AA5B1'),
        (second,'#2A938A'),(idx>cd[1],'#8897A8')]
 for mask,color in masks:
  pts=xyz[mask]
  if len(pts)>1: ax.plot(pts[:,0],pts[:,1],pts[:,2],linewidth=2.0,color=color)
 # segment joints in neutral color
 ax.plot(xyz[:,0],xyz[:,1],xyz[:,2],linewidth=.35,color='#9FA6AE',alpha=.45)
 ax.set_title(f'{tit}\npLDDT {conf.mean():.0f} | contacts <8Å: {contacts8}',fontsize=10,fontweight='medium')
 ax.set_axis_off()
 mins=xyz.min(axis=0);maxs=xyz.max(axis=0);cent=(mins+maxs)/2;span=float(np.max(maxs-mins)*.6)
 ax.set(xlim=(cent[0]-span,cent[0]+span),ylim=(cent[1]-span,cent[1]+span),zlim=(cent[2]-span,cent[2]+span))
 ax.view_init(20,-57)
df=pd.DataFrame(rows)
df.to_csv(O/'SnoaL2_five_ESMFold_models_structural_contacts.tsv',sep='\t',index=False)
fig.text(.06,.05,'Orange: region 1; teal: region 2; grey: intervening or terminal sequence. Predicted Cα backbone traces; not experimental structures.',fontsize=11,color='#4C5B66')
fig.text(.06,.025,'High pLDDT supports local residue geometry but not relative domain orientation. Inter-region contacts are model-derived.',fontsize=10,color='#4C5B66')
fig.suptitle('SnoaL-like long proteins: five predicted structural representatives',fontsize=17,y=.96,fontweight='bold')
fig.savefig(O/'SnoaL2_five_long_protein_predicted_folds.png',dpi=220,bbox_inches='tight')
fig.savefig(O/'SnoaL2_five_long_protein_predicted_folds.pdf',bbox_inches='tight')
plt.close(fig)
# N/C CA inter-region contact matrix for each representative
fig,axes=plt.subplots(2,3,figsize=(13,9),facecolor='white')
for i,(name,tit,l1,ab,l2,cd) in enumerate(models):
 r=list(parser.get_structure(name,str(O/f'{name}_ESMFold_v1.pdb')).get_residues())
 xx=np.array([x['CA'].coord for x in r if 'CA' in x]);ii=np.array([x.id[1] for x in r if 'CA' in x])
 xx1=xx[(ii>=ab[0])&(ii<=ab[1])];xx2=xx[(ii>=cd[0])&(ii<=cd[1])]
 dd=np.sqrt(((xx1[:,None,:]-xx2[None,:,:])**2).sum(axis=2))
 ax=axes.flat[i];ax.imshow(dd.T,cmap='viridis_r',origin='lower',aspect='auto',vmin=0,vmax=35);ax.set_title(tit,fontsize=10)
 ax.set_xlabel(l1,fontsize=8);ax.set_ylabel(l2,fontsize=8)
axes.flat[-1].axis('off')
fig.suptitle('Between-region Cα distance maps (Å): darker means closer',fontsize=15)
fig.tight_layout()
fig.savefig(O/'SnoaL2_predicted_region_distance_maps.pdf',bbox_inches='tight')
fig.savefig(O/'SnoaL2_predicted_region_distance_maps.png',dpi=160,bbox_inches='tight')
plt.close(fig)
print(df[['organism_label','length','mean_pLDDT_full','mean_pLDDT_region1','mean_pLDDT_region2','interregion_min_CA_A','interregion_CA_pairs_lt8A','interregion_centroid_distance_A']].round(2).to_string(index=False),flush=True)
