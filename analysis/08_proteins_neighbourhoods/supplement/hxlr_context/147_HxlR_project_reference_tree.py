#!/usr/bin/env python3
"""Project all 1,524 proteins to existing MMseqs-X5 family representatives on
unrooted 73+7 anchor ML phylogeny; bootstrap labels as reported by RAxML-NG.
No experimental ligand assignments transferred from closest branches.
"""
from pathlib import Path
import pandas as pd,numpy as np,json
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('__PROJECT_ROOT__/study13')
O=B/'16_HxlR_reference_phylogeny_20261014'
S=Path('__WORK_ROOT__/HOMD_HxlR_reference_20261014')
files=[O/'HxlR_refX5.raxml.support',S/'HxlR_refX5.raxml.support',
       O/'HxlR_refX5.raxml.bestTree',S/'HxlR_refX5.raxml.bestTree']
phylo_file=next((x for x in files if x.exists()),None)
if phylo_file is None:raise FileNotFoundError('RAxML-NG final reference tree not yet available')
tree=Phylo.read(phylo_file,'newick')
refs=sorted(t.name for t in tree.get_terminals() if t.name.startswith('REF__'))
h=pd.read_csv(O/'HxlR_1524_oral_closest_reference_cautious.tsv',sep='\t')
rep=pd.read_csv(O/'oral_73_X5_family_representative_selection.tsv',sep='\t')
m={r.MMseqs_X5_representative:r.protein_id for r in rep.itertuples()}
tip={t.name for t in tree.get_terminals()}
assert len(refs)==7 and len(tip)==80,(len(refs),len(tip))
assert set(m.values()).issubset(tip)
rec=[]
for r in rep.itertuples():
 d={ref:tree.distance(r.protein_id,ref) for ref in refs}
 ordered=sorted(d.items(),key=lambda kv:kv[1])
 rec.append({'MMseqs_X5_representative':r.MMseqs_X5_representative,
             'phylo_tip':r.protein_id,'nearest_reference_on_tree':ordered[0][0].replace('REF__',''),
             'branch_distance_to_reference':round(ordered[0][1],5),
             'second_nearest_reference':ordered[1][0].replace('REF__',''),
             'ratio_nearest_to_next_branch_distance':round(ordered[0][1]/max(ordered[1][1],1e-9),3),
             'branch_distance_second_reference':round(ordered[1][1],5)})
d=pd.DataFrame(rec)
d.to_csv(O/'HxlR_73_X5_cluster_phylogenetic_reference_map.tsv',sep='\t',index=False)
mapped=h.merge(d,on='MMseqs_X5_representative',how='left',validate='many_to_one')
mapped.to_csv(O/'HxlR_1524_projection_onto_reference_X5_tree.tsv',sep='\t',index=False)
group=mapped.groupby('nearest_reference_on_tree').agg(n_proteins=('protein_id','size'),n_genomes=('GCA','nunique'),n_HMT=('HMT','nunique'),n_genera=('Genus','nunique'),n_X5_cluster=('MMseqs_X5_representative','nunique')).reset_index()
group.to_csv(O/'HxlR_reference_tree_nearest_branch_count.tsv',sep='\t',index=False)
print('PHYLO',phylo_file,'tip',len(tip),'refs',len(refs),flush=True)
print('GROUPS\n'+group.to_string(index=False),flush=True)
genera=mapped.groupby(['MMseqs_X5_representative','Genus']).size().reset_index(name='n').sort_values('n',ascending=False).drop_duplicates('MMseqs_X5_representative')
gen_by_rep=mapped.groupby('MMseqs_X5_representative').Genus.nunique().to_dict()
labels={}
for r in rep.itertuples():
 sub=mapped[mapped.MMseqs_X5_representative==r.MMseqs_X5_representative]
 majority=sub.Genus.value_counts().idxmax()
 mult=' +' if sub.Genus.nunique()>1 else ''
 labels[r.protein_id]=f'{majority}{mult}  {r.MMseqs_X5_representative}  [n={len(sub)}]'
for ref in refs:labels[ref]=ref.replace('REF__','REFERENCE: ')
def label(c):return labels.get(c.name,'') if c.is_terminal() else None
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(14,20))
fig.patch.set_facecolor('white')
Phylo.draw(tree,axes=ax,label_func=label,do_show=False,show_confidence=False)
for item in ax.texts:
 if 'REFERENCE:' in item.get_text():
  item.set_color('#a94931');item.set_fontweight('bold');item.set_fontsize(8.3)
 else:
  item.set_fontsize(6.65);item.set_color('#263d4a')
ax.set_title('HxlR/DUF24 family: 73 oral sequence-cluster representatives + 7 functional references',loc='left',fontsize=13,fontweight='bold',pad=14)
ax.set_xlabel('Estimated amino-acid substitutions per site (LG+G)',fontsize=9)
ax.set_ylabel('')
fig.text(.12,.015,'Exploratory unrooted phylogeny displayed with arbitrary orientation; low-support clades and short domain sequences limit functional assignment. + marks clusters comprising several genera.',fontsize=8,color='#586a75')
fig.savefig(O/'HxlR_73cluster_reference_tree.pdf',bbox_inches='tight')
fig.savefig(O/'HxlR_73cluster_reference_tree.svg',bbox_inches='tight')
fig.savefig(O/'HxlR_73cluster_reference_tree.png',dpi=170,bbox_inches='tight')
plt.close(fig)
print('FIGURE_CREATED',flush=True)
