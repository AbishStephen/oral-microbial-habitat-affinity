#!/usr/bin/env python3
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr,binomtest
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
B=Path('__PROJECT_ROOT__')
E=B/'study13'
O=E/'08_validation'/'confound_audit';O.mkdir(exist_ok=True,parents=True)
d=pd.read_csv(E/'08_validation/SPARSE_DEDUP_outer_predictions.tsv',sep='\t')
f=pd.read_csv(E/'06_hmt_matrices/X2_Pfam_HMT_prevalence_TRANSFERABLE.tsv',sep='\t')
master=pd.read_csv(B/'00_metadata/HOMD_genomic_metadata.tsv',sep='\t',low_memory=False)
manifest=pd.read_csv(E/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')
covcols=['GENOME_ID','combined_size','GC','prokka_CDS','contigs','isolation_source','hmt_primary_body_site_w_abundance','MAG']
master=master.drop_duplicates(subset=["GENOME_ID"],keep="first")
m=manifest.merge(master[covcols],left_on='GCA',right_on='GENOME_ID',how='left',validate='one_to_one')
print('GENOMES',len(m),'metadata missing',m.combined_size.isna().sum())
m[covcols[1:5]] = m[covcols[1:5]].apply(pd.to_numeric,errors='coerce')
agg=m.groupby('HMT').agg(GenomeSize_Mb=('combined_size',lambda z:np.nanmedian(z)/1e6),
 GC_percent=('GC','median'), CDS=('prokka_CDS','median'),Contigs=('contigs','median'),
 n_genomes=('GCA','size'))
d=d.merge(agg,left_on='HMT',right_index=True,validate='one_to_one')
d=d.merge(f[['HMT','PF04072','PF12680','PF10604','PF08020','PF14534','PF13199','PF22815']],on='HMT',validate='one_to_one')
for key in ['Y','Y_pred','GenomeSize_Mb','GC_percent','CDS','Contigs','PF04072','PF12680','PF10604']:
 d[key+'_cen']=d[key]-d.groupby('Genus')[key].transform('mean')
def pair_auc(z,pred):
 x=z.loc[z.Ecology=='plaque',pred].values
 y=z.loc[z.Ecology=='mucosal',pred].values
 return float(np.mean([(1 if a>b else .5 if a==b else 0) for a in x for b in y]))
genus=[]
for g,z in d.groupby('Genus'):
 rho=spearmanr(z.Y,z.Y_pred).statistic if len(z)>2 else float('nan')
 genus.append([g,len(z),sum(z.Ecology=='plaque'),sum(z.Ecology=='mucosal'),rho,pair_auc(z,'Y_pred'),float(z.loc[z.Ecology=='plaque','Y_pred'].mean()-z.loc[z.Ecology=='mucosal','Y_pred'].mean())])
pg=pd.DataFrame(genus,columns=['Genus','n_HMTs','n_plaque','n_mucosal','rho','AUC_within','pred_plaque_minus_mucosa'])
pg.to_csv(O/'per_genus_prediction.tsv',sep='\t',index=False)
corrs=[]
for key in ['Y_pred','GenomeSize_Mb','GC_percent','CDS','Contigs','PF04072','PF12680','PF10604']:
 x=d[key+'_cen'];y=d.Y_cen
 corrs.append([key,float(pearsonr(x,y).statistic),float(spearmanr(x,y).statistic),float(spearmanr(d[key],d.Y).statistic)])
pd.DataFrame(corrs,columns=['variable','genus_resid_pearson','genus_resid_spearman','pooled_spearman']).to_csv(O/'genus_adjusted_correlations.tsv',sep='\t',index=False)
# Construct clean and comparable size-only nested LOGO baseline. All steps and tuning inside outer-fold.
from sklearn.pipeline import make_pipeline
G=d.Genus.to_numpy();Y=d.Y.to_numpy()
X=d[['GenomeSize_Mb','GC_percent','CDS']].values.astype(float)
# Guard potential NA, impute training-fold medians only
pred=np.zeros(len(d));fitpars=[]
for outer in np.unique(G):
 test=G==outer;tr=~test;XX=X[tr];YY=Y[tr];GG=G[tr]
 best=(-1e99,None)
 for aa in [.01,.1,1,10,100,1000]:
  pr=np.zeros(len(XX))
  for inner in np.unique(GG):
   ti=GG==inner;train=~ti
   fill=np.nanmedian(XX[train],axis=0)
   inp=np.where(np.isfinite(XX[train]),XX[train],fill)
   out=np.where(np.isfinite(XX[ti]),XX[ti],fill)
   mdl=make_pipeline(StandardScaler(),Ridge(alpha=aa))
   mdl.fit(inp,YY[train]);pr[ti]=mdl.predict(out)
  obj=spearmanr(YY,pr).statistic
  if not np.isfinite(obj):obj=-9
  if obj>best[0]:best=(obj,aa)
 aa=best[1]
 fill=np.nanmedian(XX,axis=0)
 mdl=make_pipeline(StandardScaler(),Ridge(alpha=aa))
 mdl.fit(np.where(np.isfinite(XX),XX,fill),YY)
 pred[test]=mdl.predict(np.where(np.isfinite(X[test]),X[test],fill))
 fitpars.append([outer,aa])
d['size_predict']=pred
d['size_predict_cen']=d.size_predict-d.groupby('Genus').size_predict.transform('mean')
d['Ypred_resid_cen']=d.Y_pred_cen
# Standardized within-genus pooled Pearson, Spearman, paired AUC, macro AUC
metrics=[]
for pname in ['Y_pred','size_predict']:
 p=d[pname].to_numpy()
 pooledauc=roc_auc_score((Y>0).astype(int),p)
 ckey=pname+'_cen'
 pear=pearsonr(d.Y_cen,d[ckey]).statistic
 spear=spearmanr(d.Y_cen,d[ckey]).statistic
 pair=np.array([pair_auc(z,pname) for _,z in d.groupby('Genus')])
 pair_weighted=np.average(pair,weights=[sum(z.Ecology=='plaque')*sum(z.Ecology=='mucosal') for _,z in d.groupby('Genus')])
 metrics.append([pname,float(spearmanr(Y,p).statistic),pooledauc,float(pear),float(spear),float(np.mean(pair)),float(pair_weighted),int((pair>0.5).sum())])
pd.DataFrame(metrics,columns=['model','pooled_rho','pooled_AUC','genus_centered_pearson','genus_centered_spearman','macro_within_auc','paired_within_auc','genera_with_auc_above_half']).to_csv(O/'within_genus_vs_genome_size_baseline.tsv',sep='\t',index=False)
# Domain co-associations, raw and within genus centered, to distinguish 'one profile' from three.
F=['PF04072','PF12680','PF10604','PF08020','PF14534','PF13199','PF22815']
co=[]
for i,x in enumerate(F):
 for y in F[i+1:]:
  xx=d[x]-d.groupby('Genus')[x].transform('mean')
  yy=d[y]-d.groupby('Genus')[y].transform('mean')
  co.append([x,y,float(spearmanr(d[x],d[y]).statistic),float(pearsonr(xx,yy).statistic),int(((d[x]>0)&(d[y]>0)).sum())])
pd.DataFrame(co,columns=['Pfam_1','Pfam_2','pooled_spearman','genus_centered_pearson','HMT_cooccurrences']).to_csv(O/'domain_profile_correlations.tsv',sep='\t',index=False)
d.to_csv(O/'hmt_covariates_and_predictions.tsv',sep='\t',index=False)
m[['GCA','HMT','Genus','Ecology']+covcols].to_csv(O/'genome_metadata_provenance.tsv',sep='\t',index=False)
print('GENUS_PREDICTION:\n',pg.to_string(index=False))
print('MODEL_COMPARE:\n',pd.read_csv(O/'within_genus_vs_genome_size_baseline.tsv',sep='\t').to_string(index=False))
print('DOMAIN_CORR:\n',pd.read_csv(O/'domain_profile_correlations.tsv',sep='\t').query("Pfam_1 in ['PF04072','PF12680','PF10604'] and Pfam_2 in ['PF04072','PF12680','PF10604']").to_string(index=False))
print('COVAR_ASSOC:\n',pd.read_csv(O/'genus_adjusted_correlations.tsv',sep='\t').to_string(index=False))
