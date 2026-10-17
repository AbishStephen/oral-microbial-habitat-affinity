#!/usr/bin/env python3
"""Build version-aware HOMD genome representative panel for HxlR/PF12680 challenge test.
HOMD public metadata endpoint missing a column header for abundance category;
parse as 41 fields, insert hmt_abundance after primary site.
Freeze source dated 2026-10-14. FTP protein annotations indexed from 2025-07-28.
"""
from pathlib import Path
import csv,collections,json
import pandas as pd
B=Path('__PROJECT_ROOT__/study13')
O=B/'14_HxlR_fullHOMD_controls_20261014'
f=O/'HOMD_genome_metadata_all.tsv'
with open(f) as h:
 rd=csv.reader(h,delimiter='\t');head=next(rd);head.insert(6,'hmt_abundance_category')
 vals=[]
 for row in rd:
  if len(row)!=len(head):continue
  vals.append(row)
m=pd.DataFrame(vals,columns=head)
ftp=pd.read_csv(O/'HOMD_PROKKA_GCA_ID_info_2025.csv')
print('metadata rows',len(m),'FTP rows',len(ftp),'source HMT',m.HMT_ID.nunique(),flush=True)
print('site label counts',m.hmt_primary_body_site_w_abundance.value_counts().to_dict(),flush=True)
print('selected FTP genomic overlap',len(set(m.GENOME_ID)&set(ftp['GCA-ID'])),flush=True)
for col in ['contigs','combined_size','checkM2_completeness','checkM2_contamination','GC']:
 if col in m.columns:m[col]=pd.to_numeric(m[col],errors='coerce')
ftp['GCA']=ftp['GCA-ID']
x=ftp.merge(m.drop_duplicates('GENOME_ID'),left_on='GCA',right_on='GENOME_ID',how='left')
print('matched',x.GENOME_ID.notna().sum(),'unmatched',x.GENOME_ID.isna().sum(),flush=True)
# Rank genome quality by assembly and QC, prioritizing non-MAG, completeness, contiguity, independent taxa.
x['non_MAG']=x['MAG'].isna() | (x['MAG'].str.lower()!='yes')
x['quality_score']=x.checkM2_completeness.fillna(90) - 5*x.checkM2_contamination.fillna(0)
x['contigs_rank']=x['contigs'].fillna(x.Contigs)
x['assembly_complete']=x['assembly_level'].astype(str).str.lower().isin(['complete genome','chromosome'])
x=x.sort_values(['non_MAG','assembly_complete','quality_score','contigs_rank','combined_size'],ascending=[False,False,False,True,False])
x['HMT_link']=x['HMT-ID'].fillna(x.HMT_ID)
rep=x.drop_duplicates('HMT_link').copy()
rep['main_genus']=rep.Genus
rep['main_species']=rep.Species
cols=['GCA','HMT_link','main_genus','main_species','hmt_primary_body_site_w_abundance','hmt_abundance_category','contigs','MAG','assembly_level','checkM2_completeness','checkM2_contamination','quality_score','assembly_complete','Sequence_Source']
rep[cols].to_csv(O/'HOMD_2025_FTP_taxon_representative_panel.tsv',sep='\t',index=False)
print('REPS',len(rep),'primary oral',sum(rep.hmt_primary_body_site_w_abundance=='Oral'),'main sites',rep.hmt_primary_body_site_w_abundance.value_counts().to_dict(),flush=True)
print('FTPTAXA',x.HMT_link.nunique(),'other',len(x),flush=True)
# All HMT ecological labels are broad "Oral" not plaque/mucosa. Six-source habitat metadata must be joined separately.
(O/'metadata_cohort_version_info.json').write_text(json.dumps({'download_date':'2026-10-14','ftp_annotation_date':'2025-07-28','ftp_genome_count':len(ftp),'metadata_genome_count':len(m),'FTP_HMT_count':len(rep),'HMT_primary_oral_count':int(sum(rep.hmt_primary_body_site_w_abundance=='Oral')),'warning':'Primary body site Oral does not differentiate plaque versus mucosa; do not assign niche labels from this field.'},indent=2))
