# Protein sequence clusters and genomic neighbourhood analysis

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Protein/cluster dictionaries, Pfam calls, GFF coordinates and ecology.

## Outputs and manuscript links
Protein repertoires, linked loci and supplementary candidate audit; Figs. 3–4.

## How to use this section
From the repository root: `python run.py atlas`.

`atlas` reconstructs the compact supplementary gene tracks from supplied coordinates. It does not repeat every protein, motif or structure analysis. The `supplement/` folders contain the additional HxlR, SnoaL and 14-candidate investigations; these support supplementary interpretation, not extra main figures.

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 64 | [scripts/64_domain_colocalization.py](scripts/64_domain_colocalization.py) | syntax only; raw/external computation not replayed |
| 69 | [scripts/69_coloc_all_near.py](scripts/69_coloc_all_near.py) | syntax only; raw/external computation not replayed |
| 77 | [scripts/77_colocalization_genus.py](scripts/77_colocalization_genus.py) | syntax only; raw/external computation not replayed |
| 100 | [scripts/100_crossgenus_candidate_physical_linkage.py](scripts/100_crossgenus_candidate_physical_linkage.py) | syntax only; raw/external computation not replayed |
| 100 | [scripts/100_crossgenus_candidate_physical_linkage.sh](scripts/100_crossgenus_candidate_physical_linkage.sh) | syntax only; raw/external computation not replayed |
| 101 | [scripts/101_neighborhood_hmt_ecology.py](scripts/101_neighborhood_hmt_ecology.py) | syntax only; raw/external computation not replayed |
| 102 | [scripts/102_candidate_neighborhood_orientation.py](scripts/102_candidate_neighborhood_orientation.py) | syntax only; raw/external computation not replayed |
| 103 | [scripts/103_linked_HxlR_Nterminal_motif_check.py](scripts/103_linked_HxlR_Nterminal_motif_check.py) | syntax only; raw/external computation not replayed |
| 107 | [scripts/107_protein_evidence.py](scripts/107_protein_evidence.py) | syntax only; raw/external computation not replayed |
| 109 | [scripts/109_orthology_specificity.py](scripts/109_orthology_specificity.py) | syntax only; raw/external computation not replayed |
| 114 | [supplement/hxlr_context/114_HxlR_regulatory_region_audit.py](supplement/hxlr_context/114_HxlR_regulatory_region_audit.py) | syntax only; raw/external computation not replayed |
| 115 | [supplement/hxlr_context/115_HxlR_recover_prokka_contigs.py](supplement/hxlr_context/115_HxlR_recover_prokka_contigs.py) | syntax only; raw/external computation not replayed |
| 116 | [supplement/hxlr_context/116_HxlR_intergenic_diversity.py](supplement/hxlr_context/116_HxlR_intergenic_diversity.py) | syntax only; raw/external computation not replayed |
| 117 | [supplement/hxlr_context/117_HxlR_crossgenus_exact_tracts.py](supplement/hxlr_context/117_HxlR_crossgenus_exact_tracts.py) | syntax only; raw/external computation not replayed |
| 118 | [supplement/hxlr_context/118_HxlR_all_genomes_survey.py](supplement/hxlr_context/118_HxlR_all_genomes_survey.py) | syntax only; raw/external computation not replayed |
| 119 | [supplement/hxlr_context/119_Porphyromonas_Prevotella_cassette_census.py](supplement/hxlr_context/119_Porphyromonas_Prevotella_cassette_census.py) | syntax only; raw/external computation not replayed |
| 120 | [supplement/hxlr_context/120_HxlR_fetch_missing_genomes.py](supplement/hxlr_context/120_HxlR_fetch_missing_genomes.py) | syntax only; raw/external computation not replayed |
| 121 | [supplement/hxlr_context/121_HxlR_allgenome_neighborhoods.py](supplement/hxlr_context/121_HxlR_allgenome_neighborhoods.py) | syntax only; raw/external computation not replayed |
| 122 | [supplement/hxlr_context/122_HxlR_prepare_allgenome_align.py](supplement/hxlr_context/122_HxlR_prepare_allgenome_align.py) | syntax only; raw/external computation not replayed |
| 123 | [supplement/hxlr_context/123_HxlR_homology.sbatch](supplement/hxlr_context/123_HxlR_homology.sbatch) | syntax only; raw/external computation not replayed |
| 124 | [supplement/hxlr_context/124_HxlR_family_motif_audit.py](supplement/hxlr_context/124_HxlR_family_motif_audit.py) | syntax only; raw/external computation not replayed |
| 125 | [supplement/hxlr_context/125_HxlR_gene_pair_recount.py](supplement/hxlr_context/125_HxlR_gene_pair_recount.py) | syntax only; raw/external computation not replayed |
| 126 | [supplement/hxlr_context/126_HxlR_alignments_summary.py](supplement/hxlr_context/126_HxlR_alignments_summary.py) | syntax only; raw/external computation not replayed |
| 127 | [supplement/hxlr_context/127_HxlR_assembly_provenance.py](supplement/hxlr_context/127_HxlR_assembly_provenance.py) | syntax only; raw/external computation not replayed |
| 128 | [supplement/broader_hmt_survey/128_HOMD_representative_panel.py](supplement/broader_hmt_survey/128_HOMD_representative_panel.py) | syntax only; raw/external computation not replayed |
| 129 | [supplement/broader_hmt_survey/129_HOMD_download_representatives.py](supplement/broader_hmt_survey/129_HOMD_download_representatives.py) | syntax only; raw/external computation not replayed |
| 130 | [supplement/broader_hmt_survey/130_HxlR_snoal_context_byHMT.py](supplement/broader_hmt_survey/130_HxlR_snoal_context_byHMT.py) | syntax only; raw/external computation not replayed |
| 131 | [supplement/broader_hmt_survey/131_HxlR_prepare_631_HMM.py](supplement/broader_hmt_survey/131_HxlR_prepare_631_HMM.py) | syntax only; raw/external computation not replayed |
| 132 | [supplement/broader_hmt_survey/132_HxlR_HOMD631_hmmsearch.sbatch](supplement/broader_hmt_survey/132_HxlR_HOMD631_hmmsearch.sbatch) | syntax only; raw/external computation not replayed |
| 133 | [supplement/broader_hmt_survey/133_HxlR_matched_neighbor_null.py](supplement/broader_hmt_survey/133_HxlR_matched_neighbor_null.py) | syntax only; raw/external computation not replayed |
| 134 | [supplement/broader_hmt_survey/134_HOMD631_hmm_survey.py](supplement/broader_hmt_survey/134_HOMD631_hmm_survey.py) | syntax only; raw/external computation not replayed |
| 135 | [supplement/broader_hmt_survey/135_HOMD_GFF_fetch_both.py](supplement/broader_hmt_survey/135_HOMD_GFF_fetch_both.py) | syntax only; raw/external computation not replayed |
| 136 | [supplement/broader_hmt_survey/136_HOMD631_pair_synteny.py](supplement/broader_hmt_survey/136_HOMD631_pair_synteny.py) | syntax only; raw/external computation not replayed |
| 137 | [supplement/broader_hmt_survey/137_HOMD631_ecology_check.py](supplement/broader_hmt_survey/137_HOMD631_ecology_check.py) | syntax only; raw/external computation not replayed |
| 138 | [supplement/broader_hmt_survey/138_HOMD631_niche_genus_permutations.py](supplement/broader_hmt_survey/138_HOMD631_niche_genus_permutations.py) | syntax only; raw/external computation not replayed |
| 139 | [supplement/broader_hmt_survey/139_HOMD631_threshold_sensitivity.py](supplement/broader_hmt_survey/139_HOMD631_threshold_sensitivity.py) | syntax only; raw/external computation not replayed |
| 140 | [supplement/hxlr_context/140_HxlR_regulator_target_prioritization.py](supplement/hxlr_context/140_HxlR_regulator_target_prioritization.py) | syntax only; raw/external computation not replayed |
| 141 | [supplement/hxlr_context/141_HxlR_candidate_promoter_extract.py](supplement/hxlr_context/141_HxlR_candidate_promoter_extract.py) | syntax only; raw/external computation not replayed |
| 142 | [supplement/hxlr_context/142_HxlR_curated_reference_extract.py](supplement/hxlr_context/142_HxlR_curated_reference_extract.py) | syntax only; raw/external computation not replayed |
| 143 | [supplement/hxlr_context/143_HxlR_reference_backbone.sbatch](supplement/hxlr_context/143_HxlR_reference_backbone.sbatch) | syntax only; raw/external computation not replayed |
| 144 | [supplement/hxlr_context/144_HxlR_align_all_to_anchors.py](supplement/hxlr_context/144_HxlR_align_all_to_anchors.py) | syntax only; raw/external computation not replayed |
| 145 | [supplement/hxlr_context/145_HxlR_regprecise_motifs.py](supplement/hxlr_context/145_HxlR_regprecise_motifs.py) | syntax only; raw/external computation not replayed |
| 146 | [supplement/hxlr_context/146_HxlR_RegPrecise_shuffle_motifs.py](supplement/hxlr_context/146_HxlR_RegPrecise_shuffle_motifs.py) | syntax only; raw/external computation not replayed |
| 147 | [supplement/hxlr_context/147_HxlR_project_reference_tree.py](supplement/hxlr_context/147_HxlR_project_reference_tree.py) | syntax only; raw/external computation not replayed |
| 148 | [supplement/snoal_architecture/148_SnoaL2_architecture_ecology.py](supplement/snoal_architecture/148_SnoaL2_architecture_ecology.py) | syntax only; raw/external computation not replayed |
| 149 | [supplement/snoal_architecture/149_SnoaL2_47HMT_genus_perm.py](supplement/snoal_architecture/149_SnoaL2_47HMT_genus_perm.py) | syntax only; raw/external computation not replayed |
| 150 | [supplement/snoal_architecture/150_SnoaL_HOMD631_PF07366.sbatch](supplement/snoal_architecture/150_SnoaL_HOMD631_PF07366.sbatch) | syntax only; raw/external computation not replayed |
| 151 | [supplement/snoal_architecture/151_SnoaL_HOMD631_architecture_ecology.py](supplement/snoal_architecture/151_SnoaL_HOMD631_architecture_ecology.py) | syntax only; raw/external computation not replayed |
| 152 | [supplement/snoal_architecture/152_HxlR_four_target_modules.py](supplement/snoal_architecture/152_HxlR_four_target_modules.py) | syntax only; raw/external computation not replayed |
| 153 | [supplement/snoal_architecture/153_SnoaL2_proteome_extract.py](supplement/snoal_architecture/153_SnoaL2_proteome_extract.py) | syntax only; raw/external computation not replayed |
| 154 | [supplement/snoal_architecture/154_SnoaL2_clusters.sbatch](supplement/snoal_architecture/154_SnoaL2_clusters.sbatch) | syntax only; raw/external computation not replayed |
| 155 | [supplement/snoal_architecture/155_SnoaL_integrity.py](supplement/snoal_architecture/155_SnoaL_integrity.py) | syntax only; raw/external computation not replayed |
| 156 | [supplement/snoal_architecture/156_snoal_long_reps.py](supplement/snoal_architecture/156_snoal_long_reps.py) | syntax only; raw/external computation not replayed |
| 157 | [supplement/snoal_architecture/157_SnoaL_long_domain_hmms.sbatch](supplement/snoal_architecture/157_SnoaL_long_domain_hmms.sbatch) | syntax only; raw/external computation not replayed |
| 158 | [supplement/snoal_architecture/158_snoal_resolved_layouts.py](supplement/snoal_architecture/158_snoal_resolved_layouts.py) | syntax only; raw/external computation not replayed |
| 159 | [supplement/snoal_architecture/159_SnoaL5_allPfam.sbatch](supplement/snoal_architecture/159_SnoaL5_allPfam.sbatch) | syntax only; raw/external computation not replayed |
| 160 | [supplement/snoal_architecture/160_SnoaL_actino17_Pfam.sbatch](supplement/snoal_architecture/160_SnoaL_actino17_Pfam.sbatch) | syntax only; raw/external computation not replayed |
| 161 | [supplement/snoal_architecture/161_SnoaL_ESMFold_predictions.py](supplement/snoal_architecture/161_SnoaL_ESMFold_predictions.py) | syntax only; raw/external computation not replayed |
| 162 | [supplement/snoal_architecture/162_SnoaL_ESMFold_structure_statistics.py](supplement/snoal_architecture/162_SnoaL_ESMFold_structure_statistics.py) | syntax only; raw/external computation not replayed |
| 163 | [supplement/fourteen_candidates/163_Fig2B_Pfam631.sbatch](supplement/fourteen_candidates/163_Fig2B_Pfam631.sbatch) | syntax only; raw/external computation not replayed |
| 164 | [supplement/fourteen_candidates/164_Fig2B_10Pfams_neighborhoods.py](supplement/fourteen_candidates/164_Fig2B_10Pfams_neighborhoods.py) | syntax only; raw/external computation not replayed |
| 165 | [supplement/fourteen_candidates/165_Fig2B_HOMD631_validation.py](supplement/fourteen_candidates/165_Fig2B_HOMD631_validation.py) | syntax only; raw/external computation not replayed |
| 166 | [supplement/fourteen_candidates/166_Fig2B_all14_evidence_audit.py](supplement/fourteen_candidates/166_Fig2B_all14_evidence_audit.py) | syntax only; raw/external computation not replayed |
| 167 | [supplement/fourteen_candidates/167_Fig2B_cassette_domain_check.py](supplement/fourteen_candidates/167_Fig2B_cassette_domain_check.py) | syntax only; raw/external computation not replayed |
