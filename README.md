# HOMD genomic signatures of oral habitat affinity

Code and data for a comparative genomic study of **13 genera, 47 human microbial taxa (HMTs) and 1,121 genomes**, examining plaque and mucosal habitat affinity. The HMT is the ecological observational unit; genomes estimate within-HMT genomic prevalence. The primary screen compares 13,917 distinct genomic profiles.

## Start here

1. **Follow the paper:** [Methods-to-code index](docs/METHODS_TO_CODE.md) links each Methods heading to its analysis folder, inputs, outputs and run command.
2. **Cross-check a figure:** [figure crosswalk](docs/figure_crosswalk.tsv) links the five main figures and supplementary atlas to analytical stages, source data and R scripts.
3. **Run the checks:** [execution guide](docs/RUNNING.md) explains installation and the compact reproduction commands.
4. **Review the evidence:** [verification report](docs/TEST_REPORT.md) distinguishes repeated calculations from archived results and upstream templates.

## Organisation

| Folder | Methods section |
|---|---|
| `analysis/01_study_design/` | Study design and taxon selection |
| `analysis/02_ecological_phenotype/` | Construction of the ecological habitat-affinity phenotype |
| `analysis/03_genome_retrieval_qc/` | Genome retrieval, quality control and taxonomic assignment |
| `analysis/04_genomic_profiles/` | Genome annotation and HMT-level genomic profiles |
| `analysis/05_phylogeny/` | Phylogenetic reconstruction and sister-pair definition |
| `analysis/06_feature_screen/` | Genus-balanced feature contrasts and permutation testing |
| `analysis/07_repertoire_coherence/` | Repertoire-wide tests of directional coherence |
| `analysis/08_proteins_neighbourhoods/` | Protein sequence clusters and genomic neighbourhood analysis |
| `analysis/09_nested_prediction/` | Nested leave-one-genus-out prediction |
| `analysis/10_feature_exclusion_sensitivity/` | Post hoc feature exclusion and additional ecological sensitivity |
| `analysis/11_computing_reproducibility/` | Computing and reproducibility |

Each analysis folder has a guide. Newer HxlR, SnoaL and 14-candidate investigations are grouped under `analysis/08_proteins_neighbourhoods/supplement/`. Supporting analyses are labelled separately from the main analytical route. The [script catalogue](workflows/workflow_catalog.tsv) records provenance and verification for every retained workflow file.

| Other folder | Contents |
|---|---|
| `data/frozen/` | Fixed inputs and reference outputs for compact reproduction |
| `figure_sources/` | Figure source tables and reference Python plotting code |
| `r/` | Native R scripts for figure generation and styling |
| `results/R/` | Five main figures and supplementary atlas, with plotted source tables |
| `src/`, `tests/` | Reproduction utilities and scientific integrity checks |

## Quick start

Use Python 3.12 on Linux or WSL:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python run.py check
```

This repeats the 4,096-permutation feature screen, four observed nested prediction models, seed-0 full-refitting null checks and supplementary coordinate extraction, and checks archived outputs. It also rebuilds reference Python figures. Raw genomes and a cluster connection are not required for this route. Upstream annotation and additional analyses have separate [input requirements](docs/PIPELINES.md).

## Generate and edit the figures in R

```bash
Rscript r/install_packages.R
Rscript r/verify_inputs.R
Rscript r/build_figures.R
```

Open `HOMD_genomic_signatures.Rproj` for RStudio editing. Change shared appearance in `r/style.R`, or individual layouts in `r/figure1.R` through `r/figure5.R` and `r/supplement.R`. R scripts read bundled data directly and export PDF, SVG and PNG to `results/R/`. See the [figure editing guide](docs/FIGURE_EDITING.md).

## Interpretation and release status

The 14 candidate annotations are an exploratory nominally supported shortlist; none passed all-profile family-wise correction. Sequence clusters, domain annotations and gene proximity do not establish biochemical activity or causality. The broader 631-HMT representative survey is supplementary and reuses ecological sources.

This is a private manuscript repository for review. Public repository, archive DOI, citation metadata and licences will be finalised before publication. [Release review notes](docs/RELEASE_REVIEW.md) record remaining editorial and deposition items.
