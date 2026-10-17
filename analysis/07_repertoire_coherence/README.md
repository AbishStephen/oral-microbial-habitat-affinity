# Repertoire-wide tests of directional coherence

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Deduplicated profiles, genus labels and genome covariates.

## Outputs and manuscript links
Six layers × two scaling choices; Fig. 3d.

## How to use this section
From the repository root: `python run.py prepare`.

This command only prepares runtime copies and frozen inputs. It does not run this section or download raw data. For upstream execution, inspect the scripts and provide the external data and tools listed in [pipeline requirements](../../docs/PIPELINES.md).

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 89 | [supporting_ontology_analyses/89_COG2024_functional_set_enrichment.py](supporting_ontology_analyses/89_COG2024_functional_set_enrichment.py) | syntax only; raw/external computation not replayed |
| 90 | [supporting_ontology_analyses/90_COG2024_codistribution.py](supporting_ontology_analyses/90_COG2024_codistribution.py) | syntax only; raw/external computation not replayed |
| 91 | [supporting_ontology_analyses/91_PFAM_GO_BP_set_enrichment.py](supporting_ontology_analyses/91_PFAM_GO_BP_set_enrichment.py) | syntax only; raw/external computation not replayed |
| 92 | [supporting_ontology_analyses/92_COG2024_primary_category_robustness.py](supporting_ontology_analyses/92_COG2024_primary_category_robustness.py) | syntax only; raw/external computation not replayed |
| 93 | [supporting_ontology_analyses/93_COG2024_primarycat_codistribution.py](supporting_ontology_analyses/93_COG2024_primarycat_codistribution.py) | syntax only; raw/external computation not replayed |
| 94 | [supporting_ontology_analyses/94_global_crossgenus_functional_convergence.py](supporting_ontology_analyses/94_global_crossgenus_functional_convergence.py) | syntax only; raw/external computation not replayed |
| 95 | [supporting_ontology_analyses/95_global_COG_convergence_LOO.py](supporting_ontology_analyses/95_global_COG_convergence_LOO.py) | syntax only; raw/external computation not replayed |
| 97 | [scripts/97_genus_balanced_repertoire_omnibus.py](scripts/97_genus_balanced_repertoire_omnibus.py) | syntax only; raw/external computation not replayed |
| 97 | [scripts/97_genus_balanced_repertoire_omnibus.sh](scripts/97_genus_balanced_repertoire_omnibus.sh) | syntax only; raw/external computation not replayed |
