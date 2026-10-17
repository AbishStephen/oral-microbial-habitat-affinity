# Genus-balanced feature contrasts and permutation testing

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Frozen 47-HMT matrix, habitat labels and genome covariates.

## Outputs and manuscript links
119 raw nominal, 23 jointly nominal, 14 breadth-filtered candidates; Fig. 2.

## How to use this section
From the repository root: `python run.py screen`.

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 58 | [supporting_analyses/58_sparse_dedup_benchmark.py](supporting_analyses/58_sparse_dedup_benchmark.py) | syntax only; raw/external computation not replayed |
| 60 | [supporting_analyses/60_expanded13_permutation_null.py](supporting_analyses/60_expanded13_permutation_null.py) | syntax only; raw/external computation not replayed |
| 60 | [supporting_analyses/60_expanded13_permutation_null_array.sh](supporting_analyses/60_expanded13_permutation_null_array.sh) | syntax only; raw/external computation not replayed |
| 61 | [supporting_analyses/61_expanded13_within_genus.py](supporting_analyses/61_expanded13_within_genus.py) | syntax only; raw/external computation not replayed |
| 61 | [supporting_analyses/61_expanded13_within_genus.sh](supporting_analyses/61_expanded13_within_genus.sh) | syntax only; raw/external computation not replayed |
| 62 | [supporting_analyses/62_target_domains.py](supporting_analyses/62_target_domains.py) | syntax only; raw/external computation not replayed |
| 63 | [supporting_analyses/63_confound_audit.py](supporting_analyses/63_confound_audit.py) | syntax only; raw/external computation not replayed |
| 65 | [supporting_analyses/65_within_genus_perm_null.py](supporting_analyses/65_within_genus_perm_null.py) | syntax only; raw/external computation not replayed |
| 66 | [supporting_analyses/66_adjusted_maxt.py](supporting_analyses/66_adjusted_maxt.py) | syntax only; raw/external computation not replayed |
| 70 | [supporting_analyses/70_pgls47_covariate.py](supporting_analyses/70_pgls47_covariate.py) | syntax only; raw/external computation not replayed |
| 87 | [supporting_analyses/87_broad_feature_candidate_slate.py](supporting_analyses/87_broad_feature_candidate_slate.py) | syntax only; raw/external computation not replayed |
| 88 | [supporting_analyses/88_export_broad_feature_atlas.py](supporting_analyses/88_export_broad_feature_atlas.py) | syntax only; raw/external computation not replayed |
| 96 | [scripts/96_crossgenus_balanced_permutation_scan.py](scripts/96_crossgenus_balanced_permutation_scan.py) | syntax + complete scientific replay, exact archived-result comparison |
| 96 | [scripts/96_crossgenus_balanced_permutation_scan.sh](scripts/96_crossgenus_balanced_permutation_scan.sh) | syntax only; raw/external computation not replayed |
| 98 | [scripts/98_crossgenus_integrated_modules.py](scripts/98_crossgenus_integrated_modules.py) | syntax only; raw/external computation not replayed |
