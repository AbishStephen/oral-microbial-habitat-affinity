# Nested leave-one-genus-out prediction

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Frozen matrix, both ecological outcomes and genome covariates.

## Outputs and manuscript links
Four observed models, null distributions and uncertainty; Fig. 5.

## How to use this section
From the repository root: `python run.py prediction`.

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 99 | [supporting_analyses/99_existing_nestedLOGO_per_genus.py](supporting_analyses/99_existing_nestedLOGO_per_genus.py) | syntax only; raw/external computation not replayed |
| 104 | [scripts/104_nested_genusbalanced.py](scripts/104_nested_genusbalanced.py) | syntax + complete scientific replay, exact archived-result comparison |
| 104 | [scripts/104_nested_genusbalanced.sbatch](scripts/104_nested_genusbalanced.sbatch) | syntax only; raw/external computation not replayed |
| 106 | [scripts/106_nested_permutation.py](scripts/106_nested_permutation.py) | syntax + seed-0 complete nested permutation refit, exact macro-AUC comparison |
| 106 | [scripts/106_nested_permutation.sbatch](scripts/106_nested_permutation.sbatch) | syntax only; raw/external computation not replayed |
| 111 | [scripts/111_foldwise_covariate_adjusted_prediction.py](scripts/111_foldwise_covariate_adjusted_prediction.py) | syntax + complete scientific replay, exact archived-result comparison |
| 111 | [scripts/111_foldwise_covariate_adjusted_prediction.sbatch](scripts/111_foldwise_covariate_adjusted_prediction.sbatch) | syntax only; raw/external computation not replayed |
| 112 | [scripts/112_covar_fullrefit_permutation.py](scripts/112_covar_fullrefit_permutation.py) | syntax + seed-0 complete nested permutation refit, exact macro-AUC comparison |
| 112 | [scripts/112_covar_fullrefit_permutation.sbatch](scripts/112_covar_fullrefit_permutation.sbatch) | syntax only; raw/external computation not replayed |
| 113 | [scripts/113_genusbootstrap.py](scripts/113_genusbootstrap.py) | syntax only; raw/external computation not replayed |
