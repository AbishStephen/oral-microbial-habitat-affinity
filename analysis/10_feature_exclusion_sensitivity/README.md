# Post hoc feature exclusion and additional ecological sensitivity

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Feature matrix, exclusion lists and source-specific abundance summaries.

## Outputs and manuscript links
Feature exclusion and five site-contrast sensitivity results.

## How to use this section
From the repository root: `python run.py prepare`.

This command only prepares runtime copies and frozen inputs. It does not run this section or download raw data. For upstream execution, inspect the scripts and provide the external data and tools listed in [pipeline requirements](../../docs/PIPELINES.md).

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 78 | [scripts/78_site_gradient_comparison.py](scripts/78_site_gradient_comparison.py) | syntax only; raw/external computation not replayed |
| 79 | [scripts/79_hypothesis_guided_programs.py](scripts/79_hypothesis_guided_programs.py) | syntax only; raw/external computation not replayed |
| 79 | [scripts/79_hypothesis_guided_programs.sh](scripts/79_hypothesis_guided_programs.sh) | syntax only; raw/external computation not replayed |
| 105 | [scripts/105_all14_site_sensitivity.py](scripts/105_all14_site_sensitivity.py) | syntax only; raw/external computation not replayed |
| 108 | [scripts/108_ablation.py](scripts/108_ablation.py) | syntax only; raw/external computation not replayed |
