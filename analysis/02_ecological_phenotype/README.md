# Construction of the ecological habitat-affinity phenotype

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Six abundance summaries and HMT manifest.

## Outputs and manuscript links
Original and source-family-balanced habitat scores; Figs. 1 and 5.

## How to use this section
From the repository root: `python run.py prepare`.

This command only prepares runtime copies and frozen inputs. It does not run this section or download raw data. For upstream execution, inspect the scripts and provide the external data and tools listed in [pipeline requirements](../../docs/PIPELINES.md).

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 67 | [scripts/67_source_family_sensitivity.py](scripts/67_source_family_sensitivity.py) | syntax only; raw/external computation not replayed |
| 68 | [scripts/68_source_balanced_sparse_benchmark.py](scripts/68_source_balanced_sparse_benchmark.py) | syntax only; raw/external computation not replayed |
