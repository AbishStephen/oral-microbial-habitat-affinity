# Phylogenetic reconstruction and sister-pair definition

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Representative proteomes, 18 marker alignments and species tree.

## Outputs and manuscript links
47-tip phylogeny and eight opposing-habitat sister pairs; Fig. 1 and Extended Data Fig. 1.

## How to use this section
From the repository root: `python run.py prepare`.

This command only prepares runtime copies and frozen inputs. It does not run this section or download raw data. For upstream execution, inspect the scripts and provide the external data and tools listed in [pipeline requirements](../../docs/PIPELINES.md).

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 55 | [scripts/55_build_phylogeny.py](scripts/55_build_phylogeny.py) | syntax only; raw/external computation not replayed |
| 59 | [scripts/59_expanded13_phylogeny.sh](scripts/59_expanded13_phylogeny.sh) | syntax only; raw/external computation not replayed |
| 83 | [scripts/83_terminal_sister_pairs.py](scripts/83_terminal_sister_pairs.py) | syntax only; raw/external computation not replayed |
| 84 | [scripts/84_sister8_exact_maxT.py](scripts/84_sister8_exact_maxT.py) | syntax only; raw/external computation not replayed |
| 85 | [scripts/85_genus_matched_extreme.py](scripts/85_genus_matched_extreme.py) | syntax only; raw/external computation not replayed |
| 85 | [scripts/85_genus_matched_extreme.sh](scripts/85_genus_matched_extreme.sh) | syntax only; raw/external computation not replayed |
| 86 | [scripts/86_pair_residualized_covar.py](scripts/86_pair_residualized_covar.py) | syntax only; raw/external computation not replayed |
