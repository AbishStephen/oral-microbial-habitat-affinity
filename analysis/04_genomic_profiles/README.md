# Genome annotation and HMT-level genomic profiles

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Prokka proteins/tables, Pfam models, EC mappings and MMseqs2 clusters.

## Outputs and manuscript links
HMT prevalence matrices and 13,917 distinct profiles; Figs. 2–3.

## How to use this section
From the repository root: `python run.py prepare`.

This command only prepares runtime copies and frozen inputs. It does not run this section or download raw data. For upstream execution, inspect the scripts and provide the external data and tools listed in [pipeline requirements](../../docs/PIPELINES.md).

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 46 | [scripts/46_prokka_expanded13_array.sh](scripts/46_prokka_expanded13_array.sh) | syntax only; raw/external computation not replayed |
| 47 | [scripts/47_build_protein_catalogue.py](scripts/47_build_protein_catalogue.py) | syntax only; raw/external computation not replayed |
| 47 | [scripts/47_postprokka.sh](scripts/47_postprokka.sh) | syntax only; raw/external computation not replayed |
| 47 | [scripts/47b_rebuild_fasta_chunks.sh](scripts/47b_rebuild_fasta_chunks.sh) | syntax only; raw/external computation not replayed |
| 48 | [scripts/48_build_function_matrices.py](scripts/48_build_function_matrices.py) | syntax only; raw/external computation not replayed |
| 49 | [scripts/49_mmseqs_job.sh](scripts/49_mmseqs_job.sh) | syntax only; raw/external computation not replayed |
| 49 | [scripts/49_mmseqs_linclust.sh](scripts/49_mmseqs_linclust.sh) | syntax only; raw/external computation not replayed |
| 50 | [scripts/50_build_X5_matrices.py](scripts/50_build_X5_matrices.py) | syntax only; raw/external computation not replayed |
| 50 | [scripts/50_x5_matrix_job.sh](scripts/50_x5_matrix_job.sh) | syntax only; raw/external computation not replayed |
| 51 | [scripts/51_pfam_split_job.sh](scripts/51_pfam_split_job.sh) | syntax only; raw/external computation not replayed |
| 51 | [scripts/51_split_pfam_fasta.py](scripts/51_split_pfam_fasta.py) | syntax only; raw/external computation not replayed |
| 52 | [scripts/52_pfam_scan_array.sh](scripts/52_pfam_scan_array.sh) | syntax only; raw/external computation not replayed |
| 53 | [scripts/53_build_Pfam_matrices.py](scripts/53_build_Pfam_matrices.py) | syntax only; raw/external computation not replayed |
| 53 | [scripts/53_pfam_matrix_job.sh](scripts/53_pfam_matrix_job.sh) | syntax only; raw/external computation not replayed |
| 54 | [scripts/54_build_X3_X4.py](scripts/54_build_X3_X4.py) | syntax only; raw/external computation not replayed |
| 54 | [scripts/54_x34_job.sh](scripts/54_x34_job.sh) | syntax only; raw/external computation not replayed |
| 56 | [scripts/56_integrate_benchmark.py](scripts/56_integrate_benchmark.py) | syntax only; raw/external computation not replayed |
| 56 | [scripts/56_integrate_job.sh](scripts/56_integrate_job.sh) | syntax only; raw/external computation not replayed |
| 57 | [scripts/57_deduplicate_features.py](scripts/57_deduplicate_features.py) | syntax only; raw/external computation not replayed |
