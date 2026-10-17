# Pipeline inputs and execution scope

Start with the [Methods-to-code index](METHODS_TO_CODE.md) and [run guide](RUNNING.md). Scientific sources live under `analysis/`; `workflows/workflow_catalog.tsv` maps each retained file to its Methods heading, stage number, runtime path, hash and verification coverage.

## Compact reproduction

The package includes cohort metadata, primary matrix, covariates, ecological outcomes, stored null distributions and coordinate subsets needed for `python run.py check`. It requires no HPC connection or private directories. Reference outputs are under `data/frozen/expected`; new outputs go under `build/`. The [data inventory](data_inventory.tsv) records compact inputs and reference files.

## Additional upstream requirements

| Methods section | Additional data and tools |
|---|---|
| Genome retrieval and QC | Accession-matched FASTA, NCBI Datasets; archived CheckM2 metadata |
| Genome annotation | Archived eHOMD and local Prokka outputs, predicted proteins and matching Prokka databases |
| Pfam profiles | HMMER and the archived Pfam-A model snapshot with gathering cutoffs |
| Sequence clusters | Full proteomes; MMseqs2 easy-linclust at identity 0.50, coverage 0.80, coverage mode 0 |
| Rhea / KEGG profiles | Archived EC mappings; module scores measure mapped EC coverage |
| Phylogeny | Representative proteins, marker definitions, MAFFT, RAxML-NG |
| Repertoire and supporting ontology tests | Feature matrices, covariates and matching COG/Pfam-to-GO mappings |
| Complete protein / neighbourhood analysis | Full domain hit tables, sequence clusters, Prokka GFF/FASTA collections |
| Supplementary HxlR / broader HMT survey | Representative genomes, proteomes, GFFs, reference sequences and motif resources |
| Supplementary SnoaL structures | Full sequences, extracted segments, HMMs and prediction services/models |
| Supplementary 14-candidate audit | Broader domain calls, neighbourhoods, ecological labels and EC mappings |

Configure tool environments, databases and scheduler resources for the execution site. Templates retain anonymous `__PROJECT_ROOT__`, `__WORK_ROOT__` and `__HOME_ROOT__` markers; preparation substitutes these in runtime copies. Some upstream tool/database locations require configuration. Syntax validation does not establish availability of raw data, tools or network services.

The release checks do not rerun annotation, all genome downloads, complete phylogenetic inference, all supplementary analyses or every prediction null replicate. See `TEST_REPORT.md` for actual checks. Original research archives remain separate from the curated release.
