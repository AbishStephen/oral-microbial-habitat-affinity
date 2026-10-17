# Genome retrieval, quality control and taxonomic assignment

[Methods index](../../docs/METHODS_TO_CODE.md) · [Run guide](../../docs/RUNNING.md)

## Inputs
Accession manifest, public assembly downloads and CheckM2 metadata.

## Outputs and manuscript links
Retained genome accessions and annotation inputs; Fig. 1.

## How to use this section
From the repository root: `python run.py cohort`.

The cohort command reconstructs selection from supplied metadata. Assembly downloads and annotation are separate upstream steps. CheckM2 values are archived metadata, not a newly repeated QC run.

## Script inventory

Stage numbers preserve provenance; they are not a command to run every file numerically. Shell files are execution wrappers, not additional statistical analyses. Supporting analyses and supplementary follow-ups are explicitly separated.

| Stage | Script | Verification coverage |
|---|---|---|
| 44 | [scripts/44_download_expanded13_array.sh](scripts/44_download_expanded13_array.sh) | syntax only; raw/external computation not replayed |
| 45 | [scripts/45_check_expanded13_downloads.sh](scripts/45_check_expanded13_downloads.sh) | syntax only; raw/external computation not replayed |
