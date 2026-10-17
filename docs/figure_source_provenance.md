# Figure source provenance

The five main figures use the included 13-genus study results and plotting inputs. Figures 1 and 5 use full-precision ecological scores and nested predictions from the frozen primary outputs. Figure 2 uses the complete feature-screen statistics and genus contrasts. Figure 3 summarises protein-family diversity and repertoire analyses. Figure 4 uses actual GFF-derived gene coordinates and linkage/sequence-cluster records.

The supplementary candidate-domain atlas uses included domain calls and selected-genome CDS records. Coordinate extraction, representative selection and orientation transformations are executable and checked against the distributed coordinate tables.

All plotting sources use repository-relative inputs. Study taxon, assembly, protein and contig identifiers remain available for scientific traceability.

## Native R figure implementation

`r/figure1.R` through `r/figure5.R` and `r/supplement.R` construct ggplot/patchwork objects directly from the same included scientific inputs. The five main panel sets and supplementary locus selections are retained. R-rendered PDFs, SVGs, PNG previews and plotted-data exports are in `results/R/`; the exact R runtime/package versions are recorded and locked. R JSON parsing treats Python-style non-finite missing annotations as missing; predictions and measured coordinates retain source precision.
