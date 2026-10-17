# Repository organisation

The layout borrows documented practices from published research repositories, while using this manuscript's Methods headings as its primary navigation.

- MetaModels (Marcelino et al., Nature Communications 2023): analysis workflow instructions and separately identifiable R figure scripts. https://github.com/vrmarcelino/MetaModels
- RAT_paper (Hauptfeld et al., Nature Communications 2024): explicit computational dependencies and downstream analysis/plotting source. https://github.com/thauptfeld/RAT_paper

Here, section guides connect Methods to code; the catalogue records original stage IDs, runtime destinations and SHA-256 hashes; figure_crosswalk.tsv links figures to analyses, source tables and R scripts. The compact reproduction route is separated from upstream templates that require raw genomes or external tools. A flat folder of dated scripts is not the reader-facing interface.
