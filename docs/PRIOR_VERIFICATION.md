# Verification report

Verified with Python 3.12.14 and the pinned environment in `requirements.lock.txt`. The complete `make check` command passed in the anonymised publication package and in a fresh Git checkout relocated to a path containing spaces.

- Scientific integrity and privacy checks: 11 tests passed; source checksums and syntax checked across 143 publication workflow templates.
- Primary cohort reconstruction: 13 genera, 47 HMTs and 1,121 genomes, including the documented assembly exclusions.
- Full 4,096-permutation screen: 13,917 profiles, 119 raw nominal candidates, 23 nominal in both screens, the same 14-feature shortlist and zero all-profile family-wise hits. Candidate permutation P values matched exactly.
- Four observed nested prediction runs: every HMT prediction matched the stored outputs exactly. Macro AUCs were 0.8076923076923077 and 0.8205128205128205 for unadjusted models, and 0.7243589743589743 and 0.7756410256410257 for covariate-adjusted models.
- Full-refitting permutation checks: one complete seed-0 unadjusted refit and one adjusted refit matched the stored null AUCs exactly. The stored complete 200- and 100-permutation null distributions passed integrity checks and empirical P-value recalculation.
- Supplementary coordinates: all 30 loci and 206 gene records reconstructed from the included coordinate inputs; strand and length transformations checked.
- Figure rendering: all five main figures and the supplementary gene-arrangement figure rebuilt as PDF, SVG and PNG. All PNG files passed full image decoding; SVG XML and PDF end markers were also validated.

Machine-readable evidence is in `docs/test_evidence/`. The hosted CI configuration is included but has not been executed on GitHub. Upstream raw-data and external-service stages receive syntax checks; they were not all executed in this verification. See `PIPELINES.md` for those requirements.

## Native R figures

All five main figures and the supplementary atlas rendered in R 4.2.0 with ggplot2 3.5.1, patchwork 1.2.0, jsonlite 1.8.8, ape 5.8, svglite 2.2.2 and ragg 1.3.2. The exact dependency graph is in `renv.lock`. Ten native R scientific-input checks passed. The full R build also passed after moving the project to a directory containing spaces and running it from a different working directory. The source-and-build interface used for RStudio editing was checked.

Vector exports were parsed and checked for editable SVG text; plotted data were compared against the frozen source values, including full-precision predictions and all supplementary coordinate fields. Source-table round-trip checks verify 206 gene records and 47 predictions. PNG previews were rasterized from validated R PDFs; the R scripts also export PNG directly with ragg. See `docs/test_evidence/R_figure_verification.json`. The lockfile was generated from the tested environment; a full dependency reinstall on a clean machine was not performed.
