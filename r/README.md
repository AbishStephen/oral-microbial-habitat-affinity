# Figure scripts in R

This folder contains the standalone figure-generation code. Statistical analysis pipelines are organised separately under `analysis/`.

| Output | Script |
|---|---|
| Figure 1 | `figure1.R` |
| Figure 2 | `figure2.R` |
| Figure 3 | `figure3.R` |
| Figure 4 | `figure4.R` |
| Figure 5 | `figure5.R` |
| Supplementary genomic neighbourhood atlas | `supplement.R` |

From the repository root, run `Rscript r/build_figures.R` for all figures, or add `--figure 2` for one figure. `style.R` controls shared appearance; `common.R` provides data access, themes and export helpers. Installation and scientific input checks are in `install_packages.R` and `verify_inputs.R`.

PDF, SVG and PNG exports are written to `results/R/`. Plotting inputs remain in `figure_sources/` and `data/frozen/`; the R workflow does not call Python. See [the editing guide](../docs/FIGURE_EDITING.md).

The displayed text was reconciled against the author's edited figure PDFs on 17 October 2026. Removed overall titles, selected subtitles and captions are omitted in the code. Scientific values, panel letters, axes and retained labels are preserved. Empty text strings intentionally preserve the spacing of the main-figure layout.
