# Editing figures in R

All five main figures and the supplementary candidate-domain atlas have native R plotting implementations using ggplot2 and patchwork. They load scientific source tables directly; they do not call Python or wrap an existing image.

## Setup and rendering

Open `HOMD_genomic_signatures.Rproj` in RStudio, or work from the repository root:

```bash
Rscript r/install_packages.R
Rscript r/verify_inputs.R
Rscript r/build_figures.R
```

`renv.lock` records the tested R runtime and exact dependency versions. The installer restores them into `.R_library/`, which is ignored by Git and loaded automatically by the R plotting helpers. System libraries for font rendering, Cairo, ragg and svglite may be needed when installing from source. If the required packages are already installed, the build can be run directly.

Outputs go to `results/R/`. Every figure has PDF, SVG and PNG exports. SVG text remains editable. `source_data/` contains the actual tables passed to the R panels, with numeric columns exported at 17 significant digits.

## Shared appearance

Edit `r/style.R` to change the font family, base size, title sizes, plaque/mucosal colours, phylum palette and gene colours. Choose a font installed on the rendering computer. Explicit local sizes in individual scripts override shared defaults and can be adjusted there.

| Figure | R plotting source | Main inputs |
|---|---|---|
| 1 | `r/figure1.R` | Full-precision HMT manifest and branch-length tree |
| 2 | `r/figure2.R` | Fourteen-candidate statistics and 41-feature genus effects |
| 3 | `r/figure3.R` | Domain breadth, sequence-cluster membership, co-annotation and omnibus tables |
| 4 | `r/figure4.R` | Real gene-coordinate loci, all 67 linkage distances and eight cluster-pair types |
| 5 | `r/figure5.R` | Full-precision held-out predictions, genus AUCs and uncertainty estimates |
| Supplementary atlas | `r/supplement.R` | Thirty selected loci and 206 coordinate records |

To rebuild only Figure 2:

```bash
Rscript r/build_figures.R --figure 2
```

To choose a separate output directory:

```bash
Rscript r/build_figures.R --figure 2 --output results/figure2_edited
```

In an RStudio session opened at the project root:

```r
source("r/figure2.R")
p <- build_figure2()
print(p)
write_plot(p, "Figure2", width = 19.2, height = 13.8)
```

The returned object is a ggplot or patchwork composition. Edit panel labels, scales, legends, margins and layout proportions in its figure script. The title/caption text is also editable there. Shared text defaults and palettes are in `style.R`; gene labels have a separate shared size setting.

## Scientific checks

`r/verify_inputs.R` checks the final cohort, full-precision ecological scores and predictions, reconstructed macro AUC, all-14 exploratory statistics, supplementary locus counts, strand transformations and coordinate lengths. It compares the plotting inputs against frozen primary results with a 1e-12 numeric tolerance.

R source tables and polygons preserve the same scientific data and main panel subjects. The R layouts are editable reproductions, rather than pixel-identical copies of the Python layouts. Missing annotation entries encoded as unquoted NaN in a source JSON are treated as missing values by the R reader; numerical observations are retained.

The Python plotting scripts remain available as reference implementations alongside the Python analysis pipelines. The R workflow can be used independently for figure editing and export.

## Author PDF text edits

The six R figure scripts were updated to match text removed in the author-edited PDFs on 17 October 2026. See `FIGURE_TEXT_EDITS.md` for the precise deletions. The main-figure empty title strings preserve spacing; retained scientific labels and values are unchanged. Further changes should be made in the R sources so all export formats remain consistent.
