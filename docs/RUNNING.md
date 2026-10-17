# Running the analyses

## Choose a route

| Goal | Command from repository root | What it does |
|---|---|---|
| Find code for a Methods section | `python run.py list` | Prints the section index |
| Verify inputs and source integrity | `python run.py test` | Scientific regression, provenance and privacy checks |
| Rebuild cohort and original ecology | `python run.py cohort` | Reconstructs 47 HMTs and 1,121 assemblies from metadata |
| Reproduce the main feature screen | `python run.py screen` | Repeats 4,096 permutations and compares archived results |
| Reproduce observed predictions | `python run.py prediction` | Refits four nested models and compares all HMT predictions |
| Check complete prediction null refits | `python run.py smoke` | Repeats seed 0 for adjusted/unadjusted workflows and both outcomes |
| Reconstruct supplementary tracks | `python run.py atlas` | Extracts 30 loci and 206 gene records from bundled coordinates |
| Generate figures in R | `python run.py figures` | Checks R inputs and builds five main figures and atlas |
| Generate only Fig. 2 in R | `python run.py figures --figure 2` | Checks R inputs and regenerates one figure |
| Prepare upstream scripts | `python run.py prepare` | Resolves anonymous path tokens; executes no upstream stage |
| Run compact Python checks | `python run.py check` | Integrity, screen, predictions, null smoke checks, coordinates and reference Python figures |

## Python environment

Linux is the tested environment. Use Linux or WSL for the full workflow because the integrity suite checks Bash/Slurm templates. Windows can hold and edit the repository. Native Windows execution of the full Python suite has not been verified.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock.txt
python -m pip install --no-deps -e .
python run.py check
```

Python 3.12 was used for recorded reproduction checks. `run.py` resolves its own root and can also be invoked by absolute path. Generated analytical files go into `build/`, leaving frozen inputs intact. Expected outputs are comparison references, not newly generated results.

## R figures and styling

Open `HOMD_genomic_signatures.Rproj` in RStudio; see `FIGURE_EDITING.md` for the supported source-and-build interface. From a terminal:

```bash
Rscript r/install_packages.R
Rscript r/verify_inputs.R
Rscript r/build_figures.R
```

`renv.lock` records tested R package versions. The installer uses a project-local library; system libraries required by Cairo, ragg and svglite must be available. Edit `r/style.R` for shared appearance and the individual `r/figure*.R` files for layouts. Outputs and plotted source tables go to `results/R/`. The R route does not require Python.

## Upstream execution

The repository is self-contained for compact reproduction from processed inputs. `python run.py prepare` stages the 133 retained workflow files at their expected relative runtime paths in `build/upstream/project/study13/scripts/` and writes a source map. Public code remains organised under `analysis/`. This compatibility layer preserves existing cross-script imports and numerical methods.

Do not execute every numbered script in order: some are scheduler wrappers, intermediate comparisons or supplementary investigations. Follow the section guide, inspect prepared script inputs, and supply raw data/tools listed in `PIPELINES.md`. No scheduler jobs or downloads are submitted by `prepare`. Complete 200/100 prediction null distributions are included; `smoke` reruns only seed 0 and `test` recalculates empirical P values from all stored replicates.

The delivered figure set is under `results/R/`. Python plotting scripts are retained for analytical provenance and regression checks; their generated copies are not included as a second figure set.
