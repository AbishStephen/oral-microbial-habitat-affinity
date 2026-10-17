PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/src
export OPENBLAS_NUM_THREADS := 1
export OMP_NUM_THREADS := 1
export MKL_NUM_THREADS := 1
export MPLBACKEND := Agg

.PHONY: test reproduce permutation-smoke supplement main-figures check
test:
	$(PYTHON) -m unittest discover -s tests -v

reproduce:
	$(PYTHON) -m homd_signatures.reproduce --output build/reproduction

permutation-smoke:
	$(PYTHON) -m homd_signatures.permutation_smoke

supplement:
	$(PYTHON) -m homd_signatures.neighborhoods --coordinates data/frozen/neighborhoods/candidate_coordinates.tsv.gz --adjacent data/frozen/neighborhoods/candidate_immediate_neighbors.tsv.gz --manifest data/frozen/project/study13/01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv --cds-table data/frozen/neighborhoods/selected_genome_CDS.tsv.gz --output build/neighborhoods
	$(PYTHON) -m homd_signatures.plot_neighborhoods --input build/neighborhoods --output build/reference_python_supplement

main-figures:
	$(PYTHON) -m homd_signatures.figures --output build/main_figures

check: test reproduce permutation-smoke supplement main-figures

.PHONY: r-figures r-check
r-figures:
	Rscript r/build_figures.R

r-check:
	Rscript r/verify_inputs.R
	Rscript r/build_figures.R
