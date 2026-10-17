#!/bin/bash
#SBATCH --job-name=HOMD13_INT
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1
#SBATCH --output=__PROJECT_ROOT__/study13/logs/integrate_%j.out
#SBATCH --error=__PROJECT_ROOT__/study13/logs/integrate_%j.err
set -euo pipefail
PY=python3
E=__PROJECT_ROOT__/study13
"$PY" "$E/scripts/56_integrate_benchmark.py"
"$PY" "$E/scripts/57_deduplicate_features.py"
"$PY" "$E/scripts/58_sparse_dedup_benchmark.py"
