#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
run_dir=${1:?Usage: submit_primary_replay.sh ABSOLUTE_RUN_DIRECTORY [PYTHON_EXECUTABLE]}
python_bin=${2:-python3}
mkdir -p -- "$run_dir"
run_dir=$(cd -- "$run_dir" && pwd)
sbatch --output="$run_dir/slurm-%j.log" "$repo/scripts/slurm/primary_replay.sbatch" "$repo" "$run_dir" "$python_bin"
