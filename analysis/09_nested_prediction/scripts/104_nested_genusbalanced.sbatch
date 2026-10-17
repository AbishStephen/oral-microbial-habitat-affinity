#!/bin/bash
#SBATCH -J HOMD_NESTED
#SBATCH -p computeshort
#SBATCH -t 00:55:00
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -o __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/nested_104_%j.out
#SBATCH -e __PROJECT_ROOT__/study13/11_manuscript_completion_20261013/nested_104_%j.err
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 __PROJECT_ROOT__/study13/scripts/104_nested_genusbalanced.py
