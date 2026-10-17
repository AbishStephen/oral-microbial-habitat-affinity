#!/usr/bin/env python3
"""Explicit entry points for the paper's reproducible analyses."""
from pathlib import Path
import argparse, os, subprocess, sys
ROOT=Path(__file__).resolve().parent
os.environ.update(PYTHONPATH=str(ROOT/'src'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
def run(*args):subprocess.run(list(args),cwd=ROOT,check=True)
def module(name,*args):run(sys.executable,'-m','homd_signatures.'+name,*args)
def atlas():
    module('neighborhoods','--coordinates','data/frozen/neighborhoods/candidate_coordinates.tsv.gz','--adjacent','data/frozen/neighborhoods/candidate_immediate_neighbors.tsv.gz','--manifest','data/frozen/project/study13/01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv','--cds-table','data/frozen/neighborhoods/selected_genome_CDS.tsv.gz','--output','build/neighborhoods')
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task',choices=['list','test','cohort','screen','prediction','reproduce','smoke','atlas','figures','prepare','check'])
    parser.add_argument('--figure',default='all',choices=['all','1','2','3','4','5','supplement'])
    a=parser.parse_args()
    if a.task=='list':print((ROOT/'docs/METHODS_TO_CODE.md').read_text());return
    if a.task=='test':run(sys.executable,'-m','unittest','discover','-s','tests','-v')
    elif a.task=='cohort':run(sys.executable,'analysis/01_study_design/rebuild_primary_manifest.py','--metadata','data/frozen/project/00_metadata','--exclusions','data/genome_exclusions.tsv','--output','build/cohort')
    elif a.task in ['screen','prediction','reproduce']:
        module('reproduce','--stages',{'screen':'96','prediction':'104,111','reproduce':'96,104,111'}[a.task])
    elif a.task=='smoke':module('permutation_smoke')
    elif a.task=='atlas':atlas()
    elif a.task=='prepare':module('reproduce','--prepare-only','--output','build/upstream')
    elif a.task=='figures':
        run('Rscript','r/verify_inputs.R');run('Rscript','r/build_figures.R','--figure',a.figure)
    elif a.task=='check':
        run(sys.executable,'-m','unittest','discover','-s','tests','-v')
        module('reproduce');module('permutation_smoke');atlas()
        module('figures','--output','build/reference_python_figures')
        print('Analysis checks passed. Native R figures are checked separately: python run.py figures')
if __name__=='__main__':main()
