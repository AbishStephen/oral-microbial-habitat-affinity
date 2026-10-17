#!/bin/bash
#SBATCH --job-name=SnoaL_long_domains
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=__PROJECT_ROOT__/study13/logs/SnoaL_long_domains_%j.log
set -euo pipefail
S=__WORK_ROOT__/SnoaL2_subfamilies_20261016
H=hmmsearch
mkdir -p "$S/domains"
python3 - <<'PY'
from pathlib import Path
source=Path('__WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm')
D=Path('__WORK_ROOT__/SnoaL2_subfamilies_20261016/domains')
wanted={'PF12680','PF07366','PF13577'}
active=[];acc=''
for line in source.open():
 if line.startswith('HMMER3/f'):active=[line];acc='';continue
 if not active:continue
 active.append(line)
 if line.startswith('ACC '):acc=line.split()[1].split('.')[0]
 if line.strip()=='//':
  if acc in wanted:
   (D/f'{acc}.hmm').write_text(''.join(active))
   wanted.remove(acc)
   print('EXTRACTED',acc,flush=True)
   if not wanted:break
  active=[]
assert not wanted,wanted
PY
for pf in PF12680 PF07366 PF13577; do
 "$H" --cpu 4 --domE 1e-3 --noali --domtblout "$S/domains/$pf.long524.domtblout" -o "$S/domains/$pf.full.log" "$S/domains/$pf.hmm" "$S/SnoaL2_all524_proteins.faa"
done
echo COMPLETED $(date)
