#!/bin/bash
#SBATCH --job-name=SnoaL_PF07366
#SBATCH --partition=compute
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=__PROJECT_ROOT__/study13/logs/SnoaL_PF07366_%j.log
set -euo pipefail
T=__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014
O=__PROJECT_ROOT__/study13/17_SnoaL2_independent_targets_20261016
H=hmmsearch
python3 - <<'PY'
from pathlib import Path
source=Path('__WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm')
dest=Path('__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014/PF07366.hmm')
active=[];tag=''
for line in source.open():
 if line.startswith('HMMER3/f'):active=[line];tag='';continue
 if not active:continue
 active.append(line)
 if line.startswith('ACC'):tag=line.split()[1].split('.')[0]
 if line.startswith('//'):
  if tag=='PF07366':dest.write_text(''.join(active));break
  active=[]
assert dest.stat().st_size>1000
print('EXTRACTED',dest.stat().st_size)
PY
"$H" --cpu 4 --cut_ga --noali --domtblout "$O/HOMD631_PF07366_domtblout.tsv" -o "$T/PF07366_hmmsearch_full_report.txt" "$T/PF07366.hmm" "$T/HOMD_631_HMTrepresentatives_2025_PROKKA.faa"
echo DONE $(date)
