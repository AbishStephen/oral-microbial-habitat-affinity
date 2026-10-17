#!/bin/bash
#SBATCH --job-name=Fig2B_Pfams631
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=__PROJECT_ROOT__/study13/logs/Fig2B_Pfams631_%j.log
set -euo pipefail
B=__PROJECT_ROOT__/study13
O="$B/19_Figure2B_14factor_mechanisms_20261017"
S=__WORK_ROOT__/Figure2B_14factor_20261017
F=__WORK_ROOT__/HOMD_HxlR_fullHOMD_20261014/HOMD_631_HMTrepresentatives_2025_PROKKA.faa
P=__WORK_ROOT__/HOMD_convergence_work/pfam/Pfam-A.hmm
HMMSEARCH=${HMMSEARCH:-hmmsearch}
mkdir -p "$O" "$S"
export B O S P
python3 - <<'PY'
from pathlib import Path
import os
P=Path(os.environ['P']);S=Path(os.environ['S']);profiles=S/'Fig2B_12Pfams.hmm'
wanted={'PF01638','PF04072','PF13577','PF13676','PF22166','PF01745','PF01337','PF13413','PF18075','PF14470','PF12680','PF07366'}
seen=set();buf=[];acc=''
with profiles.open('w') as out:
 for l in P.open():
  if l.startswith('HMMER3/f'):buf=[l];acc='';continue
  if not buf:continue
  buf.append(l)
  if l.startswith('ACC '):acc=l.split()[1].split('.')[0]
  if l.strip()=='//':
   if acc in wanted:
    out.write(''.join(buf));seen.add(acc)
    print('MODEL',acc,flush=True)
    if seen==wanted:break
   buf=[]
assert seen==wanted,(wanted-seen)
PY
"$H/hmmpress" -f "$S/Fig2B_12Pfams.hmm" > "$S/hmmpress.log" 2>&1
"$H/hmmscan" --cpu 8 --cut_ga --noali --domtblout "$S/HOMD631_Fig2B_12Pfams.domtblout" -o "$S/HOMD631_Fig2B_12Pfams.report" "$S/Fig2B_12Pfams.hmm" "$F"
cp "$S/HOMD631_Fig2B_12Pfams.domtblout" "$O/HOMD631_Fig2B_12Pfams.domtblout"
echo DONE $(date)
