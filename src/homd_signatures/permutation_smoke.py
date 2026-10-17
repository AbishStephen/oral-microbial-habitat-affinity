"""Run one complete unadjusted and adjusted null refit; compare with archived seed 0."""
from pathlib import Path
import json, os, subprocess, sys, time
import pandas as pd
from .reproduce import REPO, prepare

def main():
    out=REPO/'build/permutation_smoke';p=prepare(out);report=[]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for script,result,reference in [('106_nested_permutation.py','nested_perm_0000.tsv','nested_genusbalanced_200_full_refit_null_distribution.tsv'),('112_covar_fullrefit_permutation.py','covadjusted_nested_perm_0000.tsv','covadjusted_100_withingenus_perm_null_distribution.tsv')]:
        start=time.time()
        with (out/(script+'.log')).open('w') as f:subprocess.run([sys.executable,str(p/'study13/scripts'/script),'0'],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        a=pd.read_csv(p/'study13/11_manuscript_completion_20261013'/result,sep='\t').set_index('label')
        b=pd.read_csv(REPO/'data/frozen/expected'/reference,sep='\t');b=b[b.perm_id==0].set_index('label')
        delta=float((a.macro_auc-b.macro_auc).abs().max());report.append(dict(script=script,seconds=round(time.time()-start,2),max_auc_difference=delta,passed=delta<1e-12))
    (out/'permutation_smoke_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    assert all(x['passed'] for x in report)

if __name__=='__main__':main()
