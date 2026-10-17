"""Replay the archived final analysis code in a separate, relocatable workspace.

Publication workflow templates use anonymous path tokens. Tokens are resolved
in separate runtime copies, including scripts read by other scripts.
"""
from pathlib import Path
import argparse, csv, json, os, shutil, subprocess, sys, time
import numpy as np
import pandas as pd

REPO=Path(__file__).resolve().parents[2]

def prepare(output):
    output=Path(output).resolve();project=output/'project';scratch=output/'scratch'
    shutil.copytree(REPO/'data/frozen/project',project,dirs_exist_ok=True)
    scratch.mkdir(parents=True,exist_ok=True)
    replaced=[]
    with (REPO/'workflows/workflow_catalog.tsv').open() as f:
        catalog=list(csv.DictReader(f,delimiter='\t'))
    for row in catalog:
        p=REPO/row['script'];rel=Path(row['runtime_path']);dest=project/rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        text=p.read_text()
        text=text.replace('__PROJECT_ROOT__',str(project))
        text=text.replace('__WORK_ROOT__',str(scratch))
        text=text.replace('__HOME_ROOT__',str(output/'home'))
        dest.write_text(text);replaced.append(dict(source=row['script'],runtime=str(rel)))
    for d in ['10_additional_analyses_20261012','11_manuscript_completion_20261013']:
        (project/'study13'/d).mkdir(parents=True,exist_ok=True)
    (output/'runtime_source_map.json').write_text(json.dumps(dict(source='analysis',runtime_root=str(project),scripts=replaced,method='resolve anonymous project/work/home tokens in runtime copies'),indent=2)+'\n')
    return project

def macro_auc(frame):
    vals=[]
    for _,g in frame.groupby('Genus'):
        p=g.loc[g.Ecology=='plaque','Y_pred'].to_numpy()
        m=g.loc[g.Ecology=='mucosal','Y_pred'].to_numpy()
        vals.append(float(((p[:,None]>m).mean()+.5*(p[:,None]==m).mean())))
    return float(np.mean(vals))

def select_candidates(project):
    d=project/'study13/10_additional_analyses_20261012'
    raw=pd.read_csv(d/'crossgenus_allfeatures_raw_permutation_results.tsv',sep='\t').set_index('feature')
    adj=pd.read_csv(d/'crossgenus_allfeatures_genome_adjusted_permutation_results.tsv',sep='\t').set_index('feature')
    adj=adj.loc[raw.index]
    direction=np.where((raw.effect_genus_mean>0)&(adj.effect_genus_mean>0),'plaque',np.where((raw.effect_genus_mean<0)&(adj.effect_genus_mean<0),'mucosal','discordant'))
    rows=[]
    for f,di in zip(raw.index,direction):
        r=raw.loc[f];a=adj.loc[f]
        suffix='plaque' if di=='plaque' else 'mucosal'
        gs=min(r['n_'+suffix+'_genus_material'],a['n_'+suffix+'_genus_material'])
        ps=min(r['n_phyla_'+('plaque' if di=='plaque' else 'mucosa')],a['n_phyla_'+('plaque' if di=='plaque' else 'mucosa')])
        if di!='discordant' and gs>=5 and ps>=3 and r.perm_p_consistency<.05 and a.perm_p_consistency<.05:
            rows.append(dict(feature=f,direction=di,genus_support=gs,phyla_support=ps,raw_perm_p_consistency=r.perm_p_consistency,adj_perm_p_consistency=a.perm_p_consistency,raw_perm_p_t=r.perm_p_t,adj_perm_p_t=a.perm_p_t))
    result=pd.DataFrame(rows);result.to_csv(d/'reproduced_14_candidates.tsv',sep='\t',index=False)
    return result

def compare(project,stages):
    expected=REPO/'data/frozen/expected';e=project/'study13';checks=[]
    if 96 in stages:
        obs=select_candidates(project).set_index('feature')
        exp=pd.read_csv(expected/'Figure2_B_14focal_statistics.tsv',sep='\t').set_index('feature')
        same=set(obs.index)==set(exp.index)
        diffs={c:float(np.max(np.abs(obs.loc[exp.index,c]-exp[c]))) for c in ['raw_perm_p_consistency','adj_perm_p_consistency','raw_perm_p_t','adj_perm_p_t']} if same else {}
        folder=e/'10_additional_analyses_20261012'
        raw=pd.read_csv(folder/'crossgenus_allfeatures_raw_permutation_results.tsv',sep='\t').set_index('feature')
        adj=pd.read_csv(folder/'crossgenus_allfeatures_genome_adjusted_permutation_results.tsv',sep='\t').set_index('feature').loc[raw.index]
        counts=dict(n_profiles=len(raw),raw_nominal=int((raw.perm_p_consistency<.05).sum()),both_nominal=int(((raw.perm_p_consistency<.05)&(adj.perm_p_consistency<.05)).sum()),raw_FWER_hits=int((raw.maxT_all_feature_FWER_p_consistency<.05).sum()),adjusted_FWER_hits=int((adj.maxT_all_feature_FWER_p_consistency<.05).sum()))
        expected_counts=dict(n_profiles=13917,raw_nominal=119,both_nominal=23,raw_FWER_hits=0,adjusted_FWER_hits=0)
        checks.append(dict(check='4096-permutation genome-wide screen',candidate_set_matches=same,n_candidates=len(obs),screen_counts=counts,max_P_differences=diffs,passed=same and counts==expected_counts and max(diffs.values())<1e-12))
    for stage,prefix in [(104,'nested_genusbalanced_'),(111,'nested_genusbalanced_covariate_residualized_')]:
        if stage not in stages:continue
        for label in ['original','sourcebalanced']:
            name=prefix+label+('_HMT_predictions.tsv' if stage==104 else '.tsv')
            obs=pd.read_csv(e/'11_manuscript_completion_20261013'/name,sep='\t').sort_values('HMT')
            exp=pd.read_csv(expected/name,sep='\t').sort_values('HMT')
            assert obs.HMT.tolist()==exp.HMT.tolist()
            delta=float(np.max(np.abs(obs.Y_pred.to_numpy()-exp.Y_pred.to_numpy())))
            auc=macro_auc(obs);target=macro_auc(exp)
            prediction_atol = 1e-4
            checks.append(dict(check=name,macro_auc=auc,expected_auc=target,max_prediction_difference=delta,prediction_absolute_tolerance=prediction_atol,passed=abs(auc-target) < 1e-12 and delta <= prediction_atol,))
    return checks

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=REPO/'build/reproduction')
    ap.add_argument('--prepare-only',action='store_true',help='Materialise all workflow templates without executing analyses')
    ap.add_argument('--stages',default='96,104,111');ap.add_argument('--screen-permutations',type=int,default=4096)
    a=ap.parse_args();stages=[int(x) for x in a.stages.split(',')]
    allowed={96:'96_crossgenus_balanced_permutation_scan.py',104:'104_nested_genusbalanced.py',111:'111_foldwise_covariate_adjusted_prediction.py'}
    if any(x not in allowed for x in stages):ap.error('Supported replay stages: 96,104,111')
    project=prepare(a.output)
    if a.prepare_only:
        print('Workflow templates prepared. See docs/PIPELINES.md for external inputs.');return
    logs=Path(a.output)/'logs';logs.mkdir(exist_ok=True)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
    runs=[]
    for stage in stages:
        path=project/'study13/scripts'/allowed[stage]
        cmd=[sys.executable,str(path)]+([str(a.screen_permutations)] if stage==96 else [])
        start=time.time()
        with (logs/(str(stage)+'.log')).open('w') as f:
            subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        runs.append(dict(stage=stage,seconds=round(time.time()-start,2),log=str(logs/(str(stage)+'.log'))));print(json.dumps(runs[-1]),flush=True)
    checks=compare(project,stages) if a.screen_permutations==4096 else compare(project,[x for x in stages if x!=96])
    report=dict(runs=runs,checks=checks,passed=all(x['passed'] for x in checks),screen_permutations=a.screen_permutations)
    (Path(a.output)/'reproduction_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit('Archived-result comparison failed; inspect reproduction_report.json')

if __name__=='__main__':main()
