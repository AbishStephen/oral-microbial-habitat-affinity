"""Scientific regression and integrity checks against frozen study outputs."""
from pathlib import Path
import csv, hashlib, json, math, os, subprocess, sys, tempfile, unittest
import numpy as np
import pandas as pd
import ast, gzip, re
from homd_signatures.neighborhoods import read_tsv, choose_loci, extract, PFAMS
from homd_signatures.reproduce import REPO, prepare

class StudyIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=REPO/'data/frozen';cls.primary=cls.data/'project/study13'

    def test_frozen_input_checksums(self):
        for r in json.loads((self.data/'input_checksums.json').read_text()):
            p=self.data/r['path'];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r['sha256'],r['path'])

    def test_release_script_checksums_and_all_python_shell_syntax(self):
        rows=read_tsv(REPO/'workflows/workflow_catalog.tsv')
        self.assertEqual(len(rows),133)
        for r in rows:
            p=REPO/r['script'];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),r['sha256'],r['script'])
            if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
            elif p.suffix in ['.sh','.sbatch']:subprocess.run(['bash','-n',str(p)],check=True,capture_output=True)
        for p in (REPO/'scripts').rglob('*'):
            if p.suffix in ['.sh','.sbatch']:subprocess.run(['bash','-n',str(p)],check=True,capture_output=True)

    def test_release_contains_no_local_storage_paths(self):
        # Check distributable source/data, including compressed coordinate tables.
        forbidden=re.compile(r"/(?:gpfs|workspace|home|Users|mnt/data|data/home|share/apps)/|[A-Za-z]:\\Users\\|[A-Za-z0-9_.+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
        for top in ['src','scripts','analysis','workflows','docs','data','figure_sources','r']:
            for p in (REPO/top).rglob('*'):
                if not p.is_file() or '__pycache__' in p.parts:continue
                try:
                    text=gzip.open(p,'rt').read() if p.suffix=='.gz' else p.read_text()
                except UnicodeDecodeError:continue
                self.assertIsNone(forbidden.search(text),str(p.relative_to(REPO)))

    def test_primary_population_and_genome_nesting(self):
        h=pd.read_csv(self.primary/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t')
        g=pd.read_csv(self.primary/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv',sep='\t')
        self.assertEqual((len(h),h.Genus.nunique(),len(g),g.GCA.nunique()),(47,13,1121,1121))
        self.assertEqual(h.Ecology.value_counts().to_dict(),{'plaque':24,'mucosal':23})
        self.assertTrue(all(set(v.Ecology)=={'plaque','mucosal'} for _,v in h.groupby('Genus')))
        self.assertEqual(g.groupby('HMT').size().sort_index().tolist(),h.set_index('HMT').sort_index().n_genomes.tolist())
        self.assertTrue((g.CheckM2_complete>=95).all() and (g.CheckM2_contam<=5).all())

    def test_annotated_profiles_and_exploratory_statistics(self):
        x=pd.read_csv(self.primary/'handoff/X_primary_transferable_DEDUP.tsv',sep='\t')
        self.assertEqual(x.shape,(47,13922))
        f=pd.read_csv(self.data/'expected/Figure2_B_14focal_statistics.tsv',sep='\t')
        self.assertEqual(len(f),14)
        for k in ['raw','adj']:
            self.assertTrue((f[k+'_perm_p_consistency']<.05).all())
            self.assertTrue((f[k+'_maxT_all_feature_FWER_p_consistency']>=.05).all())
        self.assertEqual(sum(f.feature.str.startswith('PFAM:')),10)

    def test_primary_manifest_rebuild_with_documented_exclusions(self):
        with tempfile.TemporaryDirectory() as d:
            subprocess.run([sys.executable,str(REPO/'analysis/01_study_design/rebuild_primary_manifest.py'),'--metadata',str(self.data/'project/00_metadata'),'--exclusions',str(REPO/'data/genome_exclusions.tsv'),'--output',d],check=True,capture_output=True)
            for suffix,key in [('GENOME','GCA'),('HMT','HMT')]:
                a=pd.read_csv(Path(d)/f'PRIMARY_{suffix}_MANIFEST.tsv',sep='\t')
                b=pd.read_csv(self.primary/f'01_genome_manifests/EXPANDED_{suffix}_MANIFEST.tsv',sep='\t')
                pd.testing.assert_frame_equal(a,b,check_exact=False,atol=1e-12,rtol=0)

    def test_all_ten_domains_reextract_from_genome_coordinate_records(self):
        d=self.data/'neighborhoods'
        selected=choose_loci(read_tsv(d/'candidate_coordinates.tsv.gz'),read_tsv(d/'candidate_immediate_neighbors.tsv.gz'),read_tsv(self.primary/'01_genome_manifests/EXPANDED_GENOME_MANIFEST.tsv'))
        loci,genes=extract(selected,read_tsv(d/'selected_genome_CDS.tsv.gz'))
        self.assertEqual({r['Pfam'] for r in loci},set(PFAMS));self.assertEqual(len(loci),30)
        self.assertEqual(len({r['locus_id'] for r in loci}),30)
        expected=read_tsv(d/'locus_genes.tsv')
        key=lambda r:(r['locus_id'],r['protein_id'])
        a={key(r):r for r in genes};b={key(r):r for r in expected};self.assertEqual(set(a),set(b))
        for k in a:
            for field in ['start','end','strand','display_start','display_end','display_strand','focal']:
                self.assertEqual(str(a[k][field]),str(b[k][field]),(k,field))

    def test_reverse_complement_preserves_gene_length_and_focal_direction(self):
        d=self.data/'neighborhoods';rows=read_tsv(d/'locus_genes.tsv')
        for r in rows:
            self.assertEqual(int(r['end'])-int(r['start'])+1,int(r['display_end'])-int(r['display_start']))
            if r['focal']=='True':
                self.assertEqual(r['display_strand'],'+');self.assertEqual(int(r['display_start']),0)
        self.assertTrue(any(r['strand']=='-' and r['focal']=='True' for r in rows))

    def test_full_refitting_nulls_have_expected_size_and_reproduce_P(self):
        e=self.data/'expected'
        for dist,summary,n in [('nested_genusbalanced_200_full_refit_null_distribution.tsv','nested_genusbalanced_200_full_refit_permutation_results.tsv',200),('covadjusted_100_withingenus_perm_null_distribution.tsv','covadjusted_100_withingenus_perm_null_results.tsv',100)]:
            d=pd.read_csv(e/dist,sep='\t');s=pd.read_csv(e/summary,sep='\t')
            self.assertEqual(len(d),2*n);self.assertFalse(d.duplicated(['perm_id','label']).any())
            for r in s.to_dict('records'):
                label=r['outcome'];v=d[d.label==label];self.assertEqual(v.perm_id.nunique(),n)
                observed=r.get('observed_macro_auc',r.get('observed_auc'))
                p=(1+int((v.macro_auc>=observed).sum()))/(n+1)
                target=r.get('full_refit_withingenus_perm_p',r.get('p_full_nested_refitting'))
                self.assertAlmostEqual(p,target,places=14)

    def test_actual_archived_ranking_and_adjustment_exclude_heldout_taxa(self):
        with tempfile.TemporaryDirectory() as d:
            p=prepare(Path(d));script=p/'study13/scripts/111_foldwise_covariate_adjusted_prediction.py'
            lines=script.read_text().splitlines();cut=next(i for i,line in enumerate(lines) if line.startswith('for label,y in '));prefix='\n'.join(lines[:cut])+'\n'
            ns={'__file__':str(script)};exec(compile(prefix,str(script),'exec'),ns)
            test=np.flatnonzero(ns['gen']==ns['gs'][0]);train=np.flatnonzero(ns['gen']!=ns['gs'][0])
            before=ns['rank'](train);adjusted=ns['adjusted_X'](train)[train].copy()
            ns['X'][test]*=1000;ns['C'][test]+=1000
            np.testing.assert_array_equal(before,ns['rank'](train))
            np.testing.assert_allclose(adjusted,ns['adjusted_X'](train)[train],rtol=0,atol=1e-7)

    def test_recovered_plot_inputs_match_full_precision_primary_data(self):
        original=pd.read_csv(self.primary/'01_genome_manifests/EXPANDED_HMT_MANIFEST.tsv',sep='\t').set_index('HMT')
        plotted=pd.read_csv(REPO/'figure_sources/Figure1/HOMD47_HMT_manifest.tsv',sep='\t').set_index('HMT')
        self.assertEqual(set(original.index),set(plotted.index))
        np.testing.assert_allclose(plotted.loc[original.index,'Median_Y'],original.Median_Y,atol=1e-12,rtol=0)
        self.assertEqual(plotted.n_genomes.sum(),1121)
        d=json.loads((REPO/'figure_sources/Figure5/Figure5_model_predictions_and_sister_pairs.json').read_text())
        plotted=pd.DataFrame(d['pred']).set_index('HMT')
        expected=pd.read_csv(self.data/'expected/nested_genusbalanced_sourcebalanced_HMT_predictions.tsv',sep='\t').set_index('HMT')
        np.testing.assert_allclose(plotted.loc[expected.index,'Y_pred'],expected.Y_pred,atol=1e-12,rtol=0)

if __name__=='__main__':unittest.main(verbosity=2)
