"""Guard the Methods-based layout and relocated runtime imports."""
import csv
from pathlib import Path
import tempfile
import unittest
from homd_signatures.reproduce import REPO, prepare

class OrganisationIntegrity(unittest.TestCase):
    def test_catalog_covers_retained_sources_once(self):
        with (REPO/'workflows/workflow_catalog.tsv').open() as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        sources={r['script'] for r in rows}
        actual={str(p.relative_to(REPO)) for p in (REPO/'analysis').rglob('*') if p.suffix in ['.py','.sh','.sbatch']}
        actual.remove('analysis/01_study_design/rebuild_primary_manifest.py')
        self.assertEqual(sources,actual)
        self.assertEqual(len(rows),len({r['runtime_path'] for r in rows}))
        self.assertEqual(len(list((REPO/'analysis').glob('*/README.md'))),11)
        self.assertFalse({72,73,74,75,76,80,81,82}&{int(r['stage']) for r in rows})

    def test_prepared_sources_preserve_scientific_code(self):
        with tempfile.TemporaryDirectory(prefix='homd relocated ') as folder:
            output=Path(folder);project=prepare(output)
            with (REPO/'workflows/workflow_catalog.tsv').open() as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    original=(REPO/row['script']).read_text()
                    expected=original.replace('__PROJECT_ROOT__',str(project)).replace('__WORK_ROOT__',str(output/'scratch')).replace('__HOME_ROOT__',str(output/'home'))
                    self.assertEqual((project/row['runtime_path']).read_text(),expected,row['script'])

if __name__=='__main__':unittest.main()
