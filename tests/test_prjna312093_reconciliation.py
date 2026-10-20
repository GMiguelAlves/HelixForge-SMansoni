"""Historical sensitivity results cannot replace the reviewed RNA-seq sources."""
import hashlib
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'provenance/rnaseq/PRJNA312093/historical/nonadult_16'

class PicardReconciliationTests(unittest.TestCase):
    def test_historical_refit_is_audit_only(self):
        doc=json.loads((ARCHIVE/'archive_manifest.json').read_text())
        self.assertEqual(doc['status'],'AUDIT_ONLY')
        self.assertTrue(doc['superseded'])
        self.assertFalse(doc['atlas_eligible'])
        current=json.loads((ROOT/doc['current_approved_refit']).read_text())
        self.assertEqual(current['sample_count'],12)
        self.assertEqual([c['stage'] for c in current['contrasts']],[1,2,3])
        self.assertFalse((ROOT/'results/PRJNA312093/reanalysis/nonadult_16').exists())
        self.assertEqual(len(list((ARCHIVE/'contrasts').glob('*.tsv'))),4)

    def test_preserved_historical_files_match_archive_manifest(self):
        doc=json.loads((ARCHIVE/'archive_manifest.json').read_text())
        self.assertEqual(doc['source_commit'],'92b5d4ff6ebd44151a2d3c2df9ae5e9cde46e161')
        self.assertEqual(len(doc['files']),12)
        actual={p.relative_to(ARCHIVE).as_posix() for p in ARCHIVE.rglob('*')
                if p.is_file() and p.name not in {'README.md','archive_manifest.json'}}
        self.assertEqual(actual,{r['path'] for r in doc['files']})
        for row in doc['files']:
            path=ARCHIVE/row['path']
            self.assertEqual(path.stat().st_size,row['bytes'])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),row['sha256'])
