"""Validate actual published records, support, assignments and DB design."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / 'results/chipseq/PRJNA312093'
EXPECTED = {'female_adult': (12652, 7031, 4348), 'female_cercaria': (12946, 6933, 4461),
            'male_adult': (14088, 7128, 4719), 'male_cercaria': (14151, 7103, 4651)}

def table(path):
    with gzip.open(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def rank(matrix):
    rows = [list(map(float, row)) for row in matrix]
    pivot = 0
    for col in range(len(rows[0])):
        found = next((i for i in range(pivot, len(rows)) if abs(rows[i][col]) > 1e-10), None)
        if found is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        divisor = rows[pivot][col]
        rows[pivot] = [value / divisor for value in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot:
                factor = rows[i][col]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot

class NarrowPublicationTests(unittest.TestCase):
    def test_catalogs_and_source_assignments(self):
        groups = list((ROOT / 'narrow/chipseq/consensus').glob('*/manifest.json'))
        self.assertEqual({p.parent.name for p in groups}, set(EXPECTED))
        for path in groups:
            condition = path.parent.name
            doc = json.loads(path.read_text())
            self.assertEqual(doc['regional_min_libraries'], 1)
            self.assertEqual(doc['core_min_libraries'], 2)
            self.assertEqual(doc['helixforge_execution_commit'], '2833442eb01edbe2e090870ba945fc074bc4d564')
            for artifact in doc['artifacts'].values():
                target = ROOT / artifact['path']
                self.assertEqual(target.stat().st_size, artifact['bytes'])
                self.assertEqual(hashlib.sha256(target.read_bytes()).hexdigest(), artifact['sha256'])
            with gzip.open(path.parent / 'consolidated_peaks.bed.gz', 'rt') as handle:
                bed = [line.rstrip().split('\t') for line in handle if line.strip()]
            self.assertEqual(len(bed), EXPECTED[condition][0])
            ids = {row[3] for row in bed}
            self.assertEqual(len(ids), len(bed))
            for chrom, start, end, *_ in bed:
                self.assertTrue(chrom.startswith('SM_V10_'))
                self.assertGreater(int(end), int(start))
            original = {r['PeakID']: r for r in table(ROOT / doc['artifacts']['homer_native']['path'])}
            associated = table(ROOT / doc['artifacts']['peak_gene_associations.tsv']['path'])
            self.assertEqual({row['peak_id'] for row in associated}, ids)
            self.assertEqual(len({row['gene_id'] for row in associated}), EXPECTED[condition][1])
            counts = {'true': 0, 'false': 0, 'not_assessed': 0}
            libraries = {r['record_id'] for r in doc['library_metadata']}
            self.assertEqual(len(libraries), 2)
            for row in associated:
                raw = original[row['peak_id']]
                self.assertEqual(row['gene_id'], raw['Nearest PromoterID'])
                self.assertEqual(row['homer_annotation'], raw['Annotation'])
                self.assertEqual(row['distance_to_tss'], raw['Distance to TSS'])
                self.assertEqual(row['occupancy_status'], 'not_assessed')
                evidence = json.loads(row['library_support_json'])
                self.assertEqual({r['record_id'] for r in evidence}, libraries)
                called = {r['record_id'] for r in evidence if int(r['called_peak_count']) > 0}
                self.assertEqual(called, set(row['calling_libraries'].split(',')))
                self.assertEqual(len(called), int(row['calling_library_count']))
                self.assertEqual(row['library_support_category'], 'single_library' if len(called) == 1 else 'all_libraries')
                self.assertTrue(row['gene_assignment_ambiguity_reason'])
                counts[row['gene_assignment_ambiguous']] += 1
                alternatives = json.loads(row['alternative_gene_ids'])
                if row['gene_assignment_ambiguous'] == 'true':
                    self.assertGreater(len(alternatives), 0)
                    self.assertNotIn(row['gene_id'], alternatives)
            self.assertEqual(counts, doc['gene_assignment_ambiguity']['counts'])

    def test_idr_is_preserved_separately(self):
        for condition, (regions, _genes, idr_count) in EXPECTED.items():
            root = ROOT / 'narrow/chipseq/idr' / condition
            stats = json.loads((root / 'statistics.json').read_text())
            self.assertEqual(stats['idr_threshold'], 0.05)
            self.assertEqual(stats['rank_metric'], 'signal_value')
            self.assertEqual(stats['consolidated_peaks'], idr_count)
            with gzip.open(root / 'consolidated_peaks.bed.gz', 'rt') as handle:
                self.assertEqual(sum(bool(line.strip()) for line in handle), idr_count)
            self.assertLess(idr_count, regions)

    def test_real_qc_and_no_invented_db(self):
        qc = table(ROOT / 'narrow/qc/peak_qc_summary.tsv.gz')
        self.assertEqual(len(qc), 8)
        self.assertEqual(len({r['record_id'] for r in qc}), 8)
        self.assertTrue(all(r['unit'] == 'reads' for r in qc))
        for row in qc:
            self.assertAlmostEqual(float(row['frip']), int(row['units_in_peaks']) / int(row['total_units']))
        tracks = table(ROOT / 'narrow/tracks/tracks.tsv.gz')
        self.assertEqual(len(tracks), 20)
        self.assertEqual(sum(r['track_role'] == 'individual' for r in tracks), 16)
        self.assertEqual(sum(r['track_role'] == 'aggregate' for r in tracks), 4)
        execution = json.loads((ROOT / 'narrow/manifests/native_execution_summary.json').read_text())
        self.assertEqual(sum(execution['completed_tasks'].values()), 226)
        self.assertTrue(execution['idr_is_separate_evidence'])
        self.assertFalse(execution['public_copy_is_self_contained_reentry_package'])
        self.assertEqual(execution['differential_binding'], 'HOLD_NOT_EXECUTED')

    def test_db_decision_matches_frozen_design(self):
        review = json.loads((REPO / 'config/chipseq/PRJNA312093/db_review.json').read_text())
        self.assertEqual(review['status'], 'HOLD')
        self.assertFalse(review['differential_binding_executed'])
        self.assertEqual(review['pairing_decision'], 'unpaired')
        with (REPO / 'metadata/chipseq/PRJNA312093/design.tsv').open() as handle:
            design = [r for r in csv.DictReader(handle, delimiter='\t') if r['is_control'] == 'false']
        with (REPO / 'provenance/source_package/frozen/chipseq/PRJNA312093/ena_read_run.tsv').open() as handle:
            instrument = {r['run_accession']: r['instrument_model'] for r in csv.DictReader(handle, delimiter='\t')}
        for comparison in review['comparisons']:
            selected = [r for r in design if r['stage'] == comparison['stage'] and r['target'] == comparison['mark']]
            self.assertEqual(len(selected), 4)
            matrix = [[1, int(r['sex'] == 'male'), int(instrument[r['run_accession']] == 'Illumina HiSeq 2000')] for r in selected]
            if comparison['stage'] == 'adult':
                self.assertEqual(rank(matrix), 2)
                self.assertIn('sex_and_sequencing_instrument_completely_confounded', comparison['reasons'])
            if comparison['mark'] == 'H3K4me3' and comparison['stage'] == 'cercaria':
                self.assertEqual({instrument[r['run_accession']] for r in selected}, {'Illumina HiSeq 2500'})
            self.assertEqual(comparison['status'], 'HOLD')
        self.assertEqual(review['retention']['ip_bam_bai'], 'PROTECTED')
        self.assertEqual(review['retention']['input_bam_bai'], 'PROTECTED')

if __name__ == '__main__':
    unittest.main()
