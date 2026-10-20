"""Check published scientific results, not merely declared success flags."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import unittest

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/'results/chipseq/PRJNA1090249'

def rows(path):
    with path.open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f,delimiter='\t'))

class PublishedHp1Tests(unittest.TestCase):
    def test_qc_identity_and_observed_peak_counts_match_native_files(self):
        qc=rows(RESULT/'qc/qc_sample_summary.tsv')
        metadata=rows(ROOT/'metadata/chipseq/PRJNA1090249/helixforge.broad.tsv')
        self.assertEqual({r['run_accession'] for r in qc},{r['run_accession'] for r in metadata})
        self.assertEqual(len(qc),8)
        ips=[r for r in qc if r['is_control']=='false']
        self.assertEqual(len(ips),6)
        for q in ips:
            paths=list((RESULT/'broad/080-peak-calling').glob(f"{q['run_accession']}.*.peak_calling/peaks.broadPeak"))
            self.assertEqual(len(paths),1)
            self.assertEqual(sum(1 for line in paths[0].read_text().splitlines() if line.strip()),int(q['peak_count']))
            self.assertGreater(int(q['final_bam_reads']),0)
            self.assertTrue(0<=float(q['frip'])<=1)

    def test_consensus_counts_are_supported_by_retained_intervals(self):
        stats=list((RESULT/'broad/chipseq/consensus').glob('*/statistics.json'))
        self.assertEqual(len(stats),2)
        expected={'cercaria':1816,'sporocyst':867}
        for path in stats:
            stat=json.loads(path.read_text())
            condition=path.parent.name
            with gzip.open(path.parent/'consolidated_peaks.bed.gz','rt') as handle:
                intervals=[line for line in handle if line.strip()]
            self.assertEqual(len(intervals),expected[condition])
            self.assertEqual(stat['consolidated_peaks'],len(intervals))
            self.assertEqual(stat['regional_min_libraries'],1)
            self.assertEqual(stat['core_min_libraries'],2)

    def test_public_checksums_and_json_are_valid(self):
        for path in RESULT.rglob('*.json'):
            json.loads(path.read_text())
        for line in (RESULT/'SHA256SUMS').read_text().splitlines():
            digest,name=line.split('  ',1)
            self.assertEqual(hashlib.sha256((RESULT/name).read_bytes()).hexdigest(),digest,name)

    def test_delivery_inventory_tracks_the_sanitized_public_copy(self):
        manifest=json.loads((RESULT/'manifests/publication_validation.json').read_text())
        self.assertGreater(len(manifest['files']),100)
        for item in manifest['files']:
            group,name=item['path'].split('/',1)
            base=RESULT if group=='results' else ROOT/'provenance/chipseq/PRJNA1090249'
            path=base/name
            self.assertEqual(path.stat().st_size,item['bytes'],item['path'])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),item['sha256'],item['path'])

    def test_public_text_has_no_private_absolute_paths(self):
        forbidden=re.compile(r'/home/|/scratch/|/mnt/[a-z]/|[A-Z]:\\\\?Users\\',re.I)
        for base in (RESULT,ROOT/'provenance/chipseq/PRJNA1090249'):
            for path in base.rglob('*'):
                if path.is_file() and path.suffix not in {'.png','.pdf'}:
                    if path.suffix=='.gz':
                        with gzip.open(path,'rt') as handle:
                            for line in handle:
                                self.assertIsNone(forbidden.search(line),path.name)
                    else:
                        self.assertIsNone(forbidden.search(path.read_text()),path.name)

    def test_library_evidence_survives_gene_association(self):
        for condition,total in [('cercaria',3),('sporocyst',3)]:
            root=RESULT/'broad/chipseq'
            with gzip.open(root/f'consensus/{condition}/consolidated_peaks.tsv.gz','rt') as f:
                catalog={r['peak_id']:r for r in csv.DictReader(f,delimiter='\t')}
            observed=set()
            with gzip.open(root/f'peak_annotation/{condition}/peak_gene_associations.tsv.gz','rt') as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    observed.add(row['peak_id'])
                    original=catalog[row['peak_id']]
                    evidence=json.loads(row['library_support_json'])
                    self.assertEqual(len(evidence),total)
                    calling={r['record_id'] for r in evidence if int(r['called_peak_count'])>0}
                    self.assertEqual(calling,set(row['calling_libraries'].split(',')))
                    self.assertEqual(len(calling),int(row['calling_library_count']))
                    self.assertEqual(total,int(row['total_library_count']))
                    self.assertEqual(original['calling_libraries'],row['calling_libraries'])
                    self.assertEqual(row['gene_assignment_ambiguous'],'not_assessed')
                    self.assertEqual(row['occupancy_status'],'not_assessed')
                    category='single_library' if len(calling)==1 else 'all_libraries' if len(calling)==total else 'multiple_libraries'
                    self.assertEqual(category,row['library_support_category'])
            self.assertEqual(observed,set(catalog))

    def test_superseded_outputs_are_not_active(self):
        self.assertFalse((RESULT/'broad/chipseq/report').exists())
        self.assertFalse((RESULT/'broad/chipseq/chipseq_run_manifest.json').exists())
        self.assertFalse((RESULT/'manifests/native_terminal_validation.json').exists())
        self.assertFalse((RESULT/'broad/chipseq/peak_annotation/peak_annotation_aggregate').exists())

    def test_regional_manifests_and_native_homer_assignments(self):
        for manifest in (RESULT/'broad/chipseq/consensus').glob('*/manifest.json'):
            doc=json.loads(manifest.read_text())
            self.assertEqual(doc['regional_min_libraries'],1)
            self.assertEqual(doc['core_min_libraries'],2)
            self.assertFalse(doc['biological_assignments_changed'])
            for artifact in doc['artifacts'].values():
                path=RESULT/artifact['path']
                self.assertEqual(path.stat().st_size,artifact['bytes'])
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),artifact['sha256'])
            with gzip.open(RESULT/doc['artifacts']['homer_native']['path'],'rt') as f:
                native={r['PeakID']:r for r in csv.DictReader(f,delimiter='\t')}
            with gzip.open(RESULT/doc['artifacts']['annotated_peaks.tsv']['path'],'rt') as f:
                rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual({r['peak_id'] for r in rows},set(native))
            for row in rows:
                raw=native[row['peak_id']]
                self.assertEqual(row['gene_ids'],raw['Nearest PromoterID'])
                self.assertEqual(row['distance_to_tss'],raw['Distance to TSS'])
                self.assertEqual(row['homer_annotation'],raw['Annotation'])

    def test_accepted_atlas_inclusion_does_not_enable_inferential_db(self):
        plan=json.loads((ROOT/'config/chipseq/PRJNA1090249/db_plan.json').read_text())
        qc=json.loads((RESULT/'qc/qc_summary.json').read_text())
        self.assertEqual(plan['final_review']['atlas_inclusion'],qc['atlas_inclusion'])
        self.assertEqual(qc['atlas_inclusion'],'ACCEPTED')
        self.assertEqual(qc['differential_binding'],'NOT_ELIGIBLE_IN_CURRENT_ANALYSIS')
        self.assertFalse(plan['execute_if_eligible'])
        self.assertEqual(qc['stage_strain_confounding'],'COMPLETE')
        self.assertFalse(qc['fragment_length_experimentally_certified'])
        validation=json.loads((RESULT/'manifests/publication_validation.json').read_text())
        self.assertFalse(validation['public_copy_is_self_contained_reentry_package'])

if __name__=='__main__':unittest.main()
