"""Verify scientific payloads, layer separation and publication privacy."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import unittest

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/'results/chipseq/PRJNA602708'
EXPECTED={
    ('mapq30','narrow'):{'hm7w':4330,'rh0w':8080,'rh8w':1472,'rh10w':1722,'rh12w':7001},
    ('mapq0','narrow'):{'hm7w':2221,'rh0w':5278,'rh8w':1385,'rh10w':881,'rh12w':3653},
    ('mapq30','broad'):{'hm7w':75351,'rh0w':66549,'rh8w':73876,'rh10w':57827,'rh12w':65094},
    ('mapq0','broad'):{'hm7w':93100,'rh0w':84539,'rh8w':95945,'rh10w':52613,'rh12w':48907},
}

def rows(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


class PublishedPlasmaTests(unittest.TestCase):
    def test_twenty_native_catalogs_and_reference_identity(self):
        found=list(RESULT.glob('mapq*/*/chipseq/consensus/*/manifest.json'))
        self.assertEqual(len(found),20)
        for (policy,arm),groups in EXPECTED.items():
            for condition,n in groups.items():
                base=RESULT/policy/arm/'chipseq/consensus'/condition
                doc=json.loads((base/'manifest.json').read_text())
                self.assertEqual(doc['mapping_policy'],policy)
                self.assertEqual(doc['peak_type'],arm)
                self.assertEqual(doc['condition'],condition)
                self.assertEqual(doc['metrics']['regions'],n)
                self.assertEqual(doc['regional_min_libraries'],1)
                self.assertEqual(doc['core_min_libraries'],2)
                self.assertEqual(doc['occupancy_status'],'not_assessed')
                self.assertEqual(doc['reference_id'],'Schistosoma_mansoni_SM_V10_WBPS19')
                self.assertFalse(doc['biological_assignments_changed'])
                for artifact in doc['artifacts'].values():
                    p=RESULT/artifact['path']
                    self.assertEqual(p.stat().st_size,artifact['bytes'])
                    self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),artifact['sha256'])
                    with gzip.open(p,'rb') as f:
                        self.assertEqual(hashlib.sha256(f.read()).hexdigest(),artifact['source_sha256'])
                with gzip.open(base/'consolidated_peaks.bed.gz','rt') as f:
                    bed=[l.split() for l in f if l.strip()]
                self.assertEqual(len(bed),n)
                self.assertEqual(len({r[3] for r in bed}),n)
                self.assertTrue(all(int(r[1])<int(r[2]) for r in bed))
                self.assertEqual(sum(int(r[2])-int(r[1]) for r in bed),doc['metrics']['region_bp'])

    def test_library_membership_and_gene_associations_are_retained(self):
        for manifest in RESULT.glob('mapq*/*/chipseq/consensus/*/manifest.json'):
            doc=json.loads(manifest.read_text()); condition=manifest.parent.name
            ids={r['record_id'] for r in doc['libraries']}
            self.assertEqual(len(ids),2)
            cat={r['peak_id']:r for r in rows(manifest.parent/'consolidated_peaks.tsv.gz')}
            annotation=RESULT/doc['mapping_policy']/doc['peak_type']/'chipseq/peak_annotation'/condition
            seen=set()
            for row in rows(annotation/'annotated_peaks.tsv.gz'):
                seen.add(row['peak_id']); original=cat[row['peak_id']]
                self.assertEqual(row['calling_libraries'],original['calling_libraries'])
                evidence=json.loads(row['library_support_json'])
                self.assertEqual({e['record_id'] for e in evidence},ids)
                self.assertEqual(len(evidence),2)
                calling={e['record_id'] for e in evidence if int(e['called_peak_count'])>0}
                self.assertEqual(calling,set(row['calling_libraries'].split(',')))
                self.assertEqual(len(calling),int(row['calling_library_count']))
                self.assertEqual(int(row['total_library_count']),2)
                self.assertEqual(row['library_support_category'],'all_libraries' if len(calling)==2 else 'single_library')
                self.assertIn(row['gene_assignment_ambiguous'],{'true','false','not_assessed'})
                self.assertEqual(row['occupancy_status'],'not_assessed')
            self.assertEqual(seen,set(cat))
            self.assertTrue({r['peak_id'] for r in rows(annotation/'peak_gene_associations.tsv.gz')}<=seen)

    def test_qc_and_reports_describe_the_four_actual_layers(self):
        qc=rows(RESULT/'qc/qc_sample_summary.tsv')
        self.assertEqual(len(qc),48)
        self.assertEqual(len({r['run_accession'] for r in qc}),22)
        expected_counts={('mapq30','narrow'):22672,('mapq30','broad'):377974,
                         ('mapq0','narrow'):14532,('mapq0','broad'):429300}
        for key,count in expected_counts.items():
            selected=[r for r in qc if (r['mapping_policy'],r['peak_type'])==key]
            self.assertEqual(len(selected),12)
            ips=[r for r in selected if r['is_control']=='false']
            self.assertEqual(len(ips),10)
            self.assertEqual(sum(int(r['peak_count']) for r in ips),count)
            for r in ips:
                self.assertTrue(0<=float(r['frip'])<=1)
                self.assertGreater(int(r['final_bam_reads']),0)
            report=RESULT/key[0]/key[1]/'chipseq/report/chipseq_report.html'
            self.assertGreater(report.stat().st_size,1000)
        summary=json.loads((RESULT/'qc/qc_summary.json').read_text())
        self.assertEqual(summary['regional_groups'],20)
        self.assertEqual(summary['differential_binding'],'SUSPENDED')
        self.assertEqual(summary['idr_new_layers'],'NOT_EXECUTED')
        self.assertEqual(len(rows(RESULT/'tracks/track_inventory.tsv')),48)

    def test_checksums_privacy_and_no_heavy_binaries(self):
        forbidden=re.compile(r'/home/|/scratch/|/mnt/[a-z]/|[A-Z]:[/\\]Users[/\\]|ra\d{6}|srv-slurm|[\w.-]+@[\w.-]+\.(?:br|edu)',re.I)
        self.assertTrue((RESULT/'SHA256SUMS').is_file())
        for line in (RESULT/'SHA256SUMS').read_text().splitlines():
            digest,name=line.split('  ',1)
            self.assertEqual(hashlib.sha256((RESULT/name).read_bytes()).hexdigest(),digest,name)
        for base in [RESULT,ROOT/'provenance/chipseq/PRJNA602708']:
            for p in base.rglob('*'):
                if not p.is_file(): continue
                self.assertNotIn(p.suffix.lower(),{'.bam','.bai','.bw','.bigwig','.fastq','.fq'})
                opener=gzip.open if p.suffix=='.gz' else open
                with opener(p,'rt',encoding='utf-8') as f:
                    for line in f: self.assertIsNone(forbidden.search(line),p.name)

    def test_layer_policy_does_not_claim_truth_or_new_idr(self):
        p=json.loads((ROOT/'config/chipseq/PRJNA602708/analysis_policy.json').read_text())
        self.assertEqual(p['public_delivery_status'],'PUBLISHED_WITH_LIMITATIONS')
        self.assertTrue(p['layers_are_not_combined'])
        self.assertFalse(p['independent_occupancy_truth_validated'])
        self.assertFalse(p['exact_publication_reproduction_claimed'])
        self.assertEqual(p['idr_for_reanalysis_layers'],'NOT_EXECUTED')


if __name__=='__main__': unittest.main()
