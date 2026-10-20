"""Check published H3K27me3 data independently of declared success flags."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import unittest

ROOT=Path(__file__).resolve().parents[1]/'results/chipseq/PRJNA312093'
EXPECTED={'female_adult':(75155,6120),'female_cercaria':(29089,4326),
          'male_adult':(39845,5031),'male_cercaria':(44179,5455)}

class PublishedH3k27Tests(unittest.TestCase):
    def test_catalogs_preserve_intervals_and_library_evidence(self):
        groups=list((ROOT/'broad/chipseq/consensus').glob('*/manifest.json'))
        self.assertEqual({p.parent.name for p in groups},set(EXPECTED))
        for path in groups:
            group=path.parent.name
            doc=json.loads(path.read_text())
            self.assertEqual(doc['regional_min_libraries'],1)
            self.assertEqual(doc['core_min_libraries'],2)
            self.assertEqual(doc['target'],'H3K27me3')
            self.assertFalse(doc['biological_assignments_changed'])
            for a in doc['artifacts'].values():
                f=ROOT/a['path']
                self.assertEqual(f.stat().st_size,a['bytes'])
                self.assertEqual(hashlib.sha256(f.read_bytes()).hexdigest(),a['sha256'])
            with gzip.open(path.parent/'consolidated_peaks.bed.gz','rt') as f:
                bed=[line.rstrip().split('\t') for line in f if line.strip()]
            self.assertEqual(len(bed),EXPECTED[group][0])
            with gzip.open(path.parent/'unfiltered_peaks.bed.gz','rt') as f:
                self.assertEqual(bed,[line.rstrip().split('\t') for line in f if line.strip()])
            ids={r[3] for r in bed}
            self.assertEqual(len(ids),len(bed))
            libraries={r['record_id'] for r in doc['library_metadata']}
            self.assertEqual(len(libraries),2)
            for r in doc['library_metadata']:
                self.assertTrue(r['control_record_id'])
                self.assertNotIn(r['control_record_id'],libraries)
            for chrom,start,end,*_ in bed:
                self.assertTrue(chrom.startswith('SM_V10_'))
                self.assertGreater(int(end),int(start))
            with gzip.open(ROOT/doc['artifacts']['homer_native']['path'],'rt') as f:
                native={r['PeakID']:r for r in csv.DictReader(f,delimiter='\t')}
            self.assertEqual(set(native),ids)
            genes=set();observed=set()
            with gzip.open(ROOT/doc['artifacts']['peak_gene_associations.tsv']['path'],'rt') as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    observed.add(r['peak_id']);genes.add(r['gene_id'])
                    raw=native[r['peak_id']]
                    self.assertEqual(r['gene_id'],raw['Nearest PromoterID'])
                    self.assertEqual(r['distance_to_tss'],raw['Distance to TSS'])
                    self.assertEqual(r['homer_annotation'],raw['Annotation'])
                    self.assertEqual(r['gene_assignment_ambiguous'],'not_assessed')
                    self.assertEqual(r['occupancy_status'],'not_assessed')
                    evidence=json.loads(r['library_support_json'])
                    self.assertEqual({x['record_id'] for x in evidence},libraries)
                    self.assertEqual(len(evidence),2)
                    called={x['record_id'] for x in evidence if int(x['called_peak_count'])>0}
                    self.assertEqual(called,set(r['calling_libraries'].split(',')))
                    self.assertEqual(len(called),int(r['calling_library_count']))
                    self.assertEqual(int(r['total_library_count']),2)
                    self.assertEqual(r['library_support_category'],'single_library' if len(called)==1 else 'all_libraries')
                    for x in evidence:
                        self.assertTrue(0<=int(x['covered_bp'])<=int(x['region_bp']))
                        self.assertAlmostEqual(float(x['covered_fraction']),int(x['covered_bp'])/int(x['region_bp']),places=6)
            self.assertEqual(observed,ids)
            self.assertEqual(len(genes),EXPECTED[group][1])

    def test_integrity_and_privacy(self):
        listed=set()
        for line in (ROOT/'SHA256SUMS').read_text().splitlines():
            digest,name=line.split('  ',1)
            f=ROOT/name
            self.assertEqual(hashlib.sha256(f.read_bytes()).hexdigest(),digest,name)
            listed.add(name)
        actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS'}
        self.assertEqual(actual,listed)
        forbidden=re.compile(r'/home/|/scratch/|/mnt/[a-z]/|[A-Z]:[/\\]Users[/\\]',re.I)
        for name in actual:
            f=ROOT/name
            opener=gzip.open if f.suffix=='.gz' else open
            with opener(f,'rt',encoding='utf-8') as handle:
                for line in handle:
                    self.assertIsNone(forbidden.search(line),name)

    def test_partial_scope_and_comparison_failure_are_preserved(self):
        doc=json.loads((ROOT/'manifests/regional_publication.json').read_text())
        self.assertEqual(doc['pending_arms'],['H3K4me3'])
        self.assertEqual(doc['source_direct_equivalence'],'FAIL_DIRECT_EQUIVALENCE')
        self.assertFalse(doc['strict_gate_redefined'])
        self.assertFalse(doc['public_copy_is_self_contained_reentry_package'])
        self.assertEqual(sum(r['direct_changed_peaks'] for r in doc['source_comparison']),2)
        self.assertFalse((ROOT/'broad/chipseq/chipseq_run_manifest.json').exists())

if __name__=='__main__':unittest.main()
