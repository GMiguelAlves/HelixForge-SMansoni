import csv
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'config/chipseq'


def rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class ChipseqDbPlanningTests(unittest.TestCase):
    def test_candidates_cover_all_peak_producing_studies(self):
        registry = rows(ROOT / 'metadata/chipseq/study_registry.tsv')
        planned = rows(CONFIG / 'differential_binding_candidates.tsv')
        expected = {r['study_accession'] for r in registry if r['tier'] == 'Q1' and r['execution_scope'] == 'peaks'}
        actual = {r['study_accession'] for r in planned if r['initial_review_status'].startswith('PENDING_')}
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 5)
        self.assertEqual(len(planned), len({r['study_accession'] for r in planned}))
        for row in planned:
            self.assertTrue(row['comparison_candidate'])
            self.assertTrue(row['scope_constraint'])

    def test_restricted_studies_do_not_gain_automatic_db(self):
        policy = json.loads((CONFIG / 'differential_binding_review.json').read_text())
        restricted = {r['study_accession'] for r in rows(CONFIG / 'differential_binding_candidates.tsv')
                      if r['initial_review_status'] == 'NOT_AUTOMATICALLY_PLANNED'}
        self.assertEqual(restricted, {'PRJNA219440', 'PRJNA384372', 'PRJNA816041'})
        self.assertEqual(set(policy['no_automatic_db_studies']), restricted)

    def test_unresolved_review_protects_bams_and_metadata(self):
        policy = json.loads((CONFIG / 'differential_binding_review.json').read_text())
        self.assertIs(policy['unresolved_status_blocks_bam_cleanup'], True)
        self.assertIs(policy['inputs_are_independent_model_replicates'], False)
        self.assertEqual(set(policy['review_outcomes']), {'ELIGIBLE', 'HOLD', 'NOT_ELIGIBLE'})
        retained = set(policy['retain_until_review_and_db_resolved'])
        self.assertTrue({'final_bam_and_bai_per_biological_replicate', 'input_control_bam_and_bai',
                         'bam_producer_manifests_and_checksums', 'validated_metadata_and_control_mapping',
                         'individual_and_consolidated_peaks_and_provider_manifests'} <= retained)
        self.assertIn('estimable_design_and_contrast_with_no_unresolved_confounding', policy['review_requirements'])
        self.assertIn('count_filter_and_normalization_assumptions_review', policy['review_requirements'])
        self.assertIn('db_runtime_and_artifact_preflight', policy['review_requirements'])
        self.assertIn('documented_eligibility_decision_with_evidence_and_rationale', policy['cleanup_release_requires'])

    def test_hp1_plan_matches_six_distinct_biological_ip_samples(self):
        plan = json.loads((CONFIG / 'PRJNA1090249/db_plan.json').read_text())
        metadata = rows(ROOT / 'metadata/chipseq/PRJNA1090249/helixforge.broad.tsv')
        ips = [r for r in metadata if r['is_control'] == 'false']
        self.assertEqual(len(ips), 6)
        for condition in ('cercaria', 'sporocyst'):
            selected = [r for r in ips if r['condition'] == condition]
            with self.subTest(condition=condition):
                self.assertEqual(len(selected), plan['required_biological_replicates_per_condition'])
                self.assertEqual({r['sample_id'] for r in selected}, set(plan['expected_ip_samples'][condition]))
                self.assertEqual(len({r['biological_replicate'] for r in selected}), 3)
                self.assertEqual({r['target'] for r in selected}, {'HP1'})
                self.assertEqual({r['layout'] for r in selected}, {'single'})
        controls = [r for r in metadata if r['is_control'] == 'true']
        self.assertEqual({r['sample_id'] for r in controls}, set(plan['control_samples']))
        self.assertEqual(len(controls), 2)
        self.assertTrue({r['control_id'] for r in ips} <= {r['run_accession'] for r in controls})
        self.assertTrue({r['sample_id'] for r in ips}.isdisjoint(plan['control_samples']))

    def test_hp1_plan_is_conditional_and_does_not_change_the_initial_run(self):
        plan = json.loads((CONFIG / 'PRJNA1090249/db_plan.json').read_text())
        initial = json.loads((CONFIG / 'execution_contract.json').read_text())
        self.assertEqual(plan['status'], 'PENDING_QC_AND_MODEL_REVIEW')
        self.assertIs(plan['execute_if_eligible'], True)
        self.assertIs(plan['wait_for_other_studies'], False)
        self.assertIs(plan['executable_db_spec'], False)
        self.assertIs(plan['paired_design'], False)
        self.assertIs(plan['controls_are_model_columns'], False)
        self.assertIs(plan['loss_of_a_planned_replicate_requires_review_not_silent_reduction'], True)
        self.assertEqual(plan['contrast']['numerator'], 'cercaria')
        self.assertEqual(plan['contrast']['denominator'], 'sporocyst')
        self.assertEqual(plan['proposed_design'], '~ condition')
        self.assertEqual(plan['proposed_counting_unit'], 'reads')
        self.assertIs(initial['composition']['differential_binding'], False)
        self.assertTrue((CONFIG / initial['differential_binding_review_policy']).is_file())
        self.assertTrue((CONFIG / 'PRJNA1090249' / plan['retention_policy']).is_file())


if __name__ == '__main__':
    unittest.main()
