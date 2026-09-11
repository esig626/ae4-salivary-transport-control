"""Post-freeze completeness, chronological firewall and baseline reuse."""
from datetime import datetime
import csv
import unittest

from modern_full_model import task17_fixed_wt as task
from modern_full_model.task17_postfreeze import checked_checkpoint
from modern_full_model.validation import sha256_file


class PostFreezeTests(unittest.TestCase):
    def test_pushed_checkpoint_precedes_genotype_access(self):
        freeze,receipt=checked_checkpoint()
        result=task.read(task.OUT/'genotype_results_frozen.json')
        self.assertEqual(result['pushed_wt_checkpoint'],receipt['checkpoint_sha'])
        self.assertGreater(datetime.fromisoformat(result['created_utc']),
            datetime.fromisoformat(receipt['remote_ref_verified_utc'].replace('Z','+00:00')))
        self.assertFalse(freeze['genotype_outcomes_used'])

    def test_every_root_condition_calcium_and_no_ineligible_predictions(self):
        freeze,_=checked_checkpoint()
        expected={(r['candidate_id'],ca) for r in freeze['all_decisions'] for ca in (.10,.25,.50)}
        for name,ratio in (('ae4_5pct_results.csv','R_AE4_5pct'),('ae2_results.csv','R_AE2')):
            with (task.OUT/name).open() as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),90)
            self.assertEqual({(r['candidate_id'],float(r['calcium_uM'])) for r in rows},expected)
            self.assertEqual(sum(r['status']=='BASELINE_HASH_VALID_REUSE' for r in rows),30)
            for r in rows:
                self.assertEqual(r['new_genotype_simulation'],'False')
                if r['condition']=='inherited_baseline':
                    self.assertGreater(float(r[ratio]),0)
                    self.assertEqual(r['paired_numerical_gate_pass'],'True')
                    self.assertEqual(sha256_file(task.REPO/r['wt_source_path']),r['wt_source_sha256'])
                else:
                    self.assertEqual(r[ratio],'')
                    self.assertEqual(r['status'],'NOT_EVALUATED_FIXED_WT_INFEASIBLE')

    def test_frozen_genotype_results_and_no_new_simulations(self):
        result=task.read(task.OUT/'genotype_results_frozen.json')
        for name,digest in result['hashes'].items():
            self.assertEqual(sha256_file(task.REPO/name),digest)
        for key in ('new_genotype_simulations','new_WT_state_solves',
                    'new_genotype_state_solves','exact_AE4_zero_attempts'):
            self.assertEqual(result[key],0)
        self.assertEqual(result['baseline_pairs_reused'],30)
        self.assertEqual(result['modified_combinations_not_evaluated'],60)


if __name__=='__main__':
    unittest.main()
