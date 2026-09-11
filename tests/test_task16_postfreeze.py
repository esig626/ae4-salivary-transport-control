"""Post-freeze integrity and complete ensemble accounting; no simulations."""
import copy
import csv
from datetime import datetime
import unittest
from unittest.mock import patch

from modern_full_model import task16_postfreeze as post
from modern_full_model import task16_wt_allocation as task
from modern_full_model.validation import sha256_file


class PostFreezeTests(unittest.TestCase):
    def test_pushed_freeze_precedes_results_and_remains_unchanged(self):
        freeze, receipt = post.checked_checkpoint()
        result = task.read(task.OUT/"genotype_results_frozen.json")
        self.assertLess(datetime.fromisoformat(receipt["verified_utc"]), datetime.fromisoformat(result["created_utc"]))
        self.assertEqual(result["pushed_wt_checkpoint"], receipt["checkpoint_commit"])
        self.assertEqual(result["wt_manifest_sha256"], sha256_file(task.OUT/"frozen_candidate_manifest.json"))
        self.assertEqual(len(freeze["admissible_candidates"]), 30)
        self.assertTrue(all(r["condition"] == "inherited_baseline" for r in freeze["admissible_candidates"]))

    def test_unverified_remote_checkpoint_blocks_genotype_access(self):
        original = task.read
        def changed(path):
            value = original(path)
            if str(path).endswith("pushed_wt_checkpoint.json"):
                value = copy.deepcopy(value)
                value["remote_ref_verified"] = False
            return value
        with patch.object(task, "read", side_effect=changed):
            with self.assertRaises(AssertionError): post.checked_checkpoint()

    def test_all_ninety_combinations_and_unavailable_values_are_explicit(self):
        result = task.read(task.OUT/"genotype_results_frozen.json")
        for path, digest in result["hashes"].items():
            self.assertEqual(sha256_file(task.REPO/path), digest)
        expected = {(r["candidate_id"], ca) for r in task.contract()[0]["candidate_definitions"] for ca in task.CALCIUM}
        for name in ("ae4_5pct_results.csv", "ae2_results.csv", "allocation_vs_genotype.csv"):
            with (task.OUT/name).open() as f: rows=list(csv.DictReader(f))
            self.assertEqual(len(rows), 90)
            self.assertEqual({(r["candidate_id"],float(r["calcium_uM"])) for r in rows}, expected)
            self.assertEqual(sum(r["status"] == "BASELINE_HASH_VALID_REUSE" for r in rows), 30)
            for r in rows:
                if r["condition"] != "inherited_baseline":
                    self.assertEqual(r["status"], "NOT_EVALUATED_WT_INADMISSIBLE")
                    for k in ("R_AE4_5pct", "R_AE2", "D_AE4_5pct", "D_AE2"):
                        if k in r: self.assertEqual(r[k], "")
        for k in ("new_ae4_near_loss_runs", "new_ae2_loss_runs", "new_genotype_resting_solves", "new_exact_zero_ae4_attempts"):
            self.assertEqual(result[k], 0)

    def test_baseline_ratios_contrast_and_two_routing_families(self):
        with (task.OUT/"allocation_vs_genotype.csv").open() as f: rows=list(csv.DictReader(f))
        rows=[r for r in rows if r["condition"] == "inherited_baseline"]
        self.assertEqual({r["routing_family"] for r in rows}, {"AE4NA05", "AE4NA20"})
        self.assertEqual({float(r["calcium_uM"]) for r in rows}, set(task.CALCIUM))
        for r in rows:
            ratio = float(r["ae4_5pct_total_0_600_pL"])/float(r["wt_total_0_600_pL"])
            self.assertAlmostEqual(ratio, float(r["R_AE4_5pct"]), places=14)
            self.assertAlmostEqual(1-ratio, float(r["D_AE4_5pct"]), places=14)
            self.assertAlmostEqual(float(r["D_AE4_5pct"])-float(r["D_AE2"]), float(r["AE4_minus_AE2_reduction"]), places=14)
            self.assertGreater(ratio, 1.)
            self.assertLess(abs(float(r["D_AE2"])), .002)

    def test_modified_candidate_cannot_reuse_old_five_percent_state(self):
        fake={"wt_admissible":True,"condition":"share_0.10"}
        with self.assertRaises(AssertionError): post.baseline_pair(fake, {}, {}, {})


if __name__ == "__main__":
    unittest.main()
