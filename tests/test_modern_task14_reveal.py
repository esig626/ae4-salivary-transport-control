"""Post reveal checks use retained records only, never reopen the source ledgers."""
from pathlib import Path
import json
import unittest
from unittest.mock import patch

from modern_full_model import task14_compare as compare, task14_reveal as reveal


class TestTask14Reveal(unittest.TestCase):
    def test_one_time_source_evaluation_cannot_be_repeated(self):
        with self.assertRaises(RuntimeError):
            reveal.source_records({}, {})

    def test_replay_is_deterministic_and_does_not_change_blind_artifacts(self):
        path = reveal.OUT / "reveal_manifest.json"
        manifest = json.loads(path.read_text())
        first = reveal.evaluate(manifest)
        second = reveal.evaluate(manifest)
        self.assertEqual(first, second)
        for name, value in first[0].items():
            self.assertEqual((reveal.OUT / name).read_bytes(), value)
        for name, value in first[1].items():
            self.assertEqual((reveal.ANALYSIS / name).read_text(), value)
        compare.verify_blind_ledger(reveal.REPO, reveal.OUT / "blind_artifact_hashes.csv")

    def test_unavailable_predictions_are_not_scored_as_zero(self):
        rows = compare.read_csv(reveal.OUT / "ae4_holdout_comparison.csv")
        self.assertEqual(len(rows), 30)
        self.assertTrue(all(row["comparison_status"] == "UNAVAILABLE" for row in rows))
        for row in rows:
            for field in ("signed_error_D_total", "absolute_error_D_total", "inside_reported_interval",
                          "time_course_RMSE", "early_ratio_error", "sustained_ratio_error"):
                self.assertEqual(row[field], "")
            self.assertEqual(row["target_uncertainty_type"], "SEM")
            self.assertEqual(row["experimental_confidence_interval_available"], "False")
        final = json.loads((reveal.OUT / "final_classification.json").read_text())
        self.assertEqual(final["outcome_number"], 6)
        self.assertEqual(final["valid_AE4_pairs"], 0)
        self.assertEqual(len(final["all_30_panel_rows"]), 30)


if __name__ == "__main__": unittest.main()
