"""Regression checks for Task-13 machine-readable Stage-A artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/13_state_resolved_ae4"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(name: str) -> list[dict[str, str]]:
    with (RESULT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class StageAArtifactTests(unittest.TestCase):
    def test_original_reveal_bytes_and_post_reveal_amendment_are_auditable(self) -> None:
        marker = json.loads((RESULT / "stage_a_reveal_marker.json").read_text())
        amendment = json.loads((RESULT / "amendment_01.json").read_text())
        original = RESULT / "transporter_fit_summary_pre_reveal.csv"
        amended = RESULT / "transporter_fit_summary_amended.csv"
        current = RESULT / "transporter_fit_summary.csv"
        declared = marker["file_hashes"][
            "results/13_state_resolved_ae4/transporter_fit_summary.csv"
        ]
        self.assertEqual(sha256(original), declared)
        self.assertEqual(amendment["pre_reveal_summary_hash"], declared)
        self.assertEqual(sha256(amended), amendment["amended_summary_hash"])
        self.assertEqual(sha256(current), amendment["required_summary_current_hash"])
        self.assertEqual(sha256(amended), sha256(current))
        self.assertFalse(amendment["knockout_ionic_values_used"])
        self.assertFalse(amendment["knockout_secretion_values_used"])
        self.assertFalse(amendment["stage_a_outcome_changed"])

    def test_amended_transporter_summary_has_five_conditional_survivors(self) -> None:
        records = rows("transporter_fit_summary.csv")
        survivors = [row for row in records if row["transporter_gate_pass"] == "True"]
        self.assertEqual(len(survivors), 5)
        self.assertTrue(all(row["heldout_target_read"] == "False" for row in records))
        carbonate = [row for row in records if row["model_id"] == "SR4A_CARBONATE_111"]
        self.assertEqual(len(carbonate), 3)
        self.assertTrue(all(row["na_ec50_pred_mM"] not in {"", "nan"} for row in carbonate))
        self.assertTrue(all(row["transporter_gate_pass"] == "False" for row in carbonate))

    def test_every_transporter_pass_in_finite_gauge_grid_has_wrong_cl_direction(self) -> None:
        records = rows("barrier_gauge_source_scan.csv")
        self.assertEqual(
            sha256(RESULT / "barrier_gauge_source_scan.csv"),
            json.loads((RESULT / "amendment_01.json").read_text())["gauge_scan_hash"],
        )
        passes = [row for row in records if row["transporter_summary_pass"] == "True"]
        self.assertGreaterEqual(len(passes), 5)
        self.assertTrue(
            all(row["cl_loading_direction_matches_fixed_wt_demand"] == "False" for row in passes)
        )
        row_point_one = min(
            records, key=lambda row: abs(float(row["common_barrier_multiplier"]) - 0.1)
        )
        self.assertEqual(row_point_one["transporter_summary_pass"], "False")

    def test_no_failed_capacity_is_promoted_and_all_exact_wt_roots_fail_ph(self) -> None:
        attempts = rows("wt_capacity_attempts.csv")
        self.assertEqual(len(attempts), 5)
        self.assertTrue(all(row["calibration_valid"] == "False" for row in attempts))
        self.assertTrue(all(float(row["max_abs_raw_rhs"]) > 1e-8 for row in attempts))
        root_records = rows("stage_a_roots.csv")
        self.assertEqual(len(root_records), 5)
        self.assertTrue(all(float(row["max_abs_raw_rhs"]) <= 1e-8 for row in root_records))
        self.assertTrue(all(row["wt_physiology_gate_pass"] == "False" for row in root_records))
        self.assertTrue(all(abs(float(row["ph_z_wt"])) > 6.0 for row in root_records))

    def test_common_null_root_and_strict_secretion_non_evaluation(self) -> None:
        null = rows("stage_a_null_roots.csv")
        self.assertEqual(len(null), 1)
        row = null[0]
        self.assertLess(float(row["max_abs_raw_rhs"]), 1e-12)
        self.assertAlmostEqual(float(row["cl_i_mM"]), 48.782215027, places=8)
        self.assertAlmostEqual(float(row["pH_i"]), 7.378473911, places=8)
        self.assertGreater(float(row["cl_target_z"]), 7.6)
        self.assertGreater(float(row["ph_target_z"]), 24.0)
        predictions = rows("stage_a_predictions.csv")
        self.assertEqual(len(predictions), 5)
        for prediction in predictions:
            self.assertEqual(
                prediction["secretion_prediction_status"],
                "not_admissible_after_WT_gate_failure",
            )
            self.assertEqual(prediction["ae4_ko_endpoint_flow_ratio"], "")
            self.assertEqual(prediction["ae4_ko_integrated_flow_ratio"], "")
            self.assertEqual(prediction["strict_secretion_used_in_fit"], "false")

    def test_ctmc_snapshot_manifest_hashes_match(self) -> None:
        manifest = json.loads((RESULT / "ctmc_snapshot_manifest.json").read_text())
        self.assertEqual(sha256(RESULT / "ctmc_frozen_conditions.csv"), manifest["conditions_sha256"])
        self.assertEqual(sha256(RESULT / "ctmc_frozen_states.csv"), manifest["states_sha256"])
        self.assertEqual(sha256(RESULT / "ctmc_frozen_edges.csv"), manifest["edges_sha256"])
        self.assertEqual(manifest["condition_count"], 17)
        self.assertFalse(manifest["strict_holdout_read"])


if __name__ == "__main__":
    unittest.main()
