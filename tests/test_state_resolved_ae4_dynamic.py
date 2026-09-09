"""Regression tests for Task-13 Rounds 6--7 validation gating.

The only integrated trajectory here is the candidate-independent AE4-null
code-coordinate diagnostic.  It deliberately has no WT denominator and is
never converted into a secretion ratio or a physical delay.
"""

from __future__ import annotations

import csv
import inspect
import json
from pathlib import Path
import unittest

import numpy as np

from ae4_mechanism_reconstruction.chassis import (
    AE4_KNOCKOUT,
    STATE_ORDER,
    FixedChassis,
    ZeroAE4,
)
from state_resolved_ae4.adapter import StateResolvedAE4Adapter
from state_resolved_ae4.protocols import (
    CALCIUM_STEP_TIME,
    historical_calcium_step,
    immediate_pka_step,
    stage_a_time_grid,
)
from state_resolved_ae4.validation import (
    PredictiveEligibility,
    PredictiveGateError,
    require_physical_timing,
    secretion_ratios,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = REPOSITORY_ROOT / "results" / "13_state_resolved_ae4"
COMMON_NULL_ROOT = np.asarray(
    (
        118.56527670449147,
        5.5593461776616495,
        65.66070346681501,
        21.132282246431384,
        123.70751315363735,
        48.78221502702161,
        58.90444835892977,
    )
)


def _csv_rows(name: str) -> list[dict[str, str]]:
    with (RESULT_ROOT / name).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if any(None in row for row in rows):
        raise AssertionError(f"{name} is not rectangular")
    return rows


class TestValidationFirewall(unittest.TestCase):
    def test_validation_api_blocks_ratio_and_timing_upstream(self) -> None:
        eligibility = PredictiveEligibility(
            wt_capacity_valid=False,
            wt_physiology_gate_pass=False,
            chassis_independently_calibrated=False,
            physical_time_certified=False,
        )
        with self.assertRaises(PredictiveGateError):
            secretion_ratios(
                np.asarray((0.0, 1.0)),
                np.asarray((1.0, 1.0)),
                np.asarray((1.0, 0.65)),
                eligibility,
            )
        with self.assertRaises(PredictiveGateError):
            require_physical_timing(eligibility)

    def test_no_candidate_has_a_valid_wt_denominator(self) -> None:
        summary = json.loads((RESULT_ROOT / "stage_a_summary.json").read_text())
        self.assertEqual(summary["valid_wt_capacity_rows"], 0)
        self.assertEqual(summary["wt_fixed_chassis_pass_rows"], 0)
        candidates = [
            row
            for row in _csv_rows("dynamic_trajectory_status.csv")
            if row["record_type"] == "validation_candidate"
        ]
        self.assertEqual(len(candidates), 10)
        for row in candidates:
            self.assertEqual(row["status"], "blocked_before_dynamics")
            self.assertEqual(row["wt_capacity_valid"], "false")
            self.assertEqual(row["wt_gate2_pass"], "false")
            self.assertEqual(row["endpoint_ko_wt_ratio"], "")
            self.assertEqual(row["integral_ko_wt_ratio"], "")
            self.assertEqual(row["first_2_3_min_status"], "not_applicable")
            self.assertEqual(row["sustained_late_status"], "not_applicable")

        round_five = _csv_rows("chassis_reconstruction_scan.csv")
        attempts = [row for row in round_five if row["record_type"] == "wt_capacity_attempt"]
        self.assertEqual(len(attempts), 5)
        self.assertTrue(all(row["valid_wt_capacity_available"] == "False" for row in attempts))
        self.assertTrue(all(row["wt_gate_evaluable"] == "False" for row in attempts))
        self.assertTrue(all(row["secretion_target_read"] == "False" for row in round_five))

    def test_strict_targets_are_references_not_joined_predictions(self) -> None:
        rows = _csv_rows("dynamic_trajectory_status.csv")
        self.assertTrue(all(row["strict_target_used_in_fit"] == "false" for row in rows))
        self.assertTrue(
            all(row["strict_target_joined_to_prediction"] == "false" for row in rows)
        )
        primary = {
            row["quantity"]: row
            for row in rows
            if row["record_type"] == "heldout_target_reference"
            and row["target_id"] == "E15-02"
        }
        self.assertEqual(float(primary["ten_min_total_loss"]["value"]), 0.35)
        self.assertEqual(float(primary["ten_min_total_loss"]["uncertainty"]), 0.047)
        self.assertEqual(float(primary["ten_min_total_KO_WT_complement"]["value"]), 0.65)
        rounded = next(row for row in rows if row["record_type"] == "context_not_target")
        self.assertEqual(float(rounded["value"]), 0.30)
        self.assertEqual(rounded["status"], "not_an_experimental_endpoint")

    def test_protocol_has_no_physical_delay_or_carrier_memory(self) -> None:
        self.assertEqual(CALCIUM_STEP_TIME, 100.0)
        self.assertEqual(historical_calcium_step(100.0), 0.05)
        self.assertEqual(historical_calcium_step(100.0 + 1e-7), 0.55)
        self.assertEqual(immediate_pka_step(100.0), 0.0)
        self.assertEqual(immediate_pka_step(100.0 + 1e-7), 1.0)
        self.assertEqual(len(STATE_ORDER), 7)
        self.assertNotIn("occupancy", STATE_ORDER)
        adapter_source = inspect.getsource(StateResolvedAE4Adapter).lower()
        self.assertIn("quasi-steady", adapter_source)
        self.assertIn("stationary occupancy", adapter_source)

    def test_saved_diagnostic_is_not_a_heldout_trajectory(self) -> None:
        rows = _csv_rows("dynamic_diagnostic_trajectory.csv")
        self.assertEqual(len(rows), 24)
        self.assertEqual({row["solver_method"] for row in rows}, {"Radau", "BDF"})
        self.assertTrue(all(row["physical_time_min"] == "" for row in rows))
        self.assertTrue(all(row["wt_denominator_valid"] == "false" for row in rows))
        self.assertTrue(all(row["heldout_prediction"] == "false" for row in rows))
        state_fields = (
            "na_l_mM",
            "k_l_mM",
            "height_um",
            "na_i_mM",
            "k_i_mM",
            "cl_i_mM",
            "hco3_i_mM",
        )
        self.assertTrue(
            all(float(row[field]) > 0.0 for row in rows for field in state_fields)
        )


class TestCodeCoordinateDiagnostic(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        chassis = FixedChassis(ZeroAE4(), calcium_protocol=historical_calcium_step)
        cls.results = {
            method: chassis.simulate(
                AE4_KNOCKOUT,
                t_eval=stage_a_time_grid(),
                method=method,
                rtol=1e-8,
                atol=1e-10,
                y0=COMMON_NULL_ROOT,
            )
            for method in ("Radau", "BDF")
        }

    def test_two_solvers_reproduce_saved_unpromoted_metrics(self) -> None:
        status = _csv_rows("dynamic_trajectory_status.csv")
        recorded = {
            (row["quantity"], row["solver_status"].split("_")[0]): float(row["value"])
            for row in status
            if row["record_type"] == "code_diagnostic"
        }
        for method, result in self.results.items():
            q_total = np.asarray([row["q_total"] for row in result.diagnostics])
            post = result.time > CALCIUM_STEP_TIME
            integral = float(np.trapezoid(q_total[post], result.time[post]))
            self.assertAlmostEqual(
                float(result.endpoint["q_total"]),
                recorded[("endpoint_q_total", method)],
                delta=2e-13,
            )
            self.assertAlmostEqual(
                integral,
                recorded[("post_step_integral", method)],
                delta=2e-11,
            )
            self.assertGreater(float(np.min(result.state)), 0.0)
            conservation = max(
                abs(float(value))
                for row in result.diagnostics
                for value in row["checks"].values()
            )
            self.assertLessEqual(conservation, 5e-14)
        radau = self.results["Radau"]
        bdf = self.results["BDF"]
        self.assertLessEqual(float(np.max(np.abs(radau.state - bdf.state))), 2.5e-6)
        radau_q = np.asarray([row["q_total"] for row in radau.diagnostics])
        bdf_q = np.asarray([row["q_total"] for row in bdf.diagnostics])
        self.assertLessEqual(float(np.max(np.abs(radau_q - bdf_q))), 1.4e-9)


if __name__ == "__main__":
    unittest.main()
