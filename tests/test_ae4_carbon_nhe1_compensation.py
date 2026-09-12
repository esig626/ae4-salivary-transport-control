import csv
import json
from pathlib import Path
import unittest

from modern_full_model.run_ae4_carbon_nhe1_compensation import (
    AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S,
    BASELINE_REST_RELATIVE_TOLERANCE,
    BASELINE_TRAJECTORY_RELATIVE_TOLERANCE,
    CO2_FACTORS,
    EXPECTED_PARAMETER_PATHS,
    NHE1_FACTORS,
    OUT,
    build_controlled_rest_model,
    inputs,
    parameter_changes,
)
from modern_full_model.run_ae4_routing_re_equilibrated import _explicit_rest_model
from modern_full_model.validation import sha256_object


def _rows(name: str):
    with (OUT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class CarbonNhe1CompensationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, _, _ = inputs()
        cls.root = sorted(cls.manifest["roots"])[0]
        cls.rest = _rows("resting_states.csv")
        cls.comparisons = _rows("comparison.csv")
        cls.flux = _rows("flux_decomposition.csv")
        cls.trajectories = _rows("trajectory_summary.csv")
        cls.verification = json.loads((OUT / "verification.json").read_text())

    def test_exact_grid_and_counts(self):
        observed = {
            (float(row["co2_factor"]), float(row["nhe1_factor"]))
            for row in self.rest
        }
        self.assertEqual(observed, set((c, n) for c in CO2_FACTORS for n in NHE1_FACTORS))
        self.assertEqual(len(self.rest), 240)
        self.assertEqual(len(self.flux), 120)
        self.assertEqual(len(self.comparisons), 360)
        self.assertEqual(len(self.trajectories), 720)

    def test_exact_task24_baseline_reproduction(self):
        self.assertTrue(self.verification["baseline_control_pass"])
        self.assertLessEqual(
            self.verification["baseline_max_rest_distance_normalized"],
            BASELINE_REST_RELATIVE_TOLERANCE,
        )
        self.assertLessEqual(
            self.verification["baseline_max_wt_total_relative_error"],
            BASELINE_TRAJECTORY_RELATIVE_TOLERANCE,
        )
        self.assertLessEqual(
            self.verification["baseline_max_ae4_5pct_total_relative_error"],
            BASELINE_TRAJECTORY_RELATIVE_TOLERANCE,
        )
        self.assertLessEqual(
            self.verification["baseline_max_ratio_absolute_error"],
            BASELINE_TRAJECTORY_RELATIVE_TOLERANCE,
        )

    def test_only_co2_permeability_and_nhe1_capacity_change(self):
        baseline = _explicit_rest_model(self.manifest, self.root)
        identity = build_controlled_rest_model(self.manifest, self.root, 1.0, 1.0)
        self.assertEqual(parameter_changes(baseline.parameters, identity.parameters), {})
        self.assertEqual(
            sha256_object(baseline.parameters), sha256_object(identity.parameters)
        )
        controlled = build_controlled_rest_model(
            self.manifest, self.root, 10.0, 0.25
        )
        self.assertEqual(
            set(parameter_changes(baseline.parameters, controlled.parameters)),
            set(EXPECTED_PARAMETER_PATHS),
        )
        self.assertEqual(baseline.ae4_parameters, controlled.ae4_parameters)
        self.assertIs(baseline.ae4_evaluator, controlled.ae4_evaluator)

    def test_nkcc1_chloride_flux_uses_factor_two(self):
        for row in self.rest:
            expected = 2.0 * float(row["nkcc1_inward_cycle_flux_fmol_s"])
            self.assertAlmostEqual(
                float(row["nkcc1_chloride_influx_fmol_s"]), expected, places=15
            )
        self.assertEqual(
            self.verification[
                "max_rest_nkcc1_chloride_stoichiometry_residual_fmol_s"
            ],
            0.0,
        )
        self.assertTrue(
            all(
                float(row["max_nkcc1_chloride_stoichiometry_residual_fmol_s"])
                <= AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S
                for row in self.trajectories
            )
        )
        self.assertLessEqual(
            self.verification[
                "max_stimulation_nkcc1_chloride_stoichiometry_residual_fmol_s"
            ],
            AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S,
        )

    def test_ae4_expression_scaling(self):
        for row in self.rest:
            expression = float(row["ae4_expression"])
            residual = abs(
                float(row["ae4_cl_flux_fmol_s"])
                - expression * float(row["ae4_cl_flux_pre_expression_fmol_s"])
            )
            self.assertLessEqual(residual, AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S)
            self.assertAlmostEqual(
                float(row["ae4_flux_per_unit_expression_fmol_s"]),
                float(row["ae4_cl_flux_pre_expression_fmol_s"]),
                places=14,
            )

    def test_ae4_has_no_countercycle(self):
        self.assertTrue(
            all(row["ae4_no_opposing_cation_sources"] == "True" for row in self.rest)
        )
        self.assertTrue(
            all(int(row["opposing_cation_samples"]) == 0 for row in self.trajectories)
        )
        self.assertEqual(
            self.verification["opposing_cation_stimulation_sample_count"], 0
        )

    def test_conservation_and_complete_rest_closure(self):
        self.assertTrue(self.verification["all_rest_gates_pass"])
        self.assertTrue(self.verification["all_stimulation_gates_pass"])
        self.assertTrue(self.verification["all_rest_ae4_charge_neutral"])
        self.assertTrue(self.verification["all_rest_states_positive"])
        self.assertTrue(self.verification["all_rest_states_away_from_boundaries"])
        self.assertLessEqual(self.verification["max_rest_conservation_ratio"], 1.0)
        self.assertLessEqual(
            self.verification["max_stimulation_conservation_ratio"], 1.0
        )

    def test_resting_state_is_independent_of_later_calcium(self):
        groups = {}
        for row in self.trajectories:
            key = (
                row["root_id"],
                float(row["co2_factor"]),
                float(row["nhe1_factor"]),
                float(row["ae4_expression"]),
            )
            groups.setdefault(key, set()).add(row["resting_state_sha256"])
        self.assertEqual(len(groups), 240)
        self.assertTrue(all(len(hashes) == 1 for hashes in groups.values()))
        self.assertTrue(
            self.verification["resting_state_independent_of_subsequent_calcium"]
        )
        self.assertTrue(self.verification["all_stimulation_onsets_match_rest"])


if __name__ == "__main__":
    unittest.main()
