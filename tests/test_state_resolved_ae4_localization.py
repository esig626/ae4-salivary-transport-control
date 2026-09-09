"""Regression tests for Task-13 Stage-B missing-balance localization."""

from __future__ import annotations

import csv
import inspect
import math
from pathlib import Path
import unittest

import numpy as np

from ae4_mechanism_reconstruction.chassis import (
    AE4_KNOCKOUT,
    FixedChassis,
    HistoricalSodiumOnlyAE4,
    ZeroAE4,
)
from state_resolved_ae4 import localization


KO_ROOT = np.asarray(
    (
        118.56527670449144,
        5.5593461776616655,
        65.66070346655056,
        21.132282246430893,
        123.70751315363782,
        48.782215027024655,
        58.90444835877707,
    ),
    dtype=float,
)


class TestStageBFirewall(unittest.TestCase):
    def setUp(self) -> None:
        self.chassis = FixedChassis()

    def test_targets_remain_sealed_before_stage_a_failure(self) -> None:
        with self.assertRaises(PermissionError):
            localization.released_stage_b_manifold(
                self.chassis,
                stage_a_failed=False,
                cl_i_mM=36.5,
                pH_i=6.89,
            )
        with self.assertRaises(PermissionError):
            localization.required_balance_correction(
                self.chassis, KO_ROOT, stage_a_failed=False
            )

    def test_localization_module_contains_no_secretion_target_number(self) -> None:
        source = inspect.getsource(localization)
        self.assertNotIn("0.65", source)
        self.assertNotIn("35 +/-", source)
        self.assertNotIn("35%", source)


class TestMachineReadableArtifact(unittest.TestCase):
    def test_csv_is_rectangular_and_reproduces_correction_family_rank(self) -> None:
        path = (
            Path(__file__).resolve().parents[1]
            / "results"
            / "13_state_resolved_ae4"
            / "module_projection.csv"
        )
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertTrue(rows)
        self.assertTrue(all(None not in row for row in rows))
        jacobian_rows = [
            row
            for row in rows
            if row["record_type"] == "correction_family_jacobian_row"
        ]
        self.assertEqual(len(jacobian_rows), 12)
        jacobian = np.asarray(
            [
                [float(value) for value in row["signature_or_response"].split(";")]
                for row in jacobian_rows
            ]
        )
        singular_values = np.linalg.svd(jacobian, compute_uv=False)
        rank_row = next(
            row for row in rows if row["record_type"] == "correction_family_rank"
        )
        recorded = np.asarray(
            [float(value) for value in rank_row["signature_or_response"].split(";")]
        )
        np.testing.assert_allclose(singular_values, recorded, rtol=2e-10, atol=1e-12)
        self.assertEqual(np.linalg.matrix_rank(jacobian), int(rank_row["rank"]))


class TestTargetManifold(unittest.TestCase):
    def setUp(self) -> None:
        self.chassis = FixedChassis()
        self.manifold = localization.released_stage_b_manifold(
            self.chassis,
            stage_a_failed=True,
            cl_i_mM=36.5,
            pH_i=6.89,
        )

    def test_two_observations_leave_five_free_state_coordinates(self) -> None:
        self.assertEqual(self.manifold.dimension, 5)
        self.assertAlmostEqual(self.manifold.h_i_mM, 1000.0 * 10.0**-6.89)
        self.assertEqual(
            self.manifold.affine_coefficients,
            (-1.0, -1.0, self.chassis.parameters.impermeant_charge_amount, 1.0),
        )

    def test_state_satisfies_released_cl_ph_and_electroneutrality(self) -> None:
        state = self.manifold.state(
            na_l=KO_ROOT[0],
            k_l=KO_ROOT[1],
            height=KO_ROOT[2],
            na_i=KO_ROOT[3],
            k_i=KO_ROOT[4],
        )
        residuals = self.manifold.residuals(state)
        self.assertLess(max(abs(value) for value in residuals.values()), 2e-11)
        self.assertAlmostEqual(state[5], 36.5)
        self.assertNotAlmostEqual(state[6], KO_ROOT[6])

    def test_required_source_changes_across_unmeasured_slices(self) -> None:
        slices = localization.evaluated_target_slices(
            self.chassis,
            KO_ROOT,
            self.manifold,
            stage_a_failed=True,
        )
        first = np.asarray(slices["model_anchored"].required)
        second = np.asarray(slices["wt_reference_cations"].required)
        determined = np.asarray(slices["model_anchored"].determined)
        self.assertGreater(np.linalg.norm(first[determined] - second[determined]), 0.1)
        self.assertFalse(slices["model_anchored"].determined[5])
        self.assertFalse(slices["model_anchored"].determined[9])

    def test_required_source_family_varies_in_all_five_free_directions(self) -> None:
        response = localization.correction_family_jacobian(
            self.chassis,
            self.manifold,
            KO_ROOT[[0, 1, 2, 3, 4]],
            stage_a_failed=True,
        )
        self.assertEqual(response.numerical_rank, 5)
        self.assertEqual(len(response.singular_values), 5)
        self.assertGreater(response.singular_values[-1], 1e-4)

    def test_five_direction_nonuniqueness_survives_one_sem_target_corners(self) -> None:
        for chloride in (34.9, 38.1):
            for p_h in (6.87, 6.91):
                manifold = localization.released_stage_b_manifold(
                    self.chassis,
                    stage_a_failed=True,
                    cl_i_mM=chloride,
                    pH_i=p_h,
                )
                response = localization.correction_family_jacobian(
                    self.chassis,
                    manifold,
                    KO_ROOT[[0, 1, 2, 3, 4]],
                    stage_a_failed=True,
                )
                self.assertEqual(response.numerical_rank, 5)


class TestExactSignatures(unittest.TestCase):
    def test_apical_and_basolateral_topologies_are_cell_projection_aliases(self) -> None:
        pump_a = np.asarray(localization.MODULE_SIGNATURES["apical_nak_pump"])
        pump_b = np.asarray(localization.MODULE_SIGNATURES["basolateral_nak_pump"])
        k_a = np.asarray(localization.MODULE_SIGNATURES["apical_k_channel"])
        k_b = np.asarray(localization.MODULE_SIGNATURES["basolateral_k_channel"])
        np.testing.assert_array_equal(pump_a[:6], pump_b[:6])
        np.testing.assert_array_equal(k_a[:6], k_b[:6])
        self.assertFalse(np.array_equal(pump_a, pump_b))
        self.assertFalse(np.array_equal(k_a, k_b))

    def test_membrane_split_adds_two_independent_extended_directions(self) -> None:
        pump_difference = localization.topology_difference(
            "apical_nak_pump", "basolateral_nak_pump"
        )
        k_difference = localization.topology_difference(
            "apical_k_channel", "basolateral_k_channel"
        )
        differences = np.column_stack((pump_difference, k_difference))
        self.assertEqual(localization.exact_rank(differences.tolist()), 2)
        np.testing.assert_array_equal(pump_difference[:6], 0.0)
        np.testing.assert_array_equal(k_difference[:6], 0.0)

    def test_structural_ranks_are_exact(self) -> None:
        self.assertEqual(localization.exact_signature_rank(), 12)
        self.assertEqual(
            localization.exact_signature_rank(rows=localization.CELL_CHEMICAL_ROWS),
            6,
        )
        self.assertEqual(
            localization.exact_signature_rank(rows=localization.REDUCED_MEASURED_ROWS),
            11,
        )

    def test_fixed_row_pump_plus_nhe_repair_is_exact(self) -> None:
        result = localization.fixed_row_cation_repair(0.03217931476674193)
        np.testing.assert_allclose(result["target"], result["exact_reconstruction"])
        self.assertAlmostEqual(
            result["pump_alignment_absolute"], 5.0 / math.sqrt(26.0)
        )
        self.assertAlmostEqual(
            result["pump_residual_fraction"], 1.0 / math.sqrt(26.0)
        )


class TestNullFieldAndLocalResponse(unittest.TestCase):
    def test_ae4_null_field_is_candidate_independent(self) -> None:
        common, discrepancy = localization.candidate_independent_ae4_null_rhs(
            (
                lambda: FixedChassis(HistoricalSodiumOnlyAE4()),
                lambda: FixedChassis(ZeroAE4()),
            ),
            KO_ROOT,
        )
        expected = FixedChassis().raw_rhs(0.0, KO_ROOT, scenario=AE4_KNOCKOUT)
        np.testing.assert_allclose(common, expected, rtol=0.0, atol=1e-14)
        self.assertLessEqual(discrepancy, 1e-14)

    def test_local_response_retains_voltage_reclosure(self) -> None:
        response = localization.local_parameter_response(
            FixedChassis(), KO_ROOT, "nak_capacity"
        )
        self.assertTrue(math.isfinite(response.d_cl_d_log_parameter))
        self.assertTrue(math.isfinite(response.d_pH_d_log_parameter))
        self.assertGreater(response.jacobian_condition, 1.0)
        self.assertEqual(len(response.raw_state_response), 7)

    def test_one_direction_is_approximate_but_two_independent_directions_span(self) -> None:
        target = np.asarray((-2.0, -1.0))
        responses = {
            "one": np.asarray((1.0, 0.0)),
            "two": np.asarray((0.0, 1.0)),
        }
        rows = localization.sparse_response_projections(
            target, responses, scales=np.ones(2), maximum_modules=2
        )
        pair = next(row for row in rows if len(row.names) == 2)
        self.assertTrue(pair.exact_within_tolerance)
        self.assertEqual(pair.coefficients, (-2.0, -1.0))
        nonnegative = localization.sparse_response_projections(
            target,
            responses,
            scales=np.ones(2),
            maximum_modules=2,
            positive_only=True,
        )
        self.assertTrue(all(row.residual_fraction == 1.0 for row in nonnegative))


if __name__ == "__main__":
    unittest.main()
