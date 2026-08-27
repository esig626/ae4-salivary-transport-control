"""Independent numerical-route tests for Task 13."""

from __future__ import annotations

import math
import unittest

import numpy as np

from state_resolved_ae4.independent import (
    chemical_affinity,
    fixed_chassis_ae4_null_diagnostics,
    fixed_chassis_ae4_null_coordinate_residual,
    fixed_null_coordinates_to_state,
    fixed_null_state_to_coordinates,
    fixed_null_required_balance,
    exact_rational_rank,
    log_rate_affinity,
    ring_edge_currents,
    ring_generator,
    shared_two_state_check,
    solve_positive_root,
    source_projection,
    split_cation_ae4_null_residual,
    split_cation_null_coordinate_residual,
    split_cation_sr2_wt_residual,
    stage_b_correction_jacobian,
    stage_b_target_state,
    stationary_crosscheck,
)


class IndependentCtmcTests(unittest.TestCase):
    def test_svd_and_matrix_tree_match_on_asymmetric_ring(self) -> None:
        forward = np.asarray([2.1, 0.7, 4.2, 1.3, 0.9, 3.4])
        reverse = np.asarray([0.6, 1.8, 0.5, 2.7, 1.1, 0.8])
        generator = ring_generator(forward, reverse)
        check = stationary_crosscheck(generator)

        self.assertLess(check.max_probability_difference, 2e-14)
        self.assertLess(check.svd_residual_inf, 2e-14)
        self.assertLess(check.matrix_tree_residual_inf, 2e-14)
        self.assertAlmostEqual(float(np.sum(check.svd)), 1.0, places=14)
        self.assertGreater(float(np.min(check.svd)), 0.0)

        currents = ring_edge_currents(check.svd, forward, reverse)
        self.assertLess(float(np.ptp(currents)), 2e-14)
        self.assertEqual(np.sign(currents[0]), np.sign(log_rate_affinity(forward, reverse)))

    def test_shared_branch_formula_and_carrier_closure(self) -> None:
        rates = dict(a=2.4, b=0.8, c_na=0.4, d_na=1.2, c_k=3.1, d_k=0.3)
        check = shared_two_state_check(**rates)
        sigma = sum(rates.values())
        expected_o = (rates["b"] + rates["c_na"] + rates["c_k"]) / sigma
        expected_i = (rates["a"] + rates["d_na"] + rates["d_k"]) / sigma

        self.assertAlmostEqual(check.occupancy_o, expected_o, places=14)
        self.assertAlmostEqual(check.occupancy_i, expected_i, places=14)
        self.assertLess(abs(check.carrier_residual), 2e-15)
        self.assertGreaterEqual(check.entropy_production, -2e-14)

    def test_salivary_affinities_reproduce_frozen_thermodynamic_audit(self) -> None:
        a = {
            "cl_o": 124.6,
            "cl_i": 50.1,
            "b_o": 24.7,
            "b_i": 19.0,
            "na_o": 151.6,
            "na_i": 15.5,
            "k_o": 3.4,
            "k_i": 139.5,
        }
        na_112 = chemical_affinity(
            a,
            {"cl_o": 1, "na_i": 1, "b_i": 2},
            {"cl_i": 1, "na_o": 1, "b_o": 2},
        )
        k_112 = chemical_affinity(
            a,
            {"cl_o": 1, "k_i": 1, "b_i": 2},
            {"cl_i": 1, "k_o": 1, "b_o": 2},
        )
        sequential = chemical_affinity(
            a,
            {"cl_o": 1, "na_o": 1, "k_i": 1, "b_i": 1},
            {"cl_i": 1, "na_i": 1, "k_o": 1, "b_o": 1},
        )

        self.assertAlmostEqual(na_112, -1.89405, places=5)
        self.assertAlmostEqual(k_112, 4.10065, places=5)
        self.assertAlmostEqual(sequential, 6.64342, places=5)

    def test_positive_root_uses_log_state_and_scipy_root(self) -> None:
        def residual(x: np.ndarray) -> np.ndarray:
            return np.asarray([x[0] * x[1] - 6.0, x[0] + x[1] - 5.0])

        check = solve_positive_root(residual, [1.4, 3.6])
        self.assertTrue(check.success, check.message)
        self.assertLess(check.residual_inf, 1e-10)
        self.assertTrue(np.all(check.state > 0.0))
        self.assertAlmostEqual(float(np.min(check.state)), 2.0, places=9)
        self.assertAlmostEqual(float(np.max(check.state)), 3.0, places=9)

    def test_equal_forward_reverse_products_give_zero_ring_current(self) -> None:
        forward = np.asarray([2.0, 3.0, 5.0, 7.0])
        reverse = np.asarray([1.0, 2.0, 3.0, 35.0])
        self.assertAlmostEqual(math.prod(forward), math.prod(reverse), places=14)
        check = stationary_crosscheck(ring_generator(forward, reverse))
        currents = ring_edge_currents(check.svd, forward, reverse)
        self.assertLess(float(np.max(np.abs(currents))), 2e-14)

    def test_fixed_chassis_ae4_null_root_matches_saved_task12_result(self) -> None:
        # Start from a deliberately perturbed state while preserving positivity
        # through an explicit h_i coordinate.  The independently coded residual
        # contains no candidate-family argument, making null independence exact.
        reference = np.asarray(
            [
                118.56527670449144,
                5.5593461776616655,
                65.66070346655056,
                21.132282246430893,
                123.70751315363782,
                48.782215027024655,
                58.90444835877707,
            ],
            dtype=float,
        )
        initial = fixed_null_state_to_coordinates(reference) * np.asarray(
            [1.01, 0.99, 1.02, 0.98, 1.01, 0.97, 1.5]
        )
        check = solve_positive_root(
            fixed_chassis_ae4_null_coordinate_residual, initial, tolerance=1e-11
        )
        state = fixed_null_coordinates_to_state(check.state)
        diagnostics = fixed_chassis_ae4_null_diagnostics(state)

        self.assertTrue(check.success, check.message)
        self.assertLess(check.residual_inf, 2e-11)
        self.assertAlmostEqual(float(state[5]), 48.7822150270, places=8)
        self.assertAlmostEqual(diagnostics.pH_i, 7.3784739108, places=8)
        self.assertAlmostEqual(diagnostics.cell_volume_pL, 2.9741782058, places=8)
        self.assertLess(abs(diagnostics.current_residual_apical), 2e-12)
        self.assertLess(abs(diagnostics.current_residual_basolateral), 2e-12)

    def test_exact_rank_and_nonnegative_projection_are_independent_routes(self) -> None:
        signatures = np.asarray(
            [
                [-3.0, 1.0, 1.0],
                [2.0, -1.0, 1.0],
                [0.0, 0.0, 2.0],
            ]
        )
        target = signatures @ np.asarray([0.7, 0.2, 0.4])
        check = source_projection(target, signatures)

        self.assertEqual(exact_rational_rank(signatures.astype(int).tolist()), 3)
        self.assertEqual(check.numerical_rank, 3)
        self.assertLess(check.residual_norm_unconstrained, 2e-14)
        self.assertLess(check.residual_norm_nonnegative, 2e-14)
        np.testing.assert_allclose(check.coefficients_nonnegative, [0.7, 0.2, 0.4])

    def test_stage_b_target_family_has_rank_five_source_variation(self) -> None:
        free = np.asarray(
            [
                118.56527670449147,
                5.5593461776616495,
                65.66070346681501,
                21.132282246431384,
                123.70751315363735,
            ]
        )
        target_state = stage_b_target_state(free)
        demand = fixed_null_required_balance(target_state)
        expected_cell = np.asarray(
            [-0.1577266452, -0.08097603657, -0.1668785093, 0.4297469356, 0.5015711081]
        )
        np.testing.assert_allclose(demand[:5], expected_cell, rtol=0.0, atol=8e-11)

        jacobian = stage_b_correction_jacobian(free, relative_step=1e-3)
        singular_values = np.linalg.svd(jacobian, compute_uv=False)
        self.assertEqual(int(np.linalg.matrix_rank(jacobian)), 5)
        self.assertGreater(float(singular_values[-1]), 4e-4)
        self.assertAlmostEqual(float(singular_values[0]), 0.40223543, places=7)

    def test_round5_m0_is_exactly_nested_and_nonzero_split_root_closes(self) -> None:
        m0 = np.asarray(
            [
                118.56527670449145,
                5.559346177661645,
                65.66070346681074,
                21.132282246429515,
                123.70751315363918,
                48.78221502702284,
                58.90444835892612,
            ]
        )
        np.testing.assert_allclose(
            split_cation_ae4_null_residual(
                m0, apical_pump_fraction=0.0, apical_k_fraction=0.0
            ),
            fixed_chassis_ae4_null_coordinate_residual(
                fixed_null_state_to_coordinates(m0)
            ),
            rtol=0.0,
            atol=2e-16,
        )
        initial = fixed_null_state_to_coordinates(m0)
        solved = solve_positive_root(
            lambda coordinates: split_cation_null_coordinate_residual(
                coordinates,
                apical_pump_fraction=0.4,
                apical_k_fraction=0.2,
            ),
            initial,
            tolerance=1e-11,
        )
        state = fixed_null_coordinates_to_state(solved.state)
        self.assertTrue(solved.success, solved.message)
        self.assertLess(solved.residual_inf, 2e-10)
        self.assertAlmostEqual(float(state[5]), 47.6511637292, places=7)

    def test_round5_frozen_sr2_wt_residual_is_independently_reproduced(self) -> None:
        state = np.asarray(
            [
                118.55724378882537,
                5.5572668338680575,
                199.9999946584882,
                19.41068062118391,
                125.41973291322965,
                50.099603707889464,
                82.5333220078952,
            ]
        )
        expected = np.asarray(
            [
                9.33052871077e-11,
                4.02586818049e-13,
                -1.91527665379e-7,
                -1.25246811415e-5,
                -3.05047984524e-5,
                1.42344911995e-5,
                -2.89823885822e-6,
            ]
        )
        actual = split_cation_sr2_wt_residual(
            state,
            capacity=6.68829505457301e-6,
            apical_pump_fraction=0.0,
            apical_k_fraction=0.0,
        )
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=8e-16)


if __name__ == "__main__":
    unittest.main()
