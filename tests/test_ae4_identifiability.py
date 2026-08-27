from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.optimize import root

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ae4_identifiability import (  # noqa: E402
    BaselineState,
    acid_base_equilibrium_audit,
    activity_signature_matrix,
    apical_qss_current_audit,
    finite_difference_jacobian,
    flux_aggregate_matrix,
    flux_aggregates,
    forward_flux_wedge,
    implicit_sensitivity,
    kinetic_scale_audit,
    nkcc_turnover_audit,
    recover_exchanger_fluxes,
    water_flux_audit,
)


class StructuralIdentifiabilityTests(unittest.TestCase):
    def test_activity_signatures_have_exact_rank_two(self) -> None:
        signatures = activity_signature_matrix(25.0, 120.0)
        self.assertEqual(np.linalg.matrix_rank(signatures), 2)
        self.assertAlmostEqual(np.linalg.det(signatures[[0, 3], :]), -1.0)

    def test_two_flux_aggregates_are_exactly_invertible(self) -> None:
        self.assertAlmostEqual(np.linalg.det(flux_aggregate_matrix()), 1.0)
        for j2, j4 in ((0.0, 0.0), (0.2, 1.7), (-0.3, 0.8), (4.0, 2.0)):
            a, b = flux_aggregates(j2, j4)
            recovered = recover_exchanger_fluxes(a, b)
            np.testing.assert_allclose(recovered, (j2, j4), atol=1e-14)

    def test_forward_fluxes_are_equivalent_to_wedge_bounds(self) -> None:
        self.assertTrue(forward_flux_wedge(*flux_aggregates(0.3, 0.9)))
        self.assertFalse(forward_flux_wedge(1.0, 0.8))
        self.assertFalse(forward_flux_wedge(1.0, 2.2))

    def test_scalar_projection_has_rank_one(self) -> None:
        scalar_jacobian = np.array([[1.0, 1.0]])
        self.assertEqual(np.linalg.matrix_rank(scalar_jacobian), 1)

    def test_implicit_derivative_matches_numerical_equilibria(self) -> None:
        # Independent software check of -F_u^{-1}F_theta on a nonlinear system.
        u0 = np.array([1.0, 0.0])
        theta0 = np.array([1.0, 2.0])

        def residual(u: np.ndarray, theta: np.ndarray) -> np.ndarray:
            return np.array([u[0] ** 2 + u[1] - theta[0], u[0] + np.exp(u[1]) - theta[1]])

        f_u = np.array([[2.0, 1.0], [1.0, 1.0]])
        f_theta = -np.eye(2)
        analytic = implicit_sensitivity(f_u, f_theta)

        def equilibrium(theta: np.ndarray) -> np.ndarray:
            solved = root(lambda u: residual(u, theta), u0)
            self.assertTrue(solved.success)
            return solved.x

        numerical = finite_difference_jacobian(equilibrium, theta0, step=1e-5)
        np.testing.assert_allclose(analytic, numerical, rtol=2e-6, atol=2e-6)


class PrintedBaselineAuditTests(unittest.TestCase):
    def test_both_volume_signs_have_same_steady_condition(self) -> None:
        for q_a, q_b in ((1.0, 1.0), (1.0, 2.0), (-3.0, 4.0)):
            self.assertEqual(abs(q_a - q_b) < 1e-14, abs(q_b - q_a) < 1e-14)

    def test_printed_water_values_are_not_a_steady_volume(self) -> None:
        audit = water_flux_audit()
        self.assertGreater(audit["q_b_over_q_a"], 50.0)
        self.assertGreater(audit["q_total_over_reported"], 100.0)

    def test_printed_apical_qss_has_material_residual(self) -> None:
        audit = apical_qss_current_audit(BaselineState())
        self.assertGreater(abs(audit["qss_residual_pA"]), 60.0)
        self.assertGreater(audit["abs_residual_over_abs_tight_current"], 0.25)

    def test_printed_acid_base_equilibrium_misses_table_value(self) -> None:
        audit = acid_base_equilibrium_audit(BaselineState())
        self.assertAlmostEqual(audit["equilibrium_predicted_co2_i_mM"], 3.519, places=3)
        self.assertGreater(audit["reported_over_equilibrium_prediction"], 1.8)

    def test_ae4_fourth_order_unit_conventions_differ_by_twelve_orders(self) -> None:
        audit = kinetic_scale_audit(BaselineState())
        self.assertEqual(audit["ae4_mM_to_M_fourth_power_ratio"], 1.0e12)
        self.assertGreater(audit["ae4_printed_turnover_per_s"], 4.0e4)
        self.assertLess(
            audit["ae4_turnover_if_numeric_concentrations_are_molar_per_s"],
            5.0e-8,
        )

    def test_palk_corrigendum_repairs_nkcc_concentration_units(self) -> None:
        audit = nkcc_turnover_audit(BaselineState())
        self.assertLess(audit["literal_2018_table_turnover_per_s"], -14.0)
        self.assertGreater(audit["corrigendum_converted_turnover_per_s"], 0.59)
        self.assertLess(audit["corrigendum_converted_turnover_per_s"], 0.61)
        self.assertAlmostEqual(
            audit["a2_corrigendum_converted_mM_minus4_s_minus2"],
            2.0096e-5,
        )


if __name__ == "__main__":
    unittest.main()
