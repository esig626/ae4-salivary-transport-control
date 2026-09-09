"""Scientific and numerical checks for the structural diagnostic."""

import unittest

import numpy as np

from identifiability_discrimination import LinearizedPumpLeakModel, printed_activity_factor_check


class ReducedModelTests(unittest.TestCase):
    def setUp(self):
        self.model = LinearizedPumpLeakModel()

    def test_published_baseline_is_exact_anchor(self):
        state = self.model.steady_state((1.0, 1.0))
        np.testing.assert_allclose(state, np.ones(4), atol=1e-14, rtol=0)
        np.testing.assert_allclose(self.model.residual(state, (1.0, 1.0)), 0.0, atol=1e-14, rtol=0)

    def test_published_stoichiometric_signatures(self):
        s = self.model.stoichiometric_matrix
        self.assertEqual(s[0, 0], 0.0)  # AE2 has no cation term
        self.assertEqual(s[1, 0], 0.0)
        self.assertGreater(s[2, 0], 0.0)  # both import Cl
        self.assertGreater(s[2, 1], 0.0)
        self.assertLess(s[3, 0], 0.0)  # both export bicarbonate
        self.assertLess(s[3, 1], 0.0)
        self.assertLess(s[0, 1], 0.0)  # AE4 exports cations
        self.assertLess(s[1, 1], 0.0)

    def test_printed_activity_laws_are_transcribed(self):
        check = printed_activity_factor_check()
        self.assertAlmostEqual(check["AE2_dimensionless_bracket"], -0.002695445332692506)
        self.assertAlmostEqual(check["AE4_printed_kinetic_bracket"], 41646.261368249994)

    def test_implicit_derivative_matches_central_difference(self):
        numerical = self.model.finite_difference_sensitivity(step=1e-6)
        np.testing.assert_allclose(numerical, self.model.state_sensitivity, atol=2e-10, rtol=2e-9)

    def test_one_scalar_cannot_have_rank_two(self):
        for name in self.model.observation_rows:
            jacobian = self.model.observation_jacobian((name,))
            self.assertLessEqual(np.linalg.matrix_rank(jacobian), 1)

    def test_candidate_two_measurement_panels(self):
        # These are properties of the declared diagnostic, not of the unreproduced full model.
        for panel in (("Q_star", "Na_i"), ("Cl_i", "HCO3_i"), ("Na_i", "Cl_i")):
            self.assertEqual(self.model.panel_metrics(panel)["rank"], 2)

    def test_q_equivalence_is_exact_in_affine_diagnostic(self):
        points = self.model.q_equivalence_segment(n_points=37)
        values = [self.model.observation_values(point)["Q_star"] for point in points]
        reference = self.model.observation_values((0.5, 0.5))["Q_star"]
        np.testing.assert_allclose(values, reference, atol=2e-14, rtol=0)

    def test_consistent_equation_sign_change_leaves_sensitivity_unchanged(self):
        for signs in ((1, 1, 1, -1), (-1, 1, -1, 1), (-1, -1, -1, -1)):
            self.assertLess(self.model.row_sign_invariance_error(signs), 1e-14)

    def test_declared_grid_stays_positive(self):
        for g2 in np.linspace(0.0, 1.0, 11):
            for g4 in np.linspace(0.0, 1.0, 11):
                self.assertGreater(np.min(self.model.steady_state((g2, g4))), 0.0)


if __name__ == "__main__":
    unittest.main()
