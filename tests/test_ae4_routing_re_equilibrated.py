import unittest

import numpy as np

from modern_full_model.run_ae4_routing_only import inputs
from modern_full_model.run_ae4_routing_re_equilibrated import (
    _explicit_rest_model,
    solve_rest_case,
)
from modern_full_model.validation import sha256_object


class ReEquilibratedRoutingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, _, _ = inputs()
        cls.root = sorted(cls.manifest["roots"])[0]

    def test_explicit_rest_model_is_resting_at_all_times_and_preserves_parameters(self):
        model = _explicit_rest_model(self.manifest, self.root)
        payload = self.manifest["roots"][self.root]
        self.assertEqual(
            sha256_object(model.parameters), payload["whole_cell_parameters_sha256"]
        )
        self.assertEqual(
            sha256_object(model.ae4_parameters), payload["ae4_parameters_sha256"]
        )
        state = model.initial_state()
        first = model.evaluate(0.0, state)
        later = model.evaluate(600.0, state)
        np.testing.assert_array_equal(first.rhs, later.rhs)
        self.assertEqual(first.diagnostics.stimulus.beta_input, 0.0)
        self.assertEqual(later.diagnostics.stimulus.beta_input, 0.0)
        self.assertEqual(
            first.diagnostics.stimulus.calcium_uM,
            model.stimulus.resting_calcium_uM,
        )
        self.assertEqual(
            later.diagnostics.stimulus.calcium_uM,
            model.stimulus.resting_calcium_uM,
        )

    def test_both_prescribed_warm_starts_reach_one_accepted_root(self):
        result = solve_rest_case((self.root, 0.10, 1.0))
        self.assertTrue(result["rest_gate_pass"])
        self.assertTrue(result["both_warm_starts_accepted"])
        self.assertTrue(result["warm_starts_same_root"])
        self.assertFalse(result["rest_integration_used"])
        self.assertEqual(len(result["attempts"]), 2)
        self.assertLessEqual(
            result["warm_start_root_distance_normalized"],
            result["root_cluster_relative_tolerance"],
        )


if __name__ == "__main__":
    unittest.main()
