"""New continuation-only tests; no stationary solve or prior test is rerun."""
from dataclasses import asdict, replace
from pathlib import Path
import json
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from modern_full_model.task31_nhe1_repair import (
    Budget, BudgetStop, build_model, load_background, solve_stationary,
    task31_parameters, WT, AE4_NULL, NHE1_CHA_2009,
)


class ContinuationTests(unittest.TestCase):
    def test_only_selector_and_single_amount_change_in_both_backgrounds(self):
        for alias in ('R09', 'R10'):
            b = load_background(alias)
            p = task31_parameters(b, carrier_amount_fmol=1e-4)
            restored = replace(p, homeostasis=replace(p.homeostasis,
                nhe1_model=b.inherited_parameters.homeostasis.nhe1_model,
                nhe1_cha_carrier_amount_fmol=b.inherited_parameters.homeostasis.nhe1_cha_carrier_amount_fmol))
            self.assertEqual(asdict(restored), asdict(b.inherited_parameters))
            model = build_model(b, carrier_amount_fmol=1e-4)
            self.assertEqual(model.parameters.homeostasis.nhe1_model, NHE1_CHA_2009)
            self.assertEqual(AE4_NULL.ae4_expression, 0)

    def test_global_counters_stop_before_exceeding_limits(self):
        limits = {'stationary_solver_calls': 16, 'stationary_residual_evaluations': 24000,
                  'scalar_density_evaluations': 6, 'stimulated_integrations': 4}
        for key, limit in limits.items():
            b = Budget(**{key: limit-1})
            b.count(key, limit)
            with self.assertRaises(BudgetStop):
                b.count(key, limit)
            self.assertEqual(getattr(b, key), limit)
        with self.assertRaises(BudgetStop):
            Budget(numerical_wall_seconds=1801).check_time()

    def test_checkpoint_return_never_calls_stationary_optimizer(self):
        b = load_background('R09')
        model = build_model(b, carrier_amount_fmol=1e-4)
        budget = Budget()
        with patch('pathlib.Path.exists', return_value=True), \
             patch('pathlib.Path.read_text', return_value=json.dumps({'checkpoint': 'accepted'})), \
             patch('modern_full_model.task31_nhe1_repair.least_squares') as optimizer:
            result = solve_stationary(model, genotype=WT, start_coordinates=b.saved_coordinates,
                background=b, label='cache', budget=budget)
        self.assertEqual(result, {'checkpoint': 'accepted'})
        optimizer.assert_not_called()
        self.assertEqual(budget.stationary_solver_calls, 0)

    def test_exhausted_wall_time_can_still_be_checkpointed(self):
        budget = Budget(numerical_wall_seconds=1801)
        with patch('modern_full_model.task31_nhe1_repair._write_json') as write:
            budget.save()
        self.assertFalse(write.call_args.args[1]['limits_respected'])
        self.assertGreaterEqual(write.call_args.args[1]['numerical_wall_seconds'], 1801)


if __name__ == '__main__':
    unittest.main()
