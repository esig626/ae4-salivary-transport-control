"""Scientific invariants for the isolated AE4 cation routing intervention."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import unittest

import numpy as np

from modern_full_model.ae4_routing_only import evaluate_routed_ae4, with_routing_only
from modern_full_model.genotype_evaluation import _root_state_with_basal_regulation
from modern_full_model.model import WT
from modern_full_model.task14_blind import build_model
from modern_full_model.transporters import AE4Environment, AE4Parameters, evaluate_ae4_qss

REPO = Path(__file__).resolve().parents[1]


class RoutingTests(unittest.TestCase):
    def setUp(self):
        self.env = AE4Environment(15.5, 139.5, 50.1, 19., 151.6, 3.4, 124.6, 24.7)
        self.params = AE4Parameters.equal_branch_attempts(
            carrier_amount_fmol=.01, common_cl_attempt_rate_s=10., loaded_attempt_rate_s=10.)

    def test_forward_reverse_zero_and_asymmetric_cations(self):
        rng = np.random.default_rng(20260911)
        environments = [self.env, replace(self.env, hco3_i_mM=.01),
                        AE4Environment(140., 5., 120., 25., 140., 5., 120., 25.)]
        environments += [AE4Environment(*np.exp(rng.uniform(-3., 6., 8))) for _ in range(1000)]
        directions = set()
        for env in environments:
            old = evaluate_ae4_qss(env, self.params)
            new = evaluate_routed_ae4(env, self.params)
            j = new.j_ae4_fmol_s
            directions.add(int(np.sign(j)))
            self.assertGreaterEqual(new.j_na_fmol_s * new.j_k_fmol_s, 0.)
            for field in ('source_cl_i_fmol_s', 'source_hco3_i_fmol_s',
                          'source_tic_i_fmol_s', 'source_alkalinity_i_fmol_s'):
                self.assertEqual(getattr(old, field), getattr(new, field))
            self.assertAlmostEqual(new.source_na_i_fmol_s + new.source_k_i_fmol_s,
                                   -j, delta=2e-14 * max(1., abs(j)))
            self.assertLessEqual(abs(new.transported_charge_equivalents_fmol_s),
                                 2e-14 * max(1., abs(j)))
            side = 'i' if j >= 0. else 'o'
            expected = getattr(env, f'na_{side}_mM') / (getattr(env, f'na_{side}_mM') +
                                                       getattr(env, f'k_{side}_mM'))
            self.assertEqual(new.donor_na_fraction, expected)
        self.assertEqual(directions, {-1, 0, 1})

    def test_zero_capacity_and_zero_regulation_stop_both_cations(self):
        for p, gain in [(replace(self.params, carrier_amount_fmol=0.), 1.), (self.params, 0.)]:
            flux = evaluate_routed_ae4(self.env, p, regulation_gain=gain)
            self.assertEqual(flux.j_ae4_fmol_s, 0.)
            self.assertEqual(flux.j_na_fmol_s, 0.)
            self.assertEqual(flux.j_k_fmol_s, 0.)

    def test_only_na_and_k_amount_rhs_change_at_all_ten_frozen_roots(self):
        manifest = json.loads((REPO / 'results/13B_modern_full_model/native_dynamic_contract_manifest.json').read_text())
        self.assertEqual(len(manifest['roots']), 10)
        for root, payload in manifest['roots'].items():
            old = build_model(manifest, root, .1)
            new = with_routing_only(old)
            self.assertIs(new.parameters, old.parameters)
            self.assertIs(new.ae4_parameters, old.ae4_parameters)
            self.assertIs(new.regulatory_model, old.regulatory_model)
            self.assertIs(new.stimulus, old.stimulus)
            y = _root_state_with_basal_regulation(old, payload['core_state'])
            for t in (0., 60.):
                a, b = old.evaluate(t, y, genotype=WT), new.evaluate(t, y, genotype=WT)
                np.testing.assert_array_equal(a.rhs[2:], b.rhs[2:])
                self.assertAlmostEqual(sum(a.rhs[:2]), sum(b.rhs[:2]), places=13)
                for field in ('membranes', 'homeostasis', 'water', 'observables', 'regulatory'):
                    self.assertEqual(getattr(a.diagnostics, field), getattr(b.diagnostics, field))


if __name__ == '__main__':
    unittest.main()
