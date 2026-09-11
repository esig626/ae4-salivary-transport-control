"""Fixed-state invariants and independent source-bound certificates, no ODEs."""
from dataclasses import asdict
import hashlib
import json
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import task17_fixed_wt as task
from modern_full_model.model import ModernFullModel
from modern_full_model.run_calcium_fast_screen import load_freeze
from modern_full_model.validation import sha256_file, sha256_object


class FixedWTTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest,_=load_freeze()
        cls.payloads=[task.read(p) for p in sorted((task.OUT/'fixed_wt_payloads').glob('*.json'))]

    def test_exact_panel_and_all_187_frozen_inputs(self):
        expected={(r,c) for r in self.manifest['roots'] for c in task.CONDITIONS}
        self.assertEqual({(p['row']['root_id'],p['row']['condition']) for p in self.payloads},expected)
        self.assertEqual(len(self.payloads),30)
        self.assertEqual(len(task.check_inputs()['input_files']),187)
        self.assertEqual(task.read(task.OUT/'contract.json')['target_shares'],[.10,.30])

    def test_every_state_is_bitwise_inherited_without_state_solver(self):
        with patch.object(ModernFullModel,'solve_resting_root',side_effect=AssertionError('WT state solve forbidden')):
            for p in self.payloads:
                rid=p['row']['root_id']
                model,state,_=task.baseline(self.manifest,rid)
                self.assertFalse(state.flags.writeable)
                np.testing.assert_array_equal(p['core_state'],self.manifest['roots'][rid]['core_state'])
                np.testing.assert_array_equal(p['complete_state'],state)
                self.assertEqual(p['complete_state_bytes_sha256'],hashlib.sha256(state.tobytes()).hexdigest())
                self.assertEqual(p['core_state_sha256'],sha256_object(p['core_state']))
                self.assertEqual(p['row']['state_coordinates_changed'],0)

    def test_pathway_ledger_full_rhs_and_electrical_closure(self):
        for rid in self.manifest['roots']:
            model,state,ev=task.baseline(self.manifest,rid)
            ledger=task.source_ledger(model,ev)
            self.assertTrue(set(task.AXES).issubset(ledger))
            self.assertTrue({'AE4_Na_branch','AE4_K_branch','AE4_total','NKCC1','AE2',
                'CaCC','paracellular_Cl','water','lumen_outflow'}.issubset(ledger))
            total=sum((v for v,include in ledger.values() if include),np.zeros(14))
            np.testing.assert_allclose(total[:12],ev.rhs[:12],atol=4e-15,rtol=1e-10)
            np.testing.assert_allclose(total[12:],0,atol=5e-15)
            np.testing.assert_allclose(ledger['AE4_total'][0],
                ledger['AE4_Na_branch'][0]+ledger['AE4_K_branch'][0],atol=0,rtol=0)

    def test_complete_ae4_scaling_and_exact_positive_pool(self):
        for p in self.payloads:
            row=p['row']; model,state,ev=task.baseline(self.manifest,row['root_id'])
            ledger=task.source_ledger(model,ev)
            values,delta=task.allocation(ledger,row['target_share'])
            new=task.allocated_only(model,values)
            got=new.evaluate(0,state)
            self.assertEqual(set(values['positive_other_loaders']),{'NKCC1'})
            self.assertLess(ev.diagnostics.homeostasis.ae2_inward_fmol_s,0)
            for key in ('na_cell_fmol_s','k_cell_fmol_s','cl_cell_fmol_s',
                        'tic_cell_fmol_s','alkalinity_cell_fmol_s'):
                self.assertAlmostEqual(getattr(got.diagnostics.ae4,key),
                    getattr(ev.diagnostics.ae4,key)*values['ae4_capacity_multiplier'],places=13)
            for key in ('affinity_na','affinity_k','effective_rates','occupancy_inward','occupancy_outward'):
                self.assertEqual(got.diagnostics.ae4.diagnostics[key],ev.diagnostics.ae4.diagnostics[key])
            self.assertLess(abs(row['positive_loading_pool_change_fmol_s']),1e-14)
            if row['target_share'] is not None:
                self.assertAlmostEqual(row['realized_ae4_share'],row['target_share'],places=14)
            np.testing.assert_allclose(got.rhs[:12]-ev.rhs[:12],delta[:12],atol=5e-15,rtol=1e-10)
            self.assertEqual(got.diagnostics.membranes,ev.diagnostics.membranes)

    def test_accepted_payloads_satisfy_every_resting_balance(self):
        accepted=[p for p in self.payloads if p['row']['fixed_wt_feasible']]
        self.assertEqual(len(accepted),10)
        for p in accepted:
            model,state,ev=task.baseline(self.manifest,p['row']['root_id'])
            self.assertEqual(p['accepted_whole_cell_parameters'],asdict(model.parameters))
            self.assertEqual(p['accepted_ae4_parameters'],asdict(model.ae4_parameters))
            np.testing.assert_array_equal(p['baseline_rhs'],ev.rhs)
            v=task.residual_summary(ev)
            self.assertLessEqual(v['max_abs_scaled_independent_rhs'],task.SPEC.root_scaled_tolerance)
            self.assertLessEqual(v['max_abs_omitted_rhs'],task.SPEC.omitted_row_raw_tolerance)
            self.assertLessEqual(v['max_abs_state_charge_fmol'],task.SPEC.charge_tolerance_fmol)
            self.assertLessEqual(v['max_conservation_tolerance_ratio'],1)

    def test_rank_and_complete_null_family_against_production_rhs(self):
        for rid in self.manifest['roots']:
            model,state,ev=task.baseline(self.manifest,rid)
            matrix=task.source_matrix(model,task.source_ledger(model,ev))
            info=task.rank_info(matrix)
            self.assertEqual((info['rank'],info['nullity']),(9,1))
            self.assertEqual(len(info['smallest_full_rank_basis']),9)
            # Check the entire remaining free direction in the original model,
            # not merely a copy of the matrix assembly.
            _,_,vh=np.linalg.svd(matrix/np.asarray(info['row_scale'])[:,None])
            multipliers=1+1e-4*vh[-1]/max(abs(vh[-1]))
            new=task.with_effective_multipliers(model,multipliers)
            self.assertAlmostEqual(new.parameters.membranes.g_k_total_S/14e-9,1,places=13)
            got=new.evaluate(0,state)
            np.testing.assert_allclose(got.rhs,ev.rhs,atol=1e-14,rtol=1e-10)
            np.testing.assert_allclose([got.diagnostics.membranes.v_apical_V,got.diagnostics.membranes.v_basolateral_V],
                [ev.diagnostics.membranes.v_apical_V,ev.diagnostics.membranes.v_basolateral_V],atol=1e-14,rtol=0)

    def test_rejections_are_certified_bounds_not_solver_failures(self):
        modified=[p for p in self.payloads if p['row']['target_share'] is not None]
        self.assertEqual(len(modified),20)
        self.assertEqual(sum(p['K_certificate']['required_K_negative'] for p in modified),15)
        for p in modified:
            self.assertFalse(p['row']['fixed_wt_feasible'])
            self.assertIsNone(p['accepted_whole_cell_parameters'])
            self.assertTrue(all(s['status']==2 for s in p['lp_diagnostics']))
            self.assertGreater(p['K_certificate']['K_lower_bound_gap_fmol_s'],.02)
            self.assertLess(p['row']['max_abs_signed_source_equation_residual'],1e-10)
            # Both old channel driving forces are outward. A nonnegative
            # capacity cannot supply the required inward flux in 15 cases.
            self.assertGreater(p['K_certificate']['minimum_K_outward_at_14nS_fmol_s'],0)

    def test_no_frozen_manifest_payload_changes(self):
        path=task.OUT/'frozen_fixed_wt_manifest.json'
        if not path.exists():
            self.skipTest('Freeze follows pre-freeze tests')
        for name,digest in task.read(path)['hashes'].items():
            self.assertEqual(sha256_file(task.REPO/name),digest,name)


if __name__=='__main__':
    unittest.main()
