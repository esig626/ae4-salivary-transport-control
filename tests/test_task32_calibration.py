"""WT-only Task 32 structural and control-flow tests; no stationary solves."""
from dataclasses import asdict, replace
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from modern_full_model import task32_wt_chloride_calibration as t


class Task32Tests(unittest.TestCase):
    def test_capacity_scaling_changes_only_existing_carrier_amount(self):
        b = t.load_background('R09')
        amount = 2.339370005697548e-5
        base = t.build_model(b, carrier_amount_fmol=amount, ae4_activity=1)
        scaled = t.build_model(b, carrier_amount_fmol=amount, ae4_activity=.3)
        self.assertEqual(asdict(base.parameters), asdict(scaled.parameters))
        restored = replace(scaled.ae4_parameters,
            carrier_amount_fmol=base.ae4_parameters.carrier_amount_fmol)
        self.assertEqual(asdict(restored), asdict(base.ae4_parameters))
        self.assertEqual(scaled.ae4_parameters.carrier_amount_fmol,
                         .3*base.ae4_parameters.carrier_amount_fmol)
        cached = t.cached_wt(b)
        state = t.state_from_coordinates(base, cached['coordinates'])
        a, c = (m.evaluate(0, state, genotype=t.WT) for m in (base, scaled))
        for field in ('na_cell_fmol_s','k_cell_fmol_s','cl_cell_fmol_s','hco3_cell_fmol_s'):
            self.assertAlmostEqual(getattr(c.diagnostics.ae4,field),
                                   .3*getattr(a.diagnostics.ae4,field), places=15)
        self.assertEqual(asdict(a.diagnostics.homeostasis), asdict(c.diagnostics.homeostasis))
        self.assertEqual(asdict(a.diagnostics.membranes), asdict(c.diagnostics.membranes))
        self.assertEqual(asdict(a.diagnostics.water), asdict(c.diagnostics.water))
        self.assertEqual(t.WT.ae4_expression, 1)
        self.assertEqual(t.AE4_NULL.ae4_expression, 0)

    def test_positive_activity_required(self):
        b = t.load_background('R09')
        for value in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                t.build_model(b, carrier_amount_fmol=1e-5, ae4_activity=value)

    def test_null_gate_has_no_phenotype_targets_and_strict_volume(self):
        r = {'observables':{'na_i_mM':20,'volume_i_pL':2,'ph_i':7.5,'cl_i_mM':80},
             'numerical_pass':True,'finite_capacity':{'pass':True},'exact_ae4_zero':True}
        self.assertTrue(t.acceptance(r))
        self.assertFalse(t.acceptance(r, require_wt_intervals=True))
        r['observables']['volume_i_pL'] = 3
        self.assertFalse(t.acceptance(r))
        r['observables']['volume_i_pL'] = 2
        r['numerical_pass'] = False
        self.assertFalse(t.acceptance(r))

    def test_all_limits_stop_before_new_case(self):
        for field, limit in t.LIMITS.items():
            with self.assertRaises(t.BoundedStop):
                t.Budget(**{field:limit}).check()
        with self.assertRaises(t.BoundedStop):
            t.Budget(numerical_wall_seconds=1800).check()
        with self.assertRaises(t.BoundedStop):
            t.Budget(stationary_residual_evaluations=30000).count_residual()

    def test_scalar_interpolation_requires_actual_sign_bracket(self):
        def row(x):
            return {'label':str(x),'ae4_activity':x,'numerical_pass':True,
                    'observables':{'cl_i_mM':50.1+math.log10(x/.2)}}
        visited=[]
        def evaluate(x, seed):
            visited.append(x)
            return row(x)
        record={'stage':'synthetic'}
        result=t.bracketed_stage(row(1),evaluate,scalar='ae4',target=50.1,
            interval=(50.09999,50.10001),bounds=(1e-6,1e6),record=record)
        self.assertAlmostEqual(visited[0], .1)
        self.assertTrue(record['sign_brackets'])
        self.assertAlmostEqual(result['ae4_activity'], .2)

    def test_unclosed_endpoint_cannot_supply_bracket_sign(self):
        anchor={'label':'a','ae4_activity':1.,'numerical_pass':True,
                'observables':{'cl_i_mM':60.}}
        calls=[]
        def evaluate(x, seed):
            calls.append(x)
            return {**anchor,'label':'failed','ae4_activity':x,'numerical_pass':False}
        with self.assertRaisesRegex(t.BoundedStop,'unclosed WT scalar endpoint'):
            t.bracketed_stage(anchor,evaluate,scalar='ae4',target=50.1,
                interval=t.WT_CL_INTERVAL,bounds=t.AE4_BOUNDS,record={'stage':'test'})
        self.assertEqual(len(calls),1)

    def test_simultaneous_acceptance_skips_nhe1_and_final_ae4(self):
        start={'label':'cached','carrier_amount_fmol':2e-5,'ae4_activity':1.,
               'coordinates':[1]*10,'numerical_pass':True}
        result={**start,'label':'accepted','observables':{'ph_i':6.91,'cl_i_mM':50.1,
                'na_i_mM':12.,'volume_i_pL':1.5},'finite_capacity':{'pass':True},
                'exact_ae4_zero':True}
        with tempfile.TemporaryDirectory() as tmp, patch.object(t,'RESULTS',Path(tmp)), \
             patch.object(t,'cached_wt',return_value=start), \
             patch.object(t,'bracketed_stage',return_value=result) as stage:
            report={}
            self.assertEqual(t.calibrate(None,t.Budget(),report),result)
            self.assertEqual(stage.call_count,1)
            self.assertEqual(len(report['stages']),1)

    def test_signed_counterflux_excluded_from_positive_denominator(self):
        r={'label':'test','background':'synthetic','genotype':'WT','admissible':True,
            'fluxes':{'ae4_cl_cell_fmol_s':2,'nkcc1_cl_in_fmol_s':8,
                      'ae2_cl_in_hco3_out_fmol_s':-1}}
        a=t.allocation(r)
        self.assertEqual(a['ae4_fraction_positive_basolateral'],.2)
        self.assertEqual(a['net_basolateral_cl_fmol_s'],9)
        self.assertIsNone(a['ae4_over_positive_ae2'])

    def test_sequence_allows_only_one_nhe1_and_one_final_ae4(self):
        common={'label':'synthetic','carrier_amount_fmol':2e-5,'ae4_activity':1.,
            'coordinates':[1]*10,'numerical_pass':True,'finite_capacity':{'pass':True},
            'exact_ae4_zero':True}
        def row(ph, cl):
            return {**common,'observables':{'ph_i':ph,'cl_i_mM':cl,
                                           'na_i_mM':12.,'volume_i_pL':1.5}}
        responses=[row(7.1,50.1),row(6.91,54),row(7.1,50.1)]
        with tempfile.TemporaryDirectory() as tmp, patch.object(t,'RESULTS',Path(tmp)), \
             patch.object(t,'cached_wt',return_value=row(6.91,60)), \
             patch.object(t,'bracketed_stage',side_effect=responses) as stage:
            report={}
            with self.assertRaisesRegex(t.BoundedStop,'outside simultaneous WT acceptance'):
                t.calibrate(None,t.Budget(),report)
            self.assertEqual(stage.call_count,3)
            self.assertEqual([r['stage'] for r in report['stages']],
                             ['AE4_initial','NHE1_once','AE4_final_once'])

    def test_failed_calibration_cannot_write_a_parameter_freeze(self):
        r={'observables':{'na_i_mM':12,'volume_i_pL':1.5,'ph_i':6.91,'cl_i_mM':60},
            'numerical_pass':True,'finite_capacity':{'pass':True},'exact_ae4_zero':True}
        with patch.object(t,'write_json') as write:
            with self.assertRaises(t.BoundedStop):
                t.freeze_parameters(None,r,t.Budget())
            write.assert_not_called()

    def test_completed_execution_is_read_only_and_preserves_counters(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(t,'RESULTS',Path(tmp)), \
             patch.dict(t.os.environ,{k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}), \
             patch.object(t,'calibrate') as calibrate:
            t.write_json(Path(tmp)/'execution_summary.json',{'calibration_accepted':False})
            self.assertEqual(t.run_calibration(),{'calibration_accepted':False})
            calibrate.assert_not_called()


if __name__ == '__main__':
    unittest.main()
