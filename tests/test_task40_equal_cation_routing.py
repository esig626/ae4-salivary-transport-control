"""Focused Task 40 source tests. No stationary solve or integration."""
from pathlib import Path
import sys
import unittest
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
import math

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis/40_ae4_equal_cation_routing"))
from validation_common import (models_and_reference, verify_parent_sources,
    parameter_hashes, parent, np)
from modern_full_model.model import WT
from modern_full_model.ae4_routing_only import evaluate_routed_ae4, routing_only_adapter
from modern_full_model.ae4_equal_cation_routing import (
    EqualRoutedAE4Flux, evaluate_equal_routed_ae4, equal_routing_adapter)
from modern_full_model.transporters import AE4Environment


def reference_preflight(old_rest, equal_rest, y):
    old, new = old_rest.evaluate(0., y), equal_rest.evaluate(0., y)
    c = old.diagnostics.observables.cell_concentrations_mM
    j = old.diagnostics.ae4.cl_cell_fmol_s
    fraction = c['na']/(c['na']+c['k'])
    with localcontext() as ctx:
        ctx.prec = 50
        na, k = Decimal('11.636125748680639'), Decimal('116.76320932871538')
        jd = Decimal('0.004937067441546469')
        fd = na/(na+k)
        delta_decimal = -(Decimal('0.5')-fd)*jd
    fields = ('na_cell_fmol_s','k_cell_fmol_s','cl_cell_fmol_s',
              'tic_cell_fmol_s','alkalinity_cell_fmol_s')
    shifts = {key: getattr(new.diagnostics.ae4,key)-getattr(old.diagnostics.ae4,key)
              for key in fields}
    return {'donor_na_fraction': fraction, 'J4_fmol_s': j,
            'decimal_na_fraction': str(fd), 'decimal_delta_na_source_fmol_s': str(delta_decimal),
            'expected_delta_na_source_fmol_s': -(0.5-fraction)*j,
            'expected_delta_k_source_fmol_s': (0.5-fraction)*j,
            'actual_ae4_source_shifts_fmol_s': shifts,
            'full_rhs_shift': (np.array(new.rhs)-old.rhs).tolist(),
            'reference_old_rhs': list(map(float,old.rhs)),
            'reference_equal_rhs': list(map(float,new.rhs)),
            'charge_perturbation_fmol_s': math.fsum((shifts[fields[0]],shifts[fields[1]],
                -shifts[fields[2]],-shifts[fields[4]])),
            'exact_task31_rest_nesting_expected': False}


class EqualCationSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.saved, cls.old_rest, cls.old_stim, cls.rest, cls.stim, cls.y = models_and_reference()
        obs = cls.old_rest.evaluate(0., cls.y).diagnostics.observables
        ci, bath = obs.cell_concentrations_mM, cls.rest.parameters.bath
        cls.env = AE4Environment(ci['na'],ci['k'],ci['cl'],ci['hco3'],
            bath.na_mM,bath.k_mM,bath.cl_mM,obs.bath_acid_base.hco3_mM)
        v = asdict(cls.env)
        cls.reverse_env = AE4Environment(**{f'{ion}_{side}_mM':v[f'{ion}_{other}_mM']
            for ion in ('na','k','cl','hco3') for side,other in (('i','o'),('o','i'))})

    def test_forward_reverse_and_zero_literal_signature(self):
        for j in (0.004937067441546469, -0.004937067441546469, 0.):
            a = EqualRoutedAE4Flux(j,1.)
            self.assertEqual(list(a.intracellular_conserved_sources_fmol_s.values()),
                             [-0.5*j,-0.5*j,j,-2*j,-2*j])
            self.assertEqual(a.j_na_fmol_s,a.j_k_fmol_s)
            self.assertEqual(a.source_na_i_fmol_s+a.source_k_i_fmol_s,-j)
            self.assertEqual(a.transported_charge_equivalents_fmol_s,0.)

    def test_inherited_j4_and_anion_sources_exact_in_both_directions(self):
        for env,sign in ((self.env,1),(self.reverse_env,-1)):
            a = evaluate_equal_routed_ae4(env,self.rest.ae4_parameters)
            old = evaluate_routed_ae4(env,self.rest.ae4_parameters)
            self.assertGreater(sign*a.j_ae4_fmol_s,0)
            self.assertEqual(a.j_ae4_fmol_s,old.j_ae4_fmol_s)
            for key in ('cl','hco3','tic','alkalinity'):
                self.assertEqual(getattr(a,f'source_{key}_i_fmol_s'),
                                 getattr(old,f'source_{key}_i_fmol_s'))
            self.assertEqual(a.source_na_i_fmol_s,a.source_k_i_fmol_s)
            self.assertEqual(a.source_na_i_fmol_s,-0.5*old.j_ae4_fmol_s)

    def test_complete_cycle_regulatory_scaling_and_expression(self):
        for env in (self.env,self.reverse_env):
            base = evaluate_equal_routed_ae4(env,self.rest.ae4_parameters)
            for gain in (0.,1.,1.25):
                new = evaluate_equal_routed_ae4(env,self.rest.ae4_parameters,regulation_gain=gain)
                old = evaluate_routed_ae4(env,self.rest.ae4_parameters,regulation_gain=gain)
                self.assertEqual(new.j_ae4_fmol_s,old.j_ae4_fmol_s)
                for key,value in base.intracellular_conserved_sources_fmol_s.items():
                    # The inherited QSS law subtracts microscopic currents;
                    # linear capacity scaling therefore has floating-point
                    # rounding. Identity to the inherited J4 above is exact.
                    self.assertTrue(math.isclose(
                        new.intracellular_conserved_sources_fmol_s[key],gain*value,
                        rel_tol=32*np.finfo(float).eps,
                        abs_tol=32*np.finfo(float).eps*abs(base.j_ae4_fmol_s)))
        for expression in (1.,0.05,0.):
            genotype=replace(WT,ae4_expression=expression)
            a=self.rest.evaluate(0.,self.y,genotype=genotype).diagnostics.ae4
            old=self.old_rest.evaluate(0.,self.y,genotype=genotype).diagnostics.ae4
            self.assertEqual(a.cl_cell_fmol_s,old.cl_cell_fmol_s)
            self.assertEqual(a.na_cell_fmol_s,a.k_cell_fmol_s)
            self.assertEqual(a.na_cell_fmol_s,-0.5*a.cl_cell_fmol_s)
            self.assertEqual(a.tic_cell_fmol_s,-2*a.cl_cell_fmol_s)

    def test_reference_algebra_and_full_rhs_perturbation(self):
        p=reference_preflight(self.old_rest,self.rest,self.y)
        self.assertAlmostEqual(p['donor_na_fraction'],0.09062450161183984,delta=2e-16)
        self.assertEqual(p['J4_fmol_s'],0.004937067441546469)
        delta=-0.0020211144444590447
        self.assertAlmostEqual(float(p['decimal_delta_na_source_fmol_s']),delta,delta=2e-18)
        self.assertAlmostEqual(p['expected_delta_na_source_fmol_s'],delta,delta=2e-18)
        np.testing.assert_allclose(p['full_rhs_shift'],[delta,-delta]+[0.]*11,rtol=0,atol=2e-16)
        self.assertLess(abs(p['charge_perturbation_fmol_s']),2e-18)
        self.assertEqual(p['actual_ae4_source_shifts_fmol_s']['cl_cell_fmol_s'],0.)
        self.assertEqual(p['actual_ae4_source_shifts_fmol_s']['tic_cell_fmol_s'],0.)
        self.assertEqual(p['actual_ae4_source_shifts_fmol_s']['alkalinity_cell_fmol_s'],0.)

    def test_no_other_mechanism_or_parameter_changes(self):
        verify_parent_sources()
        for old,new in ((self.old_rest,self.rest),(self.old_stim,self.stim)):
            self.assertEqual(parameter_hashes(old),parameter_hashes(new))
            self.assertIs(old.ae4_evaluator,routing_only_adapter)
            self.assertIs(new.ae4_evaluator,equal_routing_adapter)
            for t,y in ((0.,self.y),(1e-6,np.r_[self.y[:12],0.4])):
                a,b=old.evaluate(t,y),new.evaluate(t,y)
                for key in ('homeostasis','membranes','water','co2_fluxes_fmol_s','observables','regulatory'):
                    self.assertEqual(getattr(a.diagnostics,key),getattr(b.diagnostics,key),key)
                self.assertEqual(a.diagnostics.ae4.cl_cell_fmol_s,b.diagnostics.ae4.cl_cell_fmol_s)
                np.testing.assert_array_equal(a.rhs[2:],b.rhs[2:])
        frozen=json_load(ROOT/'results/39_palk_nkcc1_full_validation/frozen_inputs.json')
        self.assertEqual(parameter_hashes(self.stim),frozen['parameter_hashes'])


def json_load(path):
    import json
    return json.loads(path.read_text())


if __name__=='__main__':
    unittest.main(verbosity=2)
