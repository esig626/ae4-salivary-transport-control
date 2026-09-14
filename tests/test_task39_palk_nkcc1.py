"""Independent Eq. 17/unit/reversal/stoichiometry tests; zero integrations."""
from pathlib import Path
import sys
import math
import unittest
from decimal import Decimal
from dataclasses import replace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis/39_palk_nkcc1_full_validation"))
from validation_common import inherited, resting_algebra, REST_CONCENTRATIONS, REST_CYCLE
from modern_full_model.nkcc1_palk2010 import (
    A1, A2_MM, A3, A4_MM, PALK_NKCC1, Nkcc1Kinetics,
    palk_shape_factor, palk_cycle_flux_fmol_s,
)
from modern_full_model.membranes import HomeostasisEnvironment, evaluate_homeostasis


class PalkSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, _, cls.rest, cls.stim, cls.y0 = inherited.models_and_state()
        e = cls.rest.evaluate(0., cls.y0)
        obs, p = e.diagnostics.observables, cls.rest.parameters
        ci = obs.cell_concentrations_mM
        cls.env = HomeostasisEnvironment(ci["na"], ci["k"], ci["cl"],
            obs.cell_acid_base.h_mM, ci["hco3"], p.bath.na_mM, p.bath.k_mM,
            p.bath.cl_mM, obs.bath_acid_base.h_mM, obs.bath_acid_base.hco3_mM)
        cls.alpha = resting_algebra()["alpha_eff_fmol_s"]
        cls.selection = Nkcc1Kinetics(PALK_NKCC1, cls.alpha)

    def test_literal_eq17_and_algebraic_rest_scale(self):
        na, k, cl = REST_CONCENTRATIONS
        literal = (157.5 - 0.000020096*na*k*cl*cl)/(1.0306 + 0.0000013852*na*k*cl*cl)
        self.assertTrue(math.isclose(palk_shape_factor(na,k,cl), literal, rel_tol=2e-15))
        audit = resting_algebra()
        for key, expected in (("x_rest_mM4",4941065.334966328),
                              ("shape_rest",7.391062769440939),
                              ("alpha_eff_fmol_s",0.017334746894096052)):
            self.assertTrue(math.isclose(audit[key],expected,rel_tol=2e-15),key)
        self.assertAlmostEqual(palk_cycle_flux_fmol_s(na,k,cl,
            alpha_eff_fmol_s=self.alpha,activity_multiplier=1),REST_CYCLE,delta=1e-15)

    def test_M_to_mM_fourth_order_conversion(self):
        self.assertEqual(Decimal(str(A2_MM)),Decimal('2.0096e7')*Decimal('1e-12'))
        self.assertEqual(Decimal(str(A4_MM)),Decimal('1.3852e6')*Decimal('1e-12'))
        self.assertTrue(math.isclose(A2_MM,2.0096e7*1e-12,rel_tol=2e-15))
        self.assertTrue(math.isclose(A4_MM,1.3852e6*1e-12,rel_tol=2e-15))
        na,k,cl = REST_CONCENTRATIONS
        xm = (na*1e-3)*(k*1e-3)*(cl*1e-3)**2
        source = (157.5-2.0096e7*xm)/(1.0306+1.3852e6*xm)
        self.assertTrue(math.isclose(source,palk_shape_factor(na,k,cl),rel_tol=2e-15))

    def test_reversal_and_unclipped_reverse(self):
        na,k,_ = REST_CONCENTRATIONS
        cl = math.sqrt((157.5/2.0096e-5)/(na*k))
        self.assertLess(abs(palk_shape_factor(na,k,cl)),1e-13)
        self.assertGreater(palk_shape_factor(na,k,cl*0.99),0)
        self.assertLess(palk_shape_factor(na,k,cl*1.01),0)

    def test_activity_scales_complete_flux_preserving_reversal(self):
        na,k,cl = REST_CONCENTRATIONS
        reversal = math.sqrt((157.5/2.0096e-5)/(na*k))
        j = palk_cycle_flux_fmol_s(na,k,cl,alpha_eff_fmol_s=self.alpha,activity_multiplier=1)
        recruited = palk_cycle_flux_fmol_s(na,k,cl,alpha_eff_fmol_s=self.alpha,activity_multiplier=1.75)
        self.assertAlmostEqual(recruited,1.75*j,delta=1e-15)
        self.assertLess(abs(palk_cycle_flux_fmol_s(na,k,reversal,
            alpha_eff_fmol_s=self.alpha,activity_multiplier=1.75)),1e-14)
        self.assertLess(palk_cycle_flux_fmol_s(na,k,reversal*1.01,
            alpha_eff_fmol_s=self.alpha,activity_multiplier=1.75),0)

    def test_actual_homeostasis_source_signature_both_directions(self):
        na,k,_ = REST_CONCENTRATIONS
        reverse_cl = math.sqrt((157.5/2.0096e-5)/(na*k))*1.01
        for env in (self.env, replace(self.env,cl_i_mM=reverse_cl)):
            h = evaluate_homeostasis(env,self.rest.parameters,nkcc1_kinetics=self.selection)
            j = h.nkcc1_inward_fmol_s
            isolated = [h.na_cell_fmol_s-h.nhe1_inward_fmol_s,h.k_cell_fmol_s,
                h.cl_cell_fmol_s-h.ae2_inward_fmol_s,
                h.tic_cell_fmol_s+h.ae2_inward_fmol_s,
                h.alkalinity_cell_fmol_s-h.nhe1_inward_fmol_s+h.ae2_inward_fmol_s]
            for actual,expected in zip(isolated,[j,j,2*j,0,0]):
                self.assertAlmostEqual(actual,expected,delta=1e-15)
            self.assertLess(abs(h.charge_source_residual_fmol_s),1e-15)

    def test_N1_inherited_rest_and_stimulus_exact(self):
        rest = self.rest.evaluate(0,self.y0).diagnostics.regulatory
        stim = self.stim.evaluate(1e-6,self.y0).diagnostics.regulatory
        self.assertEqual(rest["nkcc1_capacity_multiplier"],1)
        self.assertEqual(stim["nkcc1_capacity_multiplier"],1.75)

    def test_generic_comparator_retained_without_other_flux_changes(self):
        old = evaluate_homeostasis(self.env,self.rest.parameters)
        new = evaluate_homeostasis(self.env,self.rest.parameters,nkcc1_kinetics=self.selection)
        self.assertAlmostEqual(old.nkcc1_inward_fmol_s,REST_CYCLE,delta=1e-15)
        self.assertEqual(old.nhe1_inward_fmol_s,new.nhe1_inward_fmol_s)
        self.assertEqual(old.ae2_inward_fmol_s,new.ae2_inward_fmol_s)
        self.assertEqual(old.nhe1_model,new.nhe1_model)


if __name__ == "__main__":
    unittest.main(verbosity=2)
