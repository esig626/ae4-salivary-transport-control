"""Focused tests of new current insertion; no production trajectories."""
import json
from dataclasses import asdict
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from scipy.optimize import root
from vrac_model import *
from protocol_layer import Interventions
from modern_full_model.membranes import nernst_voltage_V
from modern_full_model.vbeta_diagnostic import evaluate_vbeta_current

class VracTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.y=np.array(json.loads((ROOT/'analysis/48_joint_experimental_constraint_reconstruction/output/rests/WT.json').read_text())['state'])
        # Synthetic 10% swelling: amounts stay fixed. This is only a test fixture.
        cls.swollen=cls.y.copy();cls.swollen[5]*=1.1

    def wrapper(self,arm='IPR_ONLY',g=1e-8,interventions=Interventions()):
        return VracModel(build_model(arm,interventions=interventions),conductance_S=g,reference_volume_pL=self.y[5])

    def test_rest_and_cch_exact_nesting(self):
        for arm,y in [('REST',self.y),('CCH_ONLY',self.swollen)]:
            m=self.wrapper(arm)
            old=m.parent.evaluate(1,y);new=m.evaluate(1,y)
            np.testing.assert_array_equal(old.rhs,new.rhs)
            self.assertEqual(new.diagnostics.regulatory['vrac_current_A'],0)
            self.assertEqual(old.diagnostics.water,new.diagnostics.water)

    def test_zero_scale_and_closed_volume_gate(self):
        for g,y in [(0,self.swollen),(1e-8,self.y)]:
            m=self.wrapper(g=g)
            np.testing.assert_array_equal(m.rhs(1,y),m.parent.rhs(1,y))
        m=self.wrapper();y=self.y.copy();y[5]*=.99
        self.assertEqual(m.evaluate(1,y).diagnostics.regulatory['vrac_effective_conductance_S'],0)

    def test_fixed_voltage_reversal_and_chloride_sign(self):
        m=self.wrapper();c=m.parameters.constants
        ec=nernst_voltage_V(50,130,valence=-1,thermal_voltage_V=c.thermal_voltage_V)
        args=dict(parameters=m.vrac_parameters,beta_input=1,cell_volume_pL=self.swollen[5],
            chloride_cell_mM=50,chloride_lumen_mM=130,
            thermal_voltage_V=c.thermal_voltage_V,faraday_C_mol=c.faraday_C_mol)
        at=evaluate_vbeta_current(v_apical_V=ec,**args)
        out=evaluate_vbeta_current(v_apical_V=ec-.01,**args)
        self.assertEqual(at.current_A,0)
        self.assertLess(out.chloride_cell_source_fmol_s,0)
        self.assertEqual(out.chloride_conservation_residual_fmol_s,0)
        self.assertAlmostEqual(out.charge_conversion_residual_fmol_s,0,places=14)

    def test_current_and_amount_conservation_with_nbc(self):
        m=self.wrapper('CCH_IPR');old=m.parent.evaluate(1,self.swollen);new=m.evaluate(1,self.swollen)
        d=new.diagnostics
        self.assertGreater(d.regulatory['vrac_chloride_export_fmol_s'],0)
        self.assertLess(max(map(abs,d.membranes.current_residuals_A.values())),1e-20)
        self.assertLess(max(map(abs,d.membranes.charge_residuals_fmol_s.values())),1e-10)
        for key in ['cell_bulk_charge_rate_fmol_s','lumen_bulk_charge_rate_minus_outflow_fmol_s','carbon_accounting_fmol_s']:
            self.assertLess(abs(d.conservation_residuals[key]),1e-10)
        delta=new.rhs-old.rhs
        # Intracellular/luminal Cl increments cancel except changed para return.
        ip=d.membranes.currents_A['para_cl']-old.diagnostics.membranes.currents_A['para_cl']
        self.assertAlmostEqual(delta[2]+delta[8],-current_to_fmol_s(ip,valence=-1,faraday_C_mol=m.parameters.constants.faraday_C_mol),places=13)
        self.assertNotEqual(d.regulatory['minimal_nbc_current_A'],old.diagnostics.regulatory['minimal_nbc_current_A'])
        self.assertEqual(d.water,old.diagnostics.water)
        self.assertEqual(d.homeostasis,old.diagnostics.homeostasis)
        self.assertEqual(d.ae4,old.diagnostics.ae4)
        np.testing.assert_array_equal(new.rhs[[5,11,12]],old.rhs[[5,11,12]])

    def test_direct_two_voltage_equations(self):
        m=self.wrapper('CCH_IPR');d=m.evaluate(1,self.swollen).diagnostics
        p=m.parameters;mp=p.membranes;c=p.constants;ci=d.observables.cell_concentrations_mM;li=d.observables.lumen_concentrations_mM
        h=hill_activation(d.stimulus.calcium_uM,mp.calcium_half_uM,mp.calcium_hill)
        ka=mp.g_k_total_S*h*mp.apical_k_fraction;kb=mp.g_k_total_S*h*(1-mp.apical_k_fraction)
        cl=mp.g_cl_apical_S*h;gv=d.regulatory['vrac_effective_conductance_S']
        ek_a=nernst_voltage_V(ci['k'],li['k'],valence=1,thermal_voltage_V=c.thermal_voltage_V)
        ek_b=nernst_voltage_V(ci['k'],p.bath.k_mM,valence=1,thermal_voltage_V=c.thermal_voltage_V)
        ec=nernst_voltage_V(ci['cl'],li['cl'],valence=-1,thermal_voltage_V=c.thermal_voltage_V)
        paths=[(mp.g_para_na_S,1,li['na'],p.bath.na_mM),(mp.g_para_k_S,1,li['k'],p.bath.k_mM),
               (mp.g_para_cl_S,-1,li['cl'],p.bath.cl_mM),(mp.g_para_hco3_S,-1,li['hco3'],d.observables.bath_acid_base.hco3_mM)]
        def direct(v):
            va,vb=v;ip=sum(g*(vb-va-nernst_voltage_V(a,b,valence=z,thermal_voltage_V=c.thermal_voltage_V)) for g,z,a,b in paths)
            nb=evaluate_minimal_nbc(na_i_mM=ci['na'],hco3_i_mM=ci['hco3'],na_o_mM=p.bath.na_mM,
                hco3_o_mM=d.observables.bath_acid_base.hco3_mM,v_basolateral_V=vb,
                thermal_voltage_V=c.thermal_voltage_V,faraday_C_mol=c.faraday_C_mol,activation_fraction=1,parameters=m.nbc_parameters)
            return np.array([ka*(va-ek_a)+(cl+gv)*(va-ec)+d.membranes.currents_A['pump_apical']-ip,
                kb*(vb-ek_b)+d.membranes.currents_A['pump_basolateral']+nb.current_A+ip])*1e9
        sol=root(direct,[-.05,-.05],tol=1e-11)
        self.assertLess(max(abs(direct(sol.x))),1e-10)
        np.testing.assert_allclose(sol.x,[d.membranes.v_apical_V,d.membranes.v_basolateral_V],atol=1e-11,rtol=0)

    def test_tmem16a_independence(self):
        m=self.wrapper(interventions=Interventions(t16ainh_a01=True))
        d=m.evaluate(1,self.swollen).diagnostics
        self.assertEqual(d.membranes.currents_A['cl_apical'],0)
        self.assertNotEqual(d.membranes.currents_A['vrac_apical'],0)

    def test_no_direct_ae4_expression_factor(self):
        m=self.wrapper();wt=m.evaluate(1,self.swollen);ko=m.evaluate(1,self.swollen,genotype=GENOTYPES['AE4_KO'])
        self.assertEqual(wt.diagnostics.regulatory['vrac_effective_conductance_S'],ko.diagnostics.regulatory['vrac_effective_conductance_S'])
        self.assertEqual(wt.diagnostics.regulatory['vrac_current_A'],ko.diagnostics.regulatory['vrac_current_A'])

    def test_unchanged_parameter_objects_and_pump(self):
        m=self.wrapper('CCH_IPR');before=asdict(m.parameters)
        old=m.parent.evaluate(1,self.swollen);new=m.evaluate(1,self.swollen)
        self.assertEqual(before,asdict(m.parameters))
        for k in ['pump_apical_fmol_s','pump_basolateral_fmol_s']:
            self.assertEqual(getattr(old.diagnostics.membranes,k),getattr(new.diagnostics.membranes,k))
        self.assertEqual(old.diagnostics.regulatory['nkcc1_capacity_multiplier'],new.diagnostics.regulatory['nkcc1_capacity_multiplier'])

    def test_integration_dispatch_uses_wrapper(self):
        m=self.wrapper()
        with patch('vrac_model.solve_ivp') as solver:
            m.solve_dynamics((0,1),initial_state=self.swollen)
            rhs=solver.call_args.args[0]
            np.testing.assert_array_equal(rhs(1,self.swollen),m.rhs(1,self.swollen))

if __name__=='__main__':unittest.main(verbosity=2)
