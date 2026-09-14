"""Recruitment boundary, closure, nesting, and inverse-equation verification."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis/41_ae4_loss_algebraic_design/candidate_07'))
from candidate_common import *
from modern_full_model.ae4_cacc_recruitment import recruitment_factor
from dataclasses import replace,asdict
import unittest


class RecruitmentTests(unittest.TestCase):
    def test_rest_wt_and_expression_law(self):
        b=read_json(OUT/'design.json')['null_stimulated_cacc_fraction']
        for e in (0.,.05,1.):
            self.assertEqual(recruitment_factor(0.,e,b),1.)
            self.assertAlmostEqual(recruitment_factor(1.,e,b),b+(1-b)*e)
        for a in (0.,.3,1.):self.assertEqual(recruitment_factor(a,1.,b),1.)

    def test_wt_exact_nesting_and_no_output_scaling(self):
        r,s,nr,ns,y,_=models_and_reference()
        for base,new,t in ((r,nr,0.),(s,ns,1e-6),(s,ns,600.)):
            a,b=base.evaluate(t,y),new.evaluate(t,y)
            np.testing.assert_array_equal(a.rhs,b.rhs)
            self.assertEqual(a.diagnostics.membranes,b.diagnostics.membranes)
            self.assertEqual(old.parameter_hashes(base),old.parameter_hashes(new))
        for expression in (0.,.05):
            g=replace(WT,ae4_expression=expression)
            a,b=s.evaluate(1e-6,y,genotype=g),ns.evaluate(1e-6,y,genotype=g)
            self.assertEqual(a.diagnostics.water,b.diagnostics.water)
            self.assertEqual(a.diagnostics.ae4,b.diagnostics.ae4)
            self.assertEqual(a.diagnostics.homeostasis,b.diagnostics.homeostasis)
            f=b.diagnostics.regulatory['ae4_dependent_cacc_recruitment_factor']
            ep=asdict(ns.effective_model(f).parameters);pp=asdict(s.parameters)
            ep['membranes']['g_cl_apical_S']=pp['membranes']['g_cl_apical_S']
            self.assertEqual(ep,pp)

    def test_inverse_rhs_current_and_chloride_identities(self):
        _,_,_,m,_,_=models_and_reference()
        for s in read_json(OUT/'design.json')['states']:
            g=WT if s['case']=='wt' else replace(WT,ae4_expression=0.)
            y=np.array(s['state_vector']);e=m.evaluate(600.,y,genotype=g)
            np.testing.assert_allclose(e.rhs,s['rhs'],rtol=0,atol=2e-15)
            row,failures,ratios=diagnose(m,600.,y,e)
            self.assertTrue(all(r<=1 for r in ratios.values()))
            P=e.diagnostics.membranes.pump_apical_fmol_s+e.diagnostics.membranes.pump_basolateral_fmol_s
            H=e.diagnostics.homeostasis.nhe1_inward_fmol_s
            self.assertAlmostEqual(row['cacc_cl_export_fmol_s'],6*P-H+2*e.rhs[0]-e.rhs[4]-e.rhs[2],delta=1e-12)
            rhs=2*e.diagnostics.homeostasis.nkcc1_inward_fmol_s+e.diagnostics.ae4.cl_cell_fmol_s+e.diagnostics.homeostasis.ae2_inward_fmol_s-row['paracellular_cl_return_fmol_s']-e.rhs[2]-e.rhs[8]
            self.assertAlmostEqual(e.diagnostics.water.lumen_outflow_pL_s*e.diagnostics.observables.lumen_concentrations_mM['cl'],rhs,delta=1e-12)


if __name__=='__main__':unittest.main(verbosity=2)
