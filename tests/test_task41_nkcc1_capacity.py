"""Whole-cycle capacity, reversal, invariance and algebraic-design checks."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis/41_ae4_loss_algebraic_design'))
from validation_common import *
import unittest
from dataclasses import replace,asdict
from modern_full_model.nkcc1_capacity_limited import capacity_limited_cycles


class CapacityTests(unittest.TestCase):
    def test_signed_cap_and_nkcc1_expression_scaling(self):
        for demand in (-1.,-.05,0.,.05,1.):
            j=capacity_limited_cycles(demand,.118)
            self.assertLessEqual(abs(j),.118)
            self.assertEqual(capacity_limited_cycles(-demand,.118),-j)
            self.assertAlmostEqual(capacity_limited_cycles(.2*demand,.118,.2),.2*j,delta=1e-16)
            self.assertEqual(capacity_limited_cycles(demand,.118,0.),0.)

    def test_whole_model_isolated_source_and_nonrouting_invariance(self):
        rest,stim,newrest,newstim,y,_=models_and_reference()
        for base,new,t in ((rest,newrest,0.),(stim,newstim,1e-6)):
            for expression in (1.,.05,0.):
                genotype=replace(WT,ae4_expression=expression)
                a,b=base.evaluate(t,y,genotype=genotype),new.evaluate(t,y,genotype=genotype)
                dn=b.diagnostics.homeostasis.nkcc1_inward_fmol_s-a.diagnostics.homeostasis.nkcc1_inward_fmol_s
                np.testing.assert_allclose(b.rhs-a.rhs,[dn,dn,2*dn]+[0.]*10,rtol=0,atol=1e-15)
                for field in ('ae4','membranes','water','co2_fluxes_fmol_s','observables'):
                    self.assertEqual(getattr(a.diagnostics,field),getattr(b.diagnostics,field))
                self.assertLess(abs(b.rhs[0]+b.rhs[1]-b.rhs[2]-b.rhs[4]),1e-10)
            self.assertEqual(old.parameter_hashes(base),old.parameter_hashes(new))

    def test_production_matches_selected_inverse_equations(self):
        _,_,_,m,_,_=models_and_reference()
        design=read_json(OUT/'capacity_design.json')
        for state in design['states']:
            g=WT if state['case']=='wt' else replace(WT,ae4_expression=0.)
            e=m.evaluate(600.,np.array(state['state_vector']),genotype=g)
            np.testing.assert_allclose(e.rhs,state['raw_rhs'],rtol=0,atol=1e-15)
            self.assertLess(np.max(np.abs(e.rhs[:12])),1e-10)
            d=e.diagnostics;h=d.homeostasis;r=d.regulatory
            # Exact transient identity: CaCC=6P-H+2dotNa-dotTA-dotCl (r=1/2).
            row=old.parent.inherited.row_for(m,600.,np.array(state['state_vector']),e)
            pump=d.membranes.pump_apical_fmol_s+d.membranes.pump_basolateral_fmol_s
            expected=6*pump-h.nhe1_inward_fmol_s+2*e.rhs[0]-e.rhs[4]-e.rhs[2]
            self.assertAlmostEqual(row['cacc_cl_export_fmol_s'],expected,delta=1e-12)


if __name__=='__main__':unittest.main(verbosity=2)
