"""Synthetic software verification only; no production trajectory or fit."""
import sys
from dataclasses import asdict
from pathlib import Path
import numpy as np
from task51_model import *
from test_vrac_model import VracTests
from modern_full_model.nkcc1_palk2010 import palk_cycle_flux_fmol_s,A1,A2_MM
from modern_full_model.membranes import current_to_fmol_s

class Task51Tests(VracTests):
    def wrapper(self,arm='IPR_ONLY',g=1e-8,interventions=Interventions()):
        return build_task51(arm,beta_only_multiplier=2.5,conductance_S=g,
            reference_volume_pL=self.y[5],source_bath=True,interventions=interventions)

    def test_complete_parent_nesting_two_tracks(self):
        for source in [False,True]:
            y=load_rest('track2' if source else 'track1')
            for arm in ['REST','CCH_ONLY']:
                parent=build_model(arm,source_bath=source)
                m=build_task51(arm,beta_only_multiplier=2.5,conductance_S=1e-8,
                    reference_volume_pL=y[5],source_bath=source)
                for genotype in [GENOTYPES['WT'],GENOTYPES['AE4_KO'],GENOTYPES['AE2_KO']]:
                    a=parent.evaluate(1,y,genotype=genotype);b=m.evaluate(1,y,genotype=genotype)
                    np.testing.assert_array_equal(a.rhs,b.rhs)
                    self.assertEqual(a.diagnostics.homeostasis,b.diagnostics.homeostasis)
                    self.assertEqual(a.diagnostics.water,b.diagnostics.water)
                    self.assertEqual(a.diagnostics.ae4,b.diagnostics.ae4)

    def test_additive_multiplier_and_exact_source_signature(self):
        for arm,total in [('IPR_ONLY',2.5),('CCH_IPR',3.25)]:
            parent=build_model(arm)
            m=self.wrapper(arm,g=0)
            a=parent.evaluate(1,self.y);b=m.evaluate(1,self.y)
            r=b.diagnostics.regulatory
            self.assertEqual(r['nkcc1_capacity_multiplier'],total)
            self.assertEqual(r['nkcc1_beta_increment'],1.5)
            self.assertEqual(b.diagnostics.ae4,a.diagnostics.ae4)
            h0=a.diagnostics.homeostasis;h1=b.diagnostics.homeostasis
            dN=h1.nkcc1_inward_fmol_s-h0.nkcc1_inward_fmol_s
            self.assertGreater(dN,0)
            delta=np.zeros_like(a.rhs);delta[:3]=[dN,dN,2*dN]
            np.testing.assert_allclose(b.rhs-a.rhs,delta,rtol=0,atol=3e-16)
            self.assertEqual(h0.nhe1_inward_fmol_s,h1.nhe1_inward_fmol_s)
            self.assertEqual(h0.ae2_inward_fmol_s,h1.ae2_inward_fmol_s)
            # The VRAC wrapper adds two labelled current keys even when off.
            # Require exact equality of every inherited membrane quantity.
            ma=asdict(a.diagnostics.membranes);mb=asdict(b.diagnostics.membranes)
            mb['currents_A']={k:mb['currents_A'][k] for k in ma['currents_A']}
            self.assertEqual(ma,mb)

    def test_palk_reversal_unchanged(self):
        alpha=self.wrapper().nkcc1_kinetics.alpha_eff_fmol_s
        cl=(A1/A2_MM/(12*115))**.5
        for activity in [1,1.75,2.5,3.25]:
            args=dict(alpha_eff_fmol_s=alpha,activity_multiplier=activity)
            self.assertLess(abs(palk_cycle_flux_fmol_s(12,115,cl,**args)),1e-14)
            self.assertGreater(palk_cycle_flux_fmol_s(12,115,cl*.99,**args),0)
            self.assertLess(palk_cycle_flux_fmol_s(12,115,cl*1.01,**args),0)

    def test_beta_disabled_exact_task49_regression(self):
        for arm in ['IPR_ONLY','CCH_IPR']:
            m=build_task51(arm,beta_only_multiplier=1,conductance_S=1e-8,
                reference_volume_pL=self.y[5],source_bath=True)
            old=VracModel(build_model(arm),conductance_S=1e-8,reference_volume_pL=self.y[5])
            for y in [self.y,self.swollen]:
                np.testing.assert_array_equal(old.rhs(1,y),m.rhs(1,y))
        y=load_rest('track2','AE4_KO');ko=GENOTYPES['AE4_KO']
        m=build_task51('IPR_ONLY',beta_only_multiplier=1,conductance_S=1e-8,
            reference_volume_pL=y[5],source_bath=True)
        r=build_model('REST').evaluate(1,y,genotype=ko)
        e=m.evaluate(1,y,genotype=ko)
        np.testing.assert_array_equal(r.rhs[:12],e.rhs[:12])
        self.assertEqual(e.diagnostics.regulatory['vrac_effective_conductance_S'],0)

    def test_vrac_disabled_keeps_regulated_transport(self):
        a=self.wrapper('CCH_IPR',g=0).evaluate(1,self.swollen)
        b=self.wrapper('CCH_IPR',g=1e-8).evaluate(1,self.swollen)
        self.assertEqual(a.diagnostics.homeostasis,b.diagnostics.homeostasis)
        self.assertEqual(a.diagnostics.ae4,b.diagnostics.ae4)
        self.assertEqual(a.diagnostics.regulatory['nkcc1_capacity_multiplier'],3.25)

    def test_bumetanide_masks_new_nkcc_arm(self):
        m=self.wrapper('CCH_IPR',interventions=Interventions(bumetanide=True))
        self.assertEqual(m.evaluate(1,self.swollen).diagnostics.homeostasis.nkcc1_inward_fmol_s,0)

    def test_total_export_balance_identity(self):
        e=self.wrapper('CCH_IPR').evaluate(1,self.swollen)
        d=e.diagnostics;c=d.membranes.currents_A;F=self.wrapper().parameters.constants.faraday_C_mol
        total=current_to_fmol_s(c['cl_apical_total'],valence=-1,faraday_C_mol=F)
        P=d.membranes.pump_apical_fmol_s+d.membranes.pump_basolateral_fmol_s
        exact=6*P-d.homeostasis.nhe1_inward_fmol_s+2*e.rhs[0]-e.rhs[4]-e.rhs[2]
        self.assertAlmostEqual(total,exact,places=12)

    def test_no_forbidden_model_imports_or_extra_states(self):
        m=self.wrapper()
        self.assertEqual(len(m.state_names),13)
        for name in ['modern_full_model.task50_effective_coupling','modern_full_model.task41_selected','modern_full_model.ae4_cacc_recruitment']:
            self.assertNotIn(name,sys.modules)
