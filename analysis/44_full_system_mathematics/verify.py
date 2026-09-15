"""Regression checks for Task 44 analysis, independent of manuscript wording."""
from core import *
from modern_full_model.acid_base import _carbon_fractions,total_alkalinity_mM
import unittest
OUT=Path(__file__).resolve().parent/'output'

def read(name):return json.loads((OUT/name).read_text())
class AnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _,s,cls.y0,_,_=reference();cls.m=clone(s);cls.sc=Scaling(cls.m,cls.y0)
    def test_scale_inverse(self):
        self.assertLess(np.max(abs(self.sc.decode(self.sc.encode(self.y0))-self.y0)),1e-12)
    def test_charge_coordinates_are_exact(self):
        self.assertLess(np.max(abs(self.sc.expand(self.sc.independent(self.y0))-self.y0)),1e-10)
        f=self.m.rhs(1.,self.y0)
        self.assertLess(abs(f[0]+f[1]-f[2]-f[4]),1e-10)
        self.assertLess(abs(f[6]+f[7]-f[8]-f[10]),1e-10)
    def test_carbonate_derivative(self):
        p=self.m.parameters.acid_base
        for ph in np.linspace(6.6,7.3,15):
            T,B=8.,20.;a=np.array(_carbon_fractions(ph,p));k=np.arange(3.)
            b=1/(1+10**(p.cell_buffer_pka-ph));H=1e3*10**(-ph);OH=1e3*10**(ph-p.water_pkw)
            deriv=np.log(10)*(T*(a@k**2-(a@k)**2)+B*b*(1-b)+H+OH)
            fun=lambda ph:total_alkalinity_mM(ph=ph,total_carbon_mM=T,buffer_total_mM=B,buffer_pka=p.cell_buffer_pka,parameters=p)
            numerical=(fun(ph+1e-5)-fun(ph-1e-5))/2e-5
            self.assertGreater(deriv,0);self.assertAlmostEqual(deriv,numerical,places=7)
    def test_dynamic_ph_formula(self):
        p=self.m.parameters.acid_base;y=self.y0;ev=self.m.evaluate(1,y);a=ev.diagnostics.observables.cell_acid_base
        H,OH,T,B=a.h_mM,a.oh_mM,y[3]/y[5],self.m.parameters.geometry.cell_buffer_total_fmol/y[5]
        frac=np.array(_carbon_fractions(a.ph,p));k=np.arange(3.);bar=frac@k;b=1/(1+10**(p.cell_buffer_pka-a.ph))
        beta=np.log(10)*(T*(frac@k**2-bar**2)+B*b*(1-b)+H+OH)
        expected=(ev.rhs[4]-bar*ev.rhs[3]-(OH-H)*ev.rhs[5])/(y[5]*beta)
        h=1e-4
        actual=(observe(self.m,1,y+h*ev.rhs)['ph_i']-observe(self.m,1,y-h*ev.rhs)['ph_i'])/(2*h)
        self.assertAlmostEqual(expected,actual,places=7)
    def test_dimensional_trajectory_matches(self):
        v=read('rescaling_verification.json')
        self.assertLess(v['dimensionless_max_scaled_error'],1e-7)
        self.assertLess(abs(v['original_sampled_cumulative_pL']-v['reproduced_sampled_cumulative_pL']),1e-8)
    def test_source_projection(self):
        self.assertLess(read('algebra_verification.json')['max_identity_error_fmol_s'],1e-11)
        expected={'C0':0,'C1':-2,'C2':-3,'C3a':-1,'C3b':1,'C4a':-1,'C4b':1}
        for r in csv.DictReader((OUT/'source_projections.csv').open()):self.assertEqual(float(r['kappa']),expected[r['mechanism']])
    def test_roots_and_jacobians(self):
        r=read('equilibria.json');self.assertEqual(len(r['roots']),26);self.assertFalse(r['failures'])
        for x in r['roots']:
            self.assertLess(x['residual_max'],1e-7);self.assertLess(x['jacobian_relative_step_error'],1e-4)
            self.assertLess(x['max_real_eigenvalue_per_s'],0)
    def test_implicit_sensitivities(self):
        for r in csv.DictReader((OUT/'expression_sensitivities.csv').open()):
            self.assertLess(float(r['validation_relative_error']),1e-4)
            self.assertAlmostEqual(float(r['chloride_compensation_gain']),-2*float(r['dN_de'])/float(r['dJ4_de']),places=10)
        for r in csv.DictReader((OUT/'parameter_sensitivity.csv').open()):self.assertLess(float(r['FD_relative_error']),1e-4)
    def test_compensation_not_relative_increase(self):
        r=read('compensation_comparison.json');self.assertAlmostEqual(r['signed_integrated_compensation'],.885070351,places=8)
        self.assertAlmostEqual(r['relative_NKCC_increase'],.231634495,places=8)
    def test_inverse_threshold_and_gates(self):
        r=read('inverse_threshold.json');self.assertLess(abs(r['tight_check']['adaptive_deficit']-.303),1e-6)
        self.assertTrue(r['tight_check']['physiological']);self.assertLess(r['tight_check']['max_conservation_ratio'],1.)
        self.assertLess(r['nearby'][0]['adaptive_deficit'],.303);self.assertGreater(r['nearby'][1]['adaptive_deficit'],.303)
        l=read('inverse_legacy_threshold.json');self.assertLess(abs(l['check']['sampled_deficit']-.303),1e-6)
    def test_late_ph_failure_is_retained(self):
        r=read('inverse_long_time_summary.json');self.assertGreater(r['ph_gate_crossing_s'],600);self.assertLess(r['ph_gate_crossing_s'],3600)
    def test_source_code_unchanged(self):
        audit=read('source_audit.json')
        for r in audit['files']:
            if r['path'].endswith('.py') or r['path'].startswith('src/'):self.assertTrue(r['matches_task40'])
if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(AnalysisTests)
    result=unittest.TextTestRunner(verbosity=1).run(suite)
    write_json(OUT/'test_summary.json',dict(tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),success=result.wasSuccessful()))
    sys.exit(0 if result.wasSuccessful() else 1)
