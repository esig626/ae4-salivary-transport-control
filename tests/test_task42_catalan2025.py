"""Source, carbonate interface and inheritance checks; no trajectory solves."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import importlib.util
import sys
import unittest
from dataclasses import asdict, replace
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
spec = importlib.util.spec_from_file_location('task42_test_common',
    ROOT / 'analysis/42_catalan_2025_ae4_mechanism_classes/validation_common.py')
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
import numpy as np
from modern_full_model.ae4_catalan2025 import (
    CLASSES, ClassContribution, forward_affinity, with_mechanism_class)
from modern_full_model.model import WT

LITERAL = {'C0': (-.5, -.5, 1, -2, -2), 'C1': (1, -1, 1, -1, -1),
           'C2': (1, -2, 1, -1, -2), 'C3a': (0, -1, 1, -2, -2),
           'C3b': (-1, 0, 1, -2, -2), 'C4a': (0, -1, 1, -1, -2),
           'C4b': (-1, 0, 1, -1, -2)}
FIELDS = ('na_cell_fmol_s', 'k_cell_fmol_s', 'cl_cell_fmol_s',
          'tic_cell_fmol_s', 'alkalinity_cell_fmol_s')


class SourceClassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rest, cls.stim, cls.y, cls.record = common.parent_models_and_state()

    def test_predeclared_panel_exact_charge_and_reversal(self):
        self.assertEqual(tuple(CLASSES), tuple(LITERAL))
        for name, c in CLASSES.items():
            self.assertEqual(c.conserved_coefficients, LITERAL[name])
            self.assertEqual(c.exact_charge, Fraction(0))
            for j in (.004933930166208064, -.004933930166208064, 0.):
                a, b = ClassContribution(c, j, {}), ClassContribution(c, -j, {})
                self.assertEqual(tuple(getattr(a, f) for f in FIELDS), tuple(v*j for v in LITERAL[name]))
                self.assertEqual(tuple(getattr(b, f) for f in FIELDS), tuple(-getattr(a, f) for f in FIELDS))
                self.assertEqual(a.charge_source_residual_fmol_s, 0.)

    def test_c0_full_rhs_is_bitwise_identical_to_task40(self):
        for old in (self.rest, self.stim):
            new = with_mechanism_class(old, 'C0')
            for t, y in ((0., self.y), (1e-6, np.r_[self.y[:12], .4])):
                for expression in (1., .05, 0.):
                    g = replace(WT, ae4_expression=expression)
                    np.testing.assert_array_equal(new.rhs(t, y, genotype=g), old.rhs(t, y, genotype=g))

    def test_every_class_changes_only_declared_cellular_sources(self):
        for old in (self.rest, self.stim):
            for name, c in CLASSES.items():
                new = with_mechanism_class(old, name)
                self.assertEqual(common.parameter_hashes(new), common.parameter_hashes(old))
                for t, y in ((0., self.y), (1e-6, np.r_[self.y[:12], .4])):
                    for expression in (1., .05, 0.):
                        g = replace(WT, ae4_expression=expression)
                        a, b = old.evaluate(t, y, genotype=g), new.evaluate(t, y, genotype=g)
                        j = a.diagnostics.ae4.cl_cell_fmol_s
                        self.assertEqual(b.diagnostics.ae4.cl_cell_fmol_s, j)
                        for field, nu in zip(FIELDS, LITERAL[name]):
                            self.assertEqual(getattr(b.diagnostics.ae4, field), nu*j)
                        expected = np.zeros(len(y));expected[:5] = (np.array(LITERAL[name])-LITERAL['C0'])*j
                        np.testing.assert_allclose(np.array(b.rhs)-a.rhs, expected, rtol=0, atol=2e-16)
                        for key in ('homeostasis', 'membranes', 'water', 'co2_fluxes_fmol_s', 'observables', 'regulatory'):
                            self.assertEqual(getattr(a.diagnostics, key), getattr(b.diagnostics, key), key)
                        self.assertLess(common.parent.inherited.bookkeeping(new, b)['cell_rhs_assembly_max_abs_fmol_s'], 1e-14)

    def test_carbonate_reaches_full_model_tic_and_ta_separately(self):
        for name in ('C2', 'C4a', 'C4b'):
            model = with_mechanism_class(self.rest, name)
            e = model.evaluate(0., self.y)
            a = e.diagnostics.ae4;j = a.cl_cell_fmol_s
            self.assertEqual(a.hco3_cell_fmol_s, 0.)
            self.assertEqual(a.co3_cell_fmol_s, -j)
            self.assertEqual(a.tic_cell_fmol_s, -j)
            self.assertEqual(a.alkalinity_cell_fmol_s, -2*j)
            self.assertNotEqual(a.tic_cell_fmol_s, a.alkalinity_cell_fmol_s)
            self.assertLess(abs(e.diagnostics.conservation_residuals['carbon_accounting_fmol_s']), 1e-14)

    def test_affinity_electrical_cancellation_and_reversal(self):
        inside = dict(na=11., k=120., cl=60., hco3=5., co3=.002)
        outside = dict(na=140., k=5., cl=120., hco3=24., co3=.03)
        for c in CLASSES.values():
            a = forward_affinity(c, inside, outside, voltage_V=-.07, thermal_voltage_V=.026)
            b = forward_affinity(c, outside, inside, voltage_V=.07, thermal_voltage_V=.026)
            zero = forward_affinity(c, inside, outside, voltage_V=0., thermal_voltage_V=.026)
            self.assertAlmostEqual(a['electrical_delta_g_over_rt'], 0., delta=2e-15)
            self.assertAlmostEqual(a['affinity_over_rt'], zero['affinity_over_rt'], delta=2e-15)
            self.assertAlmostEqual(a['affinity_over_rt'], -b['affinity_over_rt'], delta=2e-15)
        with self.assertRaises(ValueError):
            forward_affinity(CLASSES['C2'], {**inside, 'co3': 0.}, outside,
                             voltage_V=0., thermal_voltage_V=.026)

    def test_carbonate_speciation_uses_the_inherited_pka(self):
        obs = self.rest.evaluate(0., self.y).diagnostics.observables
        pka = self.rest.parameters.acid_base.carbon_pka2
        for state in (obs.cell_acid_base, obs.bath_acid_base):
            self.assertAlmostEqual(state.co3_mM / state.hco3_mM,
                                   10**(state.ph-pka), delta=1e-16)
            self.assertAlmostEqual(state.total_carbon_mM,
                                   state.co2_mM+state.hco3_mM+state.co3_mM, delta=1e-13)

    def test_reference_file_hashes(self):
        self.assertGreaterEqual(len(common.verify_parent_sources()['files']), 180)


if __name__ == '__main__':
    unittest.main(verbosity=2)
