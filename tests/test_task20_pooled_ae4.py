"""Phenotype-independent assay and production-interface tests for Task 20."""
from dataclasses import asdict, replace
from fractions import Fraction
import math
import unittest

import numpy as np

from src.modern_full_model.pooled_ae4 import (
    LEGACY_MODEL_ID, MODEL_ID, PooledAE4Environment, PooledAE4Parameters,
    evaluate_pooled_ae4, select_ae4_mechanism,
)
from src.modern_full_model.run_calcium_fast_screen import load_freeze
from src.modern_full_model.task18_fixed_wt import baseline


def assay_environment():
    # The 2016 thermodynamic illustration, not an inherited WT root.
    return PooledAE4Environment(15.5, 139.5, 50.1, 19., 151.6, 3.4, 124.6, 24.7)


def reverse_environment(e):
    return PooledAE4Environment(**{
        f"{ion}_{s}_mM": getattr(e, f"{ion}_{'o' if s == 'i' else 'i'}_mM")
        for ion in ("na", "k", "cl", "hco3") for s in ("i", "o")
    })


class PooledAssayTests(unittest.TestCase):
    def setUp(self):
        self.e = assay_environment()
        self.p = PooledAE4Parameters(0.2, 3.)

    def test_01_global_equilibrium_is_zero(self):
        e = PooledAE4Environment(140., 5., 120., 25., 140., 5., 120., 25.)
        f = evaluate_pooled_ae4(e, self.p)
        self.assertEqual(f.affinity, 0.)
        self.assertEqual(f.j_ae4_fmol_s, 0.)
        self.assertTrue(all(x == 0 for x in f.intracellular_conserved_sources_fmol_s.values()))

    def test_02_complete_gradient_reversal_reverses_cycle(self):
        f = evaluate_pooled_ae4(self.e, self.p)
        r = evaluate_pooled_ae4(reverse_environment(self.e), self.p)
        self.assertAlmostEqual(f.affinity, -r.affinity, places=14)
        self.assertAlmostEqual(f.j_ae4_fmol_s, -r.j_ae4_fmol_s, places=14)
        self.assertAlmostEqual(f.k_forward_s, r.k_reverse_s, places=14)
        self.assertEqual(r.donor_side, "extracellular")
        self.assertAlmostEqual(r.donor_na_fraction, self.e.na_i_mM / 155.)

    def test_03_positive_cycle_has_declared_sources(self):
        f = evaluate_pooled_ae4(self.e, self.p)
        self.assertGreater(f.j_ae4_fmol_s, 0)
        self.assertGreater(f.source_cl_i_fmol_s, 0)
        for source in (f.source_na_i_fmol_s, f.source_k_i_fmol_s, f.source_hco3_i_fmol_s):
            self.assertLess(source, 0)
        self.assertEqual(f.source_hco3_i_fmol_s, -2*f.j_ae4_fmol_s)
        self.assertEqual(f.source_tic_i_fmol_s, f.source_alkalinity_i_fmol_s)
        self.assertEqual(f.source_tic_i_fmol_s, f.source_hco3_i_fmol_s)
        self.assertAlmostEqual(f.donor_na_fraction, self.e.na_i_mM / 155.)

    def test_04_no_opposed_cation_sources_across_broad_states(self):
        rng = np.random.default_rng(200004)
        for values in np.exp(rng.uniform(-9., 7., (1000, 8))):
            f = evaluate_pooled_ae4(PooledAE4Environment(*values), self.p)
            self.assertGreaterEqual(f.source_na_i_fmol_s * f.source_k_i_fmol_s, 0)

    def test_05_absolute_cation_transport_equals_productive_cycle(self):
        rng = np.random.default_rng(200005)
        for values in np.exp(rng.uniform(-9., 7., (1000, 8))):
            f = evaluate_pooled_ae4(PooledAE4Environment(*values), self.p)
            total = abs(f.source_na_i_fmol_s) + abs(f.source_k_i_fmol_s)
            self.assertTrue(math.isclose(total, abs(f.j_ae4_fmol_s), rel_tol=3e-16))

    def test_06_exact_stoichiometric_charge_and_float_residual(self):
        # Exact rational identity independent of the implementation.
        for p in (Fraction(0), Fraction(1, 7), Fraction(37, 53), Fraction(1)):
            self.assertEqual(-p - (1-p) - 1 + 2, 0)
        for e in (self.e, reverse_environment(self.e)):
            f = evaluate_pooled_ae4(e, self.p)
            residual = math.fsum((f.source_na_i_fmol_s, f.source_k_i_fmol_s,
                                  -f.source_cl_i_fmol_s, -f.source_hco3_i_fmol_s))
            self.assertLessEqual(abs(residual), 5e-16 * abs(f.j_ae4_fmol_s))
            self.assertEqual(residual, f.transported_charge_equivalents_fmol_s)

    def test_07_detailed_balance_and_nonnegative_entropy(self):
        rng = np.random.default_rng(200007)
        for values in np.exp(rng.uniform(-9., 7., (1000, 8))):
            e = PooledAE4Environment(*values)
            f = evaluate_pooled_ae4(e, self.p)
            # Independent mass-action ratio, with no log-affinity call.
            ratio = (e.cl_o_mM * (e.na_i_mM+e.k_i_mM) * e.hco3_i_mM**2 /
                     (e.cl_i_mM * (e.na_o_mM+e.k_o_mM) * e.hco3_o_mM**2))
            self.assertAlmostEqual(math.log(ratio), f.affinity, places=12)
            self.assertAlmostEqual(math.log(f.k_forward_s/f.k_reverse_s), f.affinity, places=12)
            self.assertGreaterEqual(f.entropy_production_over_r_fmol_s, 0)
            self.assertGreaterEqual(f.j_ae4_fmol_s * f.affinity, 0)

    def test_08_pure_and_near_pure_na_support_transport(self):
        for minority in (0., 1e-12):
            e = replace(self.e, na_i_mM=155.-minority, k_i_mM=minority,
                        na_o_mM=155.-minority, k_o_mM=minority)
            for direction in (e, reverse_environment(e)):
                f = evaluate_pooled_ae4(direction, self.p)
                self.assertNotEqual(f.j_ae4_fmol_s, 0)
                self.assertTrue(math.isclose(f.source_na_i_fmol_s, -f.j_ae4_fmol_s, rel_tol=1e-13))

    def test_09_pure_and_near_pure_k_support_transport(self):
        for minority in (0., 1e-12):
            e = replace(self.e, na_i_mM=minority, k_i_mM=155.-minority,
                        na_o_mM=minority, k_o_mM=155.-minority)
            for direction in (e, reverse_environment(e)):
                f = evaluate_pooled_ae4(direction, self.p)
                self.assertNotEqual(f.j_ae4_fmol_s, 0)
                self.assertTrue(math.isclose(f.source_k_i_fmol_s, -f.j_ae4_fmol_s, rel_tol=1e-13))

    def test_10_cation_removal_renormalizes_donor_partition(self):
        for e in (replace(self.e, na_i_mM=0), replace(reverse_environment(self.e), na_o_mM=0)):
            f = evaluate_pooled_ae4(e, self.p)
            self.assertNotEqual(f.j_ae4_fmol_s, 0)
            self.assertEqual(f.donor_na_fraction, 0)
            self.assertEqual(f.source_na_i_fmol_s, 0)
            self.assertEqual(f.source_k_i_fmol_s, -f.j_ae4_fmol_s)

    def test_11_zero_cl_flux_cannot_hide_na_k_exchange(self):
        # Opposed individual cation gradients but equal sums and anion pools.
        e = PooledAE4Environment(5., 140., 120., 25., 140., 5., 120., 25.)
        f = evaluate_pooled_ae4(e, self.p)
        self.assertEqual(f.j_ae4_fmol_s, 0)
        self.assertEqual(f.source_na_i_fmol_s, 0)
        self.assertEqual(f.source_k_i_fmol_s, 0)

    def test_12_legacy_selection_retains_exact_production_outputs(self):
        manifest, _ = load_freeze()
        for root in manifest['roots']:
            model, state, original = baseline(manifest, root)
            pooled = select_ae4_mechanism(model, MODEL_ID)
            restored = select_ae4_mechanism(pooled, LEGACY_MODEL_ID,
                                            parameters=model.ae4_parameters)
            result = restored.evaluate(0., state)
            np.testing.assert_array_equal(result.rhs, original.rhs)
            self.assertEqual(asdict(result.diagnostics.ae4), asdict(original.diagnostics.ae4))
            self.assertEqual(asdict(model.parameters), asdict(pooled.parameters))
            self.assertEqual(model.stimulus, pooled.stimulus)
            self.assertEqual(model.regulatory_model, pooled.regulatory_model)

    def test_regulation_scales_only_active_capacity(self):
        f = evaluate_pooled_ae4(self.e, self.p)
        for gain in (0., 0.2, 1.25, 2.):
            g = evaluate_pooled_ae4(self.e, self.p, regulation_gain=gain)
            self.assertEqual(g.affinity, f.affinity)
            self.assertEqual(g.k_forward_s, f.k_forward_s)
            self.assertEqual(g.k_reverse_s, f.k_reverse_s)
            self.assertAlmostEqual(g.j_ae4_fmol_s, gain*f.j_ae4_fmol_s)

    def test_invalid_activities_and_capacities_are_rejected(self):
        for kwargs in ({'na_i_mM': -1}, {'cl_i_mM': 0}, {'na_i_mM': 0, 'k_i_mM': 0},
                       {'hco3_o_mM': math.nan}, {'k_o_mM': math.inf}):
            with self.assertRaises(ValueError):
                replace(self.e, **kwargs)
        for kwargs in ({'carrier_amount_fmol': -1}, {'attempt_rate_s': 0}, {'attempt_rate_s': math.inf}):
            with self.assertRaises(ValueError):
                replace(self.p, **kwargs)
        for gain in (-1., math.nan):
            with self.assertRaises(ValueError):
                evaluate_pooled_ae4(self.e, self.p, regulation_gain=gain)


if __name__ == '__main__':
    unittest.main()
