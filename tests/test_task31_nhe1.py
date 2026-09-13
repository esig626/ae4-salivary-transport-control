"""Task 31 algebra fixtures and separately identified source reconstruction.

No stationary solver, density calibration or integration is called.
"""
from dataclasses import replace
from fractions import Fraction
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from modern_full_model.nhe1_cha2009 import (
    NHE1_CHA_2009, ChaNHE1Kinetics,
    cha_nhe1_rate, cha_nhe1_flux_fmol_s, published_cha_kinetics,
)
from modern_full_model.membranes import (
    HomeostasisEnvironment, NHE1_VERA_SIGUENZA_2018, evaluate_homeostasis,
    evaluate_nhe1_vera_siguenza_2018,
)
from modern_full_model.parameters import FullModelParameters


def synthetic_kinetics():
    # Rational arithmetic fixture, NOT values attributed to Cha et al.
    return ChaNHE1Kinetics(
        k1_plus_per_ms=2.0, k1_minus_per_ms=3.0, k2_plus_per_ms=5.0,
        k_na_o_mM=10.0, k_na_i_mM=20.0, k_h_o_mM=1e-4,
        k_h_i_mM=2e-4, modifier_k_i_mM=1e-4, modifier_k_o_mM=4e-5,
    )


def environment(na_i=20.0, h_i=2e-4, na_o=30.0, h_o=1e-4):
    return HomeostasisEnvironment(
        na_i_mM=na_i, k_i_mM=120.0, cl_i_mM=50.0, h_i_mM=h_i,
        hco3_i_mM=10.0, na_e_mM=na_o, k_e_mM=5.0, cl_e_mM=126.0,
        h_e_mM=h_o, hco3_e_mM=22.0,
    )


class ChaEquationTests(unittest.TestCase):
    def rate(self, na_i=20.0, h_i=2e-4, na_o=30.0, h_o=1e-4):
        return cha_nhe1_rate(
            na_i_mM=na_i, h_i_mM=h_i, na_o_mM=na_o, h_o_mM=h_o,
            kinetics=synthetic_kinetics(),
        )

    def test_independent_rational_eight_state_equation(self):
        # P_o(Na only)=3/8; P_i(H only)=P_i(Na only)=1/4;
        # P_o(H only)=1/8; Eq. 6 gives k2-=10/3.
        a, b, c, d = Fraction(3, 4), Fraction(5, 4), Fraction(3, 4), Fraction(5, 12)
        exchange_ms = (a*b-c*d)/(a+b+c+d)
        modifier = Fraction(16, 23)
        result = self.rate()
        self.assertAlmostEqual(result.exchange_per_s, float(exchange_ms * 1000), places=11)
        self.assertAlmostEqual(result.modifier, float(modifier), places=14)
        self.assertAlmostEqual(result.regulated_per_s, float(exchange_ms*modifier*1000), places=11)

    def test_independent_orientation_balance(self):
        # A two-orientation balance reconstructs the same eight microstate QSS.
        # Simultaneous Na/H occupancy cannot flip and remains in the partition.
        a, b, c, d = 0.75, 1.25, 0.75, 5.0/12.0
        outside = (b+c)/(a+b+c+d)
        inside = 1.0-outside
        na_flux = a*outside-c*inside
        h_flux = b*inside-d*outside
        self.assertAlmostEqual(na_flux, h_flux, places=15)
        self.assertAlmostEqual(self.rate().exchange_per_s, 1000.0*na_flux, places=11)

    def test_microscopic_reversibility_and_zero_at_equilibrium(self):
        k = synthetic_kinetics()
        self.assertAlmostEqual(k.k2_minus_per_ms, 10.0/3.0, places=13)
        self.assertAlmostEqual(
            k.k_h_i_mM*k.k_na_o_mM/(k.k_h_o_mM*k.k_na_i_mM),
            k.k1_plus_per_ms*k.k2_plus_per_ms/(k.k1_minus_per_ms*k.k2_minus_per_ms),
            places=14,
        )
        self.assertEqual(self.rate(na_i=20, na_o=20, h_i=1e-4, h_o=1e-4).regulated_per_s, 0.0)
        self.assertAlmostEqual(self.rate(na_i=20, na_o=40, h_i=1e-4, h_o=2e-4).regulated_per_s, 0.0, places=10)

    def test_forward_and_reverse_use_same_cycle(self):
        self.assertGreater(self.rate().regulated_per_s, 0.0)
        self.assertLess(self.rate(na_i=300).regulated_per_s, 0.0)

    def test_intracellular_acid_activation_and_alkaline_suppression(self):
        # Source-supported pH range, but synthetic constants: algebra check.
        acidic = self.rate(na_i=10, na_o=140, h_i=10**(3-6.4), h_o=10**(3-7.4))
        neutral = self.rate(na_i=10, na_o=140, h_i=10**(3-6.9), h_o=10**(3-7.4))
        alkaline = self.rate(na_i=10, na_o=140, h_i=10**(3-7.5), h_o=10**(3-7.4))
        self.assertGreater(acidic.modifier, neutral.modifier)
        self.assertGreater(neutral.modifier, alkaline.modifier)
        self.assertGreater(acidic.regulated_per_s, neutral.regulated_per_s)
        self.assertGreater(neutral.regulated_per_s, alkaline.regulated_per_s)
        self.assertLess(alkaline.modifier, 0.05)

    def test_extracellular_protons_suppress_mod2(self):
        self.assertGreater(self.rate(h_o=1e-5).modifier, self.rate(h_o=1e-4).modifier)

    def test_finite_saturation_and_extreme_positive_concentrations(self):
        # Fixed numerical edge cases, not a biological parameter search.
        cases = (
            (1e-300, 1e-4, 1e300, 1e-300),
            (1e300, 1e-300, 1e-300, 1e300),
            (1e300, 1e300, 1e300, 1e300),
            (1e-300, 1e-300, 1e-300, 1e-300),
        )
        for ni, hi, no, ho in cases:
            r = self.rate(na_i=ni, h_i=hi, na_o=no, h_o=ho)
            self.assertTrue(all(math.isfinite(v) for v in r.__dict__.values()))
            self.assertLessEqual(abs(r.regulated_per_s), synthetic_kinetics().turnover_bound_per_s)
            self.assertTrue(0.0 <= r.modifier <= 1.0)
        saturated = self.rate(na_i=1e-10, na_o=1e8)
        more = self.rate(na_i=1e-10, na_o=1e12)
        self.assertAlmostEqual(saturated.regulated_per_s/more.regulated_per_s, 1.0, places=6)

    def test_amount_conversion_and_shared_activity_scale(self):
        args = dict(na_i_mM=20, h_i_mM=2e-4, na_o_mM=30, h_o_mM=1e-4,
                    kinetics=synthetic_kinetics(), carrier_amount_fmol=0.002)
        self.assertAlmostEqual(cha_nhe1_flux_fmol_s(**args), 0.002*self.rate().regulated_per_s, places=15)
        self.assertEqual(cha_nhe1_flux_fmol_s(**args, activity_scale=0), 0.0)
        self.assertEqual(cha_nhe1_flux_fmol_s(**args, activity_scale=0.5), 0.5*cha_nhe1_flux_fmol_s(**args))

    def test_na_h_stoichiometry_zero_carbon_and_charge_both_directions(self):
        base = FullModelParameters()
        p = replace(base, homeostasis=replace(
            base.homeostasis, nhe1_model=NHE1_CHA_2009,
            nhe1_cha_carrier_amount_fmol=0.002,
            nkcc1_capacity_fmol_s=0.0, ae2_capacity_fmol_s=0.0,
        ))
        # Test-only injection: no synthetic constants are installed in production.
        with patch('modern_full_model.membranes.published_cha_kinetics', return_value=synthetic_kinetics()):
            for e in (environment(), environment(na_i=300)):
                r = evaluate_homeostasis(e, p)
                self.assertEqual(r.na_cell_fmol_s, r.nhe1_inward_fmol_s)
                self.assertEqual(r.alkalinity_cell_fmol_s, r.nhe1_inward_fmol_s)
                self.assertEqual(r.tic_cell_fmol_s, 0.0)
                self.assertEqual(r.k_cell_fmol_s, 0.0)
                self.assertEqual(r.cl_cell_fmol_s, 0.0)
                self.assertEqual(r.charge_source_residual_fmol_s, 0.0)

    def test_production_uses_exact_table_s1_and_requires_salivary_amount(self):
        k = published_cha_kinetics()
        self.assertEqual(tuple(k.__dict__.values()),
                         (10.5, 0.201, 15.8, 195.0, 16.2, 1.62e-3,
                          6.05e-4, 3.07e-5, 4.8e-7, 183.0))
        self.assertEqual(k.k2_minus_per_ms, 183.0)
        base = FullModelParameters()
        p = replace(base, homeostasis=replace(base.homeostasis, nhe1_model=NHE1_CHA_2009))
        with self.assertRaisesRegex(ValueError, 'not been calibrated'):
            evaluate_homeostasis(environment(), p)

    def test_missing_density_never_defaults_to_cardiac_amount(self):
        base = FullModelParameters()
        p = replace(base, homeostasis=replace(base.homeostasis, nhe1_model=NHE1_CHA_2009))
        self.assertIsNone(p.homeostasis.nhe1_cha_carrier_amount_fmol)
        with patch('modern_full_model.membranes.published_cha_kinetics', return_value=synthetic_kinetics()):
            with self.assertRaisesRegex(ValueError, 'not been calibrated'):
                evaluate_homeostasis(environment(), p)

    def test_invalid_inputs_rejected(self):
        for bad in (0.0, -1.0, math.nan, math.inf):
            with self.assertRaises(ValueError):
                self.rate(h_i=bad)
        with self.assertRaises(ValueError):
            replace(synthetic_kinetics(), k_na_o_mM=-1)
        for bad in (-1.0, math.nan, math.inf):
            base = FullModelParameters()
            with self.assertRaises(ValueError):
                replace(base, homeostasis=replace(base.homeostasis, nhe1_cha_carrier_amount_fmol=bad))

    def test_task30_reference_remains_selectable(self):
        base = FullModelParameters()
        p = replace(base, homeostasis=replace(base.homeostasis, nhe1_model=NHE1_VERA_SIGUENZA_2018))
        self.assertEqual(evaluate_homeostasis(environment(), p).nhe1_inward_fmol_s,
                         evaluate_nhe1_vera_siguenza_2018(environment(), p))

    def test_published_turnover_and_figure_reconstruction(self):
        # Independent rational evaluation of Eq. 3 and Eq. 5, using the
        # printed Table S1 constants, not production binding functions.
        def independent(ni, hi, no, ho):
            ni, hi, no, ho = (Fraction(str(x)) for x in (ni, hi, no, ho))
            ko, ki = Fraction('195'), Fraction('16.2')
            kho, khi = Fraction('0.00162'), Fraction('0.000605')
            a = Fraction('10.5')*no/(ko+no)*kho/(kho+ho)
            b = Fraction('15.8')*hi/(khi+hi)*ki/(ki+ni)
            c = Fraction('0.201')*ni/(ki+ni)*khi/(khi+hi)
            d = Fraction('183')*ho/(kho+ho)*ko/(ko+no)
            mod = 1/(1+(1+ho/Fraction('0.00000048'))*(Fraction('0.0000307')/hi)**3)
            return float(1000*(a*b-c*d)/(a+b+c+d)*mod)

        k = published_cha_kinetics()
        # Main-text turnover constraint and Fig. 4 C/D forward/reverse
        # conditions. 1e-12 mM represents the source's zero-Na limit.
        cases = ((1e-12, 1e-3, 140, 10**-4.4),
                 (1e-12, 10**-3.5, 1, 10**-4.4),
                 (1, 10**-3.5, 1e-12, 10**-4.4))
        for ni, hi, no, ho in cases:
            actual = cha_nhe1_rate(na_i_mM=ni, h_i_mM=hi, na_o_mM=no,
                                   h_o_mM=ho, kinetics=k).regulated_per_s
            self.assertAlmostEqual(actual/independent(ni, hi, no, ho), 1.0, places=12)
        turnover_ms = independent(*cases[0])/1000
        self.assertAlmostEqual(turnover_ms, 2.521212543779516, places=12)
        # This is a documented source inconsistency, not a silently relaxed
        # 2.3--2.5 ms^-1 acceptance gate. No kinetic parameter was adjusted.
        self.assertGreater(turnover_ms, 2.5)
        self.assertGreater(independent(*cases[1]), 0)
        self.assertLess(independent(*cases[2]), 0)

    def test_printed_rounding_cycle_reversal_is_explicit(self):
        k = published_cha_kinetics()
        self.assertAlmostEqual(k.k2_minus_eq6_per_ms, 183.6074807146718, places=10)
        ratio = k.k2_minus_eq6_per_ms/k.k2_minus_per_ms
        args = dict(na_i_mM=20, na_o_mM=20, h_o_mM=1e-4, kinetics=k)
        # With the printed constants, physical chemical equilibrium has a
        # small residual; the printed cycle zero is shifted by 0.00144 pH.
        self.assertGreater(cha_nhe1_rate(**args, h_i_mM=1e-4).regulated_per_s, 0)
        self.assertAlmostEqual(cha_nhe1_rate(**args, h_i_mM=1e-4/ratio).regulated_per_s,
                               0.0, places=10)


if __name__ == '__main__':
    unittest.main()
