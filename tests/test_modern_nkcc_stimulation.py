import unittest

import numpy as np

from modern_full_model.calibration import (
    ModelVariant,
    WTCalibrationParameters,
    build_wt_model,
)
from modern_full_model.camp_pka import R1EffectiveActivation
from modern_full_model.model import ModernFullModel, WT
from modern_full_model.nkcc_stimulation import (
    N0BaselineNkcc1,
    N1AlgebraicNkcc1,
    N2EffectiveDynamicNkcc1,
    attach_stimulated_nkcc1,
    isolated_nkcc1_fractional_cl_uptake_s,
    normalized_calcium_arm,
)
from modern_full_model.validation import (
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    simulate_wt,
    solver_comparison,
)


FROZEN_CORE_ROOT = np.asarray(
    (
        19.18380554451268,
        144.43487028259776,
        65.12999999999998,
        8.05235934438533,
        20.426913621165617,
        1.3,
        14.853386850556735,
        0.4154437452481152,
        14.039093436976788,
        1.3103189911276751,
        1.2297371588280614,
        0.10164106693766255,
    ),
    dtype=float,
)


def frozen_template() -> ModernFullModel:
    variant = ModelVariant(
        "G2_BALANCED_APICAL_K_BIASED_G10",
        mixed_bath_na_attempt_fraction=0.10,
        apical_pump_fraction=0.10,
        apical_k_fraction=0.10,
    )
    calibration = WTCalibrationParameters(
        ae4_carrier_amount_fmol_at_unit_rate_gauge=0.1501264265070657,
        nhe1_capacity_fmol_s=0.007440980590414482,
        common_membrane_conductance_scale=10.0,
        cell_other_impermeant_osmoles_fmol=124.71681310732373,
    )
    resting = build_wt_model(
        variant,
        calibration,
        regulatory_model=R1EffectiveActivation(),
    )
    protocol = SecretagogueProtocol()
    return ModernFullModel(
        parameters=resting.parameters,
        stimulus=protocol,
        regulatory_model=resting.regulatory_model,
        ae4_parameters=resting.ae4_parameters,
    )


class NkccStimulationTests(unittest.TestCase):
    def test_normalized_calcium_arm_has_exact_protocol_endpoints(self) -> None:
        self.assertEqual(normalized_calcium_arm(0.058), 0.0)
        self.assertEqual(normalized_calcium_arm(0.55), 1.0)
        self.assertEqual(normalized_calcium_arm(0.0), 0.0)
        self.assertEqual(normalized_calcium_arm(5.0), 1.0)

    def test_n0_is_exact_baseline_at_rest_and_stimulation(self) -> None:
        template = frozen_template()
        wrapped = attach_stimulated_nkcc1(template, N0BaselineNkcc1())
        state = np.r_[FROZEN_CORE_ROOT, 0.0]
        for time_s in (0.0, 1.0e-6, 120.0):
            expected = template.evaluate(time_s, state, genotype=WT)
            actual = wrapped.evaluate(time_s, state, genotype=WT)
            np.testing.assert_array_equal(actual.rhs, expected.rhs)
            self.assertEqual(
                actual.diagnostics.homeostasis.nkcc1_inward_fmol_s,
                expected.diagnostics.homeostasis.nkcc1_inward_fmol_s,
            )

    def test_n1_unit_multiplier_is_exact_nested_limit_for_every_input(self) -> None:
        template = frozen_template()
        wrapped = attach_stimulated_nkcc1(
            template, N1AlgebraicNkcc1(fully_activated_multiplier=1.0)
        )
        state = np.r_[FROZEN_CORE_ROOT, 0.0]
        for time_s in (0.0, 1.0e-6, 300.0):
            np.testing.assert_array_equal(
                wrapped.evaluate(time_s, state).rhs,
                template.evaluate(time_s, state).rhs,
            )

    def test_n1_scales_only_nkcc_capacity_and_preserves_stoichiometry(self) -> None:
        template = frozen_template()
        baseline = attach_stimulated_nkcc1(template, N0BaselineNkcc1())
        regulated = attach_stimulated_nkcc1(
            template, N1AlgebraicNkcc1(fully_activated_multiplier=1.75)
        )
        state = np.r_[FROZEN_CORE_ROOT, 0.0]
        old = baseline.evaluate(1.0e-6, state)
        new = regulated.evaluate(1.0e-6, state)
        self.assertAlmostEqual(
            new.diagnostics.homeostasis.nkcc1_inward_fmol_s
            / old.diagnostics.homeostasis.nkcc1_inward_fmol_s,
            1.75,
            places=14,
        )
        delta = new.rhs - old.rhs
        self.assertAlmostEqual(delta[0], delta[1], places=14)
        self.assertAlmostEqual(delta[2], 2.0 * delta[0], places=14)
        self.assertAlmostEqual(delta[3], 0.0, places=14)
        self.assertAlmostEqual(delta[4], 0.0, places=14)
        self.assertAlmostEqual(delta[0] + delta[1] - delta[2] - delta[4], 0.0, places=14)

    def test_n2_adds_one_state_and_retains_frozen_basal_root(self) -> None:
        template = frozen_template()
        regulated = attach_stimulated_nkcc1(
            template,
            N2EffectiveDynamicNkcc1(
                fully_activated_multiplier=6.0,
                tau_activation_s=60.0,
            ),
        )
        self.assertEqual(regulated.state_names[-1], "nkcc1_effective_activation_fraction")
        self.assertEqual(len(regulated.state_names), len(template.state_names) + 1)
        state = np.r_[FROZEN_CORE_ROOT, 0.0, 0.0]
        evaluation = regulated.evaluate(0.0, state)
        self.assertLess(float(np.max(np.abs(evaluation.rhs))), 2.0e-16)
        right = regulated.evaluate(1.0e-6, state)
        self.assertAlmostEqual(right.rhs[-1], 1.0 / 60.0, places=14)
        self.assertEqual(
            right.diagnostics.regulatory["nkcc1_capacity_multiplier"], 1.0
        )

    def test_n2_conservation_and_dual_solver_agreement(self) -> None:
        template = frozen_template()
        regulated = attach_stimulated_nkcc1(
            template,
            N2EffectiveDynamicNkcc1(
                fully_activated_multiplier=6.0,
                tau_activation_s=60.0,
            ),
        )
        grid = np.r_[0.0, 1.0e-6, np.arange(5.0, 61.0, 5.0)]
        radau = simulate_wt(
            regulated,
            FROZEN_CORE_ROOT,
            family="R1+N2",
            initial_condition="baseline",
            solver=PRODUCTION_RADAU,
            time_s=grid,
        )
        bdf = simulate_wt(
            regulated,
            FROZEN_CORE_ROOT,
            family="R1+N2",
            initial_condition="baseline",
            solver=PRODUCTION_BDF,
            time_s=grid,
        )
        self.assertTrue(radau.success)
        self.assertTrue(bdf.success)
        self.assertTrue(radau.positive_core)
        self.assertLess(
            radau.max_abs_conservation_residuals["cell_bulk_charge_rate_fmol_s"],
            1.0e-12,
        )
        comparison = solver_comparison(radau, bdf)
        self.assertTrue(comparison["both_success"])
        self.assertLess(comparison["max_relative_state_difference"], 2.0e-5)

    def test_isolated_fractional_rate_is_explicit_mechanistic_diagnostic(self) -> None:
        template = frozen_template()
        baseline = attach_stimulated_nkcc1(template, N0BaselineNkcc1())
        regulated = attach_stimulated_nkcc1(
            template, N1AlgebraicNkcc1(fully_activated_multiplier=1.75)
        )
        state = np.r_[FROZEN_CORE_ROOT, 0.0]
        basal_rate = isolated_nkcc1_fractional_cl_uptake_s(
            baseline, 1.0e-6, state
        )
        regulated_rate = isolated_nkcc1_fractional_cl_uptake_s(
            regulated, 1.0e-6, state
        )
        self.assertGreater(basal_rate, 0.0)
        self.assertAlmostEqual(regulated_rate / basal_rate, 1.75, places=14)

    def test_invalid_new_parameters_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            N1AlgebraicNkcc1(fully_activated_multiplier=0.9)
        with self.assertRaises(ValueError):
            N2EffectiveDynamicNkcc1(tau_activation_s=0.0)
        with self.assertRaises(ValueError):
            normalized_calcium_arm(0.1, resting_calcium_uM=0.5, stimulated_calcium_uM=0.4)


if __name__ == "__main__":
    unittest.main()
