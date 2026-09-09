"""Focused invariant and nesting tests for the Task 13B regulator."""

from __future__ import annotations

from dataclasses import fields
import math
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from modern_full_model.camp_pka import (
    AE4Construct,
    INTEGRATION_TRIAL_TOLERANCE,
    P21_GAIN_SENSITIVITIES,
    R0Static,
    R1EffectiveActivation,
    R2CampActivation,
    R3PkaRegulatedFraction,
    R4_REJECTION,
    RegulatoryGain,
    RegulatoryOutput,
    default_regulatory_families,
)
from modern_full_model.model import ConstantStimulus, ModernFullModel, Stimulus


class TestFamilyAPI(unittest.TestCase):
    def test_registry_has_exactly_the_nested_r0_to_r3_families(self) -> None:
        families = default_regulatory_families()
        self.assertEqual(tuple(families), ("R0", "R1", "R2", "R3"))
        self.assertIsInstance(families["R0"], R0Static)
        self.assertIsInstance(families["R1"], R1EffectiveActivation)
        self.assertIsInstance(families["R2"], R2CampActivation)
        self.assertIsInstance(families["R3"], R3PkaRegulatedFraction)
        for name, family in families.items():
            self.assertEqual(family.family, name)
            self.assertTrue(hasattr(family, "state_names"))
            self.assertTrue(callable(family.initial_state))
            self.assertTrue(callable(family.rhs))
            self.assertTrue(callable(family.evaluate))

    def test_state_names_and_zero_state_static_control_are_explicit(self) -> None:
        families = default_regulatory_families()
        self.assertEqual(families["R0"].state_names, ())
        self.assertEqual(families["R0"].initial_state(0.4), ())
        self.assertEqual(families["R0"].rhs(0.0, (), 0.4), ())
        self.assertEqual(
            families["R1"].state_names,
            ("ae4_effective_activation_fraction",),
        )
        self.assertEqual(
            families["R2"].state_names,
            ("camp_effective_fraction", "ae4_activation_fraction"),
        )
        self.assertEqual(
            families["R3"].state_names,
            ("pka_effective_activity_fraction", "ae4_regulated_fraction"),
        )

    def test_output_has_the_stable_model_fields(self) -> None:
        output = R1EffectiveActivation().evaluate(2.0, (0.5,), 1.0)
        self.assertIsInstance(output, RegulatoryOutput)
        names = {item.name for item in fields(output)}
        self.assertTrue({"family", "activity", "capacity_multiplier"} <= names)
        self.assertEqual(output.family, "R1")
        self.assertEqual(output.coupling_hypothesis, "COMMON_AE4_CAPACITY_ONLY")
        self.assertEqual(output.kinetic_evidence_status, "UNMEASURED_ACUTE_KINETICS")

    def test_modern_full_model_accepts_every_family_without_adapter(self) -> None:
        for name, regulatory in default_regulatory_families().items():
            model = ModernFullModel(
                stimulus=ConstantStimulus(beta_input=0.0),
                regulatory_model=regulatory,
            )
            initial = model.initial_state()
            if regulatory.state_names:
                self.assertEqual(
                    model.state_names[-len(regulatory.state_names) :],
                    regulatory.state_names,
                )
            else:
                self.assertNotIn("ae4_effective_activation_fraction", model.state_names)
            evaluation = model.evaluate(0.0, initial)
            self.assertEqual(evaluation.diagnostics.regulatory["family"], name)
            self.assertEqual(
                evaluation.diagnostics.regulatory["capacity_multiplier"], 1.0
            )


class TestEquilibriaAndNesting(unittest.TestCase):
    def test_analytic_initial_states_are_exact_constant_input_equilibria(self) -> None:
        for beta in (0.0, 0.13, 0.5, 1.0):
            for family in default_regulatory_families().values():
                state = family.initial_state(beta)
                derivative = family.rhs(7.0, state, beta)
                np.testing.assert_allclose(derivative, np.zeros(len(state)), atol=1.0e-15)

    def test_r0_is_exact_equilibrium_output_map_of_default_dynamic_families(self) -> None:
        families = default_regulatory_families()
        for beta in np.linspace(0.0, 1.0, 11):
            reference = families["R0"].evaluate(0.0, (), beta)
            for name in ("R1", "R2", "R3"):
                family = families[name]
                output = family.evaluate(0.0, family.equilibrium(beta), beta)
                self.assertAlmostEqual(output.activity, reference.activity, places=14)
                self.assertAlmostEqual(
                    output.capacity_multiplier,
                    reference.capacity_multiplier,
                    places=14,
                )

    def test_r2_collapses_exactly_to_r1_on_constant_input_camp_manifold(self) -> None:
        beta = 0.43
        activation = 0.17
        r1 = R1EffectiveActivation(tau_activation_s=23.0)
        r2 = R2CampActivation(tau_camp_s=9.0, tau_activation_s=23.0)
        r1_rhs = r1.rhs(0.0, (activation,), beta)
        r2_rhs = r2.rhs(0.0, (beta, activation), beta)
        self.assertEqual(r2_rhs[0], 0.0)
        self.assertEqual(r2_rhs[1], r1_rhs[0])
        self.assertEqual(
            r2.evaluate(0.0, (beta, activation), beta).capacity_multiplier,
            r1.evaluate(0.0, (activation,), beta).capacity_multiplier,
        )

    def test_default_r2_and_r3_are_exactly_dynamically_isomorphic(self) -> None:
        r2 = R2CampActivation(
            tau_camp_s=11.0,
            tau_activation_s=29.0,
        )
        r3 = R3PkaRegulatedFraction(
            tau_pka_s=11.0,
            forward_regulation_rate_s=1.0 / 29.0,
            reverse_regulation_rate_s=1.0 / 29.0,
        )
        for beta, state in ((0.0, (0.2, 0.8)), (0.6, (0.1, 0.7)), (1.0, (0.9, 0.1))):
            np.testing.assert_allclose(
                r3.rhs(5.0, state, beta),
                r2.rhs(5.0, state, beta),
                rtol=0.0,
                atol=1.0e-15,
            )
            self.assertEqual(
                r3.evaluate(5.0, state, beta).capacity_multiplier,
                r2.evaluate(5.0, state, beta).capacity_multiplier,
            )

    def test_asymmetric_r3_equilibrium_is_analytic(self) -> None:
        model = R3PkaRegulatedFraction(
            tau_pka_s=13.0,
            forward_regulation_rate_s=0.2,
            reverse_regulation_rate_s=0.05,
        )
        beta = 0.3
        expected = 0.2 * beta / (0.2 * beta + 0.05 * (1.0 - beta))
        state = model.equilibrium(beta)
        self.assertAlmostEqual(state[0], beta)
        self.assertAlmostEqual(state[1], expected)
        np.testing.assert_allclose(model.rhs(0.0, state, beta), (0.0, 0.0), atol=1.0e-15)


class TestPositivityAndBounds(unittest.TestCase):
    def test_r1_vector_field_points_inward_at_fraction_boundaries(self) -> None:
        model = R1EffectiveActivation()
        self.assertGreaterEqual(model.rhs(0.0, (0.0,), 0.7)[0], 0.0)
        self.assertLessEqual(model.rhs(0.0, (1.0,), 0.2)[0], 0.0)

    def test_r2_vector_field_points_inward_on_every_face(self) -> None:
        model = R2CampActivation()
        self.assertGreaterEqual(model.rhs(0.0, (0.0, 0.4), 0.7)[0], 0.0)
        self.assertLessEqual(model.rhs(0.0, (1.0, 0.4), 0.2)[0], 0.0)
        self.assertGreaterEqual(model.rhs(0.0, (0.7, 0.0), 0.4)[1], 0.0)
        self.assertLessEqual(model.rhs(0.0, (0.3, 1.0), 0.4)[1], 0.0)

    def test_r3_vector_field_points_inward_on_every_face(self) -> None:
        model = R3PkaRegulatedFraction(
            forward_regulation_rate_s=0.2,
            reverse_regulation_rate_s=0.7,
        )
        self.assertGreaterEqual(model.rhs(0.0, (0.0, 0.4), 0.7)[0], 0.0)
        self.assertLessEqual(model.rhs(0.0, (1.0, 0.4), 0.2)[0], 0.0)
        self.assertGreaterEqual(model.rhs(0.0, (0.7, 0.0), 0.4)[1], 0.0)
        self.assertLessEqual(model.rhs(0.0, (0.3, 1.0), 0.4)[1], 0.0)

    def test_step_trajectories_remain_bounded_and_monotone(self) -> None:
        for family in default_regulatory_families().values():
            if not family.state_names:
                continue
            solution = solve_ivp(
                lambda t, y: family.rhs(t, y, 1.0),
                (0.0, 300.0),
                np.asarray(family.initial_state(0.0), dtype=float),
                method="Radau",
                rtol=1.0e-9,
                atol=1.0e-12,
                max_step=1.0,
            )
            self.assertTrue(solution.success, solution.message)
            self.assertGreaterEqual(float(np.min(solution.y)), -1.0e-12)
            self.assertLessEqual(float(np.max(solution.y)), 1.0 + 1.0e-12)
            self.assertTrue(np.all(np.diff(solution.y, axis=1) >= -1.0e-10))

    def test_implicit_solver_trial_band_is_separate_from_input_validation(self) -> None:
        model = R1EffectiveActivation()
        small_overshoot = -0.5 * INTEGRATION_TRIAL_TOLERANCE
        # A Radau internal trial state in the declared band is clamped only for
        # RHS/output evaluation; it is not accepted as a beta input.
        self.assertGreater(model.rhs(1.0, (small_overshoot,), 1.0)[0], 0.0)
        self.assertEqual(
            model.evaluate(1.0, (small_overshoot,), 1.0).activity,
            0.0,
        )
        with self.assertRaises(ValueError):
            model.rhs(1.0, (-2.0 * INTEGRATION_TRIAL_TOLERANCE,), 1.0)
        with self.assertRaises(ValueError):
            model.initial_state(1.0 + 0.5 * INTEGRATION_TRIAL_TOLERANCE)

    def test_short_full_model_radau_step_regression(self) -> None:
        class RightLimitBetaStep:
            def __call__(self, time_s: float) -> Stimulus:
                return Stimulus(calcium_uM=0.10, beta_input=0.0 if time_s <= 0.0 else 1.0)

        full = ModernFullModel(
            stimulus=RightLimitBetaStep(),
            regulatory_model=R1EffectiveActivation(tau_activation_s=30.0),
        )
        initial = full.initial_state()
        result = full.solve_dynamics(
            (1.0e-6, 0.05),
            initial_state=initial,
            method="Radau",
            rtol=1.0e-7,
            atol=1.0e-10,
            max_step_s=0.005,
        )
        self.assertTrue(result.success, result.message)
        self.assertTrue(np.all(np.isfinite(result.y)))
        self.assertGreater(result.y[-1, -1], result.y[-1, 0])
        self.assertLessEqual(result.y[-1, -1], 1.0)

    def test_invalid_domains_are_rejected(self) -> None:
        for constructor, keyword in (
            (R1EffectiveActivation, {"tau_activation_s": 0.0}),
            (R2CampActivation, {"tau_camp_s": -1.0}),
            (R3PkaRegulatedFraction, {"forward_regulation_rate_s": 0.0}),
        ):
            with self.assertRaises(ValueError):
                constructor(**keyword)
        with self.assertRaises(ValueError):
            RegulatoryGain(fully_activated_increment=-0.1)
        with self.assertRaises(ValueError):
            RegulatoryGain(
                fully_activated_increment=1.0e308,
                coupling_scale=1.0e308,
            ).capacity_multiplier(1.0)
        model = R2CampActivation()
        with self.assertRaises(ValueError):
            model.rhs(0.0, (0.1,), 0.5)
        with self.assertRaises(ValueError):
            model.rhs(0.0, (0.1, 0.2), 1.1)
        with self.assertRaises(ValueError):
            model.evaluate(math.nan, (0.1, 0.2), 0.5)
        with self.assertRaises(ValueError):
            model.evaluate(0.0, (-0.01, 0.2), 0.5)


class TestPrimaryEvidenceLogic(unittest.TestCase):
    def test_gain_discrepancy_is_preserved_as_sensitivity_not_target(self) -> None:
        self.assertEqual(
            tuple(item.fully_activated_increment for item in P21_GAIN_SENSITIVITIES),
            (0.25, 0.70),
        )
        for item in P21_GAIN_SENSITIVITIES:
            self.assertEqual(item.evidential_status, "SENSITIVITY_ONLY_NOT_POINT_TARGET")
        low = RegulatoryGain(fully_activated_increment=0.25).capacity_multiplier(1.0)
        high = RegulatoryGain(fully_activated_increment=0.70).capacity_multiplier(1.0)
        self.assertEqual(low, 1.25)
        self.assertEqual(high, 1.70)

    def test_s173_dependence_and_s273_retention_are_qualitative_switches(self) -> None:
        wt = R3PkaRegulatedFraction(construct=AE4Construct.WT)
        s173a = R3PkaRegulatedFraction(construct=AE4Construct.S173A)
        s273a = R3PkaRegulatedFraction(construct=AE4Construct.S273A)
        active_state = (1.0, 1.0)
        wt_output = wt.evaluate(0.0, active_state, 1.0)
        s173a_output = s173a.evaluate(0.0, active_state, 1.0)
        s273a_output = s273a.evaluate(0.0, active_state, 1.0)
        self.assertGreater(wt_output.capacity_multiplier, 1.0)
        self.assertEqual(s173a_output.activity, 1.0)
        self.assertEqual(s173a_output.capacity_multiplier, 1.0)
        self.assertEqual(s273a_output.capacity_multiplier, wt_output.capacity_multiplier)

    def test_zero_coupling_is_only_a_qualitative_pka_dependence_control(self) -> None:
        blocked = R1EffectiveActivation(gain=RegulatoryGain(coupling_scale=0.0))
        self.assertEqual(blocked.evaluate(0.0, (1.0,), 1.0).capacity_multiplier, 1.0)
        self.assertEqual(blocked.evaluate(0.0, (0.0,), 0.0).capacity_multiplier, 1.0)

    def test_r4_is_explicitly_rejected_not_silently_implemented(self) -> None:
        self.assertEqual(R4_REJECTION["family"], "R4")
        self.assertEqual(R4_REJECTION["status"], "REJECTED_NOT_RETAINED")
        self.assertIn("unsupported freedom", R4_REJECTION["reason"])


if __name__ == "__main__":
    unittest.main()
