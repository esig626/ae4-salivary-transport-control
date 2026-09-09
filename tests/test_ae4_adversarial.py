"""Independent adversarial gates for the Task 12 AE4 reconstruction."""

from __future__ import annotations

from dataclasses import replace
import inspect
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ae4_mechanism_reconstruction.chassis import (  # noqa: E402
    AE2_KNOCKOUT,
    AE4_KNOCKOUT,
    WT,
    AE4Balance,
    FixedChassis,
    HistoricalSodiumOnlyAE4,
    dimension_ledger,
)
from ae4_mechanism_reconstruction.comparison import (  # noqa: E402
    CAPACITY_INTERVAL_SELECTION_FRACTIONS,
    CAPACITY_PROFILE_LOG10_OFFSETS,
    CommonPkaGateMechanism,
    HISTORICAL_WT_AE4_CONTRIBUTION,
    _assay_admissible,
    calibrate_candidate,
    candidate_chassis,
    fit_nonnegative_activity_scale,
    rest_admissible_capacity_envelopes,
    solve_resting_equilibrium,
    stimulus_activation_comparisons,
)
from ae4_mechanism_reconstruction.mechanisms import (  # noqa: E402
    AE4Context,
    AE4Mechanism,
    AE4Parameters,
    CANDIDATE_SPECS,
    PKA_FOLD_ASSAY_VALUES,
    candidate_specs,
    equilibrium_context,
    evaluate_ae4,
)


CONTEXT = AE4Context(
    na_i=25.0,
    k_i=120.0,
    cl_i=50.0,
    hco3_i=10.0,
    na_e=140.2,
    k_e=5.3,
    cl_e=102.6,
    hco3_e=21.0,
)


class CandidateBalanceAudit(unittest.TestCase):
    def test_complete_source_justified_registry_is_present(self) -> None:
        identifiers = {spec.candidate_id for spec in candidate_specs()}
        required = {
            "C1",
            "C1_PKA_COMMON",
            "C1_PKA_K_RECRUIT",
            "C2",
            "C3",
            "C3_PKA_COMMON",
            "C4",
            "C5",
            "C5_GATE_ONLY_REJECTED",
            "C6_CARBONATE_UNTESTABLE",
            "C6_HCO3_SELF_EXCHANGE_UNTESTABLE",
            "C6_112",
            "C6_112_PKA_COMMON",
            "C6_123",
            "C6_123_PKA_COMMON",
            "C6_213",
            "C6_213_PKA_COMMON",
            "C7_112",
            "C7_123",
            "C7_213",
            "C8_112_MA",
            "C8_112_SAT",
            "C8_123_MA",
            "C8_123_SAT",
            "C8_213_MA",
            "C8_213_SAT",
        }
        self.assertEqual(identifiers, required)
        rejected = CANDIDATE_SPECS["C5_GATE_ONLY_REJECTED"]
        self.assertFalse(rejected.active)
        self.assertIn("direct Na and K", rejected.rejection_reason or "")
        carbonate = CANDIDATE_SPECS["C6_CARBONATE_UNTESTABLE"]
        self.assertFalse(carbonate.active)
        self.assertIsNone(carbonate.stoichiometry)
        self.assertIn("no carbonate state", carbonate.rejection_reason or "")
        self_exchange = CANDIDATE_SPECS["C6_HCO3_SELF_EXCHANGE_UNTESTABLE"]
        self.assertFalse(self_exchange.active)
        self.assertIsNone(self_exchange.stoichiometry)
        self.assertIn("unresolved stoichiometry", self_exchange.rejection_reason or "")

    def test_every_active_candidate_inserts_its_declared_balance_signature(self) -> None:
        for spec in candidate_specs(include_rejected=False):
            self.assertIsNotNone(spec.stoichiometry)
            a, b, c = spec.stoichiometry or (0, 0, 0)
            contribution = evaluate_ae4(spec.candidate_id, CONTEXT)
            scale = max(1.0, abs(contribution.cycle_flux))
            self.assertAlmostEqual(
                contribution.cl_i,
                a * contribution.cycle_flux,
                delta=2e-12 * scale,
                msg=spec.candidate_id,
            )
            self.assertAlmostEqual(
                contribution.na_i + contribution.k_i,
                -b * contribution.cycle_flux,
                delta=2e-12 * scale,
                msg=spec.candidate_id,
            )
            self.assertAlmostEqual(
                contribution.hco3_i,
                -c * contribution.cycle_flux,
                delta=2e-12 * scale,
                msg=spec.candidate_id,
            )
            self.assertAlmostEqual(
                contribution.na_i
                + contribution.k_i
                - contribution.cl_i
                - contribution.hco3_i,
                0.0,
                delta=2e-12 * scale,
                msg=spec.candidate_id,
            )
            self.assertTrue(contribution.valid, spec.candidate_id)

    def test_parallel_candidates_put_na_and_k_flux_in_separate_balances(self) -> None:
        parallel = [
            spec
            for spec in candidate_specs(include_rejected=False)
            if "parallel" in spec.cation_mode
            and spec.candidate_id != "C1_PKA_K_RECRUIT"
        ]
        for spec in parallel:
            contribution = evaluate_ae4(spec.candidate_id, CONTEXT)
            diagnostics = contribution.diagnostics
            self.assertIn("na_flux", diagnostics, spec.candidate_id)
            self.assertIn("k_flux", diagnostics, spec.candidate_id)
            b = (spec.stoichiometry or (0, 0, 0))[1]
            scale = max(1.0, abs(contribution.cycle_flux))
            self.assertAlmostEqual(
                contribution.na_i,
                -b * float(diagnostics["na_flux"]),
                delta=2e-12 * scale,
                msg=spec.candidate_id,
            )
            self.assertAlmostEqual(
                contribution.k_i,
                -b * float(diagnostics["k_flux"]),
                delta=2e-12 * scale,
                msg=spec.candidate_id,
            )

        # This explicitly novel candidate is Na-only at rest and recruits its
        # K branch only under the independent beta/cAMP/PKA input.
        recruited = evaluate_ae4(
            "C1_PKA_K_RECRUIT", CONTEXT, stimulus_activation=1.0
        )
        self.assertAlmostEqual(
            recruited.na_i, -float(recruited.diagnostics["na_flux"])
        )
        self.assertAlmostEqual(
            recruited.k_i, -float(recruited.diagnostics["k_flux"])
        )

    def test_calibrated_capacity_is_multiplied_by_genotype_not_overwritten(self) -> None:
        parameters = AE4Parameters(activity_scale=3.25, na_weight=1.0, k_weight=0.0)
        mechanism = AE4Mechanism("C7_112", parameters)
        full = mechanism.evaluate(CONTEXT, activity_scale=1.0)
        half = mechanism.evaluate(CONTEXT, activity_scale=0.5)
        knockout = mechanism.evaluate(CONTEXT, activity_scale=0.0)
        direct = evaluate_ae4("C7_112", CONTEXT, parameters)
        self.assertAlmostEqual(full.cycle_flux, direct.cycle_flux)
        self.assertAlmostEqual(half.cycle_flux, 0.5 * full.cycle_flux)
        self.assertEqual(knockout.cycle_flux, 0.0)
        self.assertEqual(knockout.na_i, 0.0)
        self.assertEqual(knockout.k_i, 0.0)
        self.assertEqual(knockout.cl_i, 0.0)
        self.assertEqual(knockout.hco3_i, 0.0)

    def test_c6_mass_action_and_c8_mass_action_are_declared_aliases(self) -> None:
        """Requested family labels must not be counted as independent evidence."""

        parameters = AE4Parameters(
            activity_scale=2.5,
            na_weight=1.5,
            k_weight=1.6,
            equilibrium_constant=1.0,
        )
        for suffix in ("112", "123", "213"):
            self.assertEqual(
                CANDIDATE_SPECS[f"C8_{suffix}_MA"].equivalent_to,
                f"C6_{suffix}",
            )
            c6 = evaluate_ae4(f"C6_{suffix}", CONTEXT, parameters)
            c8 = evaluate_ae4(f"C8_{suffix}_MA", CONTEXT, parameters)
            self.assertEqual(c6, c8, suffix)

    def test_assay_scale_k_partition_cannot_close_na_only_chassis_vector(self) -> None:
        """Comparable K capacity and exact historical Na-only closure do not overlap."""

        k_fraction = 1.6 / (1.5 + 1.6)
        signature = np.array(
            [-(1.0 - k_fraction), -k_fraction, 1.0, -2.0], dtype=float
        )
        fit = fit_nonnegative_activity_scale(signature)
        self.assertAlmostEqual(fit.max_abs_error, 0.020125768940380106)
        self.assertGreater(fit.max_abs_error, 2.0e6 * 1.0e-8)


class ThermodynamicAudit(unittest.TestCase):
    THERMODYNAMIC_IDS = tuple(
        spec.candidate_id
        for spec in candidate_specs(include_rejected=False)
        if spec.thermodynamic_sign_controlled
    )

    def test_reversible_candidates_zero_and_reverse_at_declared_equilibrium(self) -> None:
        for candidate_id in self.THERMODYNAMIC_IDS:
            stoichiometry = CANDIDATE_SPECS[candidate_id].stoichiometry
            assert stoichiometry is not None
            for cation in ("na", "k"):
                parameters = AE4Parameters(
                    na_weight=1.0 if cation == "na" else 0.0,
                    k_weight=1.0 if cation == "k" else 0.0,
                )
                equilibrium = equilibrium_context(
                    CONTEXT,
                    cation,
                    stoichiometry,
                    equilibrium_constant=parameters.equilibrium_constant,
                )
                lower = replace(equilibrium, hco3_i=equilibrium.hco3_i * 0.999)
                upper = replace(equilibrium, hco3_i=equilibrium.hco3_i * 1.001)
                at_zero = evaluate_ae4(candidate_id, equilibrium, parameters)
                below = evaluate_ae4(candidate_id, lower, parameters)
                above = evaluate_ae4(candidate_id, upper, parameters)
                comparison_scale = max(
                    1.0, abs(below.cycle_flux), abs(above.cycle_flux)
                )
                label = f"{candidate_id}:{cation}"
                self.assertLess(
                    abs(at_zero.cycle_flux),
                    2e-12 * comparison_scale,
                    label,
                )
                self.assertLess(below.cycle_flux, 0.0, label)
                self.assertGreater(above.cycle_flux, 0.0, label)
                self.assertTrue(at_zero.diagnostics["thermodynamic_valid"], label)
                self.assertTrue(below.diagnostics["thermodynamic_valid"], label)
                self.assertTrue(above.diagnostics["thermodynamic_valid"], label)

    def test_empirical_candidates_are_not_mislabelled_thermodynamic(self) -> None:
        for candidate_id in ("C1", "C2", "C3", "C4"):
            self.assertFalse(CANDIDATE_SPECS[candidate_id].thermodynamic_sign_controlled)
            result = evaluate_ae4(candidate_id, CONTEXT)
            self.assertFalse(result.diagnostics["thermodynamic_valid"])

    def test_hill_values_are_constraints_not_proof_of_transport_count(self) -> None:
        defaults = AE4Parameters()
        self.assertEqual(defaults.ec50_na, 49.0)
        self.assertEqual(defaults.ec50_k, 62.0)
        self.assertEqual(defaults.hill_na, 2.0)
        self.assertEqual(defaults.hill_k, 1.8)
        c4 = CANDIDATE_SPECS["C4"]
        self.assertIn("one cation is transported", " ".join(c4.new_assumptions))


class StimulusAudit(unittest.TestCase):
    """Source-supported PKA extensions must remain knockout-blind and local."""

    PKA_IDS = (
        "C1_PKA_COMMON",
        "C1_PKA_K_RECRUIT",
        "C3_PKA_COMMON",
        "C6_112_PKA_COMMON",
        "C6_123_PKA_COMMON",
        "C6_213_PKA_COMMON",
    )

    def test_assay_fold_envelope_is_predeclared_and_contains_no_phenotype(self) -> None:
        self.assertEqual(PKA_FOLD_ASSAY_VALUES, (1.25, 1.6, 3.0))
        self.assertNotIn(0.65, PKA_FOLD_ASSAY_VALUES)
        for candidate_id in self.PKA_IDS:
            spec = CANDIDATE_SPECS[candidate_id]
            self.assertEqual(spec.free_parameters, ("activity_scale",))
            self.assertIsNotNone(spec.core_candidate_id)
            self.assertNotEqual(spec.stimulus_mode, "none")

        parameters = AE4Parameters(na_weight=1.5, k_weight=1.6)
        self.assertFalse(
            _assay_admissible(CANDIDATE_SPECS["C1_PKA_COMMON"], parameters)
        )
        self.assertFalse(
            _assay_admissible(CANDIDATE_SPECS["C1_PKA_K_RECRUIT"], parameters)
        )
        for suffix in ("112", "123", "213"):
            self.assertTrue(
                _assay_admissible(
                    CANDIDATE_SPECS[f"C6_{suffix}_PKA_COMMON"], parameters
                )
            )

    def test_zero_pka_input_restores_each_static_core_exactly(self) -> None:
        parameters = AE4Parameters(activity_scale=2.75, pka_stimulated_fold=3.0)
        for candidate_id in self.PKA_IDS:
            spec = CANDIDATE_SPECS[candidate_id]
            assert spec.core_candidate_id is not None
            core = evaluate_ae4(spec.core_candidate_id, CONTEXT, parameters)
            resting = evaluate_ae4(
                candidate_id, CONTEXT, parameters, stimulus_activation=0.0
            )
            self.assertEqual(resting, core, candidate_id)

    def test_common_pka_gate_scales_source_without_changing_signature(self) -> None:
        parameters = AE4Parameters(activity_scale=2.75, pka_stimulated_fold=3.0)
        for candidate_id in (
            "C1_PKA_COMMON",
            "C3_PKA_COMMON",
            "C6_112_PKA_COMMON",
            "C6_123_PKA_COMMON",
            "C6_213_PKA_COMMON",
        ):
            resting = evaluate_ae4(candidate_id, CONTEXT, parameters)
            stimulated = evaluate_ae4(
                candidate_id, CONTEXT, parameters, stimulus_activation=1.0
            )
            for field in ("cycle_flux", "na_i", "k_i", "cl_i", "hco3_i"):
                self.assertAlmostEqual(
                    getattr(stimulated, field),
                    3.0 * getattr(resting, field),
                    msg=f"{candidate_id}:{field}",
                )

    def test_knockout_annihilates_basal_and_stimulus_sources(self) -> None:
        parameters = AE4Parameters(activity_scale=2.75, pka_stimulated_fold=3.0)
        for candidate_id in self.PKA_IDS:
            contribution = evaluate_ae4(
                candidate_id,
                CONTEXT,
                parameters,
                activity_scale=0.0,
                stimulus_activation=1.0,
            )
            self.assertEqual(
                (
                    contribution.cycle_flux,
                    contribution.na_i,
                    contribution.k_i,
                    contribution.cl_i,
                    contribution.hco3_i,
                ),
                (0.0, 0.0, 0.0, 0.0, 0.0),
                candidate_id,
            )

    def test_pka_input_changes_only_the_ae4_source_at_fixed_state(self) -> None:
        parameters = AE4Parameters(activity_scale=2.75, pka_stimulated_fold=3.0)
        non_ae4_diagnostics = (
            "j_nak",
            "j_nkcc1",
            "j_ae2",
            "j_nhe1",
            "j_buffer",
            "j_cl_signed",
            "j_k",
            "j_t_na",
            "j_t_k",
            "q_a",
            "q_b",
            "q_t",
            "q_total",
            "v_a",
            "v_b",
        )
        state = FixedChassis.initial_state()
        for candidate_id in self.PKA_IDS:
            chassis = FixedChassis(AE4Mechanism(candidate_id, parameters))
            raw_rest, rest = chassis.evaluate(
                state, calcium=0.55, pka_activation=0.0, scenario=WT
            )
            raw_stim, stimulated = chassis.evaluate(
                state, calcium=0.55, pka_activation=1.0, scenario=WT
            )
            self.assertEqual(
                rest["parameter_fingerprint"],
                stimulated["parameter_fingerprint"],
                candidate_id,
            )
            for name in non_ae4_diagnostics:
                self.assertEqual(rest[name], stimulated[name], f"{candidate_id}:{name}")
            # Lumen Na, lumen K, and cell-volume equations do not receive an
            # AE4 amount source directly at the same state.
            np.testing.assert_array_equal(
                raw_rest[:3], raw_stim[:3], err_msg=candidate_id
            )

            ko_rest, _ = chassis.evaluate(
                state, calcium=0.55, pka_activation=0.0, scenario=AE4_KNOCKOUT
            )
            ko_stim, _ = chassis.evaluate(
                state, calcium=0.55, pka_activation=1.0, scenario=AE4_KNOCKOUT
            )
            np.testing.assert_array_equal(ko_rest, ko_stim, err_msg=candidate_id)

    def test_pka_gate_preserves_thermodynamic_zero_and_sign(self) -> None:
        for suffix, stoichiometry in (
            ("112", (1, 1, 2)),
            ("123", (1, 2, 3)),
            ("213", (2, 1, 3)),
        ):
            candidate_id = f"C6_{suffix}_PKA_COMMON"
            for cation in ("na", "k"):
                parameters = AE4Parameters(
                    na_weight=1.0 if cation == "na" else 0.0,
                    k_weight=1.0 if cation == "k" else 0.0,
                    pka_stimulated_fold=3.0,
                )
                equilibrium = equilibrium_context(
                    CONTEXT, cation, stoichiometry
                )
                below = replace(
                    equilibrium, hco3_i=equilibrium.hco3_i * 0.999
                )
                above = replace(
                    equilibrium, hco3_i=equilibrium.hco3_i * 1.001
                )
                at_zero = evaluate_ae4(
                    candidate_id,
                    equilibrium,
                    parameters,
                    stimulus_activation=1.0,
                )
                negative = evaluate_ae4(
                    candidate_id,
                    below,
                    parameters,
                    stimulus_activation=1.0,
                )
                positive = evaluate_ae4(
                    candidate_id,
                    above,
                    parameters,
                    stimulus_activation=1.0,
                )
                label = f"{candidate_id}:{cation}"
                scale = max(1.0, abs(negative.cycle_flux), abs(positive.cycle_flux))
                self.assertLess(abs(at_zero.cycle_flux), 2e-12 * scale, label)
                self.assertLess(negative.cycle_flux, 0.0, label)
                self.assertGreater(positive.cycle_flux, 0.0, label)

    def test_generic_pka_wrapper_scales_branch_diagnostics_and_source(self) -> None:
        parameters = AE4Parameters(activity_scale=0.01)
        for candidate_id in ("C7_213", "C8_213_SAT"):
            core = AE4Mechanism(candidate_id, parameters)
            wrapper = CommonPkaGateMechanism(core, 3.0)
            basal = wrapper.evaluate(CONTEXT, stimulus_activation=0.0)
            static = core.evaluate(CONTEXT, stimulus_activation=0.0)
            stimulated = wrapper.evaluate(CONTEXT, stimulus_activation=1.0)
            knockout = wrapper.evaluate(
                CONTEXT,
                activity_scale=0.0,
                stimulus_activation=1.0,
            )
            for field in ("cycle_flux", "na_i", "k_i", "cl_i", "hco3_i"):
                self.assertEqual(
                    getattr(basal, field), getattr(static, field), candidate_id
                )
                self.assertAlmostEqual(
                    getattr(stimulated, field),
                    3.0 * getattr(basal, field),
                    msg=f"{candidate_id}:{field}",
                )
                self.assertEqual(getattr(knockout, field), 0.0, candidate_id)
            for branch in ("na_flux", "k_flux"):
                self.assertAlmostEqual(
                    float(stimulated.diagnostics[branch]),
                    3.0 * float(basal.diagnostics[branch]),
                    msg=f"{candidate_id}:{branch}",
                )
            self.assertEqual(
                stimulated.na_i,
                -(CANDIDATE_SPECS[candidate_id].stoichiometry or (0, 0, 0))[1]
                * float(stimulated.diagnostics["na_flux"]),
            )
            self.assertEqual(
                stimulated.k_i,
                -(CANDIDATE_SPECS[candidate_id].stoichiometry or (0, 0, 0))[1]
                * float(stimulated.diagnostics["k_flux"]),
            )
        with self.assertRaisesRegex(ValueError, "at least one"):
            CommonPkaGateMechanism(AE4Mechanism("C7_213", parameters), 0.9)

    def test_capacity_band_and_stimulus_selectors_have_no_phenotype_interface(self) -> None:
        self.assertEqual(
            tuple(inspect.signature(rest_admissible_capacity_envelopes).parameters),
            ("roots",),
        )
        self.assertEqual(
            tuple(inspect.signature(stimulus_activation_comparisons).parameters),
            (),
        )
        self.assertEqual(CAPACITY_PROFILE_LOG10_OFFSETS[0], -8.0)
        self.assertEqual(CAPACITY_PROFILE_LOG10_OFFSETS[-1], 4.0)
        self.assertEqual(
            CAPACITY_INTERVAL_SELECTION_FRACTIONS,
            (0.001, 0.25, 0.5, 0.75, 0.999),
        )
        for function in (
            rest_admissible_capacity_envelopes,
            stimulus_activation_comparisons,
        ):
            self.assertNotIn("0.65", inspect.getsource(function))


class FixedChassisAudit(unittest.TestCase):
    def test_historical_baseline_and_conservation_are_independently_reproduced(self) -> None:
        chassis = FixedChassis(HistoricalSodiumOnlyAE4())
        raw, diagnostics = chassis.evaluate(
            chassis.initial_state(), calcium=0.05, scenario=WT
        )
        self.assertLess(np.max(np.abs(raw)), 5.3e-12)
        self.assertAlmostEqual(diagnostics["v_a"], -50.24, delta=2e-12)
        self.assertAlmostEqual(diagnostics["v_b"], -62.8, delta=2e-12)
        self.assertAlmostEqual(diagnostics["cell_volume_pL"], 1.3, delta=2e-14)
        for name, residual in diagnostics["checks"].items():
            self.assertLess(abs(float(residual)), 6e-14, name)

    def test_non_ae4_fingerprint_does_not_depend_on_candidate_or_genotype(self) -> None:
        baseline = FixedChassis()
        fingerprint = baseline.parameter_fingerprint
        for candidate_id in ("C1", "C2", "C5", "C7_123", "C8_213_SAT"):
            chassis = FixedChassis(AE4Mechanism(candidate_id))
            self.assertEqual(chassis.parameter_fingerprint, fingerprint)
            for scenario in (WT, AE4_KNOCKOUT, AE2_KNOCKOUT):
                _, diagnostics = chassis.evaluate(
                    chassis.initial_state(), calcium=0.05, scenario=scenario
                )
                self.assertEqual(diagnostics["parameter_fingerprint"], fingerprint)

    def test_ae4_knockout_rhs_is_candidate_independent(self) -> None:
        """The zero genotype multiplier must annihilate every AE4 law exactly."""

        states = (
            FixedChassis.initial_state(),
            np.array([145.0, 4.0, 1.1, 30.0, 110.0, 45.0, 15.0]),
        )
        candidate_ids = ("C1", "C2", "C3", "C4", "C5", "C7_123", "C8_213_SAT")
        for state in states:
            reference_raw, reference_diagnostics = FixedChassis(
                AE4Mechanism(candidate_ids[0])
            ).evaluate(state, calcium=0.55, scenario=AE4_KNOCKOUT)
            self.assertEqual(reference_diagnostics["ae4_cycle_flux"], 0.0)
            for candidate_id in candidate_ids[1:]:
                raw, diagnostics = FixedChassis(AE4Mechanism(candidate_id)).evaluate(
                    state, calcium=0.55, scenario=AE4_KNOCKOUT
                )
                np.testing.assert_array_equal(raw, reference_raw, err_msg=candidate_id)
                self.assertEqual(diagnostics["ae4_cycle_flux"], 0.0, candidate_id)

    def test_primary_calibrator_has_no_knockout_target_interface(self) -> None:
        """The held-out phenotype cannot be supplied to the Round-1 calibrator."""

        signature = inspect.signature(calibrate_candidate)
        self.assertEqual(
            tuple(signature.parameters),
            ("candidate_id", "parameters", "external_hco3_mM"),
        )
        calibration = calibrate_candidate("C1")
        np.testing.assert_allclose(
            calibration.scale_fit.fitted_contribution,
            HISTORICAL_WT_AE4_CONTRIBUTION,
            rtol=0.0,
            atol=2e-14,
        )

    def test_chassis_rejects_a_non_electroneutral_ae4_source(self) -> None:
        class ChargedLaw:
            def __call__(self, _environment: object, scale: float) -> AE4Balance:
                return AE4Balance(
                    cycle_flux=scale,
                    na_i=0.0,
                    k_i=0.0,
                    cl_i=scale,
                    hco3_i=-2.0 * scale,
                )

        chassis = FixedChassis(ChargedLaw())
        with self.assertRaisesRegex(ValueError, "not electroneutral"):
            chassis.evaluate(chassis.initial_state(), calcium=0.05, scenario=WT)

    def test_chassis_has_no_runtime_dependency_on_the_immutable_archive(self) -> None:
        source = inspect.getsource(sys.modules[FixedChassis.__module__])
        self.assertNotIn("archive/", source)
        self.assertNotIn("scipy.io.loadmat", source)

    def test_unit_boundary_forbids_absolute_flow_claim(self) -> None:
        ledger = dimension_ledger()
        self.assertEqual(ledger["absolute_flow"], "not certified")
        self.assertIn("ratios", ledger["valid_primary_comparison"])


class MovingRestAudit(unittest.TestCase):
    def test_rest_solver_has_no_knockout_phenotype_interface(self) -> None:
        parameters = tuple(inspect.signature(solve_resting_equilibrium).parameters)
        self.assertEqual(
            parameters,
            (
                "chassis",
                "mechanism_id",
                "activity_scale",
                "scenario",
                "initial_states",
                "initial_start_seed",
                "calibration_mode",
            ),
        )

    def test_candidate_specific_roots_do_not_hide_rest_data_failures(self) -> None:
        # Re-solve the core moving roots directly with one neutral start.  The
        # production 20-start/C6 continuation is checked in its frozen CSV so
        # this fast audit test does not repeat hundreds of nonlinear solves.
        rows = []
        for candidate_id in ("C1", "C2", "C3"):
            calibration = calibrate_candidate(candidate_id)
            chassis = candidate_chassis(calibration)
            rows.append(
                solve_resting_equilibrium(
                    chassis,
                    mechanism_id=candidate_id,
                    activity_scale=calibration.parameters.activity_scale,
                    initial_states=(chassis.initial_state(),),
                )
            )
        self.assertEqual(
            {row.mechanism_id for row in rows},
            {"C1", "C2", "C3"},
        )
        self.assertTrue(all(row.solver_success for row in rows))
        self.assertTrue(all(row.active_bound_count == 0 for row in rows))
        self.assertTrue(all(row.initial_start_count == 1 for row in rows))
        self.assertTrue(all(row.converged_root_count >= 1 for row in rows))
        self.assertTrue(all(row.max_abs_raw_rhs <= 1e-8 for row in rows))
        self.assertTrue(all(not row.knockout_flow_target_read_by_solver for row in rows))
        self.assertEqual(
            {row.mechanism_id for row in rows if row.quantified_rest_gate_pass},
            {"C1"},
        )
        fingerprints = {row.non_ae4_chassis_fingerprint for row in rows}
        self.assertEqual(len(fingerprints), 1)

    def test_common_ae4_null_rest_contradicts_primary_cl_and_ph(self) -> None:
        """No AE4 law can repair a phenotype after its full source is zero."""

        knockout_rest = np.array(
            [
                118.56527670449144,
                5.5593461776616655,
                65.66070346655056,
                21.132282246430893,
                123.70751315363782,
                48.782215027024655,
                58.90444835877707,
            ]
        )
        chassis = candidate_chassis(calibrate_candidate("C1"))
        raw, diagnostics = chassis.evaluate(
            knockout_rest,
            calcium=0.05,
            pka_activation=0.0,
            scenario=AE4_KNOCKOUT,
        )
        self.assertLess(np.max(np.abs(raw)), 1e-8)
        cl_z = (float(knockout_rest[5]) - 36.50) / 1.60
        ph_z = (float(diagnostics["pH_i"]) - 6.89) / 0.02
        self.assertGreater(cl_z, 7.0)
        self.assertGreater(ph_z, 20.0)

        # Changing the candidate cannot change this AE4-null vector field.
        other = FixedChassis(AE4Mechanism("C8_213_SAT"))
        other_raw, _ = other.evaluate(
            knockout_rest,
            calcium=0.05,
            pka_activation=1.0,
            scenario=AE4_KNOCKOUT,
        )
        np.testing.assert_array_equal(other_raw, raw)


class FrozenArtifactAudit(unittest.TestCase):
    """Hard-gate the production freeze without repeating the costly search."""

    RESULTS = ROOT / "results/12_ae4_mechanism_reconstruction"

    @classmethod
    def _load(cls, filename: str) -> object:
        path = cls.RESULTS / filename
        if not path.is_file():
            raise AssertionError(f"required frozen artifact is missing: {path}")
        return json.loads(path.read_text(encoding="utf-8"))

    def test_capacity_points_are_ko_blind_and_rest_admissible(self) -> None:
        payload = self._load("rest_capacity_envelope.json")
        assert isinstance(payload, dict)
        self.assertFalse(payload["knockout_target_read_by_variant_selector"])
        serialized = json.dumps(payload, sort_keys=True).lower()
        self.assertNotIn("heldout", serialized)
        self.assertNotIn("knockout_flow_ratio", serialized)
        points = [
            row
            for row in payload["rows"]
            if row["profile_kind"] == "frozen_interval_prediction_point"
        ]
        self.assertEqual(len(points), 15)
        self.assertEqual(
            {row["core_candidate_id"] for row in points},
            {"C6_213", "C7_213", "C8_213_SAT"},
        )
        expected_positions = {
            "near_lower_boundary",
            "lower_quartile",
            "midpoint",
            "upper_quartile",
            "near_upper_boundary",
        }
        for core_id in ("C6_213", "C7_213", "C8_213_SAT"):
            core = [row for row in points if row["core_candidate_id"] == core_id]
            self.assertEqual(
                {row["capacity_position"] for row in core}, expected_positions
            )
            self.assertEqual(len({row["capacity_interval_lower"] for row in core}), 1)
            self.assertEqual(len({row["capacity_interval_upper"] for row in core}), 1)
        for row in points:
            self.assertTrue(row["capacity_gate_pass"])
            self.assertTrue(row["quantified_rest_gate_pass"])
            self.assertFalse(row["knockout_target_read_by_variant_selector"])
            self.assertFalse(row["knockout_flow_target_read_by_solver"])
            self.assertGreaterEqual(row["capacity_gate_margin"], -1e-10)
            self.assertLessEqual(row["max_abs_raw_rhs"], 1e-8)

    def test_pka_predictions_fail_in_the_opposite_direction_without_leakage(self) -> None:
        payload = self._load("pka_activation_envelope.json")
        assert isinstance(payload, dict)
        protocol = payload["anti_leak_protocol"]
        self.assertTrue(protocol["folds_and_capacity_points_frozen_before_heldout_join"])
        self.assertFalse(protocol["knockout_target_read_by_calibrator"])
        self.assertFalse(protocol["knockout_target_read_by_variant_selector"])
        rows = payload["rows"]
        source_rows = [
            row
            for row in rows
            if row["source_complete_under_new_volume_sensitivity"]
        ]
        self.assertEqual(len(source_rows), 30)
        self.assertEqual(
            {row["core_candidate_id"] for row in source_rows},
            {"C6_213", "C7_213", "C8_213_SAT"},
        )
        self.assertEqual({row["pka_stimulated_fold"] for row in source_rows}, {1.25, 3.0})
        fingerprint = FixedChassis().parameter_fingerprint
        for row in source_rows:
            self.assertEqual(row["simulation_status"], "pass")
            self.assertTrue(row["resting_quantified_gate_pass"])
            self.assertTrue(row["assay_admissible"])
            self.assertEqual(row["thermodynamic_status"], "pass")
            self.assertFalse(row["source_complete"])
            self.assertFalse(row["knockout_target_read_by_calibrator"])
            self.assertFalse(row["knockout_target_read_by_variant_selector"])
            self.assertNotIn("0.65", row["calibration_target_values"])
            self.assertNotIn("0.65", row["new_modelling_assumptions"])
            self.assertGreater(row["ae4_ko_endpoint_flow_ratio"], 1.0)
            self.assertGreater(row["ae4_ko_integrated_flow_ratio"], 1.0)
            self.assertLessEqual(row["max_conservation_residual"], 1e-10)
            self.assertLessEqual(row["max_thermodynamic_violation"], 1e-12)
            self.assertEqual(row["non_ae4_chassis_fingerprint"], fingerprint)
        self.assertAlmostEqual(
            min(row["ae4_ko_endpoint_flow_ratio"] for row in source_rows),
            1.0137777583543537,
            delta=5e-12,
        )
        self.assertAlmostEqual(
            min(row["ae4_ko_integrated_flow_ratio"] for row in source_rows),
            1.0209390020865694,
            delta=5e-12,
        )

    def test_genotype_specific_rest_targets_are_not_conflated(self) -> None:
        payload = self._load("genotype_resting_equilibria.json")
        assert isinstance(payload, dict)
        rows = payload["rows"]
        self.assertNotIn("COMMON_AE2_KO_REST", {row["mechanism_id"] for row in rows})
        ae4 = next(row for row in rows if row["mechanism_id"] == "COMMON_AE4_KO_REST")
        self.assertEqual(ae4["genotype_cl_target_mM"], 36.50)
        self.assertEqual(ae4["genotype_pH_target"], 6.89)
        self.assertGreater(ae4["genotype_cl_target_z"], 7.0)
        self.assertGreater(ae4["genotype_pH_target_z"], 20.0)
        self.assertFalse(ae4["genotype_rest_gate_pass"])

        c1_ae2 = next(row for row in rows if row["mechanism_id"] == "C1_AE2_KO_REST")
        self.assertEqual(c1_ae2["genotype_cl_target_mM"], 54.50)
        self.assertEqual(c1_ae2["genotype_pH_target"], 6.95)
        self.assertAlmostEqual(c1_ae2["genotype_cl_target_z"], -2.764, delta=0.002)
        self.assertFalse(c1_ae2["genotype_rest_gate_pass"])

        survivors = [
            row
            for row in rows
            if row["mechanism_id"].startswith(
                ("C6_213_", "C7_213_", "C8_213_SAT_")
            )
        ]
        self.assertEqual(len(survivors), 15)
        for row in survivors:
            self.assertEqual(row["scenario"], "ae2_knockout")
            self.assertEqual(row["genotype_cl_target_mM"], 54.50)
            self.assertEqual(row["genotype_pH_target"], 6.95)
            self.assertLess(row["genotype_cl_target_z"], -4.0)
            self.assertFalse(row["genotype_rest_gate_pass"])
            self.assertLessEqual(row["max_abs_raw_rhs"], 1e-8)

    def test_summary_and_solver_crosscheck_encode_the_scoped_exclusion(self) -> None:
        summary = self._load("summary.json")
        assert isinstance(summary, dict)
        self.assertEqual(
            summary["substantive_conclusion"],
            "AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED",
        )
        self.assertEqual(summary["heldout_AE4_KO_flow_ratio"], 0.65)
        self.assertEqual(summary["heldout_AE4_KO_flow_ratio_SEM"], 0.047)
        self.assertEqual(summary["flow_compatible_source_rest_rows"], 0)
        self.assertEqual(summary["all_hard_gate_success_rows"], 0)
        self.assertEqual(
            summary["hard_gate_failures"][
                "candidate_specific_AE2_KO_rest_passing_rows"
            ],
            0,
        )
        self.assertFalse(
            summary["hard_gate_failures"][
                "candidate_independent_AE4_KO_rest_gate_pass"
            ]
        )
        self.assertEqual(
            summary["source_rest_thermo_PKA_rows_under_new_volume_sensitivity"],
            30,
        )
        self.assertEqual(
            summary["strict_source_complete_rows_without_new_volume_sensitivity"],
            0,
        )

        solver = self._load("solver_crosscheck.json")
        assert isinstance(solver, dict)
        crosscheck = solver["independent_representative_213_BDF_Radau"]
        self.assertEqual(crosscheck["cores"], ["C6_213", "C7_213", "C8_213_SAT"])
        self.assertEqual(
            crosscheck["capacity_points_per_core"]
            * len(crosscheck["cores"])
            * len(crosscheck["folds"]),
            15,
        )
        self.assertLessEqual(
            crosscheck["reported_upper_bound_absolute_endpoint_ratio_difference"],
            3e-10,
        )
        self.assertLessEqual(
            crosscheck["reported_upper_bound_absolute_integrated_ratio_difference"],
            3e-9,
        )
        self.assertLessEqual(
            crosscheck["reported_upper_bound_conservation_residual"], 5e-14
        )
        self.assertEqual(crosscheck["maximum_thermodynamic_violation"], 0.0)
        self.assertIn("15 rows", crosscheck["provenance"])
        self.assertIn("30 positive", crosscheck["provenance"])

        source_crosscheck = solver[
            "independent_source_positive_30_BDF_vs_production_Radau"
        ]
        self.assertEqual(
            source_crosscheck["cores"],
            ["C6_213", "C7_213", "C8_213_SAT"],
        )
        self.assertEqual(source_crosscheck["capacity_points_per_core"], 5)
        self.assertEqual(source_crosscheck["folds"], [1.25, 3.0])
        self.assertEqual(source_crosscheck["row_count"], 30)
        self.assertEqual(
            source_crosscheck["row_count"],
            source_crosscheck["capacity_points_per_core"]
            * len(source_crosscheck["cores"])
            * len(source_crosscheck["folds"]),
        )
        self.assertTrue(source_crosscheck["all_BDF_endpoint_ratios_above_one"])
        self.assertTrue(source_crosscheck["all_BDF_integrated_ratios_above_one"])
        self.assertAlmostEqual(
            source_crosscheck["minimum_BDF_endpoint_ratio"],
            1.013777757953735,
            delta=1e-14,
        )
        self.assertAlmostEqual(
            source_crosscheck["minimum_BDF_integrated_ratio"],
            1.0209390004383792,
            delta=1e-14,
        )
        self.assertLessEqual(
            source_crosscheck["maximum_absolute_endpoint_ratio_difference"],
            1.36e-9,
        )
        self.assertLessEqual(
            source_crosscheck["maximum_absolute_integrated_ratio_difference"],
            3.371e-9,
        )
        self.assertLessEqual(
            source_crosscheck["maximum_production_Radau_conservation_residual"],
            5e-14,
        )
        self.assertEqual(
            source_crosscheck["maximum_thermodynamic_violation"], 0.0
        )
        self.assertIn("30 rows", source_crosscheck["provenance"])


if __name__ == "__main__":
    unittest.main()
