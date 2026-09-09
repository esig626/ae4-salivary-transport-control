"""Anti-leak and fixed-chassis tests for Task 12 model comparison."""

from __future__ import annotations

import math
import inspect
import unittest
from dataclasses import replace

import numpy as np

from ae4_mechanism_reconstruction.chassis import (
    AE4_KNOCKOUT,
    BASELINE_STATE,
    WT,
    FixedChassis,
)
from ae4_mechanism_reconstruction.comparison import (
    BathOverrideMechanism,
    CommonPkaGateMechanism,
    HELDOUT_TARGET_NAMES,
    HISTORICAL_WT_AE4_CONTRIBUTION,
    PATHWAY_SIGNATURES,
    PRIMARY_CALIBRATION_TARGET_NAMES,
    _round2_variant_definitions,
    baseline_context,
    calibrate_candidate,
    candidate_chassis,
    cosine_alignment,
    fit_nonnegative_activity_scale,
    independently_constrained_parameters,
    pathway_localization,
    projected_residual_fraction,
    stimulus_activation_comparisons,
    transported_k_fraction_bound,
)
from ae4_mechanism_reconstruction.mechanisms import AE4Mechanism
from ae4_mechanism_reconstruction.run_analysis import (
    GENOTYPE_REST_TARGETS,
    _downsample_mask,
    _genotype_rest_rows,
    _jsonable,
)


class KnockoutBlindCalibrationTests(unittest.TestCase):
    def test_calibrator_exposes_no_knockout_argument_or_target(self) -> None:
        self.assertTrue(all("KO" not in name for name in PRIMARY_CALIBRATION_TARGET_NAMES))
        self.assertTrue(all("KO" in name for name in HELDOUT_TARGET_NAMES))
        calibration = calibrate_candidate("C1")
        self.assertNotIn("knockout", repr(calibration.scale_fit).lower())

    def test_nonnegative_scale_fit_is_exact_for_historical_signature(self) -> None:
        fit = fit_nonnegative_activity_scale(HISTORICAL_WT_AE4_CONTRIBUTION)
        self.assertEqual(fit.activity_scale, 1.0)
        self.assertEqual(fit.max_abs_error, 0.0)

    def test_c1_wt_scale_and_closure_are_recovered(self) -> None:
        calibration = calibrate_candidate("C1")
        self.assertAlmostEqual(calibration.parameters.activity_scale, 18.950260722154685)
        self.assertEqual(calibration.scale_fit.max_abs_error, 0.0)
        chassis = candidate_chassis(calibration)
        raw, diagnostic = chassis.evaluate(BASELINE_STATE, calcium=0.05, scenario=WT)
        self.assertLess(np.max(np.abs(raw)), 2e-12)
        self.assertAlmostEqual(diagnostic["ae4_cycle_flux"], 0.03888333867647999)

    def test_published_partition_fails_fixed_wt_closure(self) -> None:
        calibration = calibrate_candidate("C2")
        self.assertGreater(calibration.scale_fit.max_abs_error, 0.032)
        chassis = candidate_chassis(calibration)
        raw, _ = chassis.evaluate(BASELINE_STATE, calcium=0.05, scenario=WT)
        self.assertGreater(np.max(np.abs(raw)), 1e-3)

    def test_default_reversible_family_has_nonnegative_fit_at_zero(self) -> None:
        calibration = calibrate_candidate("C6_112")
        self.assertEqual(calibration.parameters.activity_scale, 0.0)
        self.assertGreater(calibration.scale_fit.max_abs_error, 0.07)

    def test_genotype_multiplier_preserves_capacity_then_zeros_it(self) -> None:
        calibration = calibrate_candidate("C1")
        chassis = candidate_chassis(calibration)
        _, wt = chassis.evaluate(BASELINE_STATE, calcium=0.05, scenario=WT)
        _, ko = chassis.evaluate(BASELINE_STATE, calcium=0.05, scenario=AE4_KNOCKOUT)
        self.assertAlmostEqual(wt["ae4_cycle_flux"], 0.03888333867647999)
        self.assertEqual(ko["ae4_cycle_flux"], 0.0)
        self.assertEqual(ko["ae4_na_source"], 0.0)
        self.assertEqual(ko["ae4_k_source"], 0.0)

    def test_ko_rhs_is_candidate_independent(self) -> None:
        states = (
            BASELINE_STATE,
            np.array((115.0, 7.0, 25.0, 30.0, 110.0, 40.0, 15.0)),
        )
        reference = candidate_chassis(calibrate_candidate("C1"))
        comparison = candidate_chassis(calibrate_candidate("C8_123_SAT"))
        for state in states:
            np.testing.assert_array_equal(
                reference.raw_rhs(150.0, state, AE4_KNOCKOUT),
                comparison.raw_rhs(150.0, state, AE4_KNOCKOUT),
            )

    def test_pka_protocol_is_independent_of_calcium_and_ko_annihilates_it(self) -> None:
        parameters = replace(
            independently_constrained_parameters(), pka_stimulated_fold=3.0
        )
        calibration = calibrate_candidate("C1_PKA_COMMON", parameters=parameters)
        law = BathOverrideMechanism(
            AE4Mechanism("C1_PKA_COMMON", calibration.parameters), 21.0
        )
        default = FixedChassis(law)
        activated = FixedChassis(law, pka_protocol=lambda _time: 1.0)
        default_wt = default.raw_rhs(150.0, BASELINE_STATE, WT)
        activated_wt = activated.raw_rhs(150.0, BASELINE_STATE, WT)
        self.assertFalse(np.array_equal(default_wt, activated_wt))
        np.testing.assert_array_equal(
            default.raw_rhs(150.0, BASELINE_STATE, AE4_KNOCKOUT),
            activated.raw_rhs(150.0, BASELINE_STATE, AE4_KNOCKOUT),
        )

    def test_generic_pka_wrapper_scales_branch_diagnostics_and_validates_fold(self) -> None:
        parameters = independently_constrained_parameters(activity_scale=0.01)
        core = AE4Mechanism("C7_213", parameters)
        wrapper = CommonPkaGateMechanism(core, 3.0)
        context = baseline_context()
        basal = wrapper.evaluate(context, stimulus_activation=0.0)
        gated = wrapper.evaluate(context, stimulus_activation=1.0)
        self.assertAlmostEqual(gated.na_i, 3.0 * basal.na_i)
        self.assertAlmostEqual(gated.k_i, 3.0 * basal.k_i)
        self.assertAlmostEqual(
            float(gated.diagnostics["na_flux"]),
            3.0 * float(basal.diagnostics["na_flux"]),
        )
        self.assertAlmostEqual(
            float(gated.diagnostics["k_flux"]),
            3.0 * float(basal.diagnostics["k_flux"]),
        )
        with self.assertRaises(ValueError):
            CommonPkaGateMechanism(core, 0.9)


class IterativeProtocolTests(unittest.TestCase):
    def test_json_schema_preserves_boolean_type(self) -> None:
        self.assertIs(_jsonable(True), True)
        self.assertIs(_jsonable(np.bool_(False)), False)

    def test_sparse_trajectory_export_keeps_protocol_landmarks(self) -> None:
        times = np.asarray(
            (0.0, 10.0, 50.0, 100.0, 100.0 + 1e-7, 120.0, 150.0, 200.0)
        )
        np.testing.assert_array_equal(
            _downsample_mask(times),
            np.asarray((True, False, True, True, True, False, True, True)),
        )

    def test_ae2_rest_uses_direct_genotype_target_and_is_c1_specific(self) -> None:
        self.assertNotIn("COMMON_AE2_KO_REST", GENOTYPE_REST_TARGETS)
        rows = _genotype_rest_rows(
            (
                {
                    "mechanism_id": "C1_AE2_KO_REST",
                    "cl_i_mM": 49.52490516299572,
                    "pH_i": 6.91409,
                    "max_abs_raw_rhs": 5.2e-14,
                },
            )
        )
        self.assertEqual(rows[0]["genotype_cl_target_mM"], 54.50)
        self.assertEqual(rows[0]["genotype_pH_target"], 6.95)
        self.assertAlmostEqual(rows[0]["genotype_cl_target_z"], -2.7639415761)
        self.assertFalse(rows[0]["genotype_rest_gate_pass"])
        self.assertAlmostEqual(rows[0]["own_control_cl_target_z"], -2.152830465)

    def test_round2_grid_is_predeclared_without_phenotype(self) -> None:
        variants = _round2_variant_definitions()
        self.assertEqual(len(variants), 768)
        serialized = repr(variants)
        self.assertNotIn("0.65", serialized)
        self.assertNotIn("knockout", serialized.lower())

    def test_stimulus_runner_has_no_heldout_argument(self) -> None:
        signature = inspect.signature(stimulus_activation_comparisons)
        self.assertEqual(tuple(signature.parameters), ())

    def test_k_transport_fraction_bound(self) -> None:
        expected = 1e-8 / 0.03888333867647999
        self.assertAlmostEqual(transported_k_fraction_bound(), expected)
        self.assertLess(expected, 2.6e-7)

    def test_one_component_localization_limits_are_exact(self) -> None:
        correction = (-1.0, 1.0, 0.0, 0.0)
        pump = PATHWAY_SIGNATURES["NaK_ATPase"]
        self.assertAlmostEqual(cosine_alignment(correction, pump), 5 / math.sqrt(26))
        self.assertAlmostEqual(
            projected_residual_fraction(correction, pump), 1 / math.sqrt(26)
        )
        self.assertAlmostEqual(
            projected_residual_fraction(correction, PATHWAY_SIGNATURES["NHE1"]),
            1 / math.sqrt(2),
        )
        self.assertEqual(pathway_localization(correction)[0]["pathway"], "NaK_ATPase")


if __name__ == "__main__":
    unittest.main()
