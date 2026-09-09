"""Target-firewalled tests for genotype continuation and dynamics."""

from __future__ import annotations

import inspect
import itertools
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model.calibration import build_wt_model
from modern_full_model.camp_pka import R3PkaRegulatedFraction
from modern_full_model.genotype_evaluation import (
    DEFAULT_EXPRESSION_GRID,
    FROZEN_EFFECTIVE_CALCIUM_UM,
    FROZEN_NKCC1_FULLY_ACTIVATED_MULTIPLIER,
    FROZEN_REGULATORY_MEMBER_ID,
    GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
    GenotypeRoot,
    _cluster_roots,
    build_frozen_genotype_dynamic_model,
    continue_ae2_expression,
    evaluate_dual_solver_genotype_pair,
    frozen_dynamic_contract_metadata,
    genotype_with_expression,
)
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model
from modern_full_model.run_g5_cross_root_profile import load_retained_wt_roots
from modern_full_model.validation import StimulusArm, attach_basal_regulation
import modern_full_model.genotype_evaluation as genotype_evaluation


class GenotypeEvaluationFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repository = Path(__file__).resolve().parents[1]
        roots = load_retained_wt_roots(
            repository / "results" / "13B_modern_full_model"
        )
        cls.frozen_root = next(root for root in roots if root.selected_reference)
        cls.resting_model = build_wt_model(
            cls.frozen_root.variant, cls.frozen_root.calibration
        )
        cls.dynamic_model = build_frozen_genotype_dynamic_model(cls.resting_model)
        cls.ae2_continuation = continue_ae2_expression(
            cls.dynamic_model,
            cls.frozen_root.core_state,
            expression_grid=(1.0, 0.5, 0.0),
            local_start_count=1,
            max_nfev=1000,
        )


class TestFrozenGenotypeContract(GenotypeEvaluationFixture):
    def test_contract_is_exact_r3_g125_reference_timing(self) -> None:
        metadata = frozen_dynamic_contract_metadata()
        self.assertEqual(
            metadata["regulatory_member_id"], FROZEN_REGULATORY_MEMBER_ID
        )
        self.assertEqual(metadata["regulatory_family"], "R3")
        self.assertEqual(metadata["regulatory_fully_activated_multiplier"], 1.25)
        self.assertEqual(metadata["regulatory_timing_label"], "TREFERENCE")
        self.assertEqual(metadata["effective_calcium_uM"], FROZEN_EFFECTIVE_CALCIUM_UM)
        self.assertEqual(
            metadata["nkcc1_fully_activated_multiplier"],
            FROZEN_NKCC1_FULLY_ACTIVATED_MULTIPLIER,
        )
        self.assertEqual(metadata["stimulus_arm"], StimulusArm.CCH_IPR.value)

        self.assertIsInstance(self.dynamic_model, StimulatedNkcc1Model)
        ae4_regulation = (
            self.dynamic_model.regulatory_model.ae4_regulatory_model
        )
        self.assertIsInstance(ae4_regulation, R3PkaRegulatedFraction)
        self.assertEqual(ae4_regulation.tau_pka_s, 10.0)
        self.assertAlmostEqual(ae4_regulation.forward_regulation_rate_s, 1.0 / 30.0)
        self.assertAlmostEqual(ae4_regulation.reverse_regulation_rate_s, 1.0 / 30.0)

    def test_dynamic_clone_preserves_every_parameter_and_other(self) -> None:
        self.assertIs(self.dynamic_model.parameters, self.resting_model.parameters)
        self.assertEqual(
            self.dynamic_model.parameters.geometry.cell_impermeant_osmoles_fmol,
            self.resting_model.parameters.geometry.cell_impermeant_osmoles_fmol,
        )
        self.assertEqual(
            self.dynamic_model.stimulus.stimulated_calcium_uM,
            FROZEN_EFFECTIVE_CALCIUM_UM,
        )

    def test_default_expression_grid_is_predeclared_one_to_exact_zero(self) -> None:
        self.assertEqual(len(DEFAULT_EXPRESSION_GRID), 21)
        self.assertEqual(DEFAULT_EXPRESSION_GRID[0], 1.0)
        self.assertEqual(DEFAULT_EXPRESSION_GRID[-1], 0.0)
        self.assertTrue(
            all(
                right < left
                for left, right in zip(
                    DEFAULT_EXPRESSION_GRID, DEFAULT_EXPRESSION_GRID[1:]
                )
            )
        )


class TestChargeManifoldGenotypeContinuation(GenotypeEvaluationFixture):
    def test_ae2_continuation_solves_ten_for_ten_with_no_parameter_fit(self) -> None:
        report = self.ae2_continuation
        self.assertTrue(report.completed_to_exact_zero, report.branch_failure)
        self.assertEqual(report.independent_equation_count, 10)
        self.assertEqual(report.fitted_coordinate_count, 10)
        self.assertEqual(report.fitted_parameter_count, 0)
        self.assertEqual(len(report.steps), 3)
        self.assertTrue(all(step.connected_to_previous for step in report.steps))
        self.assertTrue(
            all(
                step.selected_root is not None
                and step.selected_root.normalized_jacobian_rank == 10
                and step.selected_root.normalized_jacobian_nullity == 0
                and step.selected_root.passes_numerical_gate
                for step in report.steps
            )
        )

    def test_connected_null_root_has_exact_zero_ae2_flux_and_fixed_other(self) -> None:
        report = self.ae2_continuation
        root = report.connected_null_root
        self.assertIsNotNone(root)
        assert root is not None
        self.assertEqual(root.expression, 0.0)
        self.assertTrue(root.exact_deleted_transport_zero)
        self.assertEqual(
            report.fixed_other_impermeant_osmoles_fmol,
            self.dynamic_model.parameters.geometry.cell_impermeant_osmoles_fmol,
        )
        full = attach_basal_regulation(self.dynamic_model, root.core_state)
        evaluation = self.dynamic_model.evaluate(
            0.0, full, genotype=genotype_with_expression("AE2", 0.0)
        )
        self.assertEqual(evaluation.diagnostics.homeostasis.ae2_inward_fmol_s, 0.0)

    def test_ae4_deletion_removes_every_conserved_ae4_source_exactly(self) -> None:
        full = attach_basal_regulation(
            self.dynamic_model, self.frozen_root.core_state
        )
        evaluation = self.dynamic_model.evaluate(
            0.0, full, genotype=genotype_with_expression("AE4", 0.0)
        )
        contribution = evaluation.diagnostics.ae4
        self.assertEqual(contribution.na_cell_fmol_s, 0.0)
        self.assertEqual(contribution.k_cell_fmol_s, 0.0)
        self.assertEqual(contribution.cl_cell_fmol_s, 0.0)
        self.assertEqual(contribution.hco3_cell_fmol_s, 0.0)
        self.assertEqual(contribution.transported_charge_fmol_s, 0.0)


class TestDeterministicGenotypeRootClustering(unittest.TestCase):
    @staticmethod
    def _root(label: str, coordinate: float, residual: float) -> GenotypeRoot:
        coordinates = (coordinate, *([0.0] * 9))
        return GenotypeRoot(
            root_id="",
            transporter="AE4",
            expression=0.5,
            member_start_ids=(label,),
            coordinates=coordinates,
            core_state=tuple(float(index + 1) for index in range(12)),
            max_abs_scaled_independent_rhs=residual,
            max_abs_amount_rhs_fmol_s=0.0,
            max_abs_volume_rhs_pL_s=0.0,
            max_abs_omitted_charge_rhs_fmol_equivalent_s=0.0,
            max_abs_regulatory_rhs_s_inv=0.0,
            max_abs_charge_residual_fmol=0.0,
            max_abs_current_residual_A=0.0,
            normalized_jacobian_singular_values=tuple([1.0] * 10),
            normalized_jacobian_rank=10,
            normalized_jacobian_nullity=0,
            boundary_hits=(),
            exact_deleted_transport_zero=True,
            passes_numerical_gate=True,
            gate_failures=(),
        )

    def test_chain_connected_components_and_ids_are_permutation_invariant(self) -> None:
        roots = (
            self._root("a", 0.0, 0.03),
            self._root("b", 0.075, 0.01),
            self._root("c", 0.15, 0.02),
            self._root("d", 0.8, 0.04),
        )
        signatures = set()
        for permutation in itertools.permutations(roots):
            clustered = _cluster_roots(
                permutation,
                expression=0.5,
                transporter="AE4",
                span=np.ones(10),
                relative_tolerance=0.1,
            )
            signatures.add(
                tuple(
                    (root.root_id, root.coordinates, root.member_start_ids)
                    for root in clustered
                )
            )
        self.assertEqual(len(signatures), 1)
        signature = signatures.pop()
        self.assertEqual([row[0] for row in signature], ["AE4_E0.50000000:B00", "AE4_E0.50000000:B01"])
        self.assertEqual(signature[0][2], ("a", "b", "c"))
        self.assertEqual(signature[0][1][0], 0.075)


class TestTargetFreeDynamicEvaluation(GenotypeEvaluationFixture):
    def test_short_ae2_pair_runs_with_radau_and_bdf_and_generic_ratios(self) -> None:
        root = self.ae2_continuation.connected_null_root
        self.assertIsNotNone(root)
        assert root is not None
        result = evaluate_dual_solver_genotype_pair(
            self.dynamic_model,
            self.frozen_root.core_state,
            root.core_state,
            transporter="AE2",
            time_s=np.asarray((0.0, 1.0e-6, 1.0, 5.0, 10.0)),
            landmarks_s=(1.0, 5.0, 10.0),
        )
        self.assertTrue(result.both_solver_pairs_pass)
        self.assertEqual(
            set(result.evaluations), {"production_radau", "production_bdf"}
        )
        for summary in result.evaluations.values():
            self.assertEqual(
                set(summary.cumulative_ratio_at_landmarks),
                {"1_s", "5_s", "10_s"},
            )
            self.assertEqual(
                set(summary.flow_ratio_at_landmarks),
                {"1_s", "5_s", "10_s"},
            )
            self.assertIn("cell_cl_mM", summary.genotype_endpoint)
            self.assertIn("cell_ph", summary.genotype_endpoint)
            self.assertIn("cell_volume_pL", summary.genotype_endpoint)
        for trajectory in result.genotype_trajectories.values():
            self.assertEqual(trajectory.max_abs_deleted_transport_flux_fmol_s, 0.0)
            self.assertTrue(trajectory.numerical_gate_pass)

    def test_cross_solver_disagreement_above_existing_contract_tolerance_fails(self) -> None:
        root = self.ae2_continuation.connected_null_root
        self.assertIsNotNone(root)
        assert root is not None
        with patch.object(
            genotype_evaluation,
            "_relative_difference",
            return_value=2.0 * GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        ):
            result = evaluate_dual_solver_genotype_pair(
                self.dynamic_model,
                self.frozen_root.core_state,
                root.core_state,
                transporter="AE2",
                time_s=np.asarray((0.0, 1.0e-6, 1.0, 5.0, 10.0)),
                landmarks_s=(1.0, 5.0, 10.0),
            )
        self.assertFalse(result.cross_solver_gate_pass)
        self.assertFalse(result.both_solver_pairs_pass)
        self.assertEqual(
            result.cross_solver_relative_tolerance,
            GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        )

    def test_module_has_no_target_file_reader_or_phenotype_gate(self) -> None:
        source = inspect.getsource(genotype_evaluation)
        self.assertNotIn("heldout_targets.csv", source)
        self.assertNotIn("validation_targets.csv", source)
        self.assertNotIn("phenotype_gate_pass", source)
        self.assertNotIn("target_lower", source)
        self.assertNotIn("target_upper", source)


if __name__ == "__main__":
    unittest.main()
