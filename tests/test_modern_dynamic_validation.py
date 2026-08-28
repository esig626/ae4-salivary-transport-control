"""Focused pre-reveal WT dynamic-validation tests for Task 13B."""

from __future__ import annotations

import inspect
import unittest

import numpy as np

from modern_full_model.model import ModernFullModel
from modern_full_model.camp_pka import default_regulatory_families
from modern_full_model.states import amount_charge_equivalents_fmol
from modern_full_model.transporters import AE4Parameters
from modern_full_model import validation


class TestSecretagogueSeparation(unittest.TestCase):
    def test_combined_protocol_has_a_basal_left_point_and_physical_duration(self) -> None:
        protocol = validation.SecretagogueProtocol()
        self.assertEqual(protocol.duration_s, 600.0)
        self.assertEqual(protocol(0.0).calcium_uM, protocol.resting_calcium_uM)
        self.assertEqual(protocol(0.0).beta_input, 0.0)
        self.assertEqual(protocol(1.0e-9).calcium_uM, protocol.stimulated_calcium_uM)
        self.assertEqual(protocol(1.0e-9).beta_input, 1.0)
        self.assertEqual(protocol(600.0).beta_input, 1.0)
        self.assertEqual(protocol(600.0 + 1.0e-9).beta_input, 0.0)

    def test_cch_and_ipr_enter_distinct_model_inputs(self) -> None:
        cch = validation.SecretagogueProtocol(arm=validation.StimulusArm.CCH_ONLY)
        ipr = validation.SecretagogueProtocol(arm=validation.StimulusArm.IPR_ONLY)
        cch_stimulus = cch(1.0)
        ipr_stimulus = ipr(1.0)
        self.assertGreater(cch_stimulus.calcium_uM, cch.resting_calcium_uM)
        self.assertEqual(cch_stimulus.beta_input, 0.0)
        self.assertEqual(ipr_stimulus.calcium_uM, ipr.resting_calcium_uM)
        self.assertGreater(ipr_stimulus.beta_input, 0.0)
        self.assertEqual(cch.cch_input(1.0), 0.3)
        self.assertEqual(cch.ipr_input(1.0), 0.0)
        self.assertEqual(ipr.cch_input(1.0), 0.0)
        self.assertEqual(ipr.ipr_input(1.0), 5.0)

    def test_seconds_grid_resolves_the_stimulus_right_limit(self) -> None:
        grid = validation.physical_time_grid(duration_s=10.0, sample_step_s=1.0)
        self.assertEqual(grid[0], 0.0)
        self.assertEqual(grid[1], 1.0e-6)
        self.assertEqual(grid[-1], 10.0)
        self.assertTrue(np.all(np.diff(grid) > 0.0))

    def test_solver_absolute_tolerances_follow_coordinate_units(self) -> None:
        names = ("na_i_fmol", "volume_i_pL", "regulatory_fraction")
        values = validation.PRODUCTION_RADAU.atol_vector(names)
        self.assertEqual(values[0], validation.PRODUCTION_RADAU.amount_atol_fmol)
        self.assertEqual(values[1], validation.PRODUCTION_RADAU.volume_atol_pL)
        self.assertEqual(
            values[2], validation.PRODUCTION_RADAU.regulatory_atol_fraction
        )

    def test_pre_reveal_regulatory_ensemble_covers_family_gain_and_timing_axes(self) -> None:
        ensemble = validation.pre_reveal_regulatory_ensemble()
        self.assertEqual(len(ensemble), 20)
        self.assertEqual({member.family for member in ensemble}, {"R0", "R1", "R2", "R3"})
        self.assertEqual(
            {member.fully_activated_multiplier for member in ensemble}, {1.25, 1.70}
        )
        for family in ("R1", "R2", "R3"):
            self.assertEqual(
                {member.timing_label for member in ensemble if member.family == family},
                {"TFAST", "TREFERENCE", "TSLOW"},
            )
        self.assertTrue(
            all("NOT_FIT" in member.kinetic_evidence_status for member in ensemble if member.family != "R0")
        )
        by_family_timing = {
            (member.family, member.timing_label): dict(member.timing_parameters)
            for member in ensemble
            if member.gain_label == "G125_P21_PROSE"
        }
        self.assertEqual(
            [
                by_family_timing[("R1", timing)]["tau_activation_s"]
                for timing in ("TFAST", "TREFERENCE", "TSLOW")
            ],
            [10.0, 30.0, 90.0],
        )
        self.assertEqual(
            [
                (
                    by_family_timing[("R2", timing)]["tau_camp_s"],
                    by_family_timing[("R2", timing)]["tau_activation_s"],
                )
                for timing in ("TFAST", "TREFERENCE", "TSLOW")
            ],
            [(3.0, 10.0), (10.0, 30.0), (30.0, 90.0)],
        )


class TestNearbyStatesAndSolvers(unittest.TestCase):
    def setUp(self) -> None:
        self.ae4 = AE4Parameters.equal_branch_attempts(
            carrier_amount_fmol=1.0e-4,
            common_cl_attempt_rate_s=1.0,
            loaded_attempt_rate_s=1.0,
        )
        self.model = ModernFullModel(
            stimulus=validation.SecretagogueProtocol(),
            ae4_parameters=self.ae4,
        )
        self.state = self.model.initial_state()

    def test_nearby_states_preserve_bulk_charge_exactly_to_numerical_precision(self) -> None:
        variants = validation.nearby_charge_preserving_states(self.model, self.state)
        self.assertEqual(
            set(variants),
            {
                "baseline",
                "cell_volume_down",
                "cell_volume_up",
                "cell_nacl_down",
                "cell_nacl_up",
                "lumen_volume_down",
                "lumen_volume_up",
                "lumen_nacl_down",
                "lumen_nacl_up",
            },
        )
        charges = []
        for candidate in variants.values():
            decoded = self.model.layout.decode(candidate)
            charges.append(
                amount_charge_equivalents_fmol(
                    decoded.cell,
                    fixed_anion_equivalents_fmol=(
                        self.model.parameters.geometry.fixed_cell_anion_equivalents_fmol
                    ),
                )
            )
        np.testing.assert_allclose(charges, charges[0], rtol=0.0, atol=1.0e-12)

    def test_short_physical_wt_trajectory_reproduces_with_radau_and_bdf(self) -> None:
        grid = validation.physical_time_grid(duration_s=0.05, sample_step_s=0.01)
        radau_spec = validation.SolverSpecification(
            "test_radau", "Radau", 1.0e-7, 1.0e-10, 1.0e-12, 1.0e-10, 0.01
        )
        bdf_spec = validation.SolverSpecification(
            "test_bdf", "BDF", 1.0e-7, 1.0e-10, 1.0e-12, 1.0e-10, 0.01
        )
        radau = validation.simulate_wt(
            self.model,
            self.state,
            family="R0_TEST",
            initial_condition="baseline",
            solver=radau_spec,
            time_s=grid,
        )
        bdf = validation.simulate_wt(
            self.model,
            self.state,
            family="R0_TEST",
            initial_condition="baseline",
            solver=bdf_spec,
            time_s=grid,
        )
        self.assertTrue(radau.success, radau.message)
        self.assertTrue(bdf.success, bdf.message)
        self.assertTrue(radau.positive_core)
        self.assertTrue(bdf.positive_core)
        comparison = validation.solver_comparison(radau, bdf)
        self.assertTrue(comparison["both_success"])
        self.assertLess(comparison["max_relative_state_difference"], 1.0e-5)

    def test_onset_split_integrates_every_nested_regulatory_family(self) -> None:
        grid = validation.physical_time_grid(duration_s=0.05, sample_step_s=0.01)
        specification = validation.SolverSpecification(
            "nested_smoke", "Radau", 1.0e-7, 1.0e-10, 1.0e-12, 1.0e-10, 0.01
        )
        endpoints = {}
        for family, regulatory_model in default_regulatory_families().items():
            model = ModernFullModel(
                stimulus=validation.SecretagogueProtocol(),
                regulatory_model=regulatory_model,
                ae4_parameters=self.ae4,
            )
            trajectory = validation.simulate_wt(
                model,
                model.initial_state(),
                family=family,
                initial_condition="baseline",
                solver=specification,
                time_s=grid,
            )
            self.assertTrue(trajectory.success, (family, trajectory.message))
            self.assertEqual(trajectory.capacity_multiplier[0], 1.0)
            self.assertGreater(trajectory.capacity_multiplier[-1], 1.0)
            endpoints[family] = trajectory.capacity_multiplier[-1]
        self.assertAlmostEqual(endpoints["R2"], endpoints["R3"], places=13)


class TestFlowScaleAndFirewall(unittest.TestCase):
    def test_wt_flow_scale_reports_its_implied_cell_count(self) -> None:
        grid = np.asarray((0.0, 300.0, 600.0))
        flow = np.asarray((0.1, 0.2, 0.3))
        dummy = validation.Trajectory(
            family="R0",
            initial_condition="baseline",
            solver=validation.PRODUCTION_RADAU,
            success=True,
            message="ok",
            time_s=grid,
            state_names=("x",),
            states=np.ones((1, 3)),
            flow_pL_s=flow,
            cumulative_flow_pL=np.asarray((0.0, 45.0, 120.0)),
            cell_na_mM=np.ones(3),
            cell_k_mM=np.ones(3),
            cell_cl_mM=np.ones(3),
            cell_ph=np.ones(3) * 7.0,
            cell_volume_pL=np.ones(3),
            capacity_multiplier=np.ones(3),
            max_abs_conservation_residuals={
                name: 0.0 for name in validation.CONSERVATION_RESIDUAL_UNITS
            },
        )
        scale = validation.fit_wt_flow_scale(dummy)
        self.assertAlmostEqual(scale.model_statistic_pL_s, 0.2)
        self.assertAlmostEqual(scale.convert(np.asarray((0.2,)))[0], 9.5)
        self.assertAlmostEqual(
            scale.effective_cell_count * 6.0e-5,
            scale.multiplier_uL_min_per_pL_s,
        )

    def test_flow_shape_gate_cannot_be_rescued_by_one_multiplier(self) -> None:
        grid = np.arange(0.0, 601.0, 60.0)

        def dummy(flow: np.ndarray) -> validation.Trajectory:
            size = grid.size
            return validation.Trajectory(
                family="R1",
                initial_condition="baseline",
                solver=validation.PRODUCTION_RADAU,
                success=True,
                message="ok",
                time_s=grid,
                state_names=("x",),
                states=np.ones((1, size)),
                flow_pL_s=flow,
                cumulative_flow_pL=np.r_[0.0, np.cumsum((flow[:-1] + flow[1:]) * 30.0)],
                cell_na_mM=np.ones(size),
                cell_k_mM=np.ones(size),
                cell_cl_mM=np.ones(size),
                cell_ph=np.ones(size) * 7.0,
                cell_volume_pL=np.ones(size),
                capacity_multiplier=np.ones(size),
                max_abs_conservation_residuals={
                    name: 0.0 for name in validation.CONSERVATION_RESIDUAL_UNITS
                },
            )

        flat = validation.wt_flow_envelope_feasibility(
            dummy(np.ones(grid.size) * 0.2)
        )
        declining = validation.wt_flow_envelope_feasibility(
            dummy(np.linspace(0.2, 0.1, grid.size))
        )
        self.assertTrue(flat.compatible)
        self.assertFalse(declining.compatible)
        self.assertGreater(
            declining.model_max_min_ratio, declining.envelope_max_min_ratio
        )

    def test_pre_reveal_manifest_rejects_forbidden_target_references(self) -> None:
        validation.assert_pre_reveal_payload({"status": "WT_ONLY", "seconds": 600})
        with self.assertRaises(ValueError):
            validation.assert_pre_reveal_payload(
                {"input": "results/13B_modern_full_model/heldout_targets.csv"}
            )

    def test_validation_source_contains_no_knockout_simulation(self) -> None:
        source = inspect.getsource(validation)
        self.assertNotIn("AE4_NULL", source)
        self.assertNotIn("secretion_ratios", source)
        self.assertNotIn("genotype=AE4", source)


if __name__ == "__main__":
    unittest.main()
