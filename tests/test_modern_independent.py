"""Independent numerical-transcription tests for Task 13B.

These tests contain no held-out phenotype value.  Agreement tests compare the
two source assemblies at local states; the root and trajectory tests run only
the independent equations.
"""

from __future__ import annotations

import inspect
import csv
from collections import defaultdict
from dataclasses import replace
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np

from modern_full_model import acid_base
from modern_full_model.independent import (
    FREE_CORE_INDICES,
    IndependentCalibrationDefinition,
    IndependentConstantStimulus,
    IndependentN1NkccRegulation,
    IndependentN2NkccRegulation,
    IndependentSecretagogueProtocol,
    IndependentWholeCell,
    free_charge_coordinates,
    independent_r0_regulation,
    independent_r1_regulation,
    independent_r2_regulation,
    independent_r3_regulation,
    independent_physical_time_grid,
    independent_nearby_charge_preserving_states,
    independent_speciate,
    independent_total_alkalinity_mM,
    state_from_free_charge_coordinates,
)
from modern_full_model.model import ModernFullModel
from modern_full_model.parameters import FullModelParameters
from modern_full_model.transporters import AE4Parameters
from modern_full_model.validation import SecretagogueProtocol


class IndependentFixture(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = FullModelParameters()
        self.ae4 = AE4Parameters.equal_branch_attempts(
            carrier_amount_fmol=1.0e-4,
            common_cl_attempt_rate_s=1.0,
            loaded_attempt_rate_s=1.0,
        )
        self.production = ModernFullModel(
            parameters=self.parameters,
            ae4_parameters=self.ae4,
        )
        self.reference = self.production.initial_state()
        self.independent = IndependentWholeCell(
            self.parameters,
            ae4_parameters=self.ae4,
        )


class TestIndependentAcidBase(IndependentFixture):
    def test_second_speciation_recovers_pH_and_species(self) -> None:
        ph = 7.17
        ta = independent_total_alkalinity_mM(
            ph, 19.0, 27.0, self.parameters.acid_base.cell_buffer_pka,
            self.parameters,
        )
        reproduced = independent_speciate(
            19.0, ta, 27.0, self.parameters.acid_base.cell_buffer_pka,
            self.parameters,
        )
        production = acid_base.speciate(
            total_carbon_mM=19.0,
            total_alkalinity_mM_value=ta,
            buffer_total_mM=27.0,
            buffer_pka=self.parameters.acid_base.cell_buffer_pka,
            parameters=self.parameters.acid_base,
        )
        self.assertAlmostEqual(reproduced.ph, ph, places=11)
        self.assertAlmostEqual(reproduced.ph, production.ph, places=11)
        self.assertAlmostEqual(reproduced.hco3_mM, production.hco3_mM, places=11)
        self.assertLess(abs(reproduced.alkalinity_residual_mM), 1.0e-11)


class TestIndependentEquationTranscription(IndependentFixture):
    def test_source_does_not_import_or_call_production_model(self) -> None:
        import modern_full_model.independent as independent_module

        source = inspect.getsource(independent_module)
        self.assertNotIn("from .model import", source)
        self.assertNotIn("from .model import", source)
        self.assertNotIn("heldout_targets", source)

    def test_full_rhs_matches_at_deterministic_local_states(self) -> None:
        rng = np.random.default_rng(1309)
        for _ in range(8):
            state = self.reference.copy()
            # Common amount-volume scaling preserves concentrations and pH;
            # paired cation/anion changes then exercise the source laws.
            state[:6] *= float(rng.uniform(0.90, 1.10))
            state[6:12] *= float(rng.uniform(0.90, 1.10))
            delta_i = float(rng.uniform(-0.01, 0.01)) * min(state[0], state[2])
            delta_l = float(rng.uniform(-0.01, 0.01)) * min(state[6], state[8])
            state[0] += delta_i
            state[2] += delta_i
            state[6] += delta_l
            state[8] += delta_l
            expected = self.production.evaluate(0.0, state)
            actual = self.independent.evaluate(0.0, state)
            np.testing.assert_allclose(actual.rhs, expected.rhs, rtol=3.0e-13, atol=3.0e-14)
            self.assertAlmostEqual(
                actual.cell_speciation.ph,
                expected.diagnostics.observables.cell_acid_base.ph,
                places=11,
            )
            self.assertAlmostEqual(
                actual.voltages_V["apical"],
                expected.diagnostics.membranes.v_apical_V,
                places=13,
            )

    def test_independent_current_charge_carbon_and_water_close(self) -> None:
        result = self.independent.evaluate(0.0, self.reference)
        self.assertLess(max(abs(value) for value in result.current_residuals_A.values()), 1.0e-22)
        self.assertLess(
            max(abs(value) for value in result.conservation_residuals.values()),
            1.0e-12,
        )
        self.assertGreaterEqual(result.fluxes_fmol_s["ae4_entropy_over_r"], -1.0e-14)


class TestIndependentChargeManifold(IndependentFixture):
    def test_free_coordinates_recover_exact_electroneutrality(self) -> None:
        free = free_charge_coordinates(self.reference)
        rebuilt = state_from_free_charge_coordinates(
            free,
            fixed_cell_anion_equivalents_fmol=(
                self.parameters.geometry.fixed_cell_anion_equivalents_fmol
            ),
        )
        self.assertEqual(tuple(free.shape), (len(FREE_CORE_INDICES),))
        self.assertAlmostEqual(
            rebuilt[0] + rebuilt[1] - rebuilt[2] - rebuilt[4]
            - self.parameters.geometry.fixed_cell_anion_equivalents_fmol,
            0.0,
            places=13,
        )
        self.assertAlmostEqual(rebuilt[6] + rebuilt[7] - rebuilt[8] - rebuilt[10], 0.0, places=13)

    def test_nearby_states_preserve_cell_and_lumen_bulk_charge(self) -> None:
        variants = independent_nearby_charge_preserving_states(
            self.reference, fraction=1.0e-3
        )
        self.assertEqual(len(variants), 9)
        base_cell = self.reference[0] + self.reference[1] - self.reference[2] - self.reference[4]
        base_lumen = self.reference[6] + self.reference[7] - self.reference[8] - self.reference[10]
        for label, candidate in variants.items():
            with self.subTest(label=label):
                self.assertAlmostEqual(
                    candidate[0] + candidate[1] - candidate[2] - candidate[4],
                    base_cell,
                    places=12,
                )
                self.assertAlmostEqual(
                    candidate[6] + candidate[7] - candidate[8] - candidate[10],
                    base_lumen,
                    places=12,
                )

    def test_uncalibrated_starting_parameters_fail_independent_root_gate_honestly(self) -> None:
        result = self.independent.refine_resting_root(
            self.reference,
            max_nfev=120,
            relative_bound=0.20,
        )
        self.assertFalse(result.success)
        self.assertGreater(result.max_abs_scaled_residual, 1.0e-3)
        self.assertEqual(result.jacobian_rank, 10)
        self.assertEqual(result.jacobian_nullity, 0)
        self.assertTrue(np.all(result.state > 0.0))

    def test_provisional_old_water_g2_root_is_rejected_after_unit_correction(self) -> None:
        # Provisional candidate G2_BALANCED_APICAL_K_BIASED:M00 was generated
        # before the Palk hydraulic coefficients were converted into physical
        # pL/s/mOsm units.  With the corrected WaterParameters it must not be
        # accepted.  The production build helper and production RHS are not
        # called in this independent rejection check.
        base = FullModelParameters()
        configured = replace(
            base,
            homeostasis=replace(
                base.homeostasis,
                nhe1_capacity_fmol_s=0.007531999570279085,
            ),
            membranes=replace(
                base.membranes,
                apical_pump_fraction=0.10,
                apical_k_fraction=0.10,
                g_k_total_S=base.membranes.g_k_total_S * 10.0,
                g_cl_apical_S=base.membranes.g_cl_apical_S * 10.0,
                g_basolateral_background_S=(
                    base.membranes.g_basolateral_background_S * 10.0
                ),
            ),
            geometry=replace(
                base.geometry,
                cell_impermeant_osmoles_fmol=134.04162965081832,
            ),
        )
        ae4 = SimpleNamespace(
            carrier_amount_fmol=0.15429333062794132,
            common_cl_attempt_rate_s=1.0,
            na_loaded_attempt_rate_s=0.2,
            k_loaded_attempt_rate_s=1.8,
            cooperative_gate=None,
        )

        class RestingStimulus:
            def __call__(self, _time_s: float):
                return {"calcium_uM": 0.058, "beta_input": 0.0}

        model = IndependentWholeCell(
            configured,
            ae4_parameters=ae4,
            stimulus=RestingStimulus(),
        )
        candidate = np.asarray(
            (
                19.53167330533108,
                144.06877521788118,
                65.13000000000010,
                8.031318683024280,
                20.40868631726734,
                1.300000000000001,
                15.546681498386224,
                0.4115542761886416,
                14.690123292377084,
                1.3485678704627913,
                1.2681124821977823,
                0.10142974249419931,
            )
        )
        result = model.refine_resting_root(
            candidate,
            relative_bound=0.01,
            max_nfev=200,
            row_scales=np.asarray(
                (0.05, 0.05, 0.02, 0.02, 1.0e-4, 0.02, 0.02, 0.02, 0.02, 1.0e-4)
            ),
        )
        self.assertFalse(result.success)
        self.assertGreater(result.max_abs_scaled_residual, 1.0e-3)
        self.assertEqual(result.jacobian_rank, 10)
        self.assertEqual(result.jacobian_nullity, 0)
        self.assertEqual(configured.water.apical_hydraulic_pL_s_mOsm, 4.32e-3)
        self.assertEqual(configured.water.basolateral_hydraulic_pL_s_mOsm, 5.15e-2)
        self.assertEqual(configured.water.paracellular_hydraulic_pL_s_mOsm, 2.60e-4)

    def test_every_frozen_retained_r0_root_is_independently_reproduced(self) -> None:
        repository = Path(__file__).resolve().parents[1]
        results = repository / "results" / "13B_modern_full_model"
        with (results / "wt_root_table.csv").open(newline="", encoding="utf-8") as handle:
            retained = [
                row
                for row in csv.DictReader(handle)
                if row["status"] == "WT_REST_PASS_R0_CONTROL_ONLY"
            ]
        state_records: dict[str, dict[str, float]] = defaultdict(dict)
        with (results / "wt_root_states.csv").open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                state_records[row["root_id"]][row["state_name"]] = float(row["value"])
        parameter_records: dict[str, dict[str, float]] = defaultdict(dict)
        with (results / "wt_parameter_candidates.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            for row in csv.DictReader(handle):
                parameter_records[row["root_id"]][row["parameter"]] = float(row["value"])

        self.assertEqual(len(retained), 9)
        self.assertIn(
            "G2_BALANCED_APICAL_K_BIASED_G10:M00",
            {row["root_id"] for row in retained},
        )
        row_scales = np.asarray(
            (0.05, 0.05, 0.02, 0.02, 1.0e-4, 0.02, 0.02, 0.02, 0.02, 1.0e-4)
        )
        for row in retained:
            root_id = row["root_id"]
            values = parameter_records[root_id]
            definition = IndependentCalibrationDefinition(
                mixed_bath_na_attempt_fraction=values[
                    "mixed_bath_na_attempt_fraction"
                ],
                apical_pump_fraction=values["apical_pump_fraction"],
                apical_k_fraction=values["apical_k_fraction"],
                ae4_carrier_amount_fmol_at_unit_rate_gauge=values[
                    "ae4_carrier_amount_fmol_at_unit_rate_gauge"
                ],
                nhe1_capacity_fmol_s=values["nhe1_capacity_fmol_s"],
                common_membrane_conductance_scale=values[
                    "common_membrane_conductance_scale"
                ],
                cell_other_impermeant_osmoles_fmol=values[
                    "cell_other_impermeant_osmoles_fmol"
                ],
            )
            parameters, ae4 = definition.reconstruct_parameters()
            independent = IndependentWholeCell(
                parameters,
                ae4_parameters=ae4,
                stimulus=IndependentConstantStimulus(0.058),
            )
            candidate = np.asarray(
                [state_records[root_id][name] for name in independent.state_names]
            )
            supplied = independent.evaluate(0.0, candidate)
            self.assertLess(
                max(abs(supplied.rhs[index]) for index in (0, 1, 2, 3, 4, 6, 7, 8, 9, 10)),
                5.0e-13,
                root_id,
            )
            self.assertLess(
                max(abs(supplied.rhs[index]) for index in (5, 11)),
                5.0e-12,
                root_id,
            )
            reproduced = independent.refine_resting_root(
                candidate,
                relative_bound=0.01,
                max_nfev=200,
                row_scales=row_scales,
            )
            self.assertTrue(reproduced.success, root_id)
            self.assertLess(reproduced.max_abs_scaled_residual, 1.0e-7, root_id)
            self.assertEqual(reproduced.jacobian_rank, 10, root_id)
            self.assertEqual(reproduced.jacobian_nullity, 0, root_id)
            np.testing.assert_allclose(
                reproduced.state,
                candidate,
                rtol=0.0,
                atol=1.0e-10,
                err_msg=root_id,
            )

    def test_selected_wt_candidate_manifest_has_the_frozen_digest(self) -> None:
        path = (
            Path(__file__).resolve().parents[1]
            / "results"
            / "13B_modern_full_model"
            / "wt_selected_candidate.json"
        )
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        self.assertEqual(
            digest,
            "41b91c09a69b6f50f1ea96ec49df9c56681bcc21ada8bc68b515f1b34a0eea44",
        )


class TestIndependentRegulatoryDynamics(IndependentFixture):
    def test_dynamic_capacity_and_two_stiff_solvers_reproduce(self) -> None:
        regulation = independent_r1_regulation(
            tau_activation_s=8.0,
            basal_multiplier=1.0,
            fully_activated_increment=0.7,
        )
        model = IndependentWholeCell(
            self.parameters,
            ae4_parameters=self.ae4,
            regulation=regulation,
            stimulus=SecretagogueProtocol(),
        )
        state = np.concatenate((self.reference, np.asarray(regulation.basal_state)))
        grid = np.asarray((0.0, 1.0e-6, 0.02, 0.05, 0.10))
        radau = model.integrate(
            state,
            grid,
            method="Radau",
            max_step_s=0.01,
            preserve_basal_zero=True,
        )
        bdf = model.integrate(
            state,
            grid,
            method="BDF",
            max_step_s=0.01,
            preserve_basal_zero=True,
        )
        self.assertTrue(radau.success, radau.message)
        self.assertTrue(bdf.success, bdf.message)
        self.assertEqual(
            radau.onset_handling,
            "ANALYTIC_BASAL_ZERO_THEN_POSITIVE_RIGHT_LIMIT",
        )
        self.assertGreater(radau.capacity_multiplier[-1], radau.capacity_multiplier[0])
        np.testing.assert_allclose(radau.states[:, -1], bdf.states[:, -1], rtol=3.0e-7, atol=3.0e-8)
        self.assertLess(radau.max_abs_conservation_fmol_s, 1.0e-12)
        self.assertLess(bdf.max_abs_conservation_fmol_s, 1.0e-12)
        self.assertLess(radau.max_abs_water_accounting_pL_s, 1.0e-12)
        self.assertLess(bdf.max_abs_water_accounting_pL_s, 1.0e-12)

    def test_r0_to_r3_formulas_match_the_production_regulatory_api(self) -> None:
        from modern_full_model.camp_pka import (
            R0Static,
            R1EffectiveActivation,
            R2CampActivation,
            R3PkaRegulatedFraction,
            RegulatoryGain,
        )

        gain = RegulatoryGain(
            basal_capacity_multiplier=1.1,
            fully_activated_increment=0.7,
            coupling_scale=0.8,
        )
        pairs = (
            (
                R0Static(gain=gain),
                independent_r0_regulation(
                    basal_multiplier=1.1,
                    fully_activated_increment=0.7,
                    coupling_scale=0.8,
                ),
                (),
            ),
            (
                R1EffectiveActivation(tau_activation_s=17.0, gain=gain),
                independent_r1_regulation(
                    tau_activation_s=17.0,
                    basal_multiplier=1.1,
                    fully_activated_increment=0.7,
                    coupling_scale=0.8,
                ),
                (0.37,),
            ),
            (
                R2CampActivation(
                    tau_camp_s=9.0,
                    tau_activation_s=17.0,
                    gain=gain,
                ),
                independent_r2_regulation(
                    tau_camp_s=9.0,
                    tau_activation_s=17.0,
                    basal_multiplier=1.1,
                    fully_activated_increment=0.7,
                    coupling_scale=0.8,
                ),
                (0.61, 0.37),
            ),
            (
                R3PkaRegulatedFraction(
                    tau_pka_s=9.0,
                    forward_regulation_rate_s=0.04,
                    reverse_regulation_rate_s=0.02,
                    gain=gain,
                ),
                independent_r3_regulation(
                    tau_pka_s=9.0,
                    forward_regulation_rate_s=0.04,
                    reverse_regulation_rate_s=0.02,
                    basal_multiplier=1.1,
                    fully_activated_increment=0.7,
                    coupling_scale=0.8,
                ),
                (0.61, 0.37),
            ),
        )
        for production, independent, state in pairs:
            with self.subTest(family=production.family):
                self.assertEqual(independent.state_names, production.state_names)
                for beta in (0.0, 0.23, 1.0):
                    np.testing.assert_allclose(
                        independent.rhs(3.0, state, beta),
                        production.rhs(3.0, state, beta),
                        rtol=0.0,
                        atol=2.0e-16,
                    )
                    self.assertAlmostEqual(
                        independent.capacity_multiplier(3.0, state, beta),
                        production.evaluate(3.0, state, beta).capacity_multiplier,
                        places=14,
                    )

    def test_default_r2_and_r3_are_independently_isomorphic(self) -> None:
        r2 = independent_r2_regulation(
            tau_camp_s=10.0,
            tau_activation_s=30.0,
            basal_multiplier=1.0,
            fully_activated_increment=0.25,
        )
        r3 = independent_r3_regulation(
            tau_pka_s=10.0,
            forward_regulation_rate_s=1.0 / 30.0,
            reverse_regulation_rate_s=1.0 / 30.0,
            basal_multiplier=1.0,
            fully_activated_increment=0.25,
        )
        state = (0.42, 0.19)
        np.testing.assert_allclose(
            r2.rhs(0.0, state, 0.77),
            r3.rhs(0.0, state, 0.77),
            rtol=0.0,
            atol=2.0e-17,
        )
        self.assertEqual(
            r2.capacity_multiplier(0.0, state, 0.77),
            r3.capacity_multiplier(0.0, state, 0.77),
        )


class TestFrozenIndependentR1Dynamics(unittest.TestCase):
    def test_selected_root_r1_reference_600_seconds_is_independently_reproduced(self) -> None:
        repository = Path(__file__).resolve().parents[1]
        selected_path = (
            repository
            / "results"
            / "13B_modern_full_model"
            / "wt_selected_candidate.json"
        )
        selected = json.loads(selected_path.read_text(encoding="utf-8"))
        construction = selected["construction"]
        variant = construction["variant"]
        calibration = construction["calibration"]
        definition = IndependentCalibrationDefinition(
            mixed_bath_na_attempt_fraction=variant[
                "mixed_bath_na_attempt_fraction"
            ],
            apical_pump_fraction=variant["apical_pump_fraction"],
            apical_k_fraction=variant["apical_k_fraction"],
            ae4_carrier_amount_fmol_at_unit_rate_gauge=calibration[
                "ae4_carrier_amount_fmol_at_unit_rate_gauge"
            ],
            nhe1_capacity_fmol_s=calibration["nhe1_capacity_fmol_s"],
            common_membrane_conductance_scale=calibration[
                "common_membrane_conductance_scale"
            ],
            cell_other_impermeant_osmoles_fmol=calibration[
                "cell_other_impermeant_osmoles_fmol"
            ],
        )
        parameters, ae4 = definition.reconstruct_parameters()
        regulation = independent_r1_regulation(
            tau_activation_s=30.0,
            basal_multiplier=1.0,
            fully_activated_increment=0.25,
        )
        model = IndependentWholeCell(
            parameters,
            ae4_parameters=ae4,
            regulation=regulation,
            stimulus=IndependentSecretagogueProtocol(),
        )
        initial = np.concatenate(
            (
                np.asarray(construction["root_state"], dtype=float),
                np.asarray(regulation.basal_state),
            )
        )
        grid = independent_physical_time_grid(600.0, 5.0)
        atol = np.asarray(
            [
                1.0e-10
                if name.endswith("_fmol")
                else 1.0e-12
                if name.endswith("_pL")
                else 1.0e-10
                for name in model.state_names
            ]
        )
        radau = model.integrate(
            initial,
            grid,
            method="Radau",
            rtol=1.0e-7,
            atol=atol,
            max_step_s=2.0,
            preserve_basal_zero=True,
        )
        bdf = model.integrate(
            initial,
            grid,
            method="BDF",
            rtol=1.0e-7,
            atol=atol,
            max_step_s=2.0,
            preserve_basal_zero=True,
        )
        self.assertTrue(radau.success, radau.message)
        self.assertTrue(bdf.success, bdf.message)
        self.assertTrue(np.all(radau.states[:12] > 0.0))
        self.assertTrue(np.all(bdf.states[:12] > 0.0))
        self.assertAlmostEqual(radau.cumulative_flow_pL[-1], 0.49943064273, places=8)
        self.assertAlmostEqual(radau.cell_cl_mM[-1], 20.0380117092, places=7)
        self.assertAlmostEqual(radau.cell_ph[-1], 6.81025670564, places=8)
        self.assertAlmostEqual(radau.capacity_multiplier[-1], 1.25, places=8)
        self.assertLess(
            abs(radau.cumulative_flow_pL[-1] - bdf.cumulative_flow_pL[-1]),
            5.0e-7,
        )
        denominator = np.maximum(
            np.maximum(np.abs(radau.states), np.abs(bdf.states)),
            atol[:, None] / 1.0e-7,
        )
        self.assertLess(
            np.max(np.abs(radau.states - bdf.states) / denominator),
            1.0e-5,
        )
        self.assertLess(radau.max_abs_conservation_fmol_s, 1.0e-10)
        self.assertLess(radau.max_abs_water_accounting_pL_s, 1.0e-12)

        minute_flow = np.asarray(
            [
                np.interp(time_s, radau.time_s, radau.flow_pL_s)
                for time_s in range(60, 601, 60)
            ]
        )
        # No scalar gland multiplier can flatten a shape whose max/min ratio
        # exceeds the frozen WT graphical-envelope ratio.
        self.assertGreater(
            float(np.max(minute_flow) / np.min(minute_flow)),
            10.0 / 9.0,
        )


class TestIndependentN2Sustainment(unittest.TestCase):
    @staticmethod
    def _frozen_construction():
        repository = Path(__file__).resolve().parents[1]
        selected = json.loads(
            (
                repository
                / "results"
                / "13B_modern_full_model"
                / "wt_selected_candidate.json"
            ).read_text(encoding="utf-8")
        )
        construction = selected["construction"]
        variant = construction["variant"]
        calibration = construction["calibration"]
        definition = IndependentCalibrationDefinition(
            mixed_bath_na_attempt_fraction=variant[
                "mixed_bath_na_attempt_fraction"
            ],
            apical_pump_fraction=variant["apical_pump_fraction"],
            apical_k_fraction=variant["apical_k_fraction"],
            ae4_carrier_amount_fmol_at_unit_rate_gauge=calibration[
                "ae4_carrier_amount_fmol_at_unit_rate_gauge"
            ],
            nhe1_capacity_fmol_s=calibration["nhe1_capacity_fmol_s"],
            common_membrane_conductance_scale=calibration[
                "common_membrane_conductance_scale"
            ],
            cell_other_impermeant_osmoles_fmol=calibration[
                "cell_other_impermeant_osmoles_fmol"
            ],
        )
        return construction, definition.reconstruct_parameters()

    def test_n2_formula_is_basal_nested_and_changes_only_nkcc_capacity(self) -> None:
        n2 = IndependentN2NkccRegulation(
            fully_activated_multiplier=6.0,
            tau_activation_s=60.0,
        )
        self.assertEqual(n2.normalized_calcium_arm(0.058), 0.0)
        self.assertEqual(n2.normalized_calcium_arm(0.55), 1.0)
        self.assertEqual(n2.rhs(0.0, (0.0,), 0.058), (0.0,))
        self.assertEqual(n2.capacity_multiplier((0.0,)), 1.0)
        self.assertEqual(n2.capacity_multiplier((1.0,)), 6.0)
        self.assertAlmostEqual(n2.rhs(1.0, (0.25,), 0.55)[0], 0.75 / 60.0)

    def test_n1_formula_is_an_algebraic_nested_capacity_map(self) -> None:
        n1 = IndependentN1NkccRegulation(
            fully_activated_multiplier=2.0,
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.12,
        )
        self.assertEqual(n1.rhs(0.0, (), 0.058), ())
        self.assertEqual(n1.capacity_multiplier((), 0.058), 1.0)
        self.assertEqual(n1.capacity_multiplier((), 0.12), 2.0)
        self.assertAlmostEqual(
            n1.capacity_multiplier((), 0.089),
            1.5,
        )

    def test_two_frozen_n2_profiles_are_independently_reproduced(self) -> None:
        construction, (parameters, ae4) = self._frozen_construction()
        ae4_regulation = independent_r1_regulation(
            tau_activation_s=30.0,
            basal_multiplier=1.0,
            fully_activated_increment=0.25,
        )
        protocol = IndependentSecretagogueProtocol(
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.55,
        )
        grid = independent_physical_time_grid(600.0, 5.0)
        cases = {
            "N2_M2_T120": {
                "multiplier": 2.0,
                "tau_s": 120.0,
                "cumulative_pL": 0.575729623588,
                "shape_ratio": 1.27546118817,
                "endpoint_mM": (25.7254911, 94.5848860, 31.9555554),
                "shape_pass": False,
            },
            "N2_M6_T60": {
                "multiplier": 6.0,
                "tau_s": 60.0,
                "cumulative_pL": 0.671013481820,
                "shape_ratio": 1.09918107077,
                "endpoint_mM": (35.1105033, 89.5702859, 46.6044126),
                "shape_pass": True,
            },
        }
        for member_id, case in cases.items():
            with self.subTest(member_id=member_id):
                nkcc_regulation = IndependentN2NkccRegulation(
                    fully_activated_multiplier=case["multiplier"],
                    tau_activation_s=case["tau_s"],
                )
                model = IndependentWholeCell(
                    parameters,
                    ae4_parameters=ae4,
                    regulation=ae4_regulation,
                    nkcc_regulation=nkcc_regulation,
                    stimulus=protocol,
                )
                initial = np.concatenate(
                    (
                        np.asarray(construction["root_state"], dtype=float),
                        np.asarray(ae4_regulation.basal_state),
                        np.asarray(nkcc_regulation.basal_state),
                    )
                )
                atol = np.asarray(
                    [
                        1.0e-10
                        if name.endswith("_fmol")
                        else 1.0e-12
                        if name.endswith("_pL")
                        else 1.0e-10
                        for name in model.state_names
                    ]
                )
                trajectories = {
                    method: model.integrate(
                        initial,
                        grid,
                        method=method,
                        rtol=1.0e-7,
                        atol=atol,
                        max_step_s=2.0,
                        preserve_basal_zero=True,
                    )
                    for method in ("Radau", "BDF")
                }
                radau = trajectories["Radau"]
                bdf = trajectories["BDF"]
                self.assertTrue(radau.success, radau.message)
                self.assertTrue(bdf.success, bdf.message)
                self.assertAlmostEqual(
                    radau.cumulative_flow_pL[-1], case["cumulative_pL"], places=10
                )
                minute_flow = np.asarray(
                    [
                        np.interp(time_s, radau.time_s, radau.flow_pL_s)
                        for time_s in range(60, 601, 60)
                    ]
                )
                ratio = float(np.max(minute_flow) / np.min(minute_flow))
                self.assertAlmostEqual(ratio, case["shape_ratio"], places=9)
                endpoint = radau.states[:12, -1]
                endpoint_mM = endpoint[[0, 1, 2]] / endpoint[5]
                np.testing.assert_allclose(
                    endpoint_mM,
                    case["endpoint_mM"],
                    rtol=0.0,
                    atol=8.0e-8,
                )
                self.assertEqual(ratio <= 10.0 / 9.0, case["shape_pass"])
                self.assertLess(
                    abs(
                        radau.cumulative_flow_pL[-1]
                        - bdf.cumulative_flow_pL[-1]
                    ),
                    5.0e-7,
                )
                denominator = np.maximum(
                    np.maximum(np.abs(radau.states), np.abs(bdf.states)),
                    atol[:, None] / 1.0e-7,
                )
                self.assertLess(
                    np.max(np.abs(radau.states - bdf.states) / denominator),
                    1.0e-5,
                )
                self.assertLess(radau.max_abs_conservation_fmol_s, 1.0e-10)
                self.assertLess(radau.max_abs_water_accounting_pL_s, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
