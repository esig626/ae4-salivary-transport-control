"""Thermodynamic and interface tests for the Task-13 AE4 CTMC cycles.

These tests deliberately exercise the frozen SR registry as graphs rather than
checking only one fitted output.  In particular, a shared-state multicycle is
audited by total entropy production; an individual branch is not incorrectly
required to run down its isolated affinity when another branch can drive it.
"""

from __future__ import annotations

from collections import defaultdict
import math
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from ae4_mechanism_reconstruction.chassis import (
    AE4Environment,
    AE4_KNOCKOUT,
    BASELINE_STATE,
    FixedChassis,
    ZeroAE4,
)
from state_resolved_ae4.adapter import (
    StateResolvedAE4Adapter,
    UnsupportedChassisChemistry,
    transporter_reservoirs_with_carbonate,
)
from state_resolved_ae4.chemistry import CarbonateChemistry
from state_resolved_ae4.cycles import (
    CANDIDATE_DEFINITIONS,
    CycleParameters,
    ReservoirActivities,
    branch_affinity_from_sources,
    evaluate_cycle,
    exact_charge_per_branch,
    parameters_for_mutant,
)
from state_resolved_ae4.evidence import (
    CAL_T_TARGET_IDS,
    CAL_WT_TARGET_IDS,
    EXPECTED_FREEZE_HASHES,
    STAGE_A_HELDOUT_IDS,
    STRICT_SECRETION_HOLDOUT_IDS,
    assert_calibration_targets_allowed,
    assert_stage_b_ionic_targets_only,
    verify_evidence_freeze,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def reservoirs_with_cl_ratio(cl_external: float) -> ReservoirActivities:
    """A condition in which Cl alone sets every registered cycle's sign."""

    return ReservoirActivities(
        na_i=1.0,
        k_i=1.0,
        cl_i=1.0,
        hco3_i=1.0,
        na_e=1.0,
        k_e=1.0,
        cl_e=cl_external,
        hco3_e=1.0,
        co3_i=1.0,
        co3_e=1.0,
    )


SALIVARY_RESERVOIRS = transporter_reservoirs_with_carbonate(
    na_i=15.5,
    k_i=139.5,
    cl_i=50.1,
    hco3_i=19.0,
    ph_i=6.91,
    na_e=151.6,
    k_e=3.4,
    cl_e=124.6,
    hco3_e=24.7,
    ph_e=7.4,
)


class RegisteredGraphThermodynamicsTests(unittest.TestCase):
    def test_required_sr_families_are_registered(self) -> None:
        required = {
            "SR1_NA_112",
            "SR1_K_112",
            "SR2_SHARED_112",
            "SR2_SHARED_123",
            "SR2_SHARED_213",
            "SR3_SEQUENTIAL_1111",
            "SR4A_CARBONATE_111",
            "SR4B_CARBONATE_1112",
            "SR5_COORDINATION_112",
            "SR6_PKA_EDGE_112",
            "SR6_PKA_COMMON_FLIP_112",
        }
        self.assertLessEqual(required, set(CANDIDATE_DEFINITIONS))

    def test_every_registered_branch_is_exactly_electroneutral(self) -> None:
        for model_id, definition in CANDIDATE_DEFINITIONS.items():
            with self.subTest(model=model_id):
                charges = exact_charge_per_branch(definition)
                self.assertEqual(set(charges), set(definition.branch_current_edges))
                for branch, charge in charges.items():
                    with self.subTest(branch=branch):
                        self.assertEqual(charge, 0.0)

                evaluated = evaluate_cycle(
                    definition, SALIVARY_RESERVOIRS, pka_activation=0.5
                )
                self.assertLess(abs(evaluated.transported_charge_rate), 2e-12)

    def test_local_detailed_balance_and_stationary_residual_for_all_graphs(self) -> None:
        for model_id, definition in CANDIDATE_DEFINITIONS.items():
            # Unequal symmetric barrier scalings must change kinetics without
            # changing an equilibrium ratio.
            parameters = CycleParameters(
                na_attempt_scale=1.7,
                k_attempt_scale=0.6,
                common_attempt_scale=1.2,
                edge_barrier_overrides={definition.edges[0].edge_id: 2.3},
            )
            with self.subTest(model=model_id):
                result = evaluate_cycle(
                    definition,
                    SALIVARY_RESERVOIRS,
                    parameters,
                    pka_activation=0.67,
                )
                self.assertLess(result.local_detailed_balance_residual, 5e-12)
                self.assertLess(result.generator_residual, 5e-12)
                self.assertAlmostEqual(sum(result.occupancy.values()), 1.0, places=13)
                self.assertGreater(min(result.occupancy.values()), 0.0)
                self.assertGreaterEqual(result.entropy_production, 0.0)

                for branch, affinity in result.branch_affinities.items():
                    independently_computed = branch_affinity_from_sources(
                        definition,
                        branch,
                        SALIVARY_RESERVOIRS,
                        parameters.reference_concentration_mM,
                    )
                    self.assertAlmostEqual(affinity, independently_computed, places=13)

    def test_affinity_current_sign_and_reversal_for_all_graphs(self) -> None:
        forward = reservoirs_with_cl_ratio(2.0)
        reverse = reservoirs_with_cl_ratio(0.5)
        equilibrium = reservoirs_with_cl_ratio(1.0)
        for model_id, definition in CANDIDATE_DEFINITIONS.items():
            with self.subTest(model=model_id, direction="forward"):
                result = evaluate_cycle(definition, forward, pka_activation=0.5)
                for branch in result.branch_currents:
                    self.assertGreater(result.branch_affinities[branch], 0.0)
                    self.assertGreater(result.branch_currents[branch], 0.0)
            with self.subTest(model=model_id, direction="reverse"):
                result = evaluate_cycle(definition, reverse, pka_activation=0.5)
                for branch in result.branch_currents:
                    self.assertLess(result.branch_affinities[branch], 0.0)
                    self.assertLess(result.branch_currents[branch], 0.0)
            with self.subTest(model=model_id, direction="reversal"):
                result = evaluate_cycle(definition, equilibrium, pka_activation=0.5)
                for branch in result.branch_currents:
                    self.assertAlmostEqual(result.branch_affinities[branch], 0.0, places=14)
                    self.assertLess(abs(result.branch_currents[branch]), 1e-12)


class SharedStateAndSequentialCycleTests(unittest.TestCase):
    def test_sr2_has_one_shared_pool_and_resolved_branch_currents(self) -> None:
        sr2_ids = sorted(
            model_id
            for model_id in CANDIDATE_DEFINITIONS
            if model_id.startswith("SR2_SHARED_")
        )
        self.assertEqual(len(sr2_ids), 3)
        for model_id in sr2_ids:
            definition = CANDIDATE_DEFINITIONS[model_id]
            result = evaluate_cycle(definition, SALIVARY_RESERVOIRS)
            with self.subTest(model=model_id):
                self.assertEqual(set(definition.branch_current_edges), {"na", "k"})
                self.assertEqual(definition.states.count("Eo"), 1)
                self.assertEqual(definition.states.count("Ei"), 1)
                self.assertEqual(len(definition.states), 8)
                self.assertAlmostEqual(sum(result.occupancy.values()), 1.0, places=14)

                for branch in ("na", "k"):
                    self.assertEqual(
                        definition.branch_current_edges[branch],
                        (f"{branch}_loaded_flip",),
                    )
                    self.assertAlmostEqual(
                        result.branch_currents[branch],
                        result.edge_currents[f"{branch}_loaded_flip"],
                        places=13,
                    )
                self.assertAlmostEqual(
                    result.edge_currents["common_cl_flip"],
                    sum(result.branch_currents.values()),
                    places=12,
                )

                source_from_branches: defaultdict[str, float] = defaultdict(float)
                for branch, current in result.branch_currents.items():
                    for species, coefficient in definition.intracellular_source_vectors[
                        branch
                    ].items():
                        source_from_branches[species] += coefficient * current
                for species, value in result.intracellular_sources.items():
                    self.assertAlmostEqual(value, source_from_branches[species], places=13)

    def test_sr3_1111_direction_charge_and_salivary_affinity(self) -> None:
        definition = CANDIDATE_DEFINITIONS["SR3_SEQUENTIAL_1111"]
        self.assertEqual(
            dict(definition.intracellular_source_vectors["sequential"]),
            {"na_i": 1.0, "k_i": -1.0, "cl_i": 1.0, "hco3_i": -1.0},
        )
        self.assertEqual(len(definition.states), 6)
        expected = math.log(
            SALIVARY_RESERVOIRS.cl_e
            * SALIVARY_RESERVOIRS.na_e
            * SALIVARY_RESERVOIRS.k_i
            * SALIVARY_RESERVOIRS.hco3_i
            / (
                SALIVARY_RESERVOIRS.cl_i
                * SALIVARY_RESERVOIRS.na_i
                * SALIVARY_RESERVOIRS.k_e
                * SALIVARY_RESERVOIRS.hco3_e
            )
        )
        result = evaluate_cycle(definition, SALIVARY_RESERVOIRS)
        self.assertAlmostEqual(result.branch_affinities["sequential"], expected, places=13)
        self.assertAlmostEqual(expected, 6.6434179527, places=10)
        self.assertGreater(result.branch_currents["sequential"], 0.0)

    def test_both_carbonate_variants_use_explicit_co3_and_exact_sources(self) -> None:
        shared = CANDIDATE_DEFINITIONS["SR4A_CARBONATE_111"]
        sequential = CANDIDATE_DEFINITIONS["SR4B_CARBONATE_1112"]
        for definition in (shared, sequential):
            with self.subTest(model=definition.model_id):
                required = set(definition.required_reservoir_species)
                self.assertTrue({"co3_i", "co3_e"} <= required)
                self.assertFalse({"hco3_i", "hco3_e"} & required)

        self.assertEqual(
            dict(shared.intracellular_source_vectors["na"]),
            {"na_i": -1.0, "k_i": 0.0, "cl_i": 1.0, "co3_i": -1.0},
        )
        self.assertEqual(
            dict(shared.intracellular_source_vectors["k"]),
            {"na_i": 0.0, "k_i": -1.0, "cl_i": 1.0, "co3_i": -1.0},
        )
        self.assertEqual(
            dict(sequential.intracellular_source_vectors["co3"]),
            {"na_i": 1.0, "k_i": -2.0, "cl_i": 1.0, "co3_i": -1.0},
        )

        shared_result = evaluate_cycle(shared, SALIVARY_RESERVOIRS)
        expected_na = math.log(
            SALIVARY_RESERVOIRS.cl_e
            * SALIVARY_RESERVOIRS.na_i
            * SALIVARY_RESERVOIRS.co3_i
            / (
                SALIVARY_RESERVOIRS.cl_i
                * SALIVARY_RESERVOIRS.na_e
                * SALIVARY_RESERVOIRS.co3_e
            )
        )
        expected_k = math.log(
            SALIVARY_RESERVOIRS.cl_e
            * SALIVARY_RESERVOIRS.k_i
            * SALIVARY_RESERVOIRS.co3_i
            / (
                SALIVARY_RESERVOIRS.cl_i
                * SALIVARY_RESERVOIRS.k_e
                * SALIVARY_RESERVOIRS.co3_e
            )
        )
        self.assertAlmostEqual(shared_result.branch_affinities["na"], expected_na)
        self.assertAlmostEqual(shared_result.branch_affinities["k"], expected_k)

        sequential_result = evaluate_cycle(sequential, SALIVARY_RESERVOIRS)
        expected_sequential = math.log(
            SALIVARY_RESERVOIRS.cl_e
            * SALIVARY_RESERVOIRS.na_e
            * SALIVARY_RESERVOIRS.k_i**2
            * SALIVARY_RESERVOIRS.co3_i
            / (
                SALIVARY_RESERVOIRS.cl_i
                * SALIVARY_RESERVOIRS.na_i
                * SALIVARY_RESERVOIRS.k_e**2
                * SALIVARY_RESERVOIRS.co3_e
            )
        )
        self.assertAlmostEqual(
            sequential_result.branch_affinities["co3"], expected_sequential
        )


class CoordinationAndRegulationGraphTests(unittest.TestCase):
    def test_sr5_has_explicit_catalytic_states_and_reservoir_cancels(self) -> None:
        core = CANDIDATE_DEFINITIONS["SR2_SHARED_112"]
        definition = CANDIDATE_DEFINITIONS["SR5_COORDINATION_112"]
        self.assertEqual(len(definition.states), len(core.states) + 4)
        self.assertEqual(set(definition.catalytic_bound_states), {"na", "k"})
        self.assertEqual(
            definition.intracellular_source_vectors,
            core.intracellular_source_vectors,
        )

        edge_by_id = {edge.edge_id: edge for edge in definition.edges}
        for branch in ("na", "k"):
            with self.subTest(branch=branch):
                catalytic_states = definition.catalytic_bound_states[branch]
                self.assertEqual(len(catalytic_states), 2)
                self.assertTrue(set(catalytic_states) <= set(definition.states))
                path = (
                    edge_by_id[f"{branch}_catalyst_bind"],
                    edge_by_id[f"{branch}_catalytic_flip"],
                    edge_by_id[f"{branch}_catalyst_release"],
                )
                reservoir_sum: defaultdict[str, float] = defaultdict(float)
                for edge in path:
                    for species, coefficient in edge.reservoir_stoichiometry.items():
                        reservoir_sum[species] += coefficient
                self.assertEqual(
                    {name: value for name, value in reservoir_sum.items() if value != 0.0},
                    {},
                )

        result = evaluate_cycle(definition, SALIVARY_RESERVOIRS)
        for branch in ("na", "k"):
            expected_current = (
                result.edge_currents[f"{branch}_loaded_flip"]
                + result.edge_currents[f"{branch}_catalytic_flip"]
            )
            self.assertAlmostEqual(result.branch_currents[branch], expected_current)
            self.assertGreater(
                sum(result.occupancy[state] for state in definition.catalytic_bound_states[branch]),
                0.0,
            )
            self.assertAlmostEqual(
                result.branch_affinities[branch],
                branch_affinity_from_sources(core, branch, SALIVARY_RESERVOIRS),
            )

    def test_sr6_has_explicit_u_p_layers_and_declared_pka_placement(self) -> None:
        core = CANDIDATE_DEFINITIONS["SR5_COORDINATION_112"]
        expected_edge_placement = {
            "P::na_loaded_flip",
            "P::na_catalytic_flip",
            "P::k_loaded_flip",
            "P::k_catalytic_flip",
        }
        expected_common_placement = {"P::common_cl_flip"}
        for model_id, expected_placement in (
            ("SR6_PKA_EDGE_112", expected_edge_placement),
            ("SR6_PKA_COMMON_FLIP_112", expected_common_placement),
        ):
            definition = CANDIDATE_DEFINITIONS[model_id]
            with self.subTest(model=model_id):
                self.assertEqual(len(definition.states), 2 * len(core.states))
                self.assertEqual(
                    set(definition.states),
                    {
                        f"{state}^{layer}"
                        for state in core.states
                        for layer in ("U", "P")
                    },
                )
                self.assertEqual(set(definition.pka_switch_edge_ids), {"pka_switch"})
                self.assertEqual(set(definition.pka_edge_ids), expected_placement)
                self.assertTrue(all(edge_id.startswith("P::") for edge_id in definition.pka_edge_ids))
                switch = next(edge for edge in definition.edges if edge.edge_id == "pka_switch")
                self.assertEqual((switch.source, switch.target), ("Ei^U", "Ei^P"))
                self.assertEqual(dict(switch.reservoir_stoichiometry), {})

                low = evaluate_cycle(
                    definition,
                    reservoirs_with_cl_ratio(2.0),
                    CycleParameters(pka_basal_activity=1e-6),
                    pka_activation=0.0,
                )
                high = evaluate_cycle(
                    definition,
                    reservoirs_with_cl_ratio(2.0),
                    CycleParameters(pka_basal_activity=1e-6),
                    pka_activation=1.0,
                )
                low_p = sum(value for state, value in low.occupancy.items() if state.endswith("^P"))
                high_p = sum(value for state, value in high.occupancy.items() if state.endswith("^P"))
                self.assertLess(low_p, 1e-4)
                self.assertGreater(high_p, 0.45)
                self.assertGreater(high.cycle_flux, low.cycle_flux)
                self.assertEqual(dict(high.branch_affinities), dict(low.branch_affinities))

    def test_double_mutant_collapses_na_but_retains_k_branch(self) -> None:
        almost_absent = 1e-12
        na_assay = ReservoirActivities(
            na_i=50.0,
            k_i=almost_absent,
            cl_i=4.0,
            hco3_i=25.0,
            na_e=50.0,
            k_e=almost_absent,
            cl_e=128.3,
            hco3_e=25.0,
        )
        k_assay = ReservoirActivities(
            na_i=almost_absent,
            k_i=50.0,
            cl_i=4.0,
            hco3_i=25.0,
            na_e=almost_absent,
            k_e=50.0,
            cl_e=128.3,
            hco3_e=25.0,
        )
        mutant = parameters_for_mutant(CycleParameters(), "T448I_T756A")
        self.assertEqual(mutant.mutant_na_scale, 1e-12)
        self.assertEqual(mutant.mutant_k_scale, 0.5)
        for model_id in (
            "SR2_SHARED_112",
            "SR5_COORDINATION_112",
            "SR6_PKA_EDGE_112",
            "SR6_PKA_COMMON_FLIP_112",
        ):
            definition = CANDIDATE_DEFINITIONS[model_id]
            with self.subTest(model=model_id):
                wt_na = evaluate_cycle(
                    definition, na_assay, pka_activation=0.5
                ).branch_currents["na"]
                mutant_na = evaluate_cycle(
                    definition, na_assay, mutant, pka_activation=0.5
                ).branch_currents["na"]
                wt_k = evaluate_cycle(
                    definition, k_assay, pka_activation=0.5
                ).branch_currents["k"]
                mutant_k = evaluate_cycle(
                    definition, k_assay, mutant, pka_activation=0.5
                ).branch_currents["k"]
                single_mutant_k = evaluate_cycle(
                    definition,
                    k_assay,
                    parameters_for_mutant(CycleParameters(), "T448I"),
                    pka_activation=0.5,
                ).branch_currents["k"]
                self.assertGreater(abs(wt_na), 1e-8)
                self.assertLess(abs(mutant_na / wt_na), 1e-7)
                self.assertGreater(abs(wt_k), 1e-8)
                self.assertGreater(abs(mutant_k / wt_k), 0.4)
                self.assertAlmostEqual(mutant_k, single_mutant_k, places=12)


class CarbonateChemistryAndAdapterTests(unittest.TestCase):
    def test_carbonate_chemistry_is_explicit_and_self_consistent(self) -> None:
        chemistry = CarbonateChemistry(pka2=10.33, temperature_c=25.0)
        self.assertAlmostEqual(
            chemistry.carbonate_from_bicarbonate(25.0, chemistry.pka2),
            25.0,
        )
        hco3, co3 = chemistry.speciate_total_base(25.0, chemistry.pka2)
        self.assertAlmostEqual(hco3, 12.5)
        self.assertAlmostEqual(co3, 12.5)
        self.assertAlmostEqual(hco3 + co3, 25.0)
        self.assertAlmostEqual(
            chemistry.equilibrium_residual(hco3, co3, chemistry.pka2), 0.0
        )

        reservoirs = transporter_reservoirs_with_carbonate(
            na_i=15.5,
            k_i=139.5,
            cl_i=50.1,
            hco3_i=19.0,
            ph_i=6.91,
            na_e=151.6,
            k_e=3.4,
            cl_e=124.6,
            hco3_e=24.7,
            ph_e=7.4,
            chemistry=chemistry,
        )
        self.assertEqual(reservoirs.hco3_i, 19.0)
        self.assertAlmostEqual(
            reservoirs.co3_i,
            chemistry.carbonate_from_bicarbonate(19.0, 6.91),
        )
        self.assertAlmostEqual(
            reservoirs.co3_e,
            chemistry.carbonate_from_bicarbonate(24.7, 7.4),
        )
        self.assertNotEqual(reservoirs.co3_i, reservoirs.hco3_i)
        with self.assertRaises(ValueError):
            chemistry.carbonate_from_bicarbonate(0.0, 7.4)

    def test_adapter_knockout_short_circuit_is_candidate_independent(self) -> None:
        reference = FixedChassis(ZeroAE4()).raw_rhs(
            0.0, BASELINE_STATE, scenario=AE4_KNOCKOUT
        )
        for model_id, definition in CANDIDATE_DEFINITIONS.items():
            adapter = StateResolvedAE4Adapter(definition, CycleParameters())
            with self.subTest(model=model_id):
                rhs = FixedChassis(adapter).raw_rhs(
                    0.0, BASELINE_STATE, scenario=AE4_KNOCKOUT
                )
                np.testing.assert_array_equal(rhs, reference)
                environment = FixedChassis()._ae4_environment(BASELINE_STATE)
                contribution = adapter.evaluate(environment, activity_scale=0.0)
                self.assertEqual(
                    (
                        contribution.cycle_flux,
                        contribution.na_i,
                        contribution.k_i,
                        contribution.cl_i,
                        contribution.hco3_i,
                        contribution.transported_charge,
                    ),
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                )
                self.assertTrue(contribution.diagnostics["knockout_short_circuit"])

    def test_bicarbonate_only_chassis_rejects_both_carbonate_adapters(self) -> None:
        environment = AE4Environment(
            na_i=15.5,
            k_i=139.5,
            cl_i=50.1,
            hco3_i=19.0,
            na_e=151.6,
            k_e=3.4,
            cl_e=124.6,
            hco3_e=24.7,
            pka_activation=0.5,
        )
        for model_id in ("SR4A_CARBONATE_111", "SR4B_CARBONATE_1112"):
            adapter = StateResolvedAE4Adapter(
                CANDIDATE_DEFINITIONS[model_id], CycleParameters()
            )
            with self.subTest(model=model_id):
                with self.assertRaises(UnsupportedChassisChemistry):
                    adapter.evaluate(environment, activity_scale=1.0)


class FrozenEvidenceBoundaryTests(unittest.TestCase):
    def test_repository_evidence_freeze_hashes_verify(self) -> None:
        boundary = verify_evidence_freeze(REPOSITORY_ROOT)
        self.assertEqual(boundary.status, "verified_before_calibration")
        self.assertEqual(
            set(boundary.allowed_target_ids), CAL_T_TARGET_IDS | CAL_WT_TARGET_IDS
        )
        self.assertEqual(set(boundary.forbidden_target_ids), STAGE_A_HELDOUT_IDS)
        self.assertTrue(STRICT_SECRETION_HOLDOUT_IDS <= STAGE_A_HELDOUT_IDS)

    def test_freeze_verification_fails_closed_after_document_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for relative in EXPECTED_FREEZE_HASHES:
                destination = root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(REPOSITORY_ROOT / relative, destination)
            verify_evidence_freeze(root)
            evidence = root / "analysis/13_state_resolved_ae4/evidence_freeze.md"
            evidence.write_text(
                evidence.read_text(encoding="utf-8") + "\nchanged after freeze\n",
                encoding="utf-8",
            )
            with self.assertRaises(RuntimeError):
                verify_evidence_freeze(root)

    def test_calibration_and_stage_b_guards_reject_holdout_leakage(self) -> None:
        allowed = ("E16-07", "E21-03", "E15-05")
        self.assertEqual(assert_calibration_targets_allowed(allowed), allowed)
        for target in sorted(STAGE_A_HELDOUT_IDS | {"UNKNOWN"}):
            with self.subTest(stage="calibration", target=target):
                with self.assertRaises(PermissionError):
                    assert_calibration_targets_allowed((target,))

        stage_b = ("E15-06", "E15-07_KO")
        self.assertEqual(assert_stage_b_ionic_targets_only(stage_b), stage_b)
        for target in ("E15-02", "E15-03", "E15-08_KO", "E15-16"):
            with self.subTest(stage="stage_b", target=target):
                with self.assertRaises(PermissionError):
                    assert_stage_b_ionic_targets_only((target,))


if __name__ == "__main__":
    unittest.main()
