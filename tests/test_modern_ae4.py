"""Focused invariants for the Task 13B AE4 QSS module."""

from __future__ import annotations

import math
import unittest

from modern_full_model.transporters import (
    AE4_BRANCH_SOURCE_STOICHIOMETRY,
    AE4BranchModulation,
    AE4Environment,
    AE4FastCooperativeGate,
    AE4Parameters,
    ae4_branch_affinity,
    ae4_carrier_rhs_s,
    ae4_reversal_cl_i_mM,
    evaluate_ae4_qss,
)


def salivary_environment() -> AE4Environment:
    """2016 thermodynamic-audit tuple, with K completing 155 mM per side."""

    return AE4Environment(
        na_i_mM=15.5,
        k_i_mM=139.5,
        cl_i_mM=50.1,
        hco3_i_mM=19.0,
        na_o_mM=151.6,
        k_o_mM=3.4,
        cl_o_mM=124.6,
        hco3_o_mM=24.7,
    )


def test_parameters(*, gate: AE4FastCooperativeGate | None = None) -> AE4Parameters:
    return AE4Parameters.equal_branch_attempts(
        carrier_amount_fmol=0.01,
        common_cl_attempt_rate_s=10.0,
        loaded_attempt_rate_s=10.0,
        cooperative_gate=gate,
    )


class SourceAndThermodynamicTests(unittest.TestCase):
    def test_each_branch_source_is_electroneutral_and_carbon_complete(self) -> None:
        for source in AE4_BRANCH_SOURCE_STOICHIOMETRY.values():
            explicit_charge = (
                source["na_i"]
                + source["k_i"]
                - source["cl_i"]
                - source["hco3_i"]
            )
            conserved_charge = (
                source["na_i"]
                + source["k_i"]
                - source["cl_i"]
                - source["alkalinity_i"]
            )
            self.assertEqual(explicit_charge, 0.0)
            self.assertEqual(conserved_charge, 0.0)
            self.assertEqual(source["tic_i"], source["hco3_i"])
            self.assertEqual(source["alkalinity_i"], source["hco3_i"])

    def test_salivary_affinities_reproduce_task13_audit(self) -> None:
        env = salivary_environment()
        self.assertAlmostEqual(ae4_branch_affinity(env, "na"), -1.894046379955519, places=12)
        self.assertAlmostEqual(ae4_branch_affinity(env, "k"), 4.100648238966857, places=12)
        self.assertAlmostEqual(ae4_reversal_cl_i_mM(env, "na"), 7.538133674727951, places=12)
        self.assertAlmostEqual(ae4_reversal_cl_i_mM(env, "k"), 3025.0087017055343, places=9)

    def test_reversal_has_zero_branch_affinity_and_current_at_global_equilibrium(self) -> None:
        equal = AE4Environment(
            na_i_mM=140.0,
            k_i_mM=5.0,
            cl_i_mM=120.0,
            hco3_i_mM=25.0,
            na_o_mM=140.0,
            k_o_mM=5.0,
            cl_o_mM=120.0,
            hco3_o_mM=25.0,
        )
        result = evaluate_ae4_qss(equal, test_parameters())
        self.assertEqual(result.affinity_na, 0.0)
        self.assertEqual(result.affinity_k, 0.0)
        self.assertAlmostEqual(result.j_na_fmol_s, 0.0, places=14)
        self.assertAlmostEqual(result.j_k_fmol_s, 0.0, places=14)
        self.assertAlmostEqual(result.source_cl_i_fmol_s, 0.0, places=14)


class QssAndCouplingTests(unittest.TestCase):
    def test_qss_closes_carrier_charge_and_entropy(self) -> None:
        result = evaluate_ae4_qss(salivary_environment(), test_parameters())
        self.assertEqual(result.units, "fmol/s")
        self.assertAlmostEqual(result.occupancy_outward + result.occupancy_inward, 1.0)
        self.assertLess(result.qss_residual_s, 1.0e-12)
        self.assertLess(result.local_detailed_balance_residual, 1.0e-12)
        self.assertAlmostEqual(result.transported_charge_equivalents_fmol_s, 0.0)
        self.assertAlmostEqual(result.transported_current_A, 0.0)
        self.assertGreaterEqual(result.entropy_production_over_r_fmol_s, -1.0e-14)
        self.assertEqual(
            set(result.intracellular_conserved_sources_fmol_s),
            {"na_i", "k_i", "cl_i", "tic_i", "alkalinity_i"},
        )

    def test_common_capacity_regulation_scales_flux_without_moving_reversal(self) -> None:
        env = salivary_environment()
        basal = evaluate_ae4_qss(env, test_parameters(), regulation_gain=1.0)
        active = evaluate_ae4_qss(env, test_parameters(), regulation_gain=1.7)
        self.assertAlmostEqual(active.j_na_fmol_s, 1.7 * basal.j_na_fmol_s)
        self.assertAlmostEqual(active.j_k_fmol_s, 1.7 * basal.j_k_fmol_s)
        self.assertEqual(active.affinity_na, basal.affinity_na)
        self.assertEqual(active.affinity_k, basal.affinity_k)
        self.assertEqual(active.occupancy_inward, basal.occupancy_inward)
        self.assertEqual(active.carrier_relaxation_time_s, basal.carrier_relaxation_time_s)

    def test_double_mutant_logic_collapses_na_and_retains_k(self) -> None:
        result = evaluate_ae4_qss(
            salivary_environment(),
            test_parameters(),
            branch_modulation=AE4BranchModulation.assay_variant("T448I_T756A"),
        )
        self.assertEqual(result.j_na_fmol_s, 0.0)
        self.assertNotEqual(result.j_k_fmol_s, 0.0)
        self.assertLess(result.local_detailed_balance_residual, 1.0e-12)

    def test_fast_sr5_gate_changes_conductance_not_affinity(self) -> None:
        env = salivary_environment()
        core = evaluate_ae4_qss(env, test_parameters())
        gated = evaluate_ae4_qss(
            env,
            test_parameters(gate=AE4FastCooperativeGate.from_2016_shape_sensitivity()),
        )
        self.assertEqual(core.model_id, "SR2_SHARED_112_QSS")
        self.assertEqual(gated.model_id, "SR5_FAST_GATE_112_QSS")
        self.assertEqual(gated.affinity_na, core.affinity_na)
        self.assertEqual(gated.affinity_k, core.affinity_k)
        self.assertLess(gated.local_detailed_balance_residual, 1.0e-12)
        self.assertNotEqual(gated.j_na_fmol_s, core.j_na_fmol_s)

    def test_qss_timescale_is_reported_not_assumed(self) -> None:
        result = evaluate_ae4_qss(salivary_environment(), test_parameters())
        rates = result.effective_rates
        delta = 1.0e-4
        rhs = ae4_carrier_rhs_s(result.occupancy_inward + delta, rates)
        self.assertAlmostEqual(rhs, -rates.relaxation_rate_s * delta, places=12)
        self.assertAlmostEqual(
            result.carrier_relaxation_time_s,
            1.0 / rates.relaxation_rate_s,
            places=15,
        )

    def test_mixed_bath_attempt_fraction_is_explicit(self) -> None:
        parameters = AE4Parameters(
            carrier_amount_fmol=0.01,
            common_cl_attempt_rate_s=10.0,
            na_loaded_attempt_rate_s=3.0,
            k_loaded_attempt_rate_s=7.0,
        )
        result = evaluate_ae4_qss(salivary_environment(), parameters)
        self.assertAlmostEqual(parameters.mixed_bath_na_attempt_fraction, 0.3)
        self.assertAlmostEqual(result.mixed_bath_na_attempt_fraction, 0.3)


if __name__ == "__main__":
    unittest.main()
