from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ae4_mechanism_reconstruction.structure import (  # noqa: E402
    AE4Stoichiometry,
    SOURCE_STOICHIOMETRIES,
    baseline_net_k_fraction_bound,
    balance_signature,
    c1_stimulus_source_comparison,
    delta_g_over_rt,
    equilibrium_bicarbonate_i,
    equilibrium_cation_i,
    equilibrium_chloride_i,
    forward_flux_wedge,
    historical_chloride_flux_from_flow,
    historical_flow_from_chloride_flux,
    historical_flow_slope,
    fixed_state_candidate_residual,
    mixed_cation_fixed_row_correction,
    nak_nhe1_repair_coefficients,
    parallel_112_supported_chloride_flux,
    pooled_partition_error,
    recover_ae2_ae4_fluxes,
    required_mean_gate_elasticity,
    required_wt_observable_factor,
    reversible_rate,
)
from ae4_mechanism_reconstruction.chassis import (  # noqa: E402
    AE4_KNOCKOUT,
    ChassisParameters,
    FixedChassis,
)
from ae4_mechanism_reconstruction.mechanisms import (  # noqa: E402
    AE4Mechanism,
    candidate_specs,
)


class StoichiometricStructureTests(unittest.TestCase):
    def test_all_primary_paper_possibilities_are_electroneutral(self) -> None:
        for stoichiometry in SOURCE_STOICHIOMETRIES.values():
            self.assertTrue(stoichiometry.is_electroneutral)
            self.assertEqual(
                stoichiometry.bicarbonate,
                stoichiometry.chloride + stoichiometry.cation,
            )

    def test_exact_flux_inverse_for_each_stoichiometry(self) -> None:
        for stoichiometry in SOURCE_STOICHIOMETRIES.values():
            for j2, j4 in ((0.0, 0.0), (0.3, 0.7), (-0.2, 1.1)):
                aggregate_cl = j2 + stoichiometry.chloride * j4
                aggregate_b = j2 + stoichiometry.bicarbonate * j4
                recovered = recover_ae2_ae4_fluxes(
                    aggregate_cl, aggregate_b, stoichiometry
                )
                self.assertAlmostEqual(recovered[0], j2)
                self.assertAlmostEqual(recovered[1], j4)

    def test_forward_wedges_differ_by_anion_signature(self) -> None:
        one_one_two = SOURCE_STOICHIOMETRIES["1Cl_1Cat_2B"]
        one_two_three = SOURCE_STOICHIOMETRIES["1Cl_2Cat_3B"]
        two_one_three = SOURCE_STOICHIOMETRIES["2Cl_1Cat_3B"]
        self.assertTrue(forward_flux_wedge(1.0, 1.8, one_one_two))
        self.assertTrue(forward_flux_wedge(1.0, 2.8, one_two_three))
        self.assertFalse(forward_flux_wedge(1.0, 1.8, two_one_three))

    def test_nontransported_cation_breaks_one_to_one_to_two_electroneutrality(self) -> None:
        gating_only_with_published_anion_count = AE4Stoichiometry(1, 0, 2)
        gating_only_equal_anion_count = AE4Stoichiometry(1, 0, 1)
        self.assertEqual(gating_only_with_published_anion_count.charge_residual, 1)
        self.assertFalse(gating_only_with_published_anion_count.is_electroneutral)
        self.assertTrue(gating_only_equal_anion_count.is_electroneutral)

    def test_published_partition_cannot_close_na_only_historical_state(self) -> None:
        stoichiometry = SOURCE_STOICHIOMETRIES["1Cl_1Cat_2B"]
        historical = balance_signature(stoichiometry, sodium_fraction=1.0)
        published = balance_signature(stoichiometry, sodium_fraction=25.0 / 145.0)
        j4 = 0.03888333862
        residual = fixed_state_candidate_residual(
            historical_flux=j4,
            candidate_flux=j4,
            historical_signature=historical,
            candidate_signature=published,
        )
        self.assertAlmostEqual(residual[0], 0.0)
        self.assertAlmostEqual(residual[3], 0.0)
        self.assertAlmostEqual(residual[1], 0.03217931472, places=10)
        self.assertAlmostEqual(residual[2], -0.03217931472, places=10)

        # The Cl row forces candidate_flux=j4, while exact K closure with a
        # nonzero K share forces candidate_flux=0.  No kinetic scalar can meet
        # both requirements at this fixed state.
        self.assertNotEqual(published[2], 0.0)

    def test_fixed_row_k_bound_and_exact_two_path_repair(self) -> None:
        historical_flux = 0.03888333867648
        potassium_fraction = 120.0 / 145.0
        self.assertAlmostEqual(
            baseline_net_k_fraction_bound(
                required_cycle_flux=historical_flux,
                amount_tolerance=1e-8,
            ),
            2.57179561745e-7,
            places=17,
        )

        correction = mixed_cation_fixed_row_correction(
            required_cycle_flux=historical_flux,
            potassium_fraction=potassium_fraction,
        )
        pump, nhe1 = nak_nhe1_repair_coefficients(correction)
        self.assertGreater(pump, 0.0)
        self.assertAlmostEqual(pump, nhe1, places=15)

        reconstructed = (
            -3.0 * pump + nhe1,
            2.0 * pump,
            0.0,
            0.0,
        )
        for actual, expected in zip(reconstructed, correction, strict=True):
            self.assertAlmostEqual(actual, expected, places=15)

    def test_k_recruitment_rotates_only_the_stimulated_cation_signature(self) -> None:
        basal = 0.03888333867648
        common, recruited = c1_stimulus_source_comparison(
            basal_cycle_flux=basal,
            activation_fold=3.0,
        )
        self.assertEqual(common[2:], recruited[2:])
        difference = tuple(
            k_source - common_source
            for k_source, common_source in zip(recruited, common, strict=True)
        )
        self.assertAlmostEqual(difference[0], 2.0 * basal)
        self.assertAlmostEqual(difference[1], -2.0 * basal)
        self.assertEqual(difference[2:], (0.0, 0.0))

    def test_forward_k_ae4_has_a_direct_negative_chloride_support_term(self) -> None:
        without_k = parallel_112_supported_chloride_flux(
            nkcc_cycle_flux=0.2,
            pump_cycle_flux=0.1,
            k_ae4_cycle_flux=0.0,
        )
        with_k = parallel_112_supported_chloride_flux(
            nkcc_cycle_flux=0.2,
            pump_cycle_flux=0.1,
            k_ae4_cycle_flux=0.07,
        )
        self.assertAlmostEqual(without_k - with_k, 0.07)


class ThermodynamicStructureTests(unittest.TestCase):
    concentrations = {
        "chloride_i": 50.1,
        "chloride_o": 124.6,
        "bicarbonate_i": 19.0,
        "bicarbonate_o": 24.7,
    }

    def test_na_only_candidates_have_wrong_direction_at_primary_paper_state(self) -> None:
        for stoichiometry in SOURCE_STOICHIOMETRIES.values():
            affinity = delta_g_over_rt(
                stoichiometry,
                cation_i=15.5,
                cation_o=151.6,
                **self.concentrations,
            )
            self.assertGreater(affinity, 0.0)
            self.assertLess(reversible_rate(1.0, affinity), 0.0)

    def test_nonselective_one_one_two_has_inward_direction(self) -> None:
        affinity = delta_g_over_rt(
            SOURCE_STOICHIOMETRIES["1Cl_1Cat_2B"],
            cation_i=155.0,
            cation_o=155.0,
            **self.concentrations,
        )
        self.assertAlmostEqual(affinity, -0.3863590693, places=9)
        self.assertGreater(reversible_rate(1.0, affinity), 0.0)

    def test_primary_paper_thresholds_are_reproduced(self) -> None:
        expected = {
            "1Cl_1Cat_2B": (7.5381337, 103.0162151, 48.9824548),
            "1Cl_2Cat_3B": (0.5928611, 142.4866158, 83.3774323),
            "2Cl_1Cat_3B": (26.8793809, 53.8478819, 28.7761426),
        }
        for identifier, stoichiometry in SOURCE_STOICHIOMETRIES.items():
            chloride_threshold = equilibrium_chloride_i(
                stoichiometry,
                chloride_o=124.6,
                cation_i=15.5,
                cation_o=151.6,
                bicarbonate_i=19.0,
                bicarbonate_o=24.7,
            )
            cation_threshold = equilibrium_cation_i(
                stoichiometry,
                chloride_i=50.1,
                chloride_o=124.6,
                cation_o=151.6,
                bicarbonate_i=19.0,
                bicarbonate_o=24.7,
            )
            bicarbonate_threshold = equilibrium_bicarbonate_i(
                stoichiometry,
                chloride_i=50.1,
                chloride_o=124.6,
                cation_i=15.5,
                cation_o=151.6,
                bicarbonate_o=24.7,
            )
            for actual, target in zip(
                (chloride_threshold, cation_threshold, bicarbonate_threshold),
                expected[identifier],
            ):
                self.assertAlmostEqual(actual, target, places=5)

    def test_pooled_net_flux_cannot_generally_use_intracellular_partition(self) -> None:
        error = pooled_partition_error(
            forward_factor=2.0,
            reverse_factor=0.7,
            na_i=15.5,
            k_i=139.5,
            na_o=151.6,
            k_o=3.4,
        )
        self.assertNotAlmostEqual(error, 0.0)
        self.assertAlmostEqual(
            pooled_partition_error(
                forward_factor=2.0,
                reverse_factor=0.0,
                na_i=15.5,
                k_i=139.5,
                na_o=151.6,
                k_o=3.4,
            ),
            0.0,
        )


class HistoricalFlowClosureTests(unittest.TestCase):
    def test_flow_chloride_mapping_is_exact_and_monotone(self) -> None:
        gain = 0.25
        offset = 1.2
        for chloride_flux in (0.0, 0.1, 2.0, 20.0):
            flow = historical_flow_from_chloride_flux(
                chloride_flux,
                water_gain=gain,
                fixed_osmotic_offset=offset,
            )
            recovered = historical_chloride_flux_from_flow(
                flow,
                water_gain=gain,
                fixed_osmotic_offset=offset,
            )
            self.assertAlmostEqual(recovered, chloride_flux)
            self.assertGreater(
                historical_flow_slope(
                    flow,
                    water_gain=gain,
                    fixed_osmotic_offset=offset,
                ),
                0.0,
            )

    def test_frozen_chassis_baseline_satisfies_exact_flow_map(self) -> None:
        parameters = ChassisParameters()
        chassis = FixedChassis(parameters=parameters)
        _, diagnostics = chassis.evaluate(
            chassis.initial_state(), calcium=0.05
        )
        gain = (
            parameters.water_apical * parameters.water_basolateral
            / (parameters.water_apical + parameters.water_basolateral)
            + parameters.water_paracellular
        )
        offset = (
            parameters.lumen_impermeant_osmolyte
            - parameters.external_osmotic_sum
        )
        chloride_flux = -float(diagnostics["j_cl_signed"])
        predicted = historical_flow_from_chloride_flux(
            chloride_flux,
            water_gain=gain,
            fixed_osmotic_offset=offset,
        )
        self.assertAlmostEqual(predicted, float(diagnostics["q_total"]), places=14)

    def test_frozen_chassis_baseline_has_pump_k_chloride_invariant(self) -> None:
        chassis = FixedChassis()
        _, diagnostics = chassis.evaluate(
            chassis.initial_state(), calcium=0.05
        )
        chloride_out = -float(diagnostics["j_cl_signed"])
        supported_out = float(diagnostics["j_nak"]) + float(diagnostics["j_k"])
        self.assertAlmostEqual(chloride_out, supported_out, places=14)


class KnockoutInvarianceTests(unittest.TestCase):
    def test_ae4_knockout_vector_field_is_candidate_independent(self) -> None:
        state = FixedChassis.initial_state()
        reference_raw = None
        reference_flow = None
        for spec in candidate_specs(include_rejected=False):
            chassis = FixedChassis(AE4Mechanism(spec.candidate_id))
            raw, diagnostics = chassis.evaluate(
                state,
                calcium=0.55,
                scenario=AE4_KNOCKOUT,
            )
            if reference_raw is None:
                reference_raw = raw
                reference_flow = diagnostics["q_total"]
            else:
                for actual, expected in zip(raw, reference_raw, strict=True):
                    self.assertAlmostEqual(actual, expected, places=14)
                self.assertAlmostEqual(
                    diagnostics["q_total"], reference_flow, places=14
                )

    def test_source_bounded_stimulus_gate_threshold(self) -> None:
        endpoint_ratio = 0.9985561418496379
        integral_ratio = 1.0056381632041642
        target = 0.65
        self.assertAlmostEqual(
            required_wt_observable_factor(
                current_knockout_to_wt_ratio=endpoint_ratio,
                target_knockout_to_wt_ratio=target,
            ),
            1.5362402182302122,
            places=14,
        )
        self.assertAlmostEqual(
            required_wt_observable_factor(
                current_knockout_to_wt_ratio=integral_ratio,
                target_knockout_to_wt_ratio=target,
            ),
            1.547135635698714,
            places=14,
        )
        self.assertAlmostEqual(
            required_mean_gate_elasticity(
                activation_fold=3.0,
                current_knockout_to_wt_ratio=endpoint_ratio,
                target_knockout_to_wt_ratio=target,
            ),
            0.3908003023477724,
            places=14,
        )


if __name__ == "__main__":
    unittest.main()
