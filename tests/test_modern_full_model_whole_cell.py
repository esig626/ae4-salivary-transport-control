"""Focused Task 13B invariants for the modern whole-cell architecture."""

from __future__ import annotations

from dataclasses import replace
import inspect
import unittest

import numpy as np

from modern_full_model import acid_base, membranes, model, parameters, states, water
from modern_full_model.transporters import AE4Parameters


class TestConservedCoordinates(unittest.TestCase):
    def test_fmol_per_pl_is_exactly_millimolar(self) -> None:
        compartment = states.compartment_from_concentrations(
            na_mM=12.5,
            k_mM=137.0,
            cl_mM=39.0,
            tic_mM=21.0,
            alkalinity_mM=31.0,
            volume_pL=1.7,
        )
        recovered = compartment.concentrations()
        self.assertAlmostEqual(recovered.na_mM, 12.5)
        self.assertAlmostEqual(recovered.k_mM, 137.0)
        self.assertAlmostEqual(recovered.cl_mM, 39.0)
        self.assertAlmostEqual(recovered.tic_mM, 21.0)
        self.assertAlmostEqual(recovered.alkalinity_mM, 31.0)

    def test_state_layout_has_explicit_cell_and_lumen_carbon_and_volume(self) -> None:
        required = {
            "na_i_fmol",
            "k_i_fmol",
            "cl_i_fmol",
            "tic_i_fmol",
            "alk_i_fmol",
            "volume_i_pL",
            "na_l_fmol",
            "k_l_fmol",
            "cl_l_fmol",
            "tic_l_fmol",
            "alk_l_fmol",
            "volume_l_pL",
        }
        self.assertEqual(set(states.CORE_STATE_NAMES), required)


class TestFiniteAcidBaseBlock(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = parameters.AcidBaseParameters()

    def test_tic_ta_round_trip_recovers_pH_and_explicit_species(self) -> None:
        target_ph = 7.21
        ta = acid_base.total_alkalinity_mM(
            ph=target_ph,
            total_carbon_mM=20.0,
            buffer_total_mM=30.0,
            buffer_pka=7.0,
            parameters=self.parameters,
        )
        result = acid_base.speciate(
            total_carbon_mM=20.0,
            total_alkalinity_mM_value=ta,
            buffer_total_mM=30.0,
            buffer_pka=7.0,
            parameters=self.parameters,
        )
        self.assertAlmostEqual(result.ph, target_ph, places=10)
        self.assertAlmostEqual(result.carbon_residual_mM, 0.0, places=12)
        self.assertAlmostEqual(result.alkalinity_residual_mM, 0.0, places=11)
        self.assertGreater(result.hco3_mM, 0.0)
        self.assertGreater(result.co2_mM, 0.0)
        self.assertGreater(result.buffer_h_mM, 0.0)
        self.assertGreater(result.buffer_minus_mM, 0.0)

    def test_closed_reactions_conserve_carbon_buffer_and_charge(self) -> None:
        residuals = acid_base.closed_reaction_invariant_residuals()
        for name, residual in residuals.items():
            np.testing.assert_array_equal(residual, np.zeros(2), err_msg=name)

    def test_external_species_source_has_correct_tic_ta_signature(self) -> None:
        source = acid_base.carbon_alkalinity_source(
            hco3_fmol_s=2.0,
            co3_fmol_s=3.0,
            co2_fmol_s=5.0,
            h_fmol_s=-7.0,
        )
        self.assertEqual(source.total_carbon_fmol_s, 10.0)
        self.assertEqual(source.total_alkalinity_fmol_s, 15.0)


class TestHomeostasisAndCationTopology(unittest.TestCase):
    def test_nkcc_nhe_ae2_combination_is_electroneutral(self) -> None:
        p = parameters.FullModelParameters()
        flux = membranes.evaluate_homeostasis(
            membranes.HomeostasisEnvironment(
                na_i_mM=15.0,
                k_i_mM=140.0,
                cl_i_mM=40.0,
                h_i_mM=1.0e-4,
                hco3_i_mM=18.0,
                na_e_mM=145.0,
                k_e_mM=5.0,
                cl_e_mM=126.0,
                h_e_mM=4.0e-5,
                hco3_e_mM=24.0,
            ),
            p,
        )
        self.assertLess(abs(flux.charge_source_residual_fmol_s), 1.0e-14)

    def test_each_reversible_law_vanishes_at_its_thermodynamic_reversal(self) -> None:
        p = parameters.FullModelParameters()
        environment = membranes.HomeostasisEnvironment(
            na_i_mM=10.0,
            k_i_mM=5.0,
            cl_i_mM=20.0,
            h_i_mM=1.0e-4,
            hco3_i_mM=10.0,
            na_e_mM=10.0,
            k_e_mM=5.0,
            cl_e_mM=20.0,
            h_e_mM=1.0e-4,
            hco3_e_mM=10.0,
        )
        flux = membranes.evaluate_homeostasis(environment, p)
        self.assertAlmostEqual(flux.nkcc1_inward_fmol_s, 0.0)
        self.assertAlmostEqual(flux.nhe1_inward_fmol_s, 0.0)
        self.assertAlmostEqual(flux.ae2_inward_fmol_s, 0.0)

    def test_new_apical_cation_directions_have_rank_two_and_no_acid_base_row(self) -> None:
        matrix = membranes.cation_topology_source_directions()
        self.assertEqual(np.linalg.matrix_rank(matrix), 2)
        self.assertEqual(matrix.shape, (4, 2))

    def test_zero_apical_fractions_remove_apical_pump_and_k_path(self) -> None:
        base = parameters.FullModelParameters()
        membrane_parameters = replace(
            base.membranes,
            apical_pump_fraction=0.0,
            apical_k_fraction=0.0,
            g_apical_background_S=0.0,
        )
        p = replace(base, membranes=membrane_parameters)
        full = model.ModernFullModel(parameters=p)
        diagnostics = full.evaluate(0.0, full.initial_state()).diagnostics.membranes
        self.assertEqual(diagnostics.pump_apical_fmol_s, 0.0)
        self.assertEqual(diagnostics.currents_A["k_apical"], 0.0)
        self.assertEqual(diagnostics.capacity_partition_residuals["pump_fraction"], 0.0)
        self.assertEqual(diagnostics.capacity_partition_residuals["k_fraction"], 0.0)

    def test_apical_pump_uses_local_luminal_k(self) -> None:
        full = model.ModernFullModel()
        low = full.initial_state()
        high = low.copy()
        high[7] *= 3.0
        low_flux = full.evaluate(0.0, low).diagnostics.membranes.pump_apical_fmol_s
        high_flux = full.evaluate(0.0, high).diagnostics.membranes.pump_apical_fmol_s
        self.assertGreater(high_flux, low_flux)


class TestElectricalWaterAndIntegratedBalance(unittest.TestCase):
    def setUp(self) -> None:
        self.ae4_parameters = AE4Parameters.equal_branch_attempts(
            carrier_amount_fmol=1.0e-4,
            common_cl_attempt_rate_s=1.0,
            loaded_attempt_rate_s=1.0,
        )
        self.full = model.ModernFullModel(ae4_parameters=self.ae4_parameters)
        self.state = self.full.initial_state()
        self.evaluation = self.full.evaluate(0.0, self.state)

    def test_reference_state_is_positive_electroneutral_and_recovers_target_pH(self) -> None:
        self.assertTrue(np.all(self.state > 0.0))
        diagnostics = self.evaluation.diagnostics
        self.assertLess(abs(diagnostics.state_charge_fmol["cell"]), 1.0e-12)
        self.assertLess(abs(diagnostics.state_charge_fmol["lumen"]), 1.0e-12)
        self.assertAlmostEqual(diagnostics.observables.cell_acid_base.ph, 7.2, places=10)
        self.assertAlmostEqual(diagnostics.observables.lumen_acid_base.ph, 7.4, places=10)

    def test_both_membrane_currents_close(self) -> None:
        residuals = self.evaluation.diagnostics.membranes.current_residuals_A
        self.assertLess(abs(residuals["apical"]), 1.0e-24)
        self.assertLess(abs(residuals["basolateral"]), 1.0e-24)
        charge_residuals = self.evaluation.diagnostics.membranes.charge_residuals_fmol_s
        self.assertLess(abs(charge_residuals["cell_source_minus_current"]), 1.0e-12)
        self.assertLess(abs(charge_residuals["lumen_source_minus_current"]), 1.0e-12)

    def test_integrated_carbon_charge_buffer_and_water_accounting(self) -> None:
        residuals = self.evaluation.diagnostics.conservation_residuals
        for key in (
            "cell_bulk_charge_rate_fmol_s",
            "lumen_bulk_charge_rate_minus_outflow_fmol_s",
            "carbon_accounting_fmol_s",
            "buffer_site_accounting_fmol_s",
            "water_volume_accounting_pL_s",
            "homeostasis_charge_fmol_s",
            "ae4_charge_fmol_s",
        ):
            self.assertLess(abs(residuals[key]), 1.0e-12, key)

    def test_ae4_deletion_is_an_exact_source_zero_without_other_parameter_change(self) -> None:
        wt = self.evaluation.diagnostics
        ko = self.full.evaluate(0.0, self.state, genotype=model.AE4_NULL).diagnostics
        self.assertNotEqual(wt.ae4.cl_cell_fmol_s, 0.0)
        self.assertEqual(ko.ae4.na_cell_fmol_s, 0.0)
        self.assertEqual(ko.ae4.k_cell_fmol_s, 0.0)
        self.assertEqual(ko.ae4.cl_cell_fmol_s, 0.0)
        self.assertEqual(ko.ae4.hco3_cell_fmol_s, 0.0)
        self.assertEqual(wt.membranes, ko.membranes)
        self.assertEqual(wt.homeostasis, ko.homeostasis)

    def test_short_radau_and_bdf_integrations_are_positive_and_agree(self) -> None:
        terminal = []
        for method in ("Radau", "BDF"):
            result = self.full.solve_dynamics(
                (0.0, 0.1),
                initial_state=self.state,
                method=method,
                t_eval=(0.0, 0.1),
                max_step_s=0.02,
            )
            self.assertTrue(result.success, result.message)
            self.assertTrue(np.all(result.y[: len(states.CORE_STATE_NAMES)] > 0.0))
            terminal.append(result.y[:, -1])
        np.testing.assert_allclose(terminal[0], terminal[1], rtol=2.0e-7, atol=2.0e-8)

    def test_root_attempt_is_bounded_and_does_not_claim_a_remote_solution(self) -> None:
        root = self.full.solve_resting_root(start_count=1, max_nfev=20)
        self.assertTrue(np.all(np.isfinite(root.state)))
        self.assertTrue(np.all(root.state > 0.0))
        np.testing.assert_array_less(root.state, self.state * 10.0 + 1.0e-12)
        self.assertTrue(np.isfinite(root.max_abs_scaled_rhs))
        # The executable starting parameters are not a WT calibration; a
        # failed bounded root attempt must be returned honestly.
        self.assertFalse(root.success)
        self.assertFalse(root.production_eligible)
        self.assertIn("charge manifold monitored but not parameterized", root.limitations)

    def test_no_heldout_secretion_target_is_present_in_implementation(self) -> None:
        source = "\n".join(
            inspect.getsource(module)
            for module in (acid_base, membranes, model, parameters, states, water)
        )
        self.assertNotIn("0.65", source)
        self.assertNotIn("35%", source)
        self.assertNotIn("35 +/-", source)


class TestParameterProvenance(unittest.TestCase):
    def test_every_parameter_has_a_required_provenance_category(self) -> None:
        required = {item.value for item in parameters.Provenance}
        records = parameters.parameter_records(parameters.FullModelParameters())
        self.assertTrue(records)
        for record in records:
            self.assertIn(record.provenance, required, record.name)


class TestPhysicalParameterDomains(unittest.TestCase):
    def setUp(self) -> None:
        self.base = parameters.FullModelParameters()

    def test_negative_or_nonpositive_physical_parameters_are_rejected_early(self) -> None:
        constructors = {
            "zero temperature": lambda: replace(
                self.base,
                constants=replace(self.base.constants, temperature_K=0.0),
            ),
            "zero Faraday constant": lambda: replace(
                self.base,
                constants=replace(self.base.constants, faraday_C_mol=0.0),
            ),
            "negative conductance": lambda: replace(
                self.base,
                membranes=replace(self.base.membranes, g_para_na_S=-1.0e-12),
            ),
            "negative hydraulic coefficient": lambda: replace(
                self.base,
                water=replace(
                    self.base.water, apical_hydraulic_pL_s_mOsm=-1.0e-5
                ),
            ),
            "negative outflow": lambda: replace(
                self.base,
                water=replace(self.base.water, outflow_rate_s=-0.1),
            ),
            "negative transporter capacity": lambda: replace(
                self.base,
                homeostasis=replace(
                    self.base.homeostasis, nkcc1_capacity_fmol_s=-0.1
                ),
            ),
            "negative pump capacity": lambda: replace(
                self.base,
                membranes=replace(self.base.membranes, nak_capacity_fmol_s=-0.1),
            ),
            "negative buffer amount": lambda: replace(
                self.base,
                geometry=replace(self.base.geometry, cell_buffer_total_fmol=-1.0),
            ),
            "negative impermeant amount": lambda: replace(
                self.base,
                geometry=replace(
                    self.base.geometry, cell_impermeant_osmoles_fmol=-1.0
                ),
            ),
            "zero half saturation": lambda: replace(
                self.base,
                membranes=replace(self.base.membranes, nak_na_half_mM=0.0),
            ),
        }
        for label, constructor in constructors.items():
            with self.subTest(label=label), self.assertRaises(ValueError):
                constructor()

    def test_zero_is_allowed_for_optional_paths_and_capacities(self) -> None:
        candidate = replace(
            self.base,
            geometry=replace(
                self.base.geometry,
                cell_buffer_total_fmol=0.0,
                lumen_buffer_total_fmol=0.0,
                fixed_cell_anion_equivalents_fmol=0.0,
                cell_impermeant_osmoles_fmol=0.0,
            ),
            homeostasis=replace(
                self.base.homeostasis,
                nkcc1_capacity_fmol_s=0.0,
                nhe1_capacity_fmol_s=0.0,
                ae2_capacity_fmol_s=0.0,
                co2_basolateral_permeability_fmol_s_mM=0.0,
                co2_apical_permeability_fmol_s_mM=0.0,
            ),
            membranes=replace(
                self.base.membranes,
                nak_capacity_fmol_s=0.0,
                g_k_total_S=0.0,
                g_cl_apical_S=0.0,
                g_para_na_S=0.0,
                g_para_k_S=0.0,
                g_para_cl_S=0.0,
                g_para_hco3_S=0.0,
                g_apical_background_S=0.0,
                g_basolateral_background_S=0.0,
            ),
            water=replace(
                self.base.water,
                apical_hydraulic_pL_s_mOsm=0.0,
                basolateral_hydraulic_pL_s_mOsm=0.0,
                paracellular_hydraulic_pL_s_mOsm=0.0,
                outflow_rate_s=0.0,
            ),
        )
        self.assertIsInstance(candidate, parameters.FullModelParameters)

    def test_finite_buffer_is_counted_once_as_an_osmotic_particle_pool(self) -> None:
        geometry = self.base.geometry
        total = water.cell_impermeant_osmoles_fmol(
            other_impermeant_fmol=geometry.cell_impermeant_osmoles_fmol,
            finite_buffer_fmol=geometry.cell_buffer_total_fmol,
        )
        self.assertAlmostEqual(total, 86.16159280366088)

        base_model = model.ModernFullModel(self.base)
        base_osm = base_model.evaluate(
            0.0, base_model.initial_state()
        ).diagnostics.observables.osmolarities_mOsm["cell"]
        new_buffer = geometry.cell_buffer_total_fmol + 10.0
        initial = self.base.initial
        new_ta = acid_base.total_alkalinity_mM(
            ph=initial.cell_ph,
            total_carbon_mM=initial.cell_tic_mM,
            buffer_total_mM=new_buffer / initial.cell_volume_pL,
            buffer_pka=self.base.acid_base.cell_buffer_pka,
            parameters=self.base.acid_base,
        )
        new_fixed_charge = (
            initial.cell_na_mM
            + initial.cell_k_mM
            - initial.cell_cl_mM
            - new_ta
        ) * initial.cell_volume_pL
        changed = replace(
            self.base,
            geometry=replace(
                geometry,
                cell_buffer_total_fmol=new_buffer,
                fixed_cell_anion_equivalents_fmol=new_fixed_charge,
            ),
        )
        changed_model = model.ModernFullModel(changed)
        changed_osm = changed_model.evaluate(
            0.0, changed_model.initial_state()
        ).diagnostics.observables.osmolarities_mOsm["cell"]
        self.assertAlmostEqual(changed_osm - base_osm, 10.0, places=12)

    def test_stale_derived_electroneutrality_values_are_rejected(self) -> None:
        stale_bath = replace(
            self.base,
            bath=replace(self.base.bath, na_mM=self.base.bath.na_mM + 1.0),
        )
        with self.assertRaisesRegex(ValueError, "derived electroneutrality"):
            model.ModernFullModel(stale_bath)

        stale_cell = replace(
            self.base,
            initial=replace(
                self.base.initial,
                cell_na_mM=self.base.initial.cell_na_mM + 1.0,
            ),
        )
        with self.assertRaisesRegex(ValueError, "derived electroneutrality"):
            model.ModernFullModel(stale_cell)

    def test_coordinated_nondefault_charge_tuple_is_accepted(self) -> None:
        coordinated = replace(
            self.base,
            bath=replace(
                self.base.bath,
                na_mM=self.base.bath.na_mM + 1.0,
                cl_mM=self.base.bath.cl_mM + 1.0,
            ),
            geometry=replace(
                self.base.geometry,
                fixed_cell_anion_equivalents_fmol=(
                    self.base.geometry.fixed_cell_anion_equivalents_fmol
                    + self.base.initial.cell_volume_pL
                ),
            ),
            initial=replace(
                self.base.initial,
                cell_na_mM=self.base.initial.cell_na_mM + 1.0,
            ),
        )
        full = model.ModernFullModel(coordinated)
        diagnostics = full.evaluate(0.0, full.initial_state()).diagnostics
        self.assertLess(abs(diagnostics.state_charge_fmol["cell"]), 1.0e-12)


if __name__ == "__main__":
    unittest.main()
