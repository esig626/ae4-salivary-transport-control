"""Focused tests for the single Task 30 NHE1 repair law."""
from dataclasses import replace
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from modern_full_model.membranes import (
    HomeostasisEnvironment,
    NHE1_LEGACY_TANH,
    NHE1_VERA_SIGUENZA_2018,
    evaluate_homeostasis,
    evaluate_nhe1_legacy_tanh,
    evaluate_nhe1_vera_siguenza_2018,
)
from modern_full_model.parameters import FullModelParameters, HomeostasisParameters
from modern_full_model.model import AE4_NULL, WT
from modern_full_model.task30_nhe1_repair import (
    build_task30_model,
    load_background,
    state_from_coordinates,
)


def environment(*, na_i: float, h_i: float, na_e: float, h_e: float):
    return HomeostasisEnvironment(
        na_i_mM=na_i,
        k_i_mM=120.0,
        cl_i_mM=50.0,
        h_i_mM=h_i,
        hco3_i_mM=10.0,
        na_e_mM=na_e,
        k_e_mM=5.0,
        cl_e_mM=126.0,
        h_e_mM=h_e,
        hco3_e_mM=22.0,
    )


class PublishedNhe1EquationTests(unittest.TestCase):
    def setUp(self):
        base = FullModelParameters()
        self.parameters = replace(
            base,
            homeostasis=replace(
                base.homeostasis, nhe1_model=NHE1_VERA_SIGUENZA_2018
            ),
        )

    def test_three_independently_calculated_fixed_states(self):
        cases = (
            ((15.0, 1.0e-4, 145.0, 4.0e-5), 8.121153462228204e-4),
            ((8.0, 2.0e-4, 140.0, 4.2e-5), 2.5308222676024966e-3),
            ((30.0, 3.0e-5, 150.0, 5.0e-5), -9.502367424242424e-5),
        )
        for values, expected in cases:
            actual = evaluate_nhe1_vera_siguenza_2018(
                environment(
                    na_i=values[0], h_i=values[1], na_e=values[2], h_e=values[3]
                ),
                self.parameters,
            )
            self.assertTrue(
                math.isclose(actual, expected, rel_tol=0.0, abs_tol=5.0e-19),
                (actual, expected),
            )

    def test_source_stoichiometry_and_zero_charge(self):
        homeostasis = replace(
            self.parameters.homeostasis,
            nkcc1_capacity_fmol_s=0.0,
            ae2_capacity_fmol_s=0.0,
        )
        parameters = replace(self.parameters, homeostasis=homeostasis)
        result = evaluate_homeostasis(
            environment(na_i=15.0, h_i=1.0e-4, na_e=145.0, h_e=4.0e-5),
            parameters,
        )
        flux = result.nhe1_inward_fmol_s
        self.assertEqual(result.nhe1_model, NHE1_VERA_SIGUENZA_2018)
        self.assertEqual(result.na_cell_fmol_s, flux)
        self.assertEqual(result.k_cell_fmol_s, 0.0)
        self.assertEqual(result.cl_cell_fmol_s, 0.0)
        self.assertEqual(result.tic_cell_fmol_s, 0.0)
        self.assertEqual(result.alkalinity_cell_fmol_s, flux)
        self.assertEqual(result.charge_source_residual_fmol_s, 0.0)

    def test_finite_bound_and_activity_scale(self):
        state = environment(na_i=0.01, h_i=1.0, na_e=500.0, h_e=1.0e-9)
        full = evaluate_nhe1_vera_siguenza_2018(state, self.parameters)
        half = evaluate_nhe1_vera_siguenza_2018(
            state, self.parameters, activity_scale=0.5
        )
        self.assertTrue(math.isfinite(full))
        self.assertLessEqual(abs(full), 0.0305)
        self.assertEqual(half, 0.5 * full)

    def test_acidification_increases_forward_activity_at_fixed_sodium(self):
        alkaline = environment(na_i=15.0, h_i=5.0e-5, na_e=145.0, h_e=4.0e-5)
        acidic = environment(na_i=15.0, h_i=2.0e-4, na_e=145.0, h_e=4.0e-5)
        self.assertGreater(
            evaluate_nhe1_vera_siguenza_2018(acidic, self.parameters),
            evaluate_nhe1_vera_siguenza_2018(alkaline, self.parameters),
        )

    def test_ae4_genotype_does_not_enter_nhe1_evaluator(self):
        background = load_background("R09")
        model = build_task30_model(background, stimulated=False)
        state = state_from_coordinates(model, background.saved_coordinates)
        wt = model.evaluate(0.0, state, genotype=WT)
        exact_null = model.evaluate(0.0, state, genotype=AE4_NULL)
        self.assertEqual(
            wt.diagnostics.homeostasis.nhe1_inward_fmol_s,
            exact_null.diagnostics.homeostasis.nhe1_inward_fmol_s,
        )
        self.assertEqual(
            wt.diagnostics.homeostasis.nhe1_model,
            exact_null.diagnostics.homeostasis.nhe1_model,
        )

    def test_legacy_comparator_remains_selectable(self):
        state = environment(na_i=15.0, h_i=1.0e-4, na_e=145.0, h_e=4.0e-5)
        legacy_parameters = replace(
            self.parameters,
            homeostasis=replace(
                self.parameters.homeostasis, nhe1_model=NHE1_LEGACY_TANH
            ),
        )
        result = evaluate_homeostasis(state, legacy_parameters)
        self.assertEqual(result.nhe1_model, NHE1_LEGACY_TANH)
        self.assertEqual(
            result.nhe1_inward_fmol_s,
            evaluate_nhe1_legacy_tanh(state, legacy_parameters),
        )

    def test_clean_seed_homeostasis_payload_is_backward_compatible(self):
        payload = {
            "nkcc1_capacity_fmol_s": 0.32,
            "nhe1_capacity_fmol_s": 0.0074,
            "ae2_capacity_fmol_s": 0.005,
            "thermodynamic_saturation_log_width": 2.0,
            "co2_basolateral_permeability_fmol_s_mM": 0.02,
            "co2_apical_permeability_fmol_s_mM": 0.01,
        }
        loaded = HomeostasisParameters(**payload)
        self.assertEqual(loaded.nhe1_model, NHE1_LEGACY_TANH)
        self.assertEqual(loaded.nhe1_published_g_fmol_s, 0.0305)

    def test_task30_builder_selects_only_published_nhe1_fields(self):
        for alias in ("R09", "R10"):
            background = load_background(alias)
            inherited = background.inherited_parameters
            for g_nhe1 in (0.0305, 0.123):
                model = build_task30_model(
                    background,
                    g_nhe1_fmol_s=g_nhe1,
                    stimulated=False,
                )
                active = model.parameters
                self.assertEqual(
                    active.homeostasis.nhe1_model,
                    NHE1_VERA_SIGUENZA_2018,
                )
                self.assertEqual(
                    active.homeostasis.nhe1_published_g_fmol_s,
                    g_nhe1,
                )
                self.assertEqual(active.homeostasis.nhe1_published_k_h_mM, 4.5e-4)
                self.assertEqual(active.homeostasis.nhe1_published_k_na_mM, 15.0)
                self.assertEqual(active.constants, inherited.constants)
                self.assertEqual(active.acid_base, inherited.acid_base)
                self.assertEqual(active.bath, inherited.bath)
                self.assertEqual(active.geometry, inherited.geometry)
                self.assertEqual(active.membranes, inherited.membranes)
                self.assertEqual(active.water, inherited.water)
                self.assertEqual(active.initial, inherited.initial)
                restored_homeostasis = replace(
                    active.homeostasis,
                    nhe1_model=inherited.homeostasis.nhe1_model,
                    nhe1_published_g_fmol_s=(
                        inherited.homeostasis.nhe1_published_g_fmol_s
                    ),
                    nhe1_published_k_h_mM=(
                        inherited.homeostasis.nhe1_published_k_h_mM
                    ),
                    nhe1_published_k_na_mM=(
                        inherited.homeostasis.nhe1_published_k_na_mM
                    ),
                )
                self.assertEqual(restored_homeostasis, inherited.homeostasis)


if __name__ == "__main__":
    unittest.main()
