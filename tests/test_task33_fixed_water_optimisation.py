"""Structural tests for Task 33 that do not run its numerical optimisation."""
from dataclasses import asdict
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from modern_full_model.task33_fixed_water_optimisation import (
    AE4_FIVE_PERCENT,
    AE4_ZERO,
    CONTRACT,
    CONTRACT_SHA256,
    HARD_WT_VECTOR_CEILING,
    INITIAL_MULTIPLIERS,
    INITIAL_LOG_MULTIPLIERS,
    MULTIPLIER_LOWER,
    MULTIPLIER_UPPER,
    NHE1_CARRIER_AMOUNT_FMOL,
    NHE1_CHA_2009,
    OPTIMISER_VECTOR_CEILING,
    PARAMETER_NAMES,
    REFERENCE_VALUES,
    absolute_free_parameters,
    build_model,
    changed_parameter_paths,
    multipliers_from_log,
)
from modern_full_model.validation import sha256_file


class Task33StructuralTests(unittest.TestCase):
    def test_contract_is_the_checkpointed_immutable_file(self):
        self.assertEqual(sha256_file(CONTRACT), CONTRACT_SHA256)

    def test_declared_transform_and_initial_vector(self):
        np.testing.assert_allclose(
            multipliers_from_log(INITIAL_LOG_MULTIPLIERS),
            INITIAL_MULTIPLIERS,
            rtol=5e-16,
            atol=0.0,
        )
        self.assertTrue(np.all(INITIAL_MULTIPLIERS >= MULTIPLIER_LOWER))
        self.assertTrue(np.all(INITIAL_MULTIPLIERS <= MULTIPLIER_UPPER))
        self.assertEqual(OPTIMISER_VECTOR_CEILING, 19)
        self.assertEqual(HARD_WT_VECTOR_CEILING, 20)

    def test_only_six_predeclared_magnitudes_change(self):
        trial = np.asarray((0.8, 2.0, 1.5, 1.2, 1.1, 0.9))
        self.assertEqual(
            set(changed_parameter_paths(trial)),
            {
                "parameters.homeostasis.nkcc1_capacity_fmol_s",
                "ae4_parameters.carrier_amount_fmol",
                "parameters.homeostasis.ae2_capacity_fmol_s",
                "parameters.membranes.nak_capacity_fmol_s",
                "parameters.membranes.g_k_total_S",
                "parameters.membranes.g_cl_apical_S",
            },
        )
        absolute = absolute_free_parameters(trial)
        np.testing.assert_allclose(
            np.asarray([absolute[name] for name in PARAMETER_NAMES]), REFERENCE_VALUES * trial
        )

    def test_nhe1_water_and_allocation_fractions_are_frozen(self):
        model, background = build_model(np.asarray((0.8, 2.0, 1.5, 1.2, 1.1, 0.9)), stimulated=False)
        self.assertEqual(model.parameters.homeostasis.nhe1_model, NHE1_CHA_2009)
        self.assertEqual(
            model.parameters.homeostasis.nhe1_cha_carrier_amount_fmol,
            NHE1_CARRIER_AMOUNT_FMOL,
        )
        self.assertEqual(asdict(model.parameters.water), asdict(background.inherited_parameters.water))
        self.assertEqual(
            model.parameters.membranes.apical_pump_fraction,
            background.inherited_parameters.membranes.apical_pump_fraction,
        )
        self.assertEqual(
            model.parameters.membranes.apical_k_fraction,
            background.inherited_parameters.membranes.apical_k_fraction,
        )

    def test_only_two_held_out_ae4_levels_are_declared(self):
        self.assertEqual(AE4_FIVE_PERCENT.ae4_expression, 0.05)
        self.assertEqual(AE4_ZERO.ae4_expression, 0.0)


if __name__ == "__main__":
    unittest.main()
