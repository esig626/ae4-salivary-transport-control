"""Focused dimensional invariants for the Task 13B modern full model.

These tests certify conversion identities and the declared physical coordinate
system.  They do not certify any placeholder kinetic value or a physiological
time course.
"""

from __future__ import annotations

import ast
from dataclasses import replace
import inspect
from pathlib import Path
import unittest

from modern_full_model import camp_pka, model, parameters, states, water
from modern_full_model.membranes import FMOL_PER_MOL, current_to_fmol_s


class TestAmountAndCurrentConversions(unittest.TestCase):
    def test_one_millimolar_times_one_picolitre_is_one_femtomole(self) -> None:
        # 1 mM = 1e-3 mol/L, 1 pL = 1e-12 L, and 1 mol = 1e15 fmol.
        self.assertEqual(1.0e-3 * 1.0e-12 * FMOL_PER_MOL, 1.0)
        compartment = states.compartment_from_concentrations(
            na_mM=1.0,
            k_mM=1.0,
            cl_mM=1.0,
            tic_mM=1.0,
            alkalinity_mM=1.0,
            volume_pL=1.0,
        )
        self.assertEqual(compartment.na_fmol, 1.0)
        self.assertEqual(compartment.concentrations().na_mM, 1.0)

    def test_faraday_conversion_and_anion_sign_are_exact_in_form(self) -> None:
        faraday = parameters.PhysicalConstants().faraday_C_mol
        one_fmol_monovalent_A = faraday / FMOL_PER_MOL
        self.assertAlmostEqual(
            current_to_fmol_s(
                one_fmol_monovalent_A,
                valence=+1,
                faraday_C_mol=faraday,
            ),
            1.0,
            places=15,
        )
        self.assertAlmostEqual(
            current_to_fmol_s(
                one_fmol_monovalent_A,
                valence=-1,
                faraday_C_mol=faraday,
            ),
            -1.0,
            places=15,
        )

    def test_one_picoampere_is_about_point_zero_one_fmol_per_second(self) -> None:
        faraday = parameters.PhysicalConstants().faraday_C_mol
        converted = current_to_fmol_s(
            1.0e-12,
            valence=+1,
            faraday_C_mol=faraday,
        )
        self.assertAlmostEqual(converted, 0.01036426965661773, places=15)


class TestWaterAndTimeCoordinates(unittest.TestCase):
    def test_palk_water_coefficients_use_the_exact_millimolar_conversion(self) -> None:
        # P [L^2 mol^-1 s^-1] * dc [mM] * 1e-3 [mol L^-1 mM^-1]
        # gives L/s; 1e12 pL/L therefore makes the coefficient P*1e9.
        source_L2_mol_s = {
            "apical_hydraulic_pL_s_mOsm": 4.32e-12,
            "basolateral_hydraulic_pL_s_mOsm": 5.15e-11,
            "paracellular_hydraulic_pL_s_mOsm": 2.60e-13,
        }
        expected_pL_s_mOsm = {
            "apical_hydraulic_pL_s_mOsm": 4.32e-3,
            "basolateral_hydraulic_pL_s_mOsm": 5.15e-2,
            "paracellular_hydraulic_pL_s_mOsm": 2.60e-4,
        }
        conversion = 1.0e-3 * 1.0e12
        base = parameters.FullModelParameters()
        records = {
            record.name: record
            for record in parameters.parameter_records(base)
        }
        for field_name, source_value in source_L2_mol_s.items():
            expected = expected_pL_s_mOsm[field_name]
            self.assertEqual(getattr(base.water, field_name), expected)
            self.assertAlmostEqual(source_value * conversion, expected, places=16)
            record = records[f"water.{field_name}"]
            self.assertEqual(record.units, "pL s^-1 (mOsm L^-1)^-1")
            self.assertEqual(record.provenance, "PUBLISHED_MODEL")

    def test_outflow_rate_times_volume_is_picolitres_per_second(self) -> None:
        base = parameters.FullModelParameters()
        water_parameters = replace(
            base.water,
            apical_hydraulic_pL_s_mOsm=0.0,
            basolateral_hydraulic_pL_s_mOsm=0.0,
            paracellular_hydraulic_pL_s_mOsm=0.0,
            lumen_dead_volume_pL=1.0,
            outflow_rate_s=2.0,
        )
        fluxes = water.evaluate_water_fluxes(
            water.OsmoticEnvironment(
                osm_cell_mOsm=300.0,
                osm_lumen_mOsm=300.0,
                osm_bath_mOsm=300.0,
                lumen_volume_pL=4.0,
            ),
            replace(base, water=water_parameters),
        )
        self.assertEqual(fluxes.lumen_outflow_pL_s, 6.0)
        # The same mM*pL=fmol identity makes q*c an amount rate.
        self.assertEqual(fluxes.lumen_outflow_pL_s * 3.0, 18.0)

    def test_modern_model_does_not_import_the_historical_rhs_scale(self) -> None:
        self.assertFalse(hasattr(parameters.FullModelParameters(), "rhs_scale"))
        package_directory = Path(inspect.getfile(model)).parent
        forbidden_references: list[str] = []
        for path in sorted(package_directory.glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.Name) and node.id == "rhs_scale":
                    forbidden_references.append(f"{path.name}:{node.lineno}")
                if isinstance(node, ast.Attribute) and node.attr == "rhs_scale":
                    forbidden_references.append(f"{path.name}:{node.lineno}")
        self.assertFalse(forbidden_references)

    def test_parameter_registry_has_no_unspecified_unit(self) -> None:
        records = parameters.parameter_records(parameters.FullModelParameters())
        self.assertTrue(records)
        self.assertFalse(
            [record.name for record in records if record.units == "UNSPECIFIED"]
        )


class TestRegulatoryDimensions(unittest.TestCase):
    def test_nested_regulatory_rhs_rows_are_fraction_per_second(self) -> None:
        r1 = camp_pka.R1EffectiveActivation(tau_activation_s=10.0)
        self.assertEqual(r1.rhs(2.0, (0.25,), 0.75), (0.05,))

        r2 = camp_pka.R2CampActivation(
            tau_camp_s=10.0,
            tau_activation_s=20.0,
        )
        dc_dt, dx_dt = r2.rhs(2.0, (0.20, 0.10), 0.60)
        self.assertAlmostEqual(dc_dt, 0.04, places=15)
        self.assertAlmostEqual(dx_dt, 0.005, places=15)

        r3 = camp_pka.R3PkaRegulatedFraction(
            tau_pka_s=10.0,
            forward_regulation_rate_s=0.05,
            reverse_regulation_rate_s=0.10,
        )
        dp_dt, dq_dt = r3.rhs(2.0, (0.20, 0.40), 0.60)
        self.assertAlmostEqual(dp_dt, 0.04, places=15)
        self.assertAlmostEqual(dq_dt, -0.026, places=15)

    def test_every_family_equilibrium_has_zero_regulatory_rhs(self) -> None:
        for family in camp_pka.default_regulatory_families().values():
            for beta_input in (0.0, 0.37, 1.0):
                equilibrium = family.equilibrium(beta_input)
                rates = family.rhs(0.0, equilibrium, beta_input)
                self.assertTrue(all(abs(rate) < 1.0e-15 for rate in rates))
                output = family.evaluate(0.0, equilibrium, beta_input)
                self.assertEqual(output.kinetic_evidence_status, "UNMEASURED_ACUTE_KINETICS")
                self.assertGreaterEqual(output.capacity_multiplier, 0.0)

    def test_regulatory_defaults_do_not_claim_measured_acute_kinetics(self) -> None:
        source = inspect.getsource(camp_pka)
        self.assertIn("UNMEASURED_ACUTE_KINETICS", source)
        self.assertIn("sensitivity values, not estimates", source)
        self.assertNotIn("rhs_scale", source)


if __name__ == "__main__":
    unittest.main()
