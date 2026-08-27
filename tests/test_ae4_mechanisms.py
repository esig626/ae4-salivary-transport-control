import math
from types import SimpleNamespace
import unittest

from ae4_mechanism_reconstruction.mechanisms import (
    AE4Context,
    AE4Parameters,
    CANDIDATE_SPECS,
    PKA_FOLD_ASSAY_VALUES,
    candidate_specs,
    equilibrium_context,
    evaluate_ae4,
    get_mechanism,
)


class TestAE4MechanismFamily(unittest.TestCase):
    def setUp(self):
        self.historical = AE4Context(
            na_i=25.0,
            k_i=120.0,
            cl_i=50.0,
            hco3_i=10.0,
            na_e=140.2,
            k_e=5.3,
            cl_e=102.6,
            hco3_e=40.0,
        )

    def test_registry_covers_required_families_and_stoichiometries(self):
        families = {spec.family.split()[0] for spec in candidate_specs()}
        self.assertEqual(families, {f"C{i}" for i in range(1, 9)})
        c6_stoich = {
            spec.stoichiometry
            for spec in candidate_specs()
            if spec.family.startswith("C6") and spec.active
        }
        self.assertEqual(c6_stoich, {(1, 1, 2), (1, 2, 3), (2, 1, 3)})

    def test_c8_mass_action_ids_are_declared_c6_aliases(self):
        for suffix in ("112", "123", "213"):
            with self.subTest(suffix=suffix):
                self.assertIsNone(CANDIDATE_SPECS[f"C6_{suffix}"].equivalent_to)
                self.assertEqual(
                    CANDIDATE_SPECS[f"C8_{suffix}_MA"].equivalent_to,
                    f"C6_{suffix}",
                )
                self.assertIsNone(CANDIDATE_SPECS[f"C8_{suffix}_SAT"].equivalent_to)

    def test_c3_round1_complexity_fixes_the_primary_rmax_ratio(self):
        spec = CANDIDATE_SPECS["C3"]
        self.assertEqual(spec.free_parameters, ("activity_scale",))
        self.assertTrue(
            any("1.5:1.6" in item for item in spec.independent_constraints)
        )

    def test_pka_candidates_are_predeclared_with_a_neutral_rest_input(self):
        expected = {
            "C1_PKA_COMMON": ("C1", "pka_common"),
            "C1_PKA_K_RECRUIT": ("C1", "pka_k_recruit"),
            "C3_PKA_COMMON": ("C3", "pka_common"),
            "C6_112_PKA_COMMON": ("C6_112", "pka_common"),
            "C6_123_PKA_COMMON": ("C6_123", "pka_common"),
            "C6_213_PKA_COMMON": ("C6_213", "pka_common"),
        }
        self.assertEqual(PKA_FOLD_ASSAY_VALUES, (1.25, 1.6, 3.0))
        for candidate_id, (core, mode) in expected.items():
            with self.subTest(candidate=candidate_id):
                spec = CANDIDATE_SPECS[candidate_id]
                self.assertEqual(spec.core_candidate_id, core)
                self.assertEqual(spec.stimulus_mode, mode)
                self.assertEqual(spec.investigation_round, 2)
                stimulated = evaluate_ae4(
                    candidate_id, self.historical, stimulus_activation=1.0
                )
                resting = evaluate_ae4(candidate_id, self.historical)
                self.assertNotEqual(stimulated.cycle_flux, resting.cycle_flux)
        self.assertIn("inadmissible", CANDIDATE_SPECS["C1_PKA_COMMON"].assay_status)
        self.assertIn("novel Round-2", CANDIDATE_SPECS["C1_PKA_K_RECRUIT"].assay_status)
        self.assertIn("direct-source", CANDIDATE_SPECS["C3_PKA_COMMON"].assay_status)
        self.assertIn(
            "thermodynamics-admissible",
            CANDIDATE_SPECS["C6_112_PKA_COMMON"].assay_status,
        )
        for candidate_id in ("C6_123_PKA_COMMON", "C6_213_PKA_COMMON"):
            self.assertIn(
                "thermodynamics-admissible",
                CANDIDATE_SPECS[candidate_id].assay_status,
            )

    def test_pka_common_gate_preserves_rest_and_scales_the_whole_source(self):
        parameters = AE4Parameters(pka_stimulated_fold=3.0)
        for candidate_id, core_id in (
            ("C1_PKA_COMMON", "C1"),
            ("C3_PKA_COMMON", "C3"),
            ("C6_112_PKA_COMMON", "C6_112"),
            ("C6_123_PKA_COMMON", "C6_123"),
            ("C6_213_PKA_COMMON", "C6_213"),
        ):
            with self.subTest(candidate=candidate_id):
                core = evaluate_ae4(core_id, self.historical, parameters)
                rest = evaluate_ae4(candidate_id, self.historical, parameters)
                full = evaluate_ae4(
                    candidate_id,
                    self.historical,
                    parameters,
                    stimulus_activation=1.0,
                )
                self.assertEqual(rest, core)
                for field in ("cycle_flux", "na_i", "k_i", "cl_i", "hco3_i"):
                    self.assertAlmostEqual(
                        getattr(full, field), 3.0 * getattr(core, field)
                    )
                self.assertEqual(full.diagnostics["pka_activity_multiplier"], 3.0)
                for branch in ("na", "k"):
                    key = f"{branch}_flux"
                    if key in core.diagnostics:
                        self.assertAlmostEqual(
                            full.diagnostics[key],
                            3.0 * core.diagnostics[key],
                        )

    def test_pka_recruited_k_branch_has_the_declared_balance(self):
        parameters = AE4Parameters(pka_stimulated_fold=3.0)
        basal = evaluate_ae4("C1", self.historical, parameters)
        rest = evaluate_ae4("C1_PKA_K_RECRUIT", self.historical, parameters)
        full = evaluate_ae4(
            "C1_PKA_K_RECRUIT",
            self.historical,
            parameters,
            stimulus_activation=1.0,
        )
        self.assertEqual(rest, basal)
        self.assertAlmostEqual(full.diagnostics["na_flux"], basal.cycle_flux)
        self.assertAlmostEqual(full.diagnostics["k_flux"], 2.0 * basal.cycle_flux)
        self.assertAlmostEqual(full.na_i, -basal.cycle_flux)
        self.assertAlmostEqual(full.k_i, -2.0 * basal.cycle_flux)
        self.assertAlmostEqual(full.cl_i, 3.0 * basal.cycle_flux)
        self.assertAlmostEqual(full.hco3_i, -6.0 * basal.cycle_flux)
        self.assertAlmostEqual(full.transported_charge, 0.0)
        self.assertFalse(full.diagnostics["thermodynamic_valid"])

    def test_pka_input_is_independent_and_knockout_still_annihilates_flux(self):
        context = SimpleNamespace(
            **{**self.historical.__dict__, "pka_activation": 0.5}
        )
        parameters = AE4Parameters(activity_scale=2.0, pka_stimulated_fold=3.0)
        implicit = evaluate_ae4("C1_PKA_COMMON", context, parameters)
        explicit = evaluate_ae4(
            "C1_PKA_COMMON",
            self.historical,
            parameters,
            stimulus_activation=0.5,
        )
        self.assertEqual(implicit, explicit)
        for candidate_id in (
            "C1_PKA_COMMON",
            "C1_PKA_K_RECRUIT",
            "C3_PKA_COMMON",
            "C6_112_PKA_COMMON",
            "C6_123_PKA_COMMON",
            "C6_213_PKA_COMMON",
        ):
            with self.subTest(candidate=candidate_id):
                knockout = evaluate_ae4(
                    candidate_id,
                    self.historical,
                    parameters,
                    activity_scale=0.0,
                    stimulus_activation=1.0,
                )
                self.assertEqual(knockout.cycle_flux, 0.0)
                self.assertEqual(
                    (knockout.na_i, knockout.k_i, knockout.cl_i, knockout.hco3_i),
                    (0.0, 0.0, 0.0, 0.0),
                )

    def test_pka_gate_preserves_c6_equilibrium_and_affinity_sign(self):
        parameters = AE4Parameters(
            pka_stimulated_fold=3.0,
            na_weight=1.0,
            k_weight=0.0,
        )
        for suffix, stoichiometry in (
            ("112", (1, 1, 2)),
            ("123", (1, 2, 3)),
            ("213", (2, 1, 3)),
        ):
            candidate_id = f"C6_{suffix}_PKA_COMMON"
            with self.subTest(candidate=candidate_id):
                equilibrium = equilibrium_context(
                    self.historical, "na", stoichiometry
                )
                at_rest = evaluate_ae4(
                    candidate_id,
                    equilibrium,
                    parameters,
                    stimulus_activation=0.0,
                )
                stimulated = evaluate_ae4(
                    candidate_id,
                    equilibrium,
                    parameters,
                    stimulus_activation=1.0,
                )
                self.assertTrue(
                    CANDIDATE_SPECS[candidate_id].thermodynamic_sign_controlled
                )
                self.assertAlmostEqual(at_rest.cycle_flux, 0.0, places=8)
                self.assertAlmostEqual(stimulated.cycle_flux, 0.0, places=8)

                forward = AE4Context(
                    **{
                        **equilibrium.__dict__,
                        "hco3_i": equilibrium.hco3_i * 1.01,
                    }
                )
                reverse = AE4Context(
                    **{
                        **equilibrium.__dict__,
                        "hco3_i": equilibrium.hco3_i * 0.99,
                    }
                )
                self.assertGreater(
                    evaluate_ae4(
                        candidate_id,
                        forward,
                        parameters,
                        stimulus_activation=1.0,
                    ).cycle_flux,
                    0.0,
                )
                self.assertLess(
                    evaluate_ae4(
                        candidate_id,
                        reverse,
                        parameters,
                        stimulus_activation=1.0,
                    ).cycle_flux,
                    0.0,
                )

    def test_invalid_pka_controls_fail_loudly(self):
        with self.assertRaisesRegex(ValueError, r"stimulus_activation.*\[0, 1\]"):
            evaluate_ae4(
                "C1_PKA_COMMON", self.historical, stimulus_activation=1.01
            )
        with self.assertRaisesRegex(ValueError, "pka_stimulated_fold"):
            evaluate_ae4(
                "C1_PKA_COMMON",
                self.historical,
                AE4Parameters(pka_stimulated_fold=0.99),
            )

    def test_all_active_candidates_are_electroneutral(self):
        for spec in candidate_specs(include_rejected=False):
            with self.subTest(candidate=spec.candidate_id):
                out = evaluate_ae4(spec.candidate_id, self.historical)
                self.assertTrue(out.valid)
                self.assertAlmostEqual(out.transported_charge, 0.0, places=12)
                charge_from_terms = out.na_i + out.k_i - out.cl_i - out.hco3_i
                self.assertAlmostEqual(charge_from_terms, 0.0, places=12)

    def test_historical_c1_reproduces_source_turnover_and_balance(self):
        out = evaluate_ae4("C1", self.historical)
        self.assertAlmostEqual(out.cycle_flux, 0.002051862990, places=11)
        self.assertAlmostEqual(out.na_i, -out.cycle_flux, places=14)
        self.assertEqual(out.k_i, 0.0)
        self.assertAlmostEqual(out.cl_i, out.cycle_flux, places=14)
        self.assertAlmostEqual(out.hco3_i, -2.0 * out.cycle_flux, places=14)

        scaled = evaluate_ae4(
            "C1", self.historical, AE4Parameters(activity_scale=18.9502607221547)
        )
        self.assertAlmostEqual(scaled.cycle_flux, 0.03888333862, places=9)

    def test_published_c2_uses_exact_pooled_law_and_internal_partition(self):
        out = evaluate_ae4("C2", self.historical)
        ci = self.historical.na_i + self.historical.k_i
        ce = self.historical.na_e + self.historical.k_e
        expected = (
            1.92e-2
            * self.historical.cl_e
            * self.historical.hco3_i**2
            * ci
            - 1.3e-5
            * self.historical.cl_i
            * self.historical.hco3_e**2
            * ce
        )
        self.assertAlmostEqual(out.cycle_flux, expected, places=9)
        self.assertAlmostEqual(out.na_i, -expected * self.historical.na_i / ci)
        self.assertAlmostEqual(out.k_i, -expected * self.historical.k_i / ci)
        self.assertFalse(out.diagnostics["thermodynamic_valid"])

    def test_parallel_c3_keeps_cation_branches_in_their_own_balances(self):
        out = evaluate_ae4("C3", self.historical)
        self.assertAlmostEqual(out.na_i, -out.diagnostics["na_flux"])
        self.assertAlmostEqual(out.k_i, -out.diagnostics["k_flux"])
        self.assertAlmostEqual(out.cl_i, -(out.na_i + out.k_i))
        self.assertNotEqual(out.na_i, out.k_i)

    def test_primary_hill_values_are_frozen_in_c4_and_c5_defaults(self):
        p = AE4Parameters()
        self.assertEqual((p.ec50_na, p.hill_na), (49.0, 2.0))
        self.assertEqual((p.ec50_k, p.hill_k), (62.0, 1.8))
        self.assertFalse(CANDIDATE_SPECS["C4"].thermodynamic_sign_controlled)
        self.assertTrue(CANDIDATE_SPECS["C5"].thermodynamic_sign_controlled)

    def test_pure_gating_without_cation_transport_is_rejected(self):
        spec = CANDIDATE_SPECS["C5_GATE_ONLY_REJECTED"]
        self.assertFalse(spec.active)
        self.assertIn("direct Na and K flux", spec.rejection_reason)
        with self.assertRaisesRegex(ValueError, "analytically rejected"):
            evaluate_ae4(spec.candidate_id, self.historical)

    def test_unresolved_primary_chemistries_are_recorded_not_silently_proxied(self):
        carbonate = CANDIDATE_SPECS["C6_CARBONATE_UNTESTABLE"]
        self_exchange = CANDIDATE_SPECS["C6_HCO3_SELF_EXCHANGE_UNTESTABLE"]
        self.assertFalse(carbonate.active)
        self.assertIsNone(carbonate.stoichiometry)
        self.assertIn("no carbonate state", carbonate.rejection_reason)
        self.assertFalse(self_exchange.active)
        self.assertIsNone(self_exchange.stoichiometry)
        self.assertIn("unresolved stoichiometry", self_exchange.rejection_reason)

    def test_stoichiometric_balance_signatures(self):
        for candidate, stoich in (
            ("C6_112", (1, 1, 2)),
            ("C6_123", (1, 2, 3)),
            ("C6_213", (2, 1, 3)),
        ):
            with self.subTest(candidate=candidate):
                out = evaluate_ae4(
                    candidate,
                    self.historical,
                    AE4Parameters(na_weight=1.0, k_weight=0.0),
                )
                a, b, c = stoich
                self.assertAlmostEqual(out.cl_i, a * out.cycle_flux)
                self.assertAlmostEqual(out.na_i, -b * out.cycle_flux)
                self.assertEqual(out.k_i, -0.0)
                self.assertAlmostEqual(out.hco3_i, -c * out.cycle_flux)
                self.assertEqual(c, a + b)

    def test_reversible_laws_zero_at_affinity_zero_and_reverse_sign(self):
        candidates = [
            ("C6_112", (1, 1, 2)),
            ("C6_123", (1, 2, 3)),
            ("C6_213", (2, 1, 3)),
            ("C7_112", (1, 1, 2)),
            ("C7_123", (1, 2, 3)),
            ("C7_213", (2, 1, 3)),
            ("C8_112_SAT", (1, 1, 2)),
            ("C8_123_SAT", (1, 2, 3)),
            ("C8_213_SAT", (2, 1, 3)),
        ]
        p = AE4Parameters(na_weight=1.0, k_weight=0.0)
        for candidate, stoich in candidates:
            with self.subTest(candidate=candidate):
                eq = equilibrium_context(self.historical, "na", stoich)
                at_eq = evaluate_ae4(candidate, eq, p)
                self.assertAlmostEqual(at_eq.diagnostics["affinity_na"], 0.0, places=12)
                self.assertAlmostEqual(at_eq.cycle_flux, 0.0, places=8)

                more_inside_bicarb = AE4Context(
                    **{**eq.__dict__, "hco3_i": eq.hco3_i * 1.01}
                )
                less_inside_bicarb = AE4Context(
                    **{**eq.__dict__, "hco3_i": eq.hco3_i * 0.99}
                )
                self.assertGreater(evaluate_ae4(candidate, more_inside_bicarb, p).cycle_flux, 0.0)
                self.assertLess(evaluate_ae4(candidate, less_inside_bicarb, p).cycle_flux, 0.0)

    def test_saturation_preserves_sign_and_reduces_magnitude(self):
        for suffix in ("112", "123", "213"):
            with self.subTest(stoichiometry=suffix):
                ma = evaluate_ae4(f"C8_{suffix}_MA", self.historical)
                sat = evaluate_ae4(f"C8_{suffix}_SAT", self.historical)
                self.assertEqual(math.copysign(1.0, ma.cycle_flux), math.copysign(1.0, sat.cycle_flux))
                self.assertLess(abs(sat.cycle_flux), abs(ma.cycle_flux))

    def test_genotype_multiplier_preserves_capacity_for_wt_and_zeroes_ko(self):
        mechanism = get_mechanism("C1", AE4Parameters(activity_scale=18.9502607221547))
        direct = mechanism.evaluate(self.historical)
        wt = mechanism.evaluate(self.historical, activity_scale=1.0)
        ko = mechanism.evaluate(self.historical, activity_scale=0.0)
        self.assertEqual(wt.cycle_flux, direct.cycle_flux)
        self.assertEqual(ko.cycle_flux, 0.0)
        self.assertEqual((ko.na_i, ko.k_i, ko.cl_i, ko.hco3_i), (-0.0, -0.0, 0.0, -0.0))

    def test_invalid_inputs_and_unknown_candidate_fail_loudly(self):
        with self.assertRaisesRegex(ValueError, "na_i"):
            evaluate_ae4("C1", AE4Context(0, 1, 1, 1, 1, 1, 1, 1))
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            evaluate_ae4("C1", self.historical, AE4Parameters(activity_scale=-1))
        with self.assertRaises(KeyError):
            evaluate_ae4("C99", self.historical)


if __name__ == "__main__":
    unittest.main()
