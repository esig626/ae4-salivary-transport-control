"""Regression and invariant tests for the fixed Task-12 chassis."""

from __future__ import annotations

import unittest

import numpy as np

from ae4_mechanism_reconstruction.chassis import (
    AE2_KNOCKOUT,
    AE4_KNOCKOUT,
    BASELINE_STATE,
    WT,
    AE4Balance,
    ChassisParameters,
    FixedChassis,
    dimension_ledger,
)


class TestFixedAE4Chassis(unittest.TestCase):
    def setUp(self) -> None:
        self.model = FixedChassis()

    def test_fixed_non_ae4_parameter_fingerprint(self) -> None:
        self.assertEqual(
            self.model.parameter_fingerprint,
            "55885266c44ecbbfe86be622c6ebd08b77151ee6ef7dfd5c8bbc8d4250b58e50",
        )
        self.assertEqual(self.model.parameter_fingerprint, ChassisParameters().fingerprint())

    def test_historical_baseline_closes_independently(self) -> None:
        raw, diagnostics = self.model.evaluate(BASELINE_STATE, calcium=0.05)
        self.assertLess(float(np.max(np.abs(raw))), 5.0e-12)
        self.assertAlmostEqual(diagnostics["q_total"], 0.000531415476760618, places=15)
        self.assertAlmostEqual(diagnostics["v_a"], -50.24, places=11)
        self.assertAlmostEqual(diagnostics["v_b"], -62.8, places=11)
        self.assertAlmostEqual(diagnostics["ae4_cycle_flux"], 0.0388833386248, places=12)

    def test_exact_reduction_and_conservation_checks(self) -> None:
        _, diagnostics = self.model.evaluate(BASELINE_STATE, calcium=0.05)
        checks = diagnostics["checks"]
        for name, value in checks.items():
            with self.subTest(name=name):
                self.assertLess(abs(float(value)), 1.0e-10)
        self.assertEqual(diagnostics["cl_l"], BASELINE_STATE[0] + BASELINE_STATE[1])

    def test_baseline_ae4_demand_is_held_out_free(self) -> None:
        demand = self.model.ae4_source_demand()
        source = demand["required_ae4_source"]
        self.assertEqual(demand["calibration_information"], "WT baseline closure only; no KO datum")
        self.assertAlmostEqual(source["na_i"], -0.0388833386254, places=12)
        self.assertAlmostEqual(source["k_i"], 0.0, places=12)
        self.assertAlmostEqual(source["cl_i"], 0.0388833386765, places=12)
        self.assertAlmostEqual(source["hco3_i"], -0.0777666773911, places=12)
        self.assertLess(abs(demand["required_charge_rate"]), 1.0e-10)

    def test_published_na_k_partition_breaks_only_cation_closure(self) -> None:
        mismatch = self.model.baseline_partition_mismatch()
        residual = mismatch["amount_balance_residual"]
        self.assertAlmostEqual(mismatch["na_fraction"], 25.0 / 145.0, places=15)
        self.assertAlmostEqual(mismatch["k_fraction"], 120.0 / 145.0, places=15)
        self.assertAlmostEqual(mismatch["expected_na_k_magnitude"], 0.03217931477, places=10)
        self.assertAlmostEqual(residual["na_i"], +0.03217931472, places=10)
        self.assertAlmostEqual(residual["k_i"], -0.03217931477, places=10)
        self.assertLess(abs(residual["cl_i"]), 1.0e-12)
        self.assertLess(abs(residual["hco3_i"]), 1.0e-10)

    def test_rejects_electrogenic_ae4_insertion(self) -> None:
        def invalid_law(_environment, _scale):
            return AE4Balance(
                cycle_flux=1.0,
                na_i=0.0,
                k_i=0.0,
                cl_i=1.0,
                hco3_i=-1.0,
                transported_charge=1.0,
            )

        with self.assertRaisesRegex(ValueError, "not electroneutral"):
            FixedChassis(invalid_law).evaluate(BASELINE_STATE, calcium=0.05)

    def test_radau_regresses_historical_knockout_ratios(self) -> None:
        grid = {"t_eval": (0.0, 100.0, 200.0)}
        wt = self.model.simulate(WT, method="BDF", **grid)
        ae2 = self.model.simulate(AE2_KNOCKOUT, method="BDF", **grid)
        # Radau agrees with the much slower BDF AE4-KO endpoint to better than
        # 1e-8 relative and remains solver-safe on the knockout transient.
        ae4 = self.model.simulate(AE4_KNOCKOUT, method="Radau", **grid)
        self.assertAlmostEqual(ae2.endpoint["q_total"] / wt.endpoint["q_total"], 0.999052, places=5)
        self.assertAlmostEqual(ae4.endpoint["q_total"] / wt.endpoint["q_total"], 0.998556, places=5)
        self.assertEqual(wt.parameter_fingerprint, ae2.parameter_fingerprint)
        self.assertEqual(wt.parameter_fingerprint, ae4.parameter_fingerprint)

    def test_dimensional_boundary_is_explicit(self) -> None:
        ledger = dimension_ledger()
        self.assertEqual(ledger["state_concentrations"], "mM")
        self.assertEqual(ledger["absolute_flow"], "not certified")
        self.assertIn("dimensionless WT/KO ratios", ledger["valid_primary_comparison"])


if __name__ == "__main__":
    unittest.main()
