"""Round-5 tests for the nested cation-current chassis extension."""

from __future__ import annotations

import csv
import inspect
from pathlib import Path
import unittest

import numpy as np

from ae4_mechanism_reconstruction.chassis import (
    AE4_KNOCKOUT,
    BASELINE_STATE,
    FixedChassis,
    HistoricalSodiumOnlyAE4,
    ZeroAE4,
)
from state_resolved_ae4 import chassis as extended


KO_ROOT = np.asarray(
    (
        118.56527670449144,
        5.5593461776616655,
        65.66070346655056,
        21.132282246430893,
        123.70751315363782,
        48.782215027024655,
        58.90444835877707,
    )
)

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "results" / "13_state_resolved_ae4" / "chassis_reconstruction_scan.csv"


class TestNestedCationTopology(unittest.TestCase):
    def test_zero_apical_fractions_are_algebraically_nested(self) -> None:
        base = FixedChassis(HistoricalSodiumOnlyAE4())
        split = extended.SplitCationChassis(
            HistoricalSodiumOnlyAE4(), topology=extended.CationTopology(0.0, 0.0)
        )
        for calcium in (0.05, 0.55):
            raw_base, diagnostics_base = base.evaluate(
                BASELINE_STATE, calcium=calcium, scenario=AE4_KNOCKOUT
            )
            raw_split, diagnostics_split = split.evaluate(
                BASELINE_STATE, calcium=calcium, scenario=AE4_KNOCKOUT
            )
            np.testing.assert_allclose(raw_split, raw_base, rtol=0.0, atol=1e-14)
            for key in ("v_a", "v_b", "j_nak", "j_k", "j_t_na", "j_t_k"):
                self.assertAlmostEqual(diagnostics_split[key], diagnostics_base[key])

    def test_split_recloses_both_currents_and_lumen_charge(self) -> None:
        split = extended.SplitCationChassis(
            ZeroAE4(), topology=extended.CationTopology(0.2, 0.4)
        )
        _, diagnostics = split.evaluate(
            BASELINE_STATE, calcium=0.55, scenario=AE4_KNOCKOUT
        )
        checks = diagnostics["checks"]
        for key in (
            "apical_current_residual",
            "basolateral_current_residual",
            "proton_balance_identity_residual",
            "lumen_chloride_redundancy_residual",
            "pump_capacity_partition_residual",
            "k_conductance_partition_residual",
        ):
            self.assertLessEqual(abs(checks[key]), 2e-12, key)

    def test_apical_amount_terms_are_explicit(self) -> None:
        split = extended.SplitCationChassis(
            ZeroAE4(), topology=extended.CationTopology(0.3, 0.4)
        )
        _, diagnostics = split.evaluate(
            BASELINE_STATE, calcium=0.55, scenario=AE4_KNOCKOUT
        )
        self.assertAlmostEqual(
            diagnostics["lumen_na_source"],
            diagnostics["j_t_na"] + 3.0 * diagnostics["j_nak_apical"],
        )
        self.assertAlmostEqual(
            diagnostics["lumen_k_source"],
            diagnostics["j_t_k"]
            + diagnostics["j_k_apical"]
            - 2.0 * diagnostics["j_nak_apical"],
        )
        self.assertNotEqual(diagnostics["j_nak_apical"], 0.0)
        self.assertNotEqual(diagnostics["j_k_apical"], 0.0)

    def test_local_lumen_k_changes_apical_pump_turnover(self) -> None:
        split = extended.SplitCationChassis(
            ZeroAE4(), topology=extended.CationTopology(1.0, 0.0)
        )
        low_k = BASELINE_STATE.copy()
        high_k = BASELINE_STATE.copy()
        low_k[1] = 2.0
        high_k[1] = 12.0
        _, low = split.evaluate(low_k, calcium=0.05, scenario=AE4_KNOCKOUT)
        _, high = split.evaluate(high_k, calcium=0.05, scenario=AE4_KNOCKOUT)
        self.assertGreater(high["j_nak_apical"], low["j_nak_apical"])

    def test_continuation_recovers_nonzero_split_root(self) -> None:
        split = extended.SplitCationChassis(
            ZeroAE4(), topology=extended.CationTopology(0.05, 0.0)
        )
        roots = extended.solve_split_roots(split, preferred=KO_ROOT, start_count=1)
        self.assertEqual(len(roots), 1)
        self.assertLess(roots[0].max_abs_raw_rhs, 1e-8)
        self.assertTrue(all(value > 0.0 for value in roots[0].state))


class TestAcidBaseStructure(unittest.TestCase):
    def test_closed_carbon_and_buffer_reactions_conserve_mass_and_charge(self) -> None:
        result = extended.conserved_carbon_buffer_rates(
            co2_i=1.2,
            hco3_i=12.0,
            h_i=1.0e-4,
            buffer_h=9.0,
            buffer_minus=3.0,
            hydration_forward=0.2,
            hydration_reverse=20.0,
            buffer_dissociation=0.03,
            buffer_association=4.0,
        )
        self.assertAlmostEqual(result.inorganic_carbon_residual, 0.0)
        self.assertAlmostEqual(result.buffer_site_residual, 0.0)
        self.assertAlmostEqual(result.charge_residual, 0.0)

    def test_cl_ph_projection_cannot_identify_acid_base_fluxes(self) -> None:
        audit = extended.acid_base_identifiability_audit()
        self.assertEqual(audit["full_reaction_signature_rank"], 3)
        self.assertEqual(audit["released_cl_ph_direct_rank"], 1)
        self.assertEqual(audit["reaction_flux_nullity_under_cl_ph"], 2)
        self.assertEqual(audit["split_topology_direct_acid_base_rank"], 0)


class TestStageBFirewall(unittest.TestCase):
    def test_scan_rejects_target_access_before_stage_a_failure(self) -> None:
        with self.assertRaises(PermissionError):
            extended.scan_split_topologies(
                ZeroAE4(),
                stage_a_failed=False,
                cl_target_mM=36.5,
                pH_target=6.89,
                pump_fractions=(0.0,),
                k_fractions=(0.0,),
            )

    def test_chassis_source_contains_no_secretion_target(self) -> None:
        source = inspect.getsource(extended)
        self.assertNotIn("0.65", source)
        self.assertNotIn("35%", source)
        self.assertNotIn("35 +/-", source)


class TestMachineReadableScan(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with SCAN.open(newline="") as handle:
            cls.rows = list(csv.DictReader(handle))

    def test_every_row_preserves_secretion_holdout(self) -> None:
        self.assertTrue(self.rows)
        self.assertTrue(all(row["secretion_target_read"] == "False" for row in self.rows))

    def test_wt_attempt_fields_do_not_repurpose_stage_b_columns(self) -> None:
        rows = [row for row in self.rows if row["record_type"] == "wt_capacity_attempt"]
        self.assertEqual(len(rows), 5)
        for row in rows:
            self.assertEqual(row["cl_target_residual_mM"], "")
            self.assertEqual(row["pH_target_residual"], "")
            self.assertEqual(row["normalized_ionic_distance"], "")
            self.assertNotEqual(row["wt_cl_residual_mM"], "")
            self.assertNotEqual(row["wt_pH_z"], "")
            self.assertNotEqual(row["optimizer_cost"], "")
            self.assertNotEqual(row["attempted_capacity"], "")
            self.assertEqual(row["calibration_valid"], "False")
            self.assertEqual(row["wt_gate2_pass"], "False")
            self.assertEqual(row["derived_height_boundary_hit"], "True")
            self.assertEqual(
                row["ae2_test_status"], "N/A_upstream_WT_calibration_failed"
            )

    def test_all_resolved_null_sensitivities_remain_far_from_cl_ph(self) -> None:
        rows = [row for row in self.rows if row["record_type"] == "resolved_null_root"]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(float(row["max_abs_raw_rhs"]) <= 1e-8 for row in rows))
        self.assertGreater(min(float(row["normalized_ionic_distance"]) for row in rows), 20.0)
        self.assertTrue(all(row["valid_wt_capacity_available"] == "False" for row in rows))

    def test_cation_split_cannot_change_fixed_slice_hco3_row(self) -> None:
        rows = [row for row in self.rows if row["record_type"] == "conditional_target_slice"]
        for slice_name in ("model_anchored", "wt_reference_cations"):
            hco3_rows = {
                round(float(row["raw_rhs_vector"].split(";")[6]), 14)
                for row in rows
                if row["conditional_slice"] == slice_name
            }
            self.assertEqual(len(hco3_rows), 1)
            self.assertNotEqual(next(iter(hco3_rows)), 0.0)


if __name__ == "__main__":
    unittest.main()
