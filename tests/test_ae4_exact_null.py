"""Focused Task 26 checks. These tests do not run simulations."""
import itertools
import json
import unittest

import numpy as np

from modern_full_model.calibration import WTCalibrationSpec
from modern_full_model.genotype_evaluation import _root_state_with_basal_regulation
from modern_full_model.run_ae4_exact_null import (
    CALCIUM, EXPRESSIONS, OUT, PAIR_FIELDS, REPO, SOURCE_NAMES, TASK24,
    TEXT_OUTPUTS, genotype, inputs, read_rows, rest_flux_metrics, saved_rest,
    source_values,
)
from modern_full_model.run_ae4_routing_re_equilibrated import _explicit_rest_model
from modern_full_model.validation import (
    CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU, sha256_file, sha256_object,
)


class ExactNullTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.frozen, _ = inputs()
        cls.roots = sorted(cls.manifest["roots"])
        cls.rest = read_rows(OUT / "resting_states.csv")
        cls.traces = read_rows(OUT / "trajectory_summary.csv")
        cls.pairs = read_rows(OUT / "comparison.csv")
        cls.flux = read_rows(OUT / "flux_decomposition.csv")
        cls.attempts = read_rows(OUT / "warm_start_attempts.csv")
        cls.contract = json.loads((OUT / "contract.json").read_text())
        cls.verification = json.loads((OUT / "verification.json").read_text())

    def test_complete_panel_without_five_percent_simulation(self):
        self.assertEqual(len(self.rest), 20)
        self.assertEqual(len(self.traces), 60)
        self.assertEqual(len(self.pairs), 30)
        self.assertEqual(len(self.flux), 10)
        self.assertEqual({(r["root_id"], float(r["ae4_expression"])) for r in self.rest}, set(itertools.product(self.roots, EXPRESSIONS)))
        self.assertEqual({(r["root_id"], float(r["ae4_expression"]), float(r["calcium_uM"])) for r in self.traces}, set(itertools.product(self.roots, EXPRESSIONS, CALCIUM)))
        self.assertEqual({(r["root_id"], float(r["calcium_uM"])) for r in self.pairs}, set(itertools.product(self.roots, CALCIUM)))
        with self.assertRaises(ValueError):
            genotype(0.05)

    def test_wt_rest_is_bitwise_unchanged_and_calcium_uses_one_rest(self):
        for row in self.rest:
            expression = float(row["ae4_expression"])
            if row["numerical_rest_gate_pass"] != "True":
                self.assertEqual(expression, 0.0)
                self.assertEqual(row["core_state_json"], "")
                continue
            core = json.loads(row["core_state_json"])
            self.assertEqual(sha256_object(core), row["core_state_sha256"])
            if expression == 1.0:
                saved, _ = saved_rest(row["root_id"], 1.0)
                np.testing.assert_array_equal(core, saved["core_state"])
            group = [r for r in self.traces if r["root_id"] == row["root_id"] and float(r["ae4_expression"]) == expression]
            self.assertEqual(len(group), 3)
            self.assertEqual({r["resting_state_sha256"] for r in group}, {row["core_state_sha256"]})
            self.assertTrue(all(r["onset_rest_rhs_match"] == "True" for r in group))

    def test_production_sources_and_parameters_remain_frozen(self):
        for relative, expected in self.contract["source_sha256"].items():
            self.assertEqual(sha256_file(REPO / relative), expected, relative)
        self.assertTrue(self.verification["production_sources_unchanged"])
        for row in self.rest + self.traces:
            if not row.get("whole_cell_parameters_sha256"):
                self.assertTrue(row.get("rest_status")=="unresolved_no_closed_null_rest" or row.get("run_status")=="not_run_no_closed_null_rest")
                continue
            frozen = self.manifest["roots"][row["root_id"]]
            self.assertEqual(row["whole_cell_parameters_sha256"], frozen["whole_cell_parameters_sha256"])
            self.assertEqual(row["ae4_parameters_sha256"], frozen["ae4_parameters_sha256"])

    def test_exact_zero_sources_at_wt_and_null_states_without_division_by_zero(self):
        for root in self.roots:
            model = _explicit_rest_model(self.manifest, root)
            for row in [r for r in self.rest if r["root_id"] == root]:
                if row["numerical_rest_gate_pass"] != "True":
                    continue
                core = json.loads(row["core_state_json"])
                state = _root_state_with_basal_regulation(model, core)
                evaluation = model.evaluate(0.0, state, genotype=genotype(0.0))
                self.assertEqual(set(source_values(evaluation)), set(SOURCE_NAMES))
                self.assertTrue(all(value == 0.0 for value in source_values(evaluation).values()))
                metrics = rest_flux_metrics(model, core, 0.0)
                self.assertEqual(metrics["ae4_cl_flux_fmol_s"], 0.0)
                self.assertTrue(metrics["null_ae4_sources_exact_zero"])
        for row in self.traces:
            if float(row["ae4_expression"]) == 0.0 and row["run_status"] == "completed":
                for name in SOURCE_NAMES:
                    self.assertEqual(float(row[f"max_abs_ae4_{name}"]), 0.0)

    def test_true_rest_full_closure_and_boundaries(self):
        spec = WTCalibrationSpec()
        for row in self.rest:
            if row["rest_status"] == "unresolved_no_closed_null_rest":
                self.assertEqual(row["rest_gate_pass"], "False")
                self.assertEqual(row["numerical_rest_gate_pass"], "False")
                self.assertEqual(row["core_state_json"], "")
                continue
            self.assertEqual(row["numerical_rest_gate_pass"], "True")
            self.assertEqual(row["rest_protocol"], "REST")
            self.assertEqual(float(row["rest_beta_input"]), 0.0)
            self.assertEqual(row["rest_integration_used"], "False")
            self.assertLessEqual(float(row["max_abs_scaled_independent_rhs"]), spec.root_scaled_tolerance)
            self.assertLessEqual(float(row["max_abs_omitted_rhs_fmol_s"]), spec.omitted_row_raw_tolerance)
            self.assertLessEqual(float(row["max_abs_regulatory_rhs_s_inv"]), 1e-10)
            self.assertLessEqual(float(row["max_abs_bulk_charge_fmol"]), spec.charge_tolerance_fmol)
            self.assertLessEqual(float(row["max_abs_current_residual_A"]), spec.current_tolerance_A)
            coords = np.asarray(json.loads(row["coordinates_json"]))
            distance = np.minimum(coords-spec.coordinate_lower, spec.coordinate_upper-coords)/(spec.coordinate_upper-spec.coordinate_lower)
            expected_hits = [bound.name for bound, value in zip(spec.coordinate_bounds, distance) if value <= spec.boundary_relative_tolerance]
            self.assertEqual(json.loads(row["boundary_hits_json"]), expected_hits)
            self.assertEqual(row["physiological_domain_gate_pass"] == "True", not expected_hits)
            self.assertEqual(row["rest_gate_pass"] == "True", not expected_hits)
            self.assertAlmostEqual(float(row["minimum_normalized_boundary_distance"]), float(min(distance)), places=12)
            self.assertEqual(int(float(row["independent_jacobian_rank"])), 10)
            self.assertEqual(row["positivity_gate_pass"], "True")
            residuals = json.loads(row["conservation_residuals_json"])
            self.assertEqual(set(residuals), set(CONSERVATION_RESIDUAL_TOLERANCES))
            for key, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items():
                self.assertLessEqual(abs(residuals[key]), tolerance)

    def test_both_null_warm_starts_are_preserved_and_agree(self):
        self.assertEqual(len(self.attempts), 50)
        for root in self.roots:
            initial = [r for r in self.attempts if r["root_id"] == root and r["method"] == "inherited_bounded_direct_least_squares"]
            self.assertEqual(len(initial), 2)
            self.assertTrue(all(r["optimizer_success"] == "False" for r in initial))
            self.assertTrue(all(int(float(r["nfev"])) == 3000 for r in initial))
            rows = [r for r in self.attempts if r["root_id"] == root and r["method"] == "direct_hybr_with_fixed_initial_Jacobian_preconditioning"]
            self.assertEqual({r["seed_id"] for r in rows}, {"matched_Task24_WT", "matched_Task24_5pct"})
            self.assertTrue(all(float(r["ae4_expression"]) == 0.0 for r in rows))
            selected = next(r for r in self.rest if r["root_id"] == root and float(r["ae4_expression"]) == 0.0)
            if selected["numerical_rest_gate_pass"] == "True":
                self.assertTrue(all(r["numerical_rest_gate_pass"] == "True" for r in rows))
                self.assertEqual(sum(r["selected"] == "True" for r in rows), 1)
                self.assertLessEqual(float(selected["warm_start_root_distance_normalized"]), WTCalibrationSpec().cluster_relative_tolerance)
            else:
                fallback = [r for r in self.attempts if r["root_id"]==root and r["method"]=="conditioned_reciprocal_volume_diagnostic"]
                self.assertEqual(len(fallback), 2)
                self.assertTrue(all(r["numerical_rest_gate_pass"] == "False" for r in rows+fallback))

    def test_all_stimulated_trajectories_meet_inherited_gates(self):
        for row in self.traces:
            if row["run_status"] != "completed":
                self.assertEqual(row["run_status"], "not_run_no_closed_null_rest")
                self.assertEqual(row["valid"], "False")
                self.assertEqual(row["total_0_600_pL"], "")
                continue
            for field in ("valid", "solver_success", "exact_time_grid", "positive_core", "finite_quantities", "nonnegative_flow", "onset_rest_rhs_match"):
                self.assertEqual(row[field], "True", (row["root_id"], field))
            self.assertEqual(row["solver"], PRODUCTION_RADAU.label)
            self.assertEqual(int(row["sample_count"]), len(self.frozen["time_grid_s"]))
            self.assertEqual(float(row["duration_s"]), 600.0)
            self.assertLessEqual(float(row["max_dimensionless_conservation_ratio"]), 1.0)
            for key, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items():
                self.assertLessEqual(float(row[f"max_abs_conservation_{key}"]), tolerance)
            self.assertEqual(int(row["opposing_cation_samples"]), 0)
            self.assertEqual(row["trajectory_arrays_written"], "False")
            self.assertEqual(row["admissible_trajectory"], row["matched_rest_admissible"])

    def test_nkcc1_stoichiometry_and_replacement_use_actual_chloride_flux(self):
        index = {(r["root_id"], float(r["ae4_expression"])): r for r in self.rest}
        for row in self.rest:
            if row["numerical_rest_gate_pass"] != "True":
                continue
            self.assertEqual(float(row["nkcc1_chloride_influx_fmol_s"]), 2*float(row["nkcc1_inward_cycle_flux_fmol_s"]))
        for row in self.flux:
            if row["numerical_pair_available"] != "True":
                self.assertEqual(row["nkcc1_replacement_fraction"], "")
                continue
            wt, null = index[row["root_id"], 1.0], index[row["root_id"], 0.0]
            loss = float(wt["ae4_cl_flux_fmol_s"]) - float(null["ae4_cl_flux_fmol_s"])
            gain = float(null["nkcc1_chloride_influx_fmol_s"]) - float(wt["nkcc1_chloride_influx_fmol_s"])
            self.assertGreater(loss, 0.0)
            self.assertAlmostEqual(float(row["nkcc1_replacement_fraction"]), gain/loss, places=14)
            for field in PAIR_FIELDS:
                self.assertAlmostEqual(float(row[f"change_{field}"]), float(null[field])-float(wt[field]), places=13)

    def test_matched_secretion_comparisons_and_saved_five_percent_reference(self):
        traces = {(r["root_id"], float(r["ae4_expression"]), float(r["calcium_uM"])): r for r in self.traces}
        old = {(r["root_id"], float(r["calcium_uM"])): r for r in read_rows(TASK24/"comparison.csv")}
        for row in self.pairs:
            root, ca = row["root_id"], float(row["calcium_uM"])
            wt, null = traces[root, 1.0, ca], traces[root, 0.0, ca]
            self.assertLessEqual(abs(float(wt["total_0_600_pL"])-float(old[root,ca]["wt_total_0_600_pL"]))/float(old[root,ca]["wt_total_0_600_pL"]),1e-6)
            if row["numerical_pair_available"] != "True":
                self.assertEqual(null["run_status"],"not_run_no_closed_null_rest")
                self.assertEqual(row["null_to_wt_ratio"],"")
                self.assertEqual(row["paired_gate_pass"],"False")
                continue
            ratio = float(null["total_0_600_pL"])/float(wt["total_0_600_pL"])
            self.assertAlmostEqual(float(row["null_to_wt_ratio"]), ratio, places=14)
            self.assertAlmostEqual(float(row["change_percent"]), 100*(ratio-1), places=12)
            self.assertEqual(float(row["saved_task24_5pct_to_wt_ratio"]), float(old[root,ca]["ae4_5pct_to_wt_ratio"]))
            self.assertLessEqual(float(row["wt_relative_difference_from_saved_task24"]), 1e-6)
            expected_gate = wt["matched_rest_admissible"] == "True" and null["matched_rest_admissible"] == "True"
            self.assertEqual(row["paired_gate_pass"] == "True", expected_gate)

    def test_only_declared_utf8_text_outputs_exist(self):
        self.assertEqual({path.name for path in OUT.iterdir()}, set(TEXT_OUTPUTS))
        for path in OUT.iterdir():
            self.assertTrue(path.is_file())
            content = path.read_bytes()
            self.assertNotIn(b"\x00", content)
            content.decode("utf-8")
            if path.suffix == ".json":
                json.loads(content, parse_constant=lambda value: self.fail(f"Nonfinite JSON value: {value}"))


if __name__ == "__main__":
    unittest.main()
