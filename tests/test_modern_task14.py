"""Focused Task 14 integrity and saved result tests; no new trajectories."""
from __future__ import annotations

import copy
from dataclasses import asdict
import inspect
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import task14_blind as b, task14_report as report, task14_compare as compare
from modern_full_model.validation import sha256_file, sha256_object, attach_basal_regulation


class TestTask14(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c, cls.digest = b.contract()
        cls.manifest, _ = b.load_freeze()
        cls.root = cls.c["root_ids"][0]

    def test_scientific_commit_tree_and_equations_are_frozen(self):
        self.assertEqual(self.c["frozen_scientific_commit"], "2f54e7c4f87b6da746d3a8427c7bdca40a77c3de")
        self.assertEqual(b.git("rev-parse", b.BASE + "^{tree}"), b.BASE_TREE)
        b.verify_contract()

    def test_all_ten_root_payloads_are_unchanged(self):
        self.assertEqual(len(self.c["root_ids"]), 10)
        self.assertEqual(self.c["root_ids"], sorted(self.manifest["roots"]))
        for root, payload in self.manifest["roots"].items():
            for key in ("native_root_object", "core_state", "whole_cell_parameters", "ae4_parameters"):
                self.assertEqual(sha256_object(payload[key]), self.c["root_payload_hashes"][root][key + "_sha256"])

    def test_exact_calcium_panel_and_R1_member(self):
        self.assertEqual(self.c["calcium_uM"], [0.10, 0.25, 0.50])
        self.assertEqual(self.c["regulatory_member_id"], "R1_G125_P21_PROSE_TREFERENCE")
        self.assertEqual(self.c["regulatory_member"]["family"], "R1")
        for ca in (0.15, 0.20, 0.30):
            with self.assertRaises(ValueError): b.build_model(self.manifest, self.root, ca)

    def test_genotypes_only_change_intended_expression(self):
        wild = asdict(b.WT)
        for root in self.c["root_ids"]:
            payload = self.manifest["roots"][root]
            for ca in b.CALCIUM:
                model = b.build_model(self.manifest, root, ca)
                self.assertEqual(sha256_object(model.parameters), payload["whole_cell_parameters_sha256"])
                self.assertEqual(sha256_object(model.ae4_parameters), payload["ae4_parameters_sha256"])
                for g in ("AE4", "AE2"):
                    obj = asdict(b.genotype_with_expression(g, 0))
                    self.assertEqual({k for k in wild if wild[k] != obj[k]}, {"name", g.lower() + "_expression"})
                self.assertEqual(model.regulatory_model.nkcc1_regulatory_model.stimulated_calcium_uM, 0.10)

    def test_exact_deletion_sources_are_zero_without_parameters_changing(self):
        for root in self.c["root_ids"]:
            model = b.build_model(self.manifest, root, 0.25)
            state = attach_basal_regulation(model, self.manifest["roots"][root]["core_state"])
            before = sha256_object(model.parameters)
            for g in ("AE4", "AE2"):
                ev = model.evaluate(1.0, state, genotype=b.genotype_with_expression(g, 0))
                if g == "AE4":
                    for key in ("na_cell_fmol_s", "k_cell_fmol_s", "cl_cell_fmol_s", "hco3_cell_fmol_s", "transported_charge_fmol_s", "tic_cell_fmol_s", "alkalinity_cell_fmol_s"):
                        self.assertEqual(getattr(ev.diagnostics.ae4, key), 0.0)
                else:
                    self.assertEqual(ev.diagnostics.homeostasis.ae2_inward_fmol_s, 0.0)
            self.assertEqual(before, sha256_object(model.parameters))

    def test_continuation_fits_zero_parameters_and_preserves_all_attempts(self):
        for root in self.c["root_ids"]:
            for g in ("AE4", "AE2"):
                d = json.loads(b.continuation_path(root, g).read_text())
                self.assertEqual(d["fitted_parameter_count"], 0)
                self.assertEqual(d["expression_grid"], self.c["continuation"]["expression_grid"])
                self.assertEqual(d["fixed_other_impermeant_osmoles_fmol"], self.manifest["roots"][root]["native_root_object"]["other_impermeant_osmoles_fmol"])
                self.assertTrue(all(len(s["attempts"]) == 3 for s in d["steps"]))
                if not d["completed_to_exact_zero"]:
                    self.assertIsNone(b.null_core(root, g))
                    self.assertIn("independent_audit", d)

    def test_matched_pairs_share_solver_and_exact_time_grid(self):
        for root in self.c["root_ids"]:
            for ca in b.CALCIUM:
                wt = b.load_trajectory(root, ca, "WT")["trajectory"]
                for g in ("AE4", "AE2"):
                    d = b.load_trajectory(root, ca, g)
                    if d["trajectory"] is None:
                        self.assertIsNone(b.null_core(root, g)); continue
                    self.assertEqual(wt["solver"], d["trajectory"]["solver"])
                    self.assertEqual(wt["time_s"], self.c["time_grid_s"])
                    self.assertEqual(wt["time_s"], d["trajectory"]["time_s"])

    def test_common_scale_cancels_in_all_ratio_metrics(self):
        wt = b.load_trajectory(self.root, 0.25, "WT")["trajectory"]
        null = b.load_trajectory(self.root, 0.25, "AE2")["trajectory"]
        base = b.ratios(wt, null)
        for scale in (1e-5, 7.2, 1e8):
            copies = copy.deepcopy([wt, null])
            for t in copies:
                for key in ("flow_pL_s", "cumulative_flow_pL"):
                    t[key] = (np.asarray(t[key]) * scale).tolist()
            scaled = b.ratios(*copies)
            for key in base: self.assertAlmostEqual(base[key], scaled[key], places=12)

    def test_ratio_arithmetic_and_minute_bins_match_direct_sums(self):
        wt = b.load_trajectory(self.root, 0.5, "WT")["trajectory"]
        null = b.load_trajectory(self.root, 0.5, "AE2")["trajectory"]
        row = b.ratios(wt, null)
        sums = []
        for t in (wt, null):
            intervals = np.diff(t["time_s"]) * (np.asarray(t["flow_pL_s"])[1:] + np.asarray(t["flow_pL_s"])[:-1]) / 2
            sums.append(float(sum(intervals)))
        self.assertAlmostEqual(row["R_total"], sums[1]/sums[0], places=12)
        self.assertAlmostEqual(row["delta_emergence"], row["R_late"] - row["R_early"], places=14)
        bad = copy.deepcopy(null); bad["time_s"][1] *= 2
        with self.assertRaises(ValueError): b.ratios(wt, bad)

    def test_blind_runner_has_no_experimental_reader_or_interval(self):
        source = inspect.getsource(b)
        for token in ("heldout_targets.csv", "validation_targets.csv", "task14_compare", "target_D_total", "D_total_interval"):
            self.assertNotIn(token, source)
        self.assertNotIn("read_csv", source)

    def test_representatives_and_required_independent_audits_are_target_independent(self):
        for family, representative in self.c["representatives_lexicographic"].items():
            self.assertEqual(representative, min(r for r in self.c["root_ids"] if self.c["routing_family_by_root"][r] == family))
            for g in ("AE4", "AE2"):
                d = json.loads(b.continuation_path(representative, g).read_text())
                self.assertIn("independent_audit", d)
                for ca in b.CALCIUM:
                    self.assertTrue(b.trajectory_path(representative, ca, g, "production_bdf").is_file())

    def test_reveal_writer_cannot_overwrite_blind_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); blind = root / "blind_pair_summary.csv"; blind.write_bytes(b"immutable")
            for name in ("blind_pair_summary.csv", "../blind_pair_summary.csv", "trajectories/x.json", "prediction_contract.json"):
                with self.assertRaises(ValueError): compare.write_final(root, name, b"overwrite")
            link = root / "ae4_holdout_comparison.csv"; link.symlink_to(blind)
            with self.assertRaises(ValueError): compare.write_final(root, link.name, b"overwrite")
            self.assertEqual(blind.read_bytes(), b"immutable")
        source = inspect.getsource(compare)
        self.assertNotIn("solve_dynamics", source)
        self.assertNotIn("task14_blind", source)

    def test_all_summaries_preserve_all_roots_calcium_and_failures(self):
        expected = {(r, ca) for r in self.c["root_ids"] for ca in b.CALCIUM}
        for name in ("ae4_null_blind_predictions.csv", "ae2_null_blind_predictions.csv", "blind_pair_summary.csv"):
            rows = report.read_rows(b.OUT / name)
            self.assertEqual(len(rows), 30)
            self.assertEqual({(r["root_id"], float(r["calcium_uM"])) for r in rows}, expected)
            self.assertTrue(all(r["prediction_contract_sha256"] == self.digest for r in rows))
        for name in ("ae4_holdout_comparison.csv", "ae2_validation_comparison.csv"):
            if (b.OUT / name).exists():
                rows = report.read_rows(b.OUT / name)
                self.assertEqual({(r["root_id"], float(r["calcium_uM"])) for r in rows}, expected)

    def test_saved_prediction_and_comparison_evaluation_are_deterministic(self):
        one = report.pair_rows(); two = report.pair_rows()
        self.assertEqual(sha256_object(one), sha256_object(two))
        # Artificial arithmetic record, unrelated to experimental data.
        rows = [{"root_id": "synthetic", "calcium_uM": 0.25, "routing_family": "synthetic",
                 "D_total": 0.4}, {"root_id": "unavailable", "calcium_uM": 0.5,
                 "routing_family": "synthetic", "D_total": None}]
        target = {"D_total": 0.2, "D_total_interval": None}
        left = compare.evaluate_rows(rows, target); right = compare.evaluate_rows(rows, target)
        self.assertEqual(compare.render_csv(left), compare.render_csv(right))
        self.assertAlmostEqual(left[0]["signed_error_D_total"], 0.2)
        self.assertIsNone(left[0]["inside_reported_interval"])
        self.assertEqual(left[1]["comparison_status"], "UNAVAILABLE")
        self.assertIsNone(left[1]["signed_error_D_total"])

    def test_numerical_failures_are_not_fabricated_predictions(self):
        audit = json.loads((b.OUT / "blind_numerical_audit.json").read_text())
        self.assertEqual(audit["primary_matrix_rows"], 90)
        self.assertEqual(audit["primary_trajectories_run"] + audit["primary_continuation_blocked"], 90)
        self.assertTrue(audit["all_available_solver_confirmations_pass"])
        self.assertLessEqual(audit["max_available_dynamic_conservation_ratio"], 1)
        self.assertTrue(audit["independent_all_agree_on_completion_and_shared_steps"])
        self.assertLessEqual(audit["max_dynamic_charge_fmol"], self.c["continuation"]["numerical_tolerances"]["charge_tolerance_fmol"])


if __name__ == "__main__": unittest.main()
