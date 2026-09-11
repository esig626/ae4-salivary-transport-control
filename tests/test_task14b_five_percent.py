"""Focused Task 14B checks. Tests never run an ODE or solve a resting state."""
from dataclasses import asdict, replace
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import task14b_five_percent as task
from modern_full_model.genotype_evaluation import GenotypeTrajectory
from modern_full_model.validation import PRODUCTION_RADAU, sha256_object


class FivePercentRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.wt, cls.ae2, cls.frozen = task._preflight()
        cls.root = sorted(cls.manifest["roots"])[0]

    def test_frozen_complete_panel_and_inputs(self):
        expected = {(root, calcium) for root in self.manifest["roots"]
                    for calcium in (0.10, 0.25, 0.50)}
        self.assertEqual(len(expected), 30)
        self.assertEqual(set(self.wt), expected)
        self.assertEqual(set(self.ae2), expected)
        self.assertEqual({(r["root_id"], r["calcium_uM"]) for r in self.frozen["requests"]}, expected)
        self.assertTrue(all(value > 0 for value in self.wt.values()))
        self.assertEqual(self.frozen["new_resting_solves"], 0)
        self.assertEqual(self.frozen["new_wt_trajectories"], 0)

    def test_saved_selected_states_are_used_exactly(self):
        for root in self.manifest["roots"]:
            data = json.loads(task.continuation_path(root, "AE4").read_text())
            step = next(s for s in data["steps"] if abs(s["expression"] - .05) < 1e-10)
            selected = next(r for r in step["roots"] if r["root_id"] == step["selected_root_id"])
            self.assertEqual(task._five_percent_core(root), tuple(selected["core_state"]))
            self.assertTrue(step["connected_to_previous"])
            self.assertEqual(selected["normalized_jacobian_rank"], 10)

    def test_bad_wt_hash_stops_before_any_simulation(self):
        wt_path = task.trajectory_path(self.root, .10, "WT", "production_radau")
        original = task.sha256_file
        with patch.object(task, "sha256_file", side_effect=lambda p:
                          "0" * 64 if Path(p) == wt_path else original(p)), \
                patch.object(task, "simulate_genotype") as simulate:
            with self.assertRaisesRegex(AssertionError, "Frozen Task 14 artifact changed"):
                task._preflight()
            simulate.assert_not_called()

    def test_changed_wt_summary_does_not_replace_denominator(self):
        contract, _ = task.verify_contract()
        wt_path = task.trajectory_path(self.root, .10, "WT", "production_radau")
        original = Path.read_text
        def altered(path, *args, **kwargs):
            content = original(path, *args, **kwargs)
            if path == wt_path:
                value = json.loads(content)
                value["trajectory"]["cumulative_flow_pL"][-1] *= 1.01
                return json.dumps(value)
            return content
        with patch.object(Path, "read_text", altered), patch.object(task, "simulate_genotype") as simulate:
            with self.assertRaisesRegex(AssertionError, "trajectory/summary mismatch"):
                task._wt_totals(contract, self.manifest)
            simulate.assert_not_called()

    def test_all_models_keep_parameters_regulation_and_calcium(self):
        for root, calcium in self.wt:
            model = task.build_model(self.manifest, root, calcium)
            payload = self.manifest["roots"][root]
            self.assertEqual(sha256_object(model.parameters), payload["whole_cell_parameters_sha256"])
            self.assertEqual(sha256_object(model.ae4_parameters), payload["ae4_parameters_sha256"])
            regulation = model.regulatory_model
            self.assertEqual(sha256_object(regulation.ae4_regulatory_model),
                             sha256_object(self.frozen["regulatory_member"]["regulatory_model"]))
            self.assertEqual(regulation.nkcc1_regulatory_model.stimulated_calcium_uM, .10)
            self.assertEqual(regulation.nkcc1_regulatory_model.fully_activated_multiplier, 1.75)
            self.assertEqual(regulation.stimulus.stimulated_calcium_uM, calcium)
            self.assertEqual(regulation.stimulus.duration_s, 600)
            self.assertEqual(regulation.stimulus.cch_dose_uM, .3)
            self.assertEqual(regulation.stimulus.ipr_dose_uM, 5.)

    def _replayed_trace(self):
        """Saved WT values are only a mock return to inspect worker call arguments."""
        data = json.loads(task.trajectory_path(self.root, .10, "WT", "production_radau").read_text())
        values = data["trajectory"]
        values["solver"] = PRODUCTION_RADAU
        for key, value in values.items():
            if isinstance(value, list):
                values[key] = tuple(value) if key == "state_names" else np.asarray(value)
        return GenotypeTrajectory(**values)

    def test_worker_uses_five_percent_and_original_wt_grid(self):
        with TemporaryDirectory(dir=task.REPO, prefix="task14b_test_") as temporary, \
                patch.object(task, "OUT", Path(temporary)), \
                patch.object(task, "simulate_genotype", return_value=self._replayed_trace()) as simulate:
            row = task._run_one((self.root, .10))
            args, kwargs = simulate.call_args
            self.assertEqual(args[1], task._five_percent_core(self.root))
            self.assertEqual(kwargs["transporter"], "AE4")
            self.assertEqual(kwargs["genotype"].ae4_expression, .05)
            self.assertEqual(kwargs["genotype"].ae2_expression, 1.)
            self.assertIs(kwargs["solver"], PRODUCTION_RADAU)
            self.assertEqual(kwargs["time_s"], self.frozen["time_grid_s"])
            simulate.assert_called_once()
            saved = json.loads(task._output_path(self.root, .10).read_text())
            self.assertEqual(saved["genotype"]["ae4_expression"], .05)
            self.assertEqual(saved["trajectory"]["solver"], asdict(PRODUCTION_RADAU))
            self.assertEqual(saved["trajectory"]["time_s"], self.frozen["time_grid_s"])
            self.assertTrue(row["numerical_gate_pass"])

    def test_failed_trajectory_is_saved_and_kept_without_primary_ratio(self):
        trace = replace(self._replayed_trace(), success=False, numerical_gate_pass=False,
                        message="Mock incomplete solver return")
        with TemporaryDirectory(dir=task.REPO, prefix="task14b_test_") as temporary, \
                patch.object(task, "OUT", Path(temporary)), \
                patch.object(task, "simulate_genotype", return_value=trace):
            row = task._pair(task._run_one((self.root, .10)), 1., 1.001)
            self.assertTrue(task._output_path(self.root, .10).exists())
            self.assertEqual(row["status"], "NUMERICAL_FAILURE")
            self.assertIsNone(row["ae4_5pct_to_wt_total_ratio"])
            summary = task._summarise([row])
            self.assertEqual(summary["trajectory_count"], 1)
            self.assertEqual(summary["valid_trajectory_count"], 0)
            self.assertEqual(len(summary["failed_cases"]), 1)


@unittest.skipUnless((task.OUT / "summary.json").exists(), "Production artifacts not yet generated")
class FivePercentArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.wt, cls.ae2, cls.frozen = task._preflight()
        with (task.OUT / "five_percent_ae4_flow.csv").open(newline="") as handle:
            cls.rows = list(csv.DictReader(handle))
        cls.summary = json.loads((task.OUT / "summary.json").read_text())

    def test_exactly_thirty_unique_saved_cases(self):
        keys = [(r["root_id"], float(r["calcium_uM"])) for r in self.rows]
        self.assertEqual(len(keys), 30)
        self.assertEqual(set(keys), set(self.wt))
        self.assertEqual(len(list((task.OUT / "trajectories").glob("*.json"))), 30)
        self.assertEqual({float(r["ae4_expression"]) for r in self.rows}, {.05})
        self.assertEqual({r["status"] for r in self.rows}, {"COMPLETE"})

    def test_saved_trajectories_numerics_and_pair_arithmetic(self):
        for row in self.rows:
            key = row["root_id"], float(row["calcium_uM"])
            data = json.loads(task._output_path(*key).read_text())
            trace = data["trajectory"]
            self.assertEqual(data["genotype"]["ae4_expression"], .05)
            self.assertEqual(data["genotype"]["ae2_expression"], 1.)
            self.assertEqual(trace["solver"], asdict(PRODUCTION_RADAU))
            self.assertEqual(trace["time_s"], self.frozen["time_grid_s"])
            self.assertEqual(tuple(np.asarray(trace["states"])[:12, 0]), task._five_percent_core(key[0]))
            for field in ("whole_cell_parameters_sha256", "ae4_parameters_sha256"):
                self.assertEqual(data[field], self.manifest["roots"][key[0]][field])
            self.assertTrue(trace["success"] and trace["numerical_gate_pass"])
            self.assertTrue(np.all(np.asarray(trace["states"])[:12] > 0))
            self.assertTrue(np.all(np.asarray(trace["flow_pL_s"]) >= 0))
            for field, tolerance in self.frozen["conservation_tolerances"].items():
                self.assertLessEqual(trace["max_abs_conservation_residuals"][field], tolerance)
            total = float(np.trapezoid(trace["flow_pL_s"], trace["time_s"]))
            self.assertAlmostEqual(total, trace["cumulative_flow_pL"][-1], places=13)
            self.assertAlmostEqual(total, float(row["ae4_5pct_total_0_600_pL"]), places=13)
            self.assertEqual(float(row["wt_total_0_600_pL"]), self.wt[key])
            ratio = total / self.wt[key]
            self.assertAlmostEqual(ratio, float(row["ae4_5pct_to_wt_total_ratio"]), places=13)
            self.assertAlmostEqual(100 * (1 - ratio), float(row["ae4_5pct_reduction_percent"]), places=11)
            self.assertEqual(float(row["ae2_to_wt_total_ratio"]), self.ae2[key])
            self.assertAlmostEqual(float(row["ae2_reduction_percent"]), 100 * (1 - self.ae2[key]), places=13)

    def test_ensemble_summaries_use_every_case_without_selection(self):
        groups = [(self.rows, self.summary["overall"])]
        groups += [([r for r in self.rows if float(r["calcium_uM"]) == item["calcium_uM"]], item)
                   for item in self.summary["per_calcium"]]
        groups += [([r for r in self.rows if r["routing_family"] == item["routing_family"]], item)
                   for item in self.summary["per_routing_family"]]
        groups += [([r for r in self.rows if r["routing_family"] == item["routing_family"]
                     and float(r["calcium_uM"]) == item["calcium_uM"]], item)
                   for item in self.summary["per_calcium_and_routing_family"]]
        for group, item in groups:
            self.assertEqual(item["pair_count"], len(group))
            self.assertEqual(item["valid_count"], len(group))
            for name, field in (("ratio", "ae4_5pct_to_wt_total_ratio"),
                                ("reduction_percent", "ae4_5pct_reduction_percent")):
                values = np.asarray([float(r[field]) for r in group])
                self.assertEqual(item[name + "_min"], float(values.min()))
                self.assertEqual(item[name + "_median"], float(np.median(values)))
                self.assertEqual(item[name + "_max"], float(values.max()))

    def test_postrun_frozen_inputs_and_saved_state_hashes(self):
        saved = json.loads((task.OUT / "frozen_inputs.json").read_text())
        for key in self.frozen:
            self.assertEqual(sha256_object(saved[key]), sha256_object(self.frozen[key]), key)
        post = json.loads((task.OUT / "postrun_integrity.json").read_text())
        self.assertTrue(post["all_frozen_inputs_unchanged"])
        self.assertEqual(post["new_trajectory_count"], 30)
        self.assertEqual(post["new_resting_solves"], 0)
        self.assertEqual(post["new_wt_trajectories"], 0)
        self.assertEqual(post["new_ae2_trajectories"], 0)


if __name__ == "__main__":
    unittest.main()
