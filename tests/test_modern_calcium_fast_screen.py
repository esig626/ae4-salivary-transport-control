"""Task 13C freeze and execution contracts; no new ODE integrations."""

from __future__ import annotations

import copy
import csv
from dataclasses import asdict
import inspect
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import run_calcium_fast_screen as screen
from modern_full_model import report_calcium_fast_screen as report


def replay_worker(job):
    """Real saved solver output through the unchanged complete WT worker."""
    manifest, request = job
    payload = json.loads(screen.cache_path(request).read_text())
    trajectory = screen.decode_trajectory(payload["trajectory"])
    with patch.object(screen, "simulate_wt", return_value=trajectory):
        result = screen.simulate_request(job)
    return screen.scientific_payload(result)


def reduction_worker(job):
    _, request = job
    payload = json.loads(screen.cache_path(request).read_text())
    return screen.summarise(screen.decode_trajectory(payload["trajectory"]))


class TestCalciumFreeze(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, cls.verification = screen.load_freeze()

    def test_all_ten_frozen_roots_are_hash_identical(self):
        self.assertEqual(len(self.manifest["roots"]), 10)
        for payload in self.manifest["roots"].values():
            for key in ("native_root_object", "core_state", "whole_cell_parameters", "ae4_parameters"):
                self.assertEqual(screen.sha256_object(payload[key]), payload[key + "_sha256"])
        self.assertEqual(screen.sha256_file(screen.INHERITED / "native_dynamic_contract_manifest.json"),
                         screen.MANIFEST_SHA256)

    def test_all_non_calcium_parameters_and_nkcc_reference_are_unchanged(self):
        for root_id, payload in self.manifest["roots"].items():
            hashes = []
            for ca in screen.NEW_CALCIUM:
                model = screen.build_model(self.manifest, screen.Request(root_id, ca))
                self.assertEqual(screen.sha256_object(model.parameters), payload["whole_cell_parameters_sha256"])
                self.assertEqual(screen.sha256_object(model.ae4_parameters), payload["ae4_parameters_sha256"])
                nkcc = model.regulatory_model.nkcc1_regulatory_model
                self.assertEqual(nkcc.stimulated_calcium_uM, 0.10)
                self.assertEqual(nkcc.resting_calcium_uM, 0.058)
                self.assertEqual(nkcc.fully_activated_multiplier, 1.75)
                protocol = asdict(model.stimulus)
                protocol["arm"] = protocol["arm"].value
                protocol["stimulated_calcium_uM"] = 0.10
                self.assertEqual(protocol, self.manifest["protocols"]["CCH_IPR"])
                hashes.append(screen.frozen_parameter_hash(model))
            self.assertEqual(*hashes)

    def test_primary_matrix_is_exactly_twenty_allowed_requests(self):
        requests = screen.primary_requests(self.manifest)
        self.assertEqual(len(requests), 20)
        self.assertEqual(len(set(requests)), 20)
        self.assertEqual({r.calcium_uM for r in requests}, {0.25, 0.50})
        self.assertEqual({r.solver_label for r in requests}, {screen.LOOSE_RADAU.label})

    def test_reference_and_other_calcium_values_cannot_be_simulated(self):
        for ca in (0.10, 0.15, 0.18, 0.20, 0.30, 0.35, 0.40, 0.55):
            with self.assertRaises(ValueError):
                screen.Request("any", ca)

    def test_reference_is_loaded_without_calling_model_or_integrator(self):
        with patch.object(screen, "simulate_wt", side_effect=AssertionError("Unexpected solve")), \
             patch.object(screen, "build_model", side_effect=AssertionError("Unexpected model")):
            rows = screen.load_reference(self.manifest)
        self.assertEqual(len(rows), 10)
        self.assertEqual({r["calcium_uM"] for r in rows}, {0.10})
        self.assertEqual({r["data_origin"] for r in rows}, {"FROZEN_TASK13B_NOT_RECOMPUTED"})

    def test_only_preferred_wt_costimulation_is_exposed(self):
        model = screen.build_model(self.manifest, screen.primary_requests(self.manifest)[0])
        self.assertEqual(model.stimulus.arm, screen.StimulusArm.CCH_IPR)
        self.assertEqual(screen.preferred_member().member_id, screen.MEMBER_ID)
        self.assertEqual(screen.preferred_member().family, "R1")
        self.assertEqual(model.stimulus(600).beta_input, 1)
        self.assertEqual(model.stimulus(0).calcium_uM, 0.058)

    def test_freeze_loader_and_complete_worker_do_not_access_target_paths(self):
        original = Path.open
        opened = []

        def guarded(path, *args, **kwargs):
            name = str(path)
            opened.append(name)
            if any(x in name.lower() for x in ("heldout", "held_out", "phenotype", "reveal_log", "archive")):
                raise AssertionError("Forbidden file access")
            return original(path, *args, **kwargs)

        with patch.object(Path, "open", guarded), \
             patch.object(screen.NativeRootRecord, "__repr__", return_value="FROZEN_WT"):
            manifest, _ = screen.load_freeze()
            screen.load_reference(manifest)
            request = screen.primary_requests(manifest)[0]
            replay_worker((manifest, request))
        self.assertTrue(opened)
        source = inspect.getsource(screen)
        self.assertNotIn("genotype_evaluation", source)
        self.assertNotIn("_opaque_firewall_hashes", source)
        self.assertNotIn("solve_resting_root(", source)

    def test_solver_grid_and_units_match_the_inherited_freeze(self):
        self.assertEqual(screen.physical_time_grid(600, 5).size, 122)
        self.assertEqual(screen.physical_time_grid(600, 5)[-1], 600)
        for spec in screen.SOLVERS.values():
            self.assertIn(asdict(spec), self.manifest["production_solvers"]
                          + self.manifest["tolerance_solvers_reference_r3_only"])


class TestCalciumDecision(unittest.TestCase):
    def rows(self, flow=None):
        q = 0.5 * screen.REQUIRED_FLOW if flow is None else flow
        return [dict(root_id=f"root{i:02}", calcium_uM=ca,
                     numerical_gate_pass=True, one_smg_geometry_gate_pass=False,
                     scientific_gate_pass=False, minimum_minute_flow_pL_s=q,
                     mean_flow_0_600_pL_s=q) for ca in screen.NEW_CALCIUM for i in range(10)]

    def test_clear_failures_confirm_only_one_high_endpoint(self):
        selected, decisions = screen.production_requests(self.rows())
        self.assertEqual(selected, [screen.Request("root00", 0.50, screen.PRODUCTION_RADAU.label)])
        self.assertFalse(any(d["all_roots_triggered"] for d in decisions))

    def test_ninety_percent_boundary_is_inclusive_and_frozen(self):
        rows = self.rows()
        rows[0]["minimum_minute_flow_pL_s"] = np.nextafter(0.9 * screen.REQUIRED_FLOW, 0.0)
        self.assertEqual(len(screen.production_requests(rows)[0]), 1)
        rows[0]["minimum_minute_flow_pL_s"] = 0.9 * screen.REQUIRED_FLOW
        self.assertEqual(len(screen.production_requests(rows)[0]), 11)
        rows[10]["minimum_minute_flow_pL_s"] = 0.9 * screen.REQUIRED_FLOW
        self.assertEqual(len(screen.production_requests(rows)[0]), 20)

    def test_pass_or_inconclusive_numerics_triggers_all_roots(self):
        for key, value in (("one_smg_geometry_gate_pass", True), ("numerical_gate_pass", False)):
            rows = self.rows()
            rows[0][key] = value
            self.assertEqual(len(screen.production_requests(rows)[0]), 11)

    def test_incomplete_primary_screen_cannot_be_classified(self):
        with self.assertRaises(ValueError):
            screen.production_requests(self.rows()[:-1])

    def test_best_root_uses_binding_minute_flow_not_average(self):
        rows = self.rows()
        rows[12]["mean_flow_0_600_pL_s"] *= 2
        rows[13]["minimum_minute_flow_pL_s"] *= 1.1
        self.assertEqual(screen.production_requests(rows)[0][0].root_id, "root03")

    def test_bdf_selection_is_limited_to_decisive_cases(self):
        rows = self.rows()[10:]
        self.assertEqual(len(screen.bdf_requests(rows)), 1)
        rows[2]["scientific_gate_pass"] = True
        self.assertEqual(len(screen.bdf_requests(rows)), 2)
        for row in rows:
            row["scientific_gate_pass"] = True
        self.assertEqual(len(screen.bdf_requests(rows)), 1)

    def test_available_worker_default_respects_affinity_and_cap(self):
        with patch.object(screen.os, "sched_getaffinity", return_value={0, 1}):
            self.assertEqual(screen.available_workers(), 2)
        with patch.object(screen.os, "sched_getaffinity", return_value=set(range(32))):
            self.assertEqual(screen.available_workers(), 10)


class TestCalciumSavedExecution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, _ = screen.load_freeze()
        cls.requests = screen.primary_requests(cls.manifest)

    def test_parallel_and_serial_complete_worker_replay_agree(self):
        # Both executor modes reconstruct real frozen models and use the same
        # saved solver trajectory. This checks execution and reduction identity
        # without adding an unnecessary twenty first loose integration.
        serial = screen.run_requests(self.manifest, self.requests[:2], 1, worker=replay_worker)
        parallel = screen.run_requests(self.manifest, self.requests[:2], 2, worker=replay_worker)
        self.assertEqual(serial, parallel)

    def test_all_twenty_saved_trajectory_reductions_are_serial_parallel_identical(self):
        serial = screen.run_requests(self.manifest, self.requests, 1, worker=reduction_worker)
        parallel = screen.run_requests(self.manifest, self.requests, 2, worker=reduction_worker)
        self.assertEqual(serial, parallel)

    def test_complete_screen_has_only_allowed_inputs_and_matching_parameter_hashes(self):
        rows = []
        for request in self.requests:
            saved = json.loads(screen.cache_path(request).read_text())
            self.assertEqual(saved["request"], asdict(request))
            self.assertEqual(saved["scientific_sha256"], screen.sha256_object(screen.scientific_payload(saved)))
            row = saved["row"]
            self.assertEqual((row["genotype"], row["arm"], row["regulatory_member_id"], row["solver_label"]),
                             ("WT", "CCH_IPR", screen.MEMBER_ID, screen.LOOSE_RADAU.label))
            self.assertEqual(row["frozen_parameter_sha256"],
                screen.frozen_parameter_hash(screen.build_model(self.manifest, request)))
            rows.append(row)
        self.assertEqual(len(rows), 20)
        self.assertEqual(len({(r["root_id"], r["frozen_parameter_sha256"]) for r in rows}), 10)

    def test_saved_confirmation_exactly_matches_the_predeclared_rule(self):
        saved = json.loads((screen.OUTPUT / "confirmation_decisions.json").read_text())
        rows = [json.loads(screen.cache_path(r).read_text())["row"] for r in self.requests]
        selected, decisions = screen.production_requests(rows)
        self.assertEqual(saved["decisions"], decisions)
        self.assertEqual(saved["radau_requests"], [asdict(r) for r in selected])
        production = [json.loads(screen.cache_path(r).read_text())["row"] for r in selected]
        self.assertEqual(saved["bdf_requests"], [asdict(r) for r in screen.bdf_requests(production)])
        freeze = json.loads((screen.OUTPUT / "frozen_manifest.json").read_text())
        self.assertEqual(freeze["decision_rule"], screen.DECISION_RULE)
        self.assertEqual(freeze["decision_rule_sha256"], screen.sha256_object(screen.DECISION_RULE))

    def test_completed_requests_are_reused_without_simulation(self):
        with patch.object(screen, "simulate_request", side_effect=AssertionError("Unexpected solve")):
            result = screen.run_cached(self.manifest, self.requests, 1)
        self.assertEqual(len(result), 20)

    def test_scientific_artifacts_are_deterministic_apart_from_runtime(self):
        saved = json.loads(screen.cache_path(self.requests[0]).read_text())
        modified = copy.deepcopy(saved)
        modified["row"]["wall_seconds"] += 10
        modified["row"]["integrator_wall_seconds"] += 10
        self.assertEqual(screen.scientific_payload(saved), screen.scientific_payload(modified))
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "first.csv"
            q = Path(directory) / "second.csv"
            row = screen.scientific_payload(saved["row"])
            screen.write_rows(p, [row])
            screen.write_rows(q, [dict(reversed(list(row.items())))])
            self.assertEqual(p.read_bytes(), q.read_bytes())

    def test_artifact_inventory_has_no_extra_trajectories(self):
        decision = json.loads((screen.OUTPUT / "confirmation_decisions.json").read_text())
        requests = self.requests + [screen.Request(**r) for r in decision["radau_requests"] + decision["bdf_requests"]]
        self.assertEqual(set((screen.OUTPUT / "trajectories").glob("*.json")),
                         {screen.cache_path(r) for r in requests})

    def test_report_requires_complete_screen_and_selected_confirmations(self):
        rows = report.read_rows(screen.OUTPUT / "calcium_screen.csv")
        confirmations = report.read_rows(screen.OUTPUT / "calcium_confirmatory.csv")
        report.validate_complete_results(self.manifest, rows, confirmations)
        with self.assertRaises(ValueError):
            report.validate_complete_results(self.manifest, rows[:-1], confirmations)
        with self.assertRaises(ValueError):
            report.validate_complete_results(self.manifest, rows, confirmations[:-1])

    def test_negative_classification_requires_successful_crosscheck(self):
        rows = report.read_rows(screen.OUTPUT / "calcium_screen.csv")
        confirmations = report.read_rows(screen.OUTPUT / "calcium_confirmatory.csv")
        self.assertEqual(report.classify(rows, confirmations), report.CLASSIFICATIONS[2])
        failed = copy.deepcopy(confirmations)
        next(r for r in failed if r["solver_label"] == screen.PRODUCTION_BDF.label)["solver_crosscheck_gate_pass"] = False
        self.assertEqual(report.classify(rows, failed), report.CLASSIFICATIONS[3])

    def test_reports_regenerate_deterministically_without_a_solver(self):
        paths = list(report.ANALYSIS.glob("*.md")) + [screen.OUTPUT / "assessment.json"]
        before = {p: p.read_bytes() for p in paths}
        with patch.object(screen, "simulate_wt", side_effect=AssertionError("Unexpected solve")), redirect_stdout(io.StringIO()):
            report.main()
        self.assertEqual(before, {p: p.read_bytes() for p in paths})


if __name__ == "__main__":
    unittest.main()
