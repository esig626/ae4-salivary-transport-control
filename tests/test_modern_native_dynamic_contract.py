"""WT-only native-source dynamic-contract invariants for Task 13B."""

from __future__ import annotations

import ast
import inspect
from dataclasses import replace
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import native_dynamic_contract as contract
from modern_full_model.validation import (
    CONSERVATION_RESIDUAL_UNITS,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    Trajectory,
    pre_reveal_regulatory_ensemble,
)
from modern_full_model.run_native_dynamic_contract import (
    DYNAMIC_FREEZE_ARTIFACT_FILENAMES,
    DynamicResult,
    FROZEN_EQUATION_PATHS,
    NATIVE_FREEZE_INPUT_FILENAMES,
    WT_DYNAMIC_REPORT_FILENAMES,
    _assert_native_gate_payload,
    _assert_native_pre_reveal_payload,
    _comparison_row,
    _has_exact_key_set,
    _opaque_firewall_hashes,
    _require_relative_files,
    eligible_native_roots,
    require_native_final_ready,
)
from modern_full_model import run_native_dynamic_contract as runner


class TestAbsoluteGlandGeometry(unittest.TestCase):
    def test_one_smg_cell_count_and_scale_are_independent_upper_bounds(self) -> None:
        payload = contract.geometry_contract_payload()
        self.assertEqual(
            payload["one_smg_effective_cell_count_max"],
            contract.ONE_SMG_EFFECTIVE_CELL_COUNT_MAX,
        )
        self.assertEqual(contract.ONE_SMG_EFFECTIVE_CELL_COUNT_MAX, 51_469_231)
        self.assertAlmostEqual(
            contract.ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S,
            3088.15386,
            places=8,
        )
        self.assertAlmostEqual(
            payload["minimum_per_cell_flow_for_9_uL_min_pL_s"],
            0.0029143625635,
            places=12,
        )
        self.assertAlmostEqual(
            payload["minimum_flow_per_frozen_cell_volume_for_9_uL_min_s_inv"],
            0.0022418173565,
            places=12,
        )
        self.assertEqual(
            contract.PAIRED_SMG_EFFECTIVE_CELL_COUNT_MAX,
            2 * contract.ONE_SMG_EFFECTIVE_CELL_COUNT_MAX,
        )

    def test_geometry_clips_but_never_refits_the_shared_interval(self) -> None:
        interval = contract.SharedScaleInterval(
            lower_uL_min_per_pL_s=3000.0,
            upper_uL_min_per_pL_s=3200.0,
            sample_count=40,
            minimum_flow_pL_s=0.003,
            maximum_flow_pL_s=0.0031,
        )
        production = contract.assess_geometry_scale(interval)
        self.assertTrue(production.geometry_gate_pass)
        self.assertEqual(production.feasible_lower_uL_min_per_pL_s, 3000.0)
        self.assertAlmostEqual(
            production.feasible_upper_uL_min_per_pL_s,
            contract.ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S,
        )
        failing = contract.assess_geometry_scale(
            contract.SharedScaleInterval(3100.0, 3200.0, 40, 0.003, 0.0031)
        )
        self.assertFalse(failing.geometry_gate_pass)
        paired = contract.assess_geometry_scale(
            contract.SharedScaleInterval(3100.0, 3200.0, 40, 0.003, 0.0031),
            paired_smg_sensitivity=True,
        )
        self.assertTrue(paired.geometry_gate_pass)
        self.assertIn("SENSITIVITY", paired.convention)

    def test_malformed_geometry_intervals_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            contract.SharedScaleInterval(-2.0, -1.0, 0, -3.0, -2.0)


class TestSharedDynamicContract(unittest.TestCase):
    def test_scale_interval_is_the_intersection_over_every_flow_sample(self) -> None:
        flows = np.asarray((0.0030, 0.0032, 0.0031))
        interval = contract.shared_observation_scale_interval(flows)
        self.assertEqual(interval.sample_count, 3)
        self.assertAlmostEqual(interval.lower_uL_min_per_pL_s, 3000.0)
        self.assertAlmostEqual(interval.upper_uL_min_per_pL_s, 3125.0)
        self.assertTrue(interval.compatible)

        empty = contract.shared_observation_scale_interval((0.002, 0.004))
        self.assertFalse(empty.compatible)

    def test_sustainment_uses_fixed_sixty_and_six_hundred_second_landmarks(self) -> None:
        time = np.asarray((0.0, 60.0, 300.0, 600.0))
        flow = np.asarray((0.003, 0.003, 0.0031, 0.0027))
        ratio = contract.sustainment_ratio(time, flow)
        self.assertAlmostEqual(ratio, 0.9)
        self.assertTrue(contract.sustainment_gate(ratio))
        self.assertFalse(contract.sustainment_gate(0.79))
        self.assertFalse(contract.sustainment_gate(1.21))

    def test_costimulation_total_is_paired_with_same_row_cch_control(self) -> None:
        self.assertAlmostEqual(contract.paired_arm_total_ratio(1.0, 1.05), 1.05)
        with self.assertRaises(ValueError):
            contract.paired_arm_total_ratio(0.0, 1.0)


class TestRegulatoryFreeze(unittest.TestCase):
    def test_twenty_member_panel_has_six_exact_r2_r3_pairs(self) -> None:
        ensemble = pre_reveal_regulatory_ensemble()
        contract.validate_regulatory_ensemble(ensemble)
        self.assertEqual(len(ensemble), 20)
        pairs = contract.regulatory_member_pairs(ensemble)
        self.assertEqual(len(pairs), 6)
        registry = {member.member_id: member for member in ensemble}
        for r2_id, r3_id in pairs:
            r2 = registry[r2_id].regulatory_model
            r3 = registry[r3_id].regulatory_model
            self.assertEqual(r2.tau_camp_s, r3.tau_pka_s)
            self.assertEqual(
                r2.tau_activation_s,
                1.0 / r3.forward_regulation_rate_s,
            )
            self.assertEqual(
                r3.forward_regulation_rate_s,
                r3.reverse_regulation_rate_s,
            )

    def test_duplicate_regulatory_axes_fail_closed(self) -> None:
        ensemble = list(pre_reveal_regulatory_ensemble())
        target_index = next(
            index
            for index, member in enumerate(ensemble)
            if member.family == "R2"
            and member.gain_label == "G170_P21_FIGURE_VISUAL"
            and member.timing_label == "TSLOW"
        )
        duplicate = replace(
            ensemble[target_index],
            member_id="DUPLICATE_AXIS_WITH_NEW_ID",
            gain_label="G125_P21_PROSE",
            timing_label="TFAST",
        )
        ensemble[target_index] = duplicate
        with self.assertRaises(ValueError):
            contract.validate_regulatory_ensemble(ensemble)

    def test_helper_source_is_strictly_wt_only(self) -> None:
        for module in (contract, runner):
            source = inspect.getsource(module)
            self.assertNotIn("AE4_NULL", source)
            self.assertNotIn("heldout_targets.csv", source)
            self.assertNotIn("genotype=AE4", source)

    def test_dynamic_intake_excludes_rejected_and_diagnostic_roots(self) -> None:
        def root(
            root_id: str,
            *,
            numerical: bool = True,
            wt: bool = True,
            eligibility: str = "PRODUCTION_PRE_REVEAL",
            hydraulic: float = 1.0,
            hydraulic_mode: str = "H1",
        ) -> SimpleNamespace:
            return SimpleNamespace(
                root_id=root_id,
                passes_numerical_gate=numerical,
                passes_wt_gate=wt,
                eligibility=eligibility,
                hydraulic_scale=hydraulic,
                hydraulic_mode=hydraulic_mode,
            )

        retained = eligible_native_roots(
            (
                root("keep"),
                root("numerical_fail", numerical=False),
                root("wt_fail", wt=False),
                root("hydraulic_diagnostic", hydraulic=2.5, hydraulic_mode="HW"),
                root("hwq_diagnostic", hydraulic=1.0, hydraulic_mode="HWQ"),
                root("stress", eligibility="DIAGNOSTIC_ONLY_STRESS"),
                root("production_rejected", eligibility="PRODUCTION_REJECTED"),
                SimpleNamespace(
                    root_id="string_false",
                    passes_numerical_gate="False",
                    passes_wt_gate="True",
                    eligibility="PRODUCTION_PRE_REVEAL",
                    hydraulic_scale=1.0,
                    hydraulic_mode="H1",
                ),
            )
        )
        self.assertEqual([item.root_id for item in retained], ["keep"])

        with self.assertRaises(ValueError):
            eligible_native_roots(
                (
                    SimpleNamespace(
                        root_id="malformed_bool",
                        passes_numerical_gate="yes",
                        passes_wt_gate=True,
                        eligibility="PRODUCTION_PRE_REVEAL",
                        hydraulic_scale=1.0,
                        hydraulic_mode="H1",
                    ),
                )
            )

    def test_native_summary_must_be_explicitly_post_geometry_ready(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "native_source_wt_summary.json"
            path.write_text(json.dumps({"production_confirmation_complete": True}))
            with self.assertRaises(RuntimeError):
                require_native_final_ready(root)
            ready = {
                "production_confirmation_complete": True,
                "geometry_confirmation_complete": True,
                "independent_root_reproduction_complete": True,
                "native_pre_reveal_freeze_ready": True,
            }
            path.write_text(json.dumps(ready))
            self.assertEqual(require_native_final_ready(root), ready)

    def test_freeze_dependencies_and_opaque_hashes_use_the_exact_passed_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(FileNotFoundError):
                _require_relative_files(root, ("missing.csv",))
            required = (
                "calibration_targets.csv",
                "validation_targets.csv",
                "heldout_" + "targets.csv",
                "reveal_log.json",
            )
            for index, name in enumerate(required):
                (root / name).write_bytes(f"opaque-{index}".encode())
            hashes = _opaque_firewall_hashes(root)
            self.assertEqual(len(hashes), 4)
            self.assertTrue(all(len(value) == 64 for value in hashes.values()))

    def test_freeze_registries_are_complete_explicit_and_present(self) -> None:
        self.assertEqual(
            DYNAMIC_FREEZE_ARTIFACT_FILENAMES,
            (
                *WT_DYNAMIC_REPORT_FILENAMES,
                "native_dynamic_contract_manifest.json",
                "native_dynamic_contract_gate.json",
            ),
        )
        self.assertEqual(len(DYNAMIC_FREEZE_ARTIFACT_FILENAMES), 7)
        self.assertEqual(
            len(set(DYNAMIC_FREEZE_ARTIFACT_FILENAMES)),
            len(DYNAMIC_FREEZE_ARTIFACT_FILENAMES),
        )
        self.assertEqual(len(NATIVE_FREEZE_INPUT_FILENAMES), 10)
        required_direct_dependencies = {
            "src/modern_full_model/calibration.py",
            "src/modern_full_model/independent.py",
            "src/modern_full_model/native_source_panel.py",
            "src/modern_full_model/independent_native.py",
        }
        self.assertTrue(required_direct_dependencies <= set(FROZEN_EQUATION_PATHS))
        repository_root = Path(__file__).resolve().parents[1]
        self.assertTrue(
            all((repository_root / path).is_file() for path in FROZEN_EQUATION_PATHS)
        )
        frozen = set(FROZEN_EQUATION_PATHS)
        for relative_path in FROZEN_EQUATION_PATHS:
            source_path = repository_root / relative_path
            tree = ast.parse(source_path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, ast.ImportFrom) or node.level != 1:
                    continue
                if not node.module:
                    continue
                dependency = f"src/modern_full_model/{node.module.split('.')[0]}.py"
                if (repository_root / dependency).is_file():
                    self.assertIn(
                        dependency,
                        frozen,
                        f"{relative_path} dependency is not frozen",
                    )

    def test_native_firewall_accepts_only_opaque_hash_bindings(self) -> None:
        payload = {
            "opaque_firewall_ledger_sha256": {
                "calibration_ledger_sha256": "a" * 64,
                "validation_ledger_sha256": "b" * 64,
                "sealed_phenotype_target_ledger_sha256": "c" * 64,
                "pre_reveal_log_sha256": "d" * 64,
            },
            "clean_auditor_authorization_required": True,
            "phenotype_reveal_permitted_by_this_manifest": False,
        }
        _assert_native_pre_reveal_payload(payload)
        with self.assertRaises(ValueError):
            _assert_native_pre_reveal_payload({**payload, "target_value": 1.0})
        gate = {
            "clean_auditor_authorization_required": True,
            "independent_numerical_reproduction_required": True,
            "phenotype_reveal_permitted": False,
            "final_composite_authorization_complete": False,
        }
        _assert_native_gate_payload(gate)
        with self.assertRaises(ValueError):
            _assert_native_gate_payload(
                {**gate, "phenotype_reveal_permitted": True}
            )

    def test_dependency_snapshot_detects_in_flight_source_or_input_change(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            repository = base / "repository"
            results = base / "results"
            repository.mkdir()
            results.mkdir()
            (repository / "source.py").write_text("source-v1\n", encoding="utf-8")
            (results / "input.csv").write_text("input-v1\n", encoding="utf-8")
            opaque = {"opaque": "0" * 64}
            with (
                patch.object(runner, "FROZEN_EQUATION_PATHS", ("source.py",)),
                patch.object(runner, "NATIVE_FREEZE_INPUT_FILENAMES", ("input.csv",)),
                patch.object(runner, "_opaque_firewall_hashes", return_value=opaque),
            ):
                snapshot = runner._capture_freeze_dependency_snapshot(
                    repository, results
                )
                runner._assert_freeze_snapshot_current(
                    snapshot, repository, results
                )
                (results / "input.csv").write_text("input-v2\n", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "changed during execution"):
                    runner._assert_freeze_snapshot_current(
                        snapshot, repository, results
                    )
                (results / "input.csv").write_text("input-v1\n", encoding="utf-8")
                (repository / "source.py").write_text("source-v2\n", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "changed during execution"):
                    runner._assert_freeze_snapshot_current(
                        snapshot, repository, results
                    )

    def test_staged_publication_keeps_hash_ledger_absent_on_commit_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            staging = base / "staging"
            results = base / "results"
            repository = base / "repository"
            staging.mkdir()
            results.mkdir()
            repository.mkdir()
            for filename in DYNAMIC_FREEZE_ARTIFACT_FILENAMES:
                (staging / filename).write_text(
                    f"staged {filename}\n", encoding="utf-8"
                )
            (staging / runner.HASH_FILENAME).write_text(
                "staged ledger\n", encoding="utf-8"
            )
            (results / runner.HASH_FILENAME).write_text(
                "stale ledger\n", encoding="utf-8"
            )
            snapshot = runner.FreezeDependencySnapshot({}, {}, {})
            with patch.object(
                runner,
                "_assert_freeze_snapshot_current",
                side_effect=(None, RuntimeError("changed during commit")),
            ):
                with self.assertRaisesRegex(RuntimeError, "changed during commit"):
                    runner._publish_staged_freeze(
                        staging, results, repository, snapshot
                    )
            self.assertFalse((results / runner.HASH_FILENAME).exists())
            self.assertTrue(
                all((results / name).is_file() for name in DYNAMIC_FREEZE_ARTIFACT_FILENAMES)
            )

    def test_staged_publication_commits_complete_ledger_last(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            staging = base / "staging"
            results = base / "results"
            repository = base / "repository"
            staging.mkdir()
            results.mkdir()
            repository.mkdir()
            for filename in DYNAMIC_FREEZE_ARTIFACT_FILENAMES:
                (staging / filename).write_text(
                    f"staged {filename}\n", encoding="utf-8"
                )
            (staging / runner.HASH_FILENAME).write_text(
                "staged ledger\n", encoding="utf-8"
            )
            snapshot = runner.FreezeDependencySnapshot({}, {}, {})
            with patch.object(
                runner,
                "_assert_freeze_snapshot_current",
                return_value=None,
            ) as current:
                runner._publish_staged_freeze(
                    staging, results, repository, snapshot
                )
            self.assertEqual(current.call_count, 3)
            self.assertEqual(
                (results / runner.HASH_FILENAME).read_text(encoding="utf-8"),
                "staged ledger\n",
            )
            self.assertTrue(
                all((results / name).is_file() for name in DYNAMIC_FREEZE_ARTIFACT_FILENAMES)
            )

    def test_runner_invalidates_prior_ledger_before_first_numerical_batch(self) -> None:
        source = inspect.getsource(runner.run_native_dynamic_contract)
        self.assertLess(
            source.index("_invalidate_hash_ledger"),
            source.index("_run_requests"),
        )


class TestSolverContract(unittest.TestCase):
    def test_complete_but_wrong_solver_or_isomorphism_keys_fail_closed(self) -> None:
        fields = ("family", "arm")
        expected = {("R2", "CCH_ONLY"), ("R3", "CCH_IPR")}
        valid = [
            {"family": "R2", "arm": "CCH_ONLY"},
            {"family": "R3", "arm": "CCH_IPR"},
        ]
        self.assertTrue(_has_exact_key_set(valid, fields, expected))
        wrong_same_count = [
            {"family": "R2", "arm": "CCH_ONLY"},
            {"family": "R1", "arm": "CCH_IPR"},
        ]
        self.assertFalse(_has_exact_key_set(wrong_same_count, fields, expected))
        duplicate = [valid[0], valid[0]]
        self.assertFalse(_has_exact_key_set(duplicate, fields, expected))
        self.assertFalse(_has_exact_key_set([{"family": "R2"}], fields, expected))

    def test_flow_time_course_disagreement_fails_even_when_states_and_total_match(self) -> None:
        grid = np.asarray((0.0, 300.0, 600.0))

        def trajectory(flow: np.ndarray, solver: object) -> Trajectory:
            return Trajectory(
                family="R3_REFERENCE",
                initial_condition="native_source_rest_root",
                solver=solver,
                success=True,
                message="ok",
                time_s=grid,
                state_names=("x",),
                states=np.ones((1, 3)),
                flow_pL_s=flow,
                cumulative_flow_pL=np.asarray((0.0, 30.0, 60.0)),
                cell_na_mM=np.ones(3),
                cell_k_mM=np.ones(3),
                cell_cl_mM=np.ones(3),
                cell_ph=np.ones(3) * 7.0,
                cell_volume_pL=np.ones(3),
                capacity_multiplier=np.ones(3),
                max_abs_conservation_residuals={
                    name: 0.0 for name in CONSERVATION_RESIDUAL_UNITS
                },
            )

        metadata = {
            "root_id": "ROOT",
            "regulatory_member_id": "R3_REFERENCE",
            "arm": "CCH_IPR",
            "initial_condition": "native_source_rest_root",
            "solver_label": PRODUCTION_RADAU.label,
            "numerical_gate_pass": True,
        }
        left = DynamicResult(
            metadata=metadata,
            trajectory=trajectory(np.asarray((0.1, 0.2, 0.1)), PRODUCTION_RADAU),
        )
        right = DynamicResult(
            metadata={**metadata, "solver_label": PRODUCTION_BDF.label},
            trajectory=trajectory(np.asarray((0.1, 0.1, 0.1)), PRODUCTION_BDF),
        )
        row = _comparison_row(left, right, comparison_kind="RADAU_VS_BDF")
        self.assertEqual(row["max_relative_state_difference"], 0.0)
        self.assertEqual(row["endpoint_cumulative_flow_relative_difference"], 0.0)
        self.assertGreater(row["max_relative_flow_difference"], 0.1)
        self.assertFalse(row["solver_gate_pass"])


if __name__ == "__main__":
    unittest.main()
