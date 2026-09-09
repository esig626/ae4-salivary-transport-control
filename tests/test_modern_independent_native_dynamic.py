"""Fail-closed independent reproduction of native WT dynamic classifications."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import run_independent_native_dynamic as runner


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _root_row(root_id: str) -> dict[str, object]:
    return {
        "root_id": root_id,
        "source_class": "N_ABS_NKCC",
        "source_scale": 2.0,
        "pump_capacity_scale": 1.0,
        "apical_pump_fraction": 0.075,
        "apical_k_fraction": 0.30,
        "ae4_cation_fraction": 0.10,
        "hydraulic_scale": 1.0,
        "hydraulic_mode": "H1",
        "eligibility": "PRODUCTION_PRE_REVEAL",
        "passes_numerical_gate": True,
        "passes_wt_gate": True,
        "core_state_json": json.dumps([float(index + 1) for index in range(12)]),
        "calibration_json": json.dumps(
            {"cell_other_impermeant_osmoles_fmol": 125.0}
        ),
    }


def _dynamic_row_fields(flow: float, *, arm: str) -> dict[str, object]:
    row: dict[str, object] = {
        "success": True,
        "positive_core": True,
        "all_flow_nonnegative": True,
        "max_dimensionless_conservation_ratio": 0.0,
        "numerical_gate_pass": True,
        "sustainment_ratio_q600_q60": 1.0,
        "sustainment_gate_pass": True,
        "capacity_multiplier_start": 1.0,
        "capacity_multiplier_endpoint": 1.0 if arm == "CCH_ONLY" else 1.25,
        "cumulative_flow_600_pL": 1.0 if arm == "CCH_ONLY" else 1.05,
        "regulatory_timing_parameters_json": json.dumps(
            {
                "tau_pka_s": 10.0,
                "forward_regulation_rate_s": 1.0 / 30.0,
                "reverse_regulation_rate_s": 1.0 / 30.0,
            }
        ),
        "ae4_fully_activated_multiplier": 1.25,
        "endpoint_na_mM": 10.0,
        "endpoint_k_mM": 130.0,
        "endpoint_cl_mM": 20.0,
        "endpoint_ph": 7.2,
        "endpoint_volume_pL": 1.3,
        "endpoint_flow_pL_s": flow,
    }
    row.update({f"flow_minute_{minute}_pL_s": flow for minute in range(1, 11)})
    return row


def _refresh_dynamic_ledger(results: Path) -> None:
    _write_csv(
        results / runner.DYNAMIC_HASH_INPUT,
        [
            {"artifact": filename, "sha256": _digest(results / filename)}
            for filename in runner.DYNAMIC_HASH_LEDGER_ARTIFACTS
        ],
    )


def _set_production_passing_ids(results: Path, root_ids: list[str]) -> None:
    manifest_path = results / runner.MANIFEST_INPUT
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["scientific_passing_group_count"] = len(root_ids)
    manifest["scientific_passing_root_ids"] = root_ids
    manifest["wt_dynamic_report_artifact_sha256"] = {
        filename: _digest(results / filename)
        for filename in runner.DYNAMIC_REPORT_INPUTS
    }
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")

    gate_path = results / runner.GATE_INPUT
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    gate["status"] = runner.PASS_GATE_STATUS if root_ids else runner.FAIL_GATE_STATUS
    gate["scientific_passing_group_count"] = len(root_ids)
    gate["scientific_passing_root_ids"] = root_ids
    gate["at_least_one_scientific_group_pass"] = bool(root_ids)
    gate_path.write_text(json.dumps(gate, sort_keys=True), encoding="utf-8")
    _refresh_dynamic_ledger(results)


def _fixture(directory: str) -> tuple[Path, Path]:
    repository = Path(directory) / "repository"
    results = repository / "results" / "13B_modern_full_model"
    results.mkdir(parents=True)
    for filename in set(
        (*runner.PRODUCTION_EQUATION_SOURCES, *runner.INDEPENDENT_DYNAMIC_SOURCES)
    ):
        path = repository / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"synthetic source: {filename}\n", encoding="utf-8")

    root_ids = ("ROOT_A", "ROOT_B")
    root_rows = [_root_row(root_id) for root_id in root_ids]
    _write_csv(results / runner.ROOT_INPUT, root_rows)
    lower = runner.WT_FLOW_LOWER_UL_MIN / 0.002
    upper = runner.WT_FLOW_UPPER_UL_MIN / 0.002
    group_rows = []
    for root_id in root_ids:
        group_rows.append(
            {
                "root_id": root_id,
                "hydraulic_scale": 1.0,
                "hydraulic_mode": "H1",
                "eligibility": "PRODUCTION_PRE_REVEAL",
                "all_native_rest_wt_gate_pass": True,
                "production_case_count": runner.EXPECTED_PROFILE_PER_ROOT,
                "nearby_case_count": runner.EXPECTED_NEARBY_PER_ROOT,
                "all_dynamic_numerical_gates_pass": True,
                "all_nearby_numerical_gates_pass": True,
                "all_sustainment_gates_pass": True,
                "all_nearby_sustainment_gates_pass": True,
                "solver_comparison_count": runner.EXPECTED_SOLVER_PER_ROOT,
                "all_solver_and_tolerance_gates_pass": True,
                "isomorphism_comparison_count": runner.EXPECTED_ISOMORPHISM_PER_ROOT,
                "r2_r3_isomorphism_gate_pass": True,
                "p21_regulatory_direction_contract_pass": True,
                "paired_arm_count": runner.EXPECTED_ARM_RATIOS_PER_ROOT,
                "costim_to_cch_total_ratio_gate_pass": True,
                "shared_scale_sample_count": runner.EXPECTED_SCALE_SAMPLES_PER_ROOT,
                "shared_scale_lower_uL_min_per_pL_s": lower,
                "shared_scale_upper_uL_min_per_pL_s": upper,
                "shared_scale_interval_nonempty": True,
                "one_smg_geometry_gate_pass": False,
                "one_smg_feasible_scale_lower_uL_min_per_pL_s": lower,
                "one_smg_feasible_scale_upper_uL_min_per_pL_s": (
                    runner.ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
                ),
                "native_dynamic_scientific_gate_pass": False,
            }
        )
    _write_csv(results / runner.GROUP_INPUT, group_rows)

    gains = ("G125_P21_PROSE", "G170_P21_FIGURE_VISUAL")
    timings = ("TFAST", "TREFERENCE", "TSLOW")
    r2_ids = [f"R2_{gain}_{timing}" for gain in gains for timing in timings]
    r3_ids = [f"R3_{gain}_{timing}" for gain in gains for timing in timings]
    member_ids = [f"R0_{gain}_STATIC" for gain in gains] + [
        f"R1_{gain}_{timing}" for gain in gains for timing in timings
    ] + r2_ids + r3_ids
    pairs = list(zip(r2_ids, r3_ids))
    member_payloads: dict[str, dict[str, object]] = {}
    for member_id in member_ids:
        family, remainder = member_id.split("_", 1)
        if family == "R0":
            gain_label = remainder.rsplit("_", 1)[0]
            timing_label = "STATIC"
        else:
            gain_label, timing_label = remainder.rsplit("_", 1)
        member_payloads[member_id] = {
            "member_id": member_id,
            "family": family,
            "gain_label": gain_label,
            "timing_label": timing_label,
            "fully_activated_multiplier": (
                1.25 if gain_label == runner.REFERENCE_GAIN_LABEL else 1.70
            ),
        }
    profile_rows: list[dict[str, object]] = []
    for root_id in root_ids:
        for member_id in member_ids:
            for arm in runner.ARMS:
                for solver in runner.METHOD_BY_SOLVER:
                    profile_rows.append(
                        {
                            "root_id": root_id,
                            "regulatory_member_id": member_id,
                            "arm": arm,
                            "solver_label": solver,
                            **_dynamic_row_fields(0.002, arm=arm),
                        }
                    )
    _write_csv(results / runner.PROFILE_INPUT, profile_rows)

    nearby_rows: list[dict[str, object]] = []
    for root_id in root_ids:
        for condition in sorted(runner.NEARBY_CONDITIONS):
            for arm in runner.ARMS:
                for solver in runner.METHOD_BY_SOLVER:
                    nearby_rows.append(
                        {
                            "root_id": root_id,
                            "initial_condition": condition,
                            "regulatory_member_id": runner.REFERENCE_MEMBER_ID,
                            "arm": arm,
                            "solver_label": solver,
                            **_dynamic_row_fields(0.002, arm=arm),
                        }
                    )
    _write_csv(results / runner.NEARBY_INPUT, nearby_rows)

    solver_rows: list[dict[str, object]] = []
    for root_id in root_ids:
        for member_id in member_ids:
            for arm in runner.ARMS:
                solver_rows.append(
                    {
                        "root_id": root_id,
                        "comparison_kind": "RADAU_VS_BDF",
                        "regulatory_member_id": member_id,
                        "arm": arm,
                        "initial_condition": "native_source_rest_root",
                        "reference_solver_label": "production_radau",
                        "comparison_solver_label": "production_bdf",
                        "both_success": True,
                        "reference_numerical_gate_pass": True,
                        "comparison_numerical_gate_pass": True,
                        "max_relative_state_difference": 0.0,
                        "max_relative_flow_difference": 0.0,
                        "endpoint_cumulative_flow_relative_difference": 0.0,
                        "solver_gate_pass": True,
                    }
                )
        for condition in sorted(runner.NEARBY_CONDITIONS):
            for arm in runner.ARMS:
                solver_rows.append(
                    {
                        "root_id": root_id,
                        "comparison_kind": "NEARBY_RADAU_VS_BDF",
                        "regulatory_member_id": runner.REFERENCE_MEMBER_ID,
                        "arm": arm,
                        "initial_condition": condition,
                        "reference_solver_label": "production_radau",
                        "comparison_solver_label": "production_bdf",
                        "both_success": True,
                        "reference_numerical_gate_pass": True,
                        "comparison_numerical_gate_pass": True,
                        "max_relative_state_difference": 0.0,
                        "max_relative_flow_difference": 0.0,
                        "endpoint_cumulative_flow_relative_difference": 0.0,
                        "solver_gate_pass": True,
                    }
                )
        for arm in runner.ARMS:
            for solver in ("loose_radau", "tight_radau"):
                solver_rows.append(
                    {
                        "root_id": root_id,
                        "comparison_kind": f"PRODUCTION_VS_{solver.upper()}",
                        "regulatory_member_id": runner.REFERENCE_MEMBER_ID,
                        "arm": arm,
                        "initial_condition": "native_source_rest_root",
                        "reference_solver_label": "production_radau",
                        "comparison_solver_label": solver,
                        "both_success": True,
                        "reference_numerical_gate_pass": True,
                        "comparison_numerical_gate_pass": True,
                        "max_relative_state_difference": 0.0,
                        "max_relative_flow_difference": 0.0,
                        "endpoint_cumulative_flow_relative_difference": 0.0,
                        "solver_gate_pass": True,
                    }
                )
    _write_csv(results / runner.SOLVER_INPUT, solver_rows)

    isomorphism_rows = [
        {
            "root_id": root_id,
            "r2_member_id": r2_id,
            "r3_member_id": r3_id,
            "arm": arm,
            "solver_label": solver,
            "max_relative_state_difference": 0.0,
            "max_relative_flow_difference": 0.0,
            "max_relative_capacity_difference": 0.0,
            "r2_r3_isomorphism_gate_pass": True,
        }
        for root_id in root_ids
        for r2_id, r3_id in pairs
        for arm in runner.ARMS
        for solver in runner.METHOD_BY_SOLVER
    ]
    _write_csv(results / runner.ISOMORPHISM_INPUT, isomorphism_rows)

    native_summary = {
        "production_confirmation_complete": True,
        "geometry_confirmation_complete": True,
        "independent_root_reproduction_complete": True,
        "native_pre_reveal_freeze_ready": True,
    }
    (results / "native_source_wt_summary.json").write_text(
        json.dumps(native_summary, sort_keys=True), encoding="utf-8"
    )
    for filename in runner.NATIVE_INPUT_ARTIFACTS:
        path = results / filename
        if not path.exists():
            path.write_text(f"synthetic native input: {filename}\n", encoding="utf-8")

    report_hashes = {
        filename: _digest(results / filename)
        for filename in runner.DYNAMIC_REPORT_INPUTS
    }
    native_hashes = {
        filename: _digest(results / filename)
        for filename in runner.NATIVE_INPUT_ARTIFACTS
    }
    equation_hashes = {
        filename: _digest(repository / filename)
        for filename in runner.PRODUCTION_EQUATION_SOURCES
    }
    roots_payload = {
        row["root_id"]: {
            "core_state_sha256": runner._sha256_object(
                tuple(json.loads(str(row["core_state_json"])))
            )
        }
        for row in root_rows
    }
    manifest = {
        "manifest_id": runner.EXPECTED_MANIFEST_ID,
        "status": runner.EXPECTED_MANIFEST_STATUS,
        "root_count": len(root_ids),
        "roots": roots_payload,
        "scientific_passing_group_count": 0,
        "scientific_passing_root_ids": [],
        "regulatory_members": member_payloads,
        "r2_r3_isomorphism_pairs": pairs,
        "clean_auditor_authorization_required": True,
        "independent_numerical_reproduction_required_for_composite_authorization": True,
        "phenotype_reveal_permitted_by_this_manifest": False,
        "wt_dynamic_report_artifact_sha256": report_hashes,
        "native_input_artifact_sha256": native_hashes,
        "equation_source_file_sha256": equation_hashes,
        "equation_source_tree_sha256": runner._sha256_object(equation_hashes),
        "opaque_firewall_ledger_sha256": {
            key: f"{index + 1:x}" * 64
            for index, key in enumerate(sorted(runner.FIREWALL_HASH_KEYS))
        },
        "opaque_firewall_access_mode": "SHA256_BYTES_ONLY_NO_TARGET_PARSE",
    }
    gate = {
        "gate_id": runner.EXPECTED_GATE_ID,
        "status": runner.FAIL_GATE_STATUS,
        "manifest": runner.MANIFEST_INPUT,
        "native_root_count": len(root_ids),
        "scientific_passing_group_count": 0,
        "scientific_passing_root_ids": [],
        "at_least_one_scientific_group_pass": False,
        "all_rows_retained_no_first_passer": True,
        "clean_auditor_authorization_required": True,
        "clean_auditor_authorization_recorded": False,
        "independent_numerical_reproduction_required": True,
        "independent_numerical_reproduction_recorded_by_this_runner": False,
        "final_composite_authorization_complete": False,
        "phenotype_reveal_permitted": False,
    }
    (results / runner.MANIFEST_INPUT).write_text(
        json.dumps(manifest, sort_keys=True), encoding="utf-8"
    )
    (results / runner.GATE_INPUT).write_text(
        json.dumps(gate, sort_keys=True), encoding="utf-8"
    )
    _refresh_dynamic_ledger(results)
    return repository, results


def _synthetic_reproduction(
    status: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    passed = status == "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_PASS"
    return (
        [{"root_id": "ROOT_A", "production_reproduction_pass": passed}],
        [{"root_id": "ROOT_A", "all_group_classification_agreement": passed}],
        {
            "status": status,
            "firewall": "PRE_REVEAL_WT_ONLY_NO_GENOTYPE_NO_HELDOUT_ACCESS",
            "root_count": 2,
            "root_ids": ["ROOT_A", "ROOT_B"],
            "case_count": 1,
            "expected_case_count": 8,
            "all_cases_reproduce_production": passed,
            "independent_reference_radau_bdf_gate_pass": passed,
            "all_group_classifications_reproduce_production": passed,
            "production_scientific_passing_root_ids": [],
            "independent_scientific_passing_root_ids": [],
            "scientific_passing_root_id_agreement": passed,
            "maximum_relative_production_difference": 0.0 if passed else 1.0,
            "maximum_independent_radau_bdf_state_difference": 0.0,
            "clean_auditor_authorization_required": True,
            "final_composite_authorization_complete": False,
            "phenotype_reveal_permitted": False,
        },
    )


class TestIndependentDynamicFrozenIntake(unittest.TestCase):
    def test_zero_passer_freeze_admits_every_root_and_reproduces_geometry_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            intake = runner.load_frozen_independent_dynamic_intake(repository)
            self.assertEqual(intake.root_ids, ("ROOT_A", "ROOT_B"))
            self.assertEqual(intake.production_passing_root_ids, ())
            self.assertEqual(len(intake.roots), 2)
            classifications, independent_passing = (
                runner._recompute_group_classifications(intake)
            )
            self.assertEqual(independent_passing, ())
            self.assertEqual(len(classifications), 2)
            self.assertTrue(
                all(row["all_group_classification_agreement"] for row in classifications)
            )
            self.assertTrue(
                all(not row["independent_geometry_gate_pass"] for row in classifications)
            )

    def test_mixed_passing_freeze_still_loads_every_group_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            group_rows = runner._read(results / runner.GROUP_INPUT)
            group_rows[0]["native_dynamic_scientific_gate_pass"] = "True"
            _write_csv(results / runner.GROUP_INPUT, group_rows)
            _set_production_passing_ids(results, ["ROOT_A"])

            intake = runner.load_frozen_independent_dynamic_intake(repository)
            self.assertEqual(intake.production_passing_root_ids, ("ROOT_A",))
            self.assertEqual(
                tuple(root["root_id"] for root in intake.roots),
                ("ROOT_A", "ROOT_B"),
            )

    def test_noncanonical_boolean_duplicate_ids_and_gate_mismatch_fail_closed(self) -> None:
        for malformed in ("true", "yes", "1", 1, None):
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                runner._csv_bool(malformed, field="gate")

        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            rows = runner._read(results / runner.GROUP_INPUT)
            rows[0]["native_dynamic_scientific_gate_pass"] = "yes"
            _write_csv(results / runner.GROUP_INPUT, rows)
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            rows = runner._read(results / runner.GROUP_INPUT)
            rows.append(dict(rows[0]))
            _write_csv(results / runner.GROUP_INPUT, rows)
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            gate = json.loads((results / runner.GATE_INPUT).read_text())
            gate["scientific_passing_root_ids"] = ["ROOT_A"]
            (results / runner.GATE_INPUT).write_text(json.dumps(gate), encoding="utf-8")
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

    def test_ledger_and_source_or_native_hash_drift_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            rows = runner._read(results / runner.DYNAMIC_HASH_INPUT)
            rows.append(dict(rows[0]))
            _write_csv(results / runner.DYNAMIC_HASH_INPUT, rows)
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            (results / "native_source_map.json").write_text(
                "stale native input\n", encoding="utf-8"
            )
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            (repository / runner.PRODUCTION_EQUATION_SOURCES[0]).write_text(
                "stale production source\n", encoding="utf-8"
            )
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

    def test_missing_isomorphism_registry_or_report_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            manifest_path = results / runner.MANIFEST_INPUT
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            del manifest["r2_r3_isomorphism_pairs"]
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            _refresh_dynamic_ledger(results)
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            (results / runner.ISOMORPHISM_INPUT).write_text(
                "root_id,r2_member_id,r3_member_id,arm,solver_label\n",
                encoding="utf-8",
            )
            manifest_path = results / runner.MANIFEST_INPUT
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["wt_dynamic_report_artifact_sha256"][
                runner.ISOMORPHISM_INPUT
            ] = _digest(results / runner.ISOMORPHISM_INPUT)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            _refresh_dynamic_ledger(results)
            with self.assertRaises(ValueError):
                runner.load_frozen_independent_dynamic_intake(repository)

    def test_primitive_metrics_override_tampered_report_gate_flags(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            intake = runner.load_frozen_independent_dynamic_intake(repository)
            profile = next(
                row for row in intake.profile_rows if row["root_id"] == "ROOT_A"
            )
            profile["max_dimensionless_conservation_ratio"] = "2.0"
            profile["sustainment_ratio_q600_q60"] = "2.0"
            solver = next(
                row for row in intake.solver_rows if row["root_id"] == "ROOT_A"
            )
            solver["max_relative_flow_difference"] = "0.5"
            isomorphism = next(
                row
                for row in intake.isomorphism_rows
                if row["root_id"] == "ROOT_A"
            )
            isomorphism["max_relative_capacity_difference"] = "0.5"

            rows, _ = runner._recompute_group_classifications(intake)
            root_a = next(row for row in rows if row["root_id"] == "ROOT_A")
            disagreements = json.loads(root_a["component_gate_disagreements_json"])
            self.assertIn("all_dynamic_numerical_gates_pass", disagreements)
            self.assertIn("all_sustainment_gates_pass", disagreements)
            self.assertIn("all_solver_and_tolerance_gates_pass", disagreements)
            self.assertIn("r2_r3_isomorphism_gate_pass", disagreements)
            self.assertFalse(root_a["primitive_report_gate_consistency"])
            self.assertFalse(root_a["all_group_classification_agreement"])

    def test_failed_trajectory_paired_count_matches_actual_count(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            intake = runner.load_frozen_independent_dynamic_intake(repository)
            failed = next(
                row
                for row in intake.profile_rows
                if row["root_id"] == "ROOT_A"
                and row["arm"] == "CCH_ONLY"
            )
            for field in (
                "success",
                "positive_core",
                "all_flow_nonnegative",
                "numerical_gate_pass",
                "sustainment_gate_pass",
            ):
                failed[field] = "False"
            failed["max_dimensionless_conservation_ratio"] = ""
            failed["sustainment_ratio_q600_q60"] = ""
            failed["cumulative_flow_600_pL"] = ""
            for minute in range(1, 11):
                failed[f"flow_minute_{minute}_pL_s"] = ""

            group = next(
                row for row in intake.group_rows if row["root_id"] == "ROOT_A"
            )
            group["all_dynamic_numerical_gates_pass"] = "False"
            group["all_sustainment_gates_pass"] = "False"
            group["p21_regulatory_direction_contract_pass"] = "False"
            group["costim_to_cch_total_ratio_gate_pass"] = "False"
            group["paired_arm_count"] = str(runner.EXPECTED_ARM_RATIOS_PER_ROOT - 1)
            group["shared_scale_sample_count"] = str(
                runner.EXPECTED_SCALE_SAMPLES_PER_ROOT - 10
            )

            rows, _ = runner._recompute_group_classifications(intake)
            root_a = next(row for row in rows if row["root_id"] == "ROOT_A")
            self.assertTrue(root_a["count_contract_agreement"])
            self.assertTrue(root_a["all_group_classification_agreement"])

    def test_independent_classification_disagreement_is_reported_not_hidden(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            intake = runner.load_frozen_independent_dynamic_intake(repository)
            intake.group_rows[0]["one_smg_geometry_gate_pass"] = "True"
            rows, _ = runner._recompute_group_classifications(intake)
            self.assertFalse(rows[0]["geometry_gate_agreement"])
            self.assertFalse(rows[0]["all_group_classification_agreement"])

    def test_geometry_classification_is_recomputed_from_case_flows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            intake = runner.load_frozen_independent_dynamic_intake(repository)
            for row in (*intake.profile_rows, *intake.nearby_rows):
                for minute in range(1, 11):
                    row[f"flow_minute_{minute}_pL_s"] = "0.003"
            rows, independent_passing = runner._recompute_group_classifications(
                intake
            )
            self.assertEqual(independent_passing, ("ROOT_A", "ROOT_B"))
            self.assertTrue(
                all(row["independent_geometry_gate_pass"] for row in rows)
            )
            self.assertTrue(all(not row["scientific_gate_agreement"] for row in rows))

    def test_reference_reproduction_rejects_negative_transient_flow(self) -> None:
        class FakeModel:
            state_names = tuple(f"state_{index}" for index in range(14))

            def __init__(self, beta: float) -> None:
                self.beta = beta

            def integrate(self, _initial: object, time: object, **_kwargs: object) -> object:
                time_s = np.asarray(time, dtype=float)
                state = np.ones((14, time_s.size), dtype=float)
                state[0] = 13.0
                state[1] = 169.0
                state[5] = 1.3
                flow = np.full(time_s.size, 0.002, dtype=float)
                flow[2] = -1.0e-6
                capacity = np.full(
                    time_s.size, 1.25 if self.beta else 1.0, dtype=float
                )
                cumulative = np.linspace(
                    0.0, 1.05 if self.beta else 1.0, time_s.size
                )
                return SimpleNamespace(
                    success=True,
                    states=state,
                    time_s=time_s,
                    flow_pL_s=flow,
                    cell_cl_mM=np.full(time_s.size, 20.0),
                    cell_ph=np.full(time_s.size, 7.2),
                    capacity_multiplier=capacity,
                    cumulative_flow_pL=cumulative,
                    max_abs_conservation_fmol_s=0.0,
                    max_abs_water_accounting_pL_s=0.0,
                )

        class FakeDefinition:
            def build_independent_model(
                self, *, stimulus: object, **_kwargs: object
            ) -> FakeModel:
                return FakeModel(float(stimulus.beta_occupancy_on))

        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            intake = runner.load_frozen_independent_dynamic_intake(repository)
            with patch.object(
                runner, "_definition", return_value=FakeDefinition()
            ):
                records, _, summary = runner.reproduce_reference_dynamics(intake)
            self.assertTrue(records)
            self.assertTrue(
                all(not row["independent_all_flow_nonnegative"] for row in records)
            )
            self.assertTrue(
                all(not row["production_reproduction_pass"] for row in records)
            )
            self.assertEqual(
                summary["status"],
                "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_DISCREPANCY",
            )


class TestIndependentDynamicPublication(unittest.TestCase):
    def test_runner_hashes_both_outputs_and_publishes_summary_last(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            publication_order: list[str] = []
            actual_write = runner._atomic_write_bytes

            def recording_write(path: Path, payload: bytes) -> None:
                publication_order.append(path.name)
                actual_write(path, payload)

            with patch.object(
                runner,
                "reproduce_reference_dynamics",
                return_value=_synthetic_reproduction(
                    "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_PASS"
                ),
            ), patch.object(
                runner, "_atomic_write_bytes", side_effect=recording_write
            ):
                summary = runner.run_independent_native_dynamic(repository)
            self.assertEqual(
                publication_order,
                [
                    runner.OUTPUT_FILENAME,
                    runner.GROUP_OUTPUT_FILENAME,
                    runner.OUTPUT_HASH_FILENAME,
                    runner.SUMMARY_FILENAME,
                ],
            )
            ledger = runner._read(results / runner.OUTPUT_HASH_FILENAME)
            self.assertEqual(
                [row["artifact"] for row in ledger],
                list(runner.OUTPUT_HASH_ARTIFACTS),
            )
            self.assertTrue(
                all(
                    row["sha256"] == _digest(results / row["artifact"])
                    for row in ledger
                )
            )
            self.assertEqual(
                summary["independent_output_hash_ledger_sha256"],
                _digest(results / runner.OUTPUT_HASH_FILENAME),
            )
            self.assertEqual(
                json.loads((results / runner.SUMMARY_FILENAME).read_text()),
                summary,
            )

    def test_failed_ledger_publication_cannot_create_summary_and_failure_exits_two(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            (results / runner.SUMMARY_FILENAME).write_text(
                '{"status":"STALE_PRIOR_RUN"}\n', encoding="utf-8"
            )
            actual_write = runner._atomic_write_bytes

            def fail_on_ledger(path: Path, payload: bytes) -> None:
                if path.name == runner.OUTPUT_HASH_FILENAME:
                    raise OSError("synthetic ledger failure")
                actual_write(path, payload)

            with patch.object(
                runner,
                "reproduce_reference_dynamics",
                return_value=_synthetic_reproduction(
                    "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_PASS"
                ),
            ), patch.object(
                runner, "_atomic_write_bytes", side_effect=fail_on_ledger
            ):
                with self.assertRaises(OSError):
                    runner.run_independent_native_dynamic(repository)
            self.assertFalse((results / runner.SUMMARY_FILENAME).exists())

        with tempfile.TemporaryDirectory() as directory:
            repository, _ = _fixture(directory)
            with patch.object(
                runner,
                "reproduce_reference_dynamics",
                return_value=_synthetic_reproduction(
                    "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_DISCREPANCY"
                ),
            ), patch("builtins.print"):
                self.assertEqual(
                    runner.main(["--repository", str(repository)]), 2
                )

    def test_end_of_run_input_drift_aborts_before_publication(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository, results = _fixture(directory)
            stale_summary = results / runner.SUMMARY_FILENAME
            stale_summary.write_text(
                '{"status":"VALID_PRIOR_RUN"}\n', encoding="utf-8"
            )

            def drift_after_reproduction(*_args: object) -> object:
                (results / "native_source_map.json").write_text(
                    "mid-run drift\n", encoding="utf-8"
                )
                return _synthetic_reproduction(
                    "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_PASS"
                )

            with patch.object(
                runner,
                "reproduce_reference_dynamics",
                side_effect=drift_after_reproduction,
            ):
                with self.assertRaises(ValueError):
                    runner.run_independent_native_dynamic(repository)
            self.assertEqual(
                json.loads(stale_summary.read_text(encoding="utf-8"))["status"],
                "VALID_PRIOR_RUN",
            )
            self.assertFalse((results / runner.OUTPUT_FILENAME).exists())


if __name__ == "__main__":
    unittest.main()
