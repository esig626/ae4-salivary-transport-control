"""Tests for the target-free native genotype production runner."""

from __future__ import annotations

import csv
import hashlib
import inspect
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model.genotype_evaluation import (
    DEFAULT_EXPRESSION_GRID,
    ContinuationAttempt,
    ContinuationStep,
    DualSolverGenotypeEvaluation,
    GenotypeContinuation,
    GenotypePairSummary,
    GenotypeRoot,
    GenotypeTrajectory,
    GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
)
from modern_full_model.run_native_genotype_evaluation import (
    CONTINUATION_FILENAME,
    DYNAMIC_GATE_INPUT,
    DYNAMIC_GROUP_INPUT,
    DYNAMIC_MANIFEST_INPUT,
    EXPECTED_DYNAMIC_GATE_ID,
    EXPECTED_DYNAMIC_GATE_STATUS,
    EXPECTED_DYNAMIC_MANIFEST_ID,
    EXPECTED_DYNAMIC_MANIFEST_STATUS,
    HASH_FILENAME,
    MANIFEST_FILENAME,
    NATIVE_ROOT_INPUT,
    NATIVE_SUMMARY_INPUT,
    PROFILE_FILENAME,
    SOLVER_FILENAME,
    SUMMARY_FILENAME,
    TRAJECTORY_FILENAME,
    NativeGenotypeRootResult,
    load_frozen_genotype_intake,
    match_passing_native_roots,
    passing_dynamic_group_rows,
    run_native_genotype_evaluation,
)
from modern_full_model.states import CORE_STATE_NAMES
from modern_full_model.validation import (
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    sha256_object,
)
import modern_full_model.run_native_genotype_evaluation as native_runner


REPOSITORY = Path(__file__).resolve().parents[1]


def _native_root(root_id: str) -> SimpleNamespace:
    return SimpleNamespace(
        root_id=root_id,
        panel_id=f"PANEL_{root_id}",
        topology_id="PNOM_KMID",
        routing_id="AE4NA05",
        source_class="N_ABS_NKCC",
        source_scale=4.0,
        pump_capacity_scale=1.0,
        apical_pump_fraction=0.075075,
        apical_k_fraction=0.30,
        ae4_cation_fraction=0.05,
        hydraulic_scale=1.0,
        hydraulic_mode="H1",
        eligibility="PRODUCTION_PRE_REVEAL",
        passes_numerical_gate=True,
        passes_wt_gate=True,
        core_state=tuple(float(index + 1) for index in range(12)),
        other_impermeant_osmoles_fmol=125.0,
    )


def _group_row(root_id: str, passed: bool) -> dict[str, str]:
    return {
        "root_id": root_id,
        "hydraulic_scale": "1.0",
        "hydraulic_mode": "H1",
        "all_native_rest_wt_gate_pass": "True",
        "native_dynamic_scientific_gate_pass": str(passed),
    }


def _refresh_dynamic_hash_ledger(directory: Path) -> None:
    with (directory / native_runner.DYNAMIC_HASH_INPUT).open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=("artifact", "sha256"))
        writer.writeheader()
        writer.writerows(
            {
                "artifact": filename,
                "sha256": hashlib.sha256((directory / filename).read_bytes()).hexdigest(),
            }
            for filename in native_runner.DYNAMIC_HASH_LEDGER_ARTIFACTS
        )


def _write_gate_fixture(directory: Path, roots: tuple[SimpleNamespace, ...]) -> None:
    rows = [_group_row(root.root_id, True) for root in reversed(roots)]
    rows.append(_group_row("rejected_row", False))
    with (directory / DYNAMIC_GROUP_INPUT).open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for filename in native_runner.DYNAMIC_REPORT_INPUTS:
        if filename != DYNAMIC_GROUP_INPUT:
            (directory / filename).write_text(
                f"synthetic frozen report {filename}\n", encoding="utf-8"
            )
    (directory / NATIVE_ROOT_INPUT).write_text(
        "synthetic root table; loader is mocked\n", encoding="utf-8"
    )
    (directory / NATIVE_SUMMARY_INPUT).write_text(
        json.dumps(
            {
                "production_confirmation_complete": True,
                "geometry_confirmation_complete": True,
                "independent_root_reproduction_complete": True,
                "native_pre_reveal_freeze_ready": True,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    for filename in native_runner.NATIVE_INPUT_ARTIFACTS:
        path = directory / filename
        if not path.exists():
            path.write_text(f"synthetic native input {filename}\n", encoding="utf-8")
    report_hashes = {
        filename: hashlib.sha256((directory / filename).read_bytes()).hexdigest()
        for filename in native_runner.DYNAMIC_REPORT_INPUTS
    }
    native_hashes = {
        filename: hashlib.sha256((directory / filename).read_bytes()).hexdigest()
        for filename in native_runner.NATIVE_INPUT_ARTIFACTS
    }
    equation_hashes = {
        filename: hashlib.sha256((REPOSITORY / filename).read_bytes()).hexdigest()
        for filename in native_runner.DYNAMIC_EQUATION_SOURCES
    }
    passing_ids = [root.root_id for root in roots]
    (directory / DYNAMIC_GATE_INPUT).write_text(
        json.dumps(
            {
                "gate_id": EXPECTED_DYNAMIC_GATE_ID,
                "status": EXPECTED_DYNAMIC_GATE_STATUS,
                "manifest": DYNAMIC_MANIFEST_INPUT,
                "native_root_count": len(rows),
                "at_least_one_scientific_group_pass": True,
                "all_rows_retained_no_first_passer": True,
                "scientific_passing_group_count": len(roots),
                "scientific_passing_root_ids": passing_ids,
                "clean_auditor_authorization_required": True,
                "clean_auditor_authorization_recorded": False,
                "independent_numerical_reproduction_required": True,
                "independent_numerical_reproduction_recorded_by_this_runner": False,
                "final_composite_authorization_complete": False,
                "phenotype_reveal_permitted": False,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    (directory / DYNAMIC_MANIFEST_INPUT).write_text(
        json.dumps(
            {
                "manifest_id": EXPECTED_DYNAMIC_MANIFEST_ID,
                "status": EXPECTED_DYNAMIC_MANIFEST_STATUS,
                "root_count": len(rows),
                "roots": {
                    **{
                        root.root_id: {
                            "core_state_sha256": sha256_object(root.core_state)
                        }
                        for root in roots
                    },
                    "rejected_row": {"core_state_sha256": "0" * 64},
                },
                "scientific_passing_group_count": len(roots),
                "scientific_passing_root_ids": passing_ids,
                "clean_auditor_authorization_required": True,
                "independent_numerical_reproduction_required_for_composite_authorization": True,
                "phenotype_reveal_permitted_by_this_manifest": False,
                "wt_dynamic_report_artifact_sha256": report_hashes,
                "native_input_artifact_sha256": native_hashes,
                "equation_source_file_sha256": equation_hashes,
                "equation_source_tree_sha256": sha256_object(equation_hashes),
                "opaque_firewall_ledger_sha256": {
                    key: f"{index + 1:x}" * 64
                    for index, key in enumerate(sorted(native_runner.FIREWALL_HASH_KEYS))
                },
                "opaque_firewall_access_mode": "SHA256_BYTES_ONLY_NO_TARGET_PARSE",
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    _refresh_dynamic_hash_ledger(directory)


def _attempt(expression: float) -> ContinuationAttempt:
    return ContinuationAttempt(
        start_id="previous_connected_root",
        optimizer_success=True,
        admissible_state=True,
        converged_root=True,
        max_abs_scaled_independent_rhs=1.0e-12,
        cost=0.0,
        optimality=1.0e-12,
        nfev=4,
        message="synthetic",
        coordinates=tuple(float(index + 1) for index in range(10)),
    )


def _continuation(transporter: str, *, complete: bool) -> GenotypeContinuation:
    steps: list[ContinuationStep] = []
    grid = DEFAULT_EXPRESSION_GRID
    stop = len(grid) if complete else 11
    for index, expression in enumerate(grid[:stop]):
        failing = not complete and index == stop - 1
        if failing:
            steps.append(
                ContinuationStep(
                    expression=expression,
                    attempts=(_attempt(expression),),
                    roots=(),
                    selected_root_id=None,
                    selected_distance_from_previous_normalized=None,
                    connected_to_previous=False,
                    alternate_root_count=0,
                    status="NO_ADMISSIBLE_ROOT",
                )
            )
            break
        root = GenotypeRoot(
            root_id=f"{transporter}_E{expression:.8f}:B00",
            transporter=transporter,
            expression=expression,
            member_start_ids=("previous_connected_root",),
            coordinates=tuple(float(index + 1) for index in range(10)),
            core_state=tuple(float(index + 1) for index in range(12)),
            max_abs_scaled_independent_rhs=1.0e-12,
            max_abs_amount_rhs_fmol_s=1.0e-13,
            max_abs_volume_rhs_pL_s=1.0e-14,
            max_abs_omitted_charge_rhs_fmol_equivalent_s=1.0e-13,
            max_abs_regulatory_rhs_s_inv=0.0,
            max_abs_charge_residual_fmol=1.0e-14,
            max_abs_current_residual_A=1.0e-22,
            normalized_jacobian_singular_values=tuple(
                float(10 - index) for index in range(10)
            ),
            normalized_jacobian_rank=10,
            normalized_jacobian_nullity=0,
            boundary_hits=(),
            exact_deleted_transport_zero=True,
            passes_numerical_gate=True,
            gate_failures=(),
        )
        steps.append(
            ContinuationStep(
                expression=expression,
                attempts=(_attempt(expression),),
                roots=(root,),
                selected_root_id=root.root_id,
                selected_distance_from_previous_normalized=0.01,
                connected_to_previous=True,
                alternate_root_count=0,
                status="CONNECTED",
            )
        )
    return GenotypeContinuation(
        transporter=transporter,
        expression_grid=grid,
        fixed_other_impermeant_osmoles_fmol=125.0,
        independent_equation_count=10,
        fitted_coordinate_count=10,
        fitted_parameter_count=0,
        steps=tuple(steps),
        completed_to_exact_zero=complete,
        branch_failure=None if complete else "synthetic branch loss",
    )


def _trajectory(transporter: str, solver: object, role: str) -> GenotypeTrajectory:
    time = np.asarray((0.0, 1.0e-6, 600.0))
    names = (*CORE_STATE_NAMES, "pka_fraction", "ae4_regulated_fraction")
    states = np.ones((len(names), len(time)))
    flow = np.asarray((1.0, 1.1, 1.2))
    cumulative = np.asarray((0.0, 1.0e-6, 660.0))
    base = np.asarray((10.0, 10.1, 10.2))
    return GenotypeTrajectory(
        genotype_name=role,
        transporter=transporter,
        solver=solver,
        success=True,
        message="synthetic",
        time_s=time,
        state_names=names,
        states=states,
        flow_pL_s=flow,
        cumulative_flow_pL=cumulative,
        cell_na_mM=base,
        cell_k_mM=base + 100.0,
        cell_cl_mM=base + 40.0,
        cell_ph=np.asarray((6.9, 6.91, 6.92)),
        cell_volume_pL=np.asarray((1.3, 1.3, 1.31)),
        lumen_na_mM=base + 130.0,
        lumen_k_mM=base - 5.0,
        lumen_cl_mM=base + 120.0,
        lumen_ph=np.asarray((7.2, 7.21, 7.22)),
        lumen_volume_pL=np.asarray((0.1, 0.1, 0.11)),
        ae4_capacity_multiplier=np.asarray((1.0, 1.0, 1.25)),
        max_abs_deleted_transport_flux_fmol_s=(
            0.0 if role == "EXACT_ZERO_EXPRESSION" else 1.0e-3
        ),
        max_abs_conservation_residuals={"synthetic": 0.0},
        max_dimensionless_conservation_ratio=0.01,
        positive_core=True,
        all_flow_nonnegative=True,
        numerical_gate_pass=True,
    )


def _dual(transporter: str) -> DualSolverGenotypeEvaluation:
    wt: dict[str, GenotypeTrajectory] = {}
    zero: dict[str, GenotypeTrajectory] = {}
    summaries: dict[str, GenotypePairSummary] = {}
    for solver in (PRODUCTION_RADAU, PRODUCTION_BDF):
        wt[solver.label] = _trajectory(transporter, solver, "FROZEN_WT")
        zero[solver.label] = _trajectory(
            transporter, solver, "EXACT_ZERO_EXPRESSION"
        )
        summaries[solver.label] = GenotypePairSummary(
            transporter=transporter,
            solver_label=solver.label,
            wt_success=True,
            genotype_success=True,
            paired_numerical_gate_pass=True,
            cumulative_ratio_at_landmarks={"600_s": 0.9},
            flow_ratio_at_landmarks={"600_s": 0.8},
            wt_endpoint={"cell_cl_mM": 50.0, "cell_ph": 6.9},
            genotype_endpoint={"cell_cl_mM": 48.0, "cell_ph": 6.8},
            genotype_to_wt_endpoint_ratios={
                "cell_cl_mM": 0.96,
                "cell_ph_difference": -0.1,
            },
            wt_max_dimensionless_conservation_ratio=0.01,
            genotype_max_dimensionless_conservation_ratio=0.01,
        )
    return DualSolverGenotypeEvaluation(
        transporter=transporter,
        evaluations=summaries,
        wt_trajectories=wt,
        genotype_trajectories=zero,
        cross_solver_relative_differences={
            "wt_state": 1.0e-7,
            "wt_flow": 1.0e-7,
            "wt_cumulative": 1.0e-7,
            "genotype_state": 1.0e-7,
            "genotype_flow": 1.0e-7,
            "genotype_cumulative": 1.0e-7,
        },
        cross_solver_relative_tolerance=GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        cross_solver_gate_pass=True,
        both_solver_pairs_pass=True,
    )


def _root_result(root: SimpleNamespace, success_transporter: str) -> NativeGenotypeRootResult:
    continuations = {
        transporter: _continuation(
            transporter, complete=transporter == success_transporter
        )
        for transporter in ("AE4", "AE2")
    }
    return NativeGenotypeRootResult(
        root=root,
        continuations=continuations,
        dynamics={success_transporter: _dual(success_transporter)},
    )


class TestNativeGenotypeIntake(unittest.TestCase):
    def test_runner_remains_blocked_until_all_frozen_gate_inputs_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with self.assertRaisesRegex(RuntimeError, "blocked"):
                load_frozen_genotype_intake(output, repository_root=REPOSITORY)
            for filename in (
                CONTINUATION_FILENAME,
                TRAJECTORY_FILENAME,
                PROFILE_FILENAME,
                SOLVER_FILENAME,
                SUMMARY_FILENAME,
                MANIFEST_FILENAME,
                HASH_FILENAME,
            ):
                self.assertFalse((output / filename).exists())

    def test_all_passing_rows_are_sorted_and_matched_without_first_selection(self) -> None:
        rows = (
            _group_row("root_b", True),
            _group_row("rejected", False),
            _group_row("root_a", True),
        )
        passing = passing_dynamic_group_rows(rows)
        self.assertEqual([row["root_id"] for row in passing], ["root_a", "root_b"])
        roots = (_native_root("root_b"), _native_root("unused"), _native_root("root_a"))
        matched = match_passing_native_roots(passing, roots)
        self.assertEqual([root.root_id for root in matched], ["root_a", "root_b"])

    def test_frozen_intake_crosschecks_gate_manifest_and_core_hash(self) -> None:
        roots = (_native_root("root_a"), _native_root("root_b"))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_gate_fixture(output, roots)
            with patch.object(native_runner, "load_native_source_roots", return_value=roots):
                intake = load_frozen_genotype_intake(
                    output, repository_root=REPOSITORY
                )
            self.assertEqual(
                [root.root_id for root in intake.roots], ["root_a", "root_b"]
            )
            self.assertEqual(len(intake.group_rows), 2)
            self.assertEqual(
                len(intake.input_artifact_hashes),
                len(native_runner.FROZEN_INPUT_ARTIFACTS),
            )

    def test_noncanonical_csv_boole_missing_h1_and_duplicate_rows_fail_closed(self) -> None:
        with self.subTest("lowercase boolean"):
            row = _group_row("root_a", True)
            row["native_dynamic_scientific_gate_pass"] = "true"
            with self.assertRaisesRegex(ValueError, "canonical"):
                passing_dynamic_group_rows((row,))
        with self.subTest("missing H1 mode"):
            row = _group_row("root_a", True)
            del row["hydraulic_mode"]
            with self.assertRaisesRegex(ValueError, "mode H1"):
                passing_dynamic_group_rows((row,))
        with self.subTest("duplicate group root"):
            row = _group_row("root_a", True)
            with self.assertRaisesRegex(ValueError, "duplicate root"):
                passing_dynamic_group_rows((row, dict(row)))
        with self.subTest("duplicate frozen root"):
            rows = passing_dynamic_group_rows((_group_row("root_a", True),))
            with self.assertRaisesRegex(ValueError, "duplicate root"):
                match_passing_native_roots(
                    rows, (_native_root("root_a"), _native_root("root_a"))
                )
        with self.subTest("production-prefix impostor"):
            rows = passing_dynamic_group_rows((_group_row("root_a", True),))
            rejected = _native_root("root_a")
            rejected.eligibility = "PRODUCTION_REJECTED"
            with self.assertRaisesRegex(ValueError, "exact token"):
                match_passing_native_roots(rows, (rejected,))

    def test_json_string_boole_and_recorded_authority_fail_closed(self) -> None:
        root = _native_root("root_a")
        mutations = (
            ("at_least_one_scientific_group_pass", "True"),
            ("clean_auditor_authorization_recorded", "False"),
            ("clean_auditor_authorization_recorded", True),
            ("final_composite_authorization_complete", True),
            ("phenotype_reveal_permitted", True),
        )
        for field, value in mutations:
            with self.subTest(field=field, value=value):
                with tempfile.TemporaryDirectory() as temporary:
                    output = Path(temporary)
                    _write_gate_fixture(output, (root,))
                    gate = json.loads(
                        (output / DYNAMIC_GATE_INPUT).read_text(encoding="utf-8")
                    )
                    gate[field] = value
                    (output / DYNAMIC_GATE_INPUT).write_text(
                        json.dumps(gate, sort_keys=True), encoding="utf-8"
                    )
                    with self.assertRaises((ValueError, RuntimeError)):
                        load_frozen_genotype_intake(
                            output, repository_root=REPOSITORY
                        )

    def test_composite_status_root_id_mismatch_and_stale_hash_fail_closed(self) -> None:
        root = _native_root("root_a")
        with self.subTest("composite status"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                gate = json.loads(
                    (output / DYNAMIC_GATE_INPUT).read_text(encoding="utf-8")
                )
                gate["status"] = "SCIENTIFIC_GATE_FAIL"
                (output / DYNAMIC_GATE_INPUT).write_text(
                    json.dumps(gate, sort_keys=True), encoding="utf-8"
                )
                with self.assertRaisesRegex(RuntimeError, "composite"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )
        with self.subTest("manifest passing IDs"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                manifest = json.loads(
                    (output / DYNAMIC_MANIFEST_INPUT).read_text(encoding="utf-8")
                )
                manifest["scientific_passing_root_ids"] = ["different_root"]
                (output / DYNAMIC_MANIFEST_INPUT).write_text(
                    json.dumps(manifest, sort_keys=True), encoding="utf-8"
                )
                with self.assertRaisesRegex(ValueError, "disagree"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )
        with self.subTest("stale native root hash"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                with (output / NATIVE_ROOT_INPUT).open("a", encoding="utf-8") as handle:
                    handle.write("stale mutation\n")
                with self.assertRaisesRegex(ValueError, "stale manifest hash"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )

    def test_duplicate_json_root_keys_are_rejected_before_loading(self) -> None:
        root = _native_root("root_a")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_gate_fixture(output, (root,))
            manifest = (output / DYNAMIC_MANIFEST_INPUT).read_text(encoding="utf-8")
            marker = '"roots": {'
            duplicate = (
                marker
                + '"root_a":{"core_state_sha256":"duplicate"},'
            )
            (output / DYNAMIC_MANIFEST_INPUT).write_text(
                manifest.replace(marker, duplicate, 1), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_frozen_genotype_intake(
                    output, repository_root=REPOSITORY
                )

    def test_exact_dynamic_ledger_and_complete_report_set_fail_closed(self) -> None:
        root = _native_root("root_a")
        with self.subTest("missing report"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                (output / native_runner.DYNAMIC_PROFILE_INPUT).unlink()
                with self.assertRaisesRegex(RuntimeError, "missing artifacts"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )
        with self.subTest("incomplete ledger"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                ledger = native_runner._read_csv(
                    output / native_runner.DYNAMIC_HASH_INPUT
                )[:-1]
                with (output / native_runner.DYNAMIC_HASH_INPUT).open(
                    "w", newline="", encoding="utf-8"
                ) as handle:
                    writer = csv.DictWriter(
                        handle, fieldnames=("artifact", "sha256")
                    )
                    writer.writeheader()
                    writer.writerows(ledger)
                with self.assertRaisesRegex(ValueError, "exactly seven"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )
        with self.subTest("stale ledger target"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                with (output / native_runner.DYNAMIC_PROFILE_INPUT).open(
                    "a", encoding="utf-8"
                ) as handle:
                    handle.write("stale\n")
                with self.assertRaisesRegex(ValueError, "stale manifest hash"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )

    def test_every_manifest_hash_registry_must_be_exact_and_current(self) -> None:
        root = _native_root("root_a")
        cases = (
            ("wt_dynamic_report_artifact_sha256", DYNAMIC_GROUP_INPUT),
            ("native_input_artifact_sha256", NATIVE_ROOT_INPUT),
            (
                "equation_source_file_sha256",
                native_runner.DYNAMIC_EQUATION_SOURCES[0],
            ),
            (
                "opaque_firewall_ledger_sha256",
                sorted(native_runner.FIREWALL_HASH_KEYS)[0],
            ),
        )
        for registry_name, removed_key in cases:
            with self.subTest(registry=registry_name):
                with tempfile.TemporaryDirectory() as temporary:
                    output = Path(temporary)
                    _write_gate_fixture(output, (root,))
                    manifest = json.loads(
                        (output / DYNAMIC_MANIFEST_INPUT).read_text(encoding="utf-8")
                    )
                    del manifest[registry_name][removed_key]
                    (output / DYNAMIC_MANIFEST_INPUT).write_text(
                        json.dumps(manifest, sort_keys=True), encoding="utf-8"
                    )
                    _refresh_dynamic_hash_ledger(output)
                    with self.assertRaisesRegex(ValueError, "exact"):
                        load_frozen_genotype_intake(
                            output, repository_root=REPOSITORY
                        )

    def test_all_group_root_ids_and_root_counts_close_exactly(self) -> None:
        root = _native_root("root_a")
        with self.subTest("manifest omits a nonpassing group row"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                manifest = json.loads(
                    (output / DYNAMIC_MANIFEST_INPUT).read_text(encoding="utf-8")
                )
                del manifest["roots"]["rejected_row"]
                manifest["root_count"] -= 1
                (output / DYNAMIC_MANIFEST_INPUT).write_text(
                    json.dumps(manifest, sort_keys=True), encoding="utf-8"
                )
                _refresh_dynamic_hash_ledger(output)
                with self.assertRaisesRegex(ValueError, "every group CSV root ID"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )
        with self.subTest("gate count is stale"):
            with tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary)
                _write_gate_fixture(output, (root,))
                gate = json.loads(
                    (output / DYNAMIC_GATE_INPUT).read_text(encoding="utf-8")
                )
                gate["native_root_count"] -= 1
                (output / DYNAMIC_GATE_INPUT).write_text(
                    json.dumps(gate, sort_keys=True), encoding="utf-8"
                )
                _refresh_dynamic_hash_ledger(output)
                with self.assertRaisesRegex(ValueError, "root counts"):
                    load_frozen_genotype_intake(
                        output, repository_root=REPOSITORY
                    )

    def test_genotype_source_hashes_are_complete_and_fail_closed(self) -> None:
        hashes = native_runner._source_hashes(REPOSITORY)
        self.assertEqual(
            tuple(hashes), native_runner.GENOTYPE_SOURCE_DEPENDENCIES
        )
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(FileNotFoundError, "dependency set"):
                native_runner._source_hashes(Path(temporary))


class TestNativeGenotypeArtifactRunner(unittest.TestCase):
    def test_model_build_exception_is_retained_for_both_transporter_branches(self) -> None:
        root = _native_root("root_failure")
        with patch.object(
            native_runner,
            "build_native_source_model",
            side_effect=RuntimeError("synthetic build failure"),
        ):
            result = native_runner._evaluate_root(root)
        self.assertEqual(set(result.branch_errors), {"AE4", "AE2"})
        continuation = native_runner._continuation_rows((result,))
        profile = native_runner._profile_rows((result,))
        solvers = native_runner._solver_rows((result,))
        self.assertEqual(len(continuation), 2)
        self.assertEqual(len(profile), 2)
        self.assertEqual(len(solvers), 2)
        self.assertTrue(all(row["runner_exception"] for row in continuation))
        self.assertTrue(all(row["runner_exception"] for row in profile))
        self.assertTrue(all(row["runner_exception"] for row in solvers))

    def test_dependency_snapshot_detects_input_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            repository = base / "repository"
            output = base / "output"
            repository.mkdir()
            output.mkdir()
            (repository / "source.py").write_text("source-v1\n", encoding="utf-8")
            (output / "input.csv").write_text("input-v1\n", encoding="utf-8")
            with (
                patch.object(
                    native_runner, "GENOTYPE_SOURCE_DEPENDENCIES", ("source.py",)
                ),
                patch.object(
                    native_runner, "FROZEN_INPUT_ARTIFACTS", ("input.csv",)
                ),
            ):
                snapshot = native_runner._capture_genotype_dependency_snapshot(
                    repository, output
                )
                native_runner._assert_genotype_snapshot_current(
                    snapshot, repository, output
                )
                (output / "input.csv").write_text("input-v2\n", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "changed during execution"):
                    native_runner._assert_genotype_snapshot_current(
                        snapshot, repository, output
                    )

    def test_failed_staged_commit_leaves_no_authority_ledger(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            repository = base / "repository"
            output = base / "output"
            staging = base / "staging"
            repository.mkdir()
            output.mkdir()
            staging.mkdir()
            for filename in native_runner.GENOTYPE_FREEZE_ARTIFACT_FILENAMES:
                (staging / filename).write_text(filename, encoding="utf-8")
            (staging / HASH_FILENAME).write_text("new ledger", encoding="utf-8")
            (output / HASH_FILENAME).write_text("stale ledger", encoding="utf-8")
            snapshot = native_runner.GenotypeDependencySnapshot({}, {})
            with patch.object(
                native_runner,
                "_assert_genotype_snapshot_current",
                side_effect=(None, RuntimeError("dependency changed")),
            ):
                with self.assertRaisesRegex(RuntimeError, "dependency changed"):
                    native_runner._publish_staged_genotype_freeze(
                        staging, output, repository, snapshot
                    )
            self.assertFalse((output / HASH_FILENAME).exists())

    def test_prior_ledger_is_invalidated_before_numerical_root_work(self) -> None:
        source = inspect.getsource(native_runner.run_native_genotype_evaluation)
        self.assertLess(
            source.index("_invalidate_genotype_hash_ledger"),
            source.index("_run_roots"),
        )

    def test_mocked_full_runner_writes_deterministic_complete_artifact_set(self) -> None:
        roots = (_native_root("root_a"), _native_root("root_b"))
        root_results = (
            _root_result(roots[0], "AE4"),
            _root_result(roots[1], "AE2"),
        )
        repository = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_gate_fixture(output, roots)
            with (
                patch.object(native_runner, "load_native_source_roots", return_value=roots),
                patch.object(native_runner, "_run_roots", return_value=root_results),
            ):
                summary = run_native_genotype_evaluation(repository, output)
                first_hashes = {
                    filename: hashlib.sha256((output / filename).read_bytes()).hexdigest()
                    for filename in (
                        CONTINUATION_FILENAME,
                        TRAJECTORY_FILENAME,
                        PROFILE_FILENAME,
                        SOLVER_FILENAME,
                        SUMMARY_FILENAME,
                        MANIFEST_FILENAME,
                        HASH_FILENAME,
                    )
                }
                run_native_genotype_evaluation(repository, output)
                second_hashes = {
                    filename: hashlib.sha256((output / filename).read_bytes()).hexdigest()
                    for filename in first_hashes
                }
            self.assertEqual(first_hashes, second_hashes)
            self.assertEqual(summary["evaluated_root_count"], 2)
            self.assertTrue(summary["all_passing_rows_retained"])
            self.assertEqual(summary["dynamic_pair_count"], 2)
            self.assertEqual(len(summary["branch_losses"]), 2)

            with (output / TRAJECTORY_FILENAME).open(
                newline="", encoding="utf-8"
            ) as handle:
                trajectories = list(csv.DictReader(handle))
            self.assertEqual(len(trajectories), 24)
            self.assertEqual(
                {row["genotype_role"] for row in trajectories},
                {"FROZEN_WT", "EXACT_ZERO_EXPRESSION"},
            )
            with (output / PROFILE_FILENAME).open(
                newline="", encoding="utf-8"
            ) as handle:
                profiles = list(csv.DictReader(handle))
            self.assertEqual(len(profiles), 6)
            self.assertEqual(
                sum(row["dynamic_executed"] == "True" for row in profiles), 4
            )
            with (output / HASH_FILENAME).open(
                newline="", encoding="utf-8"
            ) as handle:
                hashes = list(csv.DictReader(handle))
            self.assertEqual(len(hashes), 6)
            self.assertEqual(
                {row["artifact"] for row in hashes},
                {
                    CONTINUATION_FILENAME,
                    TRAJECTORY_FILENAME,
                    PROFILE_FILENAME,
                    SOLVER_FILENAME,
                    SUMMARY_FILENAME,
                    MANIFEST_FILENAME,
                },
            )

    def test_runner_source_has_no_experimental_target_reader_or_score(self) -> None:
        source = inspect.getsource(native_runner)
        self.assertNotIn("heldout_targets.csv", source)
        self.assertNotIn("validation_targets.csv", source)
        self.assertNotIn("phenotype_gate_pass", source)
        self.assertNotIn("target_lower", source)
        self.assertNotIn("target_upper", source)


if __name__ == "__main__":
    unittest.main()
