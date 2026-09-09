"""Target-free tests for the independent native genotype continuation."""

from __future__ import annotations

import csv
import inspect
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model.independent import (
    IndependentConstantStimulus,
    IndependentN1NkccRegulation,
    IndependentSecretagogueProtocol,
    independent_r3_regulation,
)
from modern_full_model.independent_genotype import (
    COORDINATE_LOWER,
    COORDINATE_UPPER,
    IndependentGenotypeContinuation,
    IndependentGenotypeRoot,
    _cluster_roots,
    _one_attempt,
    independent_coordinates_from_state,
)
from modern_full_model.independent_native import (
    N_ABS_NKCC,
    IndependentNativeSourceDefinition,
    state_from_independent_rest_coordinates,
)
from modern_full_model.run_independent_native_genotypes import (
    DYNAMIC_FREEZE_ARTIFACTS,
    GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
    INDEPENDENT_GENOTYPE_FREEZE_ARTIFACTS,
    NATIVE_GENOTYPE_INPUTS,
    NATIVE_GENOTYPE_SOURCE_PATHS,
    OPAQUE_FIREWALL_KEYS,
    PRODUCTION_CONTINUATION,
    PRODUCTION_GENOTYPE_FREEZE_ARTIFACTS,
    PRODUCTION_GENOTYPE_REPORTS,
    PRODUCTION_HASH,
    PRODUCTION_MANIFEST,
    PRODUCTION_PROFILE,
    PRODUCTION_SOLVER,
    PRODUCTION_SUMMARY,
    HASH_OUTPUT,
    IndependentGenotypeDependencySnapshot,
    _eligible_group_ids,
    _failed_root_rows,
    _matching_root_rows,
    _production_discrepancies,
    _run_one_root,
    _sha256_file,
    _sha256_object,
    _validate_native_genotype_freeze,
    _validate_dynamic_hash_ledger,
    _validate_exact_current_registry,
    _validate_opaque_registry,
)
import modern_full_model.run_independent_native_genotypes as independent_runner


KNOWN_N4_STATE = np.asarray(
    (
        24.375004480505822,
        138.68089788712203,
        66.94776908416362,
        6.830518062268592,
        18.04637107751941,
        1.3,
        13.534346011579268,
        2.243735106264099,
        15.265233459193988,
        0.589897675936328,
        0.5128476586493772,
        0.10495597693809071,
    ),
    dtype=float,
)


def definition() -> IndependentNativeSourceDefinition:
    return IndependentNativeSourceDefinition(
        source_class=N_ABS_NKCC,
        source_scale=4.0,
        pump_capacity_scale=1.0,
        apical_pump_fraction=0.075075,
        apical_k_fraction=0.30,
        ae4_cation_fraction=0.05,
        hydraulic_scale=1.0,
        hydraulic_mode="H1",
        cell_other_impermeant_osmoles_fmol=124.69936937127669,
    )


def dynamic_model():
    regulation = independent_r3_regulation(
        tau_pka_s=10.0,
        forward_regulation_rate_s=1.0 / 30.0,
        reverse_regulation_rate_s=1.0 / 30.0,
        basal_multiplier=1.0,
        fully_activated_increment=0.25,
    )
    return definition().build_independent_model(
        regulation=regulation,
        nkcc_regulation=IndependentN1NkccRegulation(1.75, 0.058, 0.10),
        stimulus=IndependentSecretagogueProtocol(
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.10,
            beta_occupancy_on=1.0,
        ),
    )


class TestIndependentNativeGenotype(unittest.TestCase):
    @staticmethod
    def _synthetic_root(
        label: str, coordinate: float, residual: float
    ) -> IndependentGenotypeRoot:
        coordinates = (coordinate, *([1.0] * 9))
        return IndependentGenotypeRoot(
            root_id="",
            transporter="AE4",
            expression=0.5,
            member_start_ids=(label,),
            coordinates=coordinates,
            core_state=tuple(float(index + 1) for index in range(12)),
            max_abs_scaled_residual=residual,
            max_abs_omitted_charge_rhs_fmol_s=0.0,
            max_abs_bulk_charge_residual_fmol=0.0,
            max_abs_current_residual_A=0.0,
            max_abs_regulatory_rhs_s_inv=0.0,
            jacobian_singular_values=tuple([1.0] * 10),
            jacobian_rank=10,
            jacobian_nullity=0,
            boundary_hits=(),
            exact_deleted_transport_zero=True,
            passes_numerical_gate=True,
            gate_failures=(),
        )

    def test_chain_connected_components_and_ids_are_permutation_invariant(self) -> None:
        first_span = float((COORDINATE_UPPER - COORDINATE_LOWER)[0])
        base = float(COORDINATE_LOWER[0] + 0.2 * first_span)
        roots = (
            self._synthetic_root("a", base, 0.03),
            self._synthetic_root("b", base + 1.5e-5 * first_span, 0.01),
            self._synthetic_root("c", base + 3.0e-5 * first_span, 0.02),
            self._synthetic_root("d", base + 0.1 * first_span, 0.04),
        )
        signatures = set()
        for permutation in itertools.permutations(roots):
            clustered = _cluster_roots(
                permutation,
                transporter="AE4",
                expression=0.5,
            )
            signatures.add(
                tuple(
                    (root.root_id, root.coordinates, root.member_start_ids)
                    for root in clustered
                )
            )
        self.assertEqual(len(signatures), 1)
        signature = signatures.pop()
        self.assertEqual([row[0] for row in signature], ["AE4_E0.50000000:IB00", "AE4_E0.50000000:IB01"])
        self.assertEqual(signature[0][2], ("a", "b", "c"))
        self.assertAlmostEqual(signature[0][1][0], roots[1].coordinates[0])

    def test_modules_do_not_import_production_genotype_or_target_reader(self) -> None:
        import modern_full_model.independent_genotype as equations
        import modern_full_model.run_independent_native_genotypes as runner

        source = inspect.getsource(equations) + inspect.getsource(runner)
        self.assertNotIn("from .genotype_evaluation", source)
        self.assertNotIn("from .native_source_panel", source)
        self.assertNotIn("heldout_targets", source)

    def test_numerical_solver_tolerance_matches_existing_native_contract(self) -> None:
        from modern_full_model.genotype_evaluation import (
            GENOTYPE_SOLVER_RELATIVE_TOLERANCE as production_tolerance,
        )

        self.assertEqual(GENOTYPE_SOLVER_RELATIVE_TOLERANCE, production_tolerance)

    def test_root_exception_retains_both_branch_failure_records(self) -> None:
        continuation, profile, solver = _failed_root_rows(
            {"root_id": "R1"}, RuntimeError("synthetic")
        )
        for rows in (continuation, profile, solver):
            self.assertEqual({row["transporter"] for row in rows}, {"AE4", "AE2"})
            self.assertTrue(all(row["runner_exception"] for row in rows))
        self.assertTrue(all(row["cross_solver_gate_pass"] is None for row in solver))

    def test_continuation_exception_does_not_discard_other_branch(self) -> None:
        branch_loss = IndependentGenotypeContinuation(
            transporter="AE2",
            expression_grid=(1.0,),
            steps=(),
            completed_to_exact_zero=False,
            branch_failure="synthetic branch loss",
        )

        def continue_branch(model, core, *, transporter):
            if transporter == "AE4":
                raise RuntimeError("synthetic AE4 exception")
            return branch_loss

        row = {
            "root_id": "R1",
            "core_state_json": json.dumps([1.0] * 12),
        }
        with (
            patch.object(
                independent_runner.IndependentNativeSourceDefinition,
                "from_root_row",
                return_value=object(),
            ),
            patch.object(independent_runner, "_dynamic_model", return_value=object()),
            patch.object(
                independent_runner,
                "_initial_state",
                return_value=np.ones(12, dtype=float),
            ),
            patch.object(
                independent_runner,
                "independent_physical_time_grid",
                return_value=np.asarray((0.0, 1.0)),
            ),
            patch.object(
                independent_runner,
                "continue_independent_genotype_expression",
                side_effect=continue_branch,
            ),
        ):
            continuation, profile, solver = _run_one_root(row)

        self.assertEqual(len(continuation), 1)
        self.assertEqual(continuation[0]["transporter"], "AE4")
        self.assertTrue(continuation[0]["runner_exception"])
        by_transporter = {item["transporter"]: item for item in solver}
        self.assertTrue(by_transporter["AE4"]["runner_exception"])
        self.assertFalse(by_transporter["AE2"]["runner_exception"])
        self.assertEqual(
            by_transporter["AE2"]["continuation_failure"],
            "synthetic branch loss",
        )
        self.assertEqual({item["transporter"] for item in profile}, {"AE4", "AE2"})

    def test_exact_production_comparison_passes_then_fails_on_ratio_discrepancy(self) -> None:
        continuation = [
            {
                "root_id": "R1",
                "transporter": "AE4",
                "expression_step_index": 0,
                "continuation_completed_to_exact_zero": True,
                "connected_to_previous": True,
                "runner_exception": False,
                "selected_core_state_json": json.dumps([1.0, 2.0]),
            }
        ]
        profile = [
            {
                "root_id": "R1",
                "transporter": "AE4",
                "solver_label": "production_radau",
                "dynamic_executed": True,
                "continuation_completed_to_exact_zero": True,
                "paired_numerical_gate_pass": True,
                "runner_exception": False,
                "cumulative_ratio_at_landmarks_json": json.dumps({"60_s": 0.9}),
                "flow_ratio_at_landmarks_json": json.dumps({"60_s": 0.8}),
            }
        ]
        solver = [
            {
                "root_id": "R1",
                "transporter": "AE4",
                "dynamic_executed": True,
                "continuation_completed_to_exact_zero": True,
                "cross_solver_gate_pass": True,
                "runner_exception": False,
            }
        ]

        def write_csv(path: Path, rows) -> None:
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
                writer.writeheader()
                writer.writerows(rows)

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            write_csv(output / PRODUCTION_CONTINUATION, continuation)
            write_csv(output / PRODUCTION_PROFILE, profile)
            write_csv(output / PRODUCTION_SOLVER, solver)
            passing = _production_discrepancies(
                output, continuation, profile, solver
            )
            self.assertTrue(passing["comparison_gate_pass"])
            changed_profile = [dict(profile[0])]
            changed_profile[0]["flow_ratio_at_landmarks_json"] = json.dumps(
                {"60_s": 0.7}
            )
            failing = _production_discrepancies(
                output, continuation, changed_profile, solver
            )
            self.assertFalse(failing["comparison_gate_pass"])
            self.assertTrue(failing["profile_key_discrepancies"])

    def test_physiological_coordinate_roundtrip_is_independent(self) -> None:
        model = definition().build_independent_model(
            stimulus=IndependentConstantStimulus(0.058, 0.0)
        )
        coordinates = independent_coordinates_from_state(model, KNOWN_N4_STATE)
        self.assertTrue(np.all(coordinates > COORDINATE_LOWER))
        self.assertTrue(np.all(coordinates < COORDINATE_UPPER))
        reproduced = state_from_independent_rest_coordinates(
            model.parameters, coordinates
        )
        np.testing.assert_allclose(reproduced, KNOWN_N4_STATE, rtol=0.0, atol=2.0e-13)

    def test_local_ae4_expression_attempt_retains_full_rank_root(self) -> None:
        model = dynamic_model()
        full = np.r_[KNOWN_N4_STATE, model.regulation.basal_state]
        start = independent_coordinates_from_state(model, full)
        attempt, root = _one_attempt(
            model,
            transporter="AE4",
            expression=0.95,
            start_id="known_local_start",
            start=start,
            max_nfev=500,
        )
        self.assertTrue(attempt.optimizer_success, attempt.message)
        self.assertIsNotNone(root)
        self.assertTrue(root.passes_numerical_gate)
        self.assertEqual(root.jacobian_rank, 10)
        self.assertEqual(root.jacobian_nullity, 0)

    def test_intake_retains_every_passing_h1_root(self) -> None:
        group = [
            {
                "root_id": "R2",
                "native_dynamic_scientific_gate_pass": "True",
                "all_native_rest_wt_gate_pass": "True",
                "hydraulic_scale": "1",
                "hydraulic_mode": "H1",
            },
            {
                "root_id": "R1",
                "native_dynamic_scientific_gate_pass": "True",
                "all_native_rest_wt_gate_pass": "True",
                "hydraulic_scale": "1",
                "hydraulic_mode": "H1",
            },
            {
                "root_id": "REJECTED",
                "native_dynamic_scientific_gate_pass": "False",
                "all_native_rest_wt_gate_pass": "True",
                "hydraulic_scale": "1",
                "hydraulic_mode": "H1",
            },
        ]
        ids = _eligible_group_ids(group)
        self.assertEqual(ids, ("R1", "R2"))
        roots = [
            {
                "root_id": root_id,
                "passes_numerical_gate": "True",
                "passes_wt_gate": "True",
                "eligibility": "PRODUCTION_PRE_REVEAL",
                "hydraulic_scale": "1",
                "hydraulic_mode": "H1",
            }
            for root_id in ("R2", "R1")
        ]
        matched = _matching_root_rows(ids, roots)
        self.assertEqual(tuple(row["root_id"] for row in matched), ids)

    def test_string_boolean_lookalikes_and_rejected_eligibility_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "canonical case-sensitive"):
            _eligible_group_ids(
                [
                    {
                        "root_id": "R1",
                        "native_dynamic_scientific_gate_pass": "true",
                        "all_native_rest_wt_gate_pass": "True",
                        "hydraulic_scale": "1",
                        "hydraulic_mode": "H1",
                    }
                ]
            )
        with self.assertRaisesRegex(ValueError, "not production eligible"):
            _matching_root_rows(
                ("R1",),
                [
                    {
                        "root_id": "R1",
                        "passes_numerical_gate": "True",
                        "passes_wt_gate": "True",
                        "eligibility": "PRODUCTION_REJECTED",
                        "hydraulic_scale": "1",
                        "hydraulic_mode": "H1",
                    }
                ],
            )

    def test_dynamic_hash_ledger_requires_exact_seven_current_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            for index, filename in enumerate(DYNAMIC_FREEZE_ARTIFACTS):
                (output / filename).write_text(f"artifact {index}\n", encoding="utf-8")
            ledger = output / "native_dynamic_contract_hashes.csv"

            def write_ledger(rows):
                with ledger.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=("artifact", "sha256"))
                    writer.writeheader()
                    writer.writerows(rows)

            rows = [
                {"artifact": filename, "sha256": _sha256_file(output / filename)}
                for filename in DYNAMIC_FREEZE_ARTIFACTS
            ]
            write_ledger(rows)
            _validate_dynamic_hash_ledger(output)
            write_ledger((*rows, rows[0]))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                _validate_dynamic_hash_ledger(output)
            write_ledger(rows)
            (output / DYNAMIC_FREEZE_ARTIFACTS[0]).write_text(
                "stale\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "frozen manifest hash"):
                _validate_dynamic_hash_ledger(output)

    def test_exact_registry_rejects_missing_unexpected_and_stale_entries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "a").write_text("a\n", encoding="utf-8")
            registry = {"a": _sha256_file(output / "a")}
            _validate_exact_current_registry(
                registry, ("a",), base=output, label="test_registry"
            )
            with self.assertRaisesRegex(ValueError, "unexpected"):
                _validate_exact_current_registry(
                    {**registry, "extra": "0" * 64},
                    ("a",),
                    base=output,
                    label="test_registry",
                )
            (output / "a").write_text("changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "frozen manifest hash"):
                _validate_exact_current_registry(
                    registry, ("a",), base=output, label="test_registry"
                )

    def test_native_genotype_freeze_must_be_exact_and_current(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            output = base / "output"
            repository = base / "repository"
            output.mkdir()
            repository.mkdir()
            for filename in NATIVE_GENOTYPE_INPUTS:
                (output / filename).write_text(
                    f"native input {filename}\n", encoding="utf-8"
                )
            for filename in NATIVE_GENOTYPE_SOURCE_PATHS:
                path = repository / filename
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"source {filename}\n", encoding="utf-8")
            for filename in PRODUCTION_GENOTYPE_REPORTS:
                path = output / filename
                if filename == PRODUCTION_SUMMARY:
                    path.write_text(
                        json.dumps(
                            {
                                "all_passing_rows_retained": True,
                                "target_free_numerical_reproducibility_gate_pass": True,
                                "evaluated_root_ids": ["R1"],
                            }
                        ),
                        encoding="utf-8",
                    )
                else:
                    path.write_text(f"report {filename}\n", encoding="utf-8")
            source_hashes = {
                filename: _sha256_file(repository / filename)
                for filename in NATIVE_GENOTYPE_SOURCE_PATHS
            }
            manifest = {
                "manifest_id": "TASK13B_NATIVE_GENOTYPE_TARGET_FREE_V1",
                "status": "FROZEN_TARGET_FREE_GENOTYPE_PREDICTIONS",
                "independent_genotype_reproduction_required": True,
                "clean_auditor_authorization_required": True,
                "phenotype_reveal_permitted_by_this_manifest": False,
                "experimental_genotype_comparison_loaded": False,
                "experimental_acceptance_criterion_applied": False,
                "solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
                "root_count": 1,
                "roots": {"R1": {}},
                "input_artifact_sha256": {
                    filename: _sha256_file(output / filename)
                    for filename in NATIVE_GENOTYPE_INPUTS
                },
                "report_artifact_sha256": {
                    filename: _sha256_file(output / filename)
                    for filename in PRODUCTION_GENOTYPE_REPORTS
                },
                "equation_source_file_sha256": source_hashes,
                "equation_source_tree_sha256": _sha256_object(source_hashes),
            }
            (output / PRODUCTION_MANIFEST).write_text(
                json.dumps(manifest), encoding="utf-8"
            )
            with (output / PRODUCTION_HASH).open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=("artifact", "sha256"))
                writer.writeheader()
                writer.writerows(
                    {
                        "artifact": filename,
                        "sha256": _sha256_file(output / filename),
                    }
                    for filename in PRODUCTION_GENOTYPE_FREEZE_ARTIFACTS
                )
            _validate_native_genotype_freeze(output, repository, ("R1",))
            (output / PRODUCTION_PROFILE).write_text("stale\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "frozen manifest hash"):
                _validate_native_genotype_freeze(output, repository, ("R1",))

    def test_independent_snapshot_and_failed_publish_are_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            output = base / "output"
            repository = base / "repository"
            staging = base / "staging"
            output.mkdir()
            repository.mkdir()
            staging.mkdir()
            (output / "input").write_text("input-v1\n", encoding="utf-8")
            (repository / "source").write_text("source-v1\n", encoding="utf-8")
            with (
                patch.object(
                    independent_runner, "INDEPENDENT_REQUIRED_INPUTS", ("input",)
                ),
                patch.object(
                    independent_runner, "INDEPENDENT_SOURCE_PATHS", ("source",)
                ),
            ):
                snapshot = independent_runner._capture_independent_genotype_snapshot(
                    output, repository
                )
                independent_runner._assert_independent_genotype_snapshot_current(
                    snapshot, output, repository
                )
                (output / "input").write_text("input-v2\n", encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "changed during execution"):
                    independent_runner._assert_independent_genotype_snapshot_current(
                        snapshot, output, repository
                    )

            for filename in INDEPENDENT_GENOTYPE_FREEZE_ARTIFACTS:
                (staging / filename).write_text(filename, encoding="utf-8")
            (staging / HASH_OUTPUT).write_text("new ledger", encoding="utf-8")
            (output / HASH_OUTPUT).write_text("stale ledger", encoding="utf-8")
            snapshot = IndependentGenotypeDependencySnapshot({}, {})
            with patch.object(
                independent_runner,
                "_assert_independent_genotype_snapshot_current",
                side_effect=(None, RuntimeError("changed during commit")),
            ):
                with self.assertRaisesRegex(RuntimeError, "changed during commit"):
                    independent_runner._publish_staged_independent_genotype_freeze(
                        staging, output, repository, snapshot
                    )
            self.assertFalse((output / HASH_OUTPUT).exists())

    def test_opaque_firewall_registry_is_exact_but_never_opened(self) -> None:
        registry = {key: "0" * 64 for key in OPAQUE_FIREWALL_KEYS}
        _validate_opaque_registry(registry)
        with self.assertRaisesRegex(ValueError, "exact four-entry"):
            _validate_opaque_registry({**registry, "unexpected": "0" * 64})
        malformed = dict(registry)
        malformed[OPAQUE_FIREWALL_KEYS[0]] = "not-a-digest"
        with self.assertRaisesRegex(ValueError, "malformed"):
            _validate_opaque_registry(malformed)


if __name__ == "__main__":
    unittest.main()
