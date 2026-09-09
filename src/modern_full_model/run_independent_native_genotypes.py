"""Run independent target-free AE4/AE2 deletion checks on accepted H1 roots.

The intake is every root marked scientifically passing by the frozen native
WT dynamic group gate.  The script never selects the first passer.  It uses
the independent whole-cell transcription and independent 21-step resting
continuation, runs Radau and BDF only after a connected exact-zero branch is
found, and reports model-to-model ratios without phenotype scoring.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np

from .independent import (
    CORE_SIZE,
    IndependentN1NkccRegulation,
    IndependentSecretagogueProtocol,
    independent_physical_time_grid,
    independent_r3_regulation,
)
from .independent_genotype import (
    EXPRESSION_GRID,
    TRANSPORTERS,
    IndependentGenotypeContinuation,
    continue_independent_genotype_expression,
)
from .independent_native import IndependentNativeSourceDefinition


GROUP_INPUT = "native_dynamic_contract_group_gate.csv"
GATE_INPUT = "native_dynamic_contract_gate.json"
MANIFEST_INPUT = "native_dynamic_contract_manifest.json"
DYNAMIC_HASH_INPUT = "native_dynamic_contract_hashes.csv"
ROOT_INPUT = "native_source_wt_roots.csv"
NATIVE_SUMMARY_INPUT = "native_source_wt_summary.json"
PRODUCTION_CONTINUATION = "native_genotype_continuation.csv"
PRODUCTION_TRAJECTORY = "native_genotype_trajectories.csv"
PRODUCTION_PROFILE = "native_genotype_profile.csv"
PRODUCTION_SOLVER = "native_genotype_solver_crosschecks.csv"
PRODUCTION_SUMMARY = "native_genotype_summary.json"
PRODUCTION_MANIFEST = "native_genotype_manifest.json"
PRODUCTION_HASH = "native_genotype_hashes.csv"

CONTINUATION_OUTPUT = "independent_native_genotype_continuation.csv"
PROFILE_OUTPUT = "independent_native_genotype_profile.csv"
SOLVER_OUTPUT = "independent_native_genotype_solver_crosschecks.csv"
SUMMARY_OUTPUT = "independent_native_genotype_summary.json"
MANIFEST_OUTPUT = "independent_native_genotype_manifest.json"
HASH_OUTPUT = "independent_native_genotype_hashes.csv"

PRODUCTION_GENOTYPE_REPORTS = (
    PRODUCTION_CONTINUATION,
    PRODUCTION_TRAJECTORY,
    PRODUCTION_PROFILE,
    PRODUCTION_SOLVER,
    PRODUCTION_SUMMARY,
)
PRODUCTION_GENOTYPE_FREEZE_ARTIFACTS = (
    *PRODUCTION_GENOTYPE_REPORTS,
    PRODUCTION_MANIFEST,
)
INDEPENDENT_GENOTYPE_REPORTS = (
    CONTINUATION_OUTPUT,
    PROFILE_OUTPUT,
    SOLVER_OUTPUT,
    SUMMARY_OUTPUT,
)
INDEPENDENT_GENOTYPE_FREEZE_ARTIFACTS = (
    *INDEPENDENT_GENOTYPE_REPORTS,
    MANIFEST_OUTPUT,
)
GENOTYPE_SOLVER_RELATIVE_TOLERANCE = 1.0e-4
INDEPENDENT_REPRODUCTION_RELATIVE_TOLERANCE = 1.0e-4

WT_DYNAMIC_REPORTS = (
    "native_dynamic_contract_profile.csv",
    GROUP_INPUT,
    "native_dynamic_contract_solver_crosschecks.csv",
    "native_dynamic_contract_isomorphism.csv",
    "native_dynamic_contract_nearby.csv",
)
DYNAMIC_FREEZE_ARTIFACTS = (*WT_DYNAMIC_REPORTS, MANIFEST_INPUT, GATE_INPUT)
NATIVE_FREEZE_INPUTS = (
    "native_source_map.json",
    "native_source_root_attempts.csv",
    ROOT_INPUT,
    NATIVE_SUMMARY_INPUT,
    "native_source_panel.csv",
    "native_source_scale_one_identity.csv",
    "native_source_hashes.csv",
    "independent_native_roots.csv",
    "independent_native_attempts.csv",
    "independent_native_summary.json",
)
FROZEN_EQUATION_PATHS = tuple(
    f"src/modern_full_model/{name}"
    for name in (
        "__init__.py",
        "acid_base.py",
        "camp_pka.py",
        "calibration.py",
        "independent.py",
        "independent_native.py",
        "membranes.py",
        "model.py",
        "native_dynamic_contract.py",
        "native_source_panel.py",
        "nkcc_stimulation.py",
        "parameters.py",
        "run_native_dynamic_contract.py",
        "run_native_source_panel.py",
        "run_independent_native.py",
        "states.py",
        "transporters.py",
        "validation.py",
        "water.py",
    )
)
NATIVE_GENOTYPE_SOURCE_PATHS = (
    *FROZEN_EQUATION_PATHS,
    "src/modern_full_model/genotype_evaluation.py",
    "src/modern_full_model/run_native_genotype_evaluation.py",
)
NATIVE_GENOTYPE_INPUTS = (
    *DYNAMIC_FREEZE_ARTIFACTS,
    DYNAMIC_HASH_INPUT,
    *NATIVE_FREEZE_INPUTS,
)
INDEPENDENT_SOURCE_PATHS = (
    "src/modern_full_model/independent.py",
    "src/modern_full_model/independent_native.py",
    "src/modern_full_model/independent_genotype.py",
    "src/modern_full_model/run_independent_native_genotypes.py",
)
INDEPENDENT_REQUIRED_INPUTS = (
    *NATIVE_GENOTYPE_INPUTS,
    *PRODUCTION_GENOTYPE_FREEZE_ARTIFACTS,
    PRODUCTION_HASH,
)
OPAQUE_FIREWALL_KEYS = (
    "calibration_ledger_sha256",
    "validation_ledger_sha256",
    "sealed_phenotype_target_ledger_sha256",
    "pre_reveal_log_sha256",
)


@dataclass(frozen=True)
class IndependentGenotypeDependencySnapshot:
    input_artifact_sha256: Mapping[str, str]
    equation_source_file_sha256: Mapping[str, str]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _reject_duplicate_json_keys(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"duplicate JSON key {key!r}")
        value[key] = item
    return value


def _read_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_json_keys,
    )
    if not isinstance(value, Mapping):
        raise ValueError(f"{path.name} must contain a JSON object")
    return value


def _csv_bool(value: Any, *, field: str) -> bool:
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(
        f"{field} must be bool or canonical case-sensitive True/False CSV text; "
        f"found {value!r}"
    )


def _json_bool(
    payload: Mapping[str, Any], field: str, expected: bool, *, artifact: str
) -> None:
    if field not in payload or type(payload[field]) is not bool:
        raise ValueError(f"{artifact}.{field} must be a real JSON boolean")
    if payload[field] is not expected:
        raise ValueError(f"{artifact}.{field} must be {expected}")


def _json_nonnegative_integer(
    payload: Mapping[str, Any], field: str, *, artifact: str
) -> int:
    value = payload.get(field)
    if type(value) is not int or value < 0:
        raise ValueError(f"{artifact}.{field} must be a nonnegative JSON integer")
    return value


def _json_root_ids(
    payload: Mapping[str, Any], field: str, *, artifact: str
) -> tuple[str, ...]:
    values = payload.get(field)
    if not isinstance(values, list):
        raise ValueError(f"{artifact}.{field} must be a JSON list")
    if any(type(value) is not str or not value for value in values):
        raise ValueError(f"{artifact}.{field} contains an invalid root ID")
    if len(set(values)) != len(values):
        raise ValueError(f"{artifact}.{field} contains duplicate root IDs")
    return tuple(sorted(values))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_object(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _require_hash(
    registry: Mapping[str, Any], filename: str, actual: str, *, label: str
) -> None:
    expected = registry.get(filename)
    if (
        type(expected) is not str
        or len(expected) != 64
        or any(character not in "0123456789abcdef" for character in expected)
    ):
        raise ValueError(f"{label}.{filename} is not a lowercase SHA-256 digest")
    if expected != actual:
        raise ValueError(f"{filename} does not match its frozen manifest hash")


def _validate_exact_current_registry(
    registry: Any,
    expected_names: Sequence[str],
    *,
    base: Path,
    label: str,
) -> Mapping[str, str]:
    if not isinstance(registry, Mapping):
        raise ValueError(f"{label} must be a JSON object")
    expected = set(expected_names)
    actual = set(registry)
    if actual != expected:
        raise ValueError(
            f"{label} keys differ from the exact frozen set; "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    for filename in expected_names:
        path = base / filename
        if not path.is_file():
            raise FileNotFoundError(f"{label} binds missing artifact {filename}")
        _require_hash(
            registry,
            filename,
            _sha256_file(path),
            label=label,
        )
    return {name: str(registry[name]) for name in expected_names}


def _validate_opaque_registry(registry: Any) -> None:
    """Validate opaque bindings without opening any protected ledger."""

    if not isinstance(registry, Mapping):
        raise ValueError("opaque_firewall_ledger_sha256 must be a JSON object")
    expected = set(OPAQUE_FIREWALL_KEYS)
    actual = set(registry)
    if actual != expected:
        raise ValueError(
            "opaque firewall registry keys differ from the exact four-entry set; "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    for key in OPAQUE_FIREWALL_KEYS:
        value = registry[key]
        if (
            type(value) is not str
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
        ):
            raise ValueError(f"opaque firewall digest {key} is malformed")


def _validate_dynamic_hash_ledger(output: Path) -> None:
    path = output / DYNAMIC_HASH_INPUT
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != ("artifact", "sha256"):
            raise ValueError("dynamic hash ledger must have exact artifact,sha256 columns")
        rows = list(reader)
    artifacts = [row["artifact"] for row in rows]
    if len(set(artifacts)) != len(artifacts):
        raise ValueError("dynamic hash ledger contains duplicate artifacts")
    expected = set(DYNAMIC_FREEZE_ARTIFACTS)
    actual = set(artifacts)
    if actual != expected or len(rows) != len(DYNAMIC_FREEZE_ARTIFACTS):
        raise ValueError(
            "dynamic hash ledger must contain exactly seven frozen artifacts; "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    registry = {row["artifact"]: row["sha256"] for row in rows}
    _validate_exact_current_registry(
        registry,
        DYNAMIC_FREEZE_ARTIFACTS,
        base=output,
        label=DYNAMIC_HASH_INPUT,
    )


def _validate_exact_hash_ledger(
    path: Path,
    expected_artifacts: Sequence[str],
    *,
    base: Path,
    label: str,
) -> Mapping[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != ("artifact", "sha256"):
            raise ValueError(f"{label} must have exact artifact,sha256 columns")
        rows = list(reader)
    artifacts = [row["artifact"] for row in rows]
    if len(artifacts) != len(set(artifacts)):
        raise ValueError(f"{label} contains duplicate artifacts")
    expected = set(expected_artifacts)
    actual = set(artifacts)
    if actual != expected or len(rows) != len(expected_artifacts):
        raise ValueError(
            f"{label} must bind its exact artifact set; "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    registry = {row["artifact"]: row["sha256"] for row in rows}
    return _validate_exact_current_registry(
        registry,
        expected_artifacts,
        base=base,
        label=label,
    )


def _validate_native_genotype_freeze(
    output: Path,
    repository: Path,
    group_ids: Sequence[str],
) -> Mapping[str, Any]:
    """Require the exact current target-free native genotype publication."""

    _validate_exact_hash_ledger(
        output / PRODUCTION_HASH,
        PRODUCTION_GENOTYPE_FREEZE_ARTIFACTS,
        base=output,
        label=PRODUCTION_HASH,
    )
    manifest = _read_json(output / PRODUCTION_MANIFEST)
    if manifest.get("manifest_id") != "TASK13B_NATIVE_GENOTYPE_TARGET_FREE_V1":
        raise ValueError("native genotype manifest has unexpected manifest_id")
    if manifest.get("status") != "FROZEN_TARGET_FREE_GENOTYPE_PREDICTIONS":
        raise ValueError("native genotype manifest is not a frozen prediction set")
    for field, expected in (
        ("independent_genotype_reproduction_required", True),
        ("clean_auditor_authorization_required", True),
        ("phenotype_reveal_permitted_by_this_manifest", False),
        ("experimental_genotype_comparison_loaded", False),
        ("experimental_acceptance_criterion_applied", False),
    ):
        _json_bool(manifest, field, expected, artifact=PRODUCTION_MANIFEST)
    if float(manifest.get("solver_relative_tolerance", math.nan)) != (
        GENOTYPE_SOLVER_RELATIVE_TOLERANCE
    ):
        raise ValueError("native genotype solver tolerance differs from frozen contract")
    roots = manifest.get("roots")
    if not isinstance(roots, Mapping):
        raise ValueError("native genotype manifest lacks its root registry")
    expected_ids = tuple(sorted(group_ids))
    if tuple(sorted(str(root_id) for root_id in roots)) != expected_ids:
        raise ValueError("native genotype root registry differs from dynamic passing IDs")
    if _json_nonnegative_integer(
        manifest, "root_count", artifact=PRODUCTION_MANIFEST
    ) != len(expected_ids):
        raise ValueError("native genotype manifest root_count is stale")
    _validate_exact_current_registry(
        manifest.get("input_artifact_sha256"),
        NATIVE_GENOTYPE_INPUTS,
        base=output,
        label="native_genotype.input_artifact_sha256",
    )
    _validate_exact_current_registry(
        manifest.get("report_artifact_sha256"),
        PRODUCTION_GENOTYPE_REPORTS,
        base=output,
        label="native_genotype.report_artifact_sha256",
    )
    source_hashes = _validate_exact_current_registry(
        manifest.get("equation_source_file_sha256"),
        NATIVE_GENOTYPE_SOURCE_PATHS,
        base=repository,
        label="native_genotype.equation_source_file_sha256",
    )
    if manifest.get("equation_source_tree_sha256") != _sha256_object(
        dict(source_hashes)
    ):
        raise ValueError("native genotype equation source tree hash is stale")
    summary = _read_json(output / PRODUCTION_SUMMARY)
    _json_bool(
        summary,
        "all_passing_rows_retained",
        True,
        artifact=PRODUCTION_SUMMARY,
    )
    _json_bool(
        summary,
        "target_free_numerical_reproducibility_gate_pass",
        True,
        artifact=PRODUCTION_SUMMARY,
    )
    summary_ids = _json_root_ids(
        summary, "evaluated_root_ids", artifact=PRODUCTION_SUMMARY
    )
    if summary_ids != expected_ids:
        raise ValueError("native genotype summary root IDs are stale")
    return manifest


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _eligible_group_ids(rows: Sequence[Mapping[str, str]]) -> tuple[str, ...]:
    ids: list[str] = []
    all_ids: list[str] = []
    for row in rows:
        if "native_dynamic_scientific_gate_pass" not in row:
            raise ValueError("native dynamic group table lacks scientific gate")
        root_id = str(row.get("root_id", ""))
        if not root_id:
            raise ValueError("native dynamic group row lacks root_id")
        all_ids.append(root_id)
        if not _csv_bool(
            row["native_dynamic_scientific_gate_pass"],
            field="native_dynamic_scientific_gate_pass",
        ):
            continue
        if float(row.get("hydraulic_scale", "nan")) != 1.0:
            raise ValueError(f"{root_id}: independent genotype intake is not H=1")
        if str(row.get("hydraulic_mode", "")) != "H1":
            raise ValueError(f"{root_id}: independent genotype intake is not H1 mode")
        if "all_native_rest_wt_gate_pass" not in row or not _csv_bool(
            row["all_native_rest_wt_gate_pass"],
            field="all_native_rest_wt_gate_pass",
        ):
            raise ValueError(f"{root_id}: dynamic pass lacks native WT rest gate")
        ids.append(root_id)
    if len(set(all_ids)) != len(all_ids):
        raise ValueError("dynamic group IDs are not unique")
    return tuple(sorted(ids))


def _matching_root_rows(
    group_ids: Sequence[str], root_rows: Sequence[dict[str, str]]
) -> tuple[dict[str, str], ...]:
    registry = {row["root_id"]: row for row in root_rows}
    if len(registry) != len(root_rows):
        raise ValueError("native root table contains duplicate IDs")
    missing = tuple(root_id for root_id in group_ids if root_id not in registry)
    if missing:
        raise ValueError(f"passing roots absent from native root table: {missing}")
    matched: list[dict[str, str]] = []
    for root_id in group_ids:
        row = registry[root_id]
        if not _csv_bool(row.get("passes_numerical_gate"), field="passes_numerical_gate"):
            raise ValueError(f"{root_id}: native root fails numerical gate")
        if not _csv_bool(row.get("passes_wt_gate"), field="passes_wt_gate"):
            raise ValueError(f"{root_id}: native root fails WT gate")
        if str(row.get("eligibility", "")) != "PRODUCTION_PRE_REVEAL":
            raise ValueError(f"{root_id}: native root is not production eligible")
        if float(row.get("hydraulic_scale", "nan")) != 1.0:
            raise ValueError(f"{root_id}: native root is not H=1")
        if str(row.get("hydraulic_mode", "")) != "H1":
            raise ValueError(f"{root_id}: native root is not H1 mode")
        matched.append(row)
    return tuple(matched)


def _validate_frozen_intake(
    output: Path,
    group_rows: Sequence[Mapping[str, str]],
    group_ids: Sequence[str],
    root_rows: Sequence[dict[str, str]],
) -> None:
    """Require exact root-ID and hash agreement across the frozen WT gate."""

    gate = _read_json(output / GATE_INPUT)
    manifest = _read_json(output / MANIFEST_INPUT)
    if gate.get("gate_id") != "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1":
        raise ValueError("native dynamic gate has unexpected gate_id")
    if gate.get("status") != "SCIENTIFIC_GATE_PASS_AUDITOR_AUTHORIZATION_REQUIRED":
        raise RuntimeError("native dynamic composite scientific gate is not passing")
    if gate.get("manifest") != MANIFEST_INPUT:
        raise ValueError("native dynamic gate does not bind the expected manifest")
    for field, expected in (
        ("at_least_one_scientific_group_pass", True),
        ("all_rows_retained_no_first_passer", True),
        ("clean_auditor_authorization_required", True),
        ("independent_numerical_reproduction_required", True),
        ("clean_auditor_authorization_recorded", False),
        ("independent_numerical_reproduction_recorded_by_this_runner", False),
        ("final_composite_authorization_complete", False),
        ("phenotype_reveal_permitted", False),
    ):
        _json_bool(gate, field, expected, artifact=GATE_INPUT)
    ids = tuple(sorted(group_ids))
    if _json_root_ids(
        gate, "scientific_passing_root_ids", artifact=GATE_INPUT
    ) != ids:
        raise ValueError("dynamic gate and group CSV disagree on passing root IDs")
    if _json_nonnegative_integer(
        gate, "scientific_passing_group_count", artifact=GATE_INPUT
    ) != len(ids):
        raise ValueError("dynamic gate passing count disagrees with group CSV")
    if _json_nonnegative_integer(
        gate, "native_root_count", artifact=GATE_INPUT
    ) != len(group_rows):
        raise ValueError("dynamic gate native_root_count disagrees with group CSV")

    if manifest.get("manifest_id") != "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1":
        raise ValueError("native dynamic manifest has unexpected manifest_id")
    if manifest.get("status") != "FROZEN_FOR_CLEAN_ADVERSARIAL_REVIEW":
        raise ValueError("native dynamic manifest is not frozen for review")
    for field, expected in (
        ("clean_auditor_authorization_required", True),
        (
            "independent_numerical_reproduction_required_for_composite_authorization",
            True,
        ),
        ("phenotype_reveal_permitted_by_this_manifest", False),
    ):
        _json_bool(manifest, field, expected, artifact=MANIFEST_INPUT)
    if _json_root_ids(
        manifest, "scientific_passing_root_ids", artifact=MANIFEST_INPUT
    ) != ids:
        raise ValueError("dynamic manifest and group CSV disagree on passing root IDs")
    if _json_nonnegative_integer(
        manifest, "scientific_passing_group_count", artifact=MANIFEST_INPUT
    ) != len(ids):
        raise ValueError("dynamic manifest passing count disagrees with group CSV")
    manifest_roots = manifest.get("roots")
    if not isinstance(manifest_roots, Mapping):
        raise ValueError("dynamic manifest lacks its frozen root registry")
    all_group_ids = tuple(sorted(str(row["root_id"]) for row in group_rows))
    if tuple(sorted(manifest_roots)) != all_group_ids:
        raise ValueError("dynamic manifest root registry and group CSV IDs differ")
    if _json_nonnegative_integer(
        manifest, "root_count", artifact=MANIFEST_INPUT
    ) != len(manifest_roots):
        raise ValueError("dynamic manifest root_count is stale")

    repository = Path(__file__).resolve().parents[2]
    _validate_exact_current_registry(
        manifest.get("wt_dynamic_report_artifact_sha256"),
        WT_DYNAMIC_REPORTS,
        base=output,
        label="wt_dynamic_report_artifact_sha256",
    )
    _validate_exact_current_registry(
        manifest.get("native_input_artifact_sha256"),
        NATIVE_FREEZE_INPUTS,
        base=output,
        label="native_input_artifact_sha256",
    )
    equation_hashes = _validate_exact_current_registry(
        manifest.get("equation_source_file_sha256"),
        FROZEN_EQUATION_PATHS,
        base=repository,
        label="equation_source_file_sha256",
    )
    if manifest.get("equation_source_tree_sha256") != _sha256_object(
        dict(equation_hashes)
    ):
        raise ValueError("dynamic manifest equation source tree hash is stale")
    _validate_opaque_registry(manifest.get("opaque_firewall_ledger_sha256"))
    _validate_dynamic_hash_ledger(output)
    row_registry = {row["root_id"]: row for row in root_rows}
    for root_id in ids:
        payload = manifest_roots.get(root_id)
        if not isinstance(payload, Mapping):
            raise ValueError(f"{root_id}: dynamic manifest root payload is malformed")
        core = json.loads(row_registry[root_id]["core_state_json"])
        expected = payload.get("core_state_sha256")
        if type(expected) is not str or expected != _sha256_object(core):
            raise ValueError(f"{root_id}: frozen core-state hash mismatch")


def _capture_independent_genotype_snapshot(
    output: Path, repository: Path
) -> IndependentGenotypeDependencySnapshot:
    missing_inputs = tuple(
        name for name in INDEPENDENT_REQUIRED_INPUTS
        if not (output / name).is_file()
    )
    missing_sources = tuple(
        name for name in INDEPENDENT_SOURCE_PATHS
        if not (repository / name).is_file()
    )
    if missing_inputs or missing_sources:
        raise FileNotFoundError(
            "independent genotype snapshot is incomplete; "
            f"missing_inputs={missing_inputs}, missing_sources={missing_sources}"
        )
    return IndependentGenotypeDependencySnapshot(
        input_artifact_sha256={
            name: _sha256_file(output / name)
            for name in INDEPENDENT_REQUIRED_INPUTS
        },
        equation_source_file_sha256={
            name: _sha256_file(repository / name)
            for name in INDEPENDENT_SOURCE_PATHS
        },
    )


def _assert_independent_genotype_snapshot_current(
    snapshot: IndependentGenotypeDependencySnapshot,
    output: Path,
    repository: Path,
) -> None:
    current = _capture_independent_genotype_snapshot(output, repository)
    if current != snapshot:
        changed: list[str] = []
        for label in ("input_artifact_sha256", "equation_source_file_sha256"):
            before = getattr(snapshot, label)
            after = getattr(current, label)
            changed.extend(
                f"{label}:{key}"
                for key in sorted(set(before) | set(after))
                if before.get(key) != after.get(key)
            )
        raise RuntimeError(
            "independent genotype dependencies changed during execution: "
            f"{changed}"
        )


def _dynamic_model(definition: IndependentNativeSourceDefinition) -> Any:
    regulation = independent_r3_regulation(
        tau_pka_s=10.0,
        forward_regulation_rate_s=1.0 / 30.0,
        reverse_regulation_rate_s=1.0 / 30.0,
        basal_multiplier=1.0,
        fully_activated_increment=0.25,
    )
    nkcc = IndependentN1NkccRegulation(
        fully_activated_multiplier=1.75,
        resting_calcium_uM=0.058,
        stimulated_calcium_uM=0.10,
    )
    protocol = IndependentSecretagogueProtocol(
        resting_calcium_uM=0.058,
        stimulated_calcium_uM=0.10,
        beta_occupancy_on=1.0,
    )
    return definition.build_independent_model(
        regulation=regulation,
        nkcc_regulation=nkcc,
        stimulus=protocol,
    )


def _initial_state(model: Any, core: Sequence[float]) -> np.ndarray:
    ae4 = () if model.regulation is None else model.regulation.basal_state
    nkcc = () if model.nkcc_regulation is None else model.nkcc_regulation.basal_state
    return np.r_[np.asarray(core, dtype=float), np.asarray((*ae4, *nkcc), dtype=float)]


def _atol(state_names: Sequence[str]) -> np.ndarray:
    return np.asarray(
        [
            1.0e-10
            if name.endswith("_fmol")
            else 1.0e-12
            if name.endswith("_pL")
            else 1.0e-10
            for name in state_names
        ],
        dtype=float,
    )


def _expression_kwargs(transporter: str, expression: float) -> dict[str, float]:
    return (
        {"ae4_expression": float(expression)}
        if transporter == "AE4"
        else {"ae2_expression": float(expression)}
    )


def _deleted_flux_magnitude(evaluation: Any, transporter: str) -> float:
    if transporter == "AE4":
        return float(
            max(
                abs(float(evaluation.fluxes_fmol_s["ae4_cl"])),
                abs(float(evaluation.conservation_residuals["ae4_charge_fmol_s"])),
            )
        )
    return float(abs(float(evaluation.fluxes_fmol_s["ae2"])))


def _endpoint(model: Any, trajectory: Any, transporter: str, expression: float) -> dict[str, float]:
    evaluation = model.evaluate(
        float(trajectory.time_s[-1]),
        trajectory.states[:, -1],
        **_expression_kwargs(transporter, expression),
    )
    state = trajectory.states[:, -1]
    return {
        "flow_pL_s": float(trajectory.flow_pL_s[-1]),
        "cumulative_flow_pL": float(trajectory.cumulative_flow_pL[-1]),
        "cell_na_mM": float(evaluation.cell_concentrations_mM["na"]),
        "cell_k_mM": float(evaluation.cell_concentrations_mM["k"]),
        "cell_cl_mM": float(evaluation.cell_concentrations_mM["cl"]),
        "cell_ph": float(evaluation.cell_speciation.ph),
        "cell_volume_pL": float(state[5]),
        "lumen_na_mM": float(evaluation.lumen_concentrations_mM["na"]),
        "lumen_k_mM": float(evaluation.lumen_concentrations_mM["k"]),
        "lumen_cl_mM": float(evaluation.lumen_concentrations_mM["cl"]),
        "lumen_ph": float(evaluation.lumen_speciation.ph),
        "lumen_volume_pL": float(state[11]),
    }


def _landmark_ratios(wt: Any, genotype: Any) -> tuple[dict[str, float], dict[str, float]]:
    cumulative: dict[str, float] = {}
    flow: dict[str, float] = {}
    for landmark in (60.0, 120.0, 180.0, 600.0):
        key = f"{landmark:g}_s"
        wt_cumulative = float(np.interp(landmark, wt.time_s, wt.cumulative_flow_pL))
        wt_flow = float(np.interp(landmark, wt.time_s, wt.flow_pL_s))
        if wt_cumulative <= 0.0 or wt_flow <= 0.0:
            raise ValueError(f"{key}: WT denominator must be positive")
        cumulative[key] = float(
            np.interp(landmark, genotype.time_s, genotype.cumulative_flow_pL)
            / wt_cumulative
        )
        flow[key] = float(
            np.interp(landmark, genotype.time_s, genotype.flow_pL_s) / wt_flow
        )
    return cumulative, flow


def _trajectory_gate(
    model: Any, trajectory: Any, transporter: str, expression: float
) -> tuple[bool, float]:
    deleted_max = 0.0
    for index, time_s in enumerate(trajectory.time_s):
        evaluation = model.evaluate(
            float(time_s),
            trajectory.states[:, index],
            **_expression_kwargs(transporter, expression),
        )
        if expression == 0.0:
            deleted_max = max(
                deleted_max, _deleted_flux_magnitude(evaluation, transporter)
            )
    gate = bool(
        trajectory.success
        and np.all(trajectory.states[:CORE_SIZE, :] > 0.0)
        and np.all(np.isfinite(trajectory.flow_pL_s))
        and np.all(trajectory.flow_pL_s >= 0.0)
        and trajectory.max_abs_conservation_fmol_s <= 1.0e-10
        and trajectory.max_abs_water_accounting_pL_s <= 1.0e-12
        and (expression != 0.0 or deleted_max == 0.0)
    )
    return gate, float(deleted_max)


def _relative_difference(left: np.ndarray, right: np.ndarray, floor: float) -> float:
    scale = np.maximum(np.maximum(np.abs(left), np.abs(right)), floor)
    return float(np.max(np.abs(left - right) / scale))


def _continuation_rows(
    root_id: str,
    transporter: str,
    report: IndependentGenotypeContinuation,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, step in enumerate(report.steps):
        selected = step.selected_root
        rows.append(
            {
                "root_id": root_id,
                "transporter": transporter,
                "expression_step_index": index,
                "expression": step.expression,
                "attempt_count": len(step.attempts),
                "attempts_json": _json([asdict(item) for item in step.attempts]),
                "retained_root_count": len(step.roots),
                "selected_root_id": step.selected_root_id,
                "selected_distance_from_previous_normalized": (
                    step.selected_distance_from_previous_normalized
                ),
                "connected_to_previous": step.connected_to_previous,
                "step_status": step.status,
                "selected_coordinates_json": (
                    _json(selected.coordinates) if selected else None
                ),
                "selected_core_state_json": (
                    _json(selected.core_state) if selected else None
                ),
                "selected_max_abs_scaled_residual": (
                    selected.max_abs_scaled_residual if selected else None
                ),
                "selected_max_abs_omitted_charge_rhs": (
                    selected.max_abs_omitted_charge_rhs_fmol_s if selected else None
                ),
                "selected_max_abs_charge_residual_fmol": (
                    selected.max_abs_bulk_charge_residual_fmol if selected else None
                ),
                "selected_max_abs_current_residual_A": (
                    selected.max_abs_current_residual_A if selected else None
                ),
                "selected_jacobian_rank": selected.jacobian_rank if selected else None,
                "selected_boundary_hits_json": (
                    _json(selected.boundary_hits) if selected else None
                ),
                "selected_exact_deleted_transport_zero": (
                    selected.exact_deleted_transport_zero if selected else None
                ),
                "continuation_completed_to_exact_zero": report.completed_to_exact_zero,
                "continuation_failure": report.branch_failure,
                "runner_exception": False,
            }
        )
    return rows


def _failed_continuation_row(
    root_id: str,
    transporter: str,
    failure: str,
) -> dict[str, Any]:
    """Return one schema-complete record for an exceptional continuation."""

    return {
        "root_id": root_id,
        "transporter": transporter,
        "expression_step_index": None,
        "expression": None,
        "attempt_count": 0,
        "attempts_json": _json(()),
        "retained_root_count": 0,
        "selected_root_id": None,
        "selected_distance_from_previous_normalized": None,
        "connected_to_previous": False,
        "step_status": "RUNNER_EXCEPTION",
        "selected_coordinates_json": None,
        "selected_core_state_json": None,
        "selected_max_abs_scaled_residual": None,
        "selected_max_abs_omitted_charge_rhs": None,
        "selected_max_abs_charge_residual_fmol": None,
        "selected_max_abs_current_residual_A": None,
        "selected_jacobian_rank": None,
        "selected_boundary_hits_json": None,
        "selected_exact_deleted_transport_zero": None,
        "continuation_completed_to_exact_zero": False,
        "continuation_failure": failure,
        "runner_exception": True,
    }


def _nonexecuted_dynamic_rows(
    root_id: str,
    transporter: str,
    failure: str | None,
    *,
    continuation_completed: bool,
    runner_exception: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Retain a complete profile/solver pair when dynamics cannot execute."""

    profile = {
        "root_id": root_id,
        "transporter": transporter,
        "solver_label": None,
        "solver_method": None,
        "dynamic_executed": False,
        "continuation_completed_to_exact_zero": continuation_completed,
        "continuation_failure": failure,
        "paired_numerical_gate_pass": None,
        "cumulative_ratio_at_landmarks_json": None,
        "flow_ratio_at_landmarks_json": None,
        "wt_endpoint_json": None,
        "genotype_endpoint_json": None,
        "genotype_to_wt_endpoint_json": None,
        "wt_max_abs_conservation_fmol_s": None,
        "genotype_max_abs_conservation_fmol_s": None,
        "wt_max_abs_water_accounting_pL_s": None,
        "genotype_max_abs_water_accounting_pL_s": None,
        "max_abs_deleted_transport_flux_fmol_s": None,
        "runner_exception": runner_exception,
    }
    solver = {
        "root_id": root_id,
        "transporter": transporter,
        "dynamic_executed": False,
        "continuation_completed_to_exact_zero": continuation_completed,
        "continuation_failure": failure,
        "both_solver_pairs_pass": None,
        "cross_solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        "cross_solver_gate_pass": None,
        "max_relative_difference_wt_state": None,
        "max_relative_difference_wt_flow": None,
        "max_relative_difference_wt_cumulative": None,
        "max_relative_difference_genotype_state": None,
        "max_relative_difference_genotype_flow": None,
        "max_relative_difference_genotype_cumulative": None,
        "runner_exception": runner_exception,
    }
    return profile, solver


def _run_one_root(row: dict[str, str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    root_id = row["root_id"]
    definition = IndependentNativeSourceDefinition.from_root_row(row)
    model = _dynamic_model(definition)
    wt_core = tuple(float(value) for value in json.loads(row["core_state_json"]))
    grid = independent_physical_time_grid(600.0, 1.0)
    wt_initial = _initial_state(model, wt_core)
    wt_by_solver: dict[str, Any] = {}
    continuation_rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []
    solver_rows: list[dict[str, Any]] = []

    for transporter in TRANSPORTERS:
        try:
            report = continue_independent_genotype_expression(
                model, wt_core, transporter=transporter
            )
        except Exception as exc:
            failure = f"CONTINUATION_EXCEPTION:{type(exc).__name__}: {exc}"
            continuation_rows.append(
                _failed_continuation_row(root_id, transporter, failure)
            )
            profile, solver = _nonexecuted_dynamic_rows(
                root_id,
                transporter,
                failure,
                continuation_completed=False,
                runner_exception=True,
            )
            profile_rows.append(profile)
            solver_rows.append(solver)
            continue
        continuation_rows.extend(_continuation_rows(root_id, transporter, report))
        null_root = report.connected_null_root
        if null_root is None:
            profile, solver = _nonexecuted_dynamic_rows(
                root_id,
                transporter,
                report.branch_failure,
                continuation_completed=False,
                runner_exception=False,
            )
            profile_rows.append(profile)
            solver_rows.append(solver)
            continue
        null_initial = _initial_state(model, null_root.core_state)
        solver_trajectories: dict[str, tuple[Any, Any]] = {}
        branch_profile_rows: list[dict[str, Any]] = []
        try:
            for method in ("Radau", "BDF"):
                if method not in wt_by_solver:
                    wt_by_solver[method] = model.integrate(
                        wt_initial,
                        grid,
                        method=method,
                        rtol=1.0e-7,
                        atol=_atol(model.state_names),
                        max_step_s=2.0,
                        preserve_basal_zero=True,
                    )
                wt = wt_by_solver[method]
                genotype = model.integrate(
                    null_initial,
                    grid,
                    method=method,
                    rtol=1.0e-7,
                    atol=_atol(model.state_names),
                    max_step_s=2.0,
                    preserve_basal_zero=True,
                    **_expression_kwargs(transporter, 0.0),
                )
                solver_trajectories[method] = (wt, genotype)
                wt_gate, _ = _trajectory_gate(model, wt, transporter, 1.0)
                genotype_gate, deleted_max = _trajectory_gate(
                    model, genotype, transporter, 0.0
                )
                cumulative_ratios, flow_ratios = _landmark_ratios(wt, genotype)
                wt_endpoint = _endpoint(model, wt, transporter, 1.0)
                genotype_endpoint = _endpoint(model, genotype, transporter, 0.0)
                endpoint_comparison = {
                    key: genotype_endpoint[key] / wt_endpoint[key]
                    for key in wt_endpoint
                    if key not in {"cell_ph", "lumen_ph"}
                }
                endpoint_comparison["cell_ph_difference"] = (
                    genotype_endpoint["cell_ph"] - wt_endpoint["cell_ph"]
                )
                endpoint_comparison["lumen_ph_difference"] = (
                    genotype_endpoint["lumen_ph"] - wt_endpoint["lumen_ph"]
                )
                branch_profile_rows.append(
                    {
                        "root_id": root_id,
                        "transporter": transporter,
                        "solver_label": f"production_{method.lower()}",
                        "solver_method": method,
                        "dynamic_executed": True,
                        "continuation_completed_to_exact_zero": True,
                        "continuation_failure": None,
                        "paired_numerical_gate_pass": bool(wt_gate and genotype_gate),
                        "cumulative_ratio_at_landmarks_json": _json(cumulative_ratios),
                        "flow_ratio_at_landmarks_json": _json(flow_ratios),
                        "wt_endpoint_json": _json(wt_endpoint),
                        "genotype_endpoint_json": _json(genotype_endpoint),
                        "genotype_to_wt_endpoint_json": _json(endpoint_comparison),
                        "wt_max_abs_conservation_fmol_s": (
                            wt.max_abs_conservation_fmol_s
                        ),
                        "genotype_max_abs_conservation_fmol_s": (
                            genotype.max_abs_conservation_fmol_s
                        ),
                        "wt_max_abs_water_accounting_pL_s": (
                            wt.max_abs_water_accounting_pL_s
                        ),
                        "genotype_max_abs_water_accounting_pL_s": (
                            genotype.max_abs_water_accounting_pL_s
                        ),
                        "max_abs_deleted_transport_flux_fmol_s": deleted_max,
                        "runner_exception": False,
                    }
                )
        except Exception as exc:
            failure = f"DYNAMIC_EXCEPTION:{type(exc).__name__}: {exc}"
            profile, solver = _nonexecuted_dynamic_rows(
                root_id,
                transporter,
                failure,
                continuation_completed=True,
                runner_exception=True,
            )
            profile_rows.append(profile)
            solver_rows.append(solver)
            continue
        radau_wt, radau_null = solver_trajectories["Radau"]
        bdf_wt, bdf_null = solver_trajectories["BDF"]
        differences = {
            "wt_state": _relative_difference(
                radau_wt.states, bdf_wt.states, 1.0e-10
            ),
            "wt_flow": _relative_difference(
                radau_wt.flow_pL_s, bdf_wt.flow_pL_s, 1.0e-12
            ),
            "wt_cumulative": _relative_difference(
                radau_wt.cumulative_flow_pL,
                bdf_wt.cumulative_flow_pL,
                1.0e-12,
            ),
            "genotype_state": _relative_difference(
                radau_null.states, bdf_null.states, 1.0e-10
            ),
            "genotype_flow": _relative_difference(
                radau_null.flow_pL_s, bdf_null.flow_pL_s, 1.0e-12
            ),
            "genotype_cumulative": _relative_difference(
                radau_null.cumulative_flow_pL,
                bdf_null.cumulative_flow_pL,
                1.0e-12,
            ),
        }
        cross_solver_gate = bool(
            all(
                math.isfinite(value)
                and value <= GENOTYPE_SOLVER_RELATIVE_TOLERANCE
                for value in differences.values()
            )
            and all(
                bool(row.get("paired_numerical_gate_pass"))
                for row in branch_profile_rows
            )
        )
        profile_rows.extend(branch_profile_rows)
        solver_rows.append(
            {
                "root_id": root_id,
                "transporter": transporter,
                "dynamic_executed": True,
                "continuation_completed_to_exact_zero": True,
                "continuation_failure": None,
                "cross_solver_relative_tolerance": (
                    GENOTYPE_SOLVER_RELATIVE_TOLERANCE
                ),
                "cross_solver_gate_pass": cross_solver_gate,
                "both_solver_pairs_pass": cross_solver_gate,
                **{
                    f"max_relative_difference_{name}": value
                    for name, value in differences.items()
                },
                "runner_exception": False,
            }
        )
    return continuation_rows, profile_rows, solver_rows


def _failed_root_rows(
    row: Mapping[str, str], exc: Exception
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Retain a fail-closed record for both branches after a root-level exception."""

    root_id = str(row.get("root_id", ""))
    failure = f"ROOT_RUNNER_EXCEPTION:{type(exc).__name__}: {exc}"
    continuation = [
        _failed_continuation_row(root_id, transporter, failure)
        for transporter in TRANSPORTERS
    ]
    profile: list[dict[str, Any]] = []
    solvers: list[dict[str, Any]] = []
    for transporter in TRANSPORTERS:
        branch_profile, branch_solver = _nonexecuted_dynamic_rows(
            root_id,
            transporter,
            failure,
            continuation_completed=False,
            runner_exception=True,
        )
        profile.append(branch_profile)
        solvers.append(branch_solver)
    return continuation, profile, solvers


def _production_discrepancies(
    output: Path,
    continuation_rows: Sequence[dict[str, Any]],
    profile_rows: Sequence[dict[str, Any]],
    solver_rows: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    report: dict[str, Any] = {
        "production_artifacts_present": True,
        "relative_tolerance": INDEPENDENT_REPRODUCTION_RELATIVE_TOLERANCE,
        "continuation_key_discrepancies": [],
        "profile_key_discrepancies": [],
        "solver_key_discrepancies": [],
        "max_selected_core_state_absolute_difference": None,
        "max_selected_core_state_relative_difference": None,
        "max_cumulative_ratio_absolute_difference": None,
        "max_cumulative_ratio_relative_difference": None,
        "max_flow_ratio_absolute_difference": None,
        "max_flow_ratio_relative_difference": None,
        "comparison_gate_pass": False,
    }
    continuation_path = output / PRODUCTION_CONTINUATION
    profile_path = output / PRODUCTION_PROFILE
    solver_path = output / PRODUCTION_SOLVER
    if not all(path.is_file() for path in (continuation_path, profile_path, solver_path)):
        raise FileNotFoundError("native genotype comparison artifacts are incomplete")

    def optional_bool(row: Mapping[str, Any], field: str) -> bool | None:
        value = row.get(field)
        if value is None or value == "":
            return None
        return _csv_bool(value, field=field)

    def key_registry(
        rows: Sequence[Mapping[str, Any]],
        fields: Sequence[str],
        *,
        label: str,
    ) -> dict[tuple[str, ...], Mapping[str, Any]]:
        registry: dict[tuple[str, ...], Mapping[str, Any]] = {}
        for row in rows:
            key = tuple(
                "" if row.get(field) is None else str(row.get(field, ""))
                for field in fields
            )
            if key in registry:
                raise ValueError(f"{label} contains duplicate key {key}")
            registry[key] = row
        return registry

    def relative_difference(left: float, right: float, floor: float = 1.0e-12) -> float:
        values = (float(left), float(right))
        if not all(math.isfinite(value) for value in values):
            raise ValueError("production comparison values must be finite")
        return abs(values[0] - values[1]) / max(
            abs(values[0]), abs(values[1]), floor
        )

    production_continuation = _read_csv(continuation_path)
    continuation_fields = ("root_id", "transporter", "expression_step_index")
    production_continuation_registry = key_registry(
        production_continuation,
        continuation_fields,
        label=PRODUCTION_CONTINUATION,
    )
    independent_continuation = key_registry(
        continuation_rows,
        continuation_fields,
        label=CONTINUATION_OUTPUT,
    )
    production_keys = set(production_continuation_registry)
    independent_keys = set(independent_continuation)
    for key in sorted(production_keys - independent_keys):
        report["continuation_key_discrepancies"].append(
            {"key": key, "issue": "missing independent key"}
        )
    for key in sorted(independent_keys - production_keys):
        report["continuation_key_discrepancies"].append(
            {"key": key, "issue": "unexpected independent key"}
        )
    state_differences: list[float] = []
    state_relative_differences: list[float] = []
    for key in sorted(production_keys & independent_keys):
        row = production_continuation_registry[key]
        independent = independent_continuation[key]
        for field in (
            "continuation_completed_to_exact_zero",
            "connected_to_previous",
            "runner_exception",
        ):
            if optional_bool(row, field) != optional_bool(independent, field):
                report["continuation_key_discrepancies"].append(
                    {"key": key, "issue": f"{field} differs"}
                )
        production_state_present = bool(row.get("selected_core_state_json"))
        independent_state_present = bool(independent.get("selected_core_state_json"))
        if production_state_present != independent_state_present:
            report["continuation_key_discrepancies"].append(
                {"key": key, "issue": "selected core-state presence differs"}
            )
        elif production_state_present:
            production_state = np.asarray(
                json.loads(row["selected_core_state_json"]), dtype=float
            )
            independent_state = np.asarray(
                json.loads(independent["selected_core_state_json"]), dtype=float
            )
            if (
                production_state.shape != independent_state.shape
                or np.any(~np.isfinite(production_state))
                or np.any(~np.isfinite(independent_state))
            ):
                report["continuation_key_discrepancies"].append(
                    {"key": key, "issue": "selected core states are not finite/aligned"}
                )
            else:
                absolute = float(np.max(np.abs(production_state - independent_state)))
                scale = np.maximum(
                    np.maximum(np.abs(production_state), np.abs(independent_state)),
                    1.0e-10,
                )
                relative = float(
                    np.max(np.abs(production_state - independent_state) / scale)
                )
                state_differences.append(absolute)
                state_relative_differences.append(relative)
                if relative > INDEPENDENT_REPRODUCTION_RELATIVE_TOLERANCE:
                    report["continuation_key_discrepancies"].append(
                        {"key": key, "issue": "selected core states exceed tolerance"}
                    )

    production_profile = _read_csv(profile_path)
    profile_fields = ("root_id", "transporter", "solver_label")
    production_profile_registry = key_registry(
        production_profile, profile_fields, label=PRODUCTION_PROFILE
    )
    independent_profile = key_registry(
        profile_rows, profile_fields, label=PROFILE_OUTPUT
    )
    production_profile_keys = set(production_profile_registry)
    independent_profile_keys = set(independent_profile)
    for key in sorted(production_profile_keys - independent_profile_keys):
        report["profile_key_discrepancies"].append(
            {"key": key, "issue": "missing independent key"}
        )
    for key in sorted(independent_profile_keys - production_profile_keys):
        report["profile_key_discrepancies"].append(
            {"key": key, "issue": "unexpected independent key"}
        )
    cumulative_differences: list[float] = []
    cumulative_relative_differences: list[float] = []
    flow_differences: list[float] = []
    flow_relative_differences: list[float] = []
    for key in sorted(production_profile_keys & independent_profile_keys):
        row = production_profile_registry[key]
        independent = independent_profile[key]
        for field in (
            "dynamic_executed",
            "continuation_completed_to_exact_zero",
            "paired_numerical_gate_pass",
            "runner_exception",
        ):
            if optional_bool(row, field) != optional_bool(independent, field):
                report["profile_key_discrepancies"].append(
                    {"key": key, "issue": f"{field} differs"}
                )
        for field, absolute_destination, relative_destination in (
            (
                "cumulative_ratio_at_landmarks_json",
                cumulative_differences,
                cumulative_relative_differences,
            ),
            (
                "flow_ratio_at_landmarks_json",
                flow_differences,
                flow_relative_differences,
            ),
        ):
            production_present = bool(row.get(field))
            independent_present = bool(independent.get(field))
            if production_present != independent_present:
                report["profile_key_discrepancies"].append(
                    {"key": key, "issue": f"{field} presence differs"}
                )
                continue
            if not production_present:
                continue
            production_values = json.loads(row[field])
            independent_values = json.loads(independent[field])
            if set(production_values) != set(independent_values):
                report["profile_key_discrepancies"].append(
                    {"key": key, "issue": f"{field} landmark keys differ"}
                )
                continue
            for name in production_values:
                left = float(production_values[name])
                right = float(independent_values[name])
                absolute_destination.append(abs(left - right))
                relative = relative_difference(left, right)
                relative_destination.append(relative)
                if relative > INDEPENDENT_REPRODUCTION_RELATIVE_TOLERANCE:
                    report["profile_key_discrepancies"].append(
                        {"key": key, "issue": f"{field}.{name} exceeds tolerance"}
                    )

    production_solver = _read_csv(solver_path)
    solver_fields = ("root_id", "transporter")
    production_solver_registry = key_registry(
        production_solver, solver_fields, label=PRODUCTION_SOLVER
    )
    independent_solver = key_registry(
        solver_rows, solver_fields, label=SOLVER_OUTPUT
    )
    production_solver_keys = set(production_solver_registry)
    independent_solver_keys = set(independent_solver)
    for key in sorted(production_solver_keys - independent_solver_keys):
        report["solver_key_discrepancies"].append(
            {"key": key, "issue": "missing independent key"}
        )
    for key in sorted(independent_solver_keys - production_solver_keys):
        report["solver_key_discrepancies"].append(
            {"key": key, "issue": "unexpected independent key"}
        )
    for key in sorted(production_solver_keys & independent_solver_keys):
        row = production_solver_registry[key]
        independent = independent_solver[key]
        for field in (
            "dynamic_executed",
            "continuation_completed_to_exact_zero",
            "cross_solver_gate_pass",
            "runner_exception",
        ):
            if optional_bool(row, field) != optional_bool(independent, field):
                report["solver_key_discrepancies"].append(
                    {"key": key, "issue": f"{field} differs"}
                )

    report["max_selected_core_state_absolute_difference"] = (
        max(state_differences) if state_differences else None
    )
    report["max_selected_core_state_relative_difference"] = (
        max(state_relative_differences) if state_relative_differences else None
    )
    report["max_cumulative_ratio_absolute_difference"] = (
        max(cumulative_differences) if cumulative_differences else None
    )
    report["max_cumulative_ratio_relative_difference"] = (
        max(cumulative_relative_differences)
        if cumulative_relative_differences
        else None
    )
    report["max_flow_ratio_absolute_difference"] = (
        max(flow_differences) if flow_differences else None
    )
    report["max_flow_ratio_relative_difference"] = (
        max(flow_relative_differences) if flow_relative_differences else None
    )
    report["comparison_gate_pass"] = not any(
        report[field]
        for field in (
            "continuation_key_discrepancies",
            "profile_key_discrepancies",
            "solver_key_discrepancies",
        )
    )
    return report


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    names: list[str] = []
    for row in rows:
        for name in row:
            if name not in names:
                names.append(name)
    if not names:
        names = list(fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _invalidate_independent_genotype_hash_ledger(output: Path) -> None:
    (output / HASH_OUTPUT).unlink(missing_ok=True)


def _publish_staged_independent_genotype_freeze(
    staging: Path,
    output: Path,
    repository: Path,
    snapshot: IndependentGenotypeDependencySnapshot,
) -> None:
    _invalidate_independent_genotype_hash_ledger(output)
    _assert_independent_genotype_snapshot_current(snapshot, output, repository)
    for filename in INDEPENDENT_GENOTYPE_FREEZE_ARTIFACTS:
        os.replace(staging / filename, output / filename)
    _assert_independent_genotype_snapshot_current(snapshot, output, repository)
    os.replace(staging / HASH_OUTPUT, output / HASH_OUTPUT)
    try:
        _assert_independent_genotype_snapshot_current(snapshot, output, repository)
    except Exception:
        _invalidate_independent_genotype_hash_ledger(output)
        raise


def run_independent_native_genotypes(output: Path) -> dict[str, Any]:
    repository = Path(__file__).resolve().parents[2]
    required = INDEPENDENT_REQUIRED_INPUTS
    missing = tuple(filename for filename in required if not (output / filename).is_file())
    if missing:
        raise RuntimeError(
            "independent genotype reproduction requires exact frozen native dynamic "
            f"and native genotype publications; missing {missing}"
        )
    snapshot = _capture_independent_genotype_snapshot(output, repository)
    group_rows = _read_csv(output / GROUP_INPUT)
    root_table = _read_csv(output / ROOT_INPUT)
    group_ids = _eligible_group_ids(group_rows)
    if not group_ids:
        raise RuntimeError("native dynamic gate contains no scientific-passing H1 root")
    root_rows = _matching_root_rows(group_ids, root_table)
    _validate_frozen_intake(output, group_rows, group_ids, root_table)
    _validate_native_genotype_freeze(output, repository, group_ids)
    _assert_independent_genotype_snapshot_current(snapshot, output, repository)
    _invalidate_independent_genotype_hash_ledger(output)
    continuation_rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []
    solver_rows: list[dict[str, Any]] = []
    for row in root_rows:
        try:
            continuation, profile, solvers = _run_one_root(row)
        except Exception as exc:
            continuation, profile, solvers = _failed_root_rows(row, exc)
        continuation_rows.extend(continuation)
        profile_rows.extend(profile)
        solver_rows.extend(solvers)
    completed = {
        transporter: len(
            {
                row["root_id"]
                for row in continuation_rows
                if row["transporter"] == transporter
                and row["expression"] == 0.0
                and row["continuation_completed_to_exact_zero"]
            }
        )
        for transporter in TRANSPORTERS
    }
    branch_losses = [
        {
            "root_id": row["root_id"],
            "transporter": row["transporter"],
            "failure": row["continuation_failure"],
        }
        for row in solver_rows
        if not row["continuation_completed_to_exact_zero"]
    ]
    production_comparison = _production_discrepancies(
        output, continuation_rows, profile_rows, solver_rows
    )
    runner_error_count = sum(
        bool(row.get("runner_exception")) for row in solver_rows
    )
    executed_solver_rows = [
        row for row in solver_rows if row.get("dynamic_executed")
    ]
    executed_profile_rows = [
        row for row in profile_rows if row.get("dynamic_executed")
    ]
    numerical_gate = bool(
        runner_error_count == 0
        and all(bool(row.get("cross_solver_gate_pass")) for row in executed_solver_rows)
        and all(
            bool(row.get("paired_numerical_gate_pass"))
            for row in executed_profile_rows
        )
    )
    reproduction_gate = bool(
        len(root_rows) == len(group_ids)
        and numerical_gate
        and production_comparison["comparison_gate_pass"]
    )
    summary = {
        "status": (
            "INDEPENDENT_TARGET_FREE_NATIVE_GENOTYPE_PASS"
            if reproduction_gate
            else "INDEPENDENT_TARGET_FREE_NATIVE_GENOTYPE_FAIL"
        ),
        "firewall": "NO_EXPERIMENTAL_GENOTYPE_TARGET_READ_OR_SCORED",
        "passing_dynamic_h1_root_count": len(group_ids),
        "passing_dynamic_h1_root_ids": list(group_ids),
        "all_passing_rows_retained": len(root_rows) == len(group_ids),
        "first_passer_selection_prohibited": True,
        "expression_step_count": len(EXPRESSION_GRID),
        "completed_to_exact_zero_count": completed,
        "branch_losses": branch_losses,
        "dynamic_profile_row_count": len(profile_rows),
        "solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        "independent_reproduction_relative_tolerance": (
            INDEPENDENT_REPRODUCTION_RELATIVE_TOLERANCE
        ),
        "runner_error_count": runner_error_count,
        "numerical_reproducibility_gate_pass": numerical_gate,
        "independent_reproduction_gate_pass": reproduction_gate,
        "maximum_deleted_transport_flux_fmol_s": max(
            (
                float(row["max_abs_deleted_transport_flux_fmol_s"])
                for row in profile_rows
                if row.get("dynamic_executed")
            ),
            default=None,
        ),
        "production_comparison": production_comparison,
    }
    with tempfile.TemporaryDirectory(
        prefix=".independent_native_genotype.", dir=output
    ) as temporary:
        staging = Path(temporary)
        _write_csv(
            staging / CONTINUATION_OUTPUT,
            continuation_rows,
            ("root_id", "transporter", "expression_step_index"),
        )
        _write_csv(
            staging / PROFILE_OUTPUT,
            profile_rows,
            ("root_id", "transporter", "solver_label"),
        )
        _write_csv(
            staging / SOLVER_OUTPUT,
            solver_rows,
            ("root_id", "transporter", "dynamic_executed"),
        )
        _write_json(staging / SUMMARY_OUTPUT, summary)
        manifest = {
            "manifest_id": "TASK13B_INDEPENDENT_NATIVE_GENOTYPE_PRE_REVEAL_V1",
            "status": (
                "FROZEN_TARGET_FREE_INDEPENDENT_REPRODUCTION_PASS"
                if reproduction_gate
                else "FROZEN_TARGET_FREE_INDEPENDENT_REPRODUCTION_FAIL"
            ),
            "firewall": "NO_EXPERIMENTAL_GENOTYPE_TARGET_READ_OR_SCORED",
            "root_ids": list(group_ids),
            "expression_grid": list(EXPRESSION_GRID),
            "all_passing_rows_retained": len(root_rows) == len(group_ids),
            "fitted_parameter_count": 0,
            "all_structural_parameters_and_other_fixed": True,
            "solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
            "independent_reproduction_relative_tolerance": (
                INDEPENDENT_REPRODUCTION_RELATIVE_TOLERANCE
            ),
            "independent_reproduction_gate_pass": reproduction_gate,
            "input_artifact_sha256": dict(snapshot.input_artifact_sha256),
            "equation_source_file_sha256": dict(
                snapshot.equation_source_file_sha256
            ),
            "equation_source_tree_sha256": _sha256_object(
                dict(snapshot.equation_source_file_sha256)
            ),
            "report_artifact_sha256": {
                filename: _sha256_file(staging / filename)
                for filename in INDEPENDENT_GENOTYPE_REPORTS
            },
            "production_comparison": production_comparison,
            "clean_auditor_authorization_required": True,
            "phenotype_reveal_permitted_by_this_manifest": False,
        }
        _write_json(staging / MANIFEST_OUTPUT, manifest)
        hash_rows = [
            {"artifact": filename, "sha256": _sha256_file(staging / filename)}
            for filename in INDEPENDENT_GENOTYPE_FREEZE_ARTIFACTS
        ]
        _write_csv(
            staging / HASH_OUTPUT,
            hash_rows,
            ("artifact", "sha256"),
        )
        _publish_staged_independent_genotype_freeze(
            staging, output, repository, snapshot
        )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    args = parser.parse_args(argv)
    output = args.repository / "results" / "13B_modern_full_model"
    summary = run_independent_native_genotypes(output)
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0 if summary["independent_reproduction_gate_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = (
    "CONTINUATION_OUTPUT",
    "HASH_OUTPUT",
    "MANIFEST_OUTPUT",
    "PROFILE_OUTPUT",
    "SOLVER_OUTPUT",
    "SUMMARY_OUTPUT",
    "run_independent_native_genotypes",
)
