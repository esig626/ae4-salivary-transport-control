"""Run target-free genotype continuations after the native WT gate freezes.

This runner is a strict downstream consumer of the native-source WT dynamic
contract.  It refuses to execute until the group gate, gate JSON, frozen
manifest, and native H1 root table all exist and agree.  Every group row with
``native_dynamic_scientific_gate_pass=true`` is retained; selecting the first
passing row is prohibited.

For every admitted H1 root the runner performs the complete 21-step WT-to-zero
expression continuation independently for AE4 and AE2.  Branch loss is a
reported result, not an exception and not a reason to skip later roots.  A
600-s CCh+IPR pair is run with Radau and BDF only when the corresponding
continuation reaches exact zero on a connected, full-rank, positive branch.

Outputs contain model-to-model ratios, state trajectories, solver checks, and
conservation diagnostics only.  No experimental genotype comparison is read,
scored, bounded, or used as a selection criterion.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict, dataclass, field, is_dataclass
from enum import Enum
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .genotype_evaluation import (
    DEFAULT_EXPRESSION_GRID,
    DualSolverGenotypeEvaluation,
    GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
    GenotypeContinuation,
    build_frozen_genotype_dynamic_model,
    continue_ae2_expression,
    continue_ae4_expression,
    evaluate_dual_solver_genotype_pair,
    frozen_dynamic_contract_metadata,
)
from .native_source_panel import (
    build_native_source_model,
    load_native_source_roots,
)
from .validation import (
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    assert_pre_reveal_payload,
    sha256_file,
    sha256_object,
)


DYNAMIC_GROUP_INPUT = "native_dynamic_contract_group_gate.csv"
DYNAMIC_GATE_INPUT = "native_dynamic_contract_gate.json"
DYNAMIC_MANIFEST_INPUT = "native_dynamic_contract_manifest.json"
DYNAMIC_HASH_INPUT = "native_dynamic_contract_hashes.csv"
DYNAMIC_PROFILE_INPUT = "native_dynamic_contract_profile.csv"
DYNAMIC_SOLVER_INPUT = "native_dynamic_contract_solver_crosschecks.csv"
DYNAMIC_ISOMORPHISM_INPUT = "native_dynamic_contract_isomorphism.csv"
DYNAMIC_NEARBY_INPUT = "native_dynamic_contract_nearby.csv"
NATIVE_ROOT_INPUT = "native_source_wt_roots.csv"
NATIVE_SUMMARY_INPUT = "native_source_wt_summary.json"

DYNAMIC_REPORT_INPUTS = (
    DYNAMIC_PROFILE_INPUT,
    DYNAMIC_GROUP_INPUT,
    DYNAMIC_SOLVER_INPUT,
    DYNAMIC_ISOMORPHISM_INPUT,
    DYNAMIC_NEARBY_INPUT,
)
DYNAMIC_HASH_LEDGER_ARTIFACTS = (
    *DYNAMIC_REPORT_INPUTS,
    DYNAMIC_MANIFEST_INPUT,
    DYNAMIC_GATE_INPUT,
)
NATIVE_INPUT_ARTIFACTS = (
    "native_source_map.json",
    "native_source_root_attempts.csv",
    NATIVE_ROOT_INPUT,
    NATIVE_SUMMARY_INPUT,
    "native_source_panel.csv",
    "native_source_scale_one_identity.csv",
    "native_source_hashes.csv",
    "independent_native_roots.csv",
    "independent_native_attempts.csv",
    "independent_native_summary.json",
)
DYNAMIC_EQUATION_SOURCES = (
    "src/modern_full_model/__init__.py",
    "src/modern_full_model/acid_base.py",
    "src/modern_full_model/camp_pka.py",
    "src/modern_full_model/calibration.py",
    "src/modern_full_model/independent.py",
    "src/modern_full_model/independent_native.py",
    "src/modern_full_model/membranes.py",
    "src/modern_full_model/model.py",
    "src/modern_full_model/native_dynamic_contract.py",
    "src/modern_full_model/native_source_panel.py",
    "src/modern_full_model/nkcc_stimulation.py",
    "src/modern_full_model/parameters.py",
    "src/modern_full_model/run_native_dynamic_contract.py",
    "src/modern_full_model/run_native_source_panel.py",
    "src/modern_full_model/run_independent_native.py",
    "src/modern_full_model/states.py",
    "src/modern_full_model/transporters.py",
    "src/modern_full_model/validation.py",
    "src/modern_full_model/water.py",
)
FIREWALL_HASH_KEYS = frozenset(
    {
        "calibration_ledger_sha256",
        "validation_ledger_sha256",
        "sealed_phenotype_target_ledger_sha256",
        "pre_reveal_log_sha256",
    }
)
GENOTYPE_SOURCE_DEPENDENCIES = (
    *DYNAMIC_EQUATION_SOURCES,
    "src/modern_full_model/genotype_evaluation.py",
    "src/modern_full_model/run_native_genotype_evaluation.py",
)
FROZEN_INPUT_ARTIFACTS = (
    *DYNAMIC_HASH_LEDGER_ARTIFACTS,
    DYNAMIC_HASH_INPUT,
    *NATIVE_INPUT_ARTIFACTS,
)

EXPECTED_DYNAMIC_GATE_ID = "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1"
EXPECTED_DYNAMIC_GATE_STATUS = "SCIENTIFIC_GATE_PASS_AUDITOR_AUTHORIZATION_REQUIRED"
EXPECTED_DYNAMIC_MANIFEST_ID = "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1"
EXPECTED_DYNAMIC_MANIFEST_STATUS = "FROZEN_FOR_CLEAN_ADVERSARIAL_REVIEW"

CONTINUATION_FILENAME = "native_genotype_continuation.csv"
TRAJECTORY_FILENAME = "native_genotype_trajectories.csv"
PROFILE_FILENAME = "native_genotype_profile.csv"
SOLVER_FILENAME = "native_genotype_solver_crosschecks.csv"
SUMMARY_FILENAME = "native_genotype_summary.json"
MANIFEST_FILENAME = "native_genotype_manifest.json"
HASH_FILENAME = "native_genotype_hashes.csv"
GENOTYPE_REPORT_FILENAMES = (
    CONTINUATION_FILENAME,
    TRAJECTORY_FILENAME,
    PROFILE_FILENAME,
    SOLVER_FILENAME,
    SUMMARY_FILENAME,
)
GENOTYPE_FREEZE_ARTIFACT_FILENAMES = (
    *GENOTYPE_REPORT_FILENAMES,
    MANIFEST_FILENAME,
)

TRANSPORTER_ORDER = ("AE4", "AE2")
SOLVER_ORDER = (PRODUCTION_RADAU.label, PRODUCTION_BDF.label)
GENOTYPE_ROLE_ORDER = ("FROZEN_WT", "EXACT_ZERO_EXPRESSION")


@dataclass(frozen=True)
class FrozenGenotypeIntake:
    """Cross-checked passing group rows and their exact frozen H1 roots."""

    group_rows: tuple[Mapping[str, str], ...]
    roots: tuple[Any, ...]
    dynamic_gate: Mapping[str, Any]
    dynamic_manifest: Mapping[str, Any]
    input_artifact_hashes: Mapping[str, str]


@dataclass(frozen=True)
class GenotypeDependencySnapshot:
    """Exact native-dynamic inputs and genotype sources used by one run."""

    input_artifact_sha256: Mapping[str, str]
    equation_source_file_sha256: Mapping[str, str]


@dataclass(frozen=True)
class NativeGenotypeRootResult:
    """Both transporter continuations and any licensed dynamic pairs."""

    root: Any
    continuations: Mapping[str, GenotypeContinuation]
    dynamics: Mapping[str, DualSolverGenotypeEvaluation]
    branch_errors: Mapping[str, str] = field(default_factory=dict)


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "__dict__"):
        return {
            str(key): _jsonable(item)
            for key, item in vars(value).items()
            if not str(key).startswith("_")
        }
    return repr(value)


def _json_text(value: Any) -> str:
    return json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    )


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
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    except ValueError as exc:
        raise ValueError(f"{path.name} is not a unique-key JSON object: {exc}") from exc
    if not isinstance(value, Mapping):
        raise ValueError(f"{path.name} must contain a JSON object")
    return value


def _csv_bool(value: Any, *, field: str) -> bool:
    """Accept only the canonical representation emitted by csv.DictWriter."""

    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(
        f"{field} must be bool or canonical case-sensitive True/False CSV text; "
        f"found {value!r}"
    )


def _require_json_bool(
    payload: Mapping[str, Any], field: str, expected: bool, *, artifact: str
) -> None:
    """Require a real JSON boolean; string lookalikes fail closed."""

    if field not in payload:
        raise ValueError(f"{artifact} lacks required boolean {field}")
    value = payload[field]
    if type(value) is not bool or value is not expected:
        raise ValueError(
            f"{artifact}.{field} must be JSON boolean {expected}, found {value!r}"
        )


def _require_json_integer(
    payload: Mapping[str, Any], field: str, *, artifact: str
) -> int:
    if field not in payload:
        raise ValueError(f"{artifact} lacks required integer {field}")
    value = payload[field]
    if type(value) is not int or value < 0:
        raise ValueError(
            f"{artifact}.{field} must be a nonnegative JSON integer, found {value!r}"
        )
    return value


def _unique_root_id_list(
    payload: Mapping[str, Any], field: str, *, artifact: str
) -> tuple[str, ...]:
    if field not in payload:
        raise ValueError(f"{artifact} lacks required root-ID list {field}")
    value = payload[field]
    if not isinstance(value, list):
        raise ValueError(f"{artifact}.{field} must be a JSON list")
    identifiers: list[str] = []
    for item in value:
        if type(item) is not str or not item:
            raise ValueError(f"{artifact}.{field} contains a non-string or empty root ID")
        identifiers.append(item)
    if len(set(identifiers)) != len(identifiers):
        raise ValueError(f"{artifact}.{field} contains duplicate root IDs")
    return tuple(sorted(identifiers))


def _require_hash(
    registry: Mapping[str, Any], filename: str, actual: str, *, registry_name: str
) -> None:
    if filename not in registry:
        raise ValueError(f"{registry_name} lacks SHA-256 binding for {filename}")
    expected = registry[filename]
    if (
        type(expected) is not str
        or len(expected) != 64
        or any(character not in "0123456789abcdef" for character in expected)
    ):
        raise ValueError(f"{registry_name}.{filename} is not a lowercase SHA-256 digest")
    if expected != actual:
        raise ValueError(f"stale manifest hash for {filename}")


def _require_exact_current_hash_registry(
    registry: Any,
    expected_filenames: Sequence[str],
    *,
    base: Path,
    registry_name: str,
) -> Mapping[str, str]:
    if not isinstance(registry, Mapping):
        raise ValueError(f"native dynamic manifest lacks {registry_name}")
    expected = set(expected_filenames)
    actual_keys = set(registry)
    if actual_keys != expected:
        missing = sorted(expected - actual_keys)
        extra = sorted(actual_keys - expected)
        raise ValueError(
            f"{registry_name} must bind exactly its declared artifacts; "
            f"missing={missing}, extra={extra}"
        )
    missing_files = sorted(name for name in expected if not (base / name).is_file())
    if missing_files:
        raise FileNotFoundError(
            f"{registry_name} references missing artifacts: {missing_files}"
        )
    current: dict[str, str] = {}
    for filename in expected_filenames:
        digest = sha256_file(base / filename)
        _require_hash(
            registry,
            filename,
            digest,
            registry_name=registry_name,
        )
        current[filename] = digest
    return current


def _require_exact_firewall_registry(registry: Any) -> None:
    """Validate the opaque registry only; target ledgers are never parsed here."""

    if not isinstance(registry, Mapping) or set(registry) != FIREWALL_HASH_KEYS:
        raise ValueError("opaque firewall registry must contain exactly four bindings")
    for label, digest in registry.items():
        if (
            type(digest) is not str
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ValueError(f"opaque firewall binding {label} is not lowercase SHA-256")


def _load_and_verify_dynamic_hash_ledger(
    results_directory: Path,
) -> Mapping[str, str]:
    path = results_directory / DYNAMIC_HASH_INPUT
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != ("artifact", "sha256"):
            raise ValueError(
                "native dynamic hash ledger must have exact artifact,sha256 columns"
            )
        rows = list(reader)
    if not rows:
        raise ValueError("native dynamic hash ledger is empty")
    registry: dict[str, str] = {}
    for row in rows:
        artifact = row["artifact"]
        if artifact in registry:
            raise ValueError(f"native dynamic hash ledger duplicates {artifact}")
        registry[artifact] = row["sha256"]
    expected = set(DYNAMIC_HASH_LEDGER_ARTIFACTS)
    if set(registry) != expected:
        raise ValueError(
            "native dynamic hash ledger must contain exactly seven frozen artifacts"
        )
    for artifact in DYNAMIC_HASH_LEDGER_ARTIFACTS:
        target = results_directory / artifact
        if not target.is_file():
            raise FileNotFoundError(f"native dynamic hash ledger target is missing: {artifact}")
        _require_hash(
            registry,
            artifact,
            sha256_file(target),
            registry_name=DYNAMIC_HASH_INPUT,
        )
    return registry


def _root_value(root: Any, name: str, default: Any = None) -> Any:
    if hasattr(root, name):
        return getattr(root, name)
    if isinstance(root, Mapping):
        return root.get(name, default)
    return default


def _root_id(root: Any) -> str:
    value = _root_value(root, "root_id")
    if value is None or not str(value):
        raise ValueError("native root lacks a nonempty root_id")
    return str(value)


def _root_metadata(root: Any) -> dict[str, Any]:
    names = (
        "root_id",
        "panel_id",
        "topology_id",
        "routing_id",
        "source_class",
        "source_scale",
        "pump_capacity_scale",
        "apical_pump_fraction",
        "apical_k_fraction",
        "ae4_cation_fraction",
        "hydraulic_scale",
        "hydraulic_mode",
        "eligibility",
    )
    return {name: _root_value(root, name) for name in names}


def _validate_h1_production_root(root: Any) -> None:
    root_id = _root_id(root)
    numerical = _root_value(root, "passes_numerical_gate")
    wt = _root_value(root, "passes_wt_gate")
    if type(numerical) is not bool:
        raise ValueError(f"{root_id}: numerical gate must be an in-memory bool")
    if type(wt) is not bool:
        raise ValueError(f"{root_id}: WT gate must be an in-memory bool")
    if not numerical:
        raise ValueError(f"{root_id}: frozen root does not pass its numerical gate")
    if not wt:
        raise ValueError(f"{root_id}: frozen root does not pass its WT gate")
    if _root_value(root, "eligibility") != "PRODUCTION_PRE_REVEAL":
        raise ValueError(
            f"{root_id}: eligibility must be exact token PRODUCTION_PRE_REVEAL"
        )
    if float(_root_value(root, "hydraulic_scale", math.nan)) != 1.0:
        raise ValueError(f"{root_id}: genotype runner admits H=1 roots only")
    if str(_root_value(root, "hydraulic_mode", "")) != "H1":
        raise ValueError(f"{root_id}: genotype runner admits hydraulic mode H1 only")
    core = tuple(float(value) for value in _root_value(root, "core_state", ()))
    if len(core) != 12 or not all(math.isfinite(value) and value > 0.0 for value in core):
        raise ValueError(f"{root_id}: frozen root must contain a positive 12-state core")


def passing_dynamic_group_rows(rows: Iterable[Mapping[str, str]]) -> tuple[Mapping[str, str], ...]:
    """Return every passing scientific group row in deterministic root order."""

    passing: list[Mapping[str, str]] = []
    all_ids: list[str] = []
    for row in rows:
        if "native_dynamic_scientific_gate_pass" not in row:
            raise ValueError("dynamic group table lacks its scientific gate column")
        root_id = str(row.get("root_id", ""))
        if not root_id:
            raise ValueError("dynamic group row lacks root_id")
        all_ids.append(root_id)
        if _csv_bool(
            row["native_dynamic_scientific_gate_pass"],
            field="native_dynamic_scientific_gate_pass",
        ):
            if "hydraulic_scale" not in row or not str(row["hydraulic_scale"]):
                raise ValueError(f"{root_id}: passing genotype intake row lacks H scale")
            if float(row["hydraulic_scale"]) != 1.0:
                raise ValueError(f"{root_id}: passing genotype intake row is not H=1")
            if "hydraulic_mode" not in row or row["hydraulic_mode"] != "H1":
                raise ValueError(f"{root_id}: passing genotype intake row is not mode H1")
            if "all_native_rest_wt_gate_pass" not in row:
                raise ValueError(f"{root_id}: dynamic row lacks its native WT rest gate")
            if not _csv_bool(
                row["all_native_rest_wt_gate_pass"],
                field="all_native_rest_wt_gate_pass",
            ):
                raise ValueError(f"{root_id}: dynamic pass lacks the native WT rest gate")
            passing.append(dict(row))
    if len(set(all_ids)) != len(all_ids):
        raise ValueError("dynamic group table contains duplicate root IDs")
    passing.sort(key=lambda row: str(row["root_id"]))
    return tuple(passing)


def match_passing_native_roots(
    group_rows: Sequence[Mapping[str, str]], roots: Iterable[Any]
) -> tuple[Any, ...]:
    """Match all and only passing group IDs to frozen H1 root objects."""

    root_registry: dict[str, Any] = {}
    for root in roots:
        root_id = _root_id(root)
        if root_id in root_registry:
            raise ValueError(f"native root table contains duplicate root {root_id}")
        root_registry[root_id] = root
    passing_ids = tuple(str(row["root_id"]) for row in group_rows)
    missing = tuple(root_id for root_id in passing_ids if root_id not in root_registry)
    if missing:
        raise ValueError(f"dynamic passing roots missing from frozen root table: {missing}")
    matched = tuple(root_registry[root_id] for root_id in sorted(passing_ids))
    for root in matched:
        _validate_h1_production_root(root)
    if len(matched) != len(group_rows):
        raise AssertionError("all passing dynamic rows must be retained exactly once")
    return matched


def load_frozen_genotype_intake(
    results_directory: str | Path,
    *,
    repository_root: str | Path,
) -> FrozenGenotypeIntake:
    """Load and cross-check the complete frozen WT prerequisite set."""

    output = Path(results_directory)
    repository = Path(repository_root)
    required = FROZEN_INPUT_ARTIFACTS
    missing = tuple(name for name in required if not (output / name).is_file())
    if missing:
        raise RuntimeError(
            "native genotype evaluation is blocked until the frozen native WT "
            f"dynamic gate exists; missing artifacts: {missing}"
        )
    all_group_rows = _read_csv(output / DYNAMIC_GROUP_INPUT)
    group_rows = passing_dynamic_group_rows(all_group_rows)
    if not group_rows:
        raise RuntimeError("native WT dynamic gate contains no scientific-passing H1 row")
    all_group_ids = tuple(sorted(str(row["root_id"]) for row in all_group_rows))
    dynamic_gate = _read_json(output / DYNAMIC_GATE_INPUT)
    dynamic_manifest = _read_json(output / DYNAMIC_MANIFEST_INPUT)
    if dynamic_gate.get("gate_id") != EXPECTED_DYNAMIC_GATE_ID:
        raise ValueError("native dynamic gate has an unexpected gate_id")
    if dynamic_gate.get("status") != EXPECTED_DYNAMIC_GATE_STATUS:
        raise RuntimeError(
            "native dynamic composite scientific gate is not explicitly passing"
        )
    if dynamic_gate.get("manifest") != DYNAMIC_MANIFEST_INPUT:
        raise ValueError("native dynamic gate does not bind the expected manifest")
    _require_json_bool(
        dynamic_gate,
        "at_least_one_scientific_group_pass",
        True,
        artifact=DYNAMIC_GATE_INPUT,
    )
    _require_json_bool(
        dynamic_gate,
        "all_rows_retained_no_first_passer",
        True,
        artifact=DYNAMIC_GATE_INPUT,
    )
    _require_json_bool(
        dynamic_gate,
        "clean_auditor_authorization_required",
        True,
        artifact=DYNAMIC_GATE_INPUT,
    )
    for field in (
        "clean_auditor_authorization_recorded",
        "independent_numerical_reproduction_recorded_by_this_runner",
        "final_composite_authorization_complete",
        "phenotype_reveal_permitted",
    ):
        _require_json_bool(
            dynamic_gate, field, False, artifact=DYNAMIC_GATE_INPUT
        )
    _require_json_bool(
        dynamic_gate,
        "independent_numerical_reproduction_required",
        True,
        artifact=DYNAMIC_GATE_INPUT,
    )

    group_ids = tuple(sorted(str(row["root_id"]) for row in group_rows))
    gate_ids = _unique_root_id_list(
        dynamic_gate,
        "scientific_passing_root_ids",
        artifact=DYNAMIC_GATE_INPUT,
    )
    if group_ids != gate_ids:
        raise ValueError("dynamic gate JSON and group CSV disagree on passing root IDs")
    gate_count = _require_json_integer(
        dynamic_gate,
        "scientific_passing_group_count",
        artifact=DYNAMIC_GATE_INPUT,
    )
    if gate_count != len(group_rows):
        raise ValueError("dynamic gate JSON passing count disagrees with group CSV")
    gate_native_root_count = _require_json_integer(
        dynamic_gate,
        "native_root_count",
        artifact=DYNAMIC_GATE_INPUT,
    )

    if dynamic_manifest.get("manifest_id") != EXPECTED_DYNAMIC_MANIFEST_ID:
        raise ValueError("native dynamic manifest has an unexpected manifest_id")
    if dynamic_manifest.get("status") != EXPECTED_DYNAMIC_MANIFEST_STATUS:
        raise ValueError("native dynamic manifest is not in its frozen review state")
    _require_json_bool(
        dynamic_manifest,
        "clean_auditor_authorization_required",
        True,
        artifact=DYNAMIC_MANIFEST_INPUT,
    )
    _require_json_bool(
        dynamic_manifest,
        "independent_numerical_reproduction_required_for_composite_authorization",
        True,
        artifact=DYNAMIC_MANIFEST_INPUT,
    )
    _require_json_bool(
        dynamic_manifest,
        "phenotype_reveal_permitted_by_this_manifest",
        False,
        artifact=DYNAMIC_MANIFEST_INPUT,
    )
    manifest_passing_ids = _unique_root_id_list(
        dynamic_manifest,
        "scientific_passing_root_ids",
        artifact=DYNAMIC_MANIFEST_INPUT,
    )
    if group_ids != manifest_passing_ids:
        raise ValueError(
            "dynamic manifest and group CSV disagree on passing root IDs"
        )
    manifest_passing_count = _require_json_integer(
        dynamic_manifest,
        "scientific_passing_group_count",
        artifact=DYNAMIC_MANIFEST_INPUT,
    )
    if manifest_passing_count != len(group_rows):
        raise ValueError("dynamic manifest passing count disagrees with group CSV")
    manifest_roots = dynamic_manifest.get("roots")
    if not isinstance(manifest_roots, Mapping):
        raise ValueError("native dynamic manifest lacks its frozen root registry")
    manifest_root_ids = tuple(sorted(str(root_id) for root_id in manifest_roots))
    if manifest_root_ids != all_group_ids:
        raise ValueError(
            "native dynamic manifest root IDs must exactly equal every group CSV root ID"
        )
    manifest_root_count = _require_json_integer(
        dynamic_manifest, "root_count", artifact=DYNAMIC_MANIFEST_INPUT
    )
    if (
        manifest_root_count != len(manifest_roots)
        or manifest_root_count != len(all_group_rows)
        or gate_native_root_count != manifest_root_count
    ):
        raise ValueError(
            "native root counts must agree across gate, manifest, and every group CSV row"
        )

    _load_and_verify_dynamic_hash_ledger(output)
    _require_exact_current_hash_registry(
        dynamic_manifest.get("wt_dynamic_report_artifact_sha256"),
        DYNAMIC_REPORT_INPUTS,
        base=output,
        registry_name="wt_dynamic_report_artifact_sha256",
    )
    _require_exact_current_hash_registry(
        dynamic_manifest.get("native_input_artifact_sha256"),
        NATIVE_INPUT_ARTIFACTS,
        base=output,
        registry_name="native_input_artifact_sha256",
    )
    equation_hashes = _require_exact_current_hash_registry(
        dynamic_manifest.get("equation_source_file_sha256"),
        DYNAMIC_EQUATION_SOURCES,
        base=repository,
        registry_name="equation_source_file_sha256",
    )
    equation_tree_hash = dynamic_manifest.get("equation_source_tree_sha256")
    if equation_tree_hash != sha256_object(dict(equation_hashes)):
        raise ValueError("native dynamic manifest equation source tree hash is stale")
    _require_exact_firewall_registry(
        dynamic_manifest.get("opaque_firewall_ledger_sha256")
    )
    if dynamic_manifest.get("opaque_firewall_access_mode") != (
        "SHA256_BYTES_ONLY_NO_TARGET_PARSE"
    ):
        raise ValueError("native dynamic manifest firewall access mode is not frozen")

    roots = match_passing_native_roots(
        group_rows, load_native_source_roots(output / NATIVE_ROOT_INPUT)
    )
    for root in roots:
        root_id = _root_id(root)
        payload = manifest_roots[root_id]
        if not isinstance(payload, Mapping):
            raise ValueError(f"{root_id}: frozen manifest root payload is malformed")
        expected_hash = str(payload.get("core_state_sha256", ""))
        actual_hash = sha256_object(tuple(_root_value(root, "core_state")))
        if not expected_hash or expected_hash != actual_hash:
            raise ValueError(f"{root_id}: frozen core-state hash mismatch")
    hashes = {name: sha256_file(output / name) for name in required}
    return FrozenGenotypeIntake(
        group_rows=group_rows,
        roots=roots,
        dynamic_gate=dynamic_gate,
        dynamic_manifest=dynamic_manifest,
        input_artifact_hashes=hashes,
    )


def _evaluate_root(root: Any) -> NativeGenotypeRootResult:
    """Execute both complete continuations and licensed dynamics for one root."""

    continuations: dict[str, GenotypeContinuation] = {}
    dynamics: dict[str, DualSolverGenotypeEvaluation] = {}
    errors: dict[str, str] = {}
    try:
        resting = build_native_source_model(root)
        model = build_frozen_genotype_dynamic_model(resting)
    except Exception as exc:
        message = f"MODEL_BUILD_EXCEPTION:{type(exc).__name__}: {exc}"
        return NativeGenotypeRootResult(
            root=root,
            continuations={},
            dynamics={},
            branch_errors={transporter: message for transporter in TRANSPORTER_ORDER},
        )
    for transporter in TRANSPORTER_ORDER:
        try:
            continuation = (
                continue_ae4_expression(model, _root_value(root, "core_state"))
                if transporter == "AE4"
                else continue_ae2_expression(model, _root_value(root, "core_state"))
            )
            continuations[transporter] = continuation
        except Exception as exc:
            errors[transporter] = (
                f"CONTINUATION_EXCEPTION:{type(exc).__name__}: {exc}"
            )
            continue
        null_root = continuation.connected_null_root
        if null_root is None:
            continue
        try:
            dynamics[transporter] = evaluate_dual_solver_genotype_pair(
                model,
                _root_value(root, "core_state"),
                null_root.core_state,
                transporter=transporter,
            )
        except Exception as exc:
            errors[transporter] = f"DYNAMIC_EXCEPTION:{type(exc).__name__}: {exc}"
    return NativeGenotypeRootResult(
        root=root,
        continuations=continuations,
        dynamics=dynamics,
        branch_errors=errors,
    )


def _run_roots(roots: Sequence[Any], workers: int) -> tuple[NativeGenotypeRootResult, ...]:
    if workers < 1:
        raise ValueError("workers must be at least one")
    if workers == 1:
        results = tuple(_evaluate_root(root) for root in roots)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = tuple(executor.map(_evaluate_root, roots, chunksize=1))
    return tuple(sorted(results, key=lambda result: _root_id(result.root)))


def _continuation_rows(results: Sequence[NativeGenotypeRootResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        metadata = _root_metadata(result.root)
        for transporter in TRANSPORTER_ORDER:
            report = result.continuations.get(transporter)
            if report is None:
                rows.append(
                    {
                        **metadata,
                        "transporter": transporter,
                        "expression_step_index": None,
                        "expression": None,
                        "declared_expression_step_count": len(DEFAULT_EXPRESSION_GRID),
                        "independent_equation_count": 10,
                        "fitted_coordinate_count": 10,
                        "fitted_parameter_count": 0,
                        "fixed_other_impermeant_osmoles_fmol": _root_value(
                            result.root, "other_impermeant_osmoles_fmol"
                        ),
                        "attempt_count": 0,
                        "attempts_json": _json_text(()),
                        "retained_root_count": 0,
                        "retained_roots_json": _json_text(()),
                        "selected_root_id": None,
                        "selected_distance_from_previous_normalized": None,
                        "connected_to_previous": False,
                        "alternate_root_count": 0,
                        "step_status": "RUNNER_EXCEPTION",
                        "selected_coordinates_json": None,
                        "selected_core_state_json": None,
                        "selected_max_abs_scaled_independent_rhs": None,
                        "selected_max_abs_omitted_charge_rhs": None,
                        "selected_max_abs_charge_residual_fmol": None,
                        "selected_max_abs_current_residual_A": None,
                        "selected_jacobian_rank": None,
                        "selected_boundary_hits_json": None,
                        "selected_exact_deleted_transport_zero": None,
                        "selected_numerical_gate_pass": None,
                        "selected_gate_failures_json": None,
                        "continuation_completed_to_exact_zero": False,
                        "continuation_failure": result.branch_errors.get(
                            transporter, "continuation result missing"
                        ),
                        "runner_exception": True,
                    }
                )
                continue
            for index, step in enumerate(report.steps):
                selected = step.selected_root
                rows.append(
                    {
                        **metadata,
                        "transporter": transporter,
                        "expression_step_index": index,
                        "expression": step.expression,
                        "declared_expression_step_count": len(report.expression_grid),
                        "independent_equation_count": report.independent_equation_count,
                        "fitted_coordinate_count": report.fitted_coordinate_count,
                        "fitted_parameter_count": report.fitted_parameter_count,
                        "fixed_other_impermeant_osmoles_fmol": (
                            report.fixed_other_impermeant_osmoles_fmol
                        ),
                        "attempt_count": len(step.attempts),
                        "attempts_json": _json_text(step.attempts),
                        "retained_root_count": len(step.roots),
                        "retained_roots_json": _json_text(step.roots),
                        "selected_root_id": step.selected_root_id,
                        "selected_distance_from_previous_normalized": (
                            step.selected_distance_from_previous_normalized
                        ),
                        "connected_to_previous": step.connected_to_previous,
                        "alternate_root_count": step.alternate_root_count,
                        "step_status": step.status,
                        "selected_coordinates_json": (
                            _json_text(selected.coordinates) if selected else None
                        ),
                        "selected_core_state_json": (
                            _json_text(selected.core_state) if selected else None
                        ),
                        "selected_max_abs_scaled_independent_rhs": (
                            selected.max_abs_scaled_independent_rhs if selected else None
                        ),
                        "selected_max_abs_omitted_charge_rhs": (
                            selected.max_abs_omitted_charge_rhs_fmol_equivalent_s
                            if selected
                            else None
                        ),
                        "selected_max_abs_charge_residual_fmol": (
                            selected.max_abs_charge_residual_fmol if selected else None
                        ),
                        "selected_max_abs_current_residual_A": (
                            selected.max_abs_current_residual_A if selected else None
                        ),
                        "selected_jacobian_rank": (
                            selected.normalized_jacobian_rank if selected else None
                        ),
                        "selected_boundary_hits_json": (
                            _json_text(selected.boundary_hits) if selected else None
                        ),
                        "selected_exact_deleted_transport_zero": (
                            selected.exact_deleted_transport_zero if selected else None
                        ),
                        "selected_numerical_gate_pass": (
                            selected.passes_numerical_gate if selected else None
                        ),
                        "selected_gate_failures_json": (
                            _json_text(selected.gate_failures) if selected else None
                        ),
                        "continuation_completed_to_exact_zero": (
                            report.completed_to_exact_zero
                        ),
                        "continuation_failure": report.branch_failure,
                        "runner_exception": False,
                    }
                )
    return rows


def _trajectory_rows(results: Sequence[NativeGenotypeRootResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        metadata = _root_metadata(result.root)
        for transporter in TRANSPORTER_ORDER:
            dynamic = result.dynamics.get(transporter)
            if dynamic is None:
                continue
            for solver_label in SOLVER_ORDER:
                role_trajectories = (
                    (GENOTYPE_ROLE_ORDER[0], dynamic.wt_trajectories[solver_label]),
                    (
                        GENOTYPE_ROLE_ORDER[1],
                        dynamic.genotype_trajectories[solver_label],
                    ),
                )
                for role, trajectory in role_trajectories:
                    for index, time_s in enumerate(trajectory.time_s):
                        state = trajectory.states[:, index]
                        rows.append(
                            {
                                **metadata,
                                "transporter": transporter,
                                "genotype_role": role,
                                "genotype_name": trajectory.genotype_name,
                                "solver_label": solver_label,
                                "solver_method": trajectory.solver.method,
                                "time_s": float(time_s),
                                "state_json": _json_text(
                                    dict(zip(trajectory.state_names, state))
                                ),
                                "flow_pL_s": float(trajectory.flow_pL_s[index]),
                                "cumulative_flow_pL": float(
                                    trajectory.cumulative_flow_pL[index]
                                ),
                                "cell_na_mM": float(trajectory.cell_na_mM[index]),
                                "cell_k_mM": float(trajectory.cell_k_mM[index]),
                                "cell_cl_mM": float(trajectory.cell_cl_mM[index]),
                                "cell_ph": float(trajectory.cell_ph[index]),
                                "cell_volume_pL": float(
                                    trajectory.cell_volume_pL[index]
                                ),
                                "lumen_na_mM": float(trajectory.lumen_na_mM[index]),
                                "lumen_k_mM": float(trajectory.lumen_k_mM[index]),
                                "lumen_cl_mM": float(trajectory.lumen_cl_mM[index]),
                                "lumen_ph": float(trajectory.lumen_ph[index]),
                                "lumen_volume_pL": float(
                                    trajectory.lumen_volume_pL[index]
                                ),
                                "ae4_capacity_multiplier": float(
                                    trajectory.ae4_capacity_multiplier[index]
                                ),
                                "trajectory_success": trajectory.success,
                                "positive_core": trajectory.positive_core,
                                "all_flow_nonnegative": (
                                    trajectory.all_flow_nonnegative
                                ),
                                "numerical_gate_pass": trajectory.numerical_gate_pass,
                                "max_abs_deleted_transport_flux_fmol_s": (
                                    trajectory.max_abs_deleted_transport_flux_fmol_s
                                ),
                                "max_dimensionless_conservation_ratio": (
                                    trajectory.max_dimensionless_conservation_ratio
                                ),
                            }
                        )
    rows.sort(
        key=lambda row: (
            str(row["root_id"]),
            TRANSPORTER_ORDER.index(str(row["transporter"])),
            GENOTYPE_ROLE_ORDER.index(str(row["genotype_role"])),
            SOLVER_ORDER.index(str(row["solver_label"])),
            float(row["time_s"]),
        )
    )
    return rows


def _profile_rows(results: Sequence[NativeGenotypeRootResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        metadata = _root_metadata(result.root)
        for transporter in TRANSPORTER_ORDER:
            continuation = result.continuations.get(transporter)
            dynamic = result.dynamics.get(transporter)
            if dynamic is None:
                completed = bool(
                    continuation is not None
                    and continuation.completed_to_exact_zero
                )
                failure = result.branch_errors.get(
                    transporter,
                    continuation.branch_failure if continuation is not None else None,
                )
                rows.append(
                    {
                        **metadata,
                        "transporter": transporter,
                        "solver_label": None,
                        "dynamic_executed": False,
                        "continuation_completed_to_exact_zero": completed,
                        "continuation_failure": failure,
                        "paired_numerical_gate_pass": None,
                        "cumulative_ratio_at_landmarks_json": None,
                        "flow_ratio_at_landmarks_json": None,
                        "wt_endpoint_json": None,
                        "genotype_endpoint_json": None,
                        "genotype_to_wt_endpoint_json": None,
                        "wt_max_dimensionless_conservation_ratio": None,
                        "genotype_max_dimensionless_conservation_ratio": None,
                        "max_abs_deleted_transport_flux_fmol_s": None,
                        "runner_exception": transporter in result.branch_errors,
                    }
                )
                continue
            for solver_label in SOLVER_ORDER:
                summary = dynamic.evaluations[solver_label]
                genotype_trajectory = dynamic.genotype_trajectories[solver_label]
                rows.append(
                    {
                        **metadata,
                        "transporter": transporter,
                        "solver_label": solver_label,
                        "dynamic_executed": True,
                        "continuation_completed_to_exact_zero": True,
                        "continuation_failure": None,
                        "paired_numerical_gate_pass": (
                            summary.paired_numerical_gate_pass
                        ),
                        "cumulative_ratio_at_landmarks_json": _json_text(
                            summary.cumulative_ratio_at_landmarks
                        ),
                        "flow_ratio_at_landmarks_json": _json_text(
                            summary.flow_ratio_at_landmarks
                        ),
                        "wt_endpoint_json": _json_text(summary.wt_endpoint),
                        "genotype_endpoint_json": _json_text(
                            summary.genotype_endpoint
                        ),
                        "genotype_to_wt_endpoint_json": _json_text(
                            summary.genotype_to_wt_endpoint_ratios
                        ),
                        "wt_max_dimensionless_conservation_ratio": (
                            summary.wt_max_dimensionless_conservation_ratio
                        ),
                        "genotype_max_dimensionless_conservation_ratio": (
                            summary.genotype_max_dimensionless_conservation_ratio
                        ),
                        "max_abs_deleted_transport_flux_fmol_s": (
                            genotype_trajectory.max_abs_deleted_transport_flux_fmol_s
                        ),
                        "runner_exception": False,
                    }
                )
    return rows


def _solver_rows(results: Sequence[NativeGenotypeRootResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        metadata = _root_metadata(result.root)
        for transporter in TRANSPORTER_ORDER:
            continuation = result.continuations.get(transporter)
            dynamic = result.dynamics.get(transporter)
            row: dict[str, Any] = {
                **metadata,
                "transporter": transporter,
                "dynamic_executed": dynamic is not None,
                "continuation_completed_to_exact_zero": (
                    continuation.completed_to_exact_zero
                    if continuation is not None
                    else False
                ),
                "continuation_failure": result.branch_errors.get(
                    transporter,
                    continuation.branch_failure if continuation is not None else None,
                ),
                "both_solver_pairs_pass": (
                    dynamic.both_solver_pairs_pass if dynamic else None
                ),
                "cross_solver_relative_tolerance": (
                    dynamic.cross_solver_relative_tolerance
                    if dynamic
                    else GENOTYPE_SOLVER_RELATIVE_TOLERANCE
                ),
                "cross_solver_gate_pass": (
                    dynamic.cross_solver_gate_pass if dynamic else None
                ),
                "runner_exception": transporter in result.branch_errors,
            }
            for name in (
                "wt_state",
                "wt_flow",
                "wt_cumulative",
                "genotype_state",
                "genotype_flow",
                "genotype_cumulative",
            ):
                row[f"max_relative_difference_{name}"] = (
                    dynamic.cross_solver_relative_differences[name]
                    if dynamic
                    else None
                )
            rows.append(row)
    return rows


def _write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    empty_fields: Sequence[str],
) -> None:
    fields: list[str] = []
    for row in rows:
        for name in row:
            if name not in fields:
                fields.append(name)
    if not fields:
        fields = list(empty_fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    serializable = _jsonable(payload)
    assert_pre_reveal_payload(serializable)
    path.write_text(
        json.dumps(serializable, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _summary_payload(
    intake: FrozenGenotypeIntake,
    results: Sequence[NativeGenotypeRootResult],
) -> Mapping[str, Any]:
    continuation_status: dict[str, dict[str, Any]] = {}
    losses: list[dict[str, Any]] = []
    dynamic_pairs = 0
    dynamic_numerical_pairs = 0
    runner_error_count = 0
    for transporter in TRANSPORTER_ORDER:
        reports = [
            result.continuations[transporter]
            for result in results
            if transporter in result.continuations
        ]
        completed = [report for report in reports if report.completed_to_exact_zero]
        continuation_status[transporter] = {
            "attempted_root_count": len(results),
            "reported_continuation_count": len(reports),
            "completed_to_exact_zero_count": len(completed),
            "branch_loss_count": len(results) - len(completed),
        }
    for result in results:
        for transporter in TRANSPORTER_ORDER:
            continuation = result.continuations.get(transporter)
            dynamic = result.dynamics.get(transporter)
            runner_error = result.branch_errors.get(transporter)
            runner_error_count += int(runner_error is not None)
            if continuation is None or not continuation.completed_to_exact_zero:
                losses.append(
                    {
                        "root_id": _root_id(result.root),
                        "transporter": transporter,
                        "failure": (
                            runner_error
                            or (
                                continuation.branch_failure
                                if continuation is not None
                                else "continuation result missing"
                            )
                        ),
                        "last_expression": (
                            continuation.steps[-1].expression
                            if continuation is not None and continuation.steps
                            else None
                        ),
                        "last_step_status": (
                            continuation.steps[-1].status
                            if continuation is not None and continuation.steps
                            else "RUNNER_EXCEPTION" if runner_error else "NO_STEP"
                        ),
                    }
                )
            elif runner_error is not None:
                losses.append(
                    {
                        "root_id": _root_id(result.root),
                        "transporter": transporter,
                        "failure": runner_error,
                        "last_expression": continuation.steps[-1].expression,
                        "last_step_status": "DYNAMIC_EXCEPTION",
                    }
                )
            if dynamic is not None:
                dynamic_pairs += 1
                dynamic_numerical_pairs += int(dynamic.both_solver_pairs_pass)
    numerical_gate = bool(
        runner_error_count == 0
        and dynamic_numerical_pairs == dynamic_pairs
    )
    return {
        "summary_id": "TASK13B_NATIVE_GENOTYPE_TARGET_FREE_V1",
        "status": "TARGET_FREE_GENOTYPE_EVALUATION_COMPLETE",
        "passing_dynamic_group_count": len(intake.group_rows),
        "evaluated_root_count": len(results),
        "evaluated_root_ids": [_root_id(result.root) for result in results],
        "all_passing_rows_retained": len(results) == len(intake.group_rows),
        "first_passer_selection_prohibited": True,
        "expression_step_count": len(DEFAULT_EXPRESSION_GRID),
        "continuation_status": continuation_status,
        "branch_losses": losses,
        "dynamic_pair_count": dynamic_pairs,
        "dynamic_pairs_passing_both_numerical_solvers": dynamic_numerical_pairs,
        "solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        "runner_error_count": runner_error_count,
        "target_free_numerical_reproducibility_gate_pass": numerical_gate,
        "comparison_scope": "MODEL_TO_MODEL_DESCRIPTIVE_PROFILE_ONLY",
        "experimental_genotype_comparison_loaded": False,
        "experimental_acceptance_criterion_applied": False,
    }


def _source_hashes(repository_root: Path) -> Mapping[str, str]:
    missing = tuple(
        path
        for path in GENOTYPE_SOURCE_DEPENDENCIES
        if not (repository_root / path).is_file()
    )
    if missing:
        raise FileNotFoundError(
            f"genotype scientific source dependency set is incomplete: {missing}"
        )
    return {
        path: sha256_file(repository_root / path)
        for path in GENOTYPE_SOURCE_DEPENDENCIES
    }


def _capture_genotype_dependency_snapshot(
    repository_root: Path, results_directory: Path
) -> GenotypeDependencySnapshot:
    missing_inputs = tuple(
        name for name in FROZEN_INPUT_ARTIFACTS
        if not (results_directory / name).is_file()
    )
    if missing_inputs:
        raise FileNotFoundError(
            f"genotype dependency snapshot lacks frozen inputs: {missing_inputs}"
        )
    return GenotypeDependencySnapshot(
        input_artifact_sha256={
            name: sha256_file(results_directory / name)
            for name in FROZEN_INPUT_ARTIFACTS
        },
        equation_source_file_sha256=dict(_source_hashes(repository_root)),
    )


def _assert_genotype_snapshot_current(
    snapshot: GenotypeDependencySnapshot,
    repository_root: Path,
    results_directory: Path,
) -> None:
    current = _capture_genotype_dependency_snapshot(
        repository_root, results_directory
    )
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
            "native genotype dependencies changed during execution: "
            f"{changed}"
        )


def _manifest_payload(
    repository_root: Path,
    results_directory: Path,
    intake: FrozenGenotypeIntake,
    results: Sequence[NativeGenotypeRootResult],
    report_filenames: Sequence[str],
    *,
    report_directory: Path,
    dependency_snapshot: GenotypeDependencySnapshot,
) -> Mapping[str, Any]:
    source_hashes = dict(dependency_snapshot.equation_source_file_sha256)
    roots = {
        _root_id(result.root): {
            "metadata": _root_metadata(result.root),
            "core_state_sha256": sha256_object(
                tuple(_root_value(result.root, "core_state"))
            ),
            "other_impermeant_osmoles_fmol": float(
                _root_value(result.root, "other_impermeant_osmoles_fmol")
            ),
            "continuation_completed": {
                transporter: bool(
                    result.continuations.get(transporter)
                    and result.continuations[transporter].completed_to_exact_zero
                )
                for transporter in TRANSPORTER_ORDER
            },
            "branch_errors": dict(result.branch_errors),
        }
        for result in results
    }
    return {
        "manifest_id": "TASK13B_NATIVE_GENOTYPE_TARGET_FREE_V1",
        "status": "FROZEN_TARGET_FREE_GENOTYPE_PREDICTIONS",
        "parent_native_dynamic_manifest": DYNAMIC_MANIFEST_INPUT,
        "parent_native_dynamic_gate": DYNAMIC_GATE_INPUT,
        "root_intake_rule": {
            "gate_column": "native_dynamic_scientific_gate_pass",
            "retain_all_true_rows": True,
            "first_passer_selection_prohibited": True,
            "hydraulic_scale": 1.0,
            "hydraulic_mode": "H1",
        },
        "root_count": len(results),
        "roots": roots,
        "expression_grid": list(DEFAULT_EXPRESSION_GRID),
        "continuation_contract": {
            "independent_equation_count": 10,
            "fitted_coordinate_count": 10,
            "fitted_parameter_count": 0,
            "other_fixed": True,
            "branch_selection": "NEAREST_NORMALIZED_CONTINUATION_ROOT",
            "alternate_roots_retained": True,
        },
        "dynamic_contract": frozen_dynamic_contract_metadata(),
        "solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        "production_solvers": _jsonable((PRODUCTION_RADAU, PRODUCTION_BDF)),
        "dynamic_execution_rule": (
            "ONLY_CONNECTED_FULL_GATE_EXACT_ZERO_CONTINUATIONS"
        ),
        "input_artifact_sha256": dict(dependency_snapshot.input_artifact_sha256),
        "report_artifact_sha256": {
            name: sha256_file(report_directory / name) for name in report_filenames
        },
        "equation_source_file_sha256": source_hashes,
        "equation_source_tree_sha256": sha256_object(source_hashes),
        "comparison_scope": "MODEL_TO_MODEL_DESCRIPTIVE_PROFILE_ONLY",
        "experimental_genotype_comparison_loaded": False,
        "experimental_acceptance_criterion_applied": False,
        "independent_genotype_reproduction_required": True,
        "clean_auditor_authorization_required": True,
        "phenotype_reveal_permitted_by_this_manifest": False,
    }


def _write_hashes(results_directory: Path, filenames: Sequence[str]) -> None:
    rows = [
        {"artifact": name, "sha256": sha256_file(results_directory / name)}
        for name in sorted(filenames)
    ]
    _write_csv(
        results_directory / HASH_FILENAME,
        rows,
        empty_fields=("artifact", "sha256"),
    )


def _invalidate_genotype_hash_ledger(results_directory: Path) -> None:
    (results_directory / HASH_FILENAME).unlink(missing_ok=True)


def _publish_staged_genotype_freeze(
    staging_directory: Path,
    results_directory: Path,
    repository_root: Path,
    dependency_snapshot: GenotypeDependencySnapshot,
) -> None:
    """Publish reports atomically per file and the authority ledger last."""

    _invalidate_genotype_hash_ledger(results_directory)
    _assert_genotype_snapshot_current(
        dependency_snapshot, repository_root, results_directory
    )
    for filename in GENOTYPE_FREEZE_ARTIFACT_FILENAMES:
        os.replace(staging_directory / filename, results_directory / filename)
    _assert_genotype_snapshot_current(
        dependency_snapshot, repository_root, results_directory
    )
    os.replace(
        staging_directory / HASH_FILENAME,
        results_directory / HASH_FILENAME,
    )
    try:
        _assert_genotype_snapshot_current(
            dependency_snapshot, repository_root, results_directory
        )
    except Exception:
        _invalidate_genotype_hash_ledger(results_directory)
        raise


def run_native_genotype_evaluation(
    repository_root: str | Path,
    results_directory: str | Path,
    *,
    workers: int = 1,
) -> Mapping[str, Any]:
    """Run every passing native H1 root without any first-passer selection."""

    repository = Path(repository_root).resolve()
    output = Path(results_directory).resolve()
    dependency_snapshot = _capture_genotype_dependency_snapshot(
        repository, output
    )
    intake = load_frozen_genotype_intake(output, repository_root=repository)
    _assert_genotype_snapshot_current(
        dependency_snapshot, repository, output
    )
    if dict(intake.input_artifact_hashes) != dict(
        dependency_snapshot.input_artifact_sha256
    ):
        raise RuntimeError("genotype intake hashes disagree with the pre-run snapshot")
    _invalidate_genotype_hash_ledger(output)
    results = _run_roots(intake.roots, workers)
    if tuple(_root_id(result.root) for result in results) != tuple(
        sorted(str(row["root_id"]) for row in intake.group_rows)
    ):
        raise AssertionError("genotype runner did not retain every passing root")

    continuation_rows = _continuation_rows(results)
    trajectory_rows = _trajectory_rows(results)
    profile_rows = _profile_rows(results)
    solver_rows = _solver_rows(results)
    summary = _summary_payload(intake, results)
    with tempfile.TemporaryDirectory(
        prefix=".native_genotype.", dir=output
    ) as temporary:
        staging = Path(temporary)
        _write_csv(
            staging / CONTINUATION_FILENAME,
            continuation_rows,
            empty_fields=("root_id", "transporter", "expression"),
        )
        _write_csv(
            staging / TRAJECTORY_FILENAME,
            trajectory_rows,
            empty_fields=("root_id", "transporter", "genotype_role", "time_s"),
        )
        _write_csv(
            staging / PROFILE_FILENAME,
            profile_rows,
            empty_fields=("root_id", "transporter", "dynamic_executed"),
        )
        _write_csv(
            staging / SOLVER_FILENAME,
            solver_rows,
            empty_fields=("root_id", "transporter", "dynamic_executed"),
        )
        _write_json(staging / SUMMARY_FILENAME, summary)
        manifest = _manifest_payload(
            repository,
            output,
            intake,
            results,
            GENOTYPE_REPORT_FILENAMES,
            report_directory=staging,
            dependency_snapshot=dependency_snapshot,
        )
        _write_json(staging / MANIFEST_FILENAME, manifest)
        _write_hashes(staging, GENOTYPE_FREEZE_ARTIFACT_FILENAMES)
        _publish_staged_genotype_freeze(
            staging, output, repository, dependency_snapshot
        )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument(
        "--results-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least one")
    summary = run_native_genotype_evaluation(
        args.repository_root,
        args.results_directory,
        workers=args.workers,
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = (
    "CONTINUATION_FILENAME",
    "HASH_FILENAME",
    "MANIFEST_FILENAME",
    "PROFILE_FILENAME",
    "SOLVER_FILENAME",
    "SUMMARY_FILENAME",
    "TRAJECTORY_FILENAME",
    "FrozenGenotypeIntake",
    "NativeGenotypeRootResult",
    "load_frozen_genotype_intake",
    "match_passing_native_roots",
    "passing_dynamic_group_rows",
    "run_native_genotype_evaluation",
)
