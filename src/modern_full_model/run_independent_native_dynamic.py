"""Independently reproduce and classify frozen native-source WT dynamics.

The runner is a fail-closed consumer of the complete production WT freeze.  It
reintegrates the decisive R3 reference member with the independent equation
transcription for every frozen group root, including roots whose production
scientific gate failed.  It also independently rebuilds every per-root dynamic
and geometry gate from the frozen case reports.  Agreement is evidence for the
WT failure or success classification only; this runner never authorizes a
genotype calculation or phenotype reveal.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np

from .independent import (
    IndependentN1NkccRegulation,
    IndependentSecretagogueProtocol,
    independent_physical_time_grid,
    independent_r3_regulation,
)
from .independent_native import IndependentNativeSourceDefinition


REFERENCE_MEMBER_ID = "R3_G125_P21_PROSE_TREFERENCE"
REFERENCE_GAIN_LABEL = "G125_P21_PROSE"
REFERENCE_TIMING_LABEL = "TREFERENCE"
REFERENCE_FULLY_ACTIVATED_MULTIPLIER = 1.25
REFERENCE_TAU_PKA_S = 10.0
REFERENCE_FORWARD_REGULATION_RATE_S = 1.0 / 30.0
REFERENCE_REVERSE_REGULATION_RATE_S = 1.0 / 30.0
EFFECTIVE_CALCIUM_UM = 0.10
NKCC1_FULLY_ACTIVATED_MULTIPLIER = 1.75
ARMS = ("CCH_ONLY", "CCH_IPR")
METHOD_BY_SOLVER = {"production_radau": "Radau", "production_bdf": "BDF"}
MINUTE_TIMES_S = np.arange(60.0, 601.0, 60.0)

ROOT_INPUT = "native_source_wt_roots.csv"
PROFILE_INPUT = "native_dynamic_contract_profile.csv"
GROUP_INPUT = "native_dynamic_contract_group_gate.csv"
SOLVER_INPUT = "native_dynamic_contract_solver_crosschecks.csv"
ISOMORPHISM_INPUT = "native_dynamic_contract_isomorphism.csv"
NEARBY_INPUT = "native_dynamic_contract_nearby.csv"
MANIFEST_INPUT = "native_dynamic_contract_manifest.json"
GATE_INPUT = "native_dynamic_contract_gate.json"
DYNAMIC_HASH_INPUT = "native_dynamic_contract_hashes.csv"

DYNAMIC_REPORT_INPUTS = (
    PROFILE_INPUT,
    GROUP_INPUT,
    SOLVER_INPUT,
    ISOMORPHISM_INPUT,
    NEARBY_INPUT,
)
DYNAMIC_HASH_LEDGER_ARTIFACTS = (
    *DYNAMIC_REPORT_INPUTS,
    MANIFEST_INPUT,
    GATE_INPUT,
)
NATIVE_INPUT_ARTIFACTS = (
    "native_source_map.json",
    "native_source_root_attempts.csv",
    ROOT_INPUT,
    "native_source_wt_summary.json",
    "native_source_panel.csv",
    "native_source_scale_one_identity.csv",
    "native_source_hashes.csv",
    "independent_native_roots.csv",
    "independent_native_attempts.csv",
    "independent_native_summary.json",
)
PRODUCTION_EQUATION_SOURCES = (
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
INDEPENDENT_DYNAMIC_SOURCES = (
    "src/modern_full_model/independent.py",
    "src/modern_full_model/independent_native.py",
    "src/modern_full_model/run_independent_native_dynamic.py",
)
FROZEN_RESULT_INPUTS = (
    *DYNAMIC_HASH_LEDGER_ARTIFACTS,
    DYNAMIC_HASH_INPUT,
    *NATIVE_INPUT_ARTIFACTS,
)
FIREWALL_HASH_KEYS = frozenset(
    {
        "calibration_ledger_sha256",
        "validation_ledger_sha256",
        "sealed_phenotype_target_ledger_sha256",
        "pre_reveal_log_sha256",
    }
)

EXPECTED_GATE_ID = "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1"
PASS_GATE_STATUS = "SCIENTIFIC_GATE_PASS_AUDITOR_AUTHORIZATION_REQUIRED"
FAIL_GATE_STATUS = "SCIENTIFIC_GATE_FAIL"
# Backward-compatible name used by frozen-intake fixtures representing a
# scientifically passing pre-reveal contract.  Runtime intake accepts the
# explicit PASS/FAIL constants above and does not rely on this alias.
EXPECTED_GATE_STATUS = PASS_GATE_STATUS
EXPECTED_MANIFEST_ID = "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1"
EXPECTED_MANIFEST_STATUS = "FROZEN_FOR_CLEAN_ADVERSARIAL_REVIEW"

EXPECTED_PROFILE_PER_ROOT = 80
EXPECTED_NEARBY_PER_ROOT = 32
EXPECTED_SOLVER_PER_ROOT = 60
EXPECTED_ISOMORPHISM_PER_ROOT = 24
EXPECTED_ARM_RATIOS_PER_ROOT = 56
EXPECTED_SCALE_SAMPLES_PER_ROOT = 1120
NEARBY_CONDITIONS = frozenset(
    {
        "nearby_cell_volume_down",
        "nearby_cell_volume_up",
        "nearby_cell_nacl_down",
        "nearby_cell_nacl_up",
        "nearby_lumen_volume_down",
        "nearby_lumen_volume_up",
        "nearby_lumen_nacl_down",
        "nearby_lumen_nacl_up",
    }
)
WT_FLOW_LOWER_UL_MIN = 9.0
WT_FLOW_UPPER_UL_MIN = 10.0
ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S = 51_469_231 * 6.0e-5
COSTIM_TO_CCH_TOTAL_RATIO_MAX = 1.10
SOLVER_RELATIVE_TOLERANCE = 1.0e-4
ISOMORPHISM_RELATIVE_TOLERANCE = 1.0e-10
SUSTAINMENT_RATIO_BOUNDS = (0.8, 1.2)

OUTPUT_FILENAME = "independent_native_dynamic.csv"
GROUP_OUTPUT_FILENAME = "independent_native_dynamic_group_gate.csv"
OUTPUT_HASH_FILENAME = "independent_native_dynamic_hashes.csv"
SUMMARY_FILENAME = "independent_native_dynamic_summary.json"
OUTPUT_HASH_ARTIFACTS = (OUTPUT_FILENAME, GROUP_OUTPUT_FILENAME)


@dataclass(frozen=True)
class FrozenIndependentDynamicIntake:
    repository: Path
    results_directory: Path
    root_ids: tuple[str, ...]
    production_passing_root_ids: tuple[str, ...]
    roots: tuple[Mapping[str, str], ...]
    group_rows: tuple[Mapping[str, str], ...]
    profile_rows: tuple[Mapping[str, str], ...]
    nearby_rows: tuple[Mapping[str, str], ...]
    solver_rows: tuple[Mapping[str, str], ...]
    isomorphism_rows: tuple[Mapping[str, str], ...]
    regulatory_member_ids: tuple[str, ...]
    isomorphism_pairs: tuple[tuple[str, str], ...]
    input_artifact_sha256: Mapping[str, str]
    production_equation_source_sha256: Mapping[str, str]
    production_equation_source_tree_sha256: str
    independent_source_file_sha256: Mapping[str, str]
    independent_source_tree_sha256: str
    opaque_firewall_ledger_sha256: Mapping[str, str]

    @property
    def passing_root_ids(self) -> tuple[str, ...]:
        return self.production_passing_root_ids

    @property
    def profile(self) -> tuple[Mapping[str, str], ...]:
        """Reference-member rows used by the independent reproduction."""
        passing = set(self.production_passing_root_ids)
        return tuple(
            row for row in self.profile_rows
            if row.get("root_id") in passing
            and row.get("regulatory_member_id") == REFERENCE_MEMBER_ID
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_object(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _reject_duplicate_json_keys(
    pairs: Sequence[tuple[str, Any]],
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


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
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(
        f"{field} must be canonical case-sensitive True/False; found {value!r}"
    )


def _csv_int(value: Any, *, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be an integer")
    try:
        parsed = int(str(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be an integer") from exc
    if str(value) != str(parsed) or parsed < 0:
        raise ValueError(f"{field} must be a canonical nonnegative integer")
    return parsed


def _finite_float(value: Any, *, field: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{field} must be finite")
    return parsed


def _require_json_bool(
    payload: Mapping[str, Any], field: str, expected: bool, *, artifact: str
) -> None:
    value = payload.get(field)
    if type(value) is not bool or value is not expected:
        raise ValueError(
            f"{artifact}.{field} must be JSON boolean {expected}, found {value!r}"
        )


def _require_json_integer(
    payload: Mapping[str, Any], field: str, *, artifact: str
) -> int:
    value = payload.get(field)
    if type(value) is not int or value < 0:
        raise ValueError(
            f"{artifact}.{field} must be a nonnegative JSON integer, found {value!r}"
        )
    return value


def _unique_root_ids(
    payload: Mapping[str, Any], field: str, *, artifact: str
) -> tuple[str, ...]:
    value = payload.get(field)
    if not isinstance(value, list):
        raise ValueError(f"{artifact}.{field} must be a JSON list")
    identifiers: list[str] = []
    for item in value:
        if type(item) is not str or not item:
            raise ValueError(f"{artifact}.{field} contains an invalid root ID")
        identifiers.append(item)
    if len(set(identifiers)) != len(identifiers):
        raise ValueError(f"{artifact}.{field} contains duplicate root IDs")
    return tuple(sorted(identifiers))


def _require_digest(value: Any, *, label: str) -> str:
    if (
        type(value) is not str
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _require_exact_current_hash_registry(
    registry: Any,
    filenames: Sequence[str],
    *,
    base: Path,
    label: str,
) -> dict[str, str]:
    if not isinstance(registry, Mapping):
        raise ValueError(f"{label} must be a hash registry")
    expected = set(filenames)
    actual = set(registry)
    if actual != expected:
        raise ValueError(
            f"{label} must have exact keys; "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )
    current: dict[str, str] = {}
    for filename in filenames:
        path = base / filename
        if not path.is_file():
            raise FileNotFoundError(f"{label} target is missing: {filename}")
        expected_digest = _require_digest(
            registry[filename], label=f"{label}.{filename}"
        )
        actual_digest = _sha256(path)
        if expected_digest != actual_digest:
            raise ValueError(f"{label} has a stale hash for {filename}")
        current[filename] = actual_digest
    return current


def _load_dynamic_hash_ledger(results_directory: Path) -> dict[str, str]:
    path = results_directory / DYNAMIC_HASH_INPUT
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != ("artifact", "sha256"):
            raise ValueError("dynamic hash ledger requires exact artifact,sha256 columns")
        rows = list(reader)
    registry: dict[str, str] = {}
    for row in rows:
        artifact = row["artifact"]
        if artifact in registry:
            raise ValueError(f"dynamic hash ledger duplicates {artifact}")
        registry[artifact] = row["sha256"]
    expected = set(DYNAMIC_HASH_LEDGER_ARTIFACTS)
    if set(registry) != expected or len(rows) != len(expected):
        raise ValueError("dynamic hash ledger must contain exactly seven artifacts")
    return _require_exact_current_hash_registry(
        registry,
        DYNAMIC_HASH_LEDGER_ARTIFACTS,
        base=results_directory,
        label=DYNAMIC_HASH_INPUT,
    )


def _require_firewall_registry(registry: Any) -> dict[str, str]:
    if not isinstance(registry, Mapping) or set(registry) != FIREWALL_HASH_KEYS:
        raise ValueError("opaque firewall registry must contain exactly four hashes")
    return {
        str(label): _require_digest(value, label=f"opaque_firewall.{label}")
        for label, value in registry.items()
    }


def _definition(root: Mapping[str, str]) -> IndependentNativeSourceDefinition:
    return IndependentNativeSourceDefinition.from_root_row(root)


def _group_root_ids(
    rows: Sequence[Mapping[str, str]],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    all_ids: list[str] = []
    passing_ids: list[str] = []
    for row in rows:
        root_id = str(row.get("root_id", ""))
        if not root_id:
            raise ValueError("dynamic group row lacks root_id")
        all_ids.append(root_id)
        if not _csv_bool(
            row.get("all_native_rest_wt_gate_pass"),
            field=f"{root_id}.all_native_rest_wt_gate_pass",
        ):
            raise ValueError(f"{root_id}: group root does not pass frozen WT rest")
        scientific_pass = _csv_bool(
            row.get("native_dynamic_scientific_gate_pass"),
            field=f"{root_id}.native_dynamic_scientific_gate_pass",
        )
        if row.get("hydraulic_mode") != "H1":
            raise ValueError(f"{root_id}: group root is not H1")
        if _finite_float(
            row.get("hydraulic_scale"), field=f"{root_id}.hydraulic_scale"
        ) != 1.0:
            raise ValueError(f"{root_id}: group root is not hydraulic scale one")
        if row.get("eligibility") != "PRODUCTION_PRE_REVEAL":
            raise ValueError(f"{root_id}: group eligibility is not exact")
        if scientific_pass:
            passing_ids.append(root_id)
    if len(set(all_ids)) != len(all_ids):
        raise ValueError("dynamic group table contains duplicate root IDs")
    if not all_ids:
        raise RuntimeError("dynamic group table contains no frozen roots")
    return tuple(sorted(all_ids)), tuple(sorted(passing_ids))


def _parse_isomorphism_pairs(value: Any) -> tuple[tuple[str, str], ...]:
    if not isinstance(value, list):
        raise ValueError("manifest R2/R3 pairs must be a JSON list")
    pairs: list[tuple[str, str]] = []
    for item in value:
        if (
            not isinstance(item, list)
            or len(item) != 2
            or any(type(member_id) is not str or not member_id for member_id in item)
        ):
            raise ValueError("manifest contains a malformed R2/R3 pair")
        pairs.append((item[0], item[1]))
    if len(pairs) != 6 or len(set(pairs)) != 6:
        raise ValueError("manifest must freeze six unique R2/R3 pairs")
    if len({left for left, _ in pairs}) != 6:
        raise ValueError("manifest R2 isomorphism axis must contain six unique IDs")
    if len({right for _, right in pairs}) != 6:
        raise ValueError("manifest R3 isomorphism axis must contain six unique IDs")
    return tuple(pairs)


def _validate_regulatory_members(
    members: Any,
    pairs: Sequence[tuple[str, str]],
) -> tuple[str, ...]:
    """Validate the exact frozen 2/6/6/6 gain/timing member registry."""

    if not isinstance(members, Mapping) or len(members) != 20:
        raise ValueError("production manifest must freeze exactly 20 regulatory members")
    member_ids = tuple(str(member_id) for member_id in members)
    if len(set(member_ids)) != 20 or REFERENCE_MEMBER_ID not in member_ids:
        raise ValueError("production regulatory member registry is malformed")

    family_counts: Counter[str] = Counter()
    axes: dict[str, set[tuple[str, str]]] = {
        family: set() for family in ("R0", "R1", "R2", "R3")
    }
    for member_id, payload in members.items():
        if not isinstance(payload, Mapping):
            raise ValueError(f"regulatory member {member_id} payload is malformed")
        if payload.get("member_id") != member_id:
            raise ValueError(f"regulatory member {member_id} has a stale member_id")
        family = payload.get("family")
        gain = payload.get("gain_label")
        timing = payload.get("timing_label")
        if family not in axes or type(gain) is not str or type(timing) is not str:
            raise ValueError(f"regulatory member {member_id} has malformed axes")
        family_counts[family] += 1
        axis = (gain, timing)
        if axis in axes[family]:
            raise ValueError(f"regulatory family {family} duplicates axis {axis}")
        axes[family].add(axis)

    expected_gains = {"G125_P21_PROSE", "G170_P21_FIGURE_VISUAL"}
    expected_timings = {"TFAST", "TREFERENCE", "TSLOW"}
    expected_dynamic_axes = {
        (gain, timing) for gain in expected_gains for timing in expected_timings
    }
    if family_counts != Counter({"R0": 2, "R1": 6, "R2": 6, "R3": 6}):
        raise ValueError("regulatory family composition is not exact 2/6/6/6")
    if axes["R0"] != {(gain, "STATIC") for gain in expected_gains}:
        raise ValueError("R0 does not contain the exact two static controls")
    for family in ("R1", "R2", "R3"):
        if axes[family] != expected_dynamic_axes:
            raise ValueError(f"{family} does not span the exact gain/timing grid")

    reference = members[REFERENCE_MEMBER_ID]
    if (
        reference.get("family") != "R3"
        or reference.get("gain_label") != REFERENCE_GAIN_LABEL
        or reference.get("timing_label") != REFERENCE_TIMING_LABEL
        or _finite_float(
            reference.get("fully_activated_multiplier"),
            field=f"{REFERENCE_MEMBER_ID}.fully_activated_multiplier",
        )
        != REFERENCE_FULLY_ACTIVATED_MULTIPLIER
    ):
        raise ValueError("frozen independent reference member is not exact")

    for r2_id, r3_id in pairs:
        if r2_id not in members or r3_id not in members:
            raise ValueError("R2/R3 pair references an unknown member")
        r2 = members[r2_id]
        r3 = members[r3_id]
        if (
            r2.get("family") != "R2"
            or r3.get("family") != "R3"
            or r2.get("gain_label") != r3.get("gain_label")
            or r2.get("timing_label") != r3.get("timing_label")
        ):
            raise ValueError("R2/R3 pair does not bind an identical gain/timing axis")
    return member_ids


def _current_hash_registry(
    filenames: Sequence[str], *, base: Path, label: str
) -> dict[str, str]:
    registry: dict[str, str] = {}
    for filename in filenames:
        path = base / filename
        if not path.is_file():
            raise FileNotFoundError(f"{label} target is missing: {filename}")
        registry[filename] = _sha256(path)
    return registry


def _require_unchanged_hash_registry(
    frozen: Mapping[str, str], filenames: Sequence[str], *, base: Path, label: str
) -> None:
    if set(frozen) != set(filenames):
        raise ValueError(f"{label} frozen registry keys changed")
    current = _current_hash_registry(filenames, base=base, label=label)
    if current != dict(frozen):
        changed = sorted(
            filename for filename in filenames if current.get(filename) != frozen.get(filename)
        )
        raise ValueError(f"{label} drifted during independent reproduction: {changed}")


def _exact_case_keys(
    rows: Sequence[Mapping[str, str]],
    fields: Sequence[str],
    expected: set[tuple[str, ...]],
    *,
    label: str,
) -> None:
    actual: list[tuple[str, ...]] = []
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in fields)
        if not all(key):
            raise ValueError(f"{label} contains an incomplete case key")
        actual.append(key)
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError(f"{label} does not contain its exact frozen case keys")


def load_frozen_independent_dynamic_intake(
    repository: str | Path,
) -> FrozenIndependentDynamicIntake:
    """Load all frozen roots, including an exact empty production PASS set."""

    repository_path = Path(repository).resolve()
    results_directory = repository_path / "results" / "13B_modern_full_model"
    missing = [
        filename
        for filename in FROZEN_RESULT_INPUTS
        if not (results_directory / filename).is_file()
    ]
    if missing:
        raise FileNotFoundError(
            f"independent dynamic intake lacks frozen artifacts: {missing}"
        )

    group_rows = _read(results_directory / GROUP_INPUT)
    root_ids, production_passing_ids = _group_root_ids(group_rows)
    group_registry = {str(row["root_id"]): row for row in group_rows}
    gate = _read_json(results_directory / GATE_INPUT)
    manifest = _read_json(results_directory / MANIFEST_INPUT)

    expected_status = PASS_GATE_STATUS if production_passing_ids else FAIL_GATE_STATUS
    if gate.get("gate_id") != EXPECTED_GATE_ID:
        raise ValueError("production dynamic gate_id is not frozen")
    if gate.get("status") != expected_status:
        raise ValueError("production gate status disagrees with its passing root IDs")
    if gate.get("manifest") != MANIFEST_INPUT:
        raise ValueError("production gate does not bind the expected manifest")
    for field, expected in (
        ("at_least_one_scientific_group_pass", bool(production_passing_ids)),
        ("all_rows_retained_no_first_passer", True),
        ("clean_auditor_authorization_required", True),
        ("clean_auditor_authorization_recorded", False),
        ("independent_numerical_reproduction_required", True),
        ("independent_numerical_reproduction_recorded_by_this_runner", False),
        ("final_composite_authorization_complete", False),
        ("phenotype_reveal_permitted", False),
    ):
        _require_json_bool(gate, field, expected, artifact=GATE_INPUT)
    if _unique_root_ids(
        gate, "scientific_passing_root_ids", artifact=GATE_INPUT
    ) != production_passing_ids:
        raise ValueError("gate and group table disagree on scientific passing IDs")
    if _require_json_integer(
        gate, "scientific_passing_group_count", artifact=GATE_INPUT
    ) != len(production_passing_ids):
        raise ValueError("gate and group table disagree on scientific passing count")
    gate_root_count = _require_json_integer(
        gate, "native_root_count", artifact=GATE_INPUT
    )

    if manifest.get("manifest_id") != EXPECTED_MANIFEST_ID:
        raise ValueError("production dynamic manifest_id is not frozen")
    if manifest.get("status") != EXPECTED_MANIFEST_STATUS:
        raise ValueError("production dynamic manifest is not frozen for review")
    for field, expected in (
        ("clean_auditor_authorization_required", True),
        (
            "independent_numerical_reproduction_required_for_composite_authorization",
            True,
        ),
        ("phenotype_reveal_permitted_by_this_manifest", False),
    ):
        _require_json_bool(manifest, field, expected, artifact=MANIFEST_INPUT)
    if _unique_root_ids(
        manifest, "scientific_passing_root_ids", artifact=MANIFEST_INPUT
    ) != production_passing_ids:
        raise ValueError("manifest and group table disagree on scientific passing IDs")
    if _require_json_integer(
        manifest, "scientific_passing_group_count", artifact=MANIFEST_INPUT
    ) != len(production_passing_ids):
        raise ValueError("manifest and group table disagree on scientific passing count")
    manifest_roots = manifest.get("roots")
    if not isinstance(manifest_roots, Mapping):
        raise ValueError("production dynamic manifest lacks its root registry")
    manifest_root_ids = tuple(sorted(str(root_id) for root_id in manifest_roots))
    manifest_root_count = _require_json_integer(
        manifest, "root_count", artifact=MANIFEST_INPUT
    )
    if (
        manifest_root_ids != root_ids
        or manifest_root_count != len(manifest_roots)
        or manifest_root_count != len(group_rows)
        or gate_root_count != manifest_root_count
    ):
        raise ValueError("gate, manifest, and group root registries are not identical")

    _load_dynamic_hash_ledger(results_directory)
    _require_exact_current_hash_registry(
        manifest.get("wt_dynamic_report_artifact_sha256"),
        DYNAMIC_REPORT_INPUTS,
        base=results_directory,
        label="wt_dynamic_report_artifact_sha256",
    )
    _require_exact_current_hash_registry(
        manifest.get("native_input_artifact_sha256"),
        NATIVE_INPUT_ARTIFACTS,
        base=results_directory,
        label="native_input_artifact_sha256",
    )
    equation_hashes = _require_exact_current_hash_registry(
        manifest.get("equation_source_file_sha256"),
        PRODUCTION_EQUATION_SOURCES,
        base=repository_path,
        label="equation_source_file_sha256",
    )
    tree_hash = _require_digest(
        manifest.get("equation_source_tree_sha256"),
        label="equation_source_tree_sha256",
    )
    if tree_hash != _sha256_object(equation_hashes):
        raise ValueError("production equation source tree hash is stale")
    firewall_hashes = _require_firewall_registry(
        manifest.get("opaque_firewall_ledger_sha256")
    )
    if manifest.get("opaque_firewall_access_mode") != (
        "SHA256_BYTES_ONLY_NO_TARGET_PARSE"
    ):
        raise ValueError("production manifest firewall mode is not frozen")

    native_summary = _read_json(results_directory / "native_source_wt_summary.json")
    for field in (
        "production_confirmation_complete",
        "geometry_confirmation_complete",
        "independent_root_reproduction_complete",
        "native_pre_reveal_freeze_ready",
    ):
        _require_json_bool(
            native_summary,
            field,
            True,
            artifact="native_source_wt_summary.json",
        )

    root_rows = _read(results_directory / ROOT_INPUT)
    root_registry: dict[str, Mapping[str, str]] = {}
    for row in root_rows:
        root_id = str(row.get("root_id", ""))
        if not root_id:
            raise ValueError("native root table contains an empty root ID")
        if root_id in root_registry:
            raise ValueError(f"native root table duplicates {root_id}")
        root_registry[root_id] = row
    selected_roots: list[Mapping[str, str]] = []
    # Every retained native root is evidence.  Scientific passing IDs are a
    # classification to reproduce, never a selector for independent intake.
    for root_id in root_ids:
        if root_id not in root_registry:
            raise ValueError(f"frozen group root is missing from native roots: {root_id}")
        root = root_registry[root_id]
        for field in ("passes_numerical_gate", "passes_wt_gate"):
            if not _csv_bool(root.get(field), field=f"{root_id}.{field}"):
                raise ValueError(f"{root_id}: frozen root fails {field}")
        if root.get("eligibility") != "PRODUCTION_PRE_REVEAL":
            raise ValueError(f"{root_id}: root eligibility is not exact")
        if root.get("hydraulic_mode") != "H1":
            raise ValueError(f"{root_id}: root hydraulic mode is not H1")
        if _finite_float(
            root.get("hydraulic_scale"), field=f"{root_id}.hydraulic_scale"
        ) != 1.0:
            raise ValueError(f"{root_id}: root hydraulic scale is not one")
        try:
            core = tuple(float(value) for value in json.loads(root["core_state_json"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{root_id}: malformed native root state") from exc
        if len(core) != 12 or not all(
            math.isfinite(value) and value > 0.0 for value in core
        ):
            raise ValueError(f"{root_id}: core state must have 12 positive values")
        _definition(root)
        root_payload = manifest_roots.get(root_id)
        if not isinstance(root_payload, Mapping):
            raise ValueError(f"{root_id}: manifest root payload is malformed")
        if root_payload.get("core_state_sha256") != _sha256_object(core):
            raise ValueError(f"{root_id}: core-state hash is stale")
        selected_roots.append(root)

    raw_pairs = manifest.get("r2_r3_isomorphism_pairs")
    pairs = _parse_isomorphism_pairs(raw_pairs)
    members = manifest.get("regulatory_members")
    member_ids = _validate_regulatory_members(members, pairs)

    profile_rows = _read(results_directory / PROFILE_INPUT)
    expected_profile_keys = {
        (root_id, member_id, arm, solver)
        for root_id in root_ids
        for member_id in member_ids
        for arm in ARMS
        for solver in METHOD_BY_SOLVER
    }
    _exact_case_keys(
        profile_rows,
        ("root_id", "regulatory_member_id", "arm", "solver_label"),
        expected_profile_keys,
        label="production profile",
    )
    nearby_rows = _read(results_directory / NEARBY_INPUT)
    expected_nearby_keys = {
        (root_id, condition, arm, solver)
        for root_id in root_ids
        for condition in NEARBY_CONDITIONS
        for arm in ARMS
        for solver in METHOD_BY_SOLVER
    }
    nearby_fields = ("root_id", "initial_condition", "arm", "solver_label")
    _exact_case_keys(
        nearby_rows, nearby_fields, expected_nearby_keys, label="nearby profile"
    )
    solver_rows = _read(results_directory / SOLVER_INPUT)
    expected_solver_keys = {
        (
            root_id,
            "RADAU_VS_BDF",
            member_id,
            arm,
            "native_source_rest_root",
            "production_radau",
            "production_bdf",
        )
        for root_id in root_ids
        for member_id in member_ids
        for arm in ARMS
    }
    expected_solver_keys.update(
        {
            (
                root_id,
                "NEARBY_RADAU_VS_BDF",
                REFERENCE_MEMBER_ID,
                arm,
                condition,
                "production_radau",
                "production_bdf",
            )
            for root_id in root_ids
            for condition in NEARBY_CONDITIONS
            for arm in ARMS
        }
    )
    expected_solver_keys.update(
        {
            (
                root_id,
                f"PRODUCTION_VS_{solver.upper()}",
                REFERENCE_MEMBER_ID,
                arm,
                "native_source_rest_root",
                "production_radau",
                solver,
            )
            for root_id in root_ids
            for arm in ARMS
            for solver in ("loose_radau", "tight_radau")
        }
    )
    solver_fields = (
            "root_id",
            "comparison_kind",
            "regulatory_member_id",
            "arm",
            "initial_condition",
            "reference_solver_label",
            "comparison_solver_label",
        )
    _exact_case_keys(
        solver_rows, solver_fields, expected_solver_keys,
        label="solver comparisons",
    )
    isomorphism_rows = _read(results_directory / ISOMORPHISM_INPUT)
    expected_isomorphism_keys = {
        (root_id, r2_id, r3_id, arm, solver)
        for root_id in root_ids
        for r2_id, r3_id in pairs
        for arm in ARMS
        for solver in METHOD_BY_SOLVER
    }
    isomorphism_fields = (
        "root_id", "r2_member_id", "r3_member_id", "arm", "solver_label"
    )
    _exact_case_keys(
        isomorphism_rows, isomorphism_fields, expected_isomorphism_keys,
        label="R2/R3 isomorphism comparisons",
    )

    input_hashes = {
        filename: _sha256(results_directory / filename)
        for filename in FROZEN_RESULT_INPUTS
    }
    independent_source_hashes = _current_hash_registry(
        INDEPENDENT_DYNAMIC_SOURCES,
        base=repository_path,
        label="independent_source_file_sha256",
    )
    return FrozenIndependentDynamicIntake(
        repository=repository_path,
        results_directory=results_directory,
        root_ids=root_ids,
        production_passing_root_ids=production_passing_ids,
        roots=tuple(selected_roots),
        group_rows=tuple(group_registry[root_id] for root_id in root_ids),
        profile_rows=tuple(profile_rows),
        nearby_rows=tuple(nearby_rows),
        solver_rows=tuple(solver_rows),
        isomorphism_rows=tuple(isomorphism_rows),
        regulatory_member_ids=member_ids,
        isomorphism_pairs=pairs,
        input_artifact_sha256=input_hashes,
        production_equation_source_sha256=equation_hashes,
        production_equation_source_tree_sha256=tree_hash,
        independent_source_file_sha256=independent_source_hashes,
        independent_source_tree_sha256=_sha256_object(independent_source_hashes),
        opaque_firewall_ledger_sha256=firewall_hashes,
    )


def _optional_finite(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _relative_close(left: float | None, right: float | None) -> bool:
    if left is None or right is None:
        return left is None and right is None
    return bool(
        abs(left - right) / max(abs(left), abs(right), 1.0e-12) <= 1.0e-10
    )


def _recomputed_numerical_gate(
    row: Mapping[str, Any], *, label: str
) -> tuple[bool, bool, bool]:
    """Return success, primitive numerical gate, and report-flag agreement."""

    success = _csv_bool(row.get("success"), field=f"{label}.success")
    positive = _csv_bool(
        row.get("positive_core"), field=f"{label}.positive_core"
    )
    nonnegative = _csv_bool(
        row.get("all_flow_nonnegative"),
        field=f"{label}.all_flow_nonnegative",
    )
    conservation_ratio = _optional_finite(
        row.get("max_dimensionless_conservation_ratio")
    )
    recomputed = bool(
        success
        and positive
        and nonnegative
        and conservation_ratio is not None
        and conservation_ratio <= 1.0
    )
    declared = _csv_bool(
        row.get("numerical_gate_pass"), field=f"{label}.numerical_gate_pass"
    )
    return success, recomputed, declared == recomputed


def _recomputed_sustainment_gate(
    row: Mapping[str, Any], *, label: str
) -> tuple[bool, bool]:
    ratio = _optional_finite(row.get("sustainment_ratio_q600_q60"))
    lower, upper = SUSTAINMENT_RATIO_BOUNDS
    recomputed = bool(ratio is not None and lower <= ratio <= upper)
    declared = _csv_bool(
        row.get("sustainment_gate_pass"),
        field=f"{label}.sustainment_gate_pass",
    )
    return recomputed, declared == recomputed


def _recomputed_solver_gate(
    row: Mapping[str, Any], *, label: str
) -> tuple[bool, bool]:
    both_success = _csv_bool(
        row.get("both_success"), field=f"{label}.both_success"
    )
    reference_numerical = _csv_bool(
        row.get("reference_numerical_gate_pass"),
        field=f"{label}.reference_numerical_gate_pass",
    )
    comparison_numerical = _csv_bool(
        row.get("comparison_numerical_gate_pass"),
        field=f"{label}.comparison_numerical_gate_pass",
    )
    metrics = tuple(
        _optional_finite(row.get(field))
        for field in (
            "max_relative_state_difference",
            "max_relative_flow_difference",
            "endpoint_cumulative_flow_relative_difference",
        )
    )
    recomputed = bool(
        both_success
        and reference_numerical
        and comparison_numerical
        and all(metric is not None for metric in metrics)
        and max((float(metric) for metric in metrics if metric is not None), default=math.inf)
        <= SOLVER_RELATIVE_TOLERANCE
    )
    declared = _csv_bool(
        row.get("solver_gate_pass"), field=f"{label}.solver_gate_pass"
    )
    return recomputed, declared == recomputed


def _recomputed_isomorphism_gate(
    row: Mapping[str, Any], *, label: str
) -> tuple[bool, bool]:
    metrics = tuple(
        _optional_finite(row.get(field))
        for field in (
            "max_relative_state_difference",
            "max_relative_flow_difference",
            "max_relative_capacity_difference",
        )
    )
    recomputed = bool(
        all(metric is not None for metric in metrics)
        and max((float(metric) for metric in metrics if metric is not None), default=math.inf)
        <= ISOMORPHISM_RELATIVE_TOLERANCE
    )
    declared = _csv_bool(
        row.get("r2_r3_isomorphism_gate_pass"),
        field=f"{label}.r2_r3_isomorphism_gate_pass",
    )
    return recomputed, declared == recomputed


def _recompute_group_classifications(
    intake: FrozenIndependentDynamicIntake,
) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    """Rebuild every production group/geometry gate from frozen case rows."""

    group_by_root = {str(row["root_id"]): row for row in intake.group_rows}
    profile_by_root = {
        root_id: [row for row in intake.profile_rows if row["root_id"] == root_id]
        for root_id in intake.root_ids
    }
    nearby_by_root = {
        root_id: [row for row in intake.nearby_rows if row["root_id"] == root_id]
        for root_id in intake.root_ids
    }
    solver_by_root = {
        root_id: [row for row in intake.solver_rows if row["root_id"] == root_id]
        for root_id in intake.root_ids
    }
    iso_by_root = {
        root_id: [row for row in intake.isomorphism_rows if row["root_id"] == root_id]
        for root_id in intake.root_ids
    }
    roots_by_id = {str(row["root_id"]): row for row in intake.roots}

    output: list[dict[str, Any]] = []
    independent_passing_ids: list[str] = []
    for root_id in intake.root_ids:
        group = group_by_root[root_id]
        profile = profile_by_root[root_id]
        nearby = nearby_by_root[root_id]
        solver = solver_by_root[root_id]
        isomorphism = iso_by_root[root_id]
        profile_numerical_rows = [
            _recomputed_numerical_gate(
                row, label=f"{root_id}.profile[{index}]"
            )
            for index, row in enumerate(profile)
        ]
        nearby_numerical_rows = [
            _recomputed_numerical_gate(
                row, label=f"{root_id}.nearby[{index}]"
            )
            for index, row in enumerate(nearby)
        ]
        profile_sustainment_rows = [
            _recomputed_sustainment_gate(
                row, label=f"{root_id}.profile[{index}]"
            )
            for index, row in enumerate(profile)
        ]
        nearby_sustainment_rows = [
            _recomputed_sustainment_gate(
                row, label=f"{root_id}.nearby[{index}]"
            )
            for index, row in enumerate(nearby)
        ]
        solver_gate_rows = [
            _recomputed_solver_gate(
                row, label=f"{root_id}.solver[{index}]"
            )
            for index, row in enumerate(solver)
        ]
        isomorphism_gate_rows = [
            _recomputed_isomorphism_gate(
                row, label=f"{root_id}.isomorphism[{index}]"
            )
            for index, row in enumerate(isomorphism)
        ]
        profile_success = [row[0] for row in profile_numerical_rows]
        nearby_success = [row[0] for row in nearby_numerical_rows]
        production_numerical = bool(
            len(profile) == EXPECTED_PROFILE_PER_ROOT
            and all(row[1] for row in profile_numerical_rows)
        )
        nearby_numerical = bool(
            len(nearby) == EXPECTED_NEARBY_PER_ROOT
            and all(row[1] for row in nearby_numerical_rows)
        )
        production_sustainment = bool(
            len(profile) == EXPECTED_PROFILE_PER_ROOT
            and all(row[0] for row in profile_sustainment_rows)
        )
        nearby_sustainment = bool(
            len(nearby) == EXPECTED_NEARBY_PER_ROOT
            and all(row[0] for row in nearby_sustainment_rows)
        )
        solver_gate = bool(
            len(solver) == EXPECTED_SOLVER_PER_ROOT
            and all(row[0] for row in solver_gate_rows)
        )
        isomorphism_gate = bool(
            len(isomorphism) == EXPECTED_ISOMORPHISM_PER_ROOT
            and all(row[0] for row in isomorphism_gate_rows)
        )
        report_gate_consistency = bool(
            all(row[2] for row in profile_numerical_rows)
            and all(row[2] for row in nearby_numerical_rows)
            and all(row[1] for row in profile_sustainment_rows)
            and all(row[1] for row in nearby_sustainment_rows)
            and all(row[1] for row in solver_gate_rows)
            and all(row[1] for row in isomorphism_gate_rows)
        )

        regulatory_gate = bool(
            len(profile) == EXPECTED_PROFILE_PER_ROOT
            and all(profile_success)
            and all(
                (value := _optional_finite(row.get("capacity_multiplier_start")))
                is not None
                and abs(value - 1.0) <= 1.0e-12
                for row in profile
            )
            and all(
                (value := _optional_finite(row.get("capacity_multiplier_endpoint")))
                is not None
                and abs(value - 1.0) <= 1.0e-10
                for row in profile
                if row["arm"] == "CCH_ONLY"
            )
            and all(
                (value := _optional_finite(row.get("capacity_multiplier_endpoint")))
                is not None
                and value > 1.0
                for row in profile
                if row["arm"] == "CCH_IPR"
            )
        )

        profile_lookup = {
            (row["regulatory_member_id"], row["arm"], row["solver_label"]): row
            for row in profile
        }
        nearby_lookup = {
            (row["initial_condition"], row["arm"], row["solver_label"]): row
            for row in nearby
        }
        arm_ratios: list[float] = []
        arm_complete = True
        for member_id in intake.regulatory_member_ids:
            for solver_label in METHOD_BY_SOLVER:
                cch_row = profile_lookup[(member_id, "CCH_ONLY", solver_label)]
                costim_row = profile_lookup[(member_id, "CCH_IPR", solver_label)]
                cch = _optional_finite(cch_row.get("cumulative_flow_600_pL"))
                costim = _optional_finite(costim_row.get("cumulative_flow_600_pL"))
                if (
                    not _csv_bool(
                        cch_row.get("success"),
                        field=f"{root_id}.{member_id}.CCH_ONLY.success",
                    )
                    or not _csv_bool(
                        costim_row.get("success"),
                        field=f"{root_id}.{member_id}.CCH_IPR.success",
                    )
                    or cch is None
                    or costim is None
                    or cch <= 0.0
                    or costim <= 0.0
                ):
                    arm_complete = False
                else:
                    arm_ratios.append(costim / cch)
        for condition in NEARBY_CONDITIONS:
            for solver_label in METHOD_BY_SOLVER:
                cch_row = nearby_lookup[(condition, "CCH_ONLY", solver_label)]
                costim_row = nearby_lookup[(condition, "CCH_IPR", solver_label)]
                cch = _optional_finite(cch_row.get("cumulative_flow_600_pL"))
                costim = _optional_finite(costim_row.get("cumulative_flow_600_pL"))
                if (
                    not _csv_bool(
                        cch_row.get("success"),
                        field=f"{root_id}.{condition}.CCH_ONLY.success",
                    )
                    or not _csv_bool(
                        costim_row.get("success"),
                        field=f"{root_id}.{condition}.CCH_IPR.success",
                    )
                    or cch is None
                    or costim is None
                    or cch <= 0.0
                    or costim <= 0.0
                ):
                    arm_complete = False
                else:
                    arm_ratios.append(costim / cch)
        beta_total_gate = bool(
            arm_complete
            and len(arm_ratios) == EXPECTED_ARM_RATIOS_PER_ROOT
            and max(arm_ratios, default=math.inf) <= COSTIM_TO_CCH_TOTAL_RATIO_MAX
        )

        minute_flows: list[float] = []
        flow_complete = True
        successful_rows = [
            row
            for row, success in zip(profile, profile_success)
            if success
        ] + [
            row
            for row, success in zip(nearby, nearby_success)
            if success
        ]
        for row in successful_rows:
            for minute in range(1, 11):
                value = _optional_finite(row.get(f"flow_minute_{minute}_pL_s"))
                if value is None or value <= 0.0:
                    flow_complete = False
                else:
                    minute_flows.append(value)
        scale_lower: float | None = None
        scale_upper: float | None = None
        feasible_upper: float | None = None
        shared_scale_gate = False
        geometry_gate = False
        if flow_complete and minute_flows:
            scale_lower = max(WT_FLOW_LOWER_UL_MIN / value for value in minute_flows)
            scale_upper = min(WT_FLOW_UPPER_UL_MIN / value for value in minute_flows)
            shared_scale_gate = scale_lower <= scale_upper
            feasible_upper = min(
                scale_upper, ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
            )
            geometry_gate = bool(shared_scale_gate and scale_lower <= feasible_upper)

        source_rest_gate = bool(
            _csv_bool(
                roots_by_id[root_id].get("passes_wt_gate"),
                field=f"{root_id}.passes_wt_gate",
            )
        )
        independent_components = {
            "all_native_rest_wt_gate_pass": source_rest_gate,
            "all_dynamic_numerical_gates_pass": production_numerical,
            "all_nearby_numerical_gates_pass": nearby_numerical,
            "all_sustainment_gates_pass": production_sustainment,
            "all_nearby_sustainment_gates_pass": nearby_sustainment,
            "all_solver_and_tolerance_gates_pass": solver_gate,
            "r2_r3_isomorphism_gate_pass": isomorphism_gate,
            "p21_regulatory_direction_contract_pass": regulatory_gate,
            "costim_to_cch_total_ratio_gate_pass": beta_total_gate,
            "shared_scale_interval_nonempty": shared_scale_gate,
            "one_smg_geometry_gate_pass": geometry_gate,
        }
        independent_scientific_gate = all(independent_components.values())
        production_components = {
            field: _csv_bool(
                group.get(field), field=f"{root_id}.group.{field}"
            )
            for field in independent_components
        }
        production_scientific_gate = _csv_bool(
            group.get("native_dynamic_scientific_gate_pass"),
            field=f"{root_id}.group.native_dynamic_scientific_gate_pass",
        )
        component_disagreements = sorted(
            field
            for field in independent_components
            if independent_components[field] != production_components[field]
        )
        count_agreement = bool(
            _csv_int(
                group.get("production_case_count"),
                field=f"{root_id}.production_case_count",
            )
            == len(profile)
            == EXPECTED_PROFILE_PER_ROOT
            and _csv_int(
                group.get("nearby_case_count"),
                field=f"{root_id}.nearby_case_count",
            )
            == len(nearby)
            == EXPECTED_NEARBY_PER_ROOT
            and _csv_int(
                group.get("solver_comparison_count"),
                field=f"{root_id}.solver_comparison_count",
            )
            == len(solver)
            == EXPECTED_SOLVER_PER_ROOT
            and _csv_int(
                group.get("isomorphism_comparison_count"),
                field=f"{root_id}.isomorphism_comparison_count",
            )
            == len(isomorphism)
            == EXPECTED_ISOMORPHISM_PER_ROOT
            and _csv_int(
                group.get("paired_arm_count"),
                field=f"{root_id}.paired_arm_count",
            )
            == len(arm_ratios)
            and _csv_int(
                group.get("shared_scale_sample_count"),
                field=f"{root_id}.shared_scale_sample_count",
            )
            == len(minute_flows)
        )
        scale_metric_agreement = bool(
            _relative_close(
                scale_lower,
                _optional_finite(group.get("shared_scale_lower_uL_min_per_pL_s")),
            )
            and _relative_close(
                scale_upper,
                _optional_finite(group.get("shared_scale_upper_uL_min_per_pL_s")),
            )
            and _relative_close(
                scale_lower,
                _optional_finite(
                    group.get("one_smg_feasible_scale_lower_uL_min_per_pL_s")
                ),
            )
            and _relative_close(
                feasible_upper,
                _optional_finite(
                    group.get("one_smg_feasible_scale_upper_uL_min_per_pL_s")
                ),
            )
        )
        scientific_agreement = (
            independent_scientific_gate == production_scientific_gate
        )
        classification_agreement = bool(
            not component_disagreements
            and count_agreement
            and scale_metric_agreement
            and scientific_agreement
            and report_gate_consistency
        )
        if independent_scientific_gate:
            independent_passing_ids.append(root_id)
        output.append(
            {
                "root_id": root_id,
                "production_scientific_gate_pass": production_scientific_gate,
                "independent_scientific_gate_pass": independent_scientific_gate,
                "scientific_gate_agreement": scientific_agreement,
                "production_geometry_gate_pass": production_components[
                    "one_smg_geometry_gate_pass"
                ],
                "independent_geometry_gate_pass": geometry_gate,
                "geometry_gate_agreement": (
                    production_components["one_smg_geometry_gate_pass"]
                    == geometry_gate
                ),
                "independent_shared_scale_lower_uL_min_per_pL_s": scale_lower,
                "independent_shared_scale_upper_uL_min_per_pL_s": scale_upper,
                "independent_one_smg_feasible_upper_uL_min_per_pL_s": (
                    feasible_upper
                ),
                "component_gate_disagreements_json": json.dumps(
                    component_disagreements, separators=(",", ":")
                ),
                "count_contract_agreement": count_agreement,
                "scale_metric_agreement": scale_metric_agreement,
                "primitive_report_gate_consistency": report_gate_consistency,
                "all_group_classification_agreement": classification_agreement,
            }
        )
    return output, tuple(sorted(independent_passing_ids))


def _relative(left: float, right: float, floor: float = 1.0e-12) -> float:
    return abs(left - right) / max(abs(left), abs(right), floor)


def _endpoint_concentrations(state: np.ndarray) -> tuple[float, float, float]:
    volume = float(state[5])
    return float(state[0] / volume), float(state[1] / volume), volume


def reproduce_reference_dynamics(
    intake: FrozenIndependentDynamicIntake,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Reintegrate R3 for every root and reproduce the frozen group decision."""

    roots = {row["root_id"]: row for row in intake.roots}
    profile = sorted(
        (
            row
            for row in intake.profile_rows
            if row["regulatory_member_id"] == REFERENCE_MEMBER_ID
        ),
        key=lambda row: (row["root_id"], row["arm"], row["solver_label"]),
    )
    expected = len(intake.root_ids) * len(ARMS) * len(METHOD_BY_SOLVER)
    if len(profile) != expected:
        raise ValueError(
            f"reference profile has {len(profile)} rows for "
            f"{len(intake.root_ids)} roots; expected {expected}"
        )

    records: list[dict[str, Any]] = []
    trajectories: dict[tuple[str, str, str], Any] = {}
    for production in profile:
        root = roots[production["root_id"]]
        definition = _definition(root)
        timing = json.loads(
            production["regulatory_timing_parameters_json"],
            object_pairs_hook=_reject_duplicate_json_keys,
        )
        if not isinstance(timing, Mapping) or set(timing) != {
            "tau_pka_s",
            "forward_regulation_rate_s",
            "reverse_regulation_rate_s",
        }:
            raise ValueError("reference profile timing declaration is malformed")
        declaration_pass = bool(
            _relative(
                _finite_float(timing["tau_pka_s"], field="reference.tau_pka_s"),
                REFERENCE_TAU_PKA_S,
            )
            <= 1.0e-12
            and _relative(
                _finite_float(
                    timing["forward_regulation_rate_s"],
                    field="reference.forward_regulation_rate_s",
                ),
                REFERENCE_FORWARD_REGULATION_RATE_S,
            )
            <= 1.0e-12
            and _relative(
                _finite_float(
                    timing["reverse_regulation_rate_s"],
                    field="reference.reverse_regulation_rate_s",
                ),
                REFERENCE_REVERSE_REGULATION_RATE_S,
            )
            <= 1.0e-12
            and _relative(
                _finite_float(
                    production["ae4_fully_activated_multiplier"],
                    field="reference.ae4_fully_activated_multiplier",
                ),
                REFERENCE_FULLY_ACTIVATED_MULTIPLIER,
            )
            <= 1.0e-12
        )
        regulation = independent_r3_regulation(
            tau_pka_s=REFERENCE_TAU_PKA_S,
            forward_regulation_rate_s=REFERENCE_FORWARD_REGULATION_RATE_S,
            reverse_regulation_rate_s=REFERENCE_REVERSE_REGULATION_RATE_S,
            basal_multiplier=1.0,
            fully_activated_increment=REFERENCE_FULLY_ACTIVATED_MULTIPLIER - 1.0,
        )
        beta = 0.0 if production["arm"] == "CCH_ONLY" else 1.0
        protocol = IndependentSecretagogueProtocol(
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=EFFECTIVE_CALCIUM_UM,
            beta_occupancy_on=beta,
        )
        nkcc = IndependentN1NkccRegulation(
            fully_activated_multiplier=NKCC1_FULLY_ACTIVATED_MULTIPLIER,
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=EFFECTIVE_CALCIUM_UM,
        )
        model = definition.build_independent_model(
            regulation=regulation,
            nkcc_regulation=nkcc,
            stimulus=protocol,
        )
        initial = np.r_[
            np.asarray(json.loads(root["core_state_json"]), dtype=float),
            np.asarray(regulation.basal_state, dtype=float),
        ]
        atol = np.asarray(
            [
                1.0e-10
                if name.endswith("_fmol")
                else 1.0e-12
                if name.endswith("_pL")
                else 1.0e-10
                for name in model.state_names
            ],
            dtype=float,
        )
        method = METHOD_BY_SOLVER[production["solver_label"]]
        trajectory = model.integrate(
            initial,
            independent_physical_time_grid(600.0, 5.0),
            method=method,
            rtol=1.0e-7,
            atol=atol,
            max_step_s=2.0,
            preserve_basal_zero=True,
        )
        key = (root["root_id"], production["arm"], method)
        if key in trajectories:
            raise AssertionError(f"duplicate independent trajectory {key}")
        trajectories[key] = trajectory
        endpoint_na, endpoint_k, endpoint_volume = _endpoint_concentrations(
            trajectory.states[:12, -1]
        )
        minute_flow = np.interp(
            MINUTE_TIMES_S, trajectory.time_s, trajectory.flow_pL_s
        )
        comparisons = {
            "endpoint_na": _relative(
                endpoint_na, float(production["endpoint_na_mM"])
            ),
            "endpoint_k": _relative(
                endpoint_k, float(production["endpoint_k_mM"])
            ),
            "endpoint_cl": _relative(
                float(trajectory.cell_cl_mM[-1]),
                float(production["endpoint_cl_mM"]),
            ),
            "endpoint_ph": _relative(
                float(trajectory.cell_ph[-1]), float(production["endpoint_ph"])
            ),
            "endpoint_volume": _relative(
                endpoint_volume, float(production["endpoint_volume_pL"])
            ),
            "endpoint_flow": _relative(
                float(trajectory.flow_pL_s[-1]),
                float(production["endpoint_flow_pL_s"]),
            ),
            "cumulative_flow": _relative(
                float(trajectory.cumulative_flow_pL[-1]),
                float(production["cumulative_flow_600_pL"]),
            ),
            "capacity": _relative(
                float(trajectory.capacity_multiplier[-1]),
                float(production["capacity_multiplier_endpoint"]),
            ),
        }
        minute_differences = [
            _relative(value, float(production[f"flow_minute_{minute}_pL_s"]))
            for minute, value in enumerate(minute_flow, start=1)
        ]
        maximum_difference = max(*comparisons.values(), *minute_differences)
        production_success = _csv_bool(
            production.get("success"), field="reference.production.success"
        )
        independent_all_flow_nonnegative = bool(
            np.all(trajectory.flow_pL_s >= 0.0)
        )
        records.append(
            {
                "root_id": root["root_id"],
                "source_class": root["source_class"],
                "source_scale": float(root["source_scale"]),
                "hydraulic_mode": root["hydraulic_mode"],
                "arm": production["arm"],
                "production_solver_label": production["solver_label"],
                "independent_method": method,
                "independent_success": trajectory.success,
                "independent_positive_core": bool(
                    np.all(trajectory.states[:12] > 0.0)
                ),
                "independent_all_flow_nonnegative": (
                    independent_all_flow_nonnegative
                ),
                "independent_reference_declaration_pass": declaration_pass,
                "independent_max_conservation_fmol_s": (
                    trajectory.max_abs_conservation_fmol_s
                ),
                "independent_max_water_residual_pL_s": (
                    trajectory.max_abs_water_accounting_pL_s
                ),
                "independent_endpoint_na_mM": endpoint_na,
                "independent_endpoint_k_mM": endpoint_k,
                "independent_endpoint_cl_mM": float(trajectory.cell_cl_mM[-1]),
                "independent_endpoint_ph": float(trajectory.cell_ph[-1]),
                "independent_endpoint_volume_pL": endpoint_volume,
                "independent_endpoint_flow_pL_s": float(
                    trajectory.flow_pL_s[-1]
                ),
                "independent_cumulative_flow_600_pL": float(
                    trajectory.cumulative_flow_pL[-1]
                ),
                "independent_sustainment_ratio_q600_q60": float(
                    trajectory.flow_pL_s[-1]
                    / np.interp(60.0, trajectory.time_s, trajectory.flow_pL_s)
                ),
                "maximum_relative_production_difference": maximum_difference,
                "production_reproduction_pass": bool(
                    production_success
                    and trajectory.success
                    and np.all(trajectory.states[:12] > 0.0)
                    and independent_all_flow_nonnegative
                    and declaration_pass
                    and maximum_difference <= 1.0e-5
                    and trajectory.max_abs_conservation_fmol_s <= 1.0e-8
                    and trajectory.max_abs_water_accounting_pL_s <= 1.0e-10
                ),
            }
        )

    solver_differences: list[float] = []
    for root_id in intake.root_ids:
        for arm in ARMS:
            radau = trajectories[(root_id, arm, "Radau")]
            bdf = trajectories[(root_id, arm, "BDF")]
            denominator = np.maximum(
                np.maximum(np.abs(radau.states), np.abs(bdf.states)),
                1.0e-10,
            )
            solver_differences.append(
                float(np.max(np.abs(radau.states - bdf.states) / denominator))
            )

    group_records, independent_passing_ids = _recompute_group_classifications(
        intake
    )
    independent_solver_gate_pass = bool(
        solver_differences
        and max(solver_differences) <= SOLVER_RELATIVE_TOLERANCE
    )
    equation_reproduction_pass = bool(
        records and all(row["production_reproduction_pass"] for row in records)
        and independent_solver_gate_pass
    )
    group_reproduction_pass = bool(
        group_records
        and all(row["all_group_classification_agreement"] for row in group_records)
    )
    passing_id_agreement = (
        independent_passing_ids == intake.production_passing_root_ids
    )
    overall_pass = bool(
        equation_reproduction_pass
        and group_reproduction_pass
        and passing_id_agreement
    )
    summary = {
        "status": (
            "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_PASS"
            if overall_pass
            else "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_DISCREPANCY"
        ),
        "firewall": "PRE_REVEAL_WT_ONLY_NO_GENOTYPE_NO_HELDOUT_ACCESS",
        "root_count": len(intake.root_ids),
        "root_ids": list(intake.root_ids),
        "case_count": len(records),
        "expected_case_count": 4 * len(intake.root_ids),
        "all_cases_reproduce_production": equation_reproduction_pass,
        "independent_reference_radau_bdf_gate_pass": (
            independent_solver_gate_pass
        ),
        "all_group_classifications_reproduce_production": group_reproduction_pass,
        "production_scientific_passing_root_ids": list(
            intake.production_passing_root_ids
        ),
        "independent_scientific_passing_root_ids": list(independent_passing_ids),
        "scientific_passing_root_id_agreement": passing_id_agreement,
        "maximum_relative_production_difference": (
            max(row["maximum_relative_production_difference"] for row in records)
            if records
            else None
        ),
        "maximum_independent_radau_bdf_state_difference": (
            max(solver_differences) if solver_differences else None
        ),
        "clean_auditor_authorization_required": True,
        "final_composite_authorization_complete": False,
        "phenotype_reveal_permitted": False,
    }
    return records, group_records, summary


def _csv_bytes(
    rows: Sequence[Mapping[str, Any]], *, empty_fields: Sequence[str]
) -> bytes:
    fields = list(rows[0]) if rows else list(empty_fields)
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8")


def _atomic_write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _require_intake_unchanged(intake: FrozenIndependentDynamicIntake) -> None:
    """Fail before publication if any frozen input or loaded source drifted."""

    _require_unchanged_hash_registry(
        intake.input_artifact_sha256,
        FROZEN_RESULT_INPUTS,
        base=intake.results_directory,
        label="frozen_input_artifact_sha256",
    )
    _require_unchanged_hash_registry(
        intake.production_equation_source_sha256,
        PRODUCTION_EQUATION_SOURCES,
        base=intake.repository,
        label="production_equation_source_sha256",
    )
    _require_unchanged_hash_registry(
        intake.independent_source_file_sha256,
        INDEPENDENT_DYNAMIC_SOURCES,
        base=intake.repository,
        label="independent_source_file_sha256",
    )


def _invalidate_completion_marker(path: Path) -> None:
    """Remove an old summary before replacing any member of its output set."""

    try:
        path.unlink()
    except FileNotFoundError:
        return
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def run_independent_native_dynamic(repository: str | Path) -> Mapping[str, Any]:
    intake = load_frozen_independent_dynamic_intake(repository)
    rows, group_rows, scientific_summary = reproduce_reference_dynamics(intake)
    output_path = intake.results_directory / OUTPUT_FILENAME
    group_output_path = intake.results_directory / GROUP_OUTPUT_FILENAME
    output_hash_path = intake.results_directory / OUTPUT_HASH_FILENAME
    summary_path = intake.results_directory / SUMMARY_FILENAME

    # Publication is a completion-marker protocol.  Recheck every source and
    # frozen input after the long solve, then invalidate an older completion
    # marker before replacing any output it might otherwise authenticate.
    _require_intake_unchanged(intake)
    _invalidate_completion_marker(summary_path)
    _atomic_write_bytes(
        output_path,
        _csv_bytes(
            rows, empty_fields=("root_id", "production_reproduction_pass")
        ),
    )
    _atomic_write_bytes(
        group_output_path,
        _csv_bytes(
            group_rows,
            empty_fields=("root_id", "all_group_classification_agreement"),
        ),
    )
    output_hashes = {
        filename: _sha256(intake.results_directory / filename)
        for filename in OUTPUT_HASH_ARTIFACTS
    }
    ledger_rows = [
        {"artifact": filename, "sha256": output_hashes[filename]}
        for filename in OUTPUT_HASH_ARTIFACTS
    ]
    _atomic_write_bytes(
        output_hash_path,
        _csv_bytes(ledger_rows, empty_fields=("artifact", "sha256")),
    )

    summary = {
        **scientific_summary,
        "frozen_input_artifact_sha256": dict(intake.input_artifact_sha256),
        "production_equation_source_sha256": dict(
            intake.production_equation_source_sha256
        ),
        "production_equation_source_tree_sha256": (
            intake.production_equation_source_tree_sha256
        ),
        "opaque_firewall_ledger_sha256": dict(
            intake.opaque_firewall_ledger_sha256
        ),
        "independent_source_file_sha256": dict(
            intake.independent_source_file_sha256
        ),
        "independent_source_tree_sha256": intake.independent_source_tree_sha256,
        "independent_output_artifact_sha256": output_hashes,
        "independent_output_hash_ledger": OUTPUT_HASH_FILENAME,
        "independent_output_hash_ledger_sha256": _sha256(output_hash_path),
        "publication_order": (
            "CASE_CSV_THEN_GROUP_CSV_THEN_HASH_LEDGER_THEN_SUMMARY_ATOMIC_REPLACE"
        ),
    }
    _atomic_write_bytes(
        summary_path,
        (
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8"),
    )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[2]
    )
    args = parser.parse_args(argv)
    summary = run_independent_native_dynamic(args.repository)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return (
        0
        if summary.get("status") == "INDEPENDENT_NATIVE_REFERENCE_DYNAMICS_PASS"
        else 2
    )


if __name__ == "__main__":
    raise SystemExit(main())
