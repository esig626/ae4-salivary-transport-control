"""Run and freeze the native-source WT dynamic contract before reveal.

The intake is ``native_source_wt_roots.csv`` through the public loader in
``native_source_panel``.  Every retained WT resting root is propagated through
the same fixed CA010/N1-M1.75 stimulus profile, both CCh-only and CCh+IPR
protocol arms, all 20 predeclared R0--R3 regulatory sensitivities, and both
production stiff solvers.  The reference R3/G1.25/t-reference member alone is
also checked at loose and tight tolerances.  No source scalar, regulatory
constant, observation multiplier, or protocol input is optimized here.

The resulting scientific gate never authorizes phenotype reveal.  It records
whether a clean adversarial auditor has enough frozen WT evidence to make that
separate decision.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .model import ModernFullModel
from .native_dynamic_contract import (
    assess_geometry_scale,
    geometry_contract_payload,
    paired_arm_total_ratio,
    regulatory_member_pairs,
    relative_difference,
    shared_observation_scale_interval,
    sustainment_gate,
    sustainment_ratio,
    validate_regulatory_ensemble,
)
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    LOOSE_RADAU,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    TIGHT_RADAU,
    SecretagogueProtocol,
    SolverSpecification,
    StimulusArm,
    Trajectory,
    assert_pre_reveal_payload,
    nearby_charge_preserving_states,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    sha256_file,
    sha256_object,
    simulate_wt,
    solver_comparison,
    write_json,
)


PROFILE_FILENAME = "native_dynamic_contract_profile.csv"
GROUP_FILENAME = "native_dynamic_contract_group_gate.csv"
SOLVER_FILENAME = "native_dynamic_contract_solver_crosschecks.csv"
ISOMORPHISM_FILENAME = "native_dynamic_contract_isomorphism.csv"
NEARBY_FILENAME = "native_dynamic_contract_nearby.csv"
MANIFEST_FILENAME = "native_dynamic_contract_manifest.json"
GATE_FILENAME = "native_dynamic_contract_gate.json"
HASH_FILENAME = "native_dynamic_contract_hashes.csv"

WT_DYNAMIC_REPORT_FILENAMES = (
    PROFILE_FILENAME,
    GROUP_FILENAME,
    SOLVER_FILENAME,
    ISOMORPHISM_FILENAME,
    NEARBY_FILENAME,
)
DYNAMIC_FREEZE_ARTIFACT_FILENAMES = (
    *WT_DYNAMIC_REPORT_FILENAMES,
    MANIFEST_FILENAME,
    GATE_FILENAME,
)
NATIVE_FREEZE_INPUT_FILENAMES = (
    "native_source_map.json",
    "native_source_root_attempts.csv",
    "native_source_wt_roots.csv",
    "native_source_wt_summary.json",
    "native_source_panel.csv",
    "native_source_scale_one_identity.csv",
    "native_source_hashes.csv",
    "independent_native_roots.csv",
    "independent_native_attempts.csv",
    "independent_native_summary.json",
)
FROZEN_EQUATION_PATHS = (
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

REFERENCE_MEMBER_ID = "R3_G125_P21_PROSE_TREFERENCE"
EFFECTIVE_CALCIUM_UM = 0.10
NKCC1_FULLY_ACTIVATED_MULTIPLIER = 1.75
SOLVER_RELATIVE_TOLERANCE = 1.0e-4
ISOMORPHISM_RELATIVE_TOLERANCE = 1.0e-10
COSTIM_TO_CCH_TOTAL_RATIO_MAX = 1.10
SUSTAINMENT_RATIO_BOUNDS = (0.8, 1.2)
MINUTE_TIMES_S = tuple(float(60 * minute) for minute in range(1, 11))
ALLOWED_PRODUCTION_ELIGIBILITY = frozenset({"PRODUCTION_PRE_REVEAL"})
EXPECTED_PRODUCTION_CASE_COUNT_PER_ROOT = 20 * 2 * 2
EXPECTED_NEARBY_CASE_COUNT_PER_ROOT = 8 * 2 * 2
EXPECTED_SOLVER_COMPARISON_COUNT_PER_ROOT = 60
EXPECTED_ISOMORPHISM_COUNT_PER_ROOT = 6 * 2 * 2
NATIVE_FINAL_READY_FLAGS = (
    "production_confirmation_complete",
    "geometry_confirmation_complete",
    "independent_root_reproduction_complete",
    "native_pre_reveal_freeze_ready",
)
SEALED_TARGET_FILENAME = "heldout_" + "targets.csv"


@dataclass(frozen=True)
class DynamicRequest:
    root: Any
    member_id: str
    arm: str
    solver_label: str
    initial_condition: str = "native_source_rest_root"
    initial_state: tuple[float, ...] | None = None


@dataclass(frozen=True)
class DynamicResult:
    metadata: Mapping[str, Any]
    trajectory: Trajectory | None


@dataclass(frozen=True)
class FreezeDependencySnapshot:
    """Exact pre-run bytes that the long dynamic calculation must retain."""

    equation_source_file_sha256: Mapping[str, str]
    native_input_artifact_sha256: Mapping[str, str]
    opaque_firewall_ledger_sha256: Mapping[str, str]


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, np.ndarray):
        return value.tolist()
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


def _root_value(root: Any, name: str, default: Any = None) -> Any:
    if hasattr(root, name):
        return getattr(root, name)
    if isinstance(root, Mapping):
        return root.get(name, default)
    return default


def _root_id(root: Any) -> str:
    value = _root_value(root, "root_id")
    if value is None:
        raise ValueError("native root lacks root_id")
    return str(value)


def _core_state(root: Any) -> tuple[float, ...]:
    value = _root_value(root, "core_state")
    if value is None:
        value = _root_value(root, "core_state_json")
    if isinstance(value, str):
        value = json.loads(value)
    if isinstance(value, Mapping):
        # The native loader should normally return an ordered tuple.  Mapping
        # fallback is deliberately rejected because silently sorting conserved
        # state names would be unsafe.
        raise ValueError("native core_state mapping must be decoded by its loader")
    state = tuple(float(item) for item in value)
    if not state or not all(math.isfinite(item) and item > 0.0 for item in state):
        raise ValueError("native root core_state must be finite and positive")
    return state


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
        "passes_numerical_gate",
        "passes_wt_gate",
    )
    return {name: _root_value(root, name) for name in names}


def _strict_bool(value: Any, name: str) -> bool:
    if value is True:
        return True
    if value is False:
        return False
    if isinstance(value, str) and value in {"True", "False"}:
        return value == "True"
    raise ValueError(f"{name} must be an explicit bool or canonical True/False string")


def require_native_final_ready(results_directory: Path) -> Mapping[str, Any]:
    """Require the post-confirmation native freeze marker before dynamics."""

    path = results_directory / "native_source_wt_summary.json"
    if not path.is_file():
        raise FileNotFoundError("native source WT summary is required")
    summary = json.loads(path.read_text(encoding="utf-8"))
    missing = [name for name in NATIVE_FINAL_READY_FLAGS if name not in summary]
    if missing:
        raise RuntimeError(f"native source summary lacks final-ready flags: {missing}")
    failed = [name for name in NATIVE_FINAL_READY_FLAGS if summary[name] is not True]
    if failed:
        raise RuntimeError(f"native source freeze is not final-ready: {failed}")
    return summary


def eligible_native_roots(roots: Iterable[Any]) -> tuple[Any, ...]:
    """Retain only H=1 production roots that already pass both rest gates.

    Diagnostic and rejected roots remain immutable in the native source
    artifacts; propagating them through 80 dynamic cases would neither repair
    their resting physiology nor constitute an admissible production model.
    """

    retained: list[Any] = []
    for root in roots:
        numerical = _strict_bool(
            _root_value(root, "passes_numerical_gate"),
            f"{_root_id(root)}.passes_numerical_gate",
        )
        wt = _strict_bool(
            _root_value(root, "passes_wt_gate"),
            f"{_root_id(root)}.passes_wt_gate",
        )
        eligibility = str(_root_value(root, "eligibility", ""))
        hydraulic_scale = float(_root_value(root, "hydraulic_scale", math.nan))
        hydraulic_mode = str(_root_value(root, "hydraulic_mode", ""))
        if (
            numerical
            and wt
            and eligibility in ALLOWED_PRODUCTION_ELIGIBILITY
            and hydraulic_scale == 1.0
            and hydraulic_mode == "H1"
        ):
            retained.append(root)
    retained_tuple = tuple(retained)
    ids = tuple(_root_id(root) for root in retained_tuple)
    if len(set(ids)) != len(ids):
        raise ValueError("eligible native root IDs must be unique")
    return retained_tuple


def _member_registry() -> Mapping[str, Any]:
    return {member.member_id: member for member in pre_reveal_regulatory_ensemble()}


def _solver_registry() -> Mapping[str, SolverSpecification]:
    return {
        solver.label: solver
        for solver in (PRODUCTION_RADAU, PRODUCTION_BDF, LOOSE_RADAU, TIGHT_RADAU)
    }


def _build_model(root: Any, member: Any, arm: StimulusArm) -> ModernFullModel:
    # Imported lazily so the pure dynamic-contract helpers remain testable
    # while a native root panel is being generated in parallel.
    from .native_source_panel import build_native_source_model

    protocol = SecretagogueProtocol(
        arm=arm,
        stimulated_calcium_uM=EFFECTIVE_CALCIUM_UM,
    )
    template = build_native_source_model(
        root,
        regulatory_model=member.regulatory_model,
        stimulus=protocol,
    )
    return attach_stimulated_nkcc1(
        template,
        N1AlgebraicNkcc1(
            fully_activated_multiplier=NKCC1_FULLY_ACTIVATED_MULTIPLIER,
            resting_calcium_uM=protocol.resting_calcium_uM,
            stimulated_calcium_uM=protocol.stimulated_calcium_uM,
        ),
    )


def _failed_result(request: DynamicRequest, message: str) -> DynamicResult:
    member = _member_registry()[request.member_id]
    return DynamicResult(
        metadata={
            **_root_metadata(request.root),
            "regulatory_member_id": member.member_id,
            "regulatory_family": member.family,
            "regulatory_gain_label": member.gain_label,
            "regulatory_timing_label": member.timing_label,
            "arm": request.arm,
            "solver_label": request.solver_label,
            "initial_condition": request.initial_condition,
            "success": False,
            "message": message,
            "positive_core": False,
            "all_flow_nonnegative": False,
            "numerical_gate_pass": False,
            "sustainment_ratio_q600_q60": math.nan,
            "sustainment_gate_pass": False,
        },
        trajectory=None,
    )


def _simulate_request(request: DynamicRequest) -> DynamicResult:
    members = _member_registry()
    solvers = _solver_registry()
    member = members[request.member_id]
    solver = solvers[request.solver_label]
    try:
        arm = StimulusArm(request.arm)
        model = _build_model(request.root, member, arm)
        trajectory = simulate_wt(
            model,
            (
                _core_state(request.root)
                if request.initial_state is None
                else request.initial_state
            ),
            family=member.member_id,
            initial_condition=request.initial_condition,
            solver=solver,
            time_s=physical_time_grid(duration_s=600.0, sample_step_s=5.0),
        )
        if not trajectory.success:
            return _failed_result(request, trajectory.message)
        conservation_ratios = {
            name: trajectory.max_abs_conservation_residuals[name] / tolerance
            for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
        }
        maximum_conservation_ratio = max(conservation_ratios.values())
        minute_flow = np.interp(
            np.asarray(MINUTE_TIMES_S), trajectory.time_s, trajectory.flow_pL_s
        )
        ratio = sustainment_ratio(trajectory.time_s, trajectory.flow_pL_s)
        metadata: dict[str, Any] = {
            **_root_metadata(request.root),
            "regulatory_member_id": member.member_id,
            "regulatory_family": member.family,
            "regulatory_gain_label": member.gain_label,
            "ae4_fully_activated_multiplier": member.fully_activated_multiplier,
            "regulatory_timing_label": member.timing_label,
            "regulatory_timing_parameters_json": json.dumps(
                dict(member.timing_parameters), sort_keys=True
            ),
            "regulatory_kinetic_evidence_status": member.kinetic_evidence_status,
            "arm": request.arm,
            "solver_label": request.solver_label,
            "initial_condition": request.initial_condition,
            "solver_method": solver.method,
            "success": trajectory.success,
            "message": trajectory.message,
            "positive_core": trajectory.positive_core,
            "all_flow_nonnegative": bool(np.all(trajectory.flow_pL_s >= 0.0)),
            "max_dimensionless_conservation_ratio": maximum_conservation_ratio,
            "numerical_gate_pass": bool(
                trajectory.success
                and trajectory.positive_core
                and np.all(trajectory.flow_pL_s >= 0.0)
                and maximum_conservation_ratio <= 1.0
            ),
            "sustainment_ratio_q600_q60": ratio,
            "sustainment_gate_pass": sustainment_gate(ratio),
            "capacity_multiplier_start": float(trajectory.capacity_multiplier[0]),
            "capacity_multiplier_endpoint": float(trajectory.capacity_multiplier[-1]),
            "endpoint_na_mM": float(trajectory.cell_na_mM[-1]),
            "endpoint_k_mM": float(trajectory.cell_k_mM[-1]),
            "endpoint_cl_mM": float(trajectory.cell_cl_mM[-1]),
            "endpoint_ph": float(trajectory.cell_ph[-1]),
            "endpoint_volume_pL": float(trajectory.cell_volume_pL[-1]),
            "endpoint_flow_pL_s": float(trajectory.flow_pL_s[-1]),
            "cumulative_flow_600_pL": float(trajectory.cumulative_flow_pL[-1]),
        }
        for minute, value in zip(range(1, 11), minute_flow):
            metadata[f"flow_minute_{minute}_pL_s"] = float(value)
        return DynamicResult(metadata=metadata, trajectory=trajectory)
    except Exception as exc:
        return _failed_result(request, f"{type(exc).__name__}: {exc}")


def _run_requests(
    requests: Sequence[DynamicRequest], workers: int
) -> list[DynamicResult]:
    if workers <= 1:
        return [_simulate_request(request) for request in requests]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_simulate_request, requests, chunksize=1))


def _write_rows(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    empty_fields: Sequence[str] = (),
) -> None:
    fields: list[str] = []
    for row in rows:
        for name in row:
            if name not in fields:
                fields.append(name)
    if not fields:
        fields = list(empty_fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _result_lookup(results: Iterable[DynamicResult]) -> Mapping[tuple[str, str, str, str], DynamicResult]:
    lookup: dict[tuple[str, str, str, str], DynamicResult] = {}
    for result in results:
        row = result.metadata
        key = (
            str(row["root_id"]),
            str(row["regulatory_member_id"]),
            str(row["arm"]),
            str(row["solver_label"]),
        )
        if key in lookup:
            raise AssertionError(f"duplicate dynamic case {key}")
        lookup[key] = result
    return lookup


def _comparison_row(
    reference: DynamicResult,
    comparison: DynamicResult,
    *,
    comparison_kind: str,
) -> dict[str, Any]:
    base = {
        **_root_metadata_from_row(reference.metadata),
        "regulatory_member_id": reference.metadata["regulatory_member_id"],
        "arm": reference.metadata["arm"],
        "initial_condition": reference.metadata["initial_condition"],
        "comparison_kind": comparison_kind,
        "reference_solver_label": reference.metadata["solver_label"],
        "comparison_solver_label": comparison.metadata["solver_label"],
        "reference_numerical_gate_pass": bool(
            reference.metadata.get("numerical_gate_pass", False)
        ),
        "comparison_numerical_gate_pass": bool(
            comparison.metadata.get("numerical_gate_pass", False)
        ),
    }
    if reference.trajectory is None or comparison.trajectory is None:
        return {
            **base,
            "both_success": False,
            "max_relative_state_difference": math.nan,
            "endpoint_cumulative_flow_relative_difference": math.nan,
            "solver_gate_pass": False,
        }
    metrics = dict(solver_comparison(reference.trajectory, comparison.trajectory))
    gate = bool(
        metrics["both_success"]
        and base["reference_numerical_gate_pass"]
        and base["comparison_numerical_gate_pass"]
        and float(metrics["max_relative_state_difference"])
        <= SOLVER_RELATIVE_TOLERANCE
        and float(metrics["max_relative_flow_difference"])
        <= SOLVER_RELATIVE_TOLERANCE
        and float(metrics["endpoint_cumulative_flow_relative_difference"])
        <= SOLVER_RELATIVE_TOLERANCE
    )
    return {**base, **metrics, "solver_gate_pass": gate}


def _root_metadata_from_row(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        name: row.get(name)
        for name in (
            "root_id",
            "panel_id",
            "source_class",
            "source_scale",
            "pump_capacity_scale",
            "apical_pump_fraction",
            "apical_k_fraction",
            "ae4_cation_fraction",
            "hydraulic_scale",
            "hydraulic_mode",
        )
    }


def _isomorphism_row(
    r2: DynamicResult,
    r3: DynamicResult,
) -> dict[str, Any]:
    base = {
        **_root_metadata_from_row(r2.metadata),
        "r2_member_id": r2.metadata["regulatory_member_id"],
        "r3_member_id": r3.metadata["regulatory_member_id"],
        "arm": r2.metadata["arm"],
        "solver_label": r2.metadata["solver_label"],
    }
    if r2.trajectory is None or r3.trajectory is None:
        return {
            **base,
            "max_relative_state_difference": math.nan,
            "max_relative_flow_difference": math.nan,
            "max_relative_capacity_difference": math.nan,
            "r2_r3_isomorphism_gate_pass": False,
        }
    state = relative_difference(r2.trajectory.states, r3.trajectory.states)
    flow = relative_difference(r2.trajectory.flow_pL_s, r3.trajectory.flow_pL_s)
    capacity = relative_difference(
        r2.trajectory.capacity_multiplier,
        r3.trajectory.capacity_multiplier,
    )
    return {
        **base,
        "max_relative_state_difference": state,
        "max_relative_flow_difference": flow,
        "max_relative_capacity_difference": capacity,
        "r2_r3_isomorphism_gate_pass": bool(
            max(state, flow, capacity) <= ISOMORPHISM_RELATIVE_TOLERANCE
        ),
    }


def _has_exact_key_set(
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    expected: set[tuple[str, ...]],
) -> bool:
    """Fail closed unless rows contain every frozen key exactly once."""

    try:
        keys = {
            tuple(str(row[field]) for field in fields)
            for row in rows
        }
    except (KeyError, TypeError):
        return False
    return bool(len(rows) == len(expected) and keys == expected)


def _group_rows(
    roots: Sequence[Any],
    production_results: Sequence[DynamicResult],
    nearby_results: Sequence[DynamicResult],
    solver_rows: Sequence[Mapping[str, Any]],
    isomorphism_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    members = validate_regulatory_ensemble(pre_reveal_regulatory_ensemble())
    arms = (StimulusArm.CCH_ONLY, StimulusArm.CCH_IPR)
    production_solvers = (PRODUCTION_RADAU, PRODUCTION_BDF)
    nearby_conditions_expected = frozenset(
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
    expected_solver_keys = {
        (
            "RADAU_VS_BDF",
            member.member_id,
            arm.value,
            "native_source_rest_root",
            PRODUCTION_RADAU.label,
            PRODUCTION_BDF.label,
        )
        for member in members
        for arm in arms
    }
    expected_solver_keys.update(
        {
            (
                "NEARBY_RADAU_VS_BDF",
                REFERENCE_MEMBER_ID,
                arm.value,
                condition,
                PRODUCTION_RADAU.label,
                PRODUCTION_BDF.label,
            )
            for condition in nearby_conditions_expected
            for arm in arms
        }
    )
    expected_solver_keys.update(
        {
            (
                f"PRODUCTION_VS_{solver.label.upper()}",
                REFERENCE_MEMBER_ID,
                arm.value,
                "native_source_rest_root",
                PRODUCTION_RADAU.label,
                solver.label,
            )
            for arm in arms
            for solver in (LOOSE_RADAU, TIGHT_RADAU)
        }
    )
    expected_isomorphism_keys = {
        (r2_id, r3_id, arm.value, solver.label)
        for r2_id, r3_id in regulatory_member_pairs(members)
        for arm in arms
        for solver in production_solvers
    }
    if len(expected_solver_keys) != EXPECTED_SOLVER_COMPARISON_COUNT_PER_ROOT:
        raise AssertionError("frozen solver-comparison key contract changed")
    if len(expected_isomorphism_keys) != EXPECTED_ISOMORPHISM_COUNT_PER_ROOT:
        raise AssertionError("frozen R2/R3 isomorphism key contract changed")

    by_root: dict[str, list[DynamicResult]] = {}
    for result in production_results:
        by_root.setdefault(str(result.metadata["root_id"]), []).append(result)
    nearby_by_root: dict[str, list[DynamicResult]] = {}
    for result in nearby_results:
        nearby_by_root.setdefault(str(result.metadata["root_id"]), []).append(result)
    solver_by_root: dict[str, list[Mapping[str, Any]]] = {}
    for row in solver_rows:
        solver_by_root.setdefault(str(row["root_id"]), []).append(row)
    iso_by_root: dict[str, list[Mapping[str, Any]]] = {}
    for row in isomorphism_rows:
        iso_by_root.setdefault(str(row["root_id"]), []).append(row)

    rows: list[dict[str, Any]] = []
    for root in roots:
        root_id = _root_id(root)
        cases = by_root.get(root_id, [])
        nearby_cases = nearby_by_root.get(root_id, [])
        expected_count = EXPECTED_PRODUCTION_CASE_COUNT_PER_ROOT
        successful = [case for case in cases if case.trajectory is not None]
        successful_nearby = [
            case for case in nearby_cases if case.trajectory is not None
        ]
        minute_flows = [
            float(case.metadata[f"flow_minute_{minute}_pL_s"])
            for case in (*successful, *successful_nearby)
            for minute in range(1, 11)
        ]
        if minute_flows:
            interval = shared_observation_scale_interval(minute_flows)
            one_smg = assess_geometry_scale(interval)
            paired_smg = assess_geometry_scale(
                interval, paired_smg_sensitivity=True
            )
            interval_payload = asdict(interval)
            one_payload = asdict(one_smg)
            paired_payload = asdict(paired_smg)
        else:
            interval_payload = {
                "lower_uL_min_per_pL_s": math.nan,
                "upper_uL_min_per_pL_s": math.nan,
                "sample_count": 0,
                "minimum_flow_pL_s": math.nan,
                "maximum_flow_pL_s": math.nan,
                "maximum_to_minimum_flow_ratio": math.nan,
            }
            one_payload = {
                "geometry_gate_pass": False,
                "feasible_lower_uL_min_per_pL_s": math.nan,
                "feasible_upper_uL_min_per_pL_s": math.nan,
                "implied_cell_count_at_feasible_lower": None,
            }
            paired_payload = {"geometry_gate_pass": False}

        lookup = {
            (
                str(case.metadata["regulatory_member_id"]),
                str(case.metadata["arm"]),
                str(case.metadata["solver_label"]),
            ): case
            for case in cases
        }
        arm_ratios: list[float] = []
        arm_pair_complete = True
        for member in members:
            for solver in production_solvers:
                cch = lookup.get(
                    (member.member_id, StimulusArm.CCH_ONLY.value, solver.label)
                )
                costim = lookup.get(
                    (member.member_id, StimulusArm.CCH_IPR.value, solver.label)
                )
                if (
                    cch is None
                    or costim is None
                    or cch.trajectory is None
                    or costim.trajectory is None
                ):
                    arm_pair_complete = False
                    continue
                arm_ratios.append(
                    paired_arm_total_ratio(
                        cch.trajectory.cumulative_flow_pL[-1],
                        costim.trajectory.cumulative_flow_pL[-1],
                    )
                )
        nearby_lookup = {
            (
                str(case.metadata["initial_condition"]),
                str(case.metadata["arm"]),
                str(case.metadata["solver_label"]),
            ): case
            for case in nearby_cases
        }
        nearby_conditions = sorted(
            {
                str(case.metadata["initial_condition"])
                for case in nearby_cases
            }
        )
        if set(nearby_conditions) != nearby_conditions_expected:
            arm_pair_complete = False
        for condition in nearby_conditions:
            for solver in production_solvers:
                cch = nearby_lookup.get(
                    (condition, StimulusArm.CCH_ONLY.value, solver.label)
                )
                costim = nearby_lookup.get(
                    (condition, StimulusArm.CCH_IPR.value, solver.label)
                )
                if (
                    cch is None
                    or costim is None
                    or cch.trajectory is None
                    or costim.trajectory is None
                ):
                    arm_pair_complete = False
                    continue
                arm_ratios.append(
                    paired_arm_total_ratio(
                        cch.trajectory.cumulative_flow_pL[-1],
                        costim.trajectory.cumulative_flow_pL[-1],
                    )
                )
        regulatory_contract = bool(
            len(successful) == expected_count
            and all(
                abs(float(case.metadata["capacity_multiplier_start"]) - 1.0)
                <= 1.0e-12
                for case in successful
            )
            and all(
                abs(float(case.metadata["capacity_multiplier_endpoint"]) - 1.0)
                <= 1.0e-10
                for case in successful
                if case.metadata["arm"] == StimulusArm.CCH_ONLY.value
            )
            and all(
                float(case.metadata["capacity_multiplier_endpoint"]) > 1.0
                for case in successful
                if case.metadata["arm"] == StimulusArm.CCH_IPR.value
            )
        )
        source_wt_pass = _strict_bool(
            _root_value(root, "passes_wt_gate"),
            f"{root_id}.passes_wt_gate",
        )
        numerical = bool(
            len(cases) == expected_count
            and all(bool(case.metadata.get("numerical_gate_pass")) for case in cases)
        )
        nearby_numerical = bool(
            len(nearby_cases) == EXPECTED_NEARBY_CASE_COUNT_PER_ROOT
            and all(
                bool(case.metadata.get("numerical_gate_pass"))
                for case in nearby_cases
            )
        )
        sustainment = bool(
            len(cases) == expected_count
            and all(bool(case.metadata.get("sustainment_gate_pass")) for case in cases)
        )
        nearby_sustainment = bool(
            len(nearby_cases) == EXPECTED_NEARBY_CASE_COUNT_PER_ROOT
            and all(
                bool(case.metadata.get("sustainment_gate_pass"))
                for case in nearby_cases
            )
        )
        root_solver_rows = solver_by_root.get(root_id, [])
        solver_gate = bool(
            _has_exact_key_set(
                root_solver_rows,
                (
                    "comparison_kind",
                    "regulatory_member_id",
                    "arm",
                    "initial_condition",
                    "reference_solver_label",
                    "comparison_solver_label",
                ),
                expected_solver_keys,
            )
            and all(bool(row["solver_gate_pass"]) for row in root_solver_rows)
        )
        root_iso_rows = iso_by_root.get(root_id, [])
        iso_gate = bool(
            _has_exact_key_set(
                root_iso_rows,
                ("r2_member_id", "r3_member_id", "arm", "solver_label"),
                expected_isomorphism_keys,
            )
            and all(
                bool(row["r2_r3_isomorphism_gate_pass"])
                for row in root_iso_rows
            )
        )
        beta_total_gate = bool(
            arm_pair_complete
            and len(arm_ratios) == 56
            and max(arm_ratios, default=math.inf) <= COSTIM_TO_CCH_TOTAL_RATIO_MAX
        )
        shared_interval_pass = bool(
            minute_flows
            and interval_payload["lower_uL_min_per_pL_s"]
            <= interval_payload["upper_uL_min_per_pL_s"]
        )
        scientific_gate = bool(
            source_wt_pass
            and numerical
            and nearby_numerical
            and sustainment
            and nearby_sustainment
            and solver_gate
            and iso_gate
            and regulatory_contract
            and beta_total_gate
            and shared_interval_pass
            and one_payload["geometry_gate_pass"]
        )
        rows.append(
            {
                **_root_metadata(root),
                "production_case_count": len(cases),
                "expected_production_case_count": expected_count,
                "nearby_case_count": len(nearby_cases),
                "expected_nearby_case_count": EXPECTED_NEARBY_CASE_COUNT_PER_ROOT,
                "all_native_rest_wt_gate_pass": source_wt_pass,
                "all_dynamic_numerical_gates_pass": numerical,
                "all_nearby_numerical_gates_pass": nearby_numerical,
                "all_sustainment_gates_pass": sustainment,
                "all_nearby_sustainment_gates_pass": nearby_sustainment,
                "solver_comparison_count": len(root_solver_rows),
                "expected_solver_comparison_count": (
                    EXPECTED_SOLVER_COMPARISON_COUNT_PER_ROOT
                ),
                "all_solver_and_tolerance_gates_pass": solver_gate,
                "isomorphism_comparison_count": len(root_iso_rows),
                "expected_isomorphism_comparison_count": (
                    EXPECTED_ISOMORPHISM_COUNT_PER_ROOT
                ),
                "r2_r3_isomorphism_gate_pass": iso_gate,
                "p21_regulatory_direction_contract_pass": regulatory_contract,
                "paired_arm_count": len(arm_ratios),
                "costim_to_cch_total_ratio_min": min(arm_ratios, default=math.nan),
                "costim_to_cch_total_ratio_max": max(arm_ratios, default=math.nan),
                "costim_to_cch_total_ratio_gate_pass": beta_total_gate,
                "shared_scale_sample_count": interval_payload["sample_count"],
                "shared_scale_lower_uL_min_per_pL_s": interval_payload[
                    "lower_uL_min_per_pL_s"
                ],
                "shared_scale_upper_uL_min_per_pL_s": interval_payload[
                    "upper_uL_min_per_pL_s"
                ],
                "shared_scale_interval_nonempty": shared_interval_pass,
                "ensemble_min_flow_pL_s": interval_payload["minimum_flow_pL_s"],
                "ensemble_max_flow_pL_s": interval_payload["maximum_flow_pL_s"],
                "ensemble_max_min_flow_ratio": (
                    interval_payload["maximum_flow_pL_s"]
                    / interval_payload["minimum_flow_pL_s"]
                    if minute_flows
                    else math.nan
                ),
                "one_smg_geometry_gate_pass": one_payload[
                    "geometry_gate_pass"
                ],
                "one_smg_feasible_scale_lower_uL_min_per_pL_s": one_payload[
                    "feasible_lower_uL_min_per_pL_s"
                ],
                "one_smg_feasible_scale_upper_uL_min_per_pL_s": one_payload[
                    "feasible_upper_uL_min_per_pL_s"
                ],
                "one_smg_implied_cell_count_at_lower": one_payload[
                    "implied_cell_count_at_feasible_lower"
                ],
                "paired_smg_sensitivity_gate_pass": paired_payload[
                    "geometry_gate_pass"
                ],
                "native_dynamic_scientific_gate_pass": scientific_gate,
            }
        )
    return rows


def _require_relative_files(
    base: Path, paths: Sequence[str]
) -> tuple[str, ...]:
    missing = [path for path in paths if not (base / path).is_file()]
    if missing:
        raise FileNotFoundError(f"required freeze dependencies are missing: {missing}")
    return tuple(paths)


def _opaque_firewall_hashes(results_directory: Path) -> Mapping[str, str]:
    """Hash staged ledgers as opaque bytes; do not parse target contents."""

    names = {
        "calibration_ledger_sha256": "calibration_targets.csv",
        "validation_ledger_sha256": "validation_targets.csv",
        "sealed_phenotype_target_ledger_sha256": SEALED_TARGET_FILENAME,
        "pre_reveal_log_sha256": "reveal_log.json",
    }
    missing = [filename for filename in names.values() if not (results_directory / filename).is_file()]
    if missing:
        raise FileNotFoundError(f"required staged firewall ledgers are missing: {missing}")
    return {
        label: sha256_file(results_directory / filename)
        for label, filename in names.items()
    }


def _capture_freeze_dependency_snapshot(
    repository_root: Path, results_directory: Path
) -> FreezeDependencySnapshot:
    """Hash every frozen dependency before any long numerical work begins."""

    equation_paths = _require_relative_files(
        repository_root,
        FROZEN_EQUATION_PATHS,
    )
    input_names = _require_relative_files(
        results_directory,
        NATIVE_FREEZE_INPUT_FILENAMES,
    )
    return FreezeDependencySnapshot(
        equation_source_file_sha256={
            path: sha256_file(repository_root / path) for path in equation_paths
        },
        native_input_artifact_sha256={
            name: sha256_file(results_directory / name) for name in input_names
        },
        opaque_firewall_ledger_sha256=dict(
            _opaque_firewall_hashes(results_directory)
        ),
    )


def _assert_freeze_snapshot_current(
    snapshot: FreezeDependencySnapshot,
    repository_root: Path,
    results_directory: Path,
) -> None:
    """Fail if a source, native input, or opaque firewall byte changed in flight."""

    current = _capture_freeze_dependency_snapshot(
        repository_root,
        results_directory,
    )
    if current != snapshot:
        changed: list[str] = []
        for label in (
            "equation_source_file_sha256",
            "native_input_artifact_sha256",
            "opaque_firewall_ledger_sha256",
        ):
            before = getattr(snapshot, label)
            after = getattr(current, label)
            keys = set(before) | set(after)
            changed.extend(
                f"{label}:{key}"
                for key in sorted(keys)
                if before.get(key) != after.get(key)
            )
        raise RuntimeError(
            "native dynamic freeze dependencies changed during execution: "
            f"{changed}"
        )


def _assert_native_pre_reveal_payload(payload: Mapping[str, Any]) -> None:
    """Enforce an opaque-hash-only firewall and separate auditor authority."""

    assert_pre_reveal_payload(payload)
    encoded = json.dumps(_jsonable(payload), sort_keys=True).lower()
    forbidden = (
        "ae4_null",
        "knockout_ratio",
        "target_value",
        "target_time_course",
        "heldout_",
    )
    found = [token for token in forbidden if token in encoded]
    if found:
        raise ValueError(f"native pre-reveal payload leaks forbidden schema: {found}")
    firewall_hashes = payload.get("opaque_firewall_ledger_sha256")
    if not isinstance(firewall_hashes, Mapping) or set(firewall_hashes) != {
        "calibration_ledger_sha256",
        "validation_ledger_sha256",
        "sealed_phenotype_target_ledger_sha256",
        "pre_reveal_log_sha256",
    }:
        raise ValueError("native freeze must bind exactly four staged firewall hashes")
    if not all(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
        for value in firewall_hashes.values()
    ):
        raise ValueError("staged firewall bindings must be lowercase SHA-256 digests")
    if payload.get("clean_auditor_authorization_required") is not True:
        raise ValueError("clean auditor authority must remain required")
    if payload.get("phenotype_reveal_permitted_by_this_manifest") is not False:
        raise ValueError("the WT dynamic manifest cannot authorize phenotype reveal")


def _assert_native_gate_payload(payload: Mapping[str, Any]) -> None:
    assert_pre_reveal_payload(payload)
    encoded = json.dumps(_jsonable(payload), sort_keys=True).lower()
    if any(token in encoded for token in ("ae4_null", "knockout_ratio", "target_value")):
        raise ValueError("native WT gate contains a forbidden phenotype schema")
    if payload.get("clean_auditor_authorization_required") is not True:
        raise ValueError("native WT gate must require clean auditor authorization")
    if payload.get("independent_numerical_reproduction_required") is not True:
        raise ValueError("native WT gate must require independent reproduction")
    if payload.get("phenotype_reveal_permitted") is not False:
        raise ValueError("native WT gate cannot permit phenotype reveal")
    if payload.get("final_composite_authorization_complete") is not False:
        raise ValueError("native WT runner cannot complete composite authorization")


def _freeze_manifest(
    repository_root: Path,
    results_directory: Path,
    roots: Sequence[Any],
    group_rows: Sequence[Mapping[str, Any]],
    *,
    report_directory: Path,
    dependency_snapshot: FreezeDependencySnapshot,
    native_ready_summary: Mapping[str, Any],
) -> Mapping[str, Any]:
    members = validate_regulatory_ensemble(pre_reveal_regulatory_ensemble())
    equation_hashes = dict(dependency_snapshot.equation_source_file_sha256)
    input_hashes = dict(dependency_snapshot.native_input_artifact_sha256)
    firewall_hashes = dict(dependency_snapshot.opaque_firewall_ledger_sha256)
    from .native_source_panel import build_native_source_model

    reference_member = next(
        member for member in members if member.member_id == REFERENCE_MEMBER_ID
    )
    root_payloads: dict[str, Mapping[str, Any]] = {}
    for root in roots:
        frozen_model = build_native_source_model(
            root,
            regulatory_model=reference_member.regulatory_model,
            stimulus=SecretagogueProtocol(
                arm=StimulusArm.CCH_IPR,
                stimulated_calcium_uM=EFFECTIVE_CALCIUM_UM,
            ),
        )
        parameter_payload = _jsonable(frozen_model.parameters)
        ae4_payload = _jsonable(frozen_model.ae4_parameters)
        parameter_hash = sha256_object(parameter_payload)
        ae4_hash = sha256_object(ae4_payload)
        invariant_cases: dict[str, Mapping[str, str]] = {}
        for member in members:
            for arm in (StimulusArm.CCH_ONLY, StimulusArm.CCH_IPR):
                comparison_model = build_native_source_model(
                    root,
                    regulatory_model=member.regulatory_model,
                    stimulus=SecretagogueProtocol(
                        arm=arm,
                        stimulated_calcium_uM=EFFECTIVE_CALCIUM_UM,
                    ),
                )
                comparison_parameter_hash = sha256_object(
                    _jsonable(comparison_model.parameters)
                )
                comparison_ae4_hash = sha256_object(
                    _jsonable(comparison_model.ae4_parameters)
                )
                if comparison_parameter_hash != parameter_hash:
                    raise AssertionError(
                        f"{_root_id(root)}: {member.member_id}/{arm.value} "
                        "changes whole-cell parameters"
                    )
                if comparison_ae4_hash != ae4_hash:
                    raise AssertionError(
                        f"{_root_id(root)}: {member.member_id}/{arm.value} "
                        "changes AE4 transport parameters"
                    )
                invariant_cases[f"{member.member_id}::{arm.value}"] = {
                    "whole_cell_parameters_sha256": comparison_parameter_hash,
                    "ae4_parameters_sha256": comparison_ae4_hash,
                }
        root_payloads[_root_id(root)] = {
            "metadata": _root_metadata(root),
            "core_state": list(_core_state(root)),
            "core_state_sha256": sha256_object(_core_state(root)),
            "native_root_object": _jsonable(root),
            "native_root_object_sha256": sha256_object(_jsonable(root)),
            "whole_cell_parameters": parameter_payload,
            "whole_cell_parameters_sha256": parameter_hash,
            "ae4_parameters": ae4_payload,
            "ae4_parameters_sha256": ae4_hash,
            "member_arm_parameter_invariance": invariant_cases,
        }
    regulatory_payloads = {
        member.member_id: {
            **_jsonable(member.metadata()),
            "state_names": list(member.regulatory_model.state_names),
            "model": _jsonable(member.regulatory_model),
        }
        for member in members
    }
    wt_report_hashes = {
        name: sha256_file(report_directory / name)
        for name in WT_DYNAMIC_REPORT_FILENAMES
    }
    scientific_passing_rows = [
        row
        for row in group_rows
        if row.get("native_dynamic_scientific_gate_pass") is True
    ]
    manifest = {
        "manifest_id": "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1",
        "status": "FROZEN_FOR_CLEAN_ADVERSARIAL_REVIEW",
        "supersedes": {
            "manifest": "g5_final_pre_reveal_manifest.json",
            "gate": "g5_final_pre_reveal_gate.json",
            "reason": "ABSOLUTE_GLAND_GEOMETRY_AND_NATIVE_SOURCE_MAP_SUPERSESSION",
        },
        "firewall": "WT_ONLY_NO_PHENOTYPE_TARGET_ACCESS",
        "root_intake_rule": {
            "retain_all_rows": True,
            "first_passer_selection_prohibited": True,
            "strict_panel_precedes_source_classes": True,
            "source_scalar_one_is_exact_nested_control": True,
            "root_loader": "native_source_panel.load_native_source_roots",
            "required_eligibility": sorted(ALLOWED_PRODUCTION_ELIGIBILITY),
            "required_hydraulic_scale": 1.0,
            "required_hydraulic_mode": "H1",
        },
        "native_final_ready_flags": {
            name: native_ready_summary[name] for name in NATIVE_FINAL_READY_FLAGS
        },
        "root_count": len(roots),
        "roots": root_payloads,
        "scientific_passing_group_count": len(scientific_passing_rows),
        "scientific_passing_root_ids": [
            str(row["root_id"]) for row in scientific_passing_rows
        ],
        "protocols": {
            arm.value: _jsonable(
                SecretagogueProtocol(
                    arm=arm,
                    stimulated_calcium_uM=EFFECTIVE_CALCIUM_UM,
                )
            )
            for arm in (StimulusArm.CCH_ONLY, StimulusArm.CCH_IPR)
        },
        "stimulated_nkcc1": {
            "family": "N1_ALGEBRAIC_CCH_CAPACITY",
            "effective_calcium_uM": EFFECTIVE_CALCIUM_UM,
            "fully_activated_multiplier": NKCC1_FULLY_ACTIVATED_MULTIPLIER,
            "status": "FROZEN_WT_PRIMARY_PROFILE_NO_REFIT",
        },
        "regulatory_member_count": len(members),
        "regulatory_members": regulatory_payloads,
        "r2_r3_isomorphism_pairs": regulatory_member_pairs(members),
        "p21_contract": {
            "causal_chain": "BETA_TO_ADENYLATE_CYCLASE_CAMP_TO_PKA_DEPENDENT_AE4_RESPONSE",
            "coupling": "COMMON_AE4_CAPACITY_ONLY",
            "production_minimum": "R1_DYNAMIC",
            "r2_r3_status": "UNIDENTIFIED_EXACT_ISOMORPHIC_SENSITIVITY_CLASS",
            "s173_interpretation": "RESPONSE_DEPENDENCE_NOT_DIRECT_PHOSPHORYLATION_PROOF",
            "s273_interpretation": "RESPONSE_RETAINED",
            "acute_kinetic_status": "UNMEASURED_PREDECLARED_NOT_FIT",
        },
        "production_solvers": _jsonable((PRODUCTION_RADAU, PRODUCTION_BDF)),
        "tolerance_solvers_reference_r3_only": _jsonable(
            (LOOSE_RADAU, TIGHT_RADAU)
        ),
        "solver_relative_tolerance": SOLVER_RELATIVE_TOLERANCE,
        "isomorphism_relative_tolerance": ISOMORPHISM_RELATIVE_TOLERANCE,
        "sustainment_ratio_bounds_q600_q60": SUSTAINMENT_RATIO_BOUNDS,
        "costim_to_cch_total_ratio_max": COSTIM_TO_CCH_TOTAL_RATIO_MAX,
        "minute_times_s": MINUTE_TIMES_S,
        "geometry": geometry_contract_payload(),
        "group_gate_payload_sha256": sha256_object(_jsonable(group_rows)),
        "wt_dynamic_report_artifact_sha256": wt_report_hashes,
        "equation_source_file_sha256": equation_hashes,
        "equation_source_tree_sha256": sha256_object(equation_hashes),
        "native_input_artifact_sha256": input_hashes,
        "opaque_firewall_ledger_sha256": firewall_hashes,
        "opaque_firewall_access_mode": "SHA256_BYTES_ONLY_NO_TARGET_PARSE",
        "independent_numerical_reproduction_required_for_composite_authorization": True,
        "clean_auditor_authorization_required": True,
        "phenotype_reveal_permitted_by_this_manifest": False,
    }
    _assert_native_pre_reveal_payload(manifest)
    return manifest


def _write_hash_ledger(
    results_directory: Path,
    filenames: Sequence[str],
) -> None:
    missing = [
        filename
        for filename in filenames
        if not (results_directory / filename).is_file()
    ]
    if missing:
        raise FileNotFoundError(f"cannot hash missing dynamic artifacts: {missing}")
    rows = [
        {
            "artifact": filename,
            "sha256": sha256_file(results_directory / filename),
        }
        for filename in filenames
    ]
    _write_rows(results_directory / HASH_FILENAME, rows)


def _invalidate_hash_ledger(results_directory: Path) -> None:
    """Remove prior publication authority before a new generation can mutate."""

    (results_directory / HASH_FILENAME).unlink(missing_ok=True)


def _publish_staged_freeze(
    staging_directory: Path,
    results_directory: Path,
    repository_root: Path,
    dependency_snapshot: FreezeDependencySnapshot,
) -> None:
    """Publish atomically per file with the complete hash ledger strictly last.

    Multi-file replacement cannot be globally atomic.  The absent ledger is the
    fail-closed transaction marker: downstream consumers have no authority to
    use a partially published generation.
    """

    _invalidate_hash_ledger(results_directory)
    _assert_freeze_snapshot_current(
        dependency_snapshot,
        repository_root,
        results_directory,
    )
    for filename in DYNAMIC_FREEZE_ARTIFACT_FILENAMES:
        os.replace(staging_directory / filename, results_directory / filename)
    _assert_freeze_snapshot_current(
        dependency_snapshot,
        repository_root,
        results_directory,
    )
    os.replace(
        staging_directory / HASH_FILENAME,
        results_directory / HASH_FILENAME,
    )
    try:
        _assert_freeze_snapshot_current(
            dependency_snapshot,
            repository_root,
            results_directory,
        )
    except Exception:
        _invalidate_hash_ledger(results_directory)
        raise


def run_native_dynamic_contract(
    repository_root: Path,
    results_directory: Path,
    *,
    workers: int = 1,
) -> Mapping[str, Any]:
    """Execute the complete WT-only dynamic contract on every native root."""

    from .native_source_panel import load_native_source_roots

    results_directory.mkdir(parents=True, exist_ok=True)
    dependency_snapshot = _capture_freeze_dependency_snapshot(
        repository_root,
        results_directory,
    )
    native_summary = require_native_final_ready(results_directory)
    all_roots = tuple(load_native_source_roots(results_directory))
    _assert_freeze_snapshot_current(
        dependency_snapshot,
        repository_root,
        results_directory,
    )
    roots = eligible_native_roots(all_roots)
    if not roots:
        raise RuntimeError(
            "native dynamic contract requires at least one numerical, WT-passing, "
            "H=1 production root; rejected and diagnostic roots remain in the "
            "native source artifacts"
        )
    summary_ids = tuple(str(value) for value in native_summary.get("production_root_ids", ()))
    root_ids = tuple(_root_id(root) for root in roots)
    if set(summary_ids) != set(root_ids) or len(summary_ids) != len(root_ids):
        raise RuntimeError(
            "final-ready native summary production_root_ids do not exactly match "
            "the eligible post-geometry root table"
        )
    members = validate_regulatory_ensemble(pre_reveal_regulatory_ensemble())
    arms = (StimulusArm.CCH_ONLY, StimulusArm.CCH_IPR)
    production_solvers = (PRODUCTION_RADAU, PRODUCTION_BDF)
    _invalidate_hash_ledger(results_directory)
    production_requests = [
        DynamicRequest(root, member.member_id, arm.value, solver.label)
        for root in roots
        for member in members
        for arm in arms
        for solver in production_solvers
    ]
    production_results = _run_requests(production_requests, workers)
    lookup = _result_lookup(production_results)

    member_registry = {member.member_id: member for member in members}
    reference_member = member_registry[REFERENCE_MEMBER_ID]
    nearby_requests: list[DynamicRequest] = []
    for root in roots:
        for arm in arms:
            template = _build_model(root, reference_member, arm)
            nearby = nearby_charge_preserving_states(
                template, _core_state(root), fraction=1.0e-3
            )
            if set(nearby) != {
                "baseline",
                "cell_volume_down",
                "cell_volume_up",
                "cell_nacl_down",
                "cell_nacl_up",
                "lumen_volume_down",
                "lumen_volume_up",
                "lumen_nacl_down",
                "lumen_nacl_up",
            }:
                raise AssertionError("nearby-state generator changed its frozen panel")
            for label, vector in nearby.items():
                if label == "baseline":
                    continue
                for solver in production_solvers:
                    nearby_requests.append(
                        DynamicRequest(
                            root,
                            REFERENCE_MEMBER_ID,
                            arm.value,
                            solver.label,
                            initial_condition=f"nearby_{label}",
                            initial_state=tuple(float(value) for value in vector),
                        )
                    )
    nearby_results = _run_requests(nearby_requests, workers)
    nearby_lookup = {
        (
            str(result.metadata["root_id"]),
            str(result.metadata["initial_condition"]),
            str(result.metadata["arm"]),
            str(result.metadata["solver_label"]),
        ): result
        for result in nearby_results
    }
    if len(nearby_lookup) != len(nearby_results):
        raise AssertionError("nearby dynamic case keys must be unique")

    tolerance_requests = [
        DynamicRequest(root, REFERENCE_MEMBER_ID, arm.value, solver.label)
        for root in roots
        for arm in arms
        for solver in (LOOSE_RADAU, TIGHT_RADAU)
    ]
    tolerance_results = _run_requests(tolerance_requests, workers)
    tolerance_lookup = _result_lookup(tolerance_results)

    solver_rows: list[dict[str, Any]] = []
    for root in roots:
        root_id = _root_id(root)
        for member in members:
            for arm in arms:
                radau = lookup[(root_id, member.member_id, arm.value, PRODUCTION_RADAU.label)]
                bdf = lookup[(root_id, member.member_id, arm.value, PRODUCTION_BDF.label)]
                solver_rows.append(
                    _comparison_row(radau, bdf, comparison_kind="RADAU_VS_BDF")
                )
        nearby_conditions = sorted(
            {
                str(result.metadata["initial_condition"])
                for result in nearby_results
                if str(result.metadata["root_id"]) == root_id
            }
        )
        if len(nearby_conditions) != 8:
            raise AssertionError("each root requires eight nearby initial states")
        for condition in nearby_conditions:
            for arm in arms:
                radau = nearby_lookup[
                    (root_id, condition, arm.value, PRODUCTION_RADAU.label)
                ]
                bdf = nearby_lookup[
                    (root_id, condition, arm.value, PRODUCTION_BDF.label)
                ]
                solver_rows.append(
                    _comparison_row(
                        radau,
                        bdf,
                        comparison_kind="NEARBY_RADAU_VS_BDF",
                    )
                )
        for arm in arms:
            reference = lookup[
                (root_id, REFERENCE_MEMBER_ID, arm.value, PRODUCTION_RADAU.label)
            ]
            for solver in (LOOSE_RADAU, TIGHT_RADAU):
                comparison = tolerance_lookup[
                    (root_id, REFERENCE_MEMBER_ID, arm.value, solver.label)
                ]
                solver_rows.append(
                    _comparison_row(
                        reference,
                        comparison,
                        comparison_kind=f"PRODUCTION_VS_{solver.label.upper()}",
                    )
                )

    pairs = regulatory_member_pairs(members)
    isomorphism_rows: list[dict[str, Any]] = []
    for root in roots:
        root_id = _root_id(root)
        for r2_id, r3_id in pairs:
            for arm in arms:
                for solver in production_solvers:
                    isomorphism_rows.append(
                        _isomorphism_row(
                            lookup[(root_id, r2_id, arm.value, solver.label)],
                            lookup[(root_id, r3_id, arm.value, solver.label)],
                        )
                    )

    group_rows = _group_rows(
        roots,
        production_results,
        nearby_results,
        solver_rows,
        isomorphism_rows,
    )
    profile_rows = [dict(result.metadata) for result in production_results]
    passing = [row for row in group_rows if row["native_dynamic_scientific_gate_pass"]]
    gate = {
        "gate_id": "TASK13B_NATIVE_SOURCE_DYNAMIC_PRE_REVEAL_V1",
        "status": (
            "SCIENTIFIC_GATE_PASS_AUDITOR_AUTHORIZATION_REQUIRED"
            if passing
            else "SCIENTIFIC_GATE_FAIL"
        ),
        "supersedes_old_g5": True,
        "old_g5_supersession_reason": (
            "ABSOLUTE_GLAND_GEOMETRY_SCALE_FAIL_AND_NATIVE_SOURCE_MAP_REPAIR"
        ),
        "native_root_count": len(roots),
        "scientific_passing_group_count": len(passing),
        "scientific_passing_root_ids": [row["root_id"] for row in passing],
        "all_rows_retained_no_first_passer": True,
        "regulatory_member_count": len(members),
        "dynamic_production_member_count": sum(
            member.family in {"R1", "R2", "R3"} for member in members
        ),
        "arm_count": len(arms),
        "production_solver_count": len(production_solvers),
        "production_case_count": len(production_results),
        "nearby_case_count": len(nearby_results),
        "all_scientific_group_gates_pass": bool(
            group_rows
            and all(row["native_dynamic_scientific_gate_pass"] for row in group_rows)
        ),
        "at_least_one_scientific_group_pass": bool(passing),
        "clean_auditor_authorization_required": True,
        "clean_auditor_authorization_recorded": False,
        "independent_numerical_reproduction_required": True,
        "independent_numerical_reproduction_recorded_by_this_runner": False,
        "final_composite_authorization_complete": False,
        "phenotype_reveal_permitted": False,
        "firewall": "WT_ONLY_NO_PHENOTYPE_TARGET_ACCESS",
        "manifest": MANIFEST_FILENAME,
        "artifacts": (
            PROFILE_FILENAME,
            GROUP_FILENAME,
            SOLVER_FILENAME,
            ISOMORPHISM_FILENAME,
            NEARBY_FILENAME,
        ),
    }
    _assert_native_gate_payload(gate)
    with tempfile.TemporaryDirectory(
        prefix=".native_dynamic_contract.",
        dir=results_directory,
    ) as temporary:
        staging_directory = Path(temporary)
        _write_rows(
            staging_directory / PROFILE_FILENAME,
            profile_rows,
            empty_fields=("root_id", "regulatory_member_id", "arm", "solver_label"),
        )
        _write_rows(
            staging_directory / GROUP_FILENAME,
            group_rows,
            empty_fields=("root_id", "native_dynamic_scientific_gate_pass"),
        )
        _write_rows(
            staging_directory / SOLVER_FILENAME,
            solver_rows,
            empty_fields=("root_id", "comparison_kind", "solver_gate_pass"),
        )
        _write_rows(
            staging_directory / ISOMORPHISM_FILENAME,
            isomorphism_rows,
            empty_fields=("root_id", "r2_member_id", "r3_member_id"),
        )
        _write_rows(
            staging_directory / NEARBY_FILENAME,
            [dict(result.metadata) for result in nearby_results],
            empty_fields=(
                "root_id",
                "initial_condition",
                "arm",
                "solver_label",
                "numerical_gate_pass",
            ),
        )
        manifest = _freeze_manifest(
            repository_root,
            results_directory,
            roots,
            group_rows,
            report_directory=staging_directory,
            dependency_snapshot=dependency_snapshot,
            native_ready_summary=native_summary,
        )
        write_json(staging_directory / MANIFEST_FILENAME, manifest, pre_reveal=True)
        write_json(staging_directory / GATE_FILENAME, gate, pre_reveal=True)
        _write_hash_ledger(
            staging_directory,
            DYNAMIC_FREEZE_ARTIFACT_FILENAMES,
        )
        _publish_staged_freeze(
            staging_directory,
            results_directory,
            repository_root,
            dependency_snapshot,
        )
    return gate


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
    gate = run_native_dynamic_contract(
        args.repository_root.resolve(),
        args.results_directory.resolve(),
        workers=args.workers,
    )
    print(json.dumps(gate, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
