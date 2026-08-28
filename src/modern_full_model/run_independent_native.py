"""Independently reproduce Task 13B native-source WT roots.

The runner reads only the pre-reveal native WT root table.  It never imports
the production native-source builder, opens a held-out ledger, or evaluates a
genotype perturbation.  Every supplied physical root is checked with the
second equation transcription and every production-passing root is refined
on the ten-coordinate electroneutral charge manifold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np

from .independent import IndependentConstantStimulus
from .independent_native import (
    AN_ABS_AE4_NHE,
    IndependentNativeSourceDefinition,
    N_ABS_NKCC,
    STRICT_NATIVE,
    independent_native_wt_gate_failures,
    independent_native_wt_observables,
    refine_independent_native_root,
    state_from_independent_rest_coordinates,
)


# Match the production root residual's ten independent rows (alkalinity rows
# 4/10 are omitted through the exact charge identities).  The independent
# local reroot uses the dual valid convention that eliminates chloride; both
# must vanish at a true root, but failed optimizer endpoints are comparable
# only when the same residual rows are reported.
INDEPENDENT_ROWS = (0, 1, 2, 3, 5, 6, 7, 8, 9, 11)
AMOUNT_ROWS = (0, 1, 2, 3, 4, 6, 7, 8, 9, 10)
VOLUME_ROWS = (5, 11)
ROOT_ROW_SCALES = np.asarray(
    (0.05, 0.05, 0.05, 0.02, 1.0e-4, 0.02, 0.02, 0.02, 0.02, 1.0e-4),
    dtype=float,
)

PRODUCTION_ELIGIBILITY = "PRODUCTION_PRE_REVEAL"
PRODUCTION_SOURCE_SCALES = (0.5, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0)
STRESS_SOURCE_SCALES = (12.0, 16.0, 24.0, 32.0)
HYDRAULIC_SCALES = (2.5, 3.5, 5.0, 6.75, 8.0)
PUMP_CONTEXTS = (0.5, 1.0, 2.0)
TOPOLOGIES = {
    "PLOW_KLOW": (0.056878, 0.20),
    "PLOW_KHIGH": (0.056878, 0.40),
    "PNOM_KMID": (0.075075, 0.30),
    "PHIGH_KLOW": (0.094586, 0.20),
    "PHIGH_KHIGH": (0.094586, 0.40),
}
ROUTINGS = {"AE4NA05": 0.05, "AE4NA10": 0.10, "AE4NA20": 0.20}
SCREEN_STAGES = {
    "STRICT_H1_PUMP1",
    "CONDITIONAL_PUMP_AFTER_ZERO_STRICT_ROOTS",
    "DIAGNOSTIC_HYDRAULIC_REROOT",
    "DIAGNOSTIC_HWQ_REROOT",
    "PRODUCTION_ABSOLUTE_SOURCE_GRID",
    "DIAGNOSTIC_ABSOLUTE_SOURCE_STRESS",
}
CONFIRMATION_STAGE = "SEVEN_START_CONFIRMATION"
GEOMETRY_STAGE = "THIRTY_THREE_START_ROOT_GEOMETRY"
CLUSTER_RELATIVE_TOLERANCE = 2.0e-5
ROOT_BOUNDS_LOWER = np.asarray(
    (4.0, 60.0, 1.0, 6.2, 0.5, 20.0, 0.05, 0.2, 5.0, 0.02, 20.0)
)
ROOT_BOUNDS_UPPER = np.asarray(
    (60.0, 210.0, 100.0, 8.0, 3.0, 280.0, 60.0, 100.0, 9.0, 0.80, 300.0)
)
ROOT_BOUNDS_SPAN = ROOT_BOUNDS_UPPER - ROOT_BOUNDS_LOWER


def _csv_snapshot(path: Path) -> tuple[list[dict[str, str]], str]:
    """Parse and hash one immutable byte snapshot, avoiding a second read."""

    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{path.name}: native CSV is not UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if not reader.fieldnames:
        raise ValueError(f"{path.name}: native CSV has no header")
    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError(f"{path.name}: native CSV has duplicate columns")
    return list(reader), digest


def _duplicate_values(rows: Sequence[Mapping[str, Any]], key: str) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for row in rows:
        value = str(row.get(key, ""))
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def _canonical_mode(row: Mapping[str, Any]) -> str:
    value = row.get("hydraulic_mode")
    if value:
        return str(value)
    return "H1" if float(row["hydraulic_scale"]) == 1.0 else "HW"


def _stage_start_count(stage: str) -> int | None:
    if stage in SCREEN_STAGES:
        return 3
    if stage == CONFIRMATION_STAGE:
        return 7
    if stage == GEOMETRY_STAGE:
        return 33
    return None


def _connected_root_components(
    endpoints_by_start: Mapping[str, np.ndarray],
) -> tuple[tuple[str, ...], ...]:
    """Return deterministic transitive root clusters at the declared tolerance."""

    names = tuple(sorted(str(name) for name in endpoints_by_start))
    endpoints: dict[str, np.ndarray] = {}
    for name in names:
        endpoint = np.asarray(endpoints_by_start[name], dtype=float)
        if endpoint.shape != (11,) or np.any(~np.isfinite(endpoint)):
            raise ValueError(f"{name}: root-cluster endpoint must be finite length 11")
        endpoints[name] = endpoint

    adjacency: dict[str, set[str]] = {name: set() for name in names}
    for index, left_name in enumerate(names):
        left = endpoints[left_name]
        for right_name in names[index + 1 :]:
            right = endpoints[right_name]
            distance = float(
                np.max(np.abs((left - right) / ROOT_BOUNDS_SPAN))
            )
            if distance <= CLUSTER_RELATIVE_TOLERANCE:
                adjacency[left_name].add(right_name)
                adjacency[right_name].add(left_name)

    components: list[tuple[str, ...]] = []
    unseen = set(names)
    while unseen:
        seed = min(unseen)
        pending = [seed]
        component: set[str] = set()
        while pending:
            current = pending.pop()
            if current in component:
                continue
            component.add(current)
            unseen.discard(current)
            pending.extend(sorted(adjacency[current] - component, reverse=True))
        components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def _accepted_legacy_diagnostic_start_design(
    panel: Mapping[str, Any],
    design: Mapping[str, Any],
    actual_start_ids: Sequence[str],
) -> bool:
    """Accept only the frozen diagnostic-only legacy three-start displacement."""

    ids = tuple(str(value) for value in actual_start_ids)
    try:
        hydraulic_scale = float(panel["hydraulic_scale"])
        wt_pass_count = int(panel.get("wt_passing_root_count") or 0)
    except (KeyError, TypeError, ValueError):
        return False
    mode = _canonical_mode(panel)
    return bool(
        design.get("kind")
        == "LEGACY_PRIORITY_TRUNCATION_UPPER_QUARTILE_DISPLACED"
        and "DIAGNOSTIC_ONLY" in str(panel.get("eligibility", ""))
        and hydraulic_scale > 1.0
        and mode in {"HW", "HWQ"}
        and wt_pass_count == 0
        and len(ids) == 3
        and int(design.get("count", -1)) == 3
        and "frozen_reference" in ids
        and "lower_quartile" in ids
        and "upper_quartile" not in ids
        and any(value.startswith("continuation_") for value in ids)
        and bool(design.get("actual_persisted_attempt_set", False))
        and bool(design.get("legacy_priority_truncation", False))
    )


def _panel_id(
    source_class: str,
    source_scale: float,
    topology_id: str,
    routing_id: str,
    pump_scale: float,
    hydraulic_scale: float,
    hydraulic_mode: str,
) -> str:
    hydraulic_id = (
        f"HWQ_H{hydraulic_scale:g}"
        if hydraulic_mode == "HWQ"
        else f"H{hydraulic_scale:g}"
    )
    return (
        f"{source_class}_S{source_scale:g}_{topology_id}_{routing_id}_"
        f"P{pump_scale:g}_{hydraulic_id}"
    )


def _canonical_panel_definitions() -> dict[str, dict[str, Any]]:
    definitions: dict[str, dict[str, Any]] = {}
    for pump in PUMP_CONTEXTS:
        for topology_id, (pump_fraction, k_fraction) in TOPOLOGIES.items():
            for routing_id, routing_fraction in ROUTINGS.items():
                base = {
                    "topology_id": topology_id,
                    "routing_id": routing_id,
                    "pump_capacity_scale": pump,
                    "apical_pump_fraction": pump_fraction,
                    "apical_k_fraction": k_fraction,
                    "ae4_cation_fraction": routing_fraction,
                }

                def add(
                    source_class: str,
                    source_scale: float,
                    hydraulic_scale: float,
                    hydraulic_mode: str,
                    eligibility: str,
                    kind: str,
                ) -> None:
                    panel_id = _panel_id(
                        source_class,
                        source_scale,
                        topology_id,
                        routing_id,
                        pump,
                        hydraulic_scale,
                        hydraulic_mode,
                    )
                    strict_pump_one_id = _panel_id(
                        STRICT_NATIVE,
                        1.0,
                        topology_id,
                        routing_id,
                        1.0,
                        1.0,
                        "H1",
                    )
                    strict_same_pump_id = _panel_id(
                        STRICT_NATIVE,
                        1.0,
                        topology_id,
                        routing_id,
                        pump,
                        1.0,
                        "H1",
                    )
                    if source_class != STRICT_NATIVE or hydraulic_scale != 1.0:
                        parent_panel_id = strict_same_pump_id
                    elif pump != 1.0:
                        parent_panel_id = strict_pump_one_id
                    else:
                        parent_panel_id = "ROOT_GENERATION_NATIVE_SOURCE_MAP"
                    definitions[panel_id] = {
                        **base,
                        "source_class": source_class,
                        "source_scale": source_scale,
                        "hydraulic_scale": hydraulic_scale,
                        "hydraulic_mode": hydraulic_mode,
                        "eligibility": eligibility,
                        "parent_panel_id": parent_panel_id,
                        "kind": kind,
                    }

                add(STRICT_NATIVE, 1.0, 1.0, "H1", PRODUCTION_ELIGIBILITY, "strict")
                for mode in ("HW", "HWQ"):
                    for hydraulic in HYDRAULIC_SCALES:
                        add(
                            STRICT_NATIVE,
                            1.0,
                            hydraulic,
                            mode,
                            f"DIAGNOSTIC_ONLY_{mode}_HYDRAULIC_TRANSFER",
                            mode.lower(),
                        )
                for source_class in (N_ABS_NKCC, AN_ABS_AE4_NHE):
                    for source_scale in PRODUCTION_SOURCE_SCALES:
                        add(
                            source_class,
                            source_scale,
                            1.0,
                            "H1",
                            PRODUCTION_ELIGIBILITY,
                            "absolute",
                        )
                    for source_scale in STRESS_SOURCE_SCALES:
                        add(
                            source_class,
                            source_scale,
                            1.0,
                            "H1",
                            "DIAGNOSTIC_ONLY_ABSOLUTE_SOURCE_STRESS",
                            "stress",
                        )
    return definitions


def _truth(value: Any) -> bool:
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(
        "native artifact booleans must be bool or canonical case-sensitive "
        f"True/False CSV text, found {value!r}"
    )


def _decoded(row: Mapping[str, str], key: str) -> Any:
    value = json.loads(row[key])
    return value


def reproduce_native_roots(root_table: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    supplied_rows, root_digest = _csv_snapshot(root_table)
    return _reproduce_native_root_rows(
        supplied_rows,
        root_table_name=root_table.name,
        root_table_sha256=root_digest,
    )


def _reproduce_native_root_rows(
    supplied_rows: Sequence[Mapping[str, str]],
    *,
    root_table_name: str,
    root_table_sha256: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not supplied_rows:
        raise ValueError("native WT root table is empty")
    duplicate_root_ids = _duplicate_values(supplied_rows, "root_id")

    records: list[dict[str, Any]] = []
    for row in supplied_rows:
        root_id = row["root_id"]
        production_numerical = _truth(row["passes_numerical_gate"])
        production_wt = _truth(row["passes_wt_gate"])
        definition = IndependentNativeSourceDefinition.from_root_row(row)
        model = definition.build_independent_model(
            stimulus=IndependentConstantStimulus(0.058, 0.0)
        )
        candidate = np.asarray(_decoded(row, "core_state_json"), dtype=float)
        if candidate.shape != (12,) or np.any(~np.isfinite(candidate)) or np.any(candidate <= 0.0):
            raise ValueError(f"{root_id}: supplied core state is not positive finite length 12")
        evaluated = model.evaluate(0.0, candidate)
        raw = np.asarray(evaluated.rhs, dtype=float)
        independent_scaled_max = float(
            np.max(np.abs(raw[list(INDEPENDENT_ROWS)] / ROOT_ROW_SCALES))
        )
        amount_max = float(np.max(np.abs(raw[list(AMOUNT_ROWS)])))
        volume_max = float(np.max(np.abs(raw[list(VOLUME_ROWS)])))
        current_max = float(max(abs(value) for value in evaluated.current_residuals_A.values()))
        conservation_max = float(
            max(abs(value) for value in evaluated.conservation_residuals.values())
        )
        independent_observables = independent_native_wt_observables(evaluated, candidate)
        independent_wt_failures = independent_native_wt_gate_failures(
            independent_observables
        )
        declared_observables = _decoded(row, "observables_json")
        if not isinstance(declared_observables, Mapping):
            raise ValueError(f"{root_id}: observables_json must decode to a mapping")
        observable_differences: dict[str, float] = {}
        observable_integrity_issues: list[str] = []
        for key, independent_value in independent_observables.items():
            if key not in declared_observables:
                observable_integrity_issues.append(f"{key}:missing")
                continue
            try:
                declared_value = float(declared_observables[key])
            except (TypeError, ValueError):
                observable_integrity_issues.append(f"{key}:invalid")
                continue
            difference = abs(declared_value - independent_value)
            observable_differences[key] = difference
            if not math.isfinite(declared_value) or difference > 1.0e-8:
                observable_integrity_issues.append(f"{key}:mismatch")

        refine_requested = bool(production_numerical or production_wt)
        if refine_requested:
            refined = refine_independent_native_root(
                definition,
                candidate,
                relative_bound=0.01,
                max_nfev=500,
                row_scales=ROOT_ROW_SCALES,
            )
            refine_success: bool | None = bool(refined.success)
            refine_scaled_max: float | None = float(refined.max_abs_scaled_residual)
            rank: int | None = int(refined.jacobian_rank)
            nullity: int | None = int(refined.jacobian_nullity)
            state_difference: float | None = float(
                np.max(np.abs(refined.state[:12] - candidate))
            )
        else:
            refine_success = None
            refine_scaled_max = None
            rank = None
            nullity = None
            state_difference = None

        independent_pass = bool(
            independent_scaled_max <= 1.0e-6
            and amount_max <= 1.0e-8
            and volume_max <= 1.0e-10
            and current_max <= 1.0e-18
            and conservation_max <= 1.0e-8
            and (not refine_requested or refine_success)
        )
        independent_wt_pass = bool(
            independent_pass and not independent_wt_failures
        )
        hydraulic_scale = float(row["hydraulic_scale"])
        hydraulic_mode = row.get("hydraulic_mode") or (
            "H1" if hydraulic_scale == 1.0 else "HW"
        )
        production_h1_eligible = bool(
            row.get("eligibility") == PRODUCTION_ELIGIBILITY
            and hydraulic_scale == 1.0
            and hydraulic_mode == "H1"
        )
        records.append(
            {
                "root_id": root_id,
                "panel_id": row["panel_id"],
                "source_class": row["source_class"],
                "source_scale": float(row["source_scale"]),
                "pump_capacity_scale": float(row["pump_capacity_scale"]),
                "hydraulic_scale": hydraulic_scale,
                "hydraulic_mode": hydraulic_mode,
                "eligibility": row.get("eligibility", ""),
                "production_h1_eligible_marker": production_h1_eligible,
                "production_passes_numerical_gate": production_numerical,
                "production_passes_wt_gate": production_wt,
                "independent_direct_scaled_max": independent_scaled_max,
                "independent_amount_source_max_fmol_s": amount_max,
                "independent_volume_source_max_pL_s": volume_max,
                "independent_current_residual_max_A": current_max,
                "independent_conservation_residual_max": conservation_max,
                "independent_refine_requested": refine_requested,
                "independent_refine_success": refine_success,
                "independent_refine_scaled_max": refine_scaled_max,
                "independent_jacobian_rank": rank,
                "independent_jacobian_nullity": nullity,
                "independent_refined_state_max_abs_difference": state_difference,
                "independent_pass": independent_pass,
                "independent_wt_observables_json": json.dumps(
                    independent_observables, sort_keys=True, separators=(",", ":")
                ),
                "independent_wt_gate_failures_json": json.dumps(
                    independent_wt_failures, separators=(",", ":")
                ),
                "independent_wt_pass": independent_wt_pass,
                "root_table_observable_max_abs_difference": (
                    max(observable_differences.values())
                    if observable_differences
                    else None
                ),
                "root_table_observable_integrity_issues_json": json.dumps(
                    observable_integrity_issues, separators=(",", ":")
                ),
                "production_independent_gate_agree": bool(
                    production_numerical == independent_pass
                ),
                "production_independent_wt_gate_agree": bool(
                    production_wt == independent_wt_pass
                ),
            }
        )

    passing_production = [row for row in records if row["production_passes_numerical_gate"]]
    passing_independent = [row for row in records if row["independent_pass"]]
    disagreements = [row for row in records if not row["production_independent_gate_agree"]]
    wt_disagreements = [
        row
        for row in records
        if not row["production_independent_wt_gate_agree"]
    ]
    observable_disagreements = [
        row
        for row in records
        if json.loads(row["root_table_observable_integrity_issues_json"])
    ]
    counts_by_source_and_mode: dict[str, dict[str, int]] = {}
    for source_class, hydraulic_mode in sorted(
        {(row["source_class"], row["hydraulic_mode"]) for row in records}
    ):
        selected = [
            row
            for row in records
            if row["source_class"] == source_class
            and row["hydraulic_mode"] == hydraulic_mode
        ]
        counts_by_source_and_mode[f"{source_class}::{hydraulic_mode}"] = {
            "root_count": len(selected),
            "production_numerical_pass_count": sum(
                bool(row["production_passes_numerical_gate"]) for row in selected
            ),
            "production_wt_pass_count": sum(
                bool(row["production_passes_wt_gate"]) for row in selected
            ),
            "independent_pass_count": sum(
                bool(row["independent_pass"]) for row in selected
            ),
        }
    summary = {
        "status": (
            "INDEPENDENT_NATIVE_ROOT_REPRODUCTION_PASS"
            if not disagreements
            and not wt_disagreements
            and not observable_disagreements
            and not duplicate_root_ids
            else "INDEPENDENT_NATIVE_ROOT_REPRODUCTION_DISCREPANCY"
        ),
        "firewall": "PRE_REVEAL_WT_ONLY_NO_GENOTYPE_NO_HELDOUT_ACCESS",
        "root_table": root_table_name,
        "root_table_sha256": root_table_sha256,
        "duplicate_root_ids": duplicate_root_ids,
        "row_count": len(records),
        "production_numerical_pass_count": len(passing_production),
        "production_wt_pass_count": sum(
            bool(row["production_passes_wt_gate"]) for row in records
        ),
        "independent_pass_count": len(passing_independent),
        "gate_disagreement_count": len(disagreements),
        "gate_disagreement_root_ids": [row["root_id"] for row in disagreements],
        "wt_gate_disagreement_count": len(wt_disagreements),
        "wt_gate_disagreement_root_ids": [
            row["root_id"] for row in wt_disagreements
        ],
        "observable_integrity_disagreement_count": len(observable_disagreements),
        "observable_integrity_disagreement_root_ids": [
            row["root_id"] for row in observable_disagreements
        ],
        "root_counts_by_source_class_and_hydraulic_mode": counts_by_source_and_mode,
        "maximum_direct_scaled_residual_among_production_passes": (
            max(row["independent_direct_scaled_max"] for row in passing_production)
            if passing_production
            else None
        ),
        "maximum_refined_scaled_residual_among_production_passes": (
            max(row["independent_refine_scaled_max"] for row in passing_production)
            if passing_production
            else None
        ),
        "all_production_passing_jacobians_rank_eleven_nullity_zero": bool(
            passing_production
            and all(
                row["independent_jacobian_rank"] == 11
                and row["independent_jacobian_nullity"] == 0
                for row in passing_production
            )
        ),
    }
    return records, summary


def reproduce_native_attempt_residuals(
    panel_table: Path,
    attempt_table: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Re-evaluate every production optimizer endpoint with independent equations."""

    panel_rows, panel_digest = _csv_snapshot(panel_table)
    attempts, attempt_digest = _csv_snapshot(attempt_table)
    return _reproduce_native_attempt_rows(
        panel_rows,
        attempts,
        panel_table_sha256=panel_digest,
        attempt_table_sha256=attempt_digest,
    )


def _reproduce_native_attempt_rows(
    panel_rows: Sequence[Mapping[str, str]],
    attempts: Sequence[Mapping[str, str]],
    *,
    panel_table_sha256: str,
    attempt_table_sha256: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    panels = {row["panel_id"]: row for row in panel_rows}
    duplicate_panel_ids = _duplicate_values(panel_rows, "panel_id")
    duplicate_attempt_ids = _duplicate_values(attempts, "attempt_id")
    records: list[dict[str, Any]] = []
    for row in attempts:
        if row["panel_id"] not in panels:
            raise ValueError(
                f"{row.get('attempt_id', '<unknown>')}: attempt references unknown panel"
            )
        panel = panels[row["panel_id"]]
        vector = np.asarray(json.loads(row["coordinates_and_other_json"]), dtype=float)
        if vector.shape != (11,) or np.any(~np.isfinite(vector)):
            independent_admissible = False
            scaled_max = math.inf
            amount_max = math.inf
            volume_max = math.inf
        else:
            definition = IndependentNativeSourceDefinition(
                source_class=panel["source_class"],
                source_scale=float(panel["source_scale"]),
                pump_capacity_scale=float(panel["pump_capacity_scale"]),
                apical_pump_fraction=float(panel["apical_pump_fraction"]),
                apical_k_fraction=float(panel["apical_k_fraction"]),
                ae4_cation_fraction=float(panel["ae4_cation_fraction"]),
                hydraulic_scale=float(panel["hydraulic_scale"]),
                cell_other_impermeant_osmoles_fmol=float(vector[-1]),
                hydraulic_mode=str(
                    panel.get("hydraulic_mode")
                    or ("H1" if float(panel["hydraulic_scale"]) == 1.0 else "HW")
                ),
            )
            parameters, _ae4 = definition.reconstruct()
            try:
                state = state_from_independent_rest_coordinates(
                    parameters,
                    vector[:10],
                )
                model = definition.build_independent_model(
                    stimulus=IndependentConstantStimulus(0.058, 0.0)
                )
                evaluation = model.evaluate(0.0, state)
                raw = np.asarray(evaluation.rhs, dtype=float)
                residual = np.r_[
                    raw[list(INDEPENDENT_ROWS)] / ROOT_ROW_SCALES,
                    (state[5] - 1.3) / 0.13,
                ]
                independent_admissible = bool(np.all(state > 0.0))
                scaled_max = float(np.max(np.abs(residual)))
                amount_max = float(np.max(np.abs(raw[list(AMOUNT_ROWS)])))
                volume_max = float(np.max(np.abs(raw[list(VOLUME_ROWS)])))
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
                independent_admissible = False
                scaled_max = math.inf
                amount_max = math.inf
                volume_max = math.inf
        production_scaled = float(row["max_abs_scaled_residual"])
        production_converged = _truth(row["converged_root"])
        independent_converged = bool(
            independent_admissible and scaled_max <= 1.0e-7
        )
        records.append(
            {
                "attempt_id": row["attempt_id"],
                "panel_id": row["panel_id"],
                "source_class": row["source_class"],
                "source_scale": float(row["source_scale"]),
                "pump_capacity_scale": float(row["pump_capacity_scale"]),
                "hydraulic_scale": float(row["hydraulic_scale"]),
                "hydraulic_mode": panel.get("hydraulic_mode") or (
                    "H1" if float(row["hydraulic_scale"]) == 1.0 else "HW"
                ),
                "eligibility": row["eligibility"],
                "production_admissible": _truth(row["admissible_state"]),
                "independent_admissible": independent_admissible,
                "production_converged_root": production_converged,
                "production_passes_numerical_gate": _truth(
                    row["passes_numerical_gate"]
                ),
                "independent_converged_root": independent_converged,
                "production_scaled_max": production_scaled,
                "independent_scaled_max": scaled_max,
                "scaled_max_absolute_difference": abs(production_scaled - scaled_max),
                "independent_amount_source_max_fmol_s": amount_max,
                "independent_volume_source_max_pL_s": volume_max,
                "convergence_gate_agree": bool(
                    production_converged == independent_converged
                ),
            }
        )

    disagreements = [row for row in records if not row["convergence_gate_agree"]]
    finite_differences = [
        row["scaled_max_absolute_difference"]
        for row in records
        if math.isfinite(row["scaled_max_absolute_difference"])
    ]
    by_class: dict[str, dict[str, int]] = {}
    for source_class in sorted({row["source_class"] for row in records}):
        selected = [row for row in records if row["source_class"] == source_class]
        by_class[source_class] = {
            "attempt_count": len(selected),
            "production_converged_count": sum(
                bool(row["production_converged_root"]) for row in selected
            ),
            "independent_converged_count": sum(
                bool(row["independent_converged_root"]) for row in selected
            ),
        }
    clusters_by_panel: dict[str, int] = {}
    cluster_members_by_panel: dict[str, list[list[str]]] = {}
    cluster_disagreements: list[str] = []
    for panel_id, panel in panels.items():
        endpoints = {
            str(attempt["start_id"]): np.asarray(
                json.loads(attempt["coordinates_and_other_json"]), dtype=float
            )
            for attempt in attempts
            if attempt["panel_id"] == panel_id
            and _truth(attempt["passes_numerical_gate"])
        }
        components = _connected_root_components(endpoints)
        clusters_by_panel[panel_id] = len(components)
        cluster_members_by_panel[panel_id] = [list(item) for item in components]
        if len(components) != int(panel["numerical_root_count"]):
            cluster_disagreements.append(panel_id)

    attempts_by_panel: dict[str, list[Mapping[str, str]]] = {}
    for attempt in attempts:
        attempts_by_panel.setdefault(attempt["panel_id"], []).append(attempt)
    panel_metadata_disagreements: list[dict[str, Any]] = []
    confirmation_panels: list[str] = []
    geometry_panels: list[str] = []
    accepted_legacy_diagnostic_panels: list[str] = []
    for panel_id, panel in panels.items():
        actual_attempts = attempts_by_panel.get(panel_id, [])
        actual_ids = tuple(sorted(row["start_id"] for row in actual_attempts))
        declared_attempt_count = int(panel["attempt_count"])
        design = json.loads(panel["start_design_json"])
        declared_design_ids = tuple(sorted(str(value) for value in design["start_ids"]))
        issues: list[str] = []
        legacy_accepted = _accepted_legacy_diagnostic_start_design(
            panel, design, actual_ids
        )
        if legacy_accepted:
            accepted_legacy_diagnostic_panels.append(panel_id)
        elif design.get("kind") == (
            "LEGACY_PRIORITY_TRUNCATION_UPPER_QUARTILE_DISPLACED"
        ):
            issues.append("legacy_start_design_not_eligible")
        stage = str(panel["stage"])
        expected_current_count = _stage_start_count(stage)
        if expected_current_count is None:
            issues.append("unknown_current_stage")
        elif len(actual_attempts) != expected_current_count:
            issues.append(f"current_attempt_count_not_exact_{expected_current_count}")
        if len(actual_attempts) != declared_attempt_count:
            issues.append("attempt_count")
        if len(actual_ids) != len(set(actual_ids)):
            issues.append("duplicate_start_ids")
        if int(design["count"]) != len(actual_attempts):
            issues.append("start_design_count")
        if actual_ids != declared_design_ids:
            issues.append("start_design_ids")
        if not bool(design.get("actual_persisted_attempt_set", False)):
            issues.append("actual_persisted_attempt_set")
        if not legacy_accepted and not bool(
            design.get("both_quartiles_present", False)
        ):
            issues.append("both_quartiles_present")
        if not legacy_accepted and not {
            "lower_quartile",
            "upper_quartile",
        }.issubset(actual_ids):
            issues.append("quartile_start_ids")
        if stage in {CONFIRMATION_STAGE, GEOMETRY_STAGE} and "frozen_reference" not in actual_ids:
            issues.append("frozen_reference_missing_from_deep_search")
        if stage in SCREEN_STAGES:
            primary_ids = {
                value
                for value in actual_ids
                if value in {"frozen_reference", "parent_continuation_00"}
            }
            if len(primary_ids) != 1:
                issues.append("screen_primary_start")
        screen_count = int(panel.get("screen_start_count") or 0)
        confirmation_count = int(panel.get("confirmation_start_count") or 0)
        geometry_count = int(panel.get("geometry_start_count") or 0)
        if screen_count != 3:
            issues.append("screen_start_count_not_exact_3")
        expected_confirmation = (
            7 if stage in {CONFIRMATION_STAGE, GEOMETRY_STAGE} else 0
        )
        expected_geometry = 33 if stage == GEOMETRY_STAGE else 0
        if confirmation_count != expected_confirmation:
            issues.append("confirmation_start_count_not_current_exact_7_or_0")
        if geometry_count != expected_geometry:
            issues.append("geometry_start_count_not_current_exact_33_or_0")
        if confirmation_count:
            confirmation_panels.append(panel_id)
            if len(actual_attempts) != confirmation_count and stage == CONFIRMATION_STAGE:
                issues.append("confirmation_attempts_not_retained")
        if geometry_count:
            geometry_panels.append(panel_id)
            if len(actual_attempts) != geometry_count:
                issues.append("geometry_attempts_not_retained")
        if issues:
            panel_metadata_disagreements.append(
                {"panel_id": panel_id, "issues": issues}
            )

    panel_counts_by_source_mode_stage: dict[str, int] = {}
    for panel in panels.values():
        key = "::".join(
            (
                panel["source_class"],
                panel.get("hydraulic_mode") or (
                    "H1" if float(panel["hydraulic_scale"]) == 1.0 else "HW"
                ),
                panel["stage"],
            )
        )
        panel_counts_by_source_mode_stage[key] = (
            panel_counts_by_source_mode_stage.get(key, 0) + 1
        )
    summary = {
        "attempt_count": len(records),
        "convergence_gate_disagreement_count": len(disagreements),
        "convergence_gate_disagreement_attempt_ids": [
            row["attempt_id"] for row in disagreements
        ],
        "maximum_scaled_residual_absolute_difference": (
            max(finite_differences) if finite_differences else None
        ),
        "counts_by_source_class": by_class,
        "independent_cluster_count_by_panel": clusters_by_panel,
        "independent_cluster_members_by_panel": cluster_members_by_panel,
        "root_cluster_count_disagreement_panel_ids": cluster_disagreements,
        "panel_count_by_source_mode_stage": panel_counts_by_source_mode_stage,
        "confirmation_panel_count": len(confirmation_panels),
        "confirmation_panel_ids": sorted(confirmation_panels),
        "geometry_panel_count": len(geometry_panels),
        "geometry_panel_ids": sorted(geometry_panels),
        "accepted_legacy_diagnostic_start_design_count": len(
            accepted_legacy_diagnostic_panels
        ),
        "accepted_legacy_diagnostic_start_design_panel_ids": sorted(
            accepted_legacy_diagnostic_panels
        ),
        "panel_start_design_disagreement_count": len(panel_metadata_disagreements),
        "panel_start_design_disagreements": panel_metadata_disagreements,
        "duplicate_panel_ids": duplicate_panel_ids,
        "duplicate_attempt_ids": duplicate_attempt_ids,
        "all_declared_confirmation_panels_have_exactly_7_current_starts": bool(
            confirmation_panels
            and all(
                int(panels[panel_id].get("confirmation_start_count") or 0) == 7
                and panels[panel_id]["stage"] in {CONFIRMATION_STAGE, GEOMETRY_STAGE}
                for panel_id in confirmation_panels
            )
        ),
        "all_declared_geometry_panels_have_exactly_33_current_starts": bool(
            geometry_panels
            and all(
                int(panels[panel_id].get("geometry_start_count") or 0) == 33
                and panels[panel_id]["stage"] == GEOMETRY_STAGE
                for panel_id in geometry_panels
            )
        ),
        "panel_table_sha256": panel_table_sha256,
        "attempt_table_sha256": attempt_table_sha256,
    }
    return records, summary


def audit_native_hierarchy(
    panel_rows: Sequence[Mapping[str, str]],
    root_rows: Sequence[Mapping[str, str]],
    attempt_rows: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    """Fail-closed structural audit of the current three-table root hierarchy."""

    issues: list[dict[str, str]] = []

    def issue(scope: str, item_id: str, name: str) -> None:
        issues.append({"scope": scope, "id": item_id, "issue": name})

    for key, rows, scope in (
        ("panel_id", panel_rows, "panel"),
        ("root_id", root_rows, "root"),
        ("attempt_id", attempt_rows, "attempt"),
    ):
        for duplicate in _duplicate_values(rows, key):
            issue(scope, duplicate, f"duplicate_{key}")

    panels = {str(row["panel_id"]): row for row in panel_rows}
    roots_by_panel: dict[str, list[Mapping[str, str]]] = {}
    attempts_by_panel: dict[str, list[Mapping[str, str]]] = {}
    accepted_legacy_diagnostic_panels: list[str] = []
    for row in root_rows:
        roots_by_panel.setdefault(str(row["panel_id"]), []).append(row)
    for row in attempt_rows:
        attempts_by_panel.setdefault(str(row["panel_id"]), []).append(row)

    canonical = _canonical_panel_definitions()
    numeric_fields = (
        "source_scale",
        "pump_capacity_scale",
        "apical_pump_fraction",
        "apical_k_fraction",
        "ae4_cation_fraction",
        "hydraulic_scale",
    )
    text_fields = (
        "parent_panel_id",
        "topology_id",
        "routing_id",
        "source_class",
        "eligibility",
    )
    for panel_id, panel in panels.items():
        expected = canonical.get(panel_id)
        if expected is None:
            issue("panel", panel_id, "panel_id_not_in_independent_canonical_registry")
            continue
        for field in numeric_fields:
            try:
                matches = float(panel[field]) == float(expected[field])
            except (KeyError, TypeError, ValueError):
                matches = False
            if not matches:
                issue("panel", panel_id, f"canonical_{field}_mismatch")
        for field in text_fields:
            if str(panel.get(field, "")) != str(expected[field]):
                issue("panel", panel_id, f"canonical_{field}_mismatch")
        try:
            mode = _canonical_mode(panel)
        except (KeyError, TypeError, ValueError):
            mode = "INVALID"
        if mode != expected["hydraulic_mode"]:
            issue("panel", panel_id, "canonical_hydraulic_mode_mismatch")
        expected_screen_stage = {
            "strict": (
                "STRICT_H1_PUMP1"
                if expected["pump_capacity_scale"] == 1.0
                else "CONDITIONAL_PUMP_AFTER_ZERO_STRICT_ROOTS"
            ),
            "hw": "DIAGNOSTIC_HYDRAULIC_REROOT",
            "hwq": "DIAGNOSTIC_HWQ_REROOT",
            "absolute": "PRODUCTION_ABSOLUTE_SOURCE_GRID",
            "stress": "DIAGNOSTIC_ABSOLUTE_SOURCE_STRESS",
        }[expected["kind"]]
        allowed_stages = {expected_screen_stage}
        if expected["kind"] in {"strict", "absolute"}:
            allowed_stages.update({CONFIRMATION_STAGE, GEOMETRY_STAGE})
        if str(panel.get("stage", "")) not in allowed_stages:
            issue("panel", panel_id, "stage_incompatible_with_canonical_panel_kind")

    base_strict_ids = {
        panel_id
        for panel_id, expected in canonical.items()
        if expected["kind"] == "strict"
        and expected["pump_capacity_scale"] == 1.0
    }
    for panel_id in sorted(base_strict_ids - set(panels)):
        issue("hierarchy", panel_id, "missing_strict_pump1_panel")
    base_strict_pass = any(
        int(panels[panel_id].get("wt_passing_root_count") or 0) > 0
        for panel_id in base_strict_ids & set(panels)
    )
    licensed_pumps = (1.0,) if base_strict_pass else PUMP_CONTEXTS
    required_ids = {
        panel_id
        for panel_id, expected in canonical.items()
        if expected["kind"] in {"strict", "hw", "hwq", "absolute"}
        and expected["pump_capacity_scale"] in licensed_pumps
    }
    optional_ids = {
        panel_id
        for panel_id, expected in canonical.items()
        if expected["kind"] == "stress"
        and expected["pump_capacity_scale"] in licensed_pumps
    }
    for panel_id in sorted(required_ids - set(panels)):
        issue("hierarchy", panel_id, "missing_required_panel")
    for panel_id in sorted(set(panels) - required_ids - optional_ids):
        expected = canonical.get(panel_id)
        if expected is not None and expected["pump_capacity_scale"] not in licensed_pumps:
            issue("hierarchy", panel_id, "unlicensed_pump_context_panel")
        elif expected is not None:
            issue("hierarchy", panel_id, "unexpected_panel_for_current_hierarchy")

    production_panel_ids = {
        panel_id
        for panel_id in required_ids & set(panels)
        if canonical[panel_id]["kind"] in {"strict", "absolute"}
        and canonical[panel_id]["eligibility"] == PRODUCTION_ELIGIBILITY
        and canonical[panel_id]["hydraulic_scale"] == 1.0
        and canonical[panel_id]["hydraulic_mode"] == "H1"
    }

    attempt_metadata_fields = (
        "source_class",
        "source_scale",
        "pump_capacity_scale",
        "hydraulic_scale",
        "eligibility",
    )
    root_metadata_fields = (
        *attempt_metadata_fields,
        "topology_id",
        "routing_id",
        "apical_pump_fraction",
        "apical_k_fraction",
        "ae4_cation_fraction",
    )
    for attempt in attempt_rows:
        attempt_id = str(attempt.get("attempt_id", ""))
        panel_id = str(attempt.get("panel_id", ""))
        panel = panels.get(panel_id)
        if panel is None:
            issue("attempt", attempt_id, "orphan_panel_id")
            continue
        expected_attempt_id = f"{panel_id}|{attempt.get('start_id', '')}"
        if attempt_id != expected_attempt_id:
            issue("attempt", attempt_id, "attempt_id_not_panel_plus_start_id")
        for field in attempt_metadata_fields:
            if field in {"source_scale", "pump_capacity_scale", "hydraulic_scale"}:
                try:
                    matches = float(attempt[field]) == float(panel[field])
                except (KeyError, TypeError, ValueError):
                    matches = False
            else:
                matches = str(attempt.get(field, "")) == str(panel.get(field, ""))
            if not matches:
                issue("attempt", attempt_id, f"panel_{field}_mismatch")
        if _canonical_mode(attempt) != _canonical_mode(panel):
            issue("attempt", attempt_id, "panel_hydraulic_mode_mismatch")

    for root in root_rows:
        root_id = str(root.get("root_id", ""))
        panel_id = str(root.get("panel_id", ""))
        panel = panels.get(panel_id)
        if panel is None:
            issue("root", root_id, "orphan_panel_id")
            continue
        if not root_id.startswith(f"{panel_id}:B"):
            issue("root", root_id, "root_id_not_panel_branch")
        for field in root_metadata_fields:
            if field in {
                "source_scale",
                "pump_capacity_scale",
                "hydraulic_scale",
                "apical_pump_fraction",
                "apical_k_fraction",
                "ae4_cation_fraction",
            }:
                try:
                    matches = float(root[field]) == float(panel[field])
                except (KeyError, TypeError, ValueError):
                    matches = False
            else:
                matches = str(root.get(field, "")) == str(panel.get(field, ""))
            if not matches:
                issue("root", root_id, f"panel_{field}_mismatch")
        if _canonical_mode(root) != _canonical_mode(panel):
            issue("root", root_id, "panel_hydraulic_mode_mismatch")
        if not _truth(root["passes_numerical_gate"]):
            issue("root", root_id, "root_table_contains_non_numerical_root")
        if _truth(root["passes_wt_gate"]) and not _truth(
            root["passes_numerical_gate"]
        ):
            issue("root", root_id, "wt_pass_without_numerical_pass")
        expected_count = _stage_start_count(str(panel["stage"]))
        if str(root.get("root_search_stage", "")) != str(panel["stage"]):
            issue("root", root_id, "root_search_stage_not_current_panel_stage")
        try:
            root_start_count = int(root.get("root_search_start_count") or 0)
        except (TypeError, ValueError):
            root_start_count = -1
        if expected_count is None or root_start_count != expected_count:
            issue("root", root_id, "root_search_start_count_not_current_exact_stage")
        for root_field, panel_field in (
            ("panel_confirmation_start_count", "confirmation_start_count"),
            ("panel_geometry_start_count", "geometry_start_count"),
        ):
            try:
                matches = int(root.get(root_field) or 0) == int(
                    panel.get(panel_field) or 0
                )
            except (TypeError, ValueError):
                matches = False
            if not matches:
                issue("root", root_id, f"{root_field}_mismatch")

    for panel_id, panel in panels.items():
        panel_attempts = attempts_by_panel.get(panel_id, [])
        panel_roots = roots_by_panel.get(panel_id, [])
        stage = str(panel.get("stage", ""))
        expected_current = _stage_start_count(stage)
        actual_start_ids = tuple(
            sorted(str(row.get("start_id", "")) for row in panel_attempts)
        )
        try:
            design = json.loads(str(panel["start_design_json"]))
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            design = {}
            issue("panel", panel_id, "invalid_start_design_json")
        if not isinstance(design, Mapping):
            design = {}
            issue("panel", panel_id, "start_design_not_a_mapping")
        legacy_accepted = _accepted_legacy_diagnostic_start_design(
            panel, design, actual_start_ids
        )
        if legacy_accepted:
            accepted_legacy_diagnostic_panels.append(panel_id)
        elif design.get("kind") == (
            "LEGACY_PRIORITY_TRUNCATION_UPPER_QUARTILE_DISPLACED"
        ):
            issue("panel", panel_id, "legacy_start_design_not_eligible")
        try:
            declared_design_ids = tuple(
                sorted(str(value) for value in design["start_ids"])
            )
            design_count = int(design["count"])
        except (KeyError, TypeError, ValueError):
            declared_design_ids = ()
            design_count = -1
        if declared_design_ids != actual_start_ids:
            issue("panel", panel_id, "start_design_ids_cross_table_mismatch")
        if design_count != len(panel_attempts):
            issue("panel", panel_id, "start_design_count_cross_table_mismatch")
        if not bool(design.get("actual_persisted_attempt_set", False)):
            issue("panel", panel_id, "start_design_not_actual_persisted_attempt_set")
        if not legacy_accepted:
            if not {"lower_quartile", "upper_quartile"}.issubset(
                actual_start_ids
            ):
                issue("panel", panel_id, "canonical_quartile_starts_missing")
            if not bool(design.get("both_quartiles_present", False)):
                issue("panel", panel_id, "canonical_both_quartiles_flag_false")
        try:
            declared_attempt_count = int(panel.get("attempt_count") or 0)
        except (TypeError, ValueError):
            declared_attempt_count = -1
        if expected_current is None or len(panel_attempts) != expected_current:
            issue("panel", panel_id, "current_attempt_set_not_exact_stage_count")
        if declared_attempt_count != len(panel_attempts):
            issue("panel", panel_id, "attempt_count_cross_table_mismatch")
        if len({str(row.get("start_id", "")) for row in panel_attempts}) != len(
            panel_attempts
        ):
            issue("panel", panel_id, "duplicate_start_id")
        expected_markers = (
            3,
            7 if stage in {CONFIRMATION_STAGE, GEOMETRY_STAGE} else 0,
            33 if stage == GEOMETRY_STAGE else 0,
        )
        actual_markers = tuple(
            int(panel.get(field) or 0)
            for field in (
                "screen_start_count",
                "confirmation_start_count",
                "geometry_start_count",
            )
        )
        if actual_markers != expected_markers:
            issue("panel", panel_id, "stale_or_inexact_3_7_33_lineage")
        expected_root_ids = {
            f"{panel_id}:B{index:02d}" for index in range(len(panel_roots))
        }
        actual_root_ids = {str(row["root_id"]) for row in panel_roots}
        if actual_root_ids != expected_root_ids:
            issue("panel", panel_id, "noncanonical_or_gapped_branch_ids")
        if int(panel.get("numerical_root_count") or 0) != len(panel_roots):
            issue("panel", panel_id, "numerical_root_count_cross_table_mismatch")
        wt_root_count = sum(_truth(row["passes_wt_gate"]) for row in panel_roots)
        if int(panel.get("wt_passing_root_count") or 0) != wt_root_count:
            issue("panel", panel_id, "wt_root_count_cross_table_mismatch")
        converged_attempt_count = sum(
            _truth(row["converged_root"]) for row in panel_attempts
        )
        if int(panel.get("converged_attempt_count") or 0) != converged_attempt_count:
            issue("panel", panel_id, "converged_attempt_count_cross_table_mismatch")
        if (
            panel_id in production_panel_ids
            and wt_root_count > 0
            and stage != GEOMETRY_STAGE
        ):
            issue("panel", panel_id, "passing_production_panel_not_current_33_start_geometry")

        attempts_by_start = {
            str(row.get("start_id", "")): row for row in panel_attempts
        }
        numerical_starts = {
            start_id
            for start_id, row in attempts_by_start.items()
            if _truth(row["passes_numerical_gate"])
        }
        numerical_endpoints: dict[str, np.ndarray] = {}
        for start_id in sorted(numerical_starts):
            attempt = attempts_by_start[start_id]
            try:
                endpoint = np.asarray(
                    json.loads(attempt["coordinates_and_other_json"]), dtype=float
                )
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                endpoint = np.asarray(())
            if endpoint.shape != (11,) or np.any(~np.isfinite(endpoint)):
                issue(
                    "attempt",
                    str(attempt.get("attempt_id", "")),
                    "invalid_endpoint",
                )
                continue
            numerical_endpoints[start_id] = endpoint
        expected_components = _connected_root_components(numerical_endpoints)
        claimed_starts: set[str] = set()
        actual_components: list[tuple[str, ...]] = []
        representatives: list[tuple[str, np.ndarray]] = []
        for root in panel_roots:
            root_id = str(root["root_id"])
            try:
                members = tuple(str(value) for value in json.loads(root["member_start_ids_json"]))
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                issue("root", root_id, "invalid_member_start_ids_json")
                continue
            if not members or len(members) != len(set(members)):
                issue("root", root_id, "empty_or_duplicate_member_start_ids")
            member_set = set(members)
            if claimed_starts & member_set:
                issue("root", root_id, "cluster_membership_overlaps_another_root")
            claimed_starts |= member_set
            actual_components.append(tuple(sorted(member_set)))
            if not member_set <= numerical_starts:
                issue("root", root_id, "cluster_member_not_numerical_attempt")
            endpoints = [
                numerical_endpoints[start_id]
                for start_id in members
                if start_id in numerical_endpoints
            ]
            try:
                coordinates = [
                    float(value) for value in json.loads(root["coordinates_json"])
                ]
                calibration = json.loads(root["calibration_json"])
                representative = np.asarray(
                    [
                        *coordinates,
                        float(calibration["cell_other_impermeant_osmoles_fmol"]),
                    ]
                )
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                representative = np.asarray(())
            if representative.shape != (11,) or np.any(~np.isfinite(representative)):
                issue("root", root_id, "invalid_root_representative")
            else:
                representatives.append((root_id, representative))
                if endpoints and not any(
                    np.array_equal(representative, endpoint) for endpoint in endpoints
                ):
                    issue("root", root_id, "representative_not_a_cluster_member_endpoint")
        if claimed_starts != numerical_starts:
            issue("panel", panel_id, "cluster_members_do_not_partition_numerical_attempts")
        if tuple(sorted(actual_components)) != tuple(sorted(expected_components)):
            issue(
                "panel",
                panel_id,
                "stored_clusters_do_not_match_deterministic_connected_components",
            )
        for index, (left_id, left) in enumerate(representatives):
            for right_id, right in representatives[index + 1 :]:
                separation = float(
                    np.max(np.abs((left - right) / ROOT_BOUNDS_SPAN))
                )
                if separation <= CLUSTER_RELATIVE_TOLERANCE:
                    issue(
                        "panel",
                        panel_id,
                        f"stored_root_representatives_not_separated:{left_id}:{right_id}",
                    )

    eligible_root_ids = sorted(
        str(root["root_id"])
        for root in root_rows
        if str(root.get("panel_id", "")) in production_panel_ids
    )
    eligible_wt_root_ids = sorted(
        str(root["root_id"])
        for root in root_rows
        if str(root.get("panel_id", "")) in production_panel_ids
        and _truth(root["passes_wt_gate"])
    )
    counts_by_kind = {
        kind: sum(
            1
            for panel_id in panels
            if panel_id in canonical and canonical[panel_id]["kind"] == kind
        )
        for kind in ("strict", "hw", "hwq", "absolute", "stress")
    }
    return {
        "status": "PASS" if not issues else "FAIL",
        "issue_count": len(issues),
        "issues": issues,
        "licensed_pump_contexts": list(licensed_pumps),
        "base_strict_current_wt_pass_present": base_strict_pass,
        "required_panel_count": len(required_ids),
        "observed_panel_count": len(panel_rows),
        "panel_counts_by_kind": counts_by_kind,
        "eligible_production_panel_ids": sorted(production_panel_ids),
        "eligible_production_root_ids": eligible_root_ids,
        "eligible_production_wt_root_ids": eligible_wt_root_ids,
        "accepted_legacy_diagnostic_start_design_count": len(
            accepted_legacy_diagnostic_panels
        ),
        "accepted_legacy_diagnostic_start_design_panel_ids": sorted(
            accepted_legacy_diagnostic_panels
        ),
    }


def reproduce_native_artifacts(
    root_table: Path,
    panel_table: Path,
    attempt_table: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Reproduce and audit one hash-bound snapshot of all native artifacts."""

    root_rows, root_digest = _csv_snapshot(root_table)
    panel_rows, panel_digest = _csv_snapshot(panel_table)
    attempt_rows, attempt_digest = _csv_snapshot(attempt_table)
    root_records, summary = _reproduce_native_root_rows(
        root_rows,
        root_table_name=root_table.name,
        root_table_sha256=root_digest,
    )
    attempt_records, attempt_summary = _reproduce_native_attempt_rows(
        panel_rows,
        attempt_rows,
        panel_table_sha256=panel_digest,
        attempt_table_sha256=attempt_digest,
    )
    hierarchy = audit_native_hierarchy(panel_rows, root_rows, attempt_rows)

    eligible_ids = set(hierarchy["eligible_production_root_ids"])
    eligible_records = [
        row for row in root_records if row["root_id"] in eligible_ids
    ]
    eligible_wt_records = [
        row for row in eligible_records if row["production_passes_wt_gate"]
    ]
    all_table_counts = {
        "row_count": summary["row_count"],
        "production_numerical_pass_count": summary[
            "production_numerical_pass_count"
        ],
        "production_wt_pass_count": summary["production_wt_pass_count"],
        "independent_pass_count": summary["independent_pass_count"],
    }
    summary["all_root_table_counts_before_production_filter"] = all_table_counts
    summary["production_numerical_pass_count"] = sum(
        bool(row["production_passes_numerical_gate"])
        for row in eligible_records
    )
    summary["production_wt_pass_count"] = len(eligible_wt_records)
    summary["independent_pass_count"] = sum(
        bool(row["independent_pass"]) for row in eligible_records
    )
    summary["independent_wt_pass_count"] = sum(
        bool(row["independent_wt_pass"]) for row in eligible_records
    )
    summary["eligible_production_root_ids"] = sorted(eligible_ids)
    summary["eligible_production_wt_root_ids"] = [
        row["root_id"] for row in eligible_wt_records
    ]

    artifact_hashes = {
        panel_table.name: panel_digest,
        root_table.name: root_digest,
        attempt_table.name: attempt_digest,
    }
    artifact_set_payload = json.dumps(
        artifact_hashes, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    summary["artifact_snapshot_sha256"] = artifact_hashes
    summary["artifact_snapshot_set_sha256"] = hashlib.sha256(
        artifact_set_payload
    ).hexdigest()
    summary["artifact_hash_basis"] = (
        "SHA256_OF_THE_EXACT_BYTES_PARSED_ONCE_FOR_EACH_CSV"
    )
    summary["panel_table_sha256"] = panel_digest
    summary["attempt_table_sha256"] = attempt_digest
    summary["optimizer_endpoint_reproduction"] = attempt_summary
    summary["hierarchy_audit"] = hierarchy

    nested_failures: list[str] = []
    if summary["gate_disagreement_count"]:
        nested_failures.append("root_numerical_gate_disagreement")
    if summary["wt_gate_disagreement_count"]:
        nested_failures.append("root_wt_gate_disagreement")
    if summary["observable_integrity_disagreement_count"]:
        nested_failures.append("root_observable_integrity_disagreement")
    if summary["duplicate_root_ids"]:
        nested_failures.append("duplicate_root_ids")
    if attempt_summary["convergence_gate_disagreement_count"]:
        nested_failures.append("attempt_convergence_gate_disagreement")
    if attempt_summary["root_cluster_count_disagreement_panel_ids"]:
        nested_failures.append("root_cluster_count_disagreement")
    if attempt_summary["panel_start_design_disagreement_count"]:
        nested_failures.append("panel_start_design_disagreement")
    if attempt_summary["duplicate_panel_ids"]:
        nested_failures.append("duplicate_panel_ids")
    if attempt_summary["duplicate_attempt_ids"]:
        nested_failures.append("duplicate_attempt_ids")
    if hierarchy["issue_count"]:
        nested_failures.append("cross_table_or_hierarchy_issue")
    summary["nested_gate_failure_count"] = len(nested_failures)
    summary["nested_gate_failures"] = nested_failures
    summary["status"] = (
        "INDEPENDENT_NATIVE_ROOT_REPRODUCTION_PASS"
        if not nested_failures
        else "INDEPENDENT_NATIVE_ROOT_REPRODUCTION_DISCREPANCY"
    )
    return root_records, attempt_records, summary


def _atomic_write_text(path: Path, value: str) -> None:
    """Fsync a sibling temporary file and atomically publish it."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to publish empty independent CSV {path.name}")
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    _atomic_write_text(path, buffer.getvalue())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    result_dir = args.repository / "results" / "13B_modern_full_model"
    root_table = result_dir / "native_source_wt_roots.csv"
    rows, attempt_rows, summary = reproduce_native_artifacts(
        root_table,
        result_dir / "native_source_panel.csv",
        result_dir / "native_source_root_attempts.csv",
    )
    csv_path = result_dir / "independent_native_roots.csv"
    attempt_csv_path = result_dir / "independent_native_attempts.csv"
    json_path = result_dir / "independent_native_summary.json"
    _write_csv(csv_path, rows)
    _write_csv(attempt_csv_path, attempt_rows)
    # The summary is the final commit marker: readers never see a new summary
    # pointing at partially published independent CSV reports.
    _atomic_write_text(
        json_path,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
