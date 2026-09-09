"""Run the staged WT-only native source-map and absolute-source root panel.

The runner never imports the held-out target ledger.  Stages are deliberately
separate so a completed strict screen is not rerun merely to add a later
diagnostic or absolute-source panel.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import asdict
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Sequence

from .native_source_panel import (
    AN_ABS_AE4_NHE,
    CONDITIONAL_PUMP_CAPACITY_SCALES,
    N_ABS_NKCC,
    STRICT_NATIVE,
    NativeRootSearchReport,
    NativeSourcePanelMember,
    declared_absolute_source_panel,
    declared_absolute_source_stress_panel,
    declared_hydraulic_diagnostic_panel,
    declared_native_panel,
    native_source_map,
    root_to_csv_row,
    solve_native_member,
)


ROOT_FIELDS = (
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
    "core_state_json",
    "coordinates_json",
    "calibration_json",
    "observables_json",
    "member_start_ids_json",
    "max_abs_scaled_residual",
    "max_abs_amount_rhs_fmol_s",
    "max_abs_volume_rhs_pL_s",
    "max_abs_omitted_charge_rhs_fmol_equivalent_s",
    "max_abs_charge_residual_fmol",
    "max_abs_current_residual_A",
    "normalized_jacobian_rank",
    "normalized_jacobian_nullity",
    "boundary_hits_json",
    "gate_failures_json",
    "root_search_stage",
    "root_search_start_count",
    "panel_confirmation_start_count",
    "panel_geometry_start_count",
)

COMMIT_LEDGER_NAME = "native_source_hashes.csv"
COMMIT_ARTIFACT_NAMES = (
    "native_source_panel.csv",
    "native_source_root_attempts.csv",
    "native_source_wt_roots.csv",
)
CONFIRMATION_STAGE = "SEVEN_START_CONFIRMATION"
GEOMETRY_STAGE = "THIRTY_THREE_START_ROOT_GEOMETRY"
LEGACY_HW_DIAGNOSTIC_ELIGIBILITY = "DIAGNOSTIC_ONLY_HYDRAULIC_TRANSFER"
CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY = (
    "DIAGNOSTIC_ONLY_HW_HYDRAULIC_TRANSFER"
)
LEGACY_HW_ELIGIBILITY_MIGRATION = (
    "LEGACY_DIAGNOSTIC_HW_ELIGIBILITY_TO_MODE_SPECIFIC"
)


def _expected_start_count(stage: str) -> int:
    if stage == CONFIRMATION_STAGE:
        return 7
    if stage == GEOMETRY_STAGE:
        return 33
    return 3


def _panel_seed(seed: int, panel_id: str) -> int:
    """Derive a stable per-panel seed independent of batch order."""

    digest = hashlib.sha256(f"{int(seed)}\0{panel_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


def _solve_job(
    job: tuple[NativeSourcePanelMember, int, int, int, tuple[tuple[float, ...], ...]]
) -> NativeRootSearchReport:
    member, start_count, seed, max_nfev, anchors = job
    return solve_native_member(
        member,
        start_count=start_count,
        seed=seed,
        max_nfev=max_nfev,
        anchors=anchors,
    )


def run_panel(
    members: Sequence[NativeSourcePanelMember],
    *,
    start_count: int,
    seed: int,
    max_nfev: int,
    workers: int,
    anchors_by_panel: Mapping[str, Sequence[Sequence[float]]] | None = None,
) -> tuple[NativeRootSearchReport, ...]:
    """Run fixed slices in parallel while preserving declared member order."""

    indexed = {member.panel_id: index for index, member in enumerate(members)}
    if len(indexed) != len(members):
        raise ValueError("run_panel rejects duplicate member/panel IDs")
    jobs = [
        (
            member,
            start_count,
            _panel_seed(seed, member.panel_id),
            max_nfev,
            tuple(
                tuple(float(value) for value in anchor)
                for anchor in (anchors_by_panel or {}).get(member.panel_id, ())
            ),
        )
        for index, member in enumerate(members)
    ]
    reports: list[NativeRootSearchReport] = []
    if workers == 1:
        reports = [_solve_job(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_solve_job, job): job[0].panel_id for job in jobs}
            for future in as_completed(futures):
                reports.append(future.result())
    return tuple(sorted(reports, key=lambda report: indexed[report.member.panel_id]))


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _merge_rows(
    old: Iterable[Mapping[str, Any]],
    new: Iterable[Mapping[str, Any]],
    *,
    key: str,
) -> list[Mapping[str, Any]]:
    merged: dict[str, Mapping[str, Any]] = {str(row[key]): row for row in old}
    for row in new:
        merged[str(row[key])] = row
    return [merged[item] for item in sorted(merged)]


def _canonical_hydraulic_mode(row: Mapping[str, Any]) -> str:
    value = row.get("hydraulic_mode")
    if value:
        return str(value)
    return "H1" if float(row["hydraulic_scale"]) == 1.0 else "HW"


def _replace_root_rows_atomically(
    existing: Iterable[Mapping[str, Any]],
    incoming: Iterable[Mapping[str, Any]],
    *,
    incoming_panel_ids: set[str],
) -> list[Mapping[str, Any]]:
    """Replace complete cluster sets for rerun panels; never retain stale Bxx."""

    retained = [
        row for row in existing if str(row["panel_id"]) not in incoming_panel_ids
    ]
    return _merge_rows(retained, incoming, key="root_id")


def _replace_attempt_rows_atomically(
    existing: Iterable[Mapping[str, Any]],
    incoming: Iterable[Mapping[str, Any]],
    *,
    incoming_panel_ids: set[str],
) -> list[Mapping[str, Any]]:
    """Replace the complete attempt set whenever a panel is rerun."""

    retained = [
        row for row in existing if str(row["panel_id"]) not in incoming_panel_ids
    ]
    return _merge_rows(retained, incoming, key="attempt_id")


def _actual_start_design(
    panel_row: Mapping[str, Any], attempt_rows: Sequence[Mapping[str, Any]]
) -> Mapping[str, Any]:
    """Describe the actual persisted starts, including legacy truncation."""

    start_ids = sorted(str(row["start_id"]) for row in attempt_rows)
    names = set(start_ids)
    has_both_quartiles = {"lower_quartile", "upper_quartile"}.issubset(names)
    legacy_displacement = bool(
        "frozen_reference" in names
        and any(name.startswith("continuation_") for name in names)
        and "lower_quartile" in names
        and "upper_quartile" not in names
    )
    if legacy_displacement:
        kind = "LEGACY_PRIORITY_TRUNCATION_UPPER_QUARTILE_DISPLACED"
    elif len(start_ids) == 3 and has_both_quartiles:
        kind = "one_parent_or_frozen_anchor_plus_both_quartiles"
    elif len(start_ids) >= 7 and has_both_quartiles:
        kind = "frozen_both_quartiles_unique_continuations_seeded_uniform_fill"
    else:
        kind = "ACTUAL_PERSISTED_START_SET"
    return {
        "kind": kind,
        "count": len(start_ids),
        "start_ids": start_ids,
        "actual_persisted_attempt_set": True,
        "both_quartiles_present": has_both_quartiles,
        "legacy_priority_truncation": legacy_displacement,
        "panel_stage": panel_row.get("stage"),
    }


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    """Atomically replace one CSV after fully writing and syncing a sibling temp."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary_name).replace(path)
    finally:
        if temporary_name is not None:
            temporary = Path(temporary_name)
            if temporary.exists():
                temporary.unlink()


def _write_text_atomic(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary_name).replace(path)
    finally:
        if temporary_name is not None:
            temporary = Path(temporary_name)
            if temporary.exists():
                temporary.unlink()


def _invalidate_commit(output: Path) -> None:
    """Remove the commit marker before any scientific-table mutation."""

    ledger = output / COMMIT_LEDGER_NAME
    if ledger.exists():
        ledger.unlink()


def _attempt_row(
    report: NativeRootSearchReport, attempt: Any, *, stage: str
) -> Mapping[str, Any]:
    return {
        "attempt_id": f"{report.member.panel_id}|{attempt.start_id}",
        "panel_id": report.member.panel_id,
        "source_class": report.member.source_class,
        "source_scale": report.member.source_scale,
        "pump_capacity_scale": report.member.pump_capacity_scale,
        "hydraulic_scale": report.member.hydraulic_scale,
        "hydraulic_mode": report.member.hydraulic_mode,
        "eligibility": report.member.eligibility,
        "root_search_stage": stage,
        "root_search_start_count": len(report.attempts),
        "start_id": attempt.start_id,
        "optimizer_success": attempt.optimizer_success,
        "admissible_state": attempt.admissible_state,
        "converged_root": attempt.converged_root,
        "passes_numerical_gate": attempt.passes_numerical_gate,
        "passes_wt_gate": attempt.passes_wt_gate,
        "max_abs_scaled_residual": attempt.max_abs_scaled_residual,
        "max_abs_amount_rhs_fmol_s": attempt.max_abs_amount_rhs_fmol_s,
        "max_abs_volume_rhs_pL_s": attempt.max_abs_volume_rhs_pL_s,
        "max_abs_omitted_charge_rhs_fmol_equivalent_s": (
            attempt.max_abs_omitted_charge_rhs_fmol_equivalent_s
        ),
        "normalized_jacobian_rank": attempt.normalized_jacobian_rank,
        "normalized_jacobian_nullity": attempt.normalized_jacobian_nullity,
        "boundary_hits_json": json.dumps(attempt.boundary_hits),
        "wt_gate_failures_json": json.dumps(attempt.wt_gate_failures),
        "coordinates_and_other_json": json.dumps(attempt.coordinates_and_other),
        "cost": attempt.cost,
        "optimality": attempt.optimality,
        "nfev": attempt.nfev,
        "message": attempt.message,
    }


ATTEMPT_FIELDS = (
    "attempt_id",
    "panel_id",
    "source_class",
    "source_scale",
    "pump_capacity_scale",
    "hydraulic_scale",
    "hydraulic_mode",
    "eligibility",
    "root_search_stage",
    "root_search_start_count",
    "start_id",
    "optimizer_success",
    "admissible_state",
    "converged_root",
    "passes_numerical_gate",
    "passes_wt_gate",
    "max_abs_scaled_residual",
    "max_abs_amount_rhs_fmol_s",
    "max_abs_volume_rhs_pL_s",
    "max_abs_omitted_charge_rhs_fmol_equivalent_s",
    "normalized_jacobian_rank",
    "normalized_jacobian_nullity",
    "boundary_hits_json",
    "wt_gate_failures_json",
    "coordinates_and_other_json",
    "cost",
    "optimality",
    "nfev",
    "message",
)


def _panel_row(report: NativeRootSearchReport, *, stage: str) -> Mapping[str, Any]:
    expected_start_count = _expected_start_count(stage)
    if len(report.attempts) != expected_start_count:
        raise ValueError(
            f"{stage} requires exactly {expected_start_count} starts per panel; "
            f"{report.member.panel_id} supplied {len(report.attempts)}"
        )
    numerical = [root for root in report.roots if root.passes_numerical_gate]
    passing = [root for root in numerical if root.passes_wt_gate]
    best_attempt = min(
        report.attempts,
        key=lambda attempt: attempt.max_abs_scaled_residual,
    )
    failure_sets = sorted(
        {failure for attempt in report.attempts for failure in attempt.wt_gate_failures}
    )
    member = report.member
    return {
        "panel_id": member.panel_id,
        "parent_panel_id": member.parent_panel_id,
        "stage": stage,
        "topology_id": member.topology_id,
        "routing_id": member.routing_id,
        "source_class": member.source_class,
        "source_scale": member.source_scale,
        "pump_capacity_scale": member.pump_capacity_scale,
        "apical_pump_fraction": member.apical_pump_fraction,
        "apical_k_fraction": member.apical_k_fraction,
        "ae4_cation_fraction": member.ae4_cation_fraction,
        "hydraulic_scale": member.hydraulic_scale,
        "hydraulic_mode": member.hydraulic_mode,
        "eligibility": member.eligibility,
        "attempt_count": len(report.attempts),
        "converged_attempt_count": sum(item.converged_root for item in report.attempts),
        "numerical_root_count": len(numerical),
        "wt_passing_root_count": len(passing),
        "best_max_abs_scaled_residual": best_attempt.max_abs_scaled_residual,
        "best_attempt_id": best_attempt.start_id,
        "failure_union_json": json.dumps(failure_sets),
        "start_design_json": json.dumps(report.start_design, sort_keys=True),
        "screen_start_count": (
            len(report.attempts)
            if stage not in (CONFIRMATION_STAGE, GEOMETRY_STAGE)
            else 0
        ),
        "confirmation_start_count": (
            len(report.attempts) if stage == CONFIRMATION_STAGE else 0
        ),
        "geometry_start_count": (
            len(report.attempts) if stage == GEOMETRY_STAGE else 0
        ),
    }


PANEL_FIELDS = (
    "panel_id",
    "parent_panel_id",
    "stage",
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
    "eligibility_migration",
    "attempt_count",
    "converged_attempt_count",
    "numerical_root_count",
    "wt_passing_root_count",
    "best_max_abs_scaled_residual",
    "best_attempt_id",
    "failure_union_json",
    "start_design_json",
    "screen_start_count",
    "confirmation_start_count",
    "geometry_start_count",
)


def _parent_panel_id_from_row(row: Mapping[str, Any]) -> str:
    source_class = str(row["source_class"])
    topology = str(row["topology_id"])
    routing = str(row["routing_id"])
    pump = float(row["pump_capacity_scale"])
    hydraulic = float(row["hydraulic_scale"])
    strict_pump_one = f"{STRICT_NATIVE}_S1_{topology}_{routing}_P1_H1"
    strict_same_pump = f"{STRICT_NATIVE}_S1_{topology}_{routing}_P{pump:g}_H1"
    if source_class != STRICT_NATIVE or hydraulic != 1.0:
        return strict_same_pump
    if pump != 1.0:
        return strict_pump_one
    return "ROOT_GENERATION_NATIVE_SOURCE_MAP"


def _canonical_hw_diagnostic_member(
    row: Mapping[str, Any],
) -> NativeSourcePanelMember | None:
    """Return the exact declared HW diagnostic represented by ``row``."""

    try:
        if (
            str(row.get("source_class", "")) != STRICT_NATIVE
            or _canonical_hydraulic_mode(row) != "HW"
            or float(row.get("hydraulic_scale", 1.0)) == 1.0
        ):
            return None
        pump_capacity_scale = float(row["pump_capacity_scale"])
    except (KeyError, TypeError, ValueError):
        return None
    expected_by_id = {
        member.panel_id: member
        for member in declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=(pump_capacity_scale,),
            hydraulic_modes=("HW",),
        )
    }
    member = expected_by_id.get(str(row.get("panel_id", "")))
    if member is None:
        return None
    candidate = {
        **row,
        "parent_panel_id": (
            row.get("parent_panel_id") or _parent_panel_id_from_row(row)
        ),
        "hydraulic_mode": "HW",
        "eligibility": CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY,
    }
    if _panel_row_semantic_mismatches(candidate, member):
        return None
    return member


def _legacy_hw_eligibility_migration(row: Mapping[str, Any]) -> str | None:
    """Recognize only the exact historical label on a canonical HW row."""

    if str(row.get("eligibility", "")) != LEGACY_HW_DIAGNOSTIC_ELIGIBILITY:
        return None
    if _canonical_hw_diagnostic_member(row) is None:
        return None
    return LEGACY_HW_ELIGIBILITY_MIGRATION


def _canonicalize_current_generation(
    panel_rows: Sequence[Mapping[str, Any]],
    attempt_rows: Sequence[Mapping[str, Any]],
    root_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]], list[Mapping[str, Any]]]:
    """Migrate honest current-generation provenance without inventing history."""

    attempts_by_panel: dict[str, list[Mapping[str, Any]]] = {}
    roots_by_panel: dict[str, list[Mapping[str, Any]]] = {}
    for row in attempt_rows:
        attempts_by_panel.setdefault(str(row["panel_id"]), []).append(row)
    for row in root_rows:
        roots_by_panel.setdefault(str(row["panel_id"]), []).append(row)

    migrated_panels: list[Mapping[str, Any]] = []
    generation_by_panel: dict[str, tuple[str, int, int, int, int]] = {}
    migrated_hw_panel_ids: set[str] = set()
    for original in panel_rows:
        row = dict(original)
        panel_id = str(row["panel_id"])
        eligibility_migration = _legacy_hw_eligibility_migration(row)
        if eligibility_migration is not None:
            row["eligibility"] = CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
            row["eligibility_migration"] = eligibility_migration
            migrated_hw_panel_ids.add(panel_id)
        elif (
            row.get("eligibility_migration") == LEGACY_HW_ELIGIBILITY_MIGRATION
            and row.get("eligibility") == CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
            and _canonical_hw_diagnostic_member(row) is not None
        ):
            migrated_hw_panel_ids.add(panel_id)
        else:
            row["eligibility_migration"] = row.get("eligibility_migration") or ""
        current_attempts = attempts_by_panel.get(panel_id, [])
        count = len(current_attempts)
        prior_screen = int(row.get("screen_start_count") or 0)
        prior_confirmation = int(row.get("confirmation_start_count") or 0)
        prior_stage = str(row.get("stage") or "")
        if count == 3:
            current_stage = (
                prior_stage
                if prior_stage not in (CONFIRMATION_STAGE, GEOMETRY_STAGE, "")
                else "CANONICALIZED_CURRENT_THREE_START_SCREEN"
            )
            screen_count, confirmation_count, geometry_count = 3, 0, 0
        elif count == 7:
            current_stage = CONFIRMATION_STAGE
            screen_count = 3 if prior_screen == 3 else 0
            confirmation_count, geometry_count = 7, 0
        elif count == 33:
            current_stage = GEOMETRY_STAGE
            screen_count = 3 if prior_screen == 3 else 0
            confirmation_count = 7 if prior_confirmation == 7 else 0
            geometry_count = 33
        else:
            current_stage = f"INVALID_CURRENT_{count}_START_GENERATION"
            screen_count = confirmation_count = geometry_count = 0

        current_roots = roots_by_panel.get(panel_id, [])
        wt_count = 0
        for root in current_roots:
            try:
                wt_count += int(_bool(root.get("passes_wt_gate")))
            except ValueError:
                pass
        row.update(
            {
                "stage": current_stage,
                "attempt_count": count,
                "numerical_root_count": len(current_roots),
                "wt_passing_root_count": wt_count,
                "screen_start_count": screen_count,
                "confirmation_start_count": confirmation_count,
                "geometry_start_count": geometry_count,
                "start_design_json": json.dumps(
                    _actual_start_design(row, current_attempts), sort_keys=True
                ),
            }
        )
        generation_by_panel[panel_id] = (
            current_stage,
            count,
            confirmation_count,
            geometry_count,
            screen_count,
        )
        migrated_panels.append(row)

    migrated_attempts: list[Mapping[str, Any]] = []
    for original in attempt_rows:
        row = dict(original)
        if (
            str(row.get("panel_id", "")) in migrated_hw_panel_ids
            and row.get("eligibility") == LEGACY_HW_DIAGNOSTIC_ELIGIBILITY
        ):
            row["eligibility"] = CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
        generation = generation_by_panel.get(str(row["panel_id"]))
        if generation is None:
            row["root_search_stage"] = "ORPHAN_CURRENT_GENERATION"
            row["root_search_start_count"] = 0
        else:
            row["root_search_stage"] = generation[0]
            row["root_search_start_count"] = generation[1]
        migrated_attempts.append(row)

    migrated_roots: list[Mapping[str, Any]] = []
    for original in root_rows:
        row = dict(original)
        if (
            str(row.get("panel_id", "")) in migrated_hw_panel_ids
            and row.get("eligibility") == LEGACY_HW_DIAGNOSTIC_ELIGIBILITY
        ):
            row["eligibility"] = CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
        generation = generation_by_panel.get(str(row["panel_id"]))
        if generation is None:
            row.update(
                {
                    "root_search_stage": "ORPHAN_CURRENT_GENERATION",
                    "root_search_start_count": 0,
                    "panel_confirmation_start_count": 0,
                    "panel_geometry_start_count": 0,
                }
            )
        else:
            row.update(
                {
                    "root_search_stage": generation[0],
                    "root_search_start_count": generation[1],
                    "panel_confirmation_start_count": generation[2],
                    "panel_geometry_start_count": generation[3],
                }
            )
        migrated_roots.append(row)
    return migrated_panels, migrated_attempts, migrated_roots


def persist_reports(
    output: Path,
    reports: Sequence[NativeRootSearchReport],
    *,
    stage: str,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    _invalidate_commit(output)
    panel_path = output / "native_source_panel.csv"
    attempt_path = output / "native_source_root_attempts.csv"
    root_path = output / "native_source_wt_roots.csv"
    existing_panel_rows = _read_rows(panel_path)
    existing_panel_by_id = {row["panel_id"]: row for row in existing_panel_rows}
    incoming_panel_rows: list[Mapping[str, Any]] = []
    for report in reports:
        row = dict(_panel_row(report, stage=stage))
        old = existing_panel_by_id.get(report.member.panel_id, {})
        if stage == CONFIRMATION_STAGE:
            if int(old.get("screen_start_count") or 0) != 3:
                raise RuntimeError(
                    f"{report.member.panel_id}: seven-start confirmation requires "
                    "exact three-start screen provenance"
                )
            row["screen_start_count"] = 3
            row["confirmation_start_count"] = 7
            row["geometry_start_count"] = 0
        elif stage == GEOMETRY_STAGE:
            if int(old.get("screen_start_count") or 0) != 3 or int(
                old.get("confirmation_start_count") or 0
            ) != 7:
                raise RuntimeError(
                    f"{report.member.panel_id}: 33-start geometry requires exact "
                    "three-start and seven-start provenance"
                )
            row["screen_start_count"] = 3
            row["confirmation_start_count"] = 7
            row["geometry_start_count"] = 33
        incoming_panel_rows.append(row)
    panel_rows = _merge_rows(
        existing_panel_rows,
        incoming_panel_rows,
        key="panel_id",
    )
    panel_rows = [
        {
            **row,
            "parent_panel_id": (
                row.get("parent_panel_id") or _parent_panel_id_from_row(row)
            ),
            "hydraulic_mode": _canonical_hydraulic_mode(row),
            "screen_start_count": row.get("screen_start_count") or 0,
            "confirmation_start_count": row.get("confirmation_start_count") or 0,
            "geometry_start_count": row.get("geometry_start_count") or 0,
        }
        for row in panel_rows
    ]
    incoming_panel_ids = {report.member.panel_id for report in reports}
    attempt_rows = _replace_attempt_rows_atomically(
        _read_rows(attempt_path),
        [
            _attempt_row(report, attempt, stage=stage)
            for report in reports
            for attempt in report.attempts
        ],
        incoming_panel_ids=incoming_panel_ids,
    )
    attempt_rows = [
        {**row, "hydraulic_mode": _canonical_hydraulic_mode(row)}
        for row in attempt_rows
    ]
    attempts_by_panel: dict[str, list[Mapping[str, Any]]] = {}
    for row in attempt_rows:
        attempts_by_panel.setdefault(str(row["panel_id"]), []).append(row)
    panel_rows = [
        {
            **row,
            "start_design_json": json.dumps(
                _actual_start_design(
                    row, attempts_by_panel.get(str(row["panel_id"]), [])
                ),
                sort_keys=True,
            ),
        }
        for row in panel_rows
    ]
    panel_by_id = {row["panel_id"]: row for row in panel_rows}
    incoming_root_rows = []
    for report in reports:
        panel_row = panel_by_id[report.member.panel_id]
        for root in report.roots:
            incoming_root_rows.append(
                {
                    **root_to_csv_row(root),
                    "root_search_stage": stage,
                    "root_search_start_count": len(report.attempts),
                    "panel_confirmation_start_count": panel_row[
                        "confirmation_start_count"
                    ],
                    "panel_geometry_start_count": panel_row[
                        "geometry_start_count"
                    ],
                }
            )
    root_rows = _replace_root_rows_atomically(
        _read_rows(root_path),
        incoming_root_rows,
        incoming_panel_ids=incoming_panel_ids,
    )
    root_rows = [
        {
            **row,
            "hydraulic_mode": _canonical_hydraulic_mode(row),
        }
        for row in root_rows
    ]
    panel_rows, attempt_rows, root_rows = _canonicalize_current_generation(
        panel_rows, attempt_rows, root_rows
    )
    _write_rows(panel_path, panel_rows, PANEL_FIELDS)
    _write_rows(attempt_path, attempt_rows, ATTEMPT_FIELDS)
    _write_rows(root_path, root_rows, ROOT_FIELDS)


def _bool(value: Any) -> bool:
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(
        "native artifact booleans must be bool or canonical case-sensitive "
        f"True/False CSV text, found {value!r}"
    )


def _eligible_production(row: Mapping[str, Any]) -> bool:
    return str(row.get("eligibility", "")) == "PRODUCTION_PRE_REVEAL"


def _is_production_h1(row: Mapping[str, Any]) -> bool:
    mode = row.get("hydraulic_mode") or (
        "H1" if float(row["hydraulic_scale"]) == 1.0 else "HW"
    )
    return (
        _eligible_production(row)
        and float(row["hydraulic_scale"]) == 1.0
        and mode == "H1"
    )


def _confirmation_panel_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    licensed_pump_contexts: Sequence[float] | None = None,
) -> tuple[Mapping[str, Any], ...]:
    """Production-only seven-start intake; H diagnostics are never promoted."""

    licensed = (
        None
        if licensed_pump_contexts is None
        else {float(value) for value in licensed_pump_contexts}
    )
    return tuple(
        row
        for row in rows
        if int(row["wt_passing_root_count"]) > 0 and _is_production_h1(row)
        and (
            licensed is None
            or float(row["pump_capacity_scale"]) in licensed
        )
    )


def _duplicate_values(rows: Sequence[Mapping[str, Any]], key: str) -> tuple[str, ...]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for row in rows:
        value = str(row.get(key, ""))
        if not value or value in seen:
            duplicates.add(value)
        seen.add(value)
    return tuple(sorted(duplicates))


def _current_production_lineage_status(
    panel_rows: Sequence[Mapping[str, Any]],
    attempt_rows: Sequence[Mapping[str, Any]],
    root_rows: Sequence[Mapping[str, Any]],
    production_roots: Sequence[Mapping[str, Any]],
) -> Mapping[str, Any]:
    """Validate that every retained production root is the current 33-start result."""

    issues: list[str] = []

    def checked_bool(row: Mapping[str, Any], field: str, identity: str) -> bool:
        try:
            return _bool(row.get(field))
        except ValueError:
            issues.append(f"noncanonical_boolean:{identity}:{field}")
            return False
    for rows, key, label in (
        (panel_rows, "panel_id", "panel"),
        (attempt_rows, "attempt_id", "attempt"),
        (root_rows, "root_id", "root"),
    ):
        duplicates = _duplicate_values(rows, key)
        if duplicates:
            issues.append(f"duplicate_{label}_ids:{list(duplicates)}")
    if issues:
        return {"complete": False, "issues": issues, "production_panel_ids": []}

    panels = {str(row["panel_id"]): row for row in panel_rows}
    attempts_by_panel: dict[str, list[Mapping[str, Any]]] = {}
    roots_by_panel: dict[str, list[Mapping[str, Any]]] = {}
    for row in attempt_rows:
        attempts_by_panel.setdefault(str(row["panel_id"]), []).append(row)
    for row in root_rows:
        roots_by_panel.setdefault(str(row["panel_id"]), []).append(row)
    production_panel_ids = sorted({str(row["panel_id"]) for row in production_roots})

    root_identity_fields = (
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
    float_fields = _FLOAT_MEMBER_IDENTITY_FIELDS - {"panel_id"}
    for panel_id in production_panel_ids:
        panel = panels.get(panel_id)
        if panel is None:
            issues.append(f"production_panel_missing:{panel_id}")
            continue
        if str(panel.get("stage")) != GEOMETRY_STAGE:
            issues.append(f"production_panel_not_geometry_stage:{panel_id}")
        if (
            int(panel.get("screen_start_count") or 0) != 3
            or int(panel.get("confirmation_start_count") or 0) != 7
            or int(panel.get("geometry_start_count") or 0) != 33
            or int(panel.get("attempt_count") or 0) != 33
        ):
            issues.append(f"production_panel_bad_3_7_33_provenance:{panel_id}")

        panel_attempts = attempts_by_panel.get(panel_id, [])
        if len(panel_attempts) != 33:
            issues.append(f"production_panel_attempt_count_not_33:{panel_id}")
        start_ids = [str(row.get("start_id", "")) for row in panel_attempts]
        if len(set(start_ids)) != 33:
            issues.append(f"production_panel_start_ids_not_unique:{panel_id}")
        if not {"lower_quartile", "upper_quartile"}.issubset(start_ids):
            issues.append(f"production_panel_quartile_starts_missing:{panel_id}")
        for attempt in panel_attempts:
            if (
                str(attempt.get("root_search_stage")) != GEOMETRY_STAGE
                or int(attempt.get("root_search_start_count") or 0) != 33
                or not _is_production_h1(attempt)
            ):
                issues.append(
                    f"production_attempt_not_current_geometry:{attempt.get('attempt_id')}"
                )

        panel_roots = roots_by_panel.get(panel_id, [])
        wt_roots = [
            row
            for row in panel_roots
            if checked_bool(row, "passes_wt_gate", str(row.get("root_id")))
        ]
        if int(panel.get("numerical_root_count") or 0) != len(panel_roots):
            issues.append(f"production_panel_numerical_root_count_mismatch:{panel_id}")
        if int(panel.get("wt_passing_root_count") or 0) != len(wt_roots):
            issues.append(f"production_panel_wt_root_count_mismatch:{panel_id}")
        for root in panel_roots:
            if (
                str(root.get("root_search_stage")) != GEOMETRY_STAGE
                or int(root.get("root_search_start_count") or 0) != 33
                or int(root.get("panel_confirmation_start_count") or 0) != 7
                or int(root.get("panel_geometry_start_count") or 0) != 33
                or not _is_production_h1(root)
                or not checked_bool(
                    root, "passes_numerical_gate", str(root.get("root_id"))
                )
            ):
                issues.append(f"production_root_not_current_geometry:{root.get('root_id')}")
            for field in root_identity_fields:
                left = root.get(field)
                right = panel.get(field)
                try:
                    matches = (
                        float(left) == float(right)
                        if field in float_fields
                        else str(left) == str(right)
                    )
                except (TypeError, ValueError):
                    matches = False
                if not matches:
                    issues.append(
                        f"production_root_panel_mismatch:{root.get('root_id')}:{field}"
                    )
            member_ids = tuple(json.loads(str(root.get("member_start_ids_json", "[]"))))
            if not member_ids or not set(member_ids).issubset(start_ids):
                issues.append(f"production_root_start_membership_invalid:{root.get('root_id')}")

    return {
        "complete": bool(production_panel_ids) and not issues,
        "issues": sorted(set(issues)),
        "production_panel_ids": production_panel_ids,
    }


def write_summary(output: Path, *, completed_stage: str) -> Mapping[str, Any]:
    snapshots: dict[str, bytes] = {}
    snapshot_rows: dict[str, list[dict[str, str]]] = {}
    for name in COMMIT_ARTIFACT_NAMES:
        path = output / name
        payload = path.read_bytes() if path.exists() else b""
        snapshots[name] = payload
        snapshot_rows[name] = (
            list(csv.DictReader(io.StringIO(payload.decode("utf-8"), newline="")))
            if payload
            else []
        )
    panel_rows = snapshot_rows["native_source_panel.csv"]
    attempt_rows = snapshot_rows["native_source_root_attempts.csv"]
    root_rows = snapshot_rows["native_source_wt_roots.csv"]
    current_hashes = {
        name: hashlib.sha256(payload).hexdigest()
        for name, payload in snapshots.items()
        if payload
    }
    upstream_error: str | None = None
    try:
        upstream_status = _upstream_grid_status(output)
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        upstream_error = str(exc)
        upstream_status = {
            "complete": False,
            "licensed_pump_contexts": [],
            "error": upstream_error,
        }
    licensed_pumps = {
        float(value) for value in upstream_status.get("licensed_pump_contexts", [])
    }
    boolean_parse_issues: list[str] = []

    def checked_root_bool(row: Mapping[str, Any], field: str) -> bool:
        try:
            return _bool(row.get(field))
        except ValueError:
            boolean_parse_issues.append(
                f"noncanonical_boolean:{row.get('root_id', '')}:{field}"
            )
            return False
    strict_rows = [
        row
        for row in panel_rows
        if row["source_class"] == STRICT_NATIVE
        and float(row["source_scale"]) == 1.0
        and float(row["pump_capacity_scale"]) == 1.0
        and float(row["hydraulic_scale"]) == 1.0
        and _canonical_hydraulic_mode(row) == "H1"
    ]
    eligible_numerical_rows = [
        row
        for row in root_rows
        if _eligible_production(row)
        and checked_root_bool(row, "passes_numerical_gate")
        and _is_production_h1(row)
        and float(row["pump_capacity_scale"]) in licensed_pumps
    ]
    eligible_wt_rows = [
        row
        for row in root_rows
        if _eligible_production(row) and checked_root_bool(row, "passes_wt_gate")
    ]
    production_intake_issues: list[str] = []
    production_roots: list[Mapping[str, Any]] = []
    for row in eligible_wt_rows:
        root_id = str(row.get("root_id", ""))
        if not checked_root_bool(row, "passes_numerical_gate"):
            production_intake_issues.append(f"wt_without_numerical_gate:{root_id}")
        elif not _is_production_h1(row):
            production_intake_issues.append(f"production_root_not_h1:{root_id}")
        elif float(row["pump_capacity_scale"]) not in licensed_pumps:
            production_intake_issues.append(f"production_root_unlicensed_pump:{root_id}")
        else:
            production_roots.append(row)
    strict_passing = [
        row
        for row in production_roots
        if row["source_class"] == STRICT_NATIVE
        and float(row["pump_capacity_scale"]) == 1.0
    ]
    source_passing = [
        row
        for row in production_roots
        if row["source_class"] != STRICT_NATIVE
    ]
    diagnostic_passing = [
        row
        for row in root_rows
        if not _eligible_production(row) and checked_root_bool(row, "passes_wt_gate")
    ]
    panel_by_id = {str(row["panel_id"]): row for row in panel_rows}
    production_panel_ids = sorted({row["panel_id"] for row in production_roots})
    confirmation_complete = bool(production_panel_ids) and all(
        panel_id in panel_by_id
        and int(panel_by_id[panel_id].get("screen_start_count") or 0) == 3
        and int(panel_by_id[panel_id].get("confirmation_start_count") or 0) == 7
        for panel_id in production_panel_ids
    )
    lineage_status = _current_production_lineage_status(
        panel_rows, attempt_rows, root_rows, production_roots
    )
    geometry_complete = bool(lineage_status["complete"])
    commit_status = _commit_ledger_status(output)
    if commit_status.get("complete") is True:
        committed_hashes = {
            name: data["sha256"]
            for name, data in commit_status.get("fingerprints", {}).items()
        }
        if committed_hashes != current_hashes:
            commit_status = {
                **commit_status,
                "complete": False,
                "issues": [
                    *commit_status.get("issues", []),
                    "scientific_tables_changed_during_summary_snapshot",
                ],
            }
    root_path = output / "native_source_wt_roots.csv"
    root_sha256 = current_hashes.get(root_path.name)
    independent_path = output / "independent_native_summary.json"
    independent_payload: Mapping[str, Any] = {}
    if independent_path.exists():
        independent_payload = json.loads(independent_path.read_text(encoding="utf-8"))
    independent_endpoint = independent_payload.get("optimizer_endpoint_reproduction", {})
    if not isinstance(independent_endpoint, Mapping):
        independent_endpoint = {}
    independent_panel_hash = independent_payload.get(
        "panel_table_sha256", independent_endpoint.get("panel_table_sha256")
    )
    independent_attempt_hash = independent_payload.get(
        "attempt_table_sha256", independent_endpoint.get("attempt_table_sha256")
    )
    independent_hierarchy = independent_payload.get("hierarchy_audit", {})
    if not isinstance(independent_hierarchy, Mapping):
        independent_hierarchy = {}
    expected_artifact_snapshot = {
        "native_source_panel.csv": current_hashes.get("native_source_panel.csv"),
        "native_source_wt_roots.csv": current_hashes.get(
            "native_source_wt_roots.csv"
        ),
        "native_source_root_attempts.csv": current_hashes.get(
            "native_source_root_attempts.csv"
        ),
    }
    expected_snapshot_set_sha256 = hashlib.sha256(
        json.dumps(
            expected_artifact_snapshot, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    production_root_ids = sorted(str(row["root_id"]) for row in production_roots)
    expected_eligible_numerical_root_ids = sorted(
        str(row["root_id"]) for row in eligible_numerical_rows
    )
    nested_zero = bool(
        int(independent_endpoint.get("convergence_gate_disagreement_count", -1)) == 0
        and int(independent_endpoint.get("panel_start_design_disagreement_count", -1)) == 0
        and not independent_endpoint.get("root_cluster_count_disagreement_panel_ids", ["missing"])
        and not independent_endpoint.get("duplicate_panel_ids", ["missing"])
        and not independent_endpoint.get("duplicate_attempt_ids", ["missing"])
        and independent_endpoint.get(
            "all_declared_confirmation_panels_have_exactly_7_current_starts"
        )
        is True
        and independent_endpoint.get(
            "all_declared_geometry_panels_have_exactly_33_current_starts"
        )
        is True
        and independent_hierarchy.get("status") == "PASS"
        and int(independent_hierarchy.get("issue_count", -1)) == 0
        and not independent_hierarchy.get("issues", ["missing"])
        and sorted(independent_hierarchy.get("licensed_pump_contexts", []))
        == sorted(licensed_pumps)
        and sorted(independent_hierarchy.get("eligible_production_root_ids", []))
        == expected_eligible_numerical_root_ids
        and sorted(
            independent_hierarchy.get("eligible_production_wt_root_ids", [])
        )
        == production_root_ids
    )
    independent_current = bool(
        production_roots
        and independent_payload.get("status")
        == "INDEPENDENT_NATIVE_ROOT_REPRODUCTION_PASS"
        and independent_payload.get("root_table_sha256") == root_sha256
        and independent_panel_hash == current_hashes.get("native_source_panel.csv")
        and independent_attempt_hash
        == current_hashes.get("native_source_root_attempts.csv")
        and independent_payload.get("artifact_snapshot_sha256")
        == expected_artifact_snapshot
        and independent_payload.get("artifact_snapshot_set_sha256")
        == expected_snapshot_set_sha256
        and independent_payload.get("artifact_hash_basis")
        == "SHA256_OF_THE_EXACT_BYTES_PARSED_ONCE_FOR_EACH_CSV"
        and int(independent_payload.get("row_count", -1)) == len(root_rows)
        and int(independent_payload.get("production_numerical_pass_count", -1))
        == len(expected_eligible_numerical_root_ids)
        and int(independent_payload.get("independent_pass_count", -1))
        == len(expected_eligible_numerical_root_ids)
        and int(independent_payload.get("production_wt_pass_count", -1))
        == len(production_roots)
        and int(independent_payload.get("independent_wt_pass_count", -1))
        == len(production_roots)
        and sorted(independent_payload.get("eligible_production_root_ids", []))
        == expected_eligible_numerical_root_ids
        and sorted(independent_payload.get("eligible_production_wt_root_ids", []))
        == production_root_ids
        and int(independent_payload.get("gate_disagreement_count", -1)) == 0
        and int(independent_payload.get("wt_gate_disagreement_count", -1)) == 0
        and int(
            independent_payload.get("observable_integrity_disagreement_count", -1)
        )
        == 0
        and not independent_payload.get("duplicate_root_ids", ["missing"])
        and int(independent_payload.get("nested_gate_failure_count", -1)) == 0
        and not independent_payload.get("nested_gate_failures", ["missing"])
        and nested_zero
    )
    ready_for_dynamics = bool(
        production_roots
        and upstream_status.get("complete") is True
        and not upstream_error
        and not production_intake_issues
        and not boolean_parse_issues
        and confirmation_complete
        and geometry_complete
        and commit_status.get("complete") is True
        and independent_current
    )
    summary = {
        "schema_version": 2,
        "completed_stage": completed_stage,
        "firewall": "WT_ONLY_NO_AE4_NULL_SECRETION_ACCESS",
        "strict_h1_pump1_expected_slice_count": 15,
        "strict_h1_pump1_completed_slice_count": len(strict_rows),
        "strict_h1_pump1_wt_passing_root_count": len(strict_passing),
        "production_absolute_source_wt_passing_root_count": len(source_passing),
        "all_production_wt_passing_root_count": len(production_roots),
        "diagnostic_only_wt_passing_root_count": len(diagnostic_passing),
        "production_root_ids": [row["root_id"] for row in production_roots],
        "production_numerical_root_ids": expected_eligible_numerical_root_ids,
        "diagnostic_root_ids": [row["root_id"] for row in diagnostic_passing],
        "production_panel_ids": production_panel_ids,
        "production_intake_issues": sorted(production_intake_issues),
        "boolean_parse_issues": sorted(set(boolean_parse_issues)),
        "upstream_grid_status": upstream_status,
        "exact_upstream_grids_complete": upstream_status.get("complete") is True,
        "current_production_lineage": lineage_status,
        "native_table_commit": commit_status,
        "seven_start_confirmation_complete_for_all_production_roots": (
            confirmation_complete
        ),
        "production_confirmation_complete": confirmation_complete,
        "thirty_three_start_geometry_complete_for_all_production_roots": (
            geometry_complete
        ),
        "geometry_confirmation_complete": geometry_complete,
        "independent_reproduction_current_and_passing": independent_current,
        "independent_root_reproduction_complete": independent_current,
        "native_root_table_sha256": root_sha256,
        "production_roots_ready_for_dynamics": ready_for_dynamics,
        "native_pre_reveal_freeze_ready": ready_for_dynamics,
        "native_root_freeze_status": (
            "READY_FOR_NATIVE_WT_DYNAMIC_CONTRACT"
            if ready_for_dynamics
            else "NOT_READY_FOR_NATIVE_WT_DYNAMIC_CONTRACT"
        ),
        "eligible_dynamic_root_ids": (
            [row["root_id"] for row in production_roots]
            if ready_for_dynamics
            else []
        ),
        "conditional_pump_rule": (
            "pump scales 0.5 and 2 are licensed only after the complete strict "
            "H=1,pump=1 screen has zero acceptable roots"
        ),
        "hydraulic_rule": "H>1 is diagnostic-only and cannot rescue production",
        "source_rule": (
            "all passing structural/scalar rows are retained; no first-passer selection"
        ),
        "heldout_reveal_authorized": False,
    }
    summary_path = output / "native_source_wt_summary.json"
    _write_text_atomic(
        summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def write_source_map(output: Path) -> None:
    path = output / "native_source_map.json"
    _write_text_atomic(
        path, json.dumps(native_source_map(), indent=2, sort_keys=True) + "\n"
    )


def scale_one_identity_rows(output: Path) -> tuple[Mapping[str, Any], ...]:
    """Map analytic N/AN scale-one aliases to every licensed strict slice."""

    pumps = _licensed_pump_contexts(output)
    roots_by_panel: dict[str, list[Mapping[str, Any]]] = {}
    for row in _read_rows(output / "native_source_wt_roots.csv"):
        roots_by_panel.setdefault(str(row["panel_id"]), []).append(row)
    rows: list[Mapping[str, Any]] = []
    for pump in pumps:
        for source_class in (N_ABS_NKCC, AN_ABS_AE4_NHE):
            for alias in declared_native_panel(
                source_class=source_class,
                source_scale=1.0,
                pump_capacity_scale=pump,
            ):
                strict_panel_id = alias.parent_panel_id
                strict_roots = roots_by_panel.get(strict_panel_id, [])
                if not strict_roots:
                    rows.append(
                        {
                            "alias_panel_id": alias.panel_id,
                            "alias_root_id": "",
                            "source_class": source_class,
                            "source_scale": 1.0,
                            "pump_capacity_scale": pump,
                            "topology_id": alias.topology_id,
                            "routing_id": alias.routing_id,
                            "strict_panel_id": strict_panel_id,
                            "strict_root_id": "",
                            "strict_root_status": "NO_NUMERICAL_STRICT_ROOT",
                            "exact_identity": True,
                            "identity_basis": (
                                "scale one leaves NKCC1, AE4, NHE1, all membranes, "
                                "water, and thermodynamic laws bitwise identical"
                            ),
                        }
                    )
                    continue
                for strict_root in strict_roots:
                    branch_suffix = str(strict_root["root_id"]).split(":")[-1]
                    rows.append(
                        {
                            "alias_panel_id": alias.panel_id,
                            "alias_root_id": f"{alias.panel_id}:{branch_suffix}",
                            "source_class": source_class,
                            "source_scale": 1.0,
                            "pump_capacity_scale": pump,
                            "topology_id": alias.topology_id,
                            "routing_id": alias.routing_id,
                            "strict_panel_id": strict_panel_id,
                            "strict_root_id": strict_root["root_id"],
                            "strict_root_status": "EXACT_ALIAS_OF_NUMERICAL_STRICT_ROOT",
                            "exact_identity": True,
                            "identity_basis": (
                                "scale one leaves NKCC1, AE4, NHE1, all membranes, "
                                "water, and thermodynamic laws bitwise identical"
                            ),
                        }
                    )
    return tuple(rows)


def write_scale_one_identity_controls(output: Path) -> None:
    _write_rows(
        output / "native_source_scale_one_identity.csv",
        scale_one_identity_rows(output),
        (
            "alias_panel_id",
            "alias_root_id",
            "source_class",
            "source_scale",
            "pump_capacity_scale",
            "topology_id",
            "routing_id",
            "strict_panel_id",
            "strict_root_id",
            "strict_root_status",
            "exact_identity",
            "identity_basis",
        ),
    )


def write_hashes(output: Path) -> None:
    names = (
        "native_source_map.json",
        "native_source_panel.csv",
        "native_source_root_attempts.csv",
        "native_source_wt_roots.csv",
        "native_source_wt_summary.json",
        "native_source_scale_one_identity.csv",
    )
    rows = []
    for name in names:
        path = output / name
        if path.exists():
            payload = path.read_bytes()
            rows.append(
                {
                    "artifact": name,
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    "bytes": len(payload),
                }
            )
    _write_rows(
        output / COMMIT_LEDGER_NAME,
        rows,
        ("artifact", "sha256", "bytes"),
    )


def _commit_ledger_status(output: Path) -> Mapping[str, Any]:
    """Validate the final panel/attempt/root hash ledger as a commit marker."""

    ledger_path = output / COMMIT_LEDGER_NAME
    if not ledger_path.exists():
        return {"complete": False, "issues": ["commit_ledger_missing"]}
    ledger_payload = ledger_path.read_bytes()
    rows = list(
        csv.DictReader(io.StringIO(ledger_payload.decode("utf-8"), newline=""))
    )
    names = [str(row.get("artifact", "")) for row in rows]
    issues: list[str] = []
    if len(names) != len(set(names)):
        issues.append("commit_ledger_duplicate_artifact")
    by_name = {str(row.get("artifact", "")): row for row in rows}
    fingerprints: dict[str, Mapping[str, Any]] = {}
    for name in COMMIT_ARTIFACT_NAMES:
        path = output / name
        row = by_name.get(name)
        if row is None:
            issues.append(f"commit_ledger_missing:{name}")
            continue
        if not path.exists():
            issues.append(f"commit_artifact_missing:{name}")
            continue
        payload = path.read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        size = len(payload)
        fingerprints[name] = {"sha256": digest, "bytes": size}
        if str(row.get("sha256")) != digest or int(row.get("bytes") or -1) != size:
            issues.append(f"commit_artifact_stale:{name}")
    if not ledger_path.exists() or ledger_path.read_bytes() != ledger_payload:
        issues.append("commit_ledger_changed_during_validation")
    return {
        "complete": not issues,
        "issues": issues,
        "fingerprints": fingerprints,
    }


def _existing_strict_pass_count(output: Path) -> int:
    rows = _read_rows(output / "native_source_panel.csv")
    return sum(
        int(row["wt_passing_root_count"])
        for row in rows
        if row["source_class"] == STRICT_NATIVE
        and float(row["source_scale"]) == 1.0
        and float(row["pump_capacity_scale"]) == 1.0
        and float(row["hydraulic_scale"]) == 1.0
    )


_MEMBER_IDENTITY_FIELDS = (
    "panel_id",
    "parent_panel_id",
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
_FLOAT_MEMBER_IDENTITY_FIELDS = {
    "source_scale",
    "pump_capacity_scale",
    "apical_pump_fraction",
    "apical_k_fraction",
    "ae4_cation_fraction",
    "hydraulic_scale",
}


def _panel_row_semantic_mismatches(
    row: Mapping[str, Any], member: NativeSourcePanelMember
) -> tuple[str, ...]:
    expected = asdict(member)
    mismatches: list[str] = []
    for field in _MEMBER_IDENTITY_FIELDS:
        actual_value = row.get(field)
        expected_value = expected[field]
        if field in _FLOAT_MEMBER_IDENTITY_FIELDS:
            try:
                matches = float(actual_value) == float(expected_value)
            except (TypeError, ValueError):
                matches = False
        else:
            matches = str(actual_value or "") == str(expected_value)
        if not matches:
            mismatches.append(field)
    return tuple(mismatches)


def _unique_panel_rows(rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Mapping[str, Any]]:
    by_id: dict[str, Mapping[str, Any]] = {}
    duplicates: list[str] = []
    for row in rows:
        panel_id = str(row.get("panel_id", ""))
        if not panel_id:
            raise RuntimeError("native panel table contains a blank panel_id")
        if panel_id in by_id:
            duplicates.append(panel_id)
        by_id[panel_id] = row
    if duplicates:
        raise RuntimeError(
            f"native panel table contains duplicate panel IDs: {sorted(set(duplicates))}"
        )
    return by_id


def _require_exact_members(
    rows_by_id: Mapping[str, Mapping[str, Any]],
    members: Sequence[NativeSourcePanelMember],
    *,
    context: str,
) -> tuple[Mapping[str, Any], ...]:
    expected_by_id = {member.panel_id: member for member in members}
    missing = sorted(set(expected_by_id) - set(rows_by_id))
    if missing:
        raise RuntimeError(f"{context} is missing canonical panels: {missing}")
    selected: list[Mapping[str, Any]] = []
    semantic: dict[str, tuple[str, ...]] = {}
    for panel_id, member in expected_by_id.items():
        row = rows_by_id[panel_id]
        issues = _panel_row_semantic_mismatches(row, member)
        if issues:
            semantic[panel_id] = issues
        if int(row.get("screen_start_count") or 0) != 3:
            semantic[panel_id] = tuple(
                sorted({*(semantic.get(panel_id, ())), "screen_start_count"})
            )
        selected.append(row)
    if semantic:
        raise RuntimeError(f"{context} has noncanonical rows: {semantic}")
    return tuple(selected)


def _licensed_pump_contexts(output: Path) -> tuple[float, ...]:
    """Return every pump context licensed by the completed rest hierarchy."""

    rows = _read_rows(output / "native_source_panel.csv")
    rows_by_id = _unique_panel_rows(rows)
    strict_pump_one = _require_exact_members(
        rows_by_id,
        declared_native_panel(pump_capacity_scale=1.0),
        context="strict H=1,pump=1 rest screen",
    )
    if any(int(row["wt_passing_root_count"]) for row in strict_pump_one):
        return (1.0,)
    conditional_members = tuple(
        member
        for pump in (0.5, 2.0)
        for member in declared_native_panel(pump_capacity_scale=pump)
    )
    _require_exact_members(
        rows_by_id,
        conditional_members,
        context="conditional pump rest screen",
    )
    return CONDITIONAL_PUMP_CAPACITY_SCALES


def _upstream_grid_status(output: Path) -> Mapping[str, Any]:
    rows = _read_rows(output / "native_source_panel.csv")
    rows_by_id = _unique_panel_rows(rows)
    pumps = _licensed_pump_contexts(output)
    expected_members = {
        "strict": tuple(
            member
            for pump in pumps
            for member in declared_native_panel(pump_capacity_scale=pump)
        ),
        "hw": declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=pumps, hydraulic_modes=("HW",)
        ),
        "hwq": declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=pumps, hydraulic_modes=("HWQ",)
        ),
        "absolute": declared_absolute_source_panel(
            include_identity=False, pump_capacity_scales=pumps
        ),
    }
    expected_by_category = {
        name: {member.panel_id: member for member in members}
        for name, members in expected_members.items()
    }
    expected_ids = {
        name: set(items) for name, items in expected_by_category.items()
    }

    def category(row: Mapping[str, Any]) -> str | None:
        source_class = str(row.get("source_class", ""))
        mode = _canonical_hydraulic_mode(row)
        if source_class == STRICT_NATIVE:
            if mode == "H1":
                return "strict"
            if mode == "HW":
                return "hw"
            if mode == "HWQ":
                return "hwq"
        if source_class in (N_ABS_NKCC, AN_ABS_AE4_NHE) and _eligible_production(row):
            return "absolute"
        return None

    available_ids = {name: set() for name in expected_ids}
    for row in rows:
        name = category(row)
        if name is not None:
            available_ids[name].add(str(row["panel_id"]))
    expected = {name: len(items) for name, items in expected_ids.items()}
    observed = {name: len(items) for name, items in available_ids.items()}
    missing = {
        name: sorted(expected_ids[name] - available_ids[name])
        for name in expected_ids
    }
    unexpected = {
        name: sorted(available_ids[name] - expected_ids[name])
        for name in expected_ids
    }
    semantic_mismatches: dict[str, list[str]] = {}
    screen_provenance_mismatches: list[str] = []
    eligibility_migration_panel_ids: list[str] = []
    eligibility_migration_issues: list[str] = []
    for row in rows:
        marker = str(row.get("eligibility_migration", ""))
        if not marker:
            continue
        panel_id = str(row.get("panel_id", ""))
        if (
            marker != LEGACY_HW_ELIGIBILITY_MIGRATION
            or row.get("eligibility") != CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
            or _canonical_hw_diagnostic_member(row) is None
        ):
            eligibility_migration_issues.append(panel_id)
        else:
            eligibility_migration_panel_ids.append(panel_id)
    for items in expected_by_category.values():
        for panel_id, member in items.items():
            row = rows_by_id.get(panel_id)
            if row is None:
                continue
            mismatches = _panel_row_semantic_mismatches(row, member)
            if mismatches:
                semantic_mismatches[panel_id] = list(mismatches)
            if int(row.get("screen_start_count") or 0) != 3:
                screen_provenance_mismatches.append(panel_id)
    return {
        "licensed_pump_contexts": list(pumps),
        "expected": expected,
        "observed": observed,
        "missing_panel_ids": missing,
        "unexpected_panel_ids": unexpected,
        "semantic_mismatches": semantic_mismatches,
        "screen_provenance_mismatch_panel_ids": sorted(
            screen_provenance_mismatches
        ),
        "legacy_diagnostic_eligibility_migration": {
            "from": LEGACY_HW_DIAGNOSTIC_ELIGIBILITY,
            "to": CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY,
            "scope": "EXACT_DECLARED_STRICT_SOURCE_NON_H1_HW_DIAGNOSTICS_ONLY",
            "panel_count": len(eligibility_migration_panel_ids),
            "panel_ids": sorted(eligibility_migration_panel_ids),
            "invalid_provenance_panel_ids": sorted(
                eligibility_migration_issues
            ),
        },
        "complete": all(
            not missing[name] and not unexpected[name] for name in expected_ids
        )
        and not semantic_mismatches
        and not screen_provenance_mismatches
        and not eligibility_migration_issues,
    }


def _require_complete_upstream_grids(output: Path) -> Mapping[str, Any]:
    status = _upstream_grid_status(output)
    if not status["complete"]:
        missing_counts = {
            name: len(items)
            for name, items in status["missing_panel_ids"].items()
        }
        unexpected_counts = {
            name: len(items)
            for name, items in status["unexpected_panel_ids"].items()
        }
        raise RuntimeError(
            "native hierarchy is incomplete before confirmation/geometry/stress: "
            f"observed={status['observed']} expected={status['expected']} "
            f"missing_counts={missing_counts} "
            f"unexpected_counts={unexpected_counts} "
            f"semantic_mismatch_count={len(status['semantic_mismatches'])} "
            "screen_provenance_mismatch_count="
            f"{len(status['screen_provenance_mismatch_panel_ids'])}"
        )
    return status


def _prune_unlicensed_pump_contexts(output: Path) -> tuple[str, ...]:
    """Remove stale conditional contexts after the strict licensing decision changes."""

    licensed = set(_licensed_pump_contexts(output))
    panel_path = output / "native_source_panel.csv"
    attempt_path = output / "native_source_root_attempts.csv"
    root_path = output / "native_source_wt_roots.csv"
    panel_rows = _read_rows(panel_path)
    removed = {
        str(row["panel_id"])
        for row in panel_rows
        if float(row["pump_capacity_scale"]) not in licensed
    }
    if not removed:
        return ()
    _invalidate_commit(output)
    _write_rows(
        panel_path,
        [row for row in panel_rows if str(row["panel_id"]) not in removed],
        PANEL_FIELDS,
    )
    _write_rows(
        attempt_path,
        [
            row
            for row in _read_rows(attempt_path)
            if str(row["panel_id"]) not in removed
        ],
        ATTEMPT_FIELDS,
    )
    _write_rows(
        root_path,
        [
            row
            for row in _read_rows(root_path)
            if str(row["panel_id"]) not in removed
        ],
        ROOT_FIELDS,
    )
    return tuple(sorted(removed))


def _persisted_dynamic_failure(output: Path) -> bool:
    gate_path = output / "native_dynamic_contract_gate.json"
    if not gate_path.exists():
        return False
    payload = json.loads(gate_path.read_text(encoding="utf-8"))
    return payload.get("status") == "SCIENTIFIC_GATE_FAIL"


def _member_from_panel_row(row: Mapping[str, Any]) -> NativeSourcePanelMember:
    return NativeSourcePanelMember(
        panel_id=str(row["panel_id"]),
        topology_id=str(row["topology_id"]),
        routing_id=str(row["routing_id"]),
        source_class=str(row["source_class"]),
        source_scale=float(row["source_scale"]),
        pump_capacity_scale=float(row["pump_capacity_scale"]),
        apical_pump_fraction=float(row["apical_pump_fraction"]),
        apical_k_fraction=float(row["apical_k_fraction"]),
        ae4_cation_fraction=float(row["ae4_cation_fraction"]),
        hydraulic_scale=float(row["hydraulic_scale"]),
        hydraulic_mode=(
            str(row.get("hydraulic_mode"))
            if row.get("hydraulic_mode")
            else ("H1" if float(row["hydraulic_scale"]) == 1.0 else "HW")
        ),
        eligibility=str(row["eligibility"]),
        parent_panel_id=(
            str(row.get("parent_panel_id") or _parent_panel_id_from_row(row))
        ),
    )


def _existing_root_anchors(output: Path) -> Mapping[str, tuple[tuple[float, ...], ...]]:
    anchors: dict[str, list[tuple[float, ...]]] = {}
    for row in _read_rows(output / "native_source_wt_roots.csv"):
        coordinates = [float(value) for value in json.loads(row["coordinates_json"])]
        other = float(
            json.loads(row["calibration_json"])[
                "cell_other_impermeant_osmoles_fmol"
            ]
        )
        anchors.setdefault(row["panel_id"], []).append(tuple([*coordinates, other]))
    return {key: tuple(value) for key, value in anchors.items()}


def _parent_continuation_anchors(
    output: Path, members: Sequence[NativeSourcePanelMember]
) -> Mapping[str, tuple[tuple[float, ...], ...]]:
    existing = _existing_root_anchors(output)
    return {
        member.panel_id: existing.get(member.parent_panel_id, ())
        for member in members
    }


def execute_stage(
    stage: str,
    *,
    output: Path,
    workers: int,
    start_count: int,
    max_nfev: int,
    seed: int,
    strict_dynamic_geometry_failed: bool = False,
    production_source_dynamic_failed: bool = False,
    panel_ids: Sequence[str] = (),
) -> None:
    screen_stages = {"strict", "hydraulic", "hydraulic_hwq", "absolute", "stress"}
    if stage in screen_stages and start_count != 3:
        raise RuntimeError(f"{stage} requires exactly three starts per panel")
    write_source_map(output)
    if stage == "strict":
        strict = declared_native_panel()
        reports = run_panel(
            strict,
            start_count=start_count,
            seed=seed,
            max_nfev=max_nfev,
            workers=workers,
        )
        persist_reports(output, reports, stage="STRICT_H1_PUMP1")
        if not any(root.passes_wt_gate for report in reports for root in report.roots):
            conditional: list[NativeSourcePanelMember] = []
            for pump_scale in CONDITIONAL_PUMP_CAPACITY_SCALES:
                if pump_scale == 1.0:
                    continue
                conditional.extend(
                    declared_native_panel(pump_capacity_scale=pump_scale)
                )
            conditional_reports = run_panel(
                conditional,
                start_count=start_count,
                seed=seed + 20000,
                max_nfev=max_nfev,
                workers=workers,
            )
            persist_reports(
                output,
                conditional_reports,
                stage="CONDITIONAL_PUMP_AFTER_ZERO_STRICT_ROOTS",
            )
        _prune_unlicensed_pump_contexts(output)
    elif stage == "hydraulic":
        if not strict_dynamic_geometry_failed:
            raise RuntimeError(
                "hydraulic diagnostics require explicit confirmation that the strict "
                "H=1 WT dynamic/one-SMG geometry gate failed"
            )
        pump_contexts = _licensed_pump_contexts(output)
        members = declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=pump_contexts,
            hydraulic_modes=("HW",),
        )
        reports = run_panel(
            members,
            start_count=start_count,
            seed=seed + 30000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_parent_continuation_anchors(output, members),
        )
        persist_reports(output, reports, stage="DIAGNOSTIC_HYDRAULIC_REROOT")
    elif stage == "hydraulic_hwq":
        if not strict_dynamic_geometry_failed:
            raise RuntimeError(
                "HWQ diagnostics require explicit confirmation that the strict "
                "H=1 WT dynamic/one-SMG geometry gate failed"
            )
        pump_contexts = _licensed_pump_contexts(output)
        members = declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=pump_contexts,
            hydraulic_modes=("HWQ",),
        )
        reports = run_panel(
            members,
            start_count=start_count,
            seed=seed + 35000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_parent_continuation_anchors(output, members),
        )
        persist_reports(output, reports, stage="DIAGNOSTIC_HWQ_REROOT")
    elif stage == "absolute":
        if not strict_dynamic_geometry_failed:
            raise RuntimeError(
                "absolute-source panel requires explicit confirmation that the strict "
                "H=1 WT dynamic/one-SMG geometry gate failed"
            )
        pump_contexts = _licensed_pump_contexts(output)
        # Scale-one duplicates are analytic nesting controls and are tested in
        # unit tests; they are not rerun as 30 redundant numerical slices.
        members = declared_absolute_source_panel(
            include_identity=False,
            pump_capacity_scales=pump_contexts,
        )
        reports = run_panel(
            members,
            start_count=start_count,
            seed=seed + 40000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_parent_continuation_anchors(output, members),
        )
        persist_reports(output, reports, stage="PRODUCTION_ABSOLUTE_SOURCE_GRID")
    elif stage == "stress":
        _require_complete_upstream_grids(output)
        production_source_rest_passes = sum(
            int(row["wt_passing_root_count"])
            for row in _read_rows(output / "native_source_panel.csv")
            if row["source_class"] in (N_ABS_NKCC, AN_ABS_AE4_NHE)
            and _eligible_production(row)
            and float(row["hydraulic_scale"]) == 1.0
        )
        if production_source_rest_passes and not (
            production_source_dynamic_failed and _persisted_dynamic_failure(output)
        ):
            raise RuntimeError(
                "stress rows after a rest-passing production grid require both the "
                "explicit flag and a persisted SCIENTIFIC_GATE_FAIL dynamic artifact"
            )
        pump_contexts = _licensed_pump_contexts(output)
        members = declared_absolute_source_stress_panel(
            pump_capacity_scales=pump_contexts
        )
        reports = run_panel(
            members,
            start_count=start_count,
            seed=seed + 50000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_parent_continuation_anchors(output, members),
        )
        persist_reports(output, reports, stage="DIAGNOSTIC_ABSOLUTE_SOURCE_STRESS")
    elif stage == "confirm":
        upstream = _require_complete_upstream_grids(output)
        rows = _confirmation_panel_rows(
            _read_rows(output / "native_source_panel.csv"),
            licensed_pump_contexts=upstream["licensed_pump_contexts"],
        )
        if not rows:
            raise RuntimeError("no WT-passing panel is available for seven-start confirmation")
        members = tuple(_member_from_panel_row(row) for row in rows)
        reports = run_panel(
            members,
            start_count=7,
            seed=seed + 60000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_existing_root_anchors(output),
        )
        persist_reports(output, reports, stage="SEVEN_START_CONFIRMATION")
    elif stage == "geometry":
        upstream = _require_complete_upstream_grids(output)
        if not panel_ids:
            raise RuntimeError(
                "33-start geometry requires one or more explicit retained --panel-id values"
            )
        if len(set(panel_ids)) != len(panel_ids):
            raise RuntimeError("33-start geometry rejects duplicate --panel-id values")
        rows_by_id = {
            row["panel_id"]: row
            for row in _read_rows(output / "native_source_panel.csv")
        }
        missing = sorted(set(panel_ids) - set(rows_by_id))
        if missing:
            raise RuntimeError(f"unknown retained panel ids: {missing}")
        ineligible = sorted(
            item
            for item in panel_ids
            if not _is_production_h1(rows_by_id[item])
            or float(rows_by_id[item]["pump_capacity_scale"])
            not in {float(value) for value in upstream["licensed_pump_contexts"]}
            or int(rows_by_id[item]["wt_passing_root_count"]) <= 0
            or int(rows_by_id[item].get("confirmation_start_count") or 0) != 7
        )
        if ineligible:
            raise RuntimeError(
                f"33-start production geometry rejects diagnostic/H>1 rows: {ineligible}"
            )
        members = tuple(_member_from_panel_row(rows_by_id[item]) for item in panel_ids)
        reports = run_panel(
            members,
            start_count=33,
            seed=seed + 70000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_existing_root_anchors(output),
        )
        persist_reports(output, reports, stage="THIRTY_THREE_START_ROOT_GEOMETRY")
    elif stage == "readiness":
        # No scientific rerun.  Recompute the deterministic readiness marker
        # only after the independent agent has hashed and reproduced the final
        # panel-atomic 33-start root table.
        pass
    elif stage == "canonicalize":
        # Schema/metadata migration only.  This is used after a long-running
        # pre-schema batch finishes; it never recomputes a scientific row.
        persist_reports(output, (), stage="CANONICALIZE_EXISTING_ARTIFACTS")
    else:
        raise ValueError(f"unknown stage {stage!r}")
    write_scale_one_identity_controls(output)
    write_summary(output, completed_stage=stage)
    write_hashes(output)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=(
            "strict",
            "hydraulic",
            "hydraulic_hwq",
            "absolute",
            "stress",
            "confirm",
            "geometry",
            "readiness",
            "canonicalize",
        ),
        required=True,
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    parser.add_argument(
        "--panel-id",
        action="append",
        default=[],
        help="retained production panel id for the 33-start geometry stage; repeatable",
    )
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--start-count", type=int, default=3)
    parser.add_argument("--max-nfev", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=13601)
    parser.add_argument(
        "--strict-dynamic-geometry-failed",
        action="store_true",
        help="authorize H and production source stages after the frozen strict dynamic/geometry failure",
    )
    parser.add_argument(
        "--production-source-dynamic-failed",
        action="store_true",
        help="authorize diagnostic stress rows after the production source dynamic hierarchy fails",
    )
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    execute_stage(
        args.stage,
        output=args.output_directory,
        workers=args.workers,
        start_count=args.start_count,
        max_nfev=args.max_nfev,
        seed=args.seed,
        strict_dynamic_geometry_failed=args.strict_dynamic_geometry_failed,
        production_source_dynamic_failed=args.production_source_dynamic_failed,
        panel_ids=tuple(args.panel_id),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
