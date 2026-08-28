"""Isolated WT-only absolute-source by hydraulic reroot diagnostic.

This module deliberately writes no ``native_source_*`` artifact.  It crosses
the retained fourfold NKCC source hypothesis with the two non-production
hydraulic transfer modes at H=2.5, then applies an exact 3 -> 7 -> 33 start
hierarchy.  A commit marker is published last so a reader can never mistake a
partially replaced artifact set for a complete diagnostic contract.

No held-out or genotype target is imported or evaluated here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Callable, Mapping, Sequence

from .native_source_panel import (
    CONDITIONAL_PUMP_CAPACITY_SCALES,
    N_ABS_NKCC,
    NativeRootAttempt,
    NativeRootRecord,
    NativeRootSearchReport,
    NativeSourcePanelMember,
    build_native_source_model,
    declared_native_panel,
    load_native_source_roots,
    root_to_csv_row,
)
from .run_native_source_panel import run_panel


DIAGNOSTIC_ELIGIBILITY = (
    "DIAGNOSTIC_ONLY_ABSOLUTE_SOURCE_S4_H2P5_HYDRAULIC_CROSS"
)
DIAGNOSTIC_STATUS = (
    "DIAGNOSTIC_ONLY_ABSOLUTE_SOURCE_HYDRAULIC_HIERARCHY_COMPLETE"
)
SOURCE_SCALE = 4.0
HYDRAULIC_SCALE = 2.5
HYDRAULIC_MODES = ("HW", "HWQ")
PUMP_CONTEXTS = tuple(float(value) for value in CONDITIONAL_PUMP_CAPACITY_SCALES)

SCREEN_STAGE = "ABSOLUTE_HYDRAULIC_SCREEN_3"
CONFIRMATION_STAGE = "ABSOLUTE_HYDRAULIC_CONFIRMATION_7"
GEOMETRY_STAGE = "ABSOLUTE_HYDRAULIC_GEOMETRY_33"
STAGE_START_COUNTS = {
    SCREEN_STAGE: 3,
    CONFIRMATION_STAGE: 7,
    GEOMETRY_STAGE: 33,
}
STAGE_ROOT_TAGS = {
    SCREEN_STAGE: "S03",
    CONFIRMATION_STAGE: "C07",
    GEOMETRY_STAGE: "G33",
}

SOURCE_MAP_NAME = "absolute_hydraulic_diagnostic_source_map.json"
PANEL_NAME = "absolute_hydraulic_diagnostic_panel.csv"
ATTEMPT_NAME = "absolute_hydraulic_diagnostic_attempts.csv"
ROOT_NAME = "absolute_hydraulic_diagnostic_roots.csv"
SUMMARY_NAME = "absolute_hydraulic_diagnostic_summary.json"
HASH_NAME = "absolute_hydraulic_diagnostic_hashes.csv"
COMMIT_NAME = "absolute_hydraulic_diagnostic_commit.json"
DYNAMIC_FAILURE_GATE_NAME = "native_dynamic_contract_gate.json"
HASHED_ARTIFACT_NAMES = (
    SOURCE_MAP_NAME,
    PANEL_NAME,
    ATTEMPT_NAME,
    ROOT_NAME,
    SUMMARY_NAME,
)

NATIVE_INPUT_NAMES = (
    "native_source_map.json",
    "native_source_panel.csv",
    "native_source_root_attempts.csv",
    "native_source_wt_roots.csv",
    "native_source_wt_summary.json",
    "native_source_scale_one_identity.csv",
)

PANEL_FIELDS = (
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
    "parent_anchor_count",
    "parent_anchor_root_ids_json",
    "screen_attempt_count",
    "screen_numerical_root_count",
    "screen_wt_passing_root_count",
    "confirmation_required",
    "confirmation_attempt_count",
    "confirmation_numerical_root_count",
    "confirmation_wt_passing_root_count",
    "geometry_required",
    "geometry_attempt_count",
    "geometry_numerical_root_count",
    "geometry_wt_passing_root_count",
    "final_stage",
    "final_numerical_root_count",
    "final_wt_passing_root_count",
    "final_status",
    "final_root_ids_json",
    "final_wt_root_ids_json",
    "diagnostic_only",
    "production_promotion_authorized",
    "heldout_reveal_authorized",
)

ATTEMPT_FIELDS = (
    "attempt_id",
    "panel_id",
    "root_search_stage",
    "root_search_start_count",
    "start_id",
    "source_class",
    "source_scale",
    "pump_capacity_scale",
    "hydraulic_scale",
    "hydraulic_mode",
    "eligibility",
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
    "diagnostic_only",
    "production_promotion_authorized",
    "heldout_reveal_authorized",
)

ROOT_FIELDS = (
    "root_id",
    "source_root_id",
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
    "screen_start_count",
    "confirmation_start_count",
    "geometry_start_count",
    "diagnostic_only",
    "production_promotion_authorized",
    "heldout_reveal_authorized",
)


@dataclass(frozen=True)
class AbsoluteHydraulicDiagnosticRoot:
    """Distinct wrapper that cannot enter a production native-root loader."""

    record: NativeRootRecord
    root_search_stage: str = GEOMETRY_STAGE
    diagnostic_only: bool = True
    production_promotion_authorized: bool = False
    heldout_reveal_authorized: bool = False

    def __post_init__(self) -> None:
        if self.root_search_stage != GEOMETRY_STAGE or not self.diagnostic_only:
            raise ValueError("absolute-hydraulic root must remain a geometry diagnostic")
        if self.production_promotion_authorized or self.heldout_reveal_authorized:
            raise ValueError("diagnostic root cannot authorize production or reveal")
        _validate_diagnostic_record(self.record, require_wt_pass=True)

    @property
    def root_id(self) -> str:
        return self.record.root_id

    @property
    def panel_id(self) -> str:
        return self.record.panel_id


@dataclass(frozen=True)
class AbsoluteHydraulicDiagnosticContract:
    """Fail-closed handoff for a later WT-only dynamic diagnostic."""

    roots: tuple[AbsoluteHydraulicDiagnosticRoot, ...]
    summary: Mapping[str, Any]
    status: str = DIAGNOSTIC_STATUS
    diagnostic_only: bool = True
    production_promotion_authorized: bool = False
    heldout_reveal_authorized: bool = False

    def __post_init__(self) -> None:
        if self.status != DIAGNOSTIC_STATUS or not self.diagnostic_only:
            raise ValueError("absolute-hydraulic contract must remain diagnostic-only")
        if self.production_promotion_authorized or self.heldout_reveal_authorized:
            raise ValueError("diagnostic roots cannot authorize production or reveal")
        if any(not isinstance(root, AbsoluteHydraulicDiagnosticRoot) for root in self.roots):
            raise TypeError("diagnostic contract requires distinct diagnostic root wrappers")


PanelRunner = Callable[..., tuple[NativeRootSearchReport, ...]]


def _h1_parent_panel_id(member: NativeSourcePanelMember) -> str:
    return (
        f"{N_ABS_NKCC}_S4_{member.topology_id}_{member.routing_id}_"
        f"P{member.pump_capacity_scale:g}_H1"
    )


def declared_absolute_hydraulic_diagnostic_panel(
) -> tuple[NativeSourcePanelMember, ...]:
    """Return the exact 2 modes x 3 pumps x 5 topologies x 3 routings cross."""

    members: list[NativeSourcePanelMember] = []
    for mode in HYDRAULIC_MODES:
        for pump in PUMP_CONTEXTS:
            for member in declared_native_panel(
                source_class=N_ABS_NKCC,
                source_scale=SOURCE_SCALE,
                pump_capacity_scale=pump,
                hydraulic_scale=HYDRAULIC_SCALE,
                hydraulic_mode=mode,
                eligibility=DIAGNOSTIC_ELIGIBILITY,
            ):
                members.append(
                    replace(
                        member,
                        panel_id=(
                            "ABS_HYD_DIAG_NABS_S4_H2p5_"
                            f"{mode}_{member.topology_id}_{member.routing_id}_"
                            f"P{pump:g}"
                        ),
                        parent_panel_id=_h1_parent_panel_id(member),
                    )
                )
    result = tuple(members)
    if len(result) != 90 or len({member.panel_id for member in result}) != 90:
        raise AssertionError("absolute-hydraulic diagnostic declaration must be 90 unique panels")
    return result


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _strict_bool(value: Any, *, field: str) -> bool:
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(f"{field} must use a canonical boolean")


def _native_snapshot(output: Path) -> Mapping[str, Any]:
    """Hash and validate the immutable native WT input set once."""

    ledger_path = output / "native_source_hashes.csv"
    if not ledger_path.exists():
        raise RuntimeError("native source commit ledger is missing")
    ledger_payload = ledger_path.read_bytes()
    rows = list(csv.DictReader(io.StringIO(ledger_payload.decode("utf-8"))))
    by_name = {str(row.get("artifact", "")): row for row in rows}
    if len(by_name) != len(rows) or set(by_name) != set(NATIVE_INPUT_NAMES):
        raise RuntimeError("native source commit ledger has a non-exact artifact set")
    artifacts: dict[str, Mapping[str, Any]] = {}
    for name in NATIVE_INPUT_NAMES:
        path = output / name
        if not path.exists():
            raise RuntimeError(f"native source input is missing: {name}")
        payload = path.read_bytes()
        digest = _sha256(payload)
        size = len(payload)
        row = by_name[name]
        if str(row.get("sha256")) != digest or int(row.get("bytes") or -1) != size:
            raise RuntimeError(f"native source input hash is stale: {name}")
        artifacts[name] = {"sha256": digest, "bytes": size}
    summary = json.loads((output / "native_source_wt_summary.json").read_text())
    required_true = (
        "production_confirmation_complete",
        "geometry_confirmation_complete",
        "independent_root_reproduction_complete",
        "native_pre_reveal_freeze_ready",
    )
    if any(summary.get(name) is not True for name in required_true):
        raise RuntimeError("native source WT freeze is not ready for diagnostics")
    if summary.get("heldout_reveal_authorized") is not False:
        raise RuntimeError("absolute-hydraulic diagnostic requires a pre-reveal WT freeze")
    gate_path = output / DYNAMIC_FAILURE_GATE_NAME
    if not gate_path.exists():
        raise RuntimeError("persisted WT dynamic scientific-failure gate is missing")
    gate_payload = gate_path.read_bytes()
    gate = json.loads(gate_payload)
    if (
        gate.get("status") != "SCIENTIFIC_GATE_FAIL"
        or gate.get("firewall") != "WT_ONLY_NO_PHENOTYPE_TARGET_ACCESS"
        or gate.get("phenotype_reveal_permitted") is not False
        or gate.get("at_least_one_scientific_group_pass") is not False
        or int(gate.get("scientific_passing_group_count", -1)) != 0
        or int(gate.get("native_root_count", -1))
        != int(summary.get("all_production_wt_passing_root_count", -2))
    ):
        raise RuntimeError("persisted WT dynamic failure does not license this diagnostic")
    return {
        "ledger_sha256": _sha256(ledger_payload),
        "artifacts": artifacts,
        "native_root_table_sha256": summary.get("native_root_table_sha256"),
        "dynamic_failure_gate": {
            "artifact": DYNAMIC_FAILURE_GATE_NAME,
            "sha256": _sha256(gate_payload),
            "bytes": len(gate_payload),
            "status": gate.get("status"),
            "native_root_count": gate.get("native_root_count"),
        },
    }


def _parent_anchor_maps(
    native_output: Path,
    members: Sequence[NativeSourcePanelMember],
) -> tuple[
    Mapping[str, tuple[tuple[float, ...], ...]],
    Mapping[str, tuple[str, ...]],
]:
    roots_by_panel: dict[str, list[NativeRootRecord]] = {}
    for root in load_native_source_roots(native_output):
        roots_by_panel.setdefault(root.panel_id, []).append(root)
    anchors: dict[str, tuple[tuple[float, ...], ...]] = {}
    root_ids: dict[str, tuple[str, ...]] = {}
    for member in members:
        roots = sorted(
            roots_by_panel.get(member.parent_panel_id, ()),
            key=lambda root: root.root_id,
        )
        anchors[member.panel_id] = tuple(
            (*root.coordinates, root.other_impermeant_osmoles_fmol)
            for root in roots
        )
        root_ids[member.panel_id] = tuple(root.root_id for root in roots)
    return anchors, root_ids


def _report_anchor_map(
    reports: Sequence[NativeRootSearchReport],
) -> Mapping[str, tuple[tuple[float, ...], ...]]:
    return {
        report.member.panel_id: tuple(
            (*root.coordinates, root.other_impermeant_osmoles_fmol)
            for root in sorted(report.roots, key=lambda item: item.root_id)
        )
        for report in reports
    }


def _passing_members(
    reports: Sequence[NativeRootSearchReport],
) -> tuple[NativeSourcePanelMember, ...]:
    return tuple(
        report.member
        for report in reports
        if any(root.passes_numerical_gate and root.passes_wt_gate for root in report.roots)
    )


def _validate_reports(
    reports: Sequence[NativeRootSearchReport],
    members: Sequence[NativeSourcePanelMember],
    *,
    stage: str,
) -> None:
    expected_count = STAGE_START_COUNTS[stage]
    if [report.member.panel_id for report in reports] != [
        member.panel_id for member in members
    ]:
        raise RuntimeError(f"{stage}: report/member panel order is not exact")
    for report in reports:
        member = report.member
        if member.eligibility != DIAGNOSTIC_ELIGIBILITY:
            raise RuntimeError(f"{member.panel_id}: production eligibility is forbidden")
        start_ids = [attempt.start_id for attempt in report.attempts]
        if len(start_ids) != expected_count or len(set(start_ids)) != expected_count:
            raise RuntimeError(
                f"{member.panel_id}: {stage} requires {expected_count} unique starts"
            )
        if int(report.start_design.get("count", -1)) != expected_count:
            raise RuntimeError(f"{member.panel_id}: start-design count is inconsistent")
        if set(report.start_design.get("start_ids", ())) != set(start_ids):
            raise RuntimeError(f"{member.panel_id}: start-design IDs are inconsistent")
        for root in report.roots:
            _validate_diagnostic_record(root, require_wt_pass=False)
            if root.panel_id != member.panel_id:
                raise RuntimeError(f"{member.panel_id}: cross-panel root detected")
            if not set(root.member_start_ids).issubset(start_ids):
                raise RuntimeError(f"{root.root_id}: cluster membership is not a start subset")


def _attempt_row(
    report: NativeRootSearchReport,
    attempt: NativeRootAttempt,
    *,
    stage: str,
) -> Mapping[str, Any]:
    return {
        "attempt_id": f"{stage}|{report.member.panel_id}|{attempt.start_id}",
        "panel_id": report.member.panel_id,
        "root_search_stage": stage,
        "root_search_start_count": len(report.attempts),
        "start_id": attempt.start_id,
        "source_class": report.member.source_class,
        "source_scale": report.member.source_scale,
        "pump_capacity_scale": report.member.pump_capacity_scale,
        "hydraulic_scale": report.member.hydraulic_scale,
        "hydraulic_mode": report.member.hydraulic_mode,
        "eligibility": report.member.eligibility,
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
        "diagnostic_only": True,
        "production_promotion_authorized": False,
        "heldout_reveal_authorized": False,
    }


def _stage_root_id(root: NativeRootRecord, stage: str) -> str:
    branch = root.root_id.rsplit(":", 1)[-1]
    return f"{root.panel_id}:{STAGE_ROOT_TAGS[stage]}:{branch}"


def _root_row(
    report: NativeRootSearchReport,
    root: NativeRootRecord,
    *,
    stage: str,
) -> Mapping[str, Any]:
    row = dict(root_to_csv_row(root))
    row.update(
        {
            "root_id": _stage_root_id(root, stage),
            "source_root_id": root.root_id,
            "parent_panel_id": report.member.parent_panel_id,
            "root_search_stage": stage,
            "root_search_start_count": len(report.attempts),
            "screen_start_count": 3,
            "confirmation_start_count": (
                7 if stage in (CONFIRMATION_STAGE, GEOMETRY_STAGE) else 0
            ),
            "geometry_start_count": 33 if stage == GEOMETRY_STAGE else 0,
            "diagnostic_only": True,
            "production_promotion_authorized": False,
            "heldout_reveal_authorized": False,
        }
    )
    return row


def _report_counts(report: NativeRootSearchReport | None) -> tuple[int, int, int]:
    if report is None:
        return (0, 0, 0)
    numerical = [root for root in report.roots if root.passes_numerical_gate]
    passing = [root for root in numerical if root.passes_wt_gate]
    return (len(report.attempts), len(numerical), len(passing))


def _panel_rows(
    members: Sequence[NativeSourcePanelMember],
    screen: Sequence[NativeRootSearchReport],
    confirmation: Sequence[NativeRootSearchReport],
    geometry: Sequence[NativeRootSearchReport],
    *,
    parent_anchor_ids: Mapping[str, Sequence[str]],
) -> list[Mapping[str, Any]]:
    screen_by_id = {report.member.panel_id: report for report in screen}
    confirmation_by_id = {report.member.panel_id: report for report in confirmation}
    geometry_by_id = {report.member.panel_id: report for report in geometry}
    rows: list[Mapping[str, Any]] = []
    for member in members:
        screen_report = screen_by_id[member.panel_id]
        confirmation_report = confirmation_by_id.get(member.panel_id)
        geometry_report = geometry_by_id.get(member.panel_id)
        screen_counts = _report_counts(screen_report)
        confirmation_counts = _report_counts(confirmation_report)
        geometry_counts = _report_counts(geometry_report)
        screen_pass = screen_counts[2] > 0
        confirmation_pass = confirmation_counts[2] > 0
        geometry_pass = geometry_counts[2] > 0
        final_report = geometry_report or confirmation_report or screen_report
        if geometry_pass:
            final_status = "WT_PASS_GEOMETRY_CONFIRMED_DIAGNOSTIC_ONLY"
        elif confirmation_pass:
            final_status = "WT_PASS_LOST_AT_GEOMETRY"
        elif screen_pass:
            final_status = "WT_PASS_LOST_AT_CONFIRMATION"
        else:
            final_status = "WT_FAIL_AT_SCREEN"
        final_stage = (
            GEOMETRY_STAGE
            if geometry_report is not None
            else CONFIRMATION_STAGE
            if confirmation_report is not None
            else SCREEN_STAGE
        )
        final_ids = [_stage_root_id(root, final_stage) for root in final_report.roots]
        final_wt_ids = [
            _stage_root_id(root, final_stage)
            for root in final_report.roots
            if root.passes_numerical_gate and root.passes_wt_gate
        ]
        rows.append(
            {
                **asdict(member),
                "parent_anchor_count": len(parent_anchor_ids.get(member.panel_id, ())),
                "parent_anchor_root_ids_json": json.dumps(
                    tuple(parent_anchor_ids.get(member.panel_id, ()))
                ),
                "screen_attempt_count": screen_counts[0],
                "screen_numerical_root_count": screen_counts[1],
                "screen_wt_passing_root_count": screen_counts[2],
                "confirmation_required": screen_pass,
                "confirmation_attempt_count": confirmation_counts[0],
                "confirmation_numerical_root_count": confirmation_counts[1],
                "confirmation_wt_passing_root_count": confirmation_counts[2],
                "geometry_required": confirmation_pass,
                "geometry_attempt_count": geometry_counts[0],
                "geometry_numerical_root_count": geometry_counts[1],
                "geometry_wt_passing_root_count": geometry_counts[2],
                "final_stage": final_stage,
                "final_numerical_root_count": _report_counts(final_report)[1],
                "final_wt_passing_root_count": _report_counts(final_report)[2],
                "final_status": final_status,
                "final_root_ids_json": json.dumps(final_ids),
                "final_wt_root_ids_json": json.dumps(final_wt_ids),
                "diagnostic_only": True,
                "production_promotion_authorized": False,
                "heldout_reveal_authorized": False,
            }
        )
    return rows


def _csv_payload(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode("utf-8")


def _json_payload(value: Mapping[str, Any]) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
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
        if temporary_name is not None and Path(temporary_name).exists():
            Path(temporary_name).unlink()


def _publish(
    output: Path,
    payloads: Mapping[str, bytes],
) -> None:
    if set(payloads) != set(HASHED_ARTIFACT_NAMES):
        raise ValueError("diagnostic publication has a non-exact artifact set")
    output.mkdir(parents=True, exist_ok=True)
    commit_path = output / COMMIT_NAME
    if commit_path.exists():
        commit_path.unlink()
    for name in HASHED_ARTIFACT_NAMES:
        _atomic_write(output / name, payloads[name])
    hash_rows = [
        {
            "artifact": name,
            "sha256": _sha256(payloads[name]),
            "bytes": len(payloads[name]),
        }
        for name in HASHED_ARTIFACT_NAMES
    ]
    hash_payload = _csv_payload(hash_rows, ("artifact", "sha256", "bytes"))
    _atomic_write(output / HASH_NAME, hash_payload)
    commit = {
        "schema_version": 1,
        "status": DIAGNOSTIC_STATUS,
        "diagnostic_only": True,
        "production_promotion_authorized": False,
        "heldout_reveal_authorized": False,
        "hash_ledger": HASH_NAME,
        "hash_ledger_sha256": _sha256(hash_payload),
        "hashed_artifacts": list(HASHED_ARTIFACT_NAMES),
    }
    _atomic_write(commit_path, _json_payload(commit))


def _validate_diagnostic_record(
    root: NativeRootRecord,
    *,
    require_wt_pass: bool,
) -> None:
    if (
        root.source_class != N_ABS_NKCC
        or root.source_scale != SOURCE_SCALE
        or root.pump_capacity_scale not in PUMP_CONTEXTS
        or root.hydraulic_scale != HYDRAULIC_SCALE
        or root.hydraulic_mode not in HYDRAULIC_MODES
        or root.eligibility != DIAGNOSTIC_ELIGIBILITY
    ):
        raise ValueError(f"{root.root_id}: noncanonical diagnostic root")
    if not root.passes_numerical_gate:
        raise ValueError(f"{root.root_id}: root did not pass the numerical gate")
    if require_wt_pass and not root.passes_wt_gate:
        raise ValueError(f"{root.root_id}: dynamic diagnostic requires a WT-pass root")


def diagnostic_member_from_root(
    root: AbsoluteHydraulicDiagnosticRoot,
) -> NativeSourcePanelMember:
    """Reconstruct the immutable diagnostic member without promoting it."""

    if not isinstance(root, AbsoluteHydraulicDiagnosticRoot):
        raise TypeError("a distinct absolute-hydraulic diagnostic root is required")
    record = root.record
    _validate_diagnostic_record(record, require_wt_pass=True)
    member = NativeSourcePanelMember(
        panel_id=record.panel_id,
        topology_id=record.topology_id,
        routing_id=record.routing_id,
        source_class=record.source_class,
        source_scale=record.source_scale,
        pump_capacity_scale=record.pump_capacity_scale,
        apical_pump_fraction=record.apical_pump_fraction,
        apical_k_fraction=record.apical_k_fraction,
        ae4_cation_fraction=record.ae4_cation_fraction,
        hydraulic_scale=record.hydraulic_scale,
        hydraulic_mode=record.hydraulic_mode,
        eligibility=DIAGNOSTIC_ELIGIBILITY,
        parent_panel_id=_h1_parent_panel_id(record),
    )
    expected = {item.panel_id: item for item in declared_absolute_hydraulic_diagnostic_panel()}
    if member.panel_id not in expected or asdict(member) != asdict(expected[member.panel_id]):
        raise ValueError(f"{record.root_id}: root/member semantic mismatch")
    return member


def build_absolute_hydraulic_diagnostic_model(
    root: AbsoluteHydraulicDiagnosticRoot,
    *,
    stimulus: Any = None,
    regulatory_model: Any = None,
):
    """Build a WT diagnostic model; this function grants no production status."""

    member = diagnostic_member_from_root(root)
    return build_native_source_model(
        member,
        other_impermeant_osmoles_fmol=root.record.other_impermeant_osmoles_fmol,
        stimulus=stimulus,
        regulatory_model=regulatory_model,
    )


def _artifact_payloads(
    members: Sequence[NativeSourcePanelMember],
    screen: Sequence[NativeRootSearchReport],
    confirmation: Sequence[NativeRootSearchReport],
    geometry: Sequence[NativeRootSearchReport],
    *,
    parent_anchor_ids: Mapping[str, Sequence[str]],
    native_snapshot: Mapping[str, Any],
    seed: int,
    max_nfev: int,
) -> Mapping[str, bytes]:
    stage_reports = (
        (SCREEN_STAGE, screen),
        (CONFIRMATION_STAGE, confirmation),
        (GEOMETRY_STAGE, geometry),
    )
    panel_rows = _panel_rows(
        members,
        screen,
        confirmation,
        geometry,
        parent_anchor_ids=parent_anchor_ids,
    )
    attempt_rows = [
        _attempt_row(report, attempt, stage=stage)
        for stage, reports in stage_reports
        for report in reports
        for attempt in report.attempts
    ]
    root_rows = [
        _root_row(report, root, stage=stage)
        for stage, reports in stage_reports
        for report in reports
        for root in report.roots
    ]
    screen_ids = [member.panel_id for member in _passing_members(screen)]
    confirmation_ids = [member.panel_id for member in _passing_members(confirmation)]
    geometry_ids = [member.panel_id for member in _passing_members(geometry)]
    geometry_root_ids = [
        _stage_root_id(root, GEOMETRY_STAGE)
        for report in geometry
        for root in report.roots
        if root.passes_numerical_gate and root.passes_wt_gate
    ]
    source_map = {
        "schema_version": 1,
        "status": DIAGNOSTIC_STATUS,
        "firewall": "WT_ONLY_NO_AE4_NULL_SECRETION_ACCESS",
        "source_class": N_ABS_NKCC,
        "source_scale": SOURCE_SCALE,
        "hydraulic_scale": HYDRAULIC_SCALE,
        "hydraulic_modes": list(HYDRAULIC_MODES),
        "pump_capacity_scales": list(PUMP_CONTEXTS),
        "topology_count": 5,
        "routing_count": 3,
        "expected_panel_count": 90,
        "eligibility": DIAGNOSTIC_ELIGIBILITY,
        "parent_anchor_rule": "corresponding N_ABS_NKCC S4 H1 numerical roots",
        "start_hierarchy": [3, 7, 33],
        "root_clustering": "deterministic connected components under frozen tolerance",
        "scientific_equations_gates_thresholds_and_source_values": "UNCHANGED",
        "diagnostic_only": True,
        "production_promotion_authorized": False,
        "heldout_reveal_authorized": False,
        "native_input_snapshot": native_snapshot,
    }
    summary = {
        "schema_version": 1,
        "status": DIAGNOSTIC_STATUS,
        "firewall": "WT_ONLY_NO_AE4_NULL_SECRETION_ACCESS",
        "hierarchy_complete": True,
        "panel_count": len(panel_rows),
        "screen_panel_count": len(screen),
        "screen_attempt_count": sum(len(report.attempts) for report in screen),
        "screen_wt_passing_panel_count": len(screen_ids),
        "screen_wt_passing_panel_ids": screen_ids,
        "confirmation_panel_count": len(confirmation),
        "confirmation_attempt_count": sum(
            len(report.attempts) for report in confirmation
        ),
        "confirmation_wt_passing_panel_count": len(confirmation_ids),
        "confirmation_wt_passing_panel_ids": confirmation_ids,
        "lost_at_confirmation_panel_ids": sorted(set(screen_ids) - set(confirmation_ids)),
        "geometry_panel_count": len(geometry),
        "geometry_attempt_count": sum(len(report.attempts) for report in geometry),
        "geometry_wt_passing_panel_count": len(geometry_ids),
        "geometry_wt_passing_panel_ids": geometry_ids,
        "geometry_wt_passing_root_count": len(geometry_root_ids),
        "geometry_wt_passing_root_ids": geometry_root_ids,
        "lost_at_geometry_panel_ids": sorted(set(confirmation_ids) - set(geometry_ids)),
        "attempt_history_row_count": len(attempt_rows),
        "root_history_row_count": len(root_rows),
        "parent_anchored_panel_count": sum(
            bool(parent_anchor_ids.get(member.panel_id)) for member in members
        ),
        "parent_unanchored_panel_count": sum(
            not parent_anchor_ids.get(member.panel_id) for member in members
        ),
        "seed": seed,
        "max_nfev": max_nfev,
        "diagnostic_dynamic_contract_available": bool(geometry_root_ids),
        "diagnostic_only": True,
        "production_promotion_authorized": False,
        "heldout_reveal_authorized": False,
        "native_input_snapshot": native_snapshot,
    }
    return {
        SOURCE_MAP_NAME: _json_payload(source_map),
        PANEL_NAME: _csv_payload(panel_rows, PANEL_FIELDS),
        ATTEMPT_NAME: _csv_payload(attempt_rows, ATTEMPT_FIELDS),
        ROOT_NAME: _csv_payload(root_rows, ROOT_FIELDS),
        SUMMARY_NAME: _json_payload(summary),
    }


def run_absolute_hydraulic_diagnostic(
    *,
    output: Path,
    native_output: Path,
    workers: int = 9,
    seed: int = 13601,
    max_nfev: int = 3000,
    panel_runner: PanelRunner = run_panel,
) -> Mapping[str, Any]:
    """Execute and atomically publish the full isolated WT-only hierarchy."""

    native_output = Path(native_output)
    output = Path(output)
    initial_snapshot = _native_snapshot(native_output)
    members = declared_absolute_hydraulic_diagnostic_panel()
    parent_anchors, parent_anchor_ids = _parent_anchor_maps(native_output, members)
    screen = panel_runner(
        members,
        start_count=3,
        seed=seed,
        max_nfev=max_nfev,
        workers=workers,
        anchors_by_panel=parent_anchors,
    )
    _validate_reports(screen, members, stage=SCREEN_STAGE)

    confirmation_members = _passing_members(screen)
    confirmation = (
        panel_runner(
            confirmation_members,
            start_count=7,
            seed=seed + 60000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_report_anchor_map(screen),
        )
        if confirmation_members
        else ()
    )
    _validate_reports(
        confirmation,
        confirmation_members,
        stage=CONFIRMATION_STAGE,
    )

    geometry_members = _passing_members(confirmation)
    geometry = (
        panel_runner(
            geometry_members,
            start_count=33,
            seed=seed + 70000,
            max_nfev=max_nfev,
            workers=workers,
            anchors_by_panel=_report_anchor_map(confirmation),
        )
        if geometry_members
        else ()
    )
    _validate_reports(geometry, geometry_members, stage=GEOMETRY_STAGE)

    final_snapshot = _native_snapshot(native_output)
    if final_snapshot != initial_snapshot:
        raise RuntimeError("native source inputs changed during the diagnostic")
    payloads = _artifact_payloads(
        members,
        screen,
        confirmation,
        geometry,
        parent_anchor_ids=parent_anchor_ids,
        native_snapshot=initial_snapshot,
        seed=seed,
        max_nfev=max_nfev,
    )
    _publish(output, payloads)
    return json.loads(payloads[SUMMARY_NAME].decode("utf-8"))


def _validate_publication(output: Path) -> Mapping[str, bytes]:
    commit_path = output / COMMIT_NAME
    if not commit_path.exists():
        raise RuntimeError("absolute-hydraulic diagnostic commit marker is missing")
    commit_payload = commit_path.read_bytes()
    commit = json.loads(commit_payload)
    if (
        commit.get("status") != DIAGNOSTIC_STATUS
        or commit.get("diagnostic_only") is not True
        or commit.get("production_promotion_authorized") is not False
        or commit.get("heldout_reveal_authorized") is not False
        or commit.get("hashed_artifacts") != list(HASHED_ARTIFACT_NAMES)
    ):
        raise RuntimeError("absolute-hydraulic diagnostic commit is not fail-closed")
    hash_path = output / HASH_NAME
    hash_payload = hash_path.read_bytes()
    if _sha256(hash_payload) != commit.get("hash_ledger_sha256"):
        raise RuntimeError("absolute-hydraulic diagnostic hash ledger is stale")
    rows = list(csv.DictReader(io.StringIO(hash_payload.decode("utf-8"))))
    by_name = {str(row.get("artifact", "")): row for row in rows}
    if len(by_name) != len(rows) or set(by_name) != set(HASHED_ARTIFACT_NAMES):
        raise RuntimeError("absolute-hydraulic diagnostic hash set is non-exact")
    payloads: dict[str, bytes] = {}
    for name in HASHED_ARTIFACT_NAMES:
        payload = (output / name).read_bytes()
        row = by_name[name]
        if _sha256(payload) != row.get("sha256") or len(payload) != int(
            row.get("bytes") or -1
        ):
            raise RuntimeError(f"absolute-hydraulic diagnostic artifact is stale: {name}")
        payloads[name] = payload
    if commit_path.read_bytes() != commit_payload:
        raise RuntimeError("absolute-hydraulic diagnostic commit changed during load")
    return payloads


def load_absolute_hydraulic_diagnostic_contract(
    path: str | Path,
) -> AbsoluteHydraulicDiagnosticContract:
    """Load every final WT-pass root without granting production eligibility."""

    output = Path(path)
    payloads = _validate_publication(output)
    summary = json.loads(payloads[SUMMARY_NAME])
    if (
        summary.get("status") != DIAGNOSTIC_STATUS
        or summary.get("hierarchy_complete") is not True
        or summary.get("diagnostic_only") is not True
        or summary.get("production_promotion_authorized") is not False
        or summary.get("heldout_reveal_authorized") is not False
        or int(summary.get("panel_count", -1)) != 90
        or int(summary.get("screen_panel_count", -1)) != 90
        or int(summary.get("screen_attempt_count", -1)) != 270
    ):
        raise RuntimeError("absolute-hydraulic diagnostic summary is incomplete")

    panel_rows = list(
        csv.DictReader(io.StringIO(payloads[PANEL_NAME].decode("utf-8")))
    )
    expected = {
        member.panel_id: member
        for member in declared_absolute_hydraulic_diagnostic_panel()
    }
    if len(panel_rows) != 90 or {row["panel_id"] for row in panel_rows} != set(expected):
        raise RuntimeError("absolute-hydraulic diagnostic panel set is non-exact")
    semantic_fields = (
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
    float_fields = {
        "source_scale",
        "pump_capacity_scale",
        "apical_pump_fraction",
        "apical_k_fraction",
        "ae4_cation_fraction",
        "hydraulic_scale",
    }
    for row in panel_rows:
        member_values = asdict(expected[row["panel_id"]])
        for field in semantic_fields:
            matches = (
                float(row[field]) == float(member_values[field])
                if field in float_fields
                else str(row[field]) == str(member_values[field])
            )
            if not matches:
                raise RuntimeError(f"{row['panel_id']}: diagnostic panel semantic mismatch")
        for field in (
            "diagnostic_only",
            "production_promotion_authorized",
            "heldout_reveal_authorized",
        ):
            value = _strict_bool(row[field], field=field)
            expected_value = field == "diagnostic_only"
            if value is not expected_value:
                raise RuntimeError(f"{row['panel_id']}: diagnostic authority mismatch")

    attempts = list(
        csv.DictReader(io.StringIO(payloads[ATTEMPT_NAME].decode("utf-8")))
    )
    attempt_ids = [row.get("attempt_id", "") for row in attempts]
    if len(attempt_ids) != len(set(attempt_ids)):
        raise RuntimeError("absolute-hydraulic diagnostic attempt IDs are duplicated")
    attempts_by_generation: dict[tuple[str, str], list[Mapping[str, str]]] = {}
    for row in attempts:
        panel_id = str(row.get("panel_id", ""))
        stage = str(row.get("root_search_stage", ""))
        if panel_id not in expected or stage not in STAGE_START_COUNTS:
            raise RuntimeError("absolute-hydraulic diagnostic attempt is orphaned")
        if (
            row.get("eligibility") != DIAGNOSTIC_ELIGIBILITY
            or row.get("hydraulic_mode") not in HYDRAULIC_MODES
            or _strict_bool(row.get("diagnostic_only"), field="diagnostic_only")
            is not True
            or _strict_bool(
                row.get("production_promotion_authorized"),
                field="production_promotion_authorized",
            )
            is not False
            or _strict_bool(
                row.get("heldout_reveal_authorized"),
                field="heldout_reveal_authorized",
            )
            is not False
        ):
            raise RuntimeError(f"{row.get('attempt_id')}: diagnostic attempt is unsafe")
        attempts_by_generation.setdefault((panel_id, stage), []).append(row)

    panel_by_id = {row["panel_id"]: row for row in panel_rows}
    screen_ids = list(summary.get("screen_wt_passing_panel_ids", ()))
    confirmation_ids = list(summary.get("confirmation_wt_passing_panel_ids", ()))
    geometry_ids = list(summary.get("geometry_wt_passing_panel_ids", ()))
    for name, items in (
        ("screen", screen_ids),
        ("confirmation", confirmation_ids),
        ("geometry", geometry_ids),
    ):
        if len(items) != len(set(items)) or not set(items).issubset(expected):
            raise RuntimeError(f"absolute-hydraulic {name} selection is non-exact")
    declared_screen_ids = {
        panel_id
        for panel_id, row in panel_by_id.items()
        if int(row.get("screen_wt_passing_root_count") or 0) > 0
    }
    declared_confirmation_ids = {
        panel_id
        for panel_id, row in panel_by_id.items()
        if int(row.get("confirmation_wt_passing_root_count") or 0) > 0
    }
    declared_geometry_ids = {
        panel_id
        for panel_id, row in panel_by_id.items()
        if int(row.get("geometry_wt_passing_root_count") or 0) > 0
    }
    if (
        set(screen_ids) != declared_screen_ids
        or set(confirmation_ids) != declared_confirmation_ids
        or set(geometry_ids) != declared_geometry_ids
        or not set(confirmation_ids).issubset(screen_ids)
        or not set(geometry_ids).issubset(confirmation_ids)
    ):
        raise RuntimeError("absolute-hydraulic diagnostic stage DAG is inconsistent")
    expected_stage_panels = {
        SCREEN_STAGE: set(expected),
        CONFIRMATION_STAGE: set(screen_ids),
        GEOMETRY_STAGE: set(confirmation_ids),
    }
    for stage, expected_panels in expected_stage_panels.items():
        observed_panels = {
            panel_id
            for panel_id, candidate_stage in attempts_by_generation
            if candidate_stage == stage
        }
        if observed_panels != expected_panels:
            raise RuntimeError(f"{stage}: diagnostic attempt panel selection is non-exact")
        expected_count = STAGE_START_COUNTS[stage]
        for panel_id in expected_panels:
            generation = attempts_by_generation[(panel_id, stage)]
            start_ids = [row["start_id"] for row in generation]
            if (
                len(generation) != expected_count
                or len(set(start_ids)) != expected_count
                or {"lower_quartile", "upper_quartile"} - set(start_ids)
                or any(int(row.get("root_search_start_count") or 0) != expected_count for row in generation)
            ):
                raise RuntimeError(f"{panel_id}: {stage} start contract is invalid")
    if len(attempts) != int(summary.get("attempt_history_row_count", -1)):
        raise RuntimeError("absolute-hydraulic attempt history count disagrees")

    source = output / ROOT_NAME
    all_roots = load_native_source_roots(source)
    with source.open(newline="", encoding="utf-8") as handle:
        root_rows = list(csv.DictReader(handle))
    row_by_id = {row["root_id"]: row for row in root_rows}
    if len(row_by_id) != len(root_rows) or len(all_roots) != len(root_rows):
        raise RuntimeError("absolute-hydraulic diagnostic root IDs are not unique")
    for root, row in zip(all_roots, root_rows):
        stage = str(row.get("root_search_stage", ""))
        panel_id = root.panel_id
        if panel_id not in expected or stage not in STAGE_START_COUNTS:
            raise RuntimeError(f"{root.root_id}: diagnostic root is orphaned")
        _validate_diagnostic_record(root, require_wt_pass=False)
        member = expected[panel_id]
        comparable = {
            "topology_id": root.topology_id,
            "routing_id": root.routing_id,
            "source_class": root.source_class,
            "source_scale": root.source_scale,
            "pump_capacity_scale": root.pump_capacity_scale,
            "apical_pump_fraction": root.apical_pump_fraction,
            "apical_k_fraction": root.apical_k_fraction,
            "ae4_cation_fraction": root.ae4_cation_fraction,
            "hydraulic_scale": root.hydraulic_scale,
            "hydraulic_mode": root.hydraulic_mode,
            "eligibility": root.eligibility,
        }
        if comparable != {name: getattr(member, name) for name in comparable}:
            raise RuntimeError(f"{root.root_id}: diagnostic root/member mismatch")
        start_ids = {
            attempt["start_id"]
            for attempt in attempts_by_generation.get((panel_id, stage), ())
        }
        if (
            int(row.get("root_search_start_count") or 0)
            != STAGE_START_COUNTS[stage]
            or not root.member_start_ids
            or not set(root.member_start_ids).issubset(start_ids)
            or _strict_bool(row.get("diagnostic_only"), field="diagnostic_only")
            is not True
            or _strict_bool(
                row.get("production_promotion_authorized"),
                field="production_promotion_authorized",
            )
            is not False
            or _strict_bool(
                row.get("heldout_reveal_authorized"),
                field="heldout_reveal_authorized",
            )
            is not False
        ):
            raise RuntimeError(f"{root.root_id}: diagnostic root provenance is invalid")
    if len(root_rows) != int(summary.get("root_history_row_count", -1)):
        raise RuntimeError("absolute-hydraulic root history count disagrees")

    final_ids = list(summary.get("geometry_wt_passing_root_ids", ()))
    if len(final_ids) != len(set(final_ids)):
        raise RuntimeError("absolute-hydraulic diagnostic final root IDs are duplicated")
    actual_final_ids = sorted(
        row["root_id"]
        for row in root_rows
        if row.get("root_search_stage") == GEOMETRY_STAGE
        and _strict_bool(row.get("passes_numerical_gate"), field="passes_numerical_gate")
        and _strict_bool(row.get("passes_wt_gate"), field="passes_wt_gate")
    )
    if sorted(final_ids) != actual_final_ids:
        raise RuntimeError("absolute-hydraulic final root selection is non-exact")
    final_roots: list[AbsoluteHydraulicDiagnosticRoot] = []
    roots_by_id = {root.root_id: root for root in all_roots}
    for root_id in final_ids:
        row = row_by_id.get(root_id)
        root = roots_by_id.get(root_id)
        if row is None or root is None:
            raise RuntimeError(f"missing final diagnostic root: {root_id}")
        if (
            row.get("root_search_stage") != GEOMETRY_STAGE
            or int(row.get("root_search_start_count") or 0) != 33
            or int(row.get("confirmation_start_count") or 0) != 7
            or int(row.get("geometry_start_count") or 0) != 33
            or _strict_bool(row.get("diagnostic_only"), field="diagnostic_only")
            is not True
            or _strict_bool(
                row.get("production_promotion_authorized"),
                field="production_promotion_authorized",
            )
            is not False
            or _strict_bool(
                row.get("heldout_reveal_authorized"),
                field="heldout_reveal_authorized",
            )
            is not False
        ):
            raise RuntimeError(f"{root_id}: final diagnostic root provenance is invalid")
        wrapped = AbsoluteHydraulicDiagnosticRoot(root)
        diagnostic_member_from_root(wrapped)
        final_roots.append(wrapped)
    if len(final_roots) != int(summary.get("geometry_wt_passing_root_count", -1)):
        raise RuntimeError("absolute-hydraulic diagnostic final root count disagrees")
    return AbsoluteHydraulicDiagnosticContract(tuple(final_roots), summary)


__all__ = (
    "AbsoluteHydraulicDiagnosticContract",
    "AbsoluteHydraulicDiagnosticRoot",
    "ATTEMPT_NAME",
    "COMMIT_NAME",
    "DIAGNOSTIC_ELIGIBILITY",
    "DIAGNOSTIC_STATUS",
    "HASH_NAME",
    "PANEL_NAME",
    "ROOT_NAME",
    "SOURCE_MAP_NAME",
    "SUMMARY_NAME",
    "build_absolute_hydraulic_diagnostic_model",
    "declared_absolute_hydraulic_diagnostic_panel",
    "diagnostic_member_from_root",
    "load_absolute_hydraulic_diagnostic_contract",
    "run_absolute_hydraulic_diagnostic",
)
