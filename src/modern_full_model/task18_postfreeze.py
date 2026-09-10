"""Post checkpoint WT dynamics and genotype evaluation for Task 18.

This module refuses to run unless the remotely verified WT checkpoint receipt
matches the frozen manifest.  It never changes a Task 18 capacity.  Baseline
WT, AE4 5% and AE2 results are reused only after their inherited hashes pass.
Every modified parameterisation receives new connected resting continuation
and matched production Radau trajectories.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from . import task18_fixed_wt as task
from .calibration import BoundedValue, REST_COORDINATE_NAMES, WTCalibrationSpec
from .genotype_evaluation import (
    _attempt_and_root,
    _cluster_roots,
    _core_coordinates,
    _regulatory_suffix,
    continue_genotype_expression,
    genotype_with_expression,
    simulate_genotype,
)
from .model import WT
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import build_model, routing, trajectory_path
from .task14b_five_percent import _five_percent_core, _output_path
from .validation import PRODUCTION_RADAU, physical_time_grid, sha256_file, sha256_object


EXPRESSION_AE4 = 0.05
AE4_GRID = tuple(float(value) for value in np.linspace(1.0, EXPRESSION_AE4, 20))
AE2_GRID = tuple(float(value) for value in np.linspace(1.0, 0.0, 21))
# These are numerical safety enclosures, not inherited calibration or
# physiological acceptance boxes.  Every inherited Task 13B coordinate wall
# is moved substantially outwards for the Task 18 genotype continuation.
TASK18_CONTINUATION_LIMITS = (
    (0.5, 150.0),       # intracellular Na, mM
    (5.0, 500.0),       # intracellular K, mM
    (0.05, 300.0),      # intracellular TIC, mM
    (4.5, 9.5),         # intracellular pH
    (0.05, 20.0),       # cell volume, pL
    (0.5, 700.0),       # lumen Na, mM
    (0.001, 200.0),     # lumen K, mM
    (0.01, 300.0),      # lumen TIC, mM
    (3.5, 10.5),        # lumen pH
    (0.001, 5.0),       # lumen volume, pL
)
EXPECTED_WT_CHECKPOINT = "2cdcfb3f3ab606b2c6cb84f3149f840bee1f72e0"
EXPECTED_WT_MANIFEST_BLOB = "e41e2264b2ca5bcedda7ad9e5fd26d9ac0cb9cb8"
REUSE_BLOBS = {
    "results/14_scale_free_genotype_holdout/final_artifact_hashes.csv": "876f58f5ee797a06d7b907651fe798281eaba077",
    "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv": "874cdfa2bc424f7998238eef81edc109a3202f6e",
    "results/14B_five_percent_ae4_check/artifact_hashes.csv": "b6c82d5a6e16ae6cc975edb97e24978f9d1041b1",
    "results/14B_five_percent_ae4_check/frozen_inputs.json": "56fa355d88038977152a4e6b10d297c1cd47c1a4",
}


def task18_continuation_spec() -> WTCalibrationSpec:
    """Return tolerances with noninherited numerical coordinate enclosures."""

    inherited = WTCalibrationSpec()
    bounds = tuple(
        BoundedValue(
            name=item.name,
            lower=lower,
            preferred=item.preferred,
            upper=upper,
            units=item.units,
            tier=item.tier,
            provenance="TASK18_NUMERICAL_ENCLOSURE",
            note="broad positive solver enclosure; not an acceptance interval",
        )
        for item, (lower, upper) in zip(
            inherited.coordinate_bounds, TASK18_CONTINUATION_LIMITS, strict=True
        )
    )
    return WTCalibrationSpec(
        coordinate_bounds=bounds,
        parameter_bounds=inherited.parameter_bounds,
        independent_rhs_scales=inherited.independent_rhs_scales,
        regulatory_rhs_scales_s_inv=inherited.regulatory_rhs_scales_s_inv,
        wt_cl_target_mM=inherited.wt_cl_target_mM,
        wt_cl_sem_mM=inherited.wt_cl_sem_mM,
        wt_ph_target=inherited.wt_ph_target,
        wt_ph_sem=inherited.wt_ph_sem,
        volume_lineage_target_pL=inherited.volume_lineage_target_pL,
        mean_voltage_lineage_target_mV=inherited.mean_voltage_lineage_target_mV,
        resting_calcium_uM=inherited.resting_calcium_uM,
        root_scaled_tolerance=inherited.root_scaled_tolerance,
        omitted_row_raw_tolerance=inherited.omitted_row_raw_tolerance,
        charge_tolerance_fmol=inherited.charge_tolerance_fmol,
        current_tolerance_A=inherited.current_tolerance_A,
        boundary_relative_tolerance=inherited.boundary_relative_tolerance,
        cluster_relative_tolerance=inherited.cluster_relative_tolerance,
        jacobian_rank_relative_tolerance=inherited.jacobian_rank_relative_tolerance,
    )


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def checked_checkpoint() -> tuple[dict[str, Any], dict[str, Any]]:
    receipt = task.read(task.OUT / "pushed_wt_checkpoint.json")
    frozen_path = task.OUT / "frozen_wt_manifest.json"
    if receipt["branch"] != task.BRANCH:
        raise AssertionError("Wrong WT checkpoint branch")
    if receipt["checkpoint_sha"] != EXPECTED_WT_CHECKPOINT:
        raise AssertionError("Unexpected WT checkpoint commit")
    if receipt["remote_manifest_git_blob_sha"] != EXPECTED_WT_MANIFEST_BLOB:
        raise AssertionError("Unexpected remote WT manifest blob")
    if git_blob_sha1(frozen_path) != EXPECTED_WT_MANIFEST_BLOB:
        raise AssertionError("Local WT manifest differs from the pushed checkpoint")
    if not receipt["remote_ref_verified_utc"]:
        raise AssertionError("WT checkpoint was not remotely verified")
    frozen = task.read(frozen_path)
    if frozen["genotype_outcomes_used"] or frozen["genotype_result_files_read"]:
        raise AssertionError("WT freeze phenotype firewall failed")
    if len(frozen["feasible_parameterisations"]) != 30:
        raise AssertionError("Incomplete frozen WT ensemble")
    for name, digest in frozen["hashes"].items():
        if sha256_file(task.REPO / name) != digest:
            raise AssertionError(f"Frozen WT artifact changed: {name}")
    return frozen, receipt


def _artifact_hashes(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["path"]: row["sha256"] for row in csv.DictReader(handle)}


def verify_reuse_inputs() -> dict[str, str]:
    for name, expected in REUSE_BLOBS.items():
        if git_blob_sha1(task.REPO / name) != expected:
            raise AssertionError(f"Inherited reuse anchor changed: {name}")
    task14 = _artifact_hashes(
        task.REPO / "results/14_scale_free_genotype_holdout/final_artifact_hashes.csv"
    )
    task14b = _artifact_hashes(
        task.REPO / "results/14B_five_percent_ae4_check/artifact_hashes.csv"
    )
    required: list[str] = [
        "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv"
    ]
    manifest, _ = load_freeze()
    for root_id in sorted(manifest["roots"]):
        continuation = (
            Path("results/14_scale_free_genotype_holdout/continuation")
            / f"{task.slug(root_id)}_AE4.json"
        )
        required.append(str(continuation))
        for calcium in task.CALCIUM:
            required.append(relative(trajectory_path(root_id, calcium, "WT", "production_radau")))
            required.append(relative(_output_path(root_id, calcium)))
    verified = {}
    for name in required:
        path = task.REPO / name
        expected = task14.get(name, task14b.get(name))
        if expected is None or sha256_file(path) != expected:
            raise AssertionError(f"Hash invalid baseline reuse input: {name}")
        verified[name] = expected
    return verified


def relative(path: Path) -> str:
    return str(path.relative_to(task.REPO))


def _near_ae4_continuation(model, wt_state: Sequence[float]) -> dict[str, Any]:
    """Continue only to 0.05; exact AE4 zero is neither requested nor attempted."""

    spec = task18_continuation_spec()
    lower, upper = spec.coordinate_lower.copy(), spec.coordinate_upper.copy()
    span = upper - lower
    previous = _core_coordinates(model, wt_state)
    if np.any(previous <= lower) or np.any(previous >= upper):
        raise ValueError("Frozen WT state lies outside continuation coordinates")
    suffix = _regulatory_suffix(model)
    steps = []
    selected_root = None
    for expression in AE4_GRID:
        genotype = genotype_with_expression("AE4", expression)
        interior = 1.0e-10 * span
        starts = (
            ("previous_connected_root", previous.copy()),
            ("local_minus", np.clip(previous - 1.0e-3 * span, lower + interior, upper - interior)),
            ("local_plus", np.clip(previous + 1.0e-3 * span, lower + interior, upper - interior)),
        )
        attempts = []
        roots = []
        for start_id, start in starts:
            attempt, root = _attempt_and_root(
                model=model,
                transporter="AE4",
                expression=expression,
                genotype=genotype,
                regulatory_suffix=suffix,
                spec=spec,
                start_id=start_id,
                start=start,
                lower=lower,
                upper=upper,
                max_nfev=3000,
            )
            attempts.append(attempt)
            if root is not None:
                roots.append(root)
        retained = _cluster_roots(
            roots,
            expression=expression,
            transporter="AE4",
            span=span,
            relative_tolerance=spec.cluster_relative_tolerance,
        )
        if not retained:
            steps.append(
                {
                    "expression": expression,
                    "status": "NO_ADMISSIBLE_ROOT",
                    "attempts": [asdict(value) for value in attempts],
                    "root_count": 0,
                }
            )
            return {
                "transporter": "AE4",
                "terminal_expression": EXPRESSION_AE4,
                "exact_zero_attempted": False,
                "complete": False,
                "failure": f"No admissible root at expression {expression}",
                "steps": steps,
                "selected_core_state": None,
            }
        distances = {
            root.root_id: float(
                np.max(np.abs(np.asarray(root.coordinates) - previous) / span)
            )
            for root in retained
        }
        selected_root = min(
            retained,
            key=lambda root: (distances[root.root_id], root.max_abs_scaled_independent_rhs),
        )
        distance = distances[selected_root.root_id]
        connected = distance <= 0.25
        steps.append(
            {
                "expression": expression,
                "status": "CONNECTED" if connected else "ROOT_JUMP",
                "root_count": len(retained),
                "alternate_root_count": len(retained) - 1,
                "selected_root_id": selected_root.root_id,
                "selected_distance_normalized": distance,
                "selected_max_abs_scaled_independent_rhs": selected_root.max_abs_scaled_independent_rhs,
                "selected_rank": selected_root.normalized_jacobian_rank,
                "attempts": [asdict(value) for value in attempts],
            }
        )
        if not connected:
            return {
                "transporter": "AE4",
                "terminal_expression": EXPRESSION_AE4,
                "exact_zero_attempted": False,
                "complete": False,
                "failure": f"Connected branch jump at expression {expression}",
                "steps": steps,
                "selected_core_state": None,
            }
        previous = np.asarray(selected_root.coordinates)
    assert selected_root is not None and selected_root.expression == EXPRESSION_AE4
    return {
        "transporter": "AE4",
        "terminal_expression": EXPRESSION_AE4,
        "exact_zero_attempted": False,
        "complete": True,
        "failure": None,
        "steps": steps,
        "selected_core_state": list(selected_root.core_state),
        "selected_core_state_sha256": sha256_object(selected_root.core_state),
        "selected_root": asdict(selected_root),
    }


def _continuation_summary(value) -> dict[str, Any]:
    selected = value.steps[-1].selected_root if value.steps else None
    return {
        "transporter": value.transporter,
        "expression_grid": list(value.expression_grid),
        "terminal_expression": value.expression_grid[-1],
        "exact_zero_attempted": value.expression_grid[-1] == 0.0,
        "complete": value.completed_to_exact_zero,
        "failure": value.branch_failure,
        "fixed_other_impermeant_osmoles_fmol": value.fixed_other_impermeant_osmoles_fmol,
        "fitted_parameter_count": value.fitted_parameter_count,
        "steps": [
            {
                "expression": step.expression,
                "status": step.status,
                "root_count": len(step.roots),
                "alternate_root_count": step.alternate_root_count,
                "selected_root_id": step.selected_root_id,
                "selected_distance_normalized": step.selected_distance_from_previous_normalized,
                "selected_max_abs_scaled_independent_rhs": (
                    step.selected_root.max_abs_scaled_independent_rhs
                    if step.selected_root is not None
                    else None
                ),
                "attempts": [asdict(attempt) for attempt in step.attempts],
            }
            for step in value.steps
        ],
        "selected_core_state": list(selected.core_state) if selected is not None else None,
        "selected_core_state_sha256": (
            sha256_object(selected.core_state) if selected is not None else None
        ),
        "selected_root": asdict(selected) if selected is not None else None,
    }


def _trajectory_summary(trajectory) -> dict[str, Any]:
    complete = bool(
        trajectory.success
        and trajectory.time_s.size == physical_time_grid(duration_s=600.0, sample_step_s=1.0).size
        and trajectory.time_s[-1] == 600.0
    )
    valid = bool(complete and trajectory.numerical_gate_pass)
    total = float(trajectory.cumulative_flow_pL[-1]) if complete else None
    landmarks = {}
    if complete:
        for minute in range(60, 601, 60):
            index = int(np.argmin(np.abs(trajectory.time_s - minute)))
            landmarks[str(minute)] = {
                "time_s": float(trajectory.time_s[index]),
                "flow_pL_s": float(trajectory.flow_pL_s[index]),
                "cumulative_flow_pL": float(trajectory.cumulative_flow_pL[index]),
            }
    return {
        "status": "COMPLETE" if valid else "NUMERICAL_FAILURE",
        "solver_success": bool(trajectory.success),
        "numerical_gate_pass": valid,
        "message": trajectory.message,
        "time_sample_count": int(trajectory.time_s.size),
        "total_0_600_pL": total,
        "mean_0_600_pL_s": total / 600.0 if total is not None else None,
        "positive_core": bool(trajectory.positive_core),
        "all_flow_nonnegative": bool(trajectory.all_flow_nonnegative),
        "minimum_flow_pL_s": float(np.nanmin(trajectory.flow_pL_s)),
        "maximum_flow_pL_s": float(np.nanmax(trajectory.flow_pL_s)),
        "max_dimensionless_conservation_ratio": float(
            trajectory.max_dimensionless_conservation_ratio
        ),
        "max_abs_conservation_residuals": dict(trajectory.max_abs_conservation_residuals),
        "max_abs_deleted_transport_flux_fmol_s": float(
            trajectory.max_abs_deleted_transport_flux_fmol_s
        ),
        "initial": {
            "na_mM": float(trajectory.cell_na_mM[0]),
            "k_mM": float(trajectory.cell_k_mM[0]),
            "cl_mM": float(trajectory.cell_cl_mM[0]),
            "ph": float(trajectory.cell_ph[0]),
            "cell_volume_pL": float(trajectory.cell_volume_pL[0]),
            "lumen_volume_pL": float(trajectory.lumen_volume_pL[0]),
        },
        "endpoint": {
            "na_mM": float(trajectory.cell_na_mM[-1]),
            "k_mM": float(trajectory.cell_k_mM[-1]),
            "cl_mM": float(trajectory.cell_cl_mM[-1]),
            "ph": float(trajectory.cell_ph[-1]),
            "cell_volume_pL": float(trajectory.cell_volume_pL[-1]),
            "lumen_volume_pL": float(trajectory.lumen_volume_pL[-1]),
        },
        "landmarks": landmarks,
        "state_array_sha256": sha256_object(trajectory.states.tolist()),
        "flow_array_sha256": sha256_object(trajectory.flow_pL_s.tolist()),
    }


def _failed_trajectory(error: BaseException) -> dict[str, Any]:
    return {
        "status": "EXCEPTION",
        "solver_success": False,
        "numerical_gate_pass": False,
        "message": f"{type(error).__name__}: {error}",
        "total_0_600_pL": None,
    }


def _simulate(model, core, transporter, genotype) -> dict[str, Any]:
    try:
        trajectory = simulate_genotype(
            model,
            core,
            transporter=transporter,
            genotype=genotype,
            solver=PRODUCTION_RADAU,
            time_s=physical_time_grid(duration_s=600.0, sample_step_s=1.0),
        )
        return _trajectory_summary(trajectory)
    except BaseException as error:
        return _failed_trajectory(error)


def _run_modified(candidate_id: str) -> dict[str, Any]:
    manifest, _ = load_freeze()
    payload_path = task.OUT / "wt_parameter_payloads" / f"{candidate_id}.json"
    payload = task.read(payload_path)
    row = payload["row"]
    if row["condition"] == "inherited_baseline":
        raise ValueError("Baseline must use the hash valid reuse path")
    root_id = row["root_id"]
    template, state, baseline_evaluation = task.baseline(manifest, root_id)
    model = task.model_with_capacities(template, payload["accepted_capacities"])
    verification = task.residual_summary(model.evaluate(0.0, state))
    if verification["max_abs_scaled_independent_rhs"] > task.SPEC.root_scaled_tolerance:
        raise AssertionError("Frozen WT state failed before genotype continuation")

    try:
        near = _near_ae4_continuation(model, state)
    except BaseException as error:
        near = {
            "transporter": "AE4",
            "terminal_expression": EXPRESSION_AE4,
            "exact_zero_attempted": False,
            "complete": False,
            "failure": f"{type(error).__name__}: {error}",
            "steps": [],
            "selected_core_state": None,
        }
    try:
        ae2_value = continue_genotype_expression(
            model,
            state,
            transporter="AE2",
            expression_grid=AE2_GRID,
            spec=task18_continuation_spec(),
        )
        ae2 = _continuation_summary(ae2_value)
    except BaseException as error:
        ae2 = {
            "transporter": "AE2",
            "terminal_expression": 0.0,
            "exact_zero_attempted": True,
            "complete": False,
            "failure": f"{type(error).__name__}: {error}",
            "steps": [],
            "selected_core_state": None,
        }

    rows = []
    for calcium in task.CALCIUM:
        dynamic_template = build_model(manifest, root_id, calcium)
        dynamic_model = task.model_with_capacities(
            dynamic_template, payload["accepted_capacities"]
        )
        wt = _simulate(dynamic_model, state[:12], "AE4", WT)
        ae4_result = (
            _simulate(
                dynamic_model,
                near["selected_core_state"],
                "AE4",
                genotype_with_expression("AE4", EXPRESSION_AE4),
            )
            if near["complete"]
            else {
                "status": "CONTINUATION_FAILED",
                "numerical_gate_pass": False,
                "total_0_600_pL": None,
                "message": near["failure"],
            }
        )
        ae2_result = (
            _simulate(
                dynamic_model,
                ae2["selected_core_state"],
                "AE2",
                genotype_with_expression("AE2", 0.0),
            )
            if ae2["complete"]
            else {
                "status": "CONTINUATION_FAILED",
                "numerical_gate_pass": False,
                "total_0_600_pL": None,
                "message": ae2["failure"],
            }
        )
        wt_total = wt.get("total_0_600_pL")
        ae4_total = ae4_result.get("total_0_600_pL")
        ae2_total = ae2_result.get("total_0_600_pL")
        paired_ae4 = bool(
            wt.get("numerical_gate_pass")
            and ae4_result.get("numerical_gate_pass")
            and wt_total is not None
            and wt_total > 0.0
        )
        paired_ae2 = bool(
            wt.get("numerical_gate_pass")
            and ae2_result.get("numerical_gate_pass")
            and wt_total is not None
            and wt_total > 0.0
        )
        ratio4 = ae4_total / wt_total if paired_ae4 else None
        ratio2 = ae2_total / wt_total if paired_ae2 else None
        rows.append(
            {
                "candidate_id": candidate_id,
                "root_id": root_id,
                "routing_family": row["routing_family"],
                "condition": row["condition"],
                "target_share": row["target_share"],
                "calcium_uM": calcium,
                "source": "NEW_POST_CHECKPOINT_SIMULATION",
                "capacities_refitted_after_freeze": False,
                "wt": wt,
                "ae4_5pct": ae4_result,
                "ae2_loss": ae2_result,
                "R_AE4_5pct": ratio4,
                "D_AE4_5pct": 1.0 - ratio4 if ratio4 is not None else None,
                "R_AE2": ratio2,
                "D_AE2": 1.0 - ratio2 if ratio2 is not None else None,
                "paired_ae4_numerical_gate_pass": paired_ae4,
                "paired_ae2_numerical_gate_pass": paired_ae2,
            }
        )
    result = {
        "candidate_id": candidate_id,
        "root_id": root_id,
        "condition": row["condition"],
        "target_share": row["target_share"],
        "wt_payload_path": relative(payload_path),
        "wt_payload_sha256": sha256_file(payload_path),
        "wt_checkpoint_sha": EXPECTED_WT_CHECKPOINT,
        "capacity_multipliers": payload["capacity_multipliers"],
        "wt_fixed_state_verification": verification,
        "ae4_5pct_continuation": near,
        "ae2_loss_continuation": ae2,
        "rows": rows,
        "fitted_parameter_count_after_freeze": 0,
        "exact_ae4_zero_attempts": 0,
    }
    write_json(task.OUT / "postfreeze_payloads" / f"{candidate_id}.json", result)
    return result


def _saved_summary(trace: Mapping[str, Any]) -> dict[str, Any]:
    if not trace["success"] or not trace["numerical_gate_pass"]:
        raise AssertionError("Reused trajectory does not pass its frozen numerical gate")
    total = float(trace["cumulative_flow_pL"][-1])
    if not math.isclose(
        total,
        float(np.trapezoid(trace["flow_pL_s"], trace["time_s"])),
        rel_tol=1.0e-12,
        abs_tol=1.0e-14,
    ):
        raise AssertionError("Reused trajectory integral changed")
    return {
        "status": "COMPLETE",
        "solver_success": True,
        "numerical_gate_pass": True,
        "message": trace["message"],
        "time_sample_count": len(trace["time_s"]),
        "total_0_600_pL": total,
        "mean_0_600_pL_s": total / 600.0,
        "positive_core": trace["positive_core"],
        "all_flow_nonnegative": trace["all_flow_nonnegative"],
        "minimum_flow_pL_s": min(trace["flow_pL_s"]),
        "maximum_flow_pL_s": max(trace["flow_pL_s"]),
        "max_dimensionless_conservation_ratio": trace[
            "max_dimensionless_conservation_ratio"
        ],
        "max_abs_conservation_residuals": trace["max_abs_conservation_residuals"],
        "max_abs_deleted_transport_flux_fmol_s": trace[
            "max_abs_deleted_transport_flux_fmol_s"
        ],
    }


def _baseline_rows(verified: Mapping[str, str]) -> list[dict[str, Any]]:
    manifest, _ = load_freeze()
    with (
        task.REPO / "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv"
    ).open(newline="", encoding="utf-8") as handle:
        source = list(csv.DictReader(handle))
    ae2 = {(row["root_id"], float(row["calcium_uM"])): row for row in source}
    if len(ae2) != 30:
        raise AssertionError("Incomplete inherited AE2 panel")
    rows = []
    for root_id in sorted(manifest["roots"]):
        payload_path = (
            task.OUT
            / "wt_parameter_payloads"
            / f"{task.slug(root_id)}__inherited_baseline.json"
        )
        payload = task.read(payload_path)
        _five_percent_core(root_id)
        for calcium in task.CALCIUM:
            wt_path = trajectory_path(root_id, calcium, "WT", "production_radau")
            near_path = _output_path(root_id, calcium)
            if relative(wt_path) not in verified or relative(near_path) not in verified:
                raise AssertionError("Unverified baseline reuse path")
            wt_data = task.read(wt_path)
            near_data = task.read(near_path)
            wt = _saved_summary(wt_data["trajectory"])
            near = _saved_summary(near_data["trajectory"])
            source2 = ae2[root_id, calcium]
            if source2["status"] != "COMPLETE" or source2["paired_numerical_gate_pass"] != "True":
                raise AssertionError("Inherited AE2 comparison is not valid")
            ratio4 = near["total_0_600_pL"] / wt["total_0_600_pL"]
            ratio2 = float(source2["R_total"])
            rows.append(
                {
                    "candidate_id": payload["row"]["candidate_id"],
                    "root_id": root_id,
                    "routing_family": routing(manifest, root_id),
                    "condition": "inherited_baseline",
                    "target_share": None,
                    "calcium_uM": calcium,
                    "source": "HASH_VALID_INHERITED_REUSE",
                    "capacities_refitted_after_freeze": False,
                    "wt": {**wt, "source_path": relative(wt_path), "source_sha256": verified[relative(wt_path)]},
                    "ae4_5pct": {**near, "source_path": relative(near_path), "source_sha256": verified[relative(near_path)]},
                    "ae2_loss": {
                        "status": "COMPLETE",
                        "numerical_gate_pass": True,
                        "total_0_600_pL": ratio2 * wt["total_0_600_pL"],
                        "source_path": "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv",
                        "source_sha256": verified[
                            "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv"
                        ],
                    },
                    "R_AE4_5pct": ratio4,
                    "D_AE4_5pct": 1.0 - ratio4,
                    "R_AE2": ratio2,
                    "D_AE2": 1.0 - ratio2,
                    "paired_ae4_numerical_gate_pass": True,
                    "paired_ae2_numerical_gate_pass": True,
                }
            )
    return rows


def _flat_rows(rows: Sequence[Mapping[str, Any]], kind: str) -> list[dict[str, Any]]:
    result = []
    key = {"wt": "wt", "ae4": "ae4_5pct", "ae2": "ae2_loss"}[kind]
    for row in rows:
        trajectory = row[key]
        common = {
            field: row.get(field)
            for field in (
                "candidate_id",
                "root_id",
                "routing_family",
                "condition",
                "target_share",
                "calcium_uM",
                "source",
                "capacities_refitted_after_freeze",
            )
        }
        common.update(
            status=trajectory.get("status"),
            numerical_gate_pass=trajectory.get("numerical_gate_pass"),
            total_0_600_pL=trajectory.get("total_0_600_pL"),
            mean_0_600_pL_s=trajectory.get("mean_0_600_pL_s"),
            solver_success=trajectory.get("solver_success"),
            message=trajectory.get("message"),
            minimum_flow_pL_s=trajectory.get("minimum_flow_pL_s"),
            maximum_flow_pL_s=trajectory.get("maximum_flow_pL_s"),
            max_dimensionless_conservation_ratio=trajectory.get(
                "max_dimensionless_conservation_ratio"
            ),
            source_path=trajectory.get("source_path"),
            source_sha256=trajectory.get("source_sha256"),
        )
        if kind == "ae4":
            common.update(
                ae4_expression=EXPRESSION_AE4,
                paired_numerical_gate_pass=row["paired_ae4_numerical_gate_pass"],
                secretion_ratio=row["R_AE4_5pct"],
                secretion_reduction=row["D_AE4_5pct"],
            )
        elif kind == "ae2":
            common.update(
                ae2_expression=0.0,
                paired_numerical_gate_pass=row["paired_ae2_numerical_gate_pass"],
                secretion_ratio=row["R_AE2"],
                secretion_reduction=row["D_AE2"],
            )
        result.append(common)
    return result


def evaluate(*, workers: int) -> None:
    if (task.OUT / "final_comparison.csv").exists():
        raise FileExistsError("Task 18 post checkpoint results are immutable")
    frozen, receipt = checked_checkpoint()
    verified = verify_reuse_inputs()
    rows = _baseline_rows(verified)
    modified_ids = sorted(
        item["candidate_id"]
        for item in frozen["all_decisions"]
        if item["condition"] != "inherited_baseline"
    )
    results = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_run_modified, candidate): candidate for candidate in modified_ids}
        for future in as_completed(futures):
            candidate = futures[future]
            result = future.result()
            results.append(result)
            print(
                json.dumps(
                    {
                        "candidate": candidate,
                        "ae4_continuation": result["ae4_5pct_continuation"]["complete"],
                        "ae2_continuation": result["ae2_loss_continuation"]["complete"],
                        "valid_wt": sum(row["wt"]["numerical_gate_pass"] for row in result["rows"]),
                        "valid_ae4": sum(row["paired_ae4_numerical_gate_pass"] for row in result["rows"]),
                        "valid_ae2": sum(row["paired_ae2_numerical_gate_pass"] for row in result["rows"]),
                    }
                ),
                flush=True,
            )
    for result in sorted(results, key=lambda item: item["candidate_id"]):
        rows.extend(result["rows"])
    rows.sort(key=lambda row: (row["root_id"], task.CONDITIONS.index(row["condition"]), row["calcium_uM"]))
    if len(rows) != 90:
        raise AssertionError("Incomplete Task 18 dynamic ensemble")

    write_rows(task.OUT / "wt_dynamic_results.csv", _flat_rows(rows, "wt"))
    write_rows(task.OUT / "ae4_5pct_results.csv", _flat_rows(rows, "ae4"))
    write_rows(task.OUT / "ae2_results.csv", _flat_rows(rows, "ae2"))
    decision = {item["candidate_id"]: item for item in frozen["all_decisions"]}
    comparison = []
    for row in rows:
        parent = decision[row["candidate_id"]]
        comparison.append(
            {
                "candidate_id": row["candidate_id"],
                "root_id": row["root_id"],
                "routing_family": row["routing_family"],
                "condition": row["condition"],
                "target_share": row["target_share"],
                "realized_ae4_share": parent["realized_ae4_share"],
                "calcium_uM": row["calcium_uM"],
                "largest_abs_log_fold": parent["largest_abs_log_fold"],
                "largest_fold_multiplier": parent["largest_fold_multiplier"],
                "g_k_total_nS": parent["g_k_total_nS"],
                "g_cl_apical_nS": parent["g_cl_apical_nS"],
                "v_apical_mV": parent["v_apical_mV"],
                "v_basolateral_mV": parent["v_basolateral_mV"],
                "v_transepithelial_mV": parent["v_transepithelial_mV"],
                "ae4_productive_fraction": parent["ae4_productive_fraction"],
                "ae4_cancellation_fraction": parent["ae4_cancellation_fraction"],
                "wt_total_0_600_pL": row["wt"].get("total_0_600_pL"),
                "wt_numerical_gate_pass": row["wt"].get("numerical_gate_pass"),
                "ae4_5pct_total_0_600_pL": row["ae4_5pct"].get("total_0_600_pL"),
                "R_AE4_5pct": row["R_AE4_5pct"],
                "D_AE4_5pct": row["D_AE4_5pct"],
                "ae4_paired_numerical_gate_pass": row["paired_ae4_numerical_gate_pass"],
                "ae2_total_0_600_pL": row["ae2_loss"].get("total_0_600_pL"),
                "R_AE2": row["R_AE2"],
                "D_AE2": row["D_AE2"],
                "ae2_paired_numerical_gate_pass": row["paired_ae2_numerical_gate_pass"],
                "source": row["source"],
                "capacities_refitted_after_freeze": False,
            }
        )
    write_rows(task.OUT / "final_comparison.csv", comparison)

    modified = [row for row in rows if row["condition"] != "inherited_baseline"]
    valid4 = [row for row in modified if row["paired_ae4_numerical_gate_pass"]]
    valid2 = [row for row in modified if row["paired_ae2_numerical_gate_pass"]]
    summary = {
        "created_utc": task.now(),
        "wt_checkpoint_sha": receipt["checkpoint_sha"],
        "expected_rows": 90,
        "baseline_reused_rows": 30,
        "new_modified_rows": 60,
        "new_wt_trajectories": 60,
        "new_ae4_5pct_trajectories_attempted": sum(
            3 for result in results if result["ae4_5pct_continuation"]["complete"]
        ),
        "new_ae2_loss_trajectories_attempted": sum(
            3 for result in results if result["ae2_loss_continuation"]["complete"]
        ),
        "ae4_exact_zero_attempts": 0,
        "postfreeze_parameter_refits": 0,
        "continuation_coordinate_enclosures": {
            name: [lower, upper]
            for name, (lower, upper) in zip(
                REST_COORDINATE_NAMES, TASK18_CONTINUATION_LIMITS, strict=True
            )
        },
        "modified_wt_numerically_valid": sum(row["wt"]["numerical_gate_pass"] for row in modified),
        "modified_ae4_pairs_numerically_valid": len(valid4),
        "modified_ae2_pairs_numerically_valid": len(valid2),
        "ae4_ratio_range": (
            [min(row["R_AE4_5pct"] for row in valid4), max(row["R_AE4_5pct"] for row in valid4)]
            if valid4
            else None
        ),
        "ae2_ratio_range": (
            [min(row["R_AE2"] for row in valid2), max(row["R_AE2"] for row in valid2)]
            if valid2
            else None
        ),
        "ae4_reduces_count": sum(row["R_AE4_5pct"] < 1.0 for row in valid4),
        "ae4_increases_count": sum(row["R_AE4_5pct"] > 1.0 for row in valid4),
        "ae2_reduces_count": sum(row["R_AE2"] < 1.0 for row in valid2),
        "ae2_increases_count": sum(row["R_AE2"] > 1.0 for row in valid2),
        "all_frozen_candidates_accounted": True,
        "hash_valid_reuse_files": len(verified),
        "result_hashes": {
            relative(task.OUT / name): sha256_file(task.OUT / name)
            for name in (
                "wt_dynamic_results.csv",
                "ae4_5pct_results.csv",
                "ae2_results.csv",
                "final_comparison.csv",
            )
        },
    }
    write_json(task.OUT / "postfreeze_summary.json", summary)
    print(json.dumps(summary), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    if not 1 <= args.workers <= 9:
        raise ValueError("workers must lie between one and nine")
    evaluate(workers=args.workers)


if __name__ == "__main__":
    main()
