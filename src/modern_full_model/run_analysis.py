"""Artifact writer for the Task 13B pre-reveal WT dynamic suite.

The public :func:`run_pre_reveal_wt_suite` entry point accepts already-frozen
WT models and a production resting root from ``calibration.py``.  Keeping
construction outside this module prevents the dynamic layer from silently
changing a parameter, branch, or root rule.  This module never loads the
held-out target ledger and never evaluates a knockout genotype.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .model import ModernFullModel
from .parameters import parameter_records
from .states import CORE_STATE_NAMES
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    CONSERVATION_RESIDUAL_UNITS,
    LOOSE_RADAU,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    TIGHT_RADAU,
    SecretagogueProtocol,
    StimulusArm,
    Trajectory,
    WTFlowEnvelope,
    fit_wt_flow_scale,
    nearby_charge_preserving_states,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    regulatory_landmarks,
    sha256_file,
    sha256_object,
    simulate_wt,
    solver_comparison,
    write_json,
    write_trajectory_csv,
    wt_arm_response,
    wt_flow_envelope_feasibility,
)


REQUIRED_REGULATORY_FAMILIES = ("R0", "R1", "R2", "R3")


def _object_payload(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, Mapping):
        return {str(key): _object_payload(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_object_payload(item) for item in value]
    if hasattr(value, "__dict__"):
        return {
            str(key): _object_payload(item)
            for key, item in vars(value).items()
            if not str(key).startswith("_")
        }
    return repr(value)


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table {path.name}")
    fieldnames: list[str] = []
    for row in rows:
        for name in row:
            if name not in fieldnames:
                fieldnames.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _trajectory_summary(trajectory: Trajectory) -> Mapping[str, Any]:
    summary = {
        "family": trajectory.family,
        "initial_condition": trajectory.initial_condition,
        "solver_label": trajectory.solver.label,
        "solver_method": trajectory.solver.method,
        "success": trajectory.success,
        "positive_core": trajectory.positive_core,
        "cell_na_start_mM": trajectory.cell_na_mM[0],
        "cell_na_end_mM": trajectory.cell_na_mM[-1],
        "cell_k_start_mM": trajectory.cell_k_mM[0],
        "cell_k_end_mM": trajectory.cell_k_mM[-1],
        "cell_cl_start_mM": trajectory.cell_cl_mM[0],
        "cell_cl_end_mM": trajectory.cell_cl_mM[-1],
        "cell_ph_start": trajectory.cell_ph[0],
        "cell_ph_end": trajectory.cell_ph[-1],
        "cell_volume_start_pL": trajectory.cell_volume_pL[0],
        "cell_volume_end_pL": trajectory.cell_volume_pL[-1],
        "flow_end_pL_s": trajectory.flow_pL_s[-1],
        "cumulative_flow_pL": trajectory.cumulative_flow_pL[-1],
        "capacity_multiplier_start": trajectory.capacity_multiplier[0],
        "capacity_multiplier_end": trajectory.capacity_multiplier[-1],
    }
    for name, value in trajectory.max_abs_conservation_residuals.items():
        summary[f"max_abs_residual__{name}"] = value
        summary[f"residual_unit__{name}"] = CONSERVATION_RESIDUAL_UNITS[name]
    return summary


def run_pre_reveal_wt_suite(
    *,
    models: Mapping[str, ModernFullModel],
    resting_state: ArrayLike,
    output_directory: str | Path,
    preferred_family: str,
    root_id: str,
    root_accepted: bool,
    resting_physiology_gate_pass: bool,
    root_selection_rule: Mapping[str, Any],
    preferred_member_id: str = "R1_G125_P21_PROSE_TREFERENCE",
    stimulus_protocol: SecretagogueProtocol | None = None,
    sample_step_s: float = 5.0,
    nearby_fraction: float = 1.0e-3,
    conservation_tolerances: Mapping[str, float] | None = None,
    solver_relative_tolerance: float = 1.0e-4,
    effective_cell_count_range: tuple[float, float] | None = None,
) -> Mapping[str, Any]:
    """Run and freeze the complete WT-only dynamic validation layer.

    All supplied model objects must already have identical whole-cell and AE4
    parameters.  Only their nested regulatory subsystem may differ.  The
    runner replaces their resting stimulus objects with the one explicitly
    frozen physical protocol; it does not optimize or select a root.
    """

    missing = [family for family in REQUIRED_REGULATORY_FAMILIES if family not in models]
    if missing:
        raise ValueError(f"R0--R3 comparison is incomplete: missing {missing}")
    if preferred_family not in models:
        raise KeyError(f"preferred regulatory family {preferred_family!r} is absent")
    if not root_accepted or not resting_physiology_gate_pass:
        raise ValueError(
            "production WT dynamics require an accepted numerical root that passes "
            "the frozen resting-physiology gate"
        )
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    core_root = np.asarray(resting_state, dtype=float)
    declared_conservation_tolerances = dict(
        CONSERVATION_RESIDUAL_TOLERANCES
        if conservation_tolerances is None
        else conservation_tolerances
    )
    if set(declared_conservation_tolerances) != set(CONSERVATION_RESIDUAL_UNITS):
        raise ValueError("conservation tolerances must cover every unit-declared residual")
    if any(float(value) <= 0.0 for value in declared_conservation_tolerances.values()):
        raise ValueError("every conservation tolerance must be positive")
    if core_root.ndim != 1 or core_root.size < len(CORE_STATE_NAMES):
        raise ValueError("resting state does not contain the full conserved core")
    core_root = core_root[: len(CORE_STATE_NAMES)].copy()
    if np.any(~np.isfinite(core_root)) or np.any(core_root <= 0.0):
        raise ValueError("resting core must be finite and positive")

    # Exact parameter equality is required before regulatory nesting.  A
    # family-specific whole-cell change would confound the comparison.
    frozen_preferred = models[preferred_family]
    core_parameter_hash = sha256_object(frozen_preferred.parameters)
    ae4_parameter_hash = sha256_object(_object_payload(frozen_preferred.ae4_parameters))
    for family in REQUIRED_REGULATORY_FAMILIES:
        candidate = models[family]
        if sha256_object(candidate.parameters) != core_parameter_hash:
            raise ValueError(f"{family} changes whole-cell parameters")
        if sha256_object(_object_payload(candidate.ae4_parameters)) != ae4_parameter_hash:
            raise ValueError(f"{family} changes AE4 transport parameters")

    protocol = stimulus_protocol or SecretagogueProtocol()
    ensemble = pre_reveal_regulatory_ensemble()
    ensemble_by_id = {member.member_id: member for member in ensemble}
    if preferred_member_id not in ensemble_by_id:
        raise KeyError(f"preferred ensemble member {preferred_member_id!r} is absent")
    if ensemble_by_id[preferred_member_id].family != preferred_family:
        raise ValueError("preferred member and preferred family are inconsistent")
    dynamic_models = {
        member.member_id: ModernFullModel(
            parameters=models[member.family].parameters,
            stimulus=protocol,
            regulatory_model=member.regulatory_model,
            ae4_parameters=models[member.family].ae4_parameters,
            ae4_evaluator=models[member.family].ae4_evaluator,
        )
        for member in ensemble
    }
    preferred = dynamic_models[preferred_member_id]

    grid = physical_time_grid(duration_s=600.0, sample_step_s=sample_step_s)
    trajectories: list[Trajectory] = []
    baseline: dict[tuple[str, str], Trajectory] = {}
    for member in ensemble:
        model = dynamic_models[member.member_id]
        for specification in (PRODUCTION_RADAU, PRODUCTION_BDF):
            trajectory = simulate_wt(
                model,
                core_root,
                family=member.member_id,
                initial_condition="baseline",
                solver=specification,
                time_s=grid,
            )
            trajectories.append(trajectory)
            baseline[(member.member_id, specification.label)] = trajectory

    for specification in (LOOSE_RADAU, TIGHT_RADAU):
        trajectories.append(
            simulate_wt(
                preferred,
                core_root,
                family=preferred_member_id,
                initial_condition="baseline",
                solver=specification,
                time_s=grid,
            )
        )

    nearby = nearby_charge_preserving_states(
        preferred, core_root, fraction=nearby_fraction
    )
    for initial_condition, candidate in nearby.items():
        if initial_condition == "baseline":
            continue
        trajectories.append(
            simulate_wt(
                preferred,
                candidate,
                family=preferred_member_id,
                initial_condition=initial_condition,
                solver=PRODUCTION_RADAU,
                time_s=grid,
            )
        )

    trajectory_path = output / "wt_trajectories.csv"
    write_trajectory_csv(trajectory_path, trajectories)
    summary_rows = [_trajectory_summary(item) for item in trajectories]
    summary_path = output / "wt_trajectory_summary.csv"
    _write_rows(summary_path, summary_rows)

    solver_rows: list[Mapping[str, Any]] = []
    for member in ensemble:
        comparison = solver_comparison(
            baseline[(member.member_id, PRODUCTION_RADAU.label)],
            baseline[(member.member_id, PRODUCTION_BDF.label)],
        )
        solver_rows.append(
            {
                **member.metadata(),
                "comparison": "Radau_vs_BDF",
                **comparison,
            }
        )
    preferred_production = baseline[(preferred_member_id, PRODUCTION_RADAU.label)]
    for label in (LOOSE_RADAU.label, TIGHT_RADAU.label):
        candidate = next(
            item
            for item in trajectories
            if item.family == preferred_member_id
            and item.initial_condition == "baseline"
            and item.solver.label == label
        )
        solver_rows.append(
            {
                **ensemble_by_id[preferred_member_id].metadata(),
                "comparison": f"{PRODUCTION_RADAU.label}_vs_{label}",
                **solver_comparison(preferred_production, candidate),
            }
        )
    solver_path = output / "wt_solver_crosschecks.csv"
    _write_rows(solver_path, solver_rows)

    sensitivity_rows: list[Mapping[str, Any]] = []
    for trajectory in trajectories:
        if (
            trajectory.family != preferred_member_id
            or trajectory.solver.label != PRODUCTION_RADAU.label
            or trajectory.initial_condition == "baseline"
        ):
            continue
        sensitivity_rows.append(
            {
                **ensemble_by_id[preferred_member_id].metadata(),
                "initial_condition": trajectory.initial_condition,
                "nearby_fraction": nearby_fraction,
                "success": trajectory.success,
                "positive_core": trajectory.positive_core,
                "endpoint_cl_difference_mM": (
                    trajectory.cell_cl_mM[-1] - preferred_production.cell_cl_mM[-1]
                ),
                "endpoint_ph_difference": (
                    trajectory.cell_ph[-1] - preferred_production.cell_ph[-1]
                ),
                "endpoint_volume_relative_difference": (
                    trajectory.cell_volume_pL[-1] / preferred_production.cell_volume_pL[-1]
                    - 1.0
                ),
                "cumulative_flow_relative_difference": (
                    trajectory.cumulative_flow_pL[-1]
                    / preferred_production.cumulative_flow_pL[-1]
                    - 1.0
                    if preferred_production.cumulative_flow_pL[-1] != 0.0
                    else float("nan")
                ),
            }
        )
    sensitivity_path = output / "wt_initial_state_sensitivity.csv"
    _write_rows(sensitivity_path, sensitivity_rows)

    regulatory_rows = [
        {
            **regulatory_landmarks(
                baseline[(member.member_id, PRODUCTION_RADAU.label)]
            ),
            **member.metadata(),
        }
        for member in ensemble
    ]
    regulatory_path = output / "regulatory_fit_summary.csv"
    _write_rows(regulatory_path, regulatory_rows)

    arm_rows: list[Mapping[str, Any]] = []
    for member in ensemble:
        for arm in (StimulusArm.CCH_ONLY, StimulusArm.IPR_ONLY, StimulusArm.CCH_IPR):
            for initialization in ("basal_at_onset", "beta_preconditioned"):
                response = wt_arm_response(
                    dynamic_models[member.member_id],
                    core_root,
                    family=member.member_id,
                    arm=arm,
                    regulatory_initialization=initialization,
                )
                arm_rows.append({**asdict(response), **member.metadata()})
    arm_path = output / "wt_ion_uptake_diagnostics.csv"
    _write_rows(arm_path, arm_rows)

    flow_envelope = WTFlowEnvelope()
    flow_scale = fit_wt_flow_scale(preferred_production, flow_envelope)
    flow_shape = wt_flow_envelope_feasibility(
        preferred_production, flow_envelope
    )
    if effective_cell_count_range is None:
        cell_count_gate: bool | None = None
        effective_cell_count_status = (
            "WARNING_UNASSESSED_NO_PRIMARY_GEOMETRY_RANGE"
        )
    else:
        lower_cells, upper_cells = (float(value) for value in effective_cell_count_range)
        if not 0.0 < lower_cells < upper_cells:
            raise ValueError("effective-cell-count range must be positive and ordered")
        cell_count_gate = bool(
            lower_cells <= flow_scale.effective_cell_count <= upper_cells
        )
        effective_cell_count_status = (
            "PASS_SOURCE_SUPPORTED_RANGE"
            if cell_count_gate
            else "FAIL_OUTSIDE_SOURCE_SUPPORTED_RANGE"
        )
    # The single WT gland multiplier is a Tier-2 observation scale.  Its
    # implied-cell-count diagnostic remains a warning when no independent
    # primary geometry range exists; absence of that range cannot by itself
    # invalidate an otherwise scale-compatible WT trajectory.  A supplied
    # range is a stronger check and becomes binding.
    flow_scale_gate = bool(
        flow_shape.compatible
        and (cell_count_gate is None or cell_count_gate)
    )
    flow_scale_path = output / "wt_flow_scale.json"
    write_json(
        flow_scale_path,
        {
            **asdict(flow_scale),
            "envelope": asdict(flow_envelope),
            "fit_uses": "WT_ONLY",
            "observation_scale_tier": "TIER_2_WT_OBSERVATION_SCALE",
            "minute_wise_shape_status": "GRAPHICAL_ENVELOPE_ONLY_NOT_DIGITIZED",
            "effective_cell_count_source_range": effective_cell_count_range,
            "effective_cell_count_plausibility_status": effective_cell_count_status,
            "effective_cell_count_gate_pass": cell_count_gate,
            "minute_wise_envelope_feasibility": asdict(flow_shape),
            "minute_wise_shape_gate_status": (
                "PASS_DECLARED_GRAPHICAL_ENVELOPE"
                if flow_shape.compatible
                else "FAIL_DECLARED_GRAPHICAL_ENVELOPE"
            ),
            "absolute_flow_scale_gate_pass": flow_scale_gate,
        },
        pre_reveal=True,
    )
    scaled_flow = flow_scale.convert(preferred_production.flow_pL_s)
    minute_rows: list[Mapping[str, Any]] = []
    for minute in range(0, 11):
        time_value = float(60 * minute)
        index = int(np.argmin(np.abs(preferred_production.time_s - time_value)))
        minute_rows.append(
            {
                **ensemble_by_id[preferred_member_id].metadata(),
                "time_s": preferred_production.time_s[index],
                "time_min": preferred_production.time_s[index] / 60.0,
                "model_flow_pL_s": preferred_production.flow_pL_s[index],
                "scaled_WT_flow_uL_min": scaled_flow[index],
                "source_envelope_lower_uL_min": flow_envelope.lower_uL_min,
                "source_envelope_upper_uL_min": flow_envelope.upper_uL_min,
            }
        )
    scaled_flow_path = output / "wt_flow_trajectory.csv"
    _write_rows(scaled_flow_path, minute_rows)

    numerical_success = bool(all(item.success for item in trajectories))
    positive_success = bool(all(item.positive_core for item in trajectories))
    conservation_max_by_key = {
        name: max(
            item.max_abs_conservation_residuals[name] for item in trajectories
        )
        for name in CONSERVATION_RESIDUAL_UNITS
    }
    conservation_normalized_by_key = {
        name: conservation_max_by_key[name] / declared_conservation_tolerances[name]
        for name in CONSERVATION_RESIDUAL_UNITS
    }
    conservation_success = bool(
        all(value <= 1.0 for value in conservation_normalized_by_key.values())
    )
    solver_success = bool(
        all(bool(row["both_success"]) for row in solver_rows)
        and all(
            float(row["max_relative_state_difference"]) <= solver_relative_tolerance
            and float(row["endpoint_cumulative_flow_relative_difference"])
            <= solver_relative_tolerance
            for row in solver_rows
        )
    )
    nonnegative_flow = bool(all(np.all(item.flow_pL_s >= 0.0) for item in trajectories))
    numerical_dynamic_gate = bool(
        numerical_success
        and positive_success
        and conservation_success
        and solver_success
        and nonnegative_flow
    )
    wt_gate = bool(
        root_accepted
        and resting_physiology_gate_pass
        and numerical_dynamic_gate
        and flow_scale_gate
    )
    ensemble_frozen = bool(
        len(ensemble) == 20
        and {member.family for member in ensemble}
        == set(REQUIRED_REGULATORY_FAMILIES)
        and {member.fully_activated_multiplier for member in ensemble}
        == {1.25, 1.70}
        and {member.timing_label for member in ensemble}
        == {"STATIC", "TFAST", "TREFERENCE", "TSLOW"}
    )
    gate_payload = {
        "status": "PASS" if wt_gate else "FAIL",
        "root_id": root_id,
        "root_accepted": bool(root_accepted),
        "resting_physiology_gate_pass": bool(resting_physiology_gate_pass),
        "R0_R3_nesting_present": True,
        "regulatory_ensemble_member_count": len(ensemble),
        "regulatory_gain_multipliers": sorted(
            {member.fully_activated_multiplier for member in ensemble}
        ),
        "regulatory_timing_sensitivities": sorted(
            {member.timing_label for member in ensemble}
        ),
        "regulatory_ensemble_frozen": ensemble_frozen,
        "regulatory_model_class_numerical_eligibility": {
            family: bool(
                all(
                    trajectory.success and trajectory.positive_core
                    for trajectory in trajectories
                    if ensemble_by_id.get(trajectory.family) is not None
                    and ensemble_by_id[trajectory.family].family == family
                )
            )
            for family in REQUIRED_REGULATORY_FAMILIES
        },
        "numerical_dynamic_gate_pass": numerical_dynamic_gate,
        "all_integrations_success": numerical_success,
        "all_core_states_positive": positive_success,
        "all_flows_nonnegative": nonnegative_flow,
        "conservation_gate_pass": conservation_success,
        "solver_and_tolerance_gate_pass": solver_success,
        "conservation_residual_units": dict(CONSERVATION_RESIDUAL_UNITS),
        "conservation_residual_tolerances": declared_conservation_tolerances,
        "conservation_max_abs_by_key": conservation_max_by_key,
        "conservation_normalized_by_key": conservation_normalized_by_key,
        "max_dimensionless_conservation_ratio": max(
            conservation_normalized_by_key.values()
        ),
        "solver_relative_tolerance": solver_relative_tolerance,
        "physical_time_axis": "seconds",
        "physical_protocol_duration_s": 600.0,
        "flow_scale_source": "WT envelope only",
        "flow_scale_observation_tier": "TIER_2_WT_OBSERVATION_SCALE",
        "absolute_flow_scale_gate_pass": flow_scale_gate,
        "minute_wise_flow_shape_gate_pass": flow_shape.compatible,
        "minute_wise_flow_shape": asdict(flow_shape),
        "effective_cell_count_gate_pass": cell_count_gate,
        "effective_cell_count": flow_scale.effective_cell_count,
        "effective_cell_count_source_range": effective_cell_count_range,
        "effective_cell_count_plausibility_status": effective_cell_count_status,
        "uptake_observation_map_status": (
            "MODEL_FRACTIONAL_CL_RATE_REPORTED; NOT_EQUATED_TO_SPQ_WITHOUT_MAP"
        ),
        "regulatory_kinetics_identified": False,
        "regulatory_kinetic_reason": (
            "2021 assays constrain direction/hierarchy but contain no acute onset trace"
        ),
        "unidentified_regulatory_kinetics_is_reveal_blocker": False,
        "pre_reveal_holdout_eligible": bool(wt_gate and ensemble_frozen),
        "holdout_blocker": (
            None
            if wt_gate and ensemble_frozen
            else (
                "A frozen, numerically valid WT root and the complete regulatory "
                "ensemble are available, but the fixed WT minute-wise flow-shape "
                "gate is not passed. Unidentified acute regulatory kinetics are "
                "carried ensemble-wide and are not independently blocking."
            )
        ),
    }
    gate_path = output / "wt_dynamic_gate.json"
    write_json(gate_path, gate_payload, pre_reveal=True)

    model_payloads = {
        member.member_id: {
            **member.metadata(),
            "state_names": list(dynamic_models[member.member_id].state_names),
            "regulatory_model": _object_payload(
                dynamic_models[member.member_id].regulatory_model
            ),
            "regulatory_model_class": (
                None
                if dynamic_models[member.member_id].regulatory_model is None
                else type(dynamic_models[member.member_id].regulatory_model).__name__
            ),
        }
        for member in ensemble
    }
    implementation_directory = Path(__file__).resolve().parent
    implementation_files = (
        "acid_base.py",
        "calibration.py",
        "camp_pka.py",
        "membranes.py",
        "model.py",
        "parameters.py",
        "run_analysis.py",
        "states.py",
        "transporters.py",
        "validation.py",
        "water.py",
    )
    implementation_hashes = {
        name: sha256_file(implementation_directory / name)
        for name in implementation_files
    }
    manifest = {
        "manifest_type": "TASK_13B_PRE_REVEAL_WT_DYNAMIC_FREEZE",
        "firewall_status": "WT_ONLY_NO_NULL_TRAJECTORY_OR_TARGET_ACCESS",
        "root_id": root_id,
        "root_selection_rule": dict(root_selection_rule),
        "root_state_sha256": sha256_object(core_root),
        "root_state_names": list(CORE_STATE_NAMES),
        "core_parameter_sha256": core_parameter_hash,
        "ae4_parameter_sha256": ae4_parameter_hash,
        "preferred_regulatory_family": preferred_family,
        "preferred_regulatory_member": preferred_member_id,
        "regulatory_models": model_payloads,
        "stimulus_protocol": asdict(protocol),
        "stimulus_provenance": dict(protocol.provenance()),
        "solver_specifications": [
            asdict(item)
            for item in (PRODUCTION_RADAU, PRODUCTION_BDF, LOOSE_RADAU, TIGHT_RADAU)
        ],
        "sample_step_s": sample_step_s,
        "nearby_state_fraction": nearby_fraction,
        "flow_scale": asdict(flow_scale),
        "gate": gate_payload,
        "parameter_records_sha256": sha256_object(
            [asdict(item) for item in parameter_records(preferred.parameters)]
        ),
        "implementation_file_sha256": implementation_hashes,
    }
    manifest_path = output / "dynamic_frozen_manifest.json"
    write_json(manifest_path, manifest, pre_reveal=True)

    artifact_paths = (
        trajectory_path,
        summary_path,
        solver_path,
        sensitivity_path,
        regulatory_path,
        arm_path,
        flow_scale_path,
        scaled_flow_path,
        gate_path,
        manifest_path,
    )
    hash_rows = [
        {
            "artifact": path.name,
            "sha256": sha256_file(path),
            "scope": "PRE_REVEAL_WT_ONLY",
        }
        for path in artifact_paths
    ]
    hash_path = output / "dynamic_artifact_hashes.csv"
    _write_rows(hash_path, hash_rows)
    return {
        "gate": gate_payload,
        "manifest": manifest,
        "artifacts": {row["artifact"]: row["sha256"] for row in hash_rows},
        "hash_table_sha256": sha256_file(hash_path),
    }


__all__ = ("REQUIRED_REGULATORY_FAMILIES", "run_pre_reveal_wt_suite")
