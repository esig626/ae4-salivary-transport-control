"""WT-only R0--R3 sustainment sensitivity on the provisional G2b root.

This script changes only the predeclared regulatory family, gain reading, and
unmeasured timing sensitivity.  The whole-cell parameters, AE4 carrier model,
initial core state, calcium/beta protocol, time grid, solver, and WT graphical
flow envelope are fixed.  It does not load any phenotype target file.
"""

from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from modern_full_model.calibration import (
    ModelVariant,
    WTCalibrationParameters,
    build_wt_model,
)
from modern_full_model.model import ModernFullModel
from modern_full_model.validation import (
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    WTFlowEnvelope,
    fit_wt_flow_scale,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    regulatory_landmarks,
    sha256_file,
    sha256_object,
    simulate_wt,
    solver_comparison,
    wt_flow_envelope_feasibility,
)


EXPECTED_ROOT_ID = "G2_BALANCED_APICAL_K_BIASED_G10:M00"


CSV_FIELDS = (
    "member_id",
    "family",
    "gain_label",
    "fully_activated_multiplier",
    "timing_label",
    "timing_parameters_json",
    "kinetic_evidence_status",
    "regulatory_model_class",
    "success",
    "positive_core",
    "capacity_t50_s",
    "capacity_t90_s",
    "endpoint_capacity_multiplier",
    "flow_60_pL_s",
    "flow_600_pL_s",
    "flow_600_over_60",
    "flow_decline_percent",
    "minute_min_flow_pL_s",
    "minute_max_flow_pL_s",
    "minute_max_min_ratio",
    "envelope_ratio_limit",
    "feasible_single_scale",
    "feasible_multiplier_lower_uL_min_per_pL_s",
    "feasible_multiplier_upper_uL_min_per_pL_s",
    "scaled_min_uL_min",
    "scaled_max_uL_min",
    "minutes_in_9_10_envelope",
    "mean_scale_effective_cell_count",
    "endpoint_cl_mM",
    "endpoint_ph",
    "endpoint_volume_pL",
    "endpoint_flow_pL_s",
    "cumulative_flow_pL",
)


def _load_selected_construction(
    selected_path: Path,
) -> tuple[Mapping[str, Any], ModelVariant, WTCalibrationParameters, np.ndarray]:
    selected = json.loads(selected_path.read_text(encoding="utf-8"))
    construction = selected["construction"]
    if construction["root_id"] != EXPECTED_ROOT_ID:
        raise ValueError(
            f"expected frozen staging root {EXPECTED_ROOT_ID}, got {construction['root_id']}"
        )
    return (
        selected,
        ModelVariant(**construction["variant"]),
        WTCalibrationParameters(**construction["calibration"]),
        np.asarray(construction["root_state"], dtype=float),
    )


def _trajectory(
    member: Any,
    grid: np.ndarray,
    solver: Any,
    variant: ModelVariant,
    calibration: WTCalibrationParameters,
    core_state: np.ndarray,
) -> Any:
    resting = build_wt_model(
        variant,
        calibration,
        regulatory_model=member.regulatory_model,
    )
    model = ModernFullModel(
        parameters=resting.parameters,
        stimulus=SecretagogueProtocol(),
        regulatory_model=member.regulatory_model,
        ae4_parameters=resting.ae4_parameters,
    )
    return simulate_wt(
        model,
        core_state,
        family=member.member_id,
        initial_condition="baseline",
        solver=solver,
        time_s=grid,
    )


def _one_member(
    member: Any,
    grid: np.ndarray,
    variant: ModelVariant,
    calibration: WTCalibrationParameters,
    core_state: np.ndarray,
) -> dict[str, Any]:
    trajectory = _trajectory(
        member,
        grid,
        PRODUCTION_RADAU,
        variant,
        calibration,
        core_state,
    )
    feasibility = wt_flow_envelope_feasibility(trajectory)
    scale = fit_wt_flow_scale(trajectory)
    minute_times = np.arange(60.0, 601.0, 60.0)
    minute_flow = np.interp(
        minute_times,
        trajectory.time_s,
        trajectory.flow_pL_s,
    )
    scaled_minute_flow = scale.convert(minute_flow)
    landmarks = regulatory_landmarks(trajectory)
    metadata = dict(member.metadata())
    timing_parameters = metadata.pop("timing_parameters", {})
    return {
        **metadata,
        "timing_parameters_json": json.dumps(
            timing_parameters, sort_keys=True
        ),
        "success": trajectory.success,
        "positive_core": trajectory.positive_core,
        "capacity_t50_s": landmarks["t50_s"],
        "capacity_t90_s": landmarks["t90_s"],
        "endpoint_capacity_multiplier": float(
            trajectory.capacity_multiplier[-1]
        ),
        "flow_60_pL_s": float(minute_flow[0]),
        "flow_600_pL_s": float(minute_flow[-1]),
        "flow_600_over_60": float(minute_flow[-1] / minute_flow[0]),
        "flow_decline_percent": float(
            100.0 * (1.0 - minute_flow[-1] / minute_flow[0])
        ),
        "minute_min_flow_pL_s": float(np.min(minute_flow)),
        "minute_max_flow_pL_s": float(np.max(minute_flow)),
        "minute_max_min_ratio": feasibility.model_max_min_ratio,
        "envelope_ratio_limit": feasibility.envelope_max_min_ratio,
        "feasible_single_scale": feasibility.compatible,
        "feasible_multiplier_lower_uL_min_per_pL_s": getattr(
            feasibility,
            "feasible_multiplier_lower_uL_min_per_pL_s",
        ),
        "feasible_multiplier_upper_uL_min_per_pL_s": getattr(
            feasibility,
            "feasible_multiplier_upper_uL_min_per_pL_s",
        ),
        "scaled_min_uL_min": float(np.min(scaled_minute_flow)),
        "scaled_max_uL_min": float(np.max(scaled_minute_flow)),
        "minutes_in_9_10_envelope": int(
            np.sum((scaled_minute_flow >= 9.0) & (scaled_minute_flow <= 10.0))
        ),
        "mean_scale_effective_cell_count": scale.effective_cell_count,
        "endpoint_cl_mM": float(trajectory.cell_cl_mM[-1]),
        "endpoint_ph": float(trajectory.cell_ph[-1]),
        "endpoint_volume_pL": float(trajectory.cell_volume_pL[-1]),
        "endpoint_flow_pL_s": float(trajectory.flow_pL_s[-1]),
        "cumulative_flow_pL": float(trajectory.cumulative_flow_pL[-1]),
    }


def _paired_r2_r3_difference(
    rows: list[Mapping[str, Any]], field: str
) -> float:
    maximum = 0.0
    for r2 in (row for row in rows if row["family"] == "R2"):
        partner = next(
            row
            for row in rows
            if row["family"] == "R3"
            and row["gain_label"] == r2["gain_label"]
            and row["timing_label"] == r2["timing_label"]
        )
        maximum = max(maximum, abs(float(r2[field]) - float(partner[field])))
    return maximum


def run(output_directory: Path) -> Mapping[str, Any]:
    selected_path = output_directory / "wt_selected_candidate.json"
    selected, variant, calibration, core_state = _load_selected_construction(
        selected_path
    )
    ensemble = pre_reveal_regulatory_ensemble()
    if len(ensemble) != 20:
        raise AssertionError("the frozen regulatory ensemble must contain 20 members")
    grid = physical_time_grid(duration_s=600.0, sample_step_s=5.0)
    rows = [
        _one_member(member, grid, variant, calibration, core_state)
        for member in ensemble
    ]
    dynamic_rows = [row for row in rows if row["family"] != "R0"]
    if len(dynamic_rows) != 18:
        raise AssertionError("the frozen ensemble must contain 18 dynamic members")

    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / "regulatory_sustainment_sensitivity.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    envelope = WTFlowEnvelope()
    best = min(dynamic_rows, key=lambda row: float(row["minute_max_min_ratio"]))
    worst = max(dynamic_rows, key=lambda row: float(row["minute_max_min_ratio"]))
    best_ratio = float(best["minute_max_min_ratio"])
    worst_ratio = float(worst["minute_max_min_ratio"])
    best_equivalence_class = sorted(
        str(row["member_id"])
        for row in dynamic_rows
        if abs(float(row["minute_max_min_ratio"]) - best_ratio) <= 1.0e-9
    )
    worst_equivalence_class = sorted(
        str(row["member_id"])
        for row in dynamic_rows
        if abs(float(row["minute_max_min_ratio"]) - worst_ratio) <= 1.0e-9
    )
    members_by_id = {member.member_id: member for member in ensemble}
    solver_crosschecks: dict[str, Any] = {}
    for label, row in (("best_dynamic_member", best), ("worst_dynamic_member", worst)):
        member = members_by_id[str(row["member_id"])]
        radau = _trajectory(
            member, grid, PRODUCTION_RADAU, variant, calibration, core_state
        )
        bdf = _trajectory(
            member, grid, PRODUCTION_BDF, variant, calibration, core_state
        )
        solver_crosschecks[label] = solver_comparison(radau, bdf)
    gain_summaries: dict[str, Any] = {}
    for gain_label in sorted({str(row["gain_label"]) for row in dynamic_rows}):
        subset = [row for row in dynamic_rows if row["gain_label"] == gain_label]
        gain_summaries[gain_label] = {
            "member_count": len(subset),
            "best_minute_max_min_ratio": min(
                float(row["minute_max_min_ratio"]) for row in subset
            ),
            "worst_minute_max_min_ratio": max(
                float(row["minute_max_min_ratio"]) for row in subset
            ),
            "flow_decline_percent_range": [
                min(float(row["flow_decline_percent"]) for row in subset),
                max(float(row["flow_decline_percent"]) for row in subset),
            ],
        }
    summary = {
        "analysis": "WT_ONLY_PREDECLARED_REGULATORY_SUSTAINMENT_SENSITIVITY",
        "classification": (
            "NO_PREDECLARED_R1_R3_MEMBER_REPAIRS_WT_10_MIN_FLOW_SHAPE_"
            "ON_PROVISIONAL_CORRECTED_WATER_ROOT"
        ),
        "root_label": EXPECTED_ROOT_ID,
        "root_status": selected["status"],
        "dynamic_use_status": "PROVISIONAL_G5_STAGING_NOT_PRODUCTION_FREEZE",
        "selected_candidate_sha256": sha256_file(selected_path),
        "variant": asdict(variant),
        "calibration": asdict(calibration),
        "core_state_sha256": sha256_object(core_state),
        "regulatory_ensemble_sha256": sha256_object(
            [member.metadata() for member in ensemble]
        ),
        "time_grid_sha256": sha256_object(grid),
        "protocol": asdict(SecretagogueProtocol()),
        "solver": asdict(PRODUCTION_RADAU),
        "regulatory_member_count": len(rows),
        "dynamic_member_count": len(dynamic_rows),
        "all_integrations_success": all(bool(row["success"]) for row in rows),
        "all_core_states_positive": all(bool(row["positive_core"]) for row in rows),
        "members_passing_single_scale_envelope": sum(
            bool(row["feasible_single_scale"]) for row in rows
        ),
        "dynamic_members_passing_single_scale_envelope": sum(
            bool(row["feasible_single_scale"]) for row in dynamic_rows
        ),
        "all_members_minutes_in_envelope_range": [
            min(int(row["minutes_in_9_10_envelope"]) for row in rows),
            max(int(row["minutes_in_9_10_envelope"]) for row in rows),
        ],
        "declared_envelope_max_min_ratio": (
            envelope.upper_uL_min / envelope.lower_uL_min
        ),
        "best_dynamic_member": dict(best),
        "worst_dynamic_member": dict(worst),
        "best_dynamic_equivalence_class_at_1e_9_ratio_tolerance": (
            best_equivalence_class
        ),
        "worst_dynamic_equivalence_class_at_1e_9_ratio_tolerance": (
            worst_equivalence_class
        ),
        "dynamic_ratio_range": [
            float(best["minute_max_min_ratio"]),
            float(worst["minute_max_min_ratio"]),
        ],
        "best_ratio_excess_over_envelope_fraction": (
            float(best["minute_max_min_ratio"])
            / (envelope.upper_uL_min / envelope.lower_uL_min)
            - 1.0
        ),
        "best_late_early_flow_ratio": float(best["flow_600_over_60"]),
        "minimum_late_early_ratio_required_by_envelope": (
            envelope.lower_uL_min / envelope.upper_uL_min
        ),
        "gain_summaries": gain_summaries,
        "max_abs_matched_r2_r3_ratio_difference": _paired_r2_r3_difference(
            rows, "minute_max_min_ratio"
        ),
        "max_abs_matched_r2_r3_cumulative_flow_difference_pL": (
            _paired_r2_r3_difference(rows, "cumulative_flow_pL")
        ),
        "best_and_worst_radau_bdf_crosschecks": solver_crosschecks,
        "decision": (
            "COMMON_CAPACITY_REGULATION_IS_NOT_THE_G5_SUSTAINMENT_REPAIR_FOR_"
            "THIS_PROVISIONAL_ROOT"
        ),
        "interpretation": (
            "The complete predeclared 2021-compatible gain/timing ensemble is "
            "numerically valid but every member has an empty WT minute-envelope "
            "scale intersection. Slower activation and the larger gain do not "
            "supply sustained flow; the larger gain worsens the best shape ratio. "
            "The next WT reconstruction must change a nonregulatory flux, water, "
            "outflow, or calibration axis rather than tune an unmeasured regulator."
        ),
    }
    summary_path = output_directory / "regulatory_sustainment_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary


if __name__ == "__main__":
    repository_root = Path(__file__).resolve().parents[2]
    payload = run(repository_root / "results" / "13B_modern_full_model")
    print(json.dumps(payload, indent=2, sort_keys=True))
