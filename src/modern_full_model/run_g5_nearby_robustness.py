"""Compact H23 nearby-state check for selected-root joint G5 candidates.

The check uses the preferred R1 regulatory member and both production solver
algorithms.  It perturbs cell and lumen volume plus equal Na/Cl amount pairs
by 0.1%, preserving bulk charge.  The common ensemble-wide flow scale frozen
by the joint profile is reused without refitting.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .calibration import build_wt_model
from .model import ModernFullModel
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .run_g5_cross_root_profile import load_retained_wt_roots
from .run_g5_joint_repair_profile import (
    PREFERRED_REGULATORY_MEMBER,
    declared_joint_repair_profiles,
)
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    nearby_charge_preserving_states,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    sha256_file,
    simulate_wt,
    wt_flow_envelope_feasibility,
)


PROFILE_FILENAME = "wt_g5_final_candidate_nearby_robustness.csv"
SUMMARY_FILENAME = "wt_g5_final_candidate_nearby_summary.json"
HASH_FILENAME = "wt_g5_final_candidate_nearby_hashes.csv"
NEARBY_FRACTION = 1.0e-3
OUTPUT_RELATIVE_WARNING = 1.0e-2


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for name in row:
            if name not in fields:
                fields.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _read_group_scales(results_directory: Path) -> Mapping[str, float]:
    path = results_directory / "wt_g5_joint_ca_nkcc_group_gate.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = [
        row
        for row in rows
        if row["selected_reference"] == "True"
        and row["g5_joint_group_gate_pass"] == "True"
    ]
    if not selected:
        raise ValueError("no selected-root joint G5 candidate passed before H23")
    return {
        row["profile_id"]: float(row["common_scale_multiplier_uL_min_per_pL_s"])
        for row in selected
    }


def run_nearby_robustness(results_directory: Path) -> Mapping[str, Any]:
    root = next(
        item
        for item in load_retained_wt_roots(results_directory)
        if item.selected_reference
    )
    member = next(
        item
        for item in pre_reveal_regulatory_ensemble()
        if item.member_id == PREFERRED_REGULATORY_MEMBER
    )
    profiles = {
        item.profile_id: item for item in declared_joint_repair_profiles()
    }
    scales = _read_group_scales(results_directory)
    rows: list[dict[str, Any]] = []
    for profile_id, common_scale in sorted(scales.items()):
        profile = profiles[profile_id]
        resting = build_wt_model(
            root.variant,
            root.calibration,
            regulatory_model=member.regulatory_model,
        )
        protocol = SecretagogueProtocol(
            stimulated_calcium_uM=profile.effective_calcium_uM
        )
        template = ModernFullModel(
            parameters=resting.parameters,
            stimulus=protocol,
            regulatory_model=member.regulatory_model,
            ae4_parameters=resting.ae4_parameters,
            ae4_evaluator=resting.ae4_evaluator,
        )
        model = attach_stimulated_nkcc1(
            template,
            N1AlgebraicNkcc1(
                fully_activated_multiplier=(
                    profile.nkcc1_fully_activated_multiplier
                ),
                resting_calcium_uM=protocol.resting_calcium_uM,
                stimulated_calcium_uM=protocol.stimulated_calcium_uM,
            ),
        )
        nearby = nearby_charge_preserving_states(
            model, root.core_state, fraction=NEARBY_FRACTION
        )
        for solver in (PRODUCTION_RADAU, PRODUCTION_BDF):
            trajectories = {
                label: simulate_wt(
                    model,
                    state,
                    family=f"{member.member_id}+{profile_id}",
                    initial_condition=label,
                    solver=solver,
                    time_s=physical_time_grid(600.0, 5.0),
                )
                for label, state in nearby.items()
            }
            baseline = trajectories["baseline"]
            for label, trajectory in trajectories.items():
                shape = wt_flow_envelope_feasibility(trajectory)
                minute_times = np.arange(60.0, 601.0, 60.0)
                minute_flow = np.interp(
                    minute_times, trajectory.time_s, trajectory.flow_pL_s
                )
                scaled = common_scale * minute_flow
                conservation_ratio = max(
                    trajectory.max_abs_conservation_residuals[name] / tolerance
                    for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
                )
                cumulative_relative = (
                    trajectory.cumulative_flow_pL[-1]
                    / baseline.cumulative_flow_pL[-1]
                    - 1.0
                )
                endpoint_flow_relative = (
                    trajectory.flow_pL_s[-1] / baseline.flow_pL_s[-1] - 1.0
                )
                output_warning = bool(
                    abs(cumulative_relative) > OUTPUT_RELATIVE_WARNING
                    or abs(endpoint_flow_relative) > OUTPUT_RELATIVE_WARNING
                )
                common_scale_gate = bool(
                    np.min(scaled) >= 9.0 and np.max(scaled) <= 10.0
                )
                numerical_gate = bool(
                    trajectory.success
                    and trajectory.positive_core
                    and np.all(trajectory.flow_pL_s >= 0.0)
                    and conservation_ratio <= 1.0
                )
                rows.append(
                    {
                        "root_id": root.root_id,
                        "profile_id": profile_id,
                        "regulatory_member_id": member.member_id,
                        "solver_label": solver.label,
                        "initial_condition": label,
                        "nearby_fraction": NEARBY_FRACTION,
                        "success": trajectory.success,
                        "positive_core": trajectory.positive_core,
                        "max_dimensionless_conservation_ratio": conservation_ratio,
                        "numerical_gate_pass": numerical_gate,
                        "individual_shape_gate_pass": shape.compatible,
                        "common_scale_multiplier_uL_min_per_pL_s": common_scale,
                        "common_scaled_minute_min_uL_min": float(np.min(scaled)),
                        "common_scaled_minute_max_uL_min": float(np.max(scaled)),
                        "common_scale_gate_pass": common_scale_gate,
                        "endpoint_cl_difference_mM": (
                            trajectory.cell_cl_mM[-1] - baseline.cell_cl_mM[-1]
                        ),
                        "endpoint_ph_difference": (
                            trajectory.cell_ph[-1] - baseline.cell_ph[-1]
                        ),
                        "endpoint_volume_relative_difference": (
                            trajectory.cell_volume_pL[-1]
                            / baseline.cell_volume_pL[-1]
                            - 1.0
                        ),
                        "endpoint_flow_relative_difference": endpoint_flow_relative,
                        "cumulative_flow_relative_difference": cumulative_relative,
                        "output_relative_warning_threshold": OUTPUT_RELATIVE_WARNING,
                        "output_warning": output_warning,
                        "nearby_gate_pass": bool(
                            numerical_gate
                            and shape.compatible
                            and common_scale_gate
                            and not output_warning
                        ),
                    }
                )
    profile_path = results_directory / PROFILE_FILENAME
    _write_rows(profile_path, rows)
    profile_gates = {
        profile_id: all(
            bool(row["nearby_gate_pass"])
            for row in rows
            if row["profile_id"] == profile_id
        )
        for profile_id in scales
    }
    summary = {
        "analysis": "H23_SELECTED_ROOT_CELL_AND_LUMEN_NEARBY_ROBUSTNESS",
        "root_id": root.root_id,
        "candidate_profiles": sorted(scales),
        "regulatory_member_id": member.member_id,
        "solver_algorithms": [PRODUCTION_RADAU.method, PRODUCTION_BDF.method],
        "nearby_fraction": NEARBY_FRACTION,
        "nearby_initial_conditions": sorted(
            {row["initial_condition"] for row in rows}
        ),
        "charge_preservation": (
            "ENFORCED_BY_nearby_charge_preserving_states_FOR_CELL_AND_LUMEN"
        ),
        "output_relative_warning_threshold": OUTPUT_RELATIVE_WARNING,
        "output_threshold_status": "MODELING_ROBUSTNESS_DIAGNOSTIC",
        "profile_gate_pass": profile_gates,
        "all_candidate_nearby_gates_pass": all(profile_gates.values()),
        "row_count": len(rows),
        "profile_sha256": sha256_file(profile_path),
        "input_group_gate_sha256": sha256_file(
            results_directory / "wt_g5_joint_ca_nkcc_group_gate.csv"
        ),
    }
    summary_path = results_directory / SUMMARY_FILENAME
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    hash_rows = [
        {"artifact": path.name, "sha256": sha256_file(path), "scope": "WT_ONLY"}
        for path in (profile_path, summary_path)
    ]
    hash_path = results_directory / HASH_FILENAME
    _write_rows(hash_path, hash_rows)
    return {
        "summary": summary,
        "artifact_hashes": {row["artifact"]: row["sha256"] for row in hash_rows},
        "hash_table_sha256": sha256_file(hash_path),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    args = parser.parse_args(argv)
    result = run_nearby_robustness(args.results_directory)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ("run_nearby_robustness",)
