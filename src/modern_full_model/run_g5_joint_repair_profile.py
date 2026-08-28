"""Ensemble-wide WT G5 profile for joint effective-Ca/NKCC1 repairs.

The matrix is frozen before phenotype reveal and uses only the nine retained
WT resting roots, WT flow/ion diagnostics, the historical calcium-input
lineage, P15 NKCC1 direction/scale evidence, and the mandatory P21-compatible
R0--R3 AE4 regulatory ensemble.  Each profile normalizes the algebraic NKCC1
arm to its own effective protocol calcium endpoint.  Radau and BDF are run for
every case.  One Tier-2 gland-flow scale per root/profile is fitted from the
preferred R1 member and then applied unchanged to all 20 regulatory members.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .calibration import build_wt_model
from .model import ModernFullModel, WT
from .nkcc_stimulation import (
    N1AlgebraicNkcc1,
    attach_stimulated_nkcc1,
    isolated_nkcc1_fractional_cl_uptake_s,
)
from .run_g5_cross_root_profile import FrozenRootCase, load_retained_wt_roots
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    WTFlowEnvelope,
    fit_wt_flow_scale,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    sha256_file,
    simulate_wt,
    wt_flow_envelope_feasibility,
)


PREFERRED_REGULATORY_MEMBER = "R1_G125_P21_PROSE_TREFERENCE"
PROFILE_FILENAME = "wt_g5_joint_ca_nkcc_regulatory_profile.csv"
GROUP_FILENAME = "wt_g5_joint_ca_nkcc_group_gate.csv"
SUMMARY_FILENAME = "wt_g5_joint_ca_nkcc_summary.json"
HASH_FILENAME = "wt_g5_joint_ca_nkcc_hashes.csv"


@dataclass(frozen=True)
class JointRepairProfile:
    profile_id: str
    effective_calcium_uM: float
    nkcc1_fully_activated_multiplier: float
    calcium_evidence_status: str
    nkcc_gain_evidence_status: str
    source_defensibility_status: str


@dataclass(frozen=True)
class JointSimulationRequest:
    root: FrozenRootCase
    profile: JointRepairProfile
    regulatory_member_id: str
    solver_label: str


def declared_joint_repair_profiles() -> tuple[JointRepairProfile, ...]:
    """Profiles requested from independent WT-only flow-shape screening."""

    return (
        JointRepairProfile(
            "CA010_N1_M1P75",
            0.10,
            1.75,
            "WT_FLOW_SHAPE_CALIBRATION_PROFILE; ACUTE_CA_NOT_MEASURED_IN_P15",
            "P15_DIRECTION_AND_NOMINAL_GAIN_PROFILE",
            "SOURCE_CENTERED_WT_PROFILE_REQUIRES_CA_DEFENSIBILITY_REVIEW",
        ),
        JointRepairProfile(
            "CA015_N1_M1P75",
            0.15,
            1.75,
            "WT_FLOW_SHAPE_CALIBRATION_PROFILE; ACUTE_CA_NOT_MEASURED_IN_P15",
            "P15_DIRECTION_AND_NOMINAL_GAIN_PROFILE",
            "SOURCE_CENTERED_WT_PROFILE_REQUIRES_CA_DEFENSIBILITY_REVIEW",
        ),
        JointRepairProfile(
            "CA012_N1_M2P0",
            0.12,
            2.0,
            "WT_FLOW_SHAPE_CALIBRATION_PROFILE; ACUTE_CA_NOT_MEASURED_IN_P15",
            (
                "SOURCE_CENTERED_MECHANISTIC_NKCC_DIAGNOSTIC_WITH_SPQ_CONTEXT_"
                "CAVEAT_NOT_POINT_FIT"
            ),
            (
                "CANDIDATE_ONLY_IF_EFFECTIVE_CA_AND_ASSAY_CONTEXT_CONVERSION_ARE_"
                "DEFENSIBLE"
            ),
        ),
    )


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty artifact {path.name}")
    fields: list[str] = []
    for row in rows:
        for name in row:
            if name not in fields:
                fields.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _sha256_payload(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _member_registry() -> Mapping[str, Any]:
    return {member.member_id: member for member in pre_reveal_regulatory_ensemble()}


def _request_metadata(request: JointSimulationRequest, member: Any) -> dict[str, Any]:
    timing = json.dumps(dict(member.timing_parameters), sort_keys=True)
    return {
        "root_id": request.root.root_id,
        "selected_reference": request.root.selected_reference,
        "variant_id": request.root.variant.variant_id,
        "routing_na_fraction": request.root.variant.mixed_bath_na_attempt_fraction,
        "apical_pump_fraction": request.root.variant.apical_pump_fraction,
        "apical_k_fraction": request.root.variant.apical_k_fraction,
        "profile_id": request.profile.profile_id,
        "effective_calcium_uM": request.profile.effective_calcium_uM,
        "nkcc1_fully_activated_multiplier": (
            request.profile.nkcc1_fully_activated_multiplier
        ),
        "nkcc1_normalization_rest_calcium_uM": 0.058,
        "nkcc1_normalization_stimulated_calcium_uM": (
            request.profile.effective_calcium_uM
        ),
        "calcium_evidence_status": request.profile.calcium_evidence_status,
        "nkcc_gain_evidence_status": request.profile.nkcc_gain_evidence_status,
        "source_defensibility_status": request.profile.source_defensibility_status,
        "regulatory_member_id": member.member_id,
        "regulatory_family": member.family,
        "ae4_fully_activated_multiplier": member.fully_activated_multiplier,
        "regulatory_timing_label": member.timing_label,
        "regulatory_timing_parameters_json": timing,
        "regulatory_kinetic_evidence_status": member.kinetic_evidence_status,
        "solver_label": request.solver_label,
    }


def _build_joint_model(request: JointSimulationRequest, member: Any) -> Any:
    resting = build_wt_model(
        request.root.variant,
        request.root.calibration,
        regulatory_model=member.regulatory_model,
    )
    protocol = SecretagogueProtocol(
        stimulated_calcium_uM=request.profile.effective_calcium_uM
    )
    template = ModernFullModel(
        parameters=resting.parameters,
        stimulus=protocol,
        regulatory_model=member.regulatory_model,
        ae4_parameters=resting.ae4_parameters,
        ae4_evaluator=resting.ae4_evaluator,
    )
    nkcc = N1AlgebraicNkcc1(
        fully_activated_multiplier=(
            request.profile.nkcc1_fully_activated_multiplier
        ),
        resting_calcium_uM=protocol.resting_calcium_uM,
        stimulated_calcium_uM=protocol.stimulated_calcium_uM,
    )
    return attach_stimulated_nkcc1(template, nkcc)


def _failed_row(
    request: JointSimulationRequest, member: Any, message: str
) -> dict[str, Any]:
    return {
        **_request_metadata(request, member),
        "success": False,
        "message": message,
        "positive_core": False,
        "all_flow_nonnegative": False,
        "numerical_gate_pass": False,
        "individual_shape_gate_pass": False,
    }


def _simulate_request(request: JointSimulationRequest) -> dict[str, Any]:
    registry = _member_registry()
    member = registry[request.regulatory_member_id]
    try:
        model = _build_joint_model(request, member)
        solver = (
            PRODUCTION_RADAU
            if request.solver_label == PRODUCTION_RADAU.label
            else PRODUCTION_BDF
        )
        trajectory = simulate_wt(
            model,
            request.root.core_state,
            family=f"{member.member_id}+N1_JOINT_REPAIR",
            initial_condition="baseline",
            solver=solver,
            time_s=physical_time_grid(duration_s=600.0, sample_step_s=5.0),
        )
        if not trajectory.success:
            return _failed_row(request, member, trajectory.message)
        minute_times = np.arange(60.0, 601.0, 60.0)
        minute_flow = np.interp(minute_times, trajectory.time_s, trajectory.flow_pL_s)
        shape = wt_flow_envelope_feasibility(trajectory, WTFlowEnvelope())
        scale = fit_wt_flow_scale(trajectory, WTFlowEnvelope())
        conservation_ratios = {
            name: trajectory.max_abs_conservation_residuals[name] / tolerance
            for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
        }
        maximum_conservation_ratio = max(conservation_ratios.values())
        numerical_gate = bool(
            trajectory.positive_core
            and np.all(trajectory.flow_pL_s >= 0.0)
            and maximum_conservation_ratio <= 1.0
        )
        onset_index = 1
        endpoint_index = len(trajectory.time_s) - 1
        onset = model.evaluate(
            float(trajectory.time_s[onset_index]),
            trajectory.states[:, onset_index],
            genotype=WT,
        )
        endpoint = model.evaluate(
            float(trajectory.time_s[endpoint_index]),
            trajectory.states[:, endpoint_index],
            genotype=WT,
        )
        onset_uptake = isolated_nkcc1_fractional_cl_uptake_s(
            model,
            float(trajectory.time_s[onset_index]),
            trajectory.states[:, onset_index],
            genotype=WT,
        )
        return {
            **_request_metadata(request, member),
            "success": trajectory.success,
            "message": trajectory.message,
            "positive_core": trajectory.positive_core,
            "all_flow_nonnegative": bool(np.all(trajectory.flow_pL_s >= 0.0)),
            "max_dimensionless_conservation_ratio": maximum_conservation_ratio,
            "numerical_gate_pass": numerical_gate,
            "minute_flow_60_pL_s": float(minute_flow[0]),
            "minute_flow_600_pL_s": float(minute_flow[-1]),
            "minute_flow_min_pL_s": float(np.min(minute_flow)),
            "minute_flow_max_pL_s": float(np.max(minute_flow)),
            "minute_flow_max_min_ratio": shape.model_max_min_ratio,
            "individual_shape_gate_pass": shape.compatible,
            "individual_scale_multiplier_uL_min_per_pL_s": (
                scale.multiplier_uL_min_per_pL_s
            ),
            "individual_scale_implied_effective_cell_count": (
                scale.effective_cell_count
            ),
            "implied_cell_count_status": "WARNING_NO_PRIMARY_GEOMETRY_RANGE",
            "ae4_capacity_multiplier_endpoint": float(
                trajectory.capacity_multiplier[-1]
            ),
            "nkcc1_capacity_multiplier_onset": float(
                onset.diagnostics.regulatory["nkcc1_capacity_multiplier"]
            ),
            "nkcc1_capacity_multiplier_endpoint": float(
                endpoint.diagnostics.regulatory["nkcc1_capacity_multiplier"]
            ),
            "isolated_nkcc_fractional_cl_rate_onset_s_inv": onset_uptake,
            "uptake_observation_status": (
                "MECHANISTIC_DIAGNOSTIC_NOT_EQUATED_TO_SPQ_TOTAL_WITHOUT_MAP"
            ),
            "endpoint_na_mM": float(trajectory.cell_na_mM[-1]),
            "endpoint_k_mM": float(trajectory.cell_k_mM[-1]),
            "endpoint_cl_mM": float(trajectory.cell_cl_mM[-1]),
            "endpoint_ph": float(trajectory.cell_ph[-1]),
            "endpoint_volume_pL": float(trajectory.cell_volume_pL[-1]),
            "endpoint_flow_pL_s": float(trajectory.flow_pL_s[-1]),
            "cumulative_flow_pL": float(trajectory.cumulative_flow_pL[-1]),
        }
    except Exception as exc:  # preserve every factorial case in the ledger
        return _failed_row(request, member, f"{type(exc).__name__}: {exc}")


def _run_requests(
    requests: Sequence[JointSimulationRequest], workers: int
) -> list[dict[str, Any]]:
    if workers <= 1:
        return [_simulate_request(request) for request in requests]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_simulate_request, requests, chunksize=1))


def _augment_common_scale(
    rows: list[dict[str, Any]], roots: Sequence[FrozenRootCase]
) -> list[dict[str, Any]]:
    envelope = WTFlowEnvelope()
    by_group: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        by_group.setdefault((row["root_id"], row["profile_id"]), []).append(row)
    group_rows: list[dict[str, Any]] = []
    for root in roots:
        for profile in declared_joint_repair_profiles():
            group = by_group[(root.root_id, profile.profile_id)]
            radau = [row for row in group if row["solver_label"] == PRODUCTION_RADAU.label]
            bdf = [row for row in group if row["solver_label"] == PRODUCTION_BDF.label]
            if len(radau) != 20 or len(bdf) != 20:
                raise AssertionError("each root/profile group requires 20 members per solver")
            reference = next(
                row
                for row in radau
                if row["regulatory_member_id"] == PREFERRED_REGULATORY_MEMBER
            )
            preferred_mean_scale = float(
                reference.get("individual_scale_multiplier_uL_min_per_pL_s", math.nan)
            )
            successful = [row for row in group if bool(row.get("success"))]
            if successful:
                feasible_scale_lower = max(
                    envelope.lower_uL_min / float(row["minute_flow_min_pL_s"])
                    for row in successful
                )
                feasible_scale_upper = min(
                    envelope.upper_uL_min / float(row["minute_flow_max_pL_s"])
                    for row in successful
                )
                common_scale_feasible = feasible_scale_lower <= feasible_scale_upper
            else:
                feasible_scale_lower = math.nan
                feasible_scale_upper = math.nan
                common_scale_feasible = False
            if common_scale_feasible:
                common_scale = (
                    preferred_mean_scale
                    if feasible_scale_lower
                    <= preferred_mean_scale
                    <= feasible_scale_upper
                    else 0.5 * (feasible_scale_lower + feasible_scale_upper)
                )
            else:
                common_scale = preferred_mean_scale
            for row in group:
                if not bool(row.get("success")):
                    row["common_scale_multiplier_uL_min_per_pL_s"] = common_scale
                    row["common_scaled_minute_min_uL_min"] = math.nan
                    row["common_scaled_minute_max_uL_min"] = math.nan
                    row["common_scale_member_gate_pass"] = False
                    continue
                scaled_min = common_scale * float(row["minute_flow_min_pL_s"])
                scaled_max = common_scale * float(row["minute_flow_max_pL_s"])
                row["common_scale_multiplier_uL_min_per_pL_s"] = common_scale
                row["common_scaled_minute_min_uL_min"] = scaled_min
                row["common_scaled_minute_max_uL_min"] = scaled_max
                row["common_scale_member_gate_pass"] = bool(
                    scaled_min >= envelope.lower_uL_min
                    and scaled_max <= envelope.upper_uL_min
                )
            pair_differences = []
            bdf_by_member = {row["regulatory_member_id"]: row for row in bdf}
            for radau_row in radau:
                bdf_row = bdf_by_member[radau_row["regulatory_member_id"]]
                if bool(radau_row.get("success")) and bool(bdf_row.get("success")):
                    pair_differences.append(
                        abs(
                            float(radau_row["minute_flow_max_min_ratio"])
                            - float(bdf_row["minute_flow_max_min_ratio"])
                        )
                    )
            group_pass = bool(
                all(bool(row.get("numerical_gate_pass")) for row in group)
                and all(bool(row.get("individual_shape_gate_pass")) for row in group)
                and all(bool(row.get("common_scale_member_gate_pass")) for row in group)
            )
            def finite_range(field: str) -> str:
                values = [
                    float(row[field])
                    for row in successful
                    if field in row and math.isfinite(float(row[field]))
                ]
                return json.dumps([min(values), max(values)]) if values else "[]"

            ratio_values = [
                float(row["minute_flow_max_min_ratio"])
                for row in successful
                if "minute_flow_max_min_ratio" in row
            ]
            group_rows.append(
                {
                    "root_id": root.root_id,
                    "selected_reference": root.selected_reference,
                    "profile_id": profile.profile_id,
                    "effective_calcium_uM": profile.effective_calcium_uM,
                    "nkcc1_fully_activated_multiplier": (
                        profile.nkcc1_fully_activated_multiplier
                    ),
                    "source_defensibility_status": profile.source_defensibility_status,
                    "regulatory_member_count": 20,
                    "solver_count": 2,
                    "all_integrations_success": all(
                        bool(row.get("success")) for row in group
                    ),
                    "all_numerical_gates_pass": all(
                        bool(row.get("numerical_gate_pass")) for row in group
                    ),
                    "all_individual_shape_gates_pass": all(
                        bool(row.get("individual_shape_gate_pass")) for row in group
                    ),
                    "all_common_scale_member_gates_pass": all(
                        bool(row.get("common_scale_member_gate_pass")) for row in group
                    ),
                    "common_scale_feasible_intersection": common_scale_feasible,
                    "common_scale_feasible_lower_uL_min_per_pL_s": (
                        feasible_scale_lower
                    ),
                    "common_scale_feasible_upper_uL_min_per_pL_s": (
                        feasible_scale_upper
                    ),
                    "preferred_member_mean_scale_uL_min_per_pL_s": (
                        preferred_mean_scale
                    ),
                    "common_scale_selection_rule": (
                        "USE_PREFERRED_MEMBER_MEAN_SCALE_IF_INSIDE_ENSEMBLE_"
                        "INTERSECTION_ELSE_USE_INTERSECTION_MIDPOINT"
                    ),
                    "common_scale_multiplier_uL_min_per_pL_s": common_scale,
                    "common_scale_implied_effective_cell_count": common_scale / 6.0e-5,
                    "implied_cell_count_status": "WARNING_NO_PRIMARY_GEOMETRY_RANGE",
                    "ensemble_min_ratio": min(ratio_values) if ratio_values else math.nan,
                    "ensemble_max_ratio": max(ratio_values) if ratio_values else math.nan,
                    "ensemble_endpoint_cl_range_mM": finite_range("endpoint_cl_mM"),
                    "ensemble_endpoint_ph_range": finite_range("endpoint_ph"),
                    "ensemble_endpoint_volume_range_pL": finite_range(
                        "endpoint_volume_pL"
                    ),
                    "max_abs_radau_bdf_ratio_difference": (
                        max(pair_differences) if pair_differences else math.nan
                    ),
                    "g5_joint_group_gate_pass": group_pass,
                }
            )
    return group_rows


def run_joint_repair_profile(
    results_directory: Path, *, workers: int = 1
) -> Mapping[str, Any]:
    roots = load_retained_wt_roots(results_directory)
    profiles = declared_joint_repair_profiles()
    members = pre_reveal_regulatory_ensemble()
    requests = tuple(
        JointSimulationRequest(root, profile, member.member_id, solver.label)
        for root in roots
        for profile in profiles
        for member in members
        for solver in (PRODUCTION_RADAU, PRODUCTION_BDF)
    )
    rows = _run_requests(requests, workers)
    group_rows = _augment_common_scale(rows, roots)

    profile_path = results_directory / PROFILE_FILENAME
    group_path = results_directory / GROUP_FILENAME
    _write_rows(profile_path, rows)
    _write_rows(group_path, group_rows)
    passing = [row for row in group_rows if bool(row["g5_joint_group_gate_pass"])]
    selected_passing = [
        row for row in passing if bool(row["selected_reference"])
    ]
    input_artifacts = (
        "wt_root_table.csv",
        "wt_root_states.csv",
        "wt_parameter_candidates.csv",
        "wt_selected_candidate.json",
        "dynamic_frozen_manifest.json",
    )
    summary = {
        "analysis": "WT_ONLY_G5_JOINT_EFFECTIVE_CA_NKCC1_REGULATORY_ENSEMBLE",
        "classification": (
            "AT_LEAST_ONE_JOINT_PROFILE_PASSES_FULL_WT_ENSEMBLE"
            if passing
            else "NO_JOINT_PROFILE_PASSES_FULL_WT_ENSEMBLE"
        ),
        "retained_root_count": len(roots),
        "joint_profile_count": len(profiles),
        "regulatory_member_count": len(members),
        "solver_count": 2,
        "total_simulation_count": len(rows),
        "selected_reference_root": next(
            root.root_id for root in roots if root.selected_reference
        ),
        "passing_group_count": len(passing),
        "passing_groups": passing,
        "selected_reference_passing_groups": selected_passing,
        "candidate_selection_status": (
            "EVIDENCE_DEFENSIBILITY_REVIEW_REQUIRED_BEFORE_FINAL_G5_FREEZE"
        ),
        "common_scale_rule": (
            "INTERSECT_9_TO_10_UL_MIN_MULTIPLIER_INTERVALS_ACROSS_ALL_20_"
            "MEMBERS_AND_BOTH_SOLVERS; USE_SELECTED_R1_MEAN_SCALE_IF_INSIDE_"
            "THE_INTERSECTION_OTHERWISE_ITS_MIDPOINT; APPLY_ONCE_UNCHANGED"
        ),
        "flow_scale_tier": "TIER_2_WT_OBSERVATION_SCALE",
        "implied_cell_count_policy": (
            "WARNING_ONLY_UNLESS_PRIMARY_GEOMETRY_RANGE_IS_SUPPLIED"
        ),
        "regulatory_kinetics_identified": False,
        "regulatory_kinetic_policy": (
            "FULL_FROZEN_ENSEMBLE_CARRIED_FORWARD; NONIDENTIFICATION_NOT_A_"
            "SEPARATE_BLOCKER"
        ),
        "nkcc_law": "N1_CCH_ALGEBRAIC_CAPACITY_ONLY",
        "nkcc_kinetic_parameter_added": False,
        "profile_definitions": [asdict(profile) for profile in profiles],
        "profile_definitions_sha256": _sha256_payload(
            [asdict(profile) for profile in profiles]
        ),
        "regulatory_ensemble_sha256": _sha256_payload(
            [member.metadata() for member in members]
        ),
        "input_artifact_sha256": {
            name: sha256_file(results_directory / name) for name in input_artifacts
        },
        "profile_sha256": sha256_file(profile_path),
        "group_gate_sha256": sha256_file(group_path),
        "firewall": "WT_ONLY_NO_PHENOTYPE_TARGET_ACCESS",
    }
    summary_path = results_directory / SUMMARY_FILENAME
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    hash_rows = [
        {"artifact": path.name, "sha256": sha256_file(path), "scope": "WT_ONLY"}
        for path in (profile_path, group_path, summary_path)
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
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    result = run_joint_repair_profile(args.results_directory, workers=args.workers)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = (
    "JointRepairProfile",
    "declared_joint_repair_profiles",
    "run_joint_repair_profile",
)
