"""Independent P15 AE2-null validation for WT-frozen joint G5 profiles.

P15 reports no detected difference in Ae2-line secretion kinetics or 10-min
total, but provides no equivalence margin.  This runner therefore records the
exact WT-versus-AE2-null effects and applies only a conspicuous-effect warning
threshold (5%, a modeling diagnostic, not a source uncertainty interval).
It never treats numerical proximity as proof of biological equivalence.

Every retained WT root and every joint effective-Ca/NKCC1 profile is tested
with the preferred R1 AE4 regulatory member under both Radau and BDF.  WT and
AE2-null start from the same frozen WT core to isolate the dynamic deletion
effect; this is not presented as a matched Ae2-line resting calibration.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .calibration import build_wt_model
from .model import AE2_NULL, ModernFullModel, WT
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .run_g5_cross_root_profile import FrozenRootCase, load_retained_wt_roots
from .run_g5_joint_repair_profile import (
    PREFERRED_REGULATORY_MEMBER,
    JointRepairProfile,
    declared_joint_repair_profiles,
)
from .states import CORE_STATE_NAMES
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    attach_basal_regulation,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    sha256_file,
)


PROFILE_FILENAME = "wt_g5_ae2_independent_validation.csv"
SUMMARY_FILENAME = "wt_g5_ae2_independent_summary.json"
HASH_FILENAME = "wt_g5_ae2_independent_hashes.csv"
LARGE_EFFECT_WARNING_FRACTION = 0.05


@dataclass(frozen=True)
class AE2ValidationRequest:
    root: FrozenRootCase
    profile: JointRepairProfile
    solver_label: str


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


def _build_model(request: AE2ValidationRequest) -> Any:
    member = next(
        item
        for item in pre_reveal_regulatory_ensemble()
        if item.member_id == PREFERRED_REGULATORY_MEMBER
    )
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
    return attach_stimulated_nkcc1(
        template,
        N1AlgebraicNkcc1(
            fully_activated_multiplier=(
                request.profile.nkcc1_fully_activated_multiplier
            ),
            resting_calcium_uM=protocol.resting_calcium_uM,
            stimulated_calcium_uM=protocol.stimulated_calcium_uM,
        ),
    )


def _simulate_genotype(
    model: Any,
    core_state: Sequence[float],
    *,
    genotype: Any,
    solver_label: str,
) -> dict[str, Any]:
    solver = (
        PRODUCTION_RADAU if solver_label == PRODUCTION_RADAU.label else PRODUCTION_BDF
    )
    grid = physical_time_grid(duration_s=600.0, sample_step_s=5.0)
    y0 = attach_basal_regulation(model, core_state)
    integration_grid = grid[1:]
    result = model.solve_dynamics(
        (float(integration_grid[0]), float(integration_grid[-1])),
        initial_state=y0,
        genotype=genotype,
        method=solver.method,
        rtol=solver.rtol,
        atol=solver.atol_vector(model.state_names),
        t_eval=integration_grid,
        max_step_s=solver.max_step_s,
    )
    complete = result.y.shape[1] == integration_grid.size
    if not complete:
        return {
            "success": False,
            "message": result.message,
            "positive_core": False,
            "numerical_gate_pass": False,
        }
    states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    flow = np.empty(grid.size, dtype=float)
    cl = np.empty(grid.size, dtype=float)
    ph = np.empty(grid.size, dtype=float)
    volume = states[5].copy()
    conservation_max = {name: 0.0 for name in CONSERVATION_RESIDUAL_TOLERANCES}
    for index, (time_s, state) in enumerate(zip(grid, states.T)):
        evaluation = model.evaluate(float(time_s), state, genotype=genotype)
        diagnostics = evaluation.diagnostics
        flow[index] = diagnostics.water.lumen_outflow_pL_s
        cl[index] = diagnostics.observables.cell_concentrations_mM["cl"]
        ph[index] = diagnostics.observables.cell_acid_base.ph
        for name, value in diagnostics.conservation_residuals.items():
            conservation_max[name] = max(conservation_max[name], abs(float(value)))
    conservation_ratios = {
        name: conservation_max[name] / tolerance
        for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
    }
    minute_times = np.arange(60.0, 601.0, 60.0)
    minute_flow = np.interp(minute_times, grid, flow)
    cumulative = np.concatenate((np.asarray((0.0,)), cumulative_trapezoid(flow, grid)))
    numerical_gate = bool(
        result.success
        and np.all(states[: len(CORE_STATE_NAMES)] > 0.0)
        and np.all(flow >= 0.0)
        and max(conservation_ratios.values()) <= 1.0
    )
    return {
        "success": bool(result.success),
        "message": str(result.message),
        "positive_core": bool(np.all(states[: len(CORE_STATE_NAMES)] > 0.0)),
        "all_flow_nonnegative": bool(np.all(flow >= 0.0)),
        "max_dimensionless_conservation_ratio": max(conservation_ratios.values()),
        "numerical_gate_pass": numerical_gate,
        "minute_flow_pL_s": minute_flow,
        "cumulative_flow_pL": float(cumulative[-1]),
        "endpoint_flow_pL_s": float(flow[-1]),
        "endpoint_cl_mM": float(cl[-1]),
        "endpoint_ph": float(ph[-1]),
        "endpoint_volume_pL": float(volume[-1]),
    }


def _relative_difference(comparison: float, reference: float) -> float:
    return float((comparison - reference) / reference)


def _simulate_request(request: AE2ValidationRequest) -> dict[str, Any]:
    try:
        model = _build_model(request)
        wt = _simulate_genotype(
            model,
            request.root.core_state,
            genotype=WT,
            solver_label=request.solver_label,
        )
        ae2 = _simulate_genotype(
            model,
            request.root.core_state,
            genotype=AE2_NULL,
            solver_label=request.solver_label,
        )
        if not bool(wt.get("success")) or not bool(ae2.get("success")):
            return {
                "root_id": request.root.root_id,
                "selected_reference": request.root.selected_reference,
                "profile_id": request.profile.profile_id,
                "solver_label": request.solver_label,
                "wt_success": wt.get("success", False),
                "ae2_null_success": ae2.get("success", False),
                "paired_numerical_gate_pass": False,
                "message": f"WT: {wt.get('message')}; AE2: {ae2.get('message')}",
            }
        wt_minute = np.asarray(wt.pop("minute_flow_pL_s"), dtype=float)
        ae2_minute = np.asarray(ae2.pop("minute_flow_pL_s"), dtype=float)
        cumulative_relative = _relative_difference(
            float(ae2["cumulative_flow_pL"]), float(wt["cumulative_flow_pL"])
        )
        minute_relative = (ae2_minute - wt_minute) / wt_minute
        conspicuous = bool(
            abs(cumulative_relative) > LARGE_EFFECT_WARNING_FRACTION
            or np.max(np.abs(minute_relative)) > LARGE_EFFECT_WARNING_FRACTION
        )
        return {
            "root_id": request.root.root_id,
            "selected_reference": request.root.selected_reference,
            "variant_id": request.root.variant.variant_id,
            "profile_id": request.profile.profile_id,
            "effective_calcium_uM": request.profile.effective_calcium_uM,
            "nkcc1_fully_activated_multiplier": (
                request.profile.nkcc1_fully_activated_multiplier
            ),
            "regulatory_member_id": PREFERRED_REGULATORY_MEMBER,
            "solver_label": request.solver_label,
            "wt_success": wt["success"],
            "ae2_null_success": ae2["success"],
            "wt_numerical_gate_pass": wt["numerical_gate_pass"],
            "ae2_null_numerical_gate_pass": ae2["numerical_gate_pass"],
            "paired_numerical_gate_pass": bool(
                wt["numerical_gate_pass"] and ae2["numerical_gate_pass"]
            ),
            "wt_cumulative_flow_pL": wt["cumulative_flow_pL"],
            "ae2_null_cumulative_flow_pL": ae2["cumulative_flow_pL"],
            "ae2_minus_wt_cumulative_relative": cumulative_relative,
            "max_abs_minute_flow_relative_difference": float(
                np.max(np.abs(minute_relative))
            ),
            "wt_endpoint_flow_pL_s": wt["endpoint_flow_pL_s"],
            "ae2_null_endpoint_flow_pL_s": ae2["endpoint_flow_pL_s"],
            "ae2_minus_wt_endpoint_flow_relative": _relative_difference(
                float(ae2["endpoint_flow_pL_s"]), float(wt["endpoint_flow_pL_s"])
            ),
            "wt_endpoint_cl_mM": wt["endpoint_cl_mM"],
            "ae2_null_endpoint_cl_mM": ae2["endpoint_cl_mM"],
            "ae2_minus_wt_endpoint_cl_mM": (
                float(ae2["endpoint_cl_mM"]) - float(wt["endpoint_cl_mM"])
            ),
            "wt_endpoint_ph": wt["endpoint_ph"],
            "ae2_null_endpoint_ph": ae2["endpoint_ph"],
            "ae2_minus_wt_endpoint_ph": (
                float(ae2["endpoint_ph"]) - float(wt["endpoint_ph"])
            ),
            "wt_endpoint_volume_pL": wt["endpoint_volume_pL"],
            "ae2_null_endpoint_volume_pL": ae2["endpoint_volume_pL"],
            "ae2_minus_wt_endpoint_volume_pL": (
                float(ae2["endpoint_volume_pL"])
                - float(wt["endpoint_volume_pL"])
            ),
            "large_secretion_effect_warning": conspicuous,
            "large_effect_warning_threshold_fraction": (
                LARGE_EFFECT_WARNING_FRACTION
            ),
            "warning_threshold_status": (
                "MODELING_DIAGNOSTIC_NOT_SOURCE_EQUIVALENCE_MARGIN"
            ),
            "p15_qualitative_consistency_status": (
                "NO_CONSPICUOUS_MODELED_EFFECT"
                if not conspicuous
                else "CONSPICUOUS_MODELED_EFFECT_REQUIRES_REJECTION_REVIEW"
            ),
            "initial_state_status": (
                "SAME_FROZEN_WT_ROOT_DYNAMIC_DELETION_NOT_AE2_LINE_REST_FIT"
            ),
        }
    except Exception as exc:
        return {
            "root_id": request.root.root_id,
            "selected_reference": request.root.selected_reference,
            "profile_id": request.profile.profile_id,
            "solver_label": request.solver_label,
            "wt_success": False,
            "ae2_null_success": False,
            "paired_numerical_gate_pass": False,
            "message": f"{type(exc).__name__}: {exc}",
        }


def _run_requests(
    requests: Sequence[AE2ValidationRequest], workers: int
) -> list[dict[str, Any]]:
    if workers <= 1:
        return [_simulate_request(request) for request in requests]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_simulate_request, requests, chunksize=1))


def run_ae2_validation(results_directory: Path, *, workers: int = 1) -> Mapping[str, Any]:
    roots = load_retained_wt_roots(results_directory)
    profiles = declared_joint_repair_profiles()
    requests = tuple(
        AE2ValidationRequest(root, profile, solver.label)
        for root in roots
        for profile in profiles
        for solver in (PRODUCTION_RADAU, PRODUCTION_BDF)
    )
    rows = _run_requests(requests, workers)
    profile_path = results_directory / PROFILE_FILENAME
    _write_rows(profile_path, rows)

    selected_rows = [row for row in rows if bool(row["selected_reference"])]
    selected_m2 = [
        row for row in selected_rows if row["profile_id"] == "CA012_N1_M2P0"
    ]
    summary = {
        "analysis": "P15_AE2_NULL_INDEPENDENT_VALIDATION_OF_JOINT_G5_PROFILES",
        "classification": (
            "ALL_PROFILES_NO_CONSPICUOUS_AE2_SECRETION_EFFECT"
            if all(not bool(row.get("large_secretion_effect_warning", True)) for row in rows)
            else "AT_LEAST_ONE_PROFILE_HAS_CONSPICUOUS_AE2_SECRETION_EFFECT"
        ),
        "retained_root_count": len(roots),
        "profile_count": len(profiles),
        "solver_count": 2,
        "simulation_pair_count": len(rows),
        "all_paired_numerical_gates_pass": all(
            bool(row["paired_numerical_gate_pass"]) for row in rows
        ),
        "large_effect_warning_threshold_fraction": LARGE_EFFECT_WARNING_FRACTION,
        "warning_threshold_status": (
            "MODELING_DIAGNOSTIC_NOT_SOURCE_EQUIVALENCE_MARGIN"
        ),
        "p15_target": (
            "NO_DETECTED_AE2_LINE_DIFFERENCE_IN_KINETICS_OR_10_MIN_TOTAL; "
            "NO_EQUIVALENCE_MARGIN_REPORTED"
        ),
        "selected_reference_ca012_m2_rows": selected_m2,
        "initial_state_scope": (
            "DYNAMIC_DELETION_FROM_EACH_FROZEN_WT_ROOT; NOT_A_MATCHED_AE2_LINE_"
            "RESTING_CALIBRATION"
        ),
        "profile_definitions_sha256": _sha256_payload(
            [profile.__dict__ for profile in profiles]
        ),
        "input_joint_group_gate_sha256": sha256_file(
            results_directory / "wt_g5_joint_ca_nkcc_group_gate.csv"
        ),
        "profile_sha256": sha256_file(profile_path),
        "firewall": "NO_AE4_NULL_TARGET_OR_TRAJECTORY_USED",
    }
    summary_path = results_directory / SUMMARY_FILENAME
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    hash_rows = [
        {"artifact": path.name, "sha256": sha256_file(path), "scope": "P15_AE2_ONLY"}
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
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")
    result = run_ae2_validation(args.results_directory, workers=args.workers)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ("run_ae2_validation",)
