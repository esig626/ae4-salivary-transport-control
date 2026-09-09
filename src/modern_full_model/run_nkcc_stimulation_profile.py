"""Profile the WT-only NKCC1 sustainment repair on the frozen selected root.

This script reads only ``wt_selected_candidate.json``.  It does not import or
open the strict held-out target ledger and it never evaluates a knockout.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .calibration import ModelVariant, WTCalibrationParameters, build_wt_model
from .camp_pka import R1EffectiveActivation
from .model import ModernFullModel
from .nkcc_stimulation import (
    N0BaselineNkcc1,
    N1AlgebraicNkcc1,
    N2EffectiveDynamicNkcc1,
    attach_stimulated_nkcc1,
    isolated_nkcc1_fractional_cl_uptake_s,
)
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    WTFlowEnvelope,
    physical_time_grid,
    simulate_wt,
    solver_comparison,
    wt_flow_envelope_feasibility,
)


# P15 Table 1 WT-only native-acinar measurements.  The uptake values are an
# SPQ net observation under a depletion/repletion protocol, not an NKCC1
# turnover.  They center a sensitivity scale but do not supply a hard capacity
# bound without the missing observation map.
P15_CCH_UPTAKE_S = 2.18e-3
P15_CCH_UPTAKE_SEM_S = 0.20e-3
P15_COMBINED_UPTAKE_S = 2.02e-3
P15_COMBINED_UPTAKE_SEM_S = 0.10e-3
SOURCE_CENTERED_MAX_MULTIPLIER = 2.0


def _write_rows(path: Path, rows: list[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError("profile table cannot be empty")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _load_frozen_template(selected_path: Path) -> tuple[ModernFullModel, np.ndarray, str]:
    selected = json.loads(selected_path.read_text(encoding="utf-8"))
    construction = selected["construction"]
    variant = ModelVariant(**construction["variant"])
    calibration = WTCalibrationParameters(**construction["calibration"])
    resting = build_wt_model(
        variant,
        calibration,
        regulatory_model=R1EffectiveActivation(),
    )
    protocol = SecretagogueProtocol()
    template = ModernFullModel(
        parameters=resting.parameters,
        stimulus=protocol,
        regulatory_model=resting.regulatory_model,
        ae4_parameters=resting.ae4_parameters,
        ae4_evaluator=resting.ae4_evaluator,
    )
    state = np.asarray(construction["root_state"], dtype=float)
    if state.shape != (12,) or np.any(~np.isfinite(state)) or np.any(state <= 0.0):
        raise ValueError("selected candidate does not contain a positive 12-state core root")
    return template, state, str(construction["root_id"])


def _law_panel() -> list[tuple[str, Any, str]]:
    panel: list[tuple[str, Any, str]] = [
        ("N0_M1", N0BaselineNkcc1(), "EXACT_BASELINE_CONTROL")
    ]
    multipliers = (1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 4.00, 6.00, 8.00)
    for multiplier in multipliers:
        scale_status = (
            "P15_CENTERED_CONTEXT_MISMATCH_SENSITIVITY"
            if multiplier <= SOURCE_CENTERED_MAX_MULTIPLIER
            else "UNBOUNDED_FLOW_SHAPE_PROBE"
        )
        panel.append(
            (
                f"N1_M{multiplier:g}",
                N1AlgebraicNkcc1(fully_activated_multiplier=multiplier),
                scale_status,
            )
        )
    for tau_s in (30.0, 60.0, 120.0, 240.0, 480.0):
        for multiplier in (1.75, 2.00, 3.00, 4.00, 6.00, 8.00):
            scale_status = (
                "P15_CENTERED_CONTEXT_MISMATCH_SENSITIVITY"
                if multiplier <= SOURCE_CENTERED_MAX_MULTIPLIER
                else "UNBOUNDED_FLOW_SHAPE_PROBE"
            )
            panel.append(
                (
                    f"N2_M{multiplier:g}_T{tau_s:g}",
                    N2EffectiveDynamicNkcc1(
                        fully_activated_multiplier=multiplier,
                        tau_activation_s=tau_s,
                    ),
                    scale_status,
                )
            )
    return panel


def _profile_row(
    *,
    law_id: str,
    law: Any,
    scale_status: str,
    model: Any,
    core_root: np.ndarray,
    trajectory: Any,
    root_id: str,
) -> Mapping[str, Any]:
    flow = wt_flow_envelope_feasibility(trajectory, WTFlowEnvelope())
    full_initial = np.r_[core_root, model.initial_state()[12:]]
    basal_isolated = isolated_nkcc1_fractional_cl_uptake_s(
        model, 0.0, full_initial
    )
    index_60 = int(np.argmin(np.abs(trajectory.time_s - 60.0)))
    index_600 = int(np.argmin(np.abs(trajectory.time_s - 600.0)))
    isolated_60 = isolated_nkcc1_fractional_cl_uptake_s(
        model,
        float(trajectory.time_s[index_60]),
        trajectory.states[:, index_60],
    )
    isolated_600 = isolated_nkcc1_fractional_cl_uptake_s(
        model,
        float(trajectory.time_s[index_600]),
        trajectory.states[:, index_600],
    )
    regulatory_60 = model.evaluate(
        float(trajectory.time_s[index_60]), trajectory.states[:, index_60]
    ).diagnostics.regulatory
    regulatory_600 = model.evaluate(
        float(trajectory.time_s[index_600]), trajectory.states[:, index_600]
    ).diagnostics.regulatory
    ion_gate = bool(
        np.min(trajectory.cell_na_mM) >= 4.0
        and np.max(trajectory.cell_na_mM) <= 60.0
        and np.min(trajectory.cell_k_mM) >= 60.0
        and np.max(trajectory.cell_k_mM) <= 210.0
        and np.min(trajectory.cell_cl_mM) >= 5.0
        and np.max(trajectory.cell_cl_mM) <= 90.0
        and np.min(trajectory.cell_ph) >= 6.2
        and np.max(trajectory.cell_ph) <= 8.0
        and np.min(trajectory.cell_volume_pL) >= 0.5
        and np.max(trajectory.cell_volume_pL) <= 3.0
    )
    conservation_gate = all(
        trajectory.max_abs_conservation_residuals[name] <= tolerance
        for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
    )
    source_centered = bool(
        getattr(law, "fully_activated_multiplier", 1.0)
        <= SOURCE_CENTERED_MAX_MULTIPLIER
    )
    return {
        "root_id": root_id,
        "law_id": law_id,
        "family": law.family,
        "fully_activated_multiplier": getattr(law, "fully_activated_multiplier", 1.0),
        "tau_activation_s": getattr(law, "tau_activation_s", "NA"),
        "gain_provenance_status": scale_status,
        "tau_provenance_status": (
            "UNMEASURED_WT_ONLY_SENSITIVITY"
            if hasattr(law, "tau_activation_s")
            else "NOT_APPLICABLE"
        ),
        "solver_success": trajectory.success,
        "positive_core": trajectory.positive_core,
        "conservation_gate_pass": conservation_gate,
        "broad_dynamic_ion_gate_pass": ion_gate,
        "WT_flow_shape_gate_pass": flow.compatible,
        "WT_flow_max_min_ratio": flow.model_max_min_ratio,
        "WT_flow_allowed_max_min_ratio": flow.envelope_max_min_ratio,
        "nkcc1_multiplier_60_s": regulatory_60["nkcc1_capacity_multiplier"],
        "nkcc1_multiplier_600_s": regulatory_600["nkcc1_capacity_multiplier"],
        "isolated_NKCC_fractional_Cl_rate_rest_s_inv": basal_isolated,
        "isolated_NKCC_fractional_Cl_rate_60_s_s_inv": isolated_60,
        "isolated_NKCC_fractional_Cl_rate_600_s_s_inv": isolated_600,
        "P15_CCh_total_SPQ_uptake_s_inv": P15_CCH_UPTAKE_S,
        "P15_CCh_total_SPQ_SEM_s_inv": P15_CCH_UPTAKE_SEM_S,
        "uptake_observation_map_status": "MISSING_DEPLETION_STATE_AND_SPQ_MAP_NOT_POINT_FIT",
        "cell_na_endpoint_mM": trajectory.cell_na_mM[-1],
        "cell_k_endpoint_mM": trajectory.cell_k_mM[-1],
        "cell_cl_endpoint_mM": trajectory.cell_cl_mM[-1],
        "cell_ph_endpoint": trajectory.cell_ph[-1],
        "cell_volume_endpoint_pL": trajectory.cell_volume_pL[-1],
        "flow_60_s_pL_s": trajectory.flow_pL_s[index_60],
        "flow_600_s_pL_s": trajectory.flow_pL_s[index_600],
        "cumulative_flow_600_s_pL": trajectory.cumulative_flow_pL[-1],
        "source_centered_scale": source_centered,
        "source_disciplined_acceptance": bool(
            trajectory.success
            and trajectory.positive_core
            and conservation_gate
            and ion_gate
            and flow.compatible
            and source_centered
        ),
    }


def run_profile(output_directory: Path, selected_path: Path) -> Mapping[str, Any]:
    output_directory.mkdir(parents=True, exist_ok=True)
    template, core_root, root_id = _load_frozen_template(selected_path)
    grid = physical_time_grid(duration_s=600.0, sample_step_s=5.0)
    rows: list[Mapping[str, Any]] = []
    trajectories: dict[str, tuple[Any, Any]] = {}
    for law_id, law, scale_status in _law_panel():
        model = attach_stimulated_nkcc1(template, law)
        trajectory = simulate_wt(
            model,
            core_root,
            family=f"R1+{law.family}",
            initial_condition="baseline",
            solver=PRODUCTION_RADAU,
            time_s=grid,
        )
        row = _profile_row(
            law_id=law_id,
            law=law,
            scale_status=scale_status,
            model=model,
            core_root=core_root,
            trajectory=trajectory,
            root_id=root_id,
        )
        rows.append(row)
        trajectories[law_id] = (model, trajectory)
        print(
            f"{law_id} ratio={row['WT_flow_max_min_ratio']:.6g} "
            f"shape={row['WT_flow_shape_gate_pass']} accepted={row['source_disciplined_acceptance']}",
            flush=True,
        )

    _write_rows(output_directory / "wt_nkcc_stimulation_profile.csv", rows)
    source_rows = [row for row in rows if row["source_centered_scale"]]
    passing_shape = [row for row in rows if row["WT_flow_shape_gate_pass"]]
    source_accepted = [row for row in rows if row["source_disciplined_acceptance"]]
    best_source = min(source_rows, key=lambda row: row["WT_flow_max_min_ratio"])
    first_shape = (
        min(
            passing_shape,
            key=lambda row: (
                row["fully_activated_multiplier"],
                float(row["tau_activation_s"])
                if row["tau_activation_s"] != "NA"
                else 0.0,
                row["WT_flow_max_min_ratio"],
            ),
        )
        if passing_shape
        else None
    )

    crosscheck_ids = [best_source["law_id"]]
    if first_shape is not None and first_shape["law_id"] not in crosscheck_ids:
        crosscheck_ids.append(first_shape["law_id"])
    solver_payload: dict[str, Any] = {}
    for law_id in crosscheck_ids:
        model, radau = trajectories[law_id]
        bdf = simulate_wt(
            model,
            core_root,
            family=radau.family,
            initial_condition="baseline",
            solver=PRODUCTION_BDF,
            time_s=grid,
        )
        solver_payload[law_id] = solver_comparison(radau, bdf)
    (output_directory / "wt_nkcc_solver_crosscheck.json").write_text(
        json.dumps(solver_payload, indent=2, sort_keys=True), encoding="utf-8"
    )

    summary = {
        "generation_id": "G5_NKCC1_STIMULATION_REPAIR_PROFILE_SELECTED_ROOT",
        "parent_root_id": root_id,
        "heldout_access": "NONE",
        "genotypes_evaluated": ["WT"],
        "evidence": {
            "P15_CCh_total_SPQ_uptake_s_inv": P15_CCH_UPTAKE_S,
            "P15_CCh_total_SPQ_SEM_s_inv": P15_CCH_UPTAKE_SEM_S,
            "P15_combined_total_SPQ_uptake_s_inv": P15_COMBINED_UPTAKE_S,
            "P15_combined_total_SPQ_SEM_s_inv": P15_COMBINED_UPTAKE_SEM_S,
            "P15_bumetanide_isolated_reduction_percent": 95.6,
            "P15_bumetanide_isolated_reduction_SEM_percent": 0.9,
            "observation_map_limit": (
                "The SPQ target uses a low/high-Cl depleted-state protocol and is a net "
                "whole-cell slope; it is not fitted as an NKCC1 capacity constant."
            ),
        },
        "laws": {
            "N0": "existing reversible NKCC1 law, exact control",
            "N1": "algebraic CCh/Ca capacity multiplier; one gain",
            "N2": "one effective CCh/Ca activation state; one gain plus unmeasured tau",
        },
        "profile_count": len(rows),
        "source_centered_profile_count": len(source_rows),
        "flow_shape_passing_count": len(passing_shape),
        "source_disciplined_accepted_count": len(source_accepted),
        "best_source_centered_row": dict(best_source),
        "lowest_gain_flow_shape_passing_row": (
            dict(first_shape) if first_shape is not None else None
        ),
        "classification": (
            "SOURCE_SUPPORTED_NKCC1_ONLY_REPAIR_FOUND"
            if source_accepted
            else "NKCC1_ONLY_REPAIR_NOT_SOURCE_CONSTRAINED"
        ),
        "interpretation": (
            "The nested law can be a numerical sustainment repair only if its WT-only "
            "gain/time profile passes. A passing unbounded probe is not accepted as a "
            "physiological calibration because P15 does not identify the native "
            "stimulated NKCC1 capacity or activation time course."
        ),
        "decision_critical_measurement": (
            "Native SMG-acinar bumetanide-sensitive NKCC1 Na/K/Cl flux versus physical "
            "time during the matched 0.3 uM CCh + 5 uM IPR protocol, with simultaneous "
            "intracellular Na, K, Cl and volume."
        ),
        "solver_crosscheck_ids": crosscheck_ids,
    }
    (output_directory / "wt_nkcc_stimulation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    parser.add_argument(
        "--selected",
        type=Path,
        default=Path("results/13B_modern_full_model/wt_selected_candidate.json"),
    )
    args = parser.parse_args()
    print(json.dumps(run_profile(args.output, args.selected), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
