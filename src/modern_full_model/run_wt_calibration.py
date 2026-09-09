"""Generate Task-13B WT-only calibration artifacts.

The strict AE4-null secretion ledger is neither imported nor opened here.
This runner records the pre-audit water-scale root only as rejected, then
profiles the corrected Palk water scale across mixed-bath routing, apical
topology, and conductance-scale axes.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Any, Mapping

from .camp_pka import R1EffectiveActivation
from .calibration import (
    CALIBRATED_PARAMETER_NAMES,
    ModelVariant,
    RootSearchReport,
    WTCalibrationParameters,
    WTCalibrationSpec,
    build_wt_model,
    default_variant_panel,
    fit_wt_source_compatible_variant,
    profile_wt_parameter_ablations,
    report_to_dict,
    solve_wt_resting_branches,
)
from .parameters import FullModelParameters, parameter_records


STATE_UNITS = {
    "na_i_fmol": "fmol",
    "k_i_fmol": "fmol",
    "cl_i_fmol": "fmol",
    "tic_i_fmol": "fmol",
    "alk_i_fmol": "fmol-equivalent",
    "volume_i_pL": "pL",
    "na_l_fmol": "fmol",
    "k_l_fmol": "fmol",
    "cl_l_fmol": "fmol",
    "tic_l_fmol": "fmol",
    "alk_l_fmol": "fmol-equivalent",
    "volume_l_pL": "pL",
}


def _variant_with_scale(variant: ModelVariant, scale: float) -> ModelVariant:
    label = str(int(scale)) if float(scale).is_integer() else str(scale).replace(".", "p")
    return replace(variant, variant_id=f"{variant.variant_id}_G{label}")


def _fit_profile_task(
    task: tuple[ModelVariant, float, int, int, str],
) -> RootSearchReport:
    """Process-safe unit of work for the topology/routing profile."""

    variant, scale, starts, seed, stage = task
    report = fit_wt_source_compatible_variant(
        variant,
        fixed_conductance_scale=scale,
        start_count=starts,
        seed=seed,
    )
    return replace(
        report,
        start_design={**report.start_design, "profile_stage": stage},
    )


def _run_profile_tasks(
    tasks: list[tuple[ModelVariant, float, int, int, str]],
    *,
    workers: int,
) -> list[RootSearchReport]:
    if workers <= 1:
        reports: list[RootSearchReport] = []
        for task in tasks:
            print(f"WT_{task[4].upper()} {task[0].variant_id}", flush=True)
            reports.append(_fit_profile_task(task))
        return reports

    reports_by_index: dict[int, RootSearchReport] = {}
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(_fit_profile_task, task): index
            for index, task in enumerate(tasks)
        }
        for future in as_completed(futures):
            index = futures[future]
            report = future.result()
            reports_by_index[index] = report
            print(f"WT_{tasks[index][4].upper()}_DONE {report.variant.variant_id}", flush=True)
    return [reports_by_index[index] for index in range(len(tasks))]


def _write_attempts(path: Path, reports: list[RootSearchReport]) -> None:
    fields = (
        "variant_id", "conductance_scale", "start_id", "optimizer_success",
        "admissible_state", "converged_root", "max_abs_scaled_residual",
        "max_abs_amount_rhs_fmol_s", "max_abs_volume_rhs_pL_s",
        "max_abs_regulatory_rhs_s_inv",
        "max_abs_omitted_charge_rhs_fmol_equivalent_s", "cost", "optimality",
        "nfev", "boundary_hits", "message", "profile_stage",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for report in reports:
            scale = report.start_design["fixed_conductance_scale"]
            for attempt in report.attempts:
                writer.writerow(
                    {
                        "variant_id": report.variant.variant_id,
                        "conductance_scale": scale,
                        "start_id": attempt.start_id,
                        "optimizer_success": attempt.optimizer_success,
                        "admissible_state": attempt.admissible_state,
                        "converged_root": attempt.converged_root,
                        "max_abs_scaled_residual": attempt.max_abs_scaled_independent_rhs,
                        "max_abs_amount_rhs_fmol_s": attempt.max_abs_amount_rhs_fmol_s,
                        "max_abs_volume_rhs_pL_s": attempt.max_abs_volume_rhs_pL_s,
                        "max_abs_regulatory_rhs_s_inv": attempt.max_abs_regulatory_rhs_s_inv,
                        "max_abs_omitted_charge_rhs_fmol_equivalent_s": (
                            attempt.max_abs_omitted_charge_rhs_fmol_equivalent_s
                        ),
                        "cost": attempt.cost,
                        "optimality": attempt.optimality,
                        "nfev": attempt.nfev,
                        "boundary_hits": ";".join(attempt.boundary_hits),
                        "message": attempt.message,
                        "profile_stage": report.start_design.get("profile_stage", "screen"),
                    }
                )


def _write_roots(path: Path, reports: list[RootSearchReport]) -> None:
    fields = (
        "variant_id", "root_id", "status", "routing_na_fraction",
        "apical_pump_fraction", "apical_k_fraction", "conductance_scale",
        "converged_starts", "distinct_root_count", "passes_numerical_gate",
        "passes_wt_gate", "boundary_hits", "source_supported_boundary_hits",
        "disallowed_boundary_hits", "gate_failures",
        "max_abs_scaled_independent_rhs", "max_abs_amount_rhs_fmol_s",
        "max_abs_volume_rhs_pL_s", "max_abs_regulatory_rhs_s_inv",
        "max_abs_omitted_charge_rhs_fmol_equivalent_s",
        "jacobian_rank", "jacobian_nullity",
        "jacobian_smallest_singular_value", "na_i_mM", "k_i_mM", "cl_i_mM",
        "tic_i_mM", "hco3_i_mM", "ph_i", "volume_i_pL", "v_apical_mV",
        "v_basolateral_mV", "na_l_mM", "k_l_mM", "cl_l_mM", "tic_l_mM",
        "ph_l", "volume_l_pL", "rest_outflow_pL_s",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for report in reports:
            common = {
                "variant_id": report.variant.variant_id,
                "routing_na_fraction": report.variant.mixed_bath_na_attempt_fraction,
                "apical_pump_fraction": report.variant.apical_pump_fraction,
                "apical_k_fraction": report.variant.apical_k_fraction,
                "conductance_scale": report.start_design["fixed_conductance_scale"],
                "converged_starts": sum(item.converged_root for item in report.attempts),
                "distinct_root_count": len(report.branches),
            }
            if not report.branches:
                best = min(report.attempts, key=lambda item: item.max_abs_scaled_independent_rhs)
                writer.writerow(
                    {
                        **common,
                        "root_id": "",
                        "status": "NO_EXACT_INTERIOR_ROOT",
                        "passes_numerical_gate": False,
                        "passes_wt_gate": False,
                        "boundary_hits": ";".join(best.boundary_hits),
                        "gate_failures": "no exact non-boundary charge-manifold root",
                        "max_abs_scaled_independent_rhs": best.max_abs_scaled_independent_rhs,
                        "max_abs_amount_rhs_fmol_s": best.max_abs_amount_rhs_fmol_s,
                        "max_abs_volume_rhs_pL_s": best.max_abs_volume_rhs_pL_s,
                        "max_abs_regulatory_rhs_s_inv": best.max_abs_regulatory_rhs_s_inv,
                        "max_abs_omitted_charge_rhs_fmol_equivalent_s": (
                            best.max_abs_omitted_charge_rhs_fmol_equivalent_s
                        ),
                    }
                )
                continue
            for branch in report.branches:
                writer.writerow(
                    {
                        **common,
                        "root_id": branch.branch_id,
                        "status": (
                            "WT_REST_PASS_R0_CONTROL_ONLY"
                            if branch.passes_wt_gate
                            else "ROOT_FOUND_WT_GATE_FAIL"
                        ),
                        "passes_numerical_gate": branch.passes_numerical_gate,
                        "passes_wt_gate": branch.passes_wt_gate,
                        "boundary_hits": ";".join(branch.boundary_hits),
                        "source_supported_boundary_hits": ";".join(
                            branch.source_supported_boundary_hits
                        ),
                        "disallowed_boundary_hits": ";".join(
                            branch.disallowed_boundary_hits
                        ),
                        "gate_failures": ";".join(branch.gate_failures),
                        "max_abs_scaled_independent_rhs": max(
                            abs(value) for value in branch.scaled_independent_rhs
                        ),
                        "max_abs_amount_rhs_fmol_s": max(
                            abs(branch.raw_rhs[index])
                            for index in (0, 1, 2, 3, 6, 7, 8, 9)
                        ),
                        "max_abs_volume_rhs_pL_s": max(
                            abs(branch.raw_rhs[index]) for index in (5, 11)
                        ),
                        "max_abs_regulatory_rhs_s_inv": (
                            max(abs(value) for value in branch.raw_rhs[12:])
                            if len(branch.raw_rhs) > 12 else 0.0
                        ),
                        "max_abs_omitted_charge_rhs_fmol_equivalent_s": max(
                            abs(value) for value in branch.omitted_raw_rhs
                        ),
                        "jacobian_rank": branch.independent_jacobian_rank,
                        "jacobian_nullity": branch.independent_jacobian_nullity,
                        "jacobian_smallest_singular_value": min(
                            branch.independent_jacobian_singular_values
                        ),
                        **branch.observables,
                    }
                )


def _write_states_and_parameters(
    state_path: Path,
    parameter_path: Path,
    reports: list[RootSearchReport],
) -> None:
    with state_path.open("w", newline="", encoding="utf-8") as state_handle, parameter_path.open(
        "w", newline="", encoding="utf-8"
    ) as parameter_handle:
        state_writer = csv.DictWriter(
            state_handle,
            fieldnames=("variant_id", "root_id", "state_name", "value", "unit"),
            lineterminator="\n",
        )
        parameter_writer = csv.DictWriter(
            parameter_handle,
            fieldnames=(
                "variant_id", "root_id", "parameter", "value", "unit", "tier",
                "provenance", "fit_status", "note",
            ),
            lineterminator="\n",
        )
        state_writer.writeheader()
        parameter_writer.writeheader()
        metadata = {
            CALIBRATED_PARAMETER_NAMES[0]: (
                "fmol carrier", 2,
                "WT-rest composite; attempt rates fixed at the declared 1 s^-1 gauge; amount/rate split unidentifiable"
            ),
            CALIBRATED_PARAMETER_NAMES[1]: (
                "fmol s^-1", 2, "WT-rest acid-base closure only"
            ),
            CALIBRATED_PARAMETER_NAMES[2]: (
                "dimensionless", 3, "profiled fixed voltage-gate axis; not identified"
            ),
            CALIBRATED_PARAMETER_NAMES[3]: (
                "fmol osmoles", 3, "OTHER pool; excludes explicit finite buffer"
            ),
        }
        for report in reports:
            for branch in report.branches:
                names = (
                    "na_i_fmol", "k_i_fmol", "cl_i_fmol", "tic_i_fmol",
                    "alk_i_fmol", "volume_i_pL", "na_l_fmol", "k_l_fmol",
                    "cl_l_fmol", "tic_l_fmol", "alk_l_fmol", "volume_l_pL",
                )
                for name, value in zip(names, branch.state):
                    state_writer.writerow(
                        {
                            "variant_id": report.variant.variant_id,
                            "root_id": branch.branch_id,
                            "state_name": name,
                            "value": value,
                            "unit": STATE_UNITS.get(name, "dimensionless"),
                        }
                    )
                for name, value in branch.parameters.items():
                    unit, tier, note = metadata[name]
                    parameter_writer.writerow(
                        {
                            "variant_id": report.variant.variant_id,
                            "root_id": branch.branch_id,
                            "parameter": name,
                            "value": value,
                            "unit": unit,
                            "tier": tier,
                            "provenance": "NEW_MODELING_DECISION",
                            "fit_status": (
                                "PROFILED_FIXED_AXIS"
                                if name == CALIBRATED_PARAMETER_NAMES[2]
                                else "FITTED_WT_ONLY"
                            ),
                            "note": note,
                        }
                    )
                for name, value, note in (
                    ("mixed_bath_na_attempt_fraction", report.variant.mixed_bath_na_attempt_fraction, "unmeasured native mixed-bath routing"),
                    ("apical_pump_fraction", report.variant.apical_pump_fraction, "localization supported; fraction unmeasured"),
                    ("apical_k_fraction", report.variant.apical_k_fraction, "apical current supported; fraction unmeasured"),
                ):
                    parameter_writer.writerow(
                        {
                            "variant_id": report.variant.variant_id,
                            "root_id": branch.branch_id,
                            "parameter": name,
                            "value": value,
                            "unit": "fraction",
                            "tier": 3,
                            "provenance": "NEW_MODELING_DECISION",
                            "fit_status": "PROFILED_FIXED_AXIS",
                            "note": note,
                        }
                    )


def _upgrade_root_table_boundary_schema(path: Path) -> None:
    """Add explicit allowed/disallowed-boundary columns to older R0 output."""

    if not path.exists():
        return
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or ())
        rows = list(reader)
    additions = ("source_supported_boundary_hits", "disallowed_boundary_hits")
    if all(name in fieldnames for name in additions):
        return
    insertion = fieldnames.index("boundary_hits") + 1
    for name in reversed(additions):
        if name not in fieldnames:
            fieldnames.insert(insertion, name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            row.setdefault("source_supported_boundary_hits", "")
            row.setdefault("disallowed_boundary_hits", row.get("boundary_hits", ""))
            writer.writerow(row)


def _upgrade_selected_geometry_provenance(output_directory: Path) -> None:
    geometry_path = output_directory / "wt_selected_root_geometry.json"
    selected_path = output_directory / "wt_selected_candidate.json"
    if not geometry_path.exists() or not selected_path.exists():
        return
    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    selected = json.loads(selected_path.read_text(encoding="utf-8"))
    geometry.setdefault("start_design", {}).setdefault(
        "fixed_calibration_parameters",
        selected["construction"]["calibration"],
    )
    geometry_path.write_text(
        json.dumps(geometry, indent=2, sort_keys=True), encoding="utf-8"
    )


def run_dynamic_rest_confirmation(
    output_directory: Path,
    *,
    starts: int = 7,
) -> dict[str, Any]:
    """Attach mandatory R1 beta/cAMP/PKA dynamics to the selected WT rest root.

    R1 has one effective beta/PKA/AE4 activation fraction.  At beta=0 its
    analytic basal equilibrium is exactly zero, a declared physical boundary
    rather than a fitted capacity bound.  Its 30-s time constant and the
    matching 1/30 s^-1 residual scale are sensitivity/modeling choices, not
    measured kinetics.
    """

    if starts < 3:
        raise ValueError("dynamic resting calibration requires at least three starts")
    selected_path = output_directory / "wt_selected_candidate.json"
    if not selected_path.exists():
        raise FileNotFoundError("R0 selected-candidate artifact is required first")
    selected = json.loads(selected_path.read_text(encoding="utf-8"))
    construction = selected["construction"]
    r0_variant = ModelVariant(**construction["variant"])
    calibration = WTCalibrationParameters(**construction["calibration"])
    dynamic_variant = replace(
        r0_variant,
        variant_id=r0_variant.variant_id.replace("G2_", "G3_R1_", 1),
        regulatory_family="R1_EFFECTIVE_BETA_PKA_AE4",
    )
    regulator = R1EffectiveActivation()
    spec = replace(
        WTCalibrationSpec(),
        regulatory_rhs_scales_s_inv=(1.0 / regulator.tau_activation_s,),
    )
    report = fit_wt_source_compatible_variant(
        dynamic_variant,
        spec=spec,
        regulatory_model=regulator,
        fixed_conductance_scale=calibration.common_membrane_conductance_scale,
        fixed_other_impermeant_start_fmol=(
            calibration.cell_other_impermeant_osmoles_fmol
        ),
        start_count=starts,
        seed=13131,
    )
    (output_directory / "wt_dynamic_rest_calibration.json").write_text(
        json.dumps(report_to_dict(report), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    passing = [branch for branch in report.branches if branch.passes_wt_gate]
    reference = (
        min(
            passing,
            key=lambda branch: (
                -len(branch.member_start_ids),
                max(abs(value) for value in branch.scaled_independent_rhs),
                branch.branch_id,
            ),
        )
        if passing
        else None
    )
    dynamic_candidate = {
        "status": (
            "R1_DYNAMIC_WT_RESTING_REFERENCE_NOT_FULLY_VALIDATED"
            if reference is not None
            else "R1_DYNAMIC_WT_RESTING_ROOT_FAILED"
        ),
        "regulatory_family": "R1 effective beta/PKA/AE4 activation",
        "kinetic_identification": "UNMEASURED_30_S_SENSITIVITY_NOT_FITTED",
        "regulatory_rhs_scale_s_inv": 1.0 / regulator.tau_activation_s,
        "source_supported_boundary_rule": (
            "The beta-free activation fraction is fixed by the analytic R1 equilibrium at zero; "
            "it is reported as a boundary but is not an optimizer-driven bound failure."
        ),
        "passing_root_ids": [branch.branch_id for branch in passing],
        "all_branches": [asdict(branch) for branch in report.branches],
        "construction": (
            {
                "variant": asdict(dynamic_variant),
                "calibration": dict(reference.parameters),
                "regulatory_model": {
                    "class": type(regulator).__name__,
                    "tau_activation_s": regulator.tau_activation_s,
                    "fully_activated_increment": regulator.gain.fully_activated_increment,
                },
                "resting_calcium_uM": spec.resting_calcium_uM,
                "root_id": reference.branch_id,
                "root_state": list(reference.state),
                "root_coordinates": list(reference.coordinates),
            }
            if reference is not None else None
        ),
        "full_model_validated": False,
        "reason": "A dynamic resting root does not establish WT stimulation timing or flow shape.",
    }
    (output_directory / "wt_dynamic_reference_candidate.json").write_text(
        json.dumps(dynamic_candidate, indent=2, sort_keys=True), encoding="utf-8"
    )

    generation_path = output_directory / "wt_model_generations.json"
    if generation_path.exists():
        generations = json.loads(generation_path.read_text(encoding="utf-8"))
        generations["generations"] = [
            item
            for item in generations.get("generations", [])
            if item.get("generation_id") != "G3_R1_DYNAMIC_REGULATION"
        ]
        generations["generations"].append(
            {
                "generation_id": "G3_R1_DYNAMIC_REGULATION",
                "parent": "G2b_CORRECTED_PALK_WATER_WT_REST",
                "defect_addressed": "mandatory dynamic beta/cAMP/PKA regulation of AE4",
                "changed_equations": [
                    "d(ae4_effective_activation_fraction)/dt=(beta-input-activation)/tau"
                ],
                "new_parameters": [
                    "tau_activation_s=30 sensitivity, not identified",
                    "fully_activated_increment=0.25 sensitivity, not fitted to secretion",
                ],
                "calibration_data": ["WT resting core only at beta=0"],
                "prohibited_data_used": [],
                "passing_root_ids": [branch.branch_id for branch in passing],
                "decision": (
                    "RETAIN_DYNAMIC_RESTING_SCAFFOLD_STIMULATION_VALIDATION_REQUIRED"
                    if passing else "REJECT_NO_VALID_DYNAMIC_RESTING_ROOT"
                ),
            }
        )
        generation_path.write_text(
            json.dumps(generations, indent=2, sort_keys=True), encoding="utf-8"
        )

    summary_path = output_directory / "wt_calibration_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["dynamic_R1_rest"] = {
            "passing_root_count": len(passing),
            "passing_root_ids": [branch.branch_id for branch in passing],
            "reference_root_id": reference.branch_id if reference is not None else None,
        }
        if reference is not None:
            summary["classification"] = (
                "WT_RESTING_ROOTS_FOUND_WITH_R1_DYNAMIC_REGULATION_"
                "WT_STIMULATION_VALIDATION_PENDING"
            )
            summary["selected_dynamic_reference_root_id"] = reference.branch_id
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
        )
    _upgrade_root_table_boundary_schema(output_directory / "wt_root_table.csv")
    _upgrade_selected_geometry_provenance(output_directory)
    return dynamic_candidate


def run_parameter_ablations(
    output_directory: Path,
    *,
    starts: int = 3,
) -> list[Mapping[str, Any]]:
    """Write the generation-local three-parameter minimality profile."""

    selected_path = output_directory / "wt_selected_candidate.json"
    if not selected_path.exists():
        raise FileNotFoundError("R0 selected-candidate artifact is required first")
    selected = json.loads(selected_path.read_text(encoding="utf-8"))
    construction = selected["construction"]
    variant = ModelVariant(**construction["variant"])
    calibration = WTCalibrationParameters(**construction["calibration"])
    rows = list(
        profile_wt_parameter_ablations(
            variant,
            calibration,
            construction["root_coordinates"],
            start_count=starts,
        )
    )
    (output_directory / "wt_parameter_ablations.json").write_text(
        json.dumps(rows, indent=2, sort_keys=True), encoding="utf-8"
    )
    fields = (
        "ablation_id", "removed_fit_parameter", "held_fixed_value",
        "held_fixed_unit", "active_fitted_parameters", "optimizer_success",
        "exact_root", "max_abs_dimensionless_residual", "boundary_hits",
        "normalized_jacobian_rank", "normalized_jacobian_column_count",
        "nfev", "cost", "representative_start_id", "interpretation_scope",
    )
    with (output_directory / "wt_parameter_ablations.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            flattened = {name: row.get(name, "") for name in fields}
            flattened["active_fitted_parameters"] = ";".join(
                row["active_fitted_parameters"]
            )
            flattened["boundary_hits"] = ";".join(row["boundary_hits"])
            writer.writerow(flattened)
    return rows


def synchronize_wt_dynamic_gate_summary(output_directory: Path) -> dict[str, Any]:
    """Fold independently produced WT-dynamic gate status into WT summary."""

    summary_path = output_directory / "wt_calibration_summary.json"
    gate_path = output_directory / "wt_dynamic_gate.json"
    if not summary_path.exists() or not gate_path.exists():
        raise FileNotFoundError("WT calibration summary and dynamic gate are required")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    summary["wt_dynamic_gate"] = {
        "status": gate["status"],
        "numerical_dynamic_gate_pass": gate["numerical_dynamic_gate_pass"],
        "minute_wise_flow_shape_gate_pass": gate[
            "minute_wise_flow_shape_gate_pass"
        ],
        "pre_reveal_holdout_eligible": gate["pre_reveal_holdout_eligible"],
        "effective_cell_count": gate["effective_cell_count"],
        "effective_cell_count_gate_pass": gate["effective_cell_count_gate_pass"],
        "artifact": "wt_dynamic_gate.json",
    }
    if not gate["minute_wise_flow_shape_gate_pass"]:
        summary["classification"] = (
            "WT_RESTING_ROOTS_FOUND_WITH_DYNAMIC_REGULATION_"
            "WT_STIMULATION_FLOW_SHAPE_FAILED"
        )
        summary["reason"] = (
            "The frozen 20-member R0-R3 ensemble passes numerical dynamics but "
            "fails the predeclared scale-free WT minute-flow-shape gate; heldout "
            "reveal remains ineligible."
        )
    ablation_path = output_directory / "wt_parameter_ablations.json"
    if ablation_path.exists():
        ablations = json.loads(ablation_path.read_text(encoding="utf-8"))
        summary["three_parameter_ablation_exact_root_count"] = sum(
            bool(row["exact_root"]) for row in ablations
        )
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )

    generation_path = output_directory / "wt_model_generations.json"
    if generation_path.exists():
        generations = json.loads(generation_path.read_text(encoding="utf-8"))
        generations["generations"] = [
            item
            for item in generations.get("generations", [])
            if item.get("generation_id") != "G5_FROZEN_WT_DYNAMIC_BASELINE"
        ]
        generations["generations"].append(
            {
                "generation_id": "G5_FROZEN_WT_DYNAMIC_BASELINE",
                "parent": "G3_R1_DYNAMIC_REGULATION",
                "defect_addressed": "WT CCh/IPR stimulation and physical-seconds dynamics",
                "regulatory_ensemble_member_count": gate[
                    "regulatory_ensemble_member_count"
                ],
                "numerical_dynamic_gate_pass": gate[
                    "numerical_dynamic_gate_pass"
                ],
                "minute_wise_flow_shape_gate_pass": gate[
                    "minute_wise_flow_shape_gate_pass"
                ],
                "pre_reveal_holdout_eligible": gate[
                    "pre_reveal_holdout_eligible"
                ],
                "prohibited_data_used": [],
                "decision": "REJECT_BASELINE_SHAPE; RETAIN_AS_NEGATIVE_CONTROL_FOR_WT_ONLY_REPAIR",
            }
        )
        generation_path.write_text(
            json.dumps(generations, indent=2, sort_keys=True), encoding="utf-8"
        )
    return summary


def run(
    output_directory: Path,
    *,
    screen_starts: int = 3,
    confirm_starts: int = 7,
    geometry_starts: int = 32,
    workers: int = 4,
) -> dict[str, Any]:
    if min(screen_starts, confirm_starts, geometry_starts) < 3:
        raise ValueError("every multistart stage requires at least three starts")
    if workers < 1:
        raise ValueError("workers must be positive")
    output_directory.mkdir(parents=True, exist_ok=True)
    screen_tasks: list[tuple[ModelVariant, float, int, int, str]] = []
    for scale in (10.0, 30.0, 100.0):
        for index, base_variant in enumerate(default_variant_panel()):
            variant = _variant_with_scale(base_variant, scale)
            screen_tasks.append(
                (
                    variant,
                    scale,
                    screen_starts,
                    13023 + index + int(scale),
                    "screen",
                )
            )
    reports = _run_profile_tasks(screen_tasks, workers=workers)

    # Seven-start confirmation is reserved for variants that produce a
    # physiological interior root in the complete three-start screen.  Failed
    # conductance/routing/topology slices remain explicit in the artifact.
    confirm_indices = [
        index
        for index, report in enumerate(reports)
        if any(branch.passes_wt_gate for branch in report.branches)
    ]
    confirm_tasks = [
        (
            reports[index].variant,
            float(reports[index].start_design["fixed_conductance_scale"]),
            confirm_starts,
            int(reports[index].start_design["seed"]),
            "confirm",
        )
        for index in confirm_indices
    ]
    for index, confirmed in zip(
        confirm_indices,
        _run_profile_tasks(confirm_tasks, workers=workers),
    ):
        reports[index] = confirmed

    _write_attempts(output_directory / "wt_root_attempts.csv", reports)
    _write_roots(output_directory / "wt_root_table.csv", reports)
    _write_states_and_parameters(
        output_directory / "wt_root_states.csv",
        output_directory / "wt_parameter_candidates.csv",
        reports,
    )
    (output_directory / "wt_solver_diagnostics.json").write_text(
        json.dumps([report_to_dict(report) for report in reports], indent=2, sort_keys=True),
        encoding="utf-8",
    )
    parameter_ledger = [asdict(record) for record in parameter_records(FullModelParameters())]
    (output_directory / "parameter_provenance.json").write_text(
        json.dumps(parameter_ledger, indent=2, sort_keys=True), encoding="utf-8"
    )
    passing_pairs = [
        (report, branch)
        for report in reports
        for branch in report.branches
        if branch.passes_wt_gate
    ]
    passing = [branch for _, branch in passing_pairs]

    selected_payload: dict[str, Any] | None = None
    geometry_report: RootSearchReport | None = None
    if passing_pairs:
        preferred_pairs = [
            pair
            for pair in passing_pairs
            if pair[0].variant.variant_id == "G2_BALANCED_APICAL_K_BIASED_G10"
        ]
        candidates = preferred_pairs or passing_pairs
        selected_report, selected_branch = min(
            candidates,
            key=lambda pair: (
                -len(pair[1].member_start_ids),
                max(abs(value) for value in pair[1].scaled_independent_rhs),
                pair[0].variant.variant_id,
                pair[1].branch_id,
            ),
        )
        selected_calibration = WTCalibrationParameters(**selected_branch.parameters)
        selected_model = build_wt_model(selected_report.variant, selected_calibration)
        geometry_report = solve_wt_resting_branches(
            selected_model,
            variant=selected_report.variant,
            calibration=selected_calibration,
            start_count=geometry_starts,
            seed=13091,
            reference_coordinates=selected_branch.coordinates,
        )
        (output_directory / "wt_selected_root_geometry.json").write_text(
            json.dumps(report_to_dict(geometry_report), indent=2, sort_keys=True),
            encoding="utf-8",
        )
        selected_payload = {
            "status": "R0_WT_RESTING_REFERENCE_ONLY_NOT_FINAL_FULL_MODEL",
            "selection_rule": (
                "Prefer the predeclared balanced-apical, K-biased, G=10 slice; "
                "then maximize converged-start membership and minimize the dimensionless residual."
            ),
            "construction": {
                "variant": asdict(selected_report.variant),
                "calibration": asdict(selected_calibration),
                "resting_calcium_uM": WTCalibrationSpec().resting_calcium_uM,
                "root_id": selected_branch.branch_id,
                "root_state": list(selected_branch.state),
                "root_coordinates": list(selected_branch.coordinates),
            },
            "calibration_branch": asdict(selected_branch),
            "fixed_parameter_root_geometry": {
                "start_count": len(geometry_report.attempts),
                "distinct_root_count": len(geometry_report.branches),
                "passing_root_ids": [
                    branch.branch_id
                    for branch in geometry_report.branches
                    if branch.passes_wt_gate
                ],
            },
            "unresolved_profile_axes": [
                "native mixed-bath Na/K routing",
                "apical pump/K topology fractions",
                "common membrane conductance scale",
            ],
        }
        (output_directory / "wt_selected_candidate.json").write_text(
            json.dumps(selected_payload, indent=2, sort_keys=True), encoding="utf-8"
        )
    generation = {
        "generations": [
            {
                "generation_id": "G2a_PRE_AUDIT_WATER_PLACEHOLDERS",
                "parent": "G1_CONSERVATION_EXPLICIT_CORE",
                "defect_addressed": "initial WT rest closure",
                "water_coefficients_pL_s_mOsm": [2.0e-5, 2.0e-5, 1.0e-5],
                "decision": "REJECT",
                "reason": "dimensional audit found coefficients inconsistent with Palk/2018 conversion",
            },
            {
                "generation_id": "G2b_CORRECTED_PALK_WATER_WT_REST",
                "parent": "G2a_PRE_AUDIT_WATER_PLACEHOLDERS",
                "defect_addressed": "water-permeability unit conversion and WT resting closure",
                "water_coefficients_pL_s_mOsm": [0.00432, 0.0515, 0.000260],
                "changed_equations": [],
                "changed_parameters": [
                    "water.apical_hydraulic_pL_s_mOsm",
                    "water.basolateral_hydraulic_pL_s_mOsm",
                    "water.paracellular_hydraulic_pL_s_mOsm",
                ],
                "calibration_data": [
                    "matched WT resting Cl", "matched WT resting pH",
                    "historical-lineage 1.3-pL volume convention",
                ],
                "prohibited_data_used": [],
                "passing_R0_root_ids": [branch.branch_id for branch in passing],
                "decision": "RETAIN_R0_REST_SCAFFOLD_ONLY_DYNAMIC_REGULATION_PENDING",
            },
        ]
    }
    (output_directory / "wt_model_generations.json").write_text(
        json.dumps(generation, indent=2, sort_keys=True), encoding="utf-8"
    )
    summary = {
        "classification": (
            "WT_RESTING_ROOTS_FOUND_R0_ONLY_DYNAMIC_REGULATION_PENDING"
            if passing
            else "NO_VALID_WT_RESTING_ROOT"
        ),
        "variant_profiles": len(reports),
        "screen_starts_per_profile": screen_starts,
        "confirm_starts_per_passing_profile": confirm_starts,
        "selected_fixed_parameter_geometry_starts": (
            len(geometry_report.attempts) if geometry_report is not None else 0
        ),
        "parallel_workers": workers,
        "passing_root_count": len(passing),
        "passing_root_ids": [branch.branch_id for branch in passing],
        "selected_reference_root_id": (
            selected_payload["construction"]["root_id"]
            if selected_payload is not None else None
        ),
        "heldout_access": "NONE",
        "final_full_model_validated": False,
        "reason": "Dynamic beta/cAMP/PKA integration and WT stimulation remain pending.",
    }
    (output_directory / "wt_calibration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    if selected_payload is not None:
        ablations = run_parameter_ablations(output_directory, starts=screen_starts)
        summary["three_parameter_ablation_exact_root_count"] = sum(
            bool(row["exact_root"]) for row in ablations
        )
        (output_directory / "wt_calibration_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
        )
        dynamic_candidate = run_dynamic_rest_confirmation(
            output_directory, starts=confirm_starts
        )
        summary = json.loads(
            (output_directory / "wt_calibration_summary.json").read_text(
                encoding="utf-8"
            )
        )
        summary["dynamic_R1_rest_status"] = dynamic_candidate["status"]
        (output_directory / "wt_calibration_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
        )
    if (output_directory / "wt_dynamic_gate.json").exists():
        summary = synchronize_wt_dynamic_gate_summary(output_directory)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    parser.add_argument("--screen-starts", type=int, default=3)
    parser.add_argument("--confirm-starts", type=int, default=7)
    parser.add_argument("--geometry-starts", type=int, default=32)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--dynamic-only", action="store_true")
    parser.add_argument("--ablations-only", action="store_true")
    args = parser.parse_args()
    if args.ablations_only:
        print(
            json.dumps(
                run_parameter_ablations(args.output, starts=args.screen_starts),
                indent=2,
                sort_keys=True,
            )
        )
        return
    if args.dynamic_only:
        print(
            json.dumps(
                run_dynamic_rest_confirmation(
                    args.output, starts=args.confirm_starts
                ),
                indent=2,
                sort_keys=True,
            )
        )
        return
    print(
        json.dumps(
            run(
                args.output,
                screen_starts=args.screen_starts,
                confirm_starts=args.confirm_starts,
                geometry_starts=args.geometry_starts,
                workers=args.workers,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
