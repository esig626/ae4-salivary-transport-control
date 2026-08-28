"""Independent native-source AE2 continuation and dynamic validation.

Only production-eligible H1 WT roots are admitted.  AE2 expression is
continued from one to exact zero through 21 fixed steps while every parameter
and fitted OTHER osmole remains frozen.  The connected null root is then
propagated through the frozen R3/G1.25/reference-time, CA010, N1/M1.75
CCh+IPR protocol with Radau and BDF.  No AE4-null phenotype target is read.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .independent import (
    IndependentConstantStimulus,
    IndependentN1NkccRegulation,
    IndependentSecretagogueProtocol,
    independent_physical_time_grid,
    independent_r3_regulation,
)
from .independent_native import IndependentNativeSourceDefinition


EXPRESSION_GRID = np.linspace(1.0, 0.0, 21)
INDEPENDENT_ROW_SCALES = np.asarray(
    (0.05, 0.05, 0.05, 0.02, 1.0e-4, 0.02, 0.02, 0.02, 0.02, 1.0e-4)
)


def _truth(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _eligible(row: dict[str, str]) -> bool:
    mode = row.get("hydraulic_mode") or (
        "H1" if float(row["hydraulic_scale"]) == 1.0 else "HW"
    )
    return bool(
        _truth(row["passes_numerical_gate"])
        and _truth(row["passes_wt_gate"])
        and row["eligibility"] == "PRODUCTION_PRE_REVEAL"
        and float(row["hydraulic_scale"]) == 1.0
        and mode == "H1"
    )


def _definition(row: dict[str, str]) -> IndependentNativeSourceDefinition:
    return IndependentNativeSourceDefinition.from_root_row(row)


def _atol(state_names: tuple[str, ...]) -> np.ndarray:
    return np.asarray(
        [
            1.0e-10
            if name.endswith("_fmol")
            else 1.0e-12
            if name.endswith("_pL")
            else 1.0e-10
            for name in state_names
        ]
    )


def run_independent_ae2(
    roots_path: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    with roots_path.open(newline="", encoding="utf-8") as handle:
        roots = [row for row in csv.DictReader(handle) if _eligible(row)]

    continuation_rows: list[dict[str, Any]] = []
    dynamic_rows: list[dict[str, Any]] = []
    for row in roots:
        definition = _definition(row)
        resting_model = definition.build_independent_model(
            stimulus=IndependentConstantStimulus(0.058, 0.0)
        )
        state = np.asarray(json.loads(row["core_state_json"]), dtype=float)
        completed = True
        for expression in EXPRESSION_GRID:
            result = resting_model.refine_resting_root(
                state,
                relative_bound=0.05,
                max_nfev=500,
                row_scales=INDEPENDENT_ROW_SCALES,
                ae2_expression=float(expression),
            )
            evaluation = resting_model.evaluate(
                0.0,
                result.state,
                ae2_expression=float(expression),
            )
            completed = completed and result.success
            continuation_rows.append(
                {
                    "root_id": row["root_id"],
                    "source_class": row["source_class"],
                    "source_scale": float(row["source_scale"]),
                    "ae2_expression": float(expression),
                    "independent_root_success": result.success,
                    "independent_scaled_residual": result.max_abs_scaled_residual,
                    "independent_jacobian_rank": result.jacobian_rank,
                    "independent_jacobian_nullity": result.jacobian_nullity,
                    "ae2_flux_fmol_s": evaluation.fluxes_fmol_s["ae2"],
                    "cell_na_mM": evaluation.cell_concentrations_mM["na"],
                    "cell_k_mM": evaluation.cell_concentrations_mM["k"],
                    "cell_cl_mM": evaluation.cell_concentrations_mM["cl"],
                    "cell_ph": evaluation.cell_speciation.ph,
                    "cell_volume_pL": float(result.state[5]),
                }
            )
            if not result.success:
                completed = False
                break
            state = np.asarray(result.state[:12], dtype=float)

        if not completed or len(
            [item for item in continuation_rows if item["root_id"] == row["root_id"]]
        ) != len(EXPRESSION_GRID):
            continue
        null_state = state
        regulation = independent_r3_regulation(
            tau_pka_s=10.0,
            forward_regulation_rate_s=1.0 / 30.0,
            reverse_regulation_rate_s=1.0 / 30.0,
            basal_multiplier=1.0,
            fully_activated_increment=0.25,
        )
        nkcc = IndependentN1NkccRegulation(
            fully_activated_multiplier=1.75,
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.10,
        )
        protocol = IndependentSecretagogueProtocol(
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.10,
            beta_occupancy_on=1.0,
        )
        dynamic_model = definition.build_independent_model(
            regulation=regulation,
            nkcc_regulation=nkcc,
            stimulus=protocol,
        )
        initial_wt = np.r_[
            np.asarray(json.loads(row["core_state_json"]), dtype=float),
            regulation.basal_state,
        ]
        initial_ae2 = np.r_[null_state, regulation.basal_state]
        for method in ("Radau", "BDF"):
            trajectories = {
                "WT": dynamic_model.integrate(
                    initial_wt,
                    independent_physical_time_grid(600.0, 5.0),
                    method=method,
                    rtol=1.0e-7,
                    atol=_atol(dynamic_model.state_names),
                    max_step_s=2.0,
                    preserve_basal_zero=True,
                    ae2_expression=1.0,
                ),
                "AE2_NULL": dynamic_model.integrate(
                    initial_ae2,
                    independent_physical_time_grid(600.0, 5.0),
                    method=method,
                    rtol=1.0e-7,
                    atol=_atol(dynamic_model.state_names),
                    max_step_s=2.0,
                    preserve_basal_zero=True,
                    ae2_expression=0.0,
                ),
            }
            wt = trajectories["WT"]
            null = trajectories["AE2_NULL"]
            deleted_flux_max = max(
                abs(
                    dynamic_model.evaluate(
                        float(time_s),
                        null.states[:, index],
                        ae2_expression=0.0,
                    ).fluxes_fmol_s["ae2"]
                )
                for index, time_s in enumerate(null.time_s)
            )
            dynamic_rows.append(
                {
                    "root_id": row["root_id"],
                    "source_class": row["source_class"],
                    "source_scale": float(row["source_scale"]),
                    "solver_method": method,
                    "wt_success": wt.success,
                    "ae2_null_success": null.success,
                    "wt_cumulative_flow_600_pL": float(wt.cumulative_flow_pL[-1]),
                    "ae2_null_cumulative_flow_600_pL": float(
                        null.cumulative_flow_pL[-1]
                    ),
                    "ae2_to_wt_cumulative_flow_ratio": float(
                        null.cumulative_flow_pL[-1] / wt.cumulative_flow_pL[-1]
                    ),
                    "ae2_to_wt_endpoint_flow_ratio": float(
                        null.flow_pL_s[-1] / wt.flow_pL_s[-1]
                    ),
                    "wt_endpoint_cl_mM": float(wt.cell_cl_mM[-1]),
                    "ae2_null_endpoint_cl_mM": float(null.cell_cl_mM[-1]),
                    "wt_endpoint_ph": float(wt.cell_ph[-1]),
                    "ae2_null_endpoint_ph": float(null.cell_ph[-1]),
                    "maximum_deleted_ae2_flux_fmol_s": deleted_flux_max,
                    "wt_max_conservation_fmol_s": wt.max_abs_conservation_fmol_s,
                    "ae2_null_max_conservation_fmol_s": (
                        null.max_abs_conservation_fmol_s
                    ),
                }
            )

    complete_root_ids = {
        row["root_id"]
        for row in continuation_rows
        if row["ae2_expression"] == 0.0
        and row["independent_root_success"]
        and row["ae2_flux_fmol_s"] == 0.0
    }
    dynamic_by_root: dict[str, list[dict[str, Any]]] = {}
    for row in dynamic_rows:
        dynamic_by_root.setdefault(row["root_id"], []).append(row)
    solver_differences = []
    for root_id, cases in dynamic_by_root.items():
        if len(cases) == 2:
            values = {case["solver_method"]: case for case in cases}
            solver_differences.append(
                abs(
                    values["Radau"]["ae2_to_wt_cumulative_flow_ratio"]
                    - values["BDF"]["ae2_to_wt_cumulative_flow_ratio"]
                )
            )
    summary = {
        "status": (
            "INDEPENDENT_NATIVE_AE2_PASS"
            if roots
            and len(complete_root_ids) == len(roots)
            and len(dynamic_rows) == 2 * len(roots)
            and all(row["maximum_deleted_ae2_flux_fmol_s"] == 0.0 for row in dynamic_rows)
            else "NO_ELIGIBLE_NATIVE_WT_ROOTS"
            if not roots
            else "INDEPENDENT_NATIVE_AE2_DISCREPANCY"
        ),
        "firewall": "PRE_REVEAL_AE2_VALIDATION_NO_AE4_NULL_TARGET_ACCESS",
        "eligible_wt_root_count": len(roots),
        "completed_ae2_null_continuation_count": len(complete_root_ids),
        "continuation_step_count": len(continuation_rows),
        "dynamic_case_count": len(dynamic_rows),
        "maximum_exact_deleted_ae2_flux_fmol_s": (
            max(
                (row["maximum_deleted_ae2_flux_fmol_s"] for row in dynamic_rows),
                default=None,
            )
        ),
        "ae2_to_wt_cumulative_flow_ratio_range": (
            [
                min(row["ae2_to_wt_cumulative_flow_ratio"] for row in dynamic_rows),
                max(row["ae2_to_wt_cumulative_flow_ratio"] for row in dynamic_rows),
            ]
            if dynamic_rows
            else None
        ),
        "maximum_radau_bdf_ratio_absolute_difference": (
            max(solver_differences) if solver_differences else None
        ),
    }
    return continuation_rows, dynamic_rows, summary


def _write_csv(path: Path, rows: list[dict[str, Any]], empty_fields: tuple[str, ...]) -> None:
    fields = list(rows[0]) if rows else list(empty_fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    output = args.repository / "results" / "13B_modern_full_model"
    continuation, dynamics, summary = run_independent_ae2(
        output / "native_source_wt_roots.csv"
    )
    _write_csv(
        output / "independent_native_ae2_continuation.csv",
        continuation,
        ("root_id", "ae2_expression", "independent_root_success"),
    )
    _write_csv(
        output / "independent_native_ae2_dynamics.csv",
        dynamics,
        ("root_id", "solver_method", "ae2_to_wt_cumulative_flow_ratio"),
    )
    (output / "independent_native_ae2_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
