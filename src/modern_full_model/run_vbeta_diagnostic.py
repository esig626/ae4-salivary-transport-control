"""Run the exact-rest Vbeta anti-deadlock audit on retained native WT roots.

This is a strictly pre-reveal WT/IPR-only diagnostic.  It does not evaluate a
genotype perturbation and it does not read a validation or held-out target
ledger.  V1 is evaluated only as a gate on trajectories generated *without*
V1, preventing the candidate pathway from creating its own activation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .camp_pka import R1EffectiveActivation
from .native_source_panel import (
    build_native_source_model,
    load_native_source_roots,
)
from .validation import (
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    StimulusArm,
)
from .vbeta_diagnostic import assess_pre_vbeta_trajectory


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_rows(path: Path, rows: list[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError("Vbeta diagnostic produced no rows")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_native_vbeta_diagnostic(results_directory: str | Path) -> Mapping[str, Any]:
    """Audit every numerical native root with both frozen stiff solvers.

    WT-failing resting roots remain useful only for testing whether the V1
    gate algebra deadlocks.  They are labelled ``REST_WT_FAIL`` in every row
    and can never authorize production or held-out reveal.
    """

    output = Path(results_directory)
    roots_path = output / "native_source_wt_roots.csv"
    roots = tuple(
        root
        for root in load_native_source_roots(roots_path)
        if root.passes_numerical_gate
        and root.eligibility == "PRODUCTION_PRE_REVEAL"
    )
    if not roots:
        raise ValueError("native-source table contains no numerical production-panel roots")

    regulator = R1EffectiveActivation(tau_activation_s=30.0)
    protocol = SecretagogueProtocol(
        arm=StimulusArm.IPR_ONLY,
        resting_calcium_uM=0.058,
        # In the IPR-only arm the stimulated calcium value is unused; retaining
        # the declared reference makes the protocol identity explicit.
        stimulated_calcium_uM=0.55,
    )
    time_s = np.linspace(0.0, 600.0, 301)
    rows: list[Mapping[str, Any]] = []
    for root in roots:
        model = build_native_source_model(
            root,
            regulatory_model=regulator,
            stimulus=protocol,
        )
        initial = np.r_[
            np.asarray(root.core_state, dtype=float),
            np.asarray(regulator.initial_state(beta_input=0.0), dtype=float),
        ]
        if abs(initial[5] - root.coordinates[4]) > 1.0e-12:
            raise ValueError(f"root {root.root_id} has inconsistent cell volumes")
        for solver in (PRODUCTION_RADAU, PRODUCTION_BDF):
            solution = model.solve_dynamics(
                (0.0, 600.0),
                initial_state=initial,
                method=solver.method,
                rtol=solver.rtol,
                atol=solver.atol_vector(model.state_names),
                t_eval=time_s,
                max_step_s=solver.max_step_s,
            )
            beta = np.asarray([protocol.beta_input(t) for t in solution.t])
            assessment = assess_pre_vbeta_trajectory(
                solution.y[5],
                beta,
                reference_cell_volume_pL=float(initial[5]),
            )
            rows.append(
                {
                    "root_id": root.root_id,
                    "rest_wt_gate_status": (
                        "WT_PASS_RETAINED" if root.passes_wt_gate else "REST_WT_FAIL"
                    ),
                    "root_production_eligible_before_vbeta": bool(root.passes_wt_gate),
                    "source_class": root.source_class,
                    "source_scale": root.source_scale,
                    "topology_id": root.topology_id,
                    "routing_id": root.routing_id,
                    "solver": solver.method,
                    "solver_success": bool(solution.success),
                    "initial_volume_pL": float(solution.y[5, 0]),
                    "final_volume_pL": float(solution.y[5, -1]),
                    **asdict(assessment),
                    "rejection_reasons": json.dumps(assessment.rejection_reasons),
                }
            )

    csv_path = output / "vbeta_native_exact_rest_diagnostic.csv"
    _write_rows(csv_path, rows)
    no_gate_rows = [row for row in rows if not row["activation_detected"]]
    all_solvers_succeeded = all(bool(row["solver_success"]) for row in rows)
    all_no_gate = len(no_gate_rows) == len(rows)
    rows_by_root: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        rows_by_root.setdefault(str(row["root_id"]), []).append(row)
    if any(len(group) != 2 for group in rows_by_root.values()):
        raise ValueError("each native root must have exactly two solver rows")
    solver_gate_classification_agreement = all(
        bool(group[0]["activation_detected"])
        == bool(group[1]["activation_detected"])
        for group in rows_by_root.values()
    )
    maximum_solver_swelling_difference = max(
        abs(
            float(group[0]["maximum_relative_swelling"])
            - float(group[1]["maximum_relative_swelling"])
        )
        for group in rows_by_root.values()
    )
    maximum_solver_final_volume_relative_difference = max(
        abs(float(group[0]["final_volume_pL"]) - float(group[1]["final_volume_pL"]))
        / max(
            abs(float(group[0]["final_volume_pL"])),
            abs(float(group[1]["final_volume_pL"])),
            1.0e-12,
        )
        for group in rows_by_root.values()
    )
    wt_passing_root_count = sum(bool(root.passes_wt_gate) for root in roots)
    rest_wt_failing_root_count = len(roots) - wt_passing_root_count
    activated_root_ids = sorted(
        {
            str(row["root_id"])
            for row in rows
            if bool(row["activation_detected"])
        }
    )
    wt_passing_activated_root_ids = sorted(
        {
            str(row["root_id"])
            for row in rows
            if bool(row["activation_detected"])
            and row["rest_wt_gate_status"] == "WT_PASS_RETAINED"
        }
    )
    if all_solvers_succeeded and all_no_gate and wt_passing_root_count:
        status = "V1_REJECTED_EXACT_REST_ANTI_DEADLOCK_ALL_WT_PASSING_ROOTS"
    elif all_solvers_succeeded and all_no_gate:
        status = (
            "V1_DEADLOCKS_ON_NUMERICAL_ROOTS_BUT_NO_WT_PASSING_ROOT_EXISTS"
        )
    elif all_solvers_succeeded and not wt_passing_root_count:
        status = "V1_ACTIVATES_ONLY_ON_REST_WT_FAIL_ROOTS_NO_PRODUCTION_CONCLUSION"
    else:
        status = "V1_GATE_BEHAVIOR_MIXED_CAPACITY_REMAINS_UNLICENSED"
    summary = {
        "status": status,
        "protocol": {
            "arm": "IPR_ONLY",
            "calcium_uM": 0.058,
            "beta_input": "0 at exact t=0; 1 for 0<t<=600 s",
            "regulation": "R1_EFFECTIVE_TAU_30_S_GAIN_1P25",
        },
        "topology": "gV * beta_input * max(Vi/Vrest - 1, 0)",
        "capacity_status": "UNLICENSED_NO_WHOLE_CELL_CURRENT_MAP",
        "production_eligible": False,
        "root_count": len(roots),
        "wt_passing_root_count": wt_passing_root_count,
        "rest_wt_failing_root_count": rest_wt_failing_root_count,
        "activated_root_count": len(activated_root_ids),
        "activated_root_ids": activated_root_ids,
        "wt_passing_activated_root_count": len(wt_passing_activated_root_ids),
        "wt_passing_activated_root_ids": wt_passing_activated_root_ids,
        "row_count": len(rows),
        "source_classes": sorted({root.source_class for root in roots}),
        "all_solvers_succeeded": all_solvers_succeeded,
        "solver_gate_classification_agreement": (
            solver_gate_classification_agreement
        ),
        "maximum_solver_swelling_absolute_difference": (
            maximum_solver_swelling_difference
        ),
        "maximum_solver_final_volume_relative_difference": (
            maximum_solver_final_volume_relative_difference
        ),
        "all_rows_no_gate_activation": all_no_gate,
        "minimum_relative_swelling": min(
            float(row["minimum_relative_swelling"]) for row in rows
        ),
        "maximum_relative_swelling": max(
            float(row["maximum_relative_swelling"]) for row in rows
        ),
        "maximum_rectified_gate": max(
            float(row["maximum_rectified_gate"]) for row in rows
        ),
        "firewall": "WT_IPR_ONLY_NO_GENOTYPE_PERTURBATION_NO_HELDOUT_READ",
        "input_hashes": {roots_path.name: _sha256(roots_path)},
    }
    summary_path = output / "vbeta_native_exact_rest_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    args = parser.parse_args()
    print(
        json.dumps(
            run_native_vbeta_diagnostic(args.results_directory),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
