"""Reproduce transporter constraint and fixed-chassis Stage-A artifacts."""

from __future__ import annotations

from dataclasses import asdict
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping

from ae4_mechanism_reconstruction.chassis import AE4_KNOCKOUT, FixedChassis, ZeroAE4

from .cycles import CANDIDATE_DEFINITIONS, exact_charge_per_branch
from .evidence import verify_evidence_freeze
from .fitting import TransporterFitResult, transporter_fit_grid
from .solvers import (
    build_fixed_chassis,
    calibrate_wt_capacity,
    parameters_from_fit,
    solve_fixed_capacity_roots,
)


RESULT_DIR = Path("results/13_state_resolved_ae4")


def _write_csv(path: Path, rows: Iterable[Mapping[str, object]], fieldnames: list[str] | None = None) -> None:
    records = list(rows)
    if not records and fieldnames is None:
        raise ValueError(f"empty CSV needs explicit fields: {path}")
    fields = fieldnames or list(records[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(repository_root: str | Path = ".") -> Mapping[str, object]:
    root = Path(repository_root).resolve()
    boundary = verify_evidence_freeze(root)
    result_dir = root / RESULT_DIR
    result_dir.mkdir(parents=True, exist_ok=True)

    fits = transporter_fit_grid()
    fit_path = result_dir / "transporter_fit_summary.csv"
    _write_csv(fit_path, (row.as_dict() for row in fits))

    registry = []
    for model_id, definition in CANDIDATE_DEFINITIONS.items():
        registry.append(
            {
                "model_id": model_id,
                "family": definition.family,
                "state_count": len(definition.states),
                "edge_count": len(definition.edges),
                "branches": ";".join(definition.branch_current_edges),
                "branch_charge_equivalents": json.dumps(dict(exact_charge_per_branch(definition)), sort_keys=True),
                "pka_layered": bool(definition.pka_switch_edge_ids),
                "explicit_carbonate": "co3_i" in definition.required_reservoir_species,
                "evidence_status": definition.evidence_status,
                "new_assumptions": ";".join(definition.new_assumptions),
                "stage_a_holdout_read_before_freeze": False,
            }
        )
    _write_csv(result_dir / "model_registry.csv", registry)

    calibration_rows = [
        {"target_id": "E16-03", "stage": "CAL-T", "quantity": "direct Na transport", "used": True},
        {"target_id": "E16-04", "stage": "CAL-T", "quantity": "electroneutrality/voltage independence", "used": True},
        {"target_id": "E16-05", "stage": "CAL-T", "quantity": "direct K transport", "used": True},
        {"target_id": "E16-07", "stage": "CAL-T", "quantity": "Na/K EC50, Hill, Rmin/Rmax", "used": True},
        {"target_id": "E21-03", "stage": "CAL-T", "quantity": "approximately 25% recombinant forskolin response", "used": True},
        {"target_id": "E25-04", "stage": "CAL-T", "quantity": "double-mutant Na collapse/K retention", "used": True},
        {"target_id": "E15-05", "stage": "CAL-WT", "quantity": "WT resting Cl", "used": True},
        {"target_id": "E15-07_WT", "stage": "CAL-WT", "quantity": "WT resting pH gate", "used": False},
    ]
    _write_csv(result_dir / "calibration_targets.csv", calibration_rows)
    heldout_rows = [
        {"target_id": "E15-06", "stage": "VAL-A/LOC-B", "quantity": "AE4-null resting Cl", "read_before_model_freeze": False},
        {"target_id": "E15-07_KO", "stage": "VAL-A/LOC-B", "quantity": "AE4-null resting pH", "read_before_model_freeze": False},
        {"target_id": "E15-02", "stage": "STRICT", "quantity": "10-min secretion ratio", "read_before_model_freeze": False},
        {"target_id": "E15-03", "stage": "STRICT", "quantity": "secretion time course", "read_before_model_freeze": False},
    ]
    _write_csv(result_dir / "heldout_targets.csv", heldout_rows)

    provenance = [
        {"parameter": "EC50_Na_effective", "value": 49.0, "unit": "mM", "source": "E16-07", "status": "CAL-T regression summary"},
        {"parameter": "EC50_K_effective", "value": 62.0, "unit": "mM", "source": "E16-07", "status": "CAL-T regression summary"},
        {"parameter": "Hill_Na_effective", "value": 2.0, "unit": "dimensionless", "source": "E16-07", "status": "CAL-T regression summary; not site count"},
        {"parameter": "Hill_K_effective", "value": 1.8, "unit": "dimensionless", "source": "E16-07", "status": "CAL-T regression summary; not site count"},
        {"parameter": "assay_Rmin_Na", "value": 0.3, "unit": "1e-3 s^-1", "source": "E16-07", "status": "fixed assay observation-map nuisance; no zero-cation point"},
        {"parameter": "assay_Rmin_K", "value": 0.4, "unit": "1e-3 s^-1", "source": "E16-07", "status": "fixed assay observation-map nuisance; no zero-cation point"},
        {"parameter": "PKA_recombinant_fold", "value": 1.25, "unit": "fold", "source": "E21-03", "status": "assay-specific approximate CAL-T"},
        {"parameter": "common_barrier_gauge", "value": "0.1;1;10", "unit": "multiplier", "source": "NEW sensitivity", "status": "unidentified microscopic gauge; not fitted"},
        {"parameter": "whole_cell_capacity", "value": "per-row attempted", "unit": "chassis-native", "source": "E15-05", "status": "valid only if exact interior closure succeeds"},
    ]
    _write_csv(result_dir / "parameter_provenance.csv", provenance)

    survivors = [row for row in fits if row.transporter_gate_pass]
    calibration_attempts = []
    root_rows = []
    survivor_parameter_rows = []
    for fit in survivors:
        definition = CANDIDATE_DEFINITIONS[fit.model_id]
        attempt = calibrate_wt_capacity(definition, fit)
        calibration_attempts.append(attempt.as_dict())
        survivor_parameter_rows.append(
            {
                "model_id": fit.model_id,
                "gauge_id": fit.gauge_id,
                "capacity_attempt": attempt.capacity,
                "capacity_calibration_valid": attempt.calibration_valid,
                "na_loaded_energy_lump_mM": fit.na_loaded_energy_lump_mM,
                "k_loaded_energy_lump_mM": fit.k_loaded_energy_lump_mM,
                "k_attempt_scale": fit.k_attempt_scale,
                "pka_barrier_fold": fit.pka_barrier_fold,
                "common_barrier_multiplier": fit.common_barrier_multiplier,
                "heldout_target_read": False,
            }
        )
        parameters = parameters_from_fit(fit, capacity=attempt.capacity)
        chassis = build_fixed_chassis(definition, parameters)
        roots = solve_fixed_capacity_roots(
            chassis,
            model_id=fit.model_id,
            gauge_id=fit.gauge_id,
            capacity=attempt.capacity,
            preferred=attempt.state,
            start_count=10,
        )
        root_rows.extend(root.as_dict() for root in roots)
    _write_csv(result_dir / "wt_capacity_attempts.csv", calibration_attempts)
    _write_csv(result_dir / "frozen_survivor_parameters.csv", survivor_parameter_rows)
    _write_csv(
        result_dir / "stage_a_roots.csv",
        root_rows,
        fieldnames=list(root_rows[0]) if root_rows else [
            "model_id", "gauge_id", "scenario", "capacity", "root_id", "solver_success",
            "nfev", "cost", "max_abs_raw_rhs", "na_l_mM", "k_l_mM", "height_um",
            "na_i_mM", "k_i_mM", "cl_i_mM", "hco3_i_mM", "pH_i", "volume_pL",
            "cl_z_wt", "ph_z_wt", "positivity_pass", "wt_physiology_gate_pass",
            "alternate_root_count", "initial_start_count", "calibration_target_ids",
            "heldout_target_read", "failure_reason",
        ],
    )

    # Freeze hashes before any AE4-null phenotype value is joined.
    pre_reveal_files = [
        fit_path,
        result_dir / "wt_capacity_attempts.csv",
        result_dir / "frozen_survivor_parameters.csv",
        result_dir / "stage_a_roots.csv",
    ]
    reveal_marker = {
        "status": "transporter_and_WT_candidates_frozen_before_stage_a_reveal",
        "evidence_hash": boundary.evidence_hash,
        "question_tree_hash": boundary.question_tree_hash,
        "file_hashes": {str(path.relative_to(root)): _hash(path) for path in pre_reveal_files},
        "strict_secretion_used_in_any_fit": False,
    }
    (result_dir / "stage_a_reveal_marker.json").write_text(
        json.dumps(reveal_marker, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # Stage-A reveal: the AE4-null source is exactly zero for every candidate.
    # The root solve uses no phenotype target; target values are joined only in
    # the output comparison row below.
    zero_chassis = FixedChassis(ZeroAE4())
    null_roots = solve_fixed_capacity_roots(
        zero_chassis,
        model_id="AE4_NULL_COMMON",
        gauge_id="candidate_independent",
        capacity=0.0,
        scenario=AE4_KNOCKOUT,
        start_count=20,
    )
    null_rows = []
    for row in null_roots:
        record = row.as_dict()
        record.update(
            {
                "target_cl_mM": 36.50,
                "target_cl_sem_mM": 1.60,
                "target_ph": 6.89,
                "target_ph_sem": 0.02,
                "cl_target_z": (row.cl_i_mM - 36.50) / 1.60,
                "ph_target_z": (row.pH_i - 6.89) / 0.02,
                "target_join_stage": "Stage A reveal after marker hash",
            }
        )
        null_rows.append(record)
    _write_csv(
        result_dir / "stage_a_null_roots.csv",
        null_rows,
        fieldnames=list(null_rows[0]) if null_rows else ["model_id", "failure_reason"],
    )

    model_scores = []
    rejected = []
    attempts_by_key = {(row["model_id"], row["gauge_id"]): row for row in calibration_attempts}
    roots_by_key: dict[tuple[str, str], list[Mapping[str, object]]] = {}
    for row in root_rows:
        roots_by_key.setdefault((str(row["model_id"]), str(row["gauge_id"])), []).append(row)
    for fit in fits:
        key = (fit.model_id, fit.gauge_id)
        attempt = attempts_by_key.get(key)
        roots = roots_by_key.get(key, [])
        wt_pass = any(bool(row["wt_physiology_gate_pass"]) for row in roots)
        stage_status = (
            "transporter_rejected"
            if not fit.transporter_gate_pass
            else "WT_fixed_chassis_pass" if wt_pass else "WT_fixed_chassis_fail"
        )
        model_scores.append(
            {
                "model_id": fit.model_id,
                "gauge_id": fit.gauge_id,
                "transporter_gate_pass": fit.transporter_gate_pass,
                "wt_capacity_calibration_valid": (attempt or {}).get("calibration_valid", False),
                "wt_root_gate_pass": wt_pass,
                "stage_a_status": stage_status,
                "strict_secretion_used_in_fit": False,
                "comparison_type": "conditional compatibility; not AIC",
            }
        )
        if stage_status != "WT_fixed_chassis_pass":
            rejected.append(
                {
                    "model_id": fit.model_id,
                    "gauge_id": fit.gauge_id,
                    "rejection_stage": stage_status,
                    "reason": fit.failure_reason if not fit.transporter_gate_pass else (
                        (attempt or {}).get("failure_reason", "no rest-compatible root")
                        + ("; exact roots fail WT physiology" if roots else "; no bounded exact root")
                    ),
                    "heldout_secretion_used": False,
                }
            )
    _write_csv(result_dir / "model_scores.csv", model_scores)
    _write_csv(result_dir / "rejected_models.csv", rejected)

    summary = {
        "freeze": boundary.as_dict(),
        "transporter_fit_rows": len(fits),
        "transporter_survivor_rows": len(survivors),
        "valid_wt_capacity_rows": sum(row["calibration_valid"] for row in calibration_attempts),
        "wt_fixed_chassis_pass_rows": sum(row["wt_root_gate_pass"] for row in model_scores),
        "bounded_common_null_root_count": len(null_roots),
        "stage_a_all_source_supported_families_fail_fixed_chassis": not any(
            row["wt_root_gate_pass"] for row in model_scores
        ),
        "strict_secretion_used_in_any_fit": False,
    }
    (result_dir / "stage_a_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
