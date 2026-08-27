"""Generate the frozen Round-5 chassis sensitivity ledger.

Run from the repository root with ``PYTHONPATH=src python -m
state_resolved_ae4.run_chassis_reconstruction``.  The script reads only the
Stage-A gate/root artifacts and the two released Stage-B ionic targets encoded
in the frozen evidence.  It never reads the held-out-target ledger.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.optimize import least_squares, root

from ae4_mechanism_reconstruction.chassis import AE4_KNOCKOUT, FixedChassis, ZeroAE4

from . import chassis as split
from . import solvers as fixed_solvers
from .adapter import StateResolvedAE4Adapter
from .cycles import CANDIDATE_DEFINITIONS, CycleParameters
from .localization import evaluated_target_slices, released_stage_b_manifold


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "13_state_resolved_ae4"
OUTPUT = RESULTS / "chassis_reconstruction_scan.csv"

# Released only after Stage-A failure (E15-06 and E15-07_KO).
CL_TARGET_MM = 36.50
CL_TARGET_SEM_MM = 1.60
PH_TARGET = 6.89
PH_TARGET_SEM = 0.02


def _stage_a_gate_and_root() -> np.ndarray:
    summary = json.loads((RESULTS / "stage_a_summary.json").read_text())
    if not summary["stage_a_all_source_supported_families_fail_fixed_chassis"]:
        raise PermissionError("Round-5 Stage-B analysis requires frozen Stage-A failure")
    if summary["strict_secretion_used_in_any_fit"]:
        raise RuntimeError("Stage-A ledger reports forbidden secretion target use")
    with (RESULTS / "stage_a_null_roots.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    return np.asarray(
        tuple(
            float(row[name])
            for name in (
                "na_l_mM",
                "k_l_mM",
                "height_um",
                "na_i_mM",
                "k_i_mM",
                "cl_i_mM",
                "hco3_i_mM",
            )
        )
    )


def _solve_step(
    fp: float, fk: float, log_start: np.ndarray
) -> tuple[bool, np.ndarray, np.ndarray, dict[str, Any], float, str]:
    model = split.SplitCationChassis(
        ZeroAE4(), topology=split.CationTopology(fp, fk)
    )

    def residual(log_state: np.ndarray) -> np.ndarray:
        return split._split_root_residual(  # structural analysis entry point
            log_state, model, AE4_KNOCKOUT, 0.05
        )

    solution = root(
        residual,
        log_start,
        method="hybr",
        options={"xtol": 1e-10, "maxfev": 3000},
    )
    state = np.exp(solution.x)
    try:
        raw, diagnostics = model.evaluate(
            state, calcium=0.05, pka_activation=0.0, scenario=AE4_KNOCKOUT
        )
        error = float(np.max(np.abs(raw)))
    except (ValueError, FloatingPointError, OverflowError) as exc:
        return False, solution.x, state, {}, math.inf, str(exc)
    in_domain = bool(
        np.all(state >= split.SPLIT_STATE_LOWER)
        and np.all(state <= split.SPLIT_STATE_UPPER)
    )
    accepted = bool(error <= 1e-8 and in_domain)
    reason = "" if accepted else (
        f"continuation did not close: max_abs_raw_rhs={error:.9g}; "
        f"solver_message={solution.message}; in_declared_domain={in_domain}"
    )
    return accepted, solution.x, state, dict(diagnostics), error, reason


def _trace(
    pairs: Iterable[tuple[float, float]], start_state: np.ndarray
) -> tuple[
    dict[tuple[float, float], tuple[np.ndarray, dict[str, Any], float]],
    tuple[float, float, str] | None,
]:
    log_state = np.log(start_state)
    saved: dict[tuple[float, float], tuple[np.ndarray, dict[str, Any], float]] = {}
    for fp, fk in pairs:
        accepted, candidate, state, diagnostics, error, reason = _solve_step(
            fp, fk, log_state
        )
        if not accepted:
            return saved, (fp, fk, reason)
        log_state = candidate
        saved[(fp, fk)] = (state, diagnostics, error)
    return saved, None


def _line(endpoint: float, kind: str) -> list[tuple[float, float]]:
    values = [float(round(value, 3)) for value in np.arange(0.0, endpoint + 0.0025, 0.005)]
    if values[-1] != endpoint:
        values.append(endpoint)
    if kind == "pump":
        return [(value, 0.0) for value in values]
    if kind == "k":
        return [(0.0, value) for value in values]
    if kind == "diagonal":
        return [(value, value) for value in values]
    raise KeyError(kind)


def _cross(
    start: tuple[float, float], endpoint: tuple[float, float]
) -> list[tuple[float, float]]:
    fp0, fk0 = start
    fp1, fk1 = endpoint
    pairs: list[tuple[float, float]] = []
    if fp1 != fp0:
        direction = 1.0 if fp1 > fp0 else -1.0
        values = np.arange(fp0 + direction * 0.005, fp1 + direction * 0.0025, direction * 0.005)
        pairs.extend((float(round(value, 3)), fk0) for value in values)
        if not pairs or pairs[-1][0] != fp1:
            pairs.append((fp1, fk0))
    if fk1 != fk0:
        direction = 1.0 if fk1 > fk0 else -1.0
        values = np.arange(fk0 + direction * 0.005, fk1 + direction * 0.0025, direction * 0.005)
        pairs.extend((fp1, float(round(value, 3))) for value in values)
        if not pairs or pairs[-1][1] != fk1:
            pairs.append((fp1, fk1))
    return pairs


def _prediction_row(
    label: str,
    fp: float,
    fk: float,
    state: np.ndarray,
    diagnostics: dict[str, Any],
    error: float,
    path: str,
) -> dict[str, Any]:
    d_cl = float(state[5] - CL_TARGET_MM)
    d_ph = float(diagnostics["pH_i"] - PH_TARGET)
    return {
        "record_type": "resolved_null_root",
        "label": label,
        "scenario": "ae4_knockout",
        "topology_model_id": split.CationTopology(fp, fk).model_id,
        "apical_pump_fraction": fp,
        "apical_k_fraction": fk,
        "fraction_provenance": (
            "predeclared parotid/model sensitivity; uncalibrated for mouse submandibular acini"
        ),
        "continuation_path": path,
        "root_success": True,
        "failure_at_fraction": "",
        "failure_reason": "",
        "max_abs_raw_rhs": error,
        "na_l_mM": state[0],
        "k_l_mM": state[1],
        "height_um": state[2],
        "na_i_mM": state[3],
        "k_i_mM": state[4],
        "cl_i_mM": state[5],
        "hco3_i_mM": state[6],
        "pH_i": diagnostics["pH_i"],
        "cl_target_residual_mM": d_cl,
        "pH_target_residual": d_ph,
        "normalized_ionic_distance": math.hypot(
            d_cl / CL_TARGET_SEM_MM, d_ph / PH_TARGET_SEM
        ),
        "target_stage": "LOC-B",
        "conditional_slice": "re-solved model root; not a target-manifold slice",
        "raw_rhs_vector": ";".join(
            f"{value:.12g}" for value in diagnostics["raw_rhs"].values()
        ),
        "acid_base_full_rank": "",
        "acid_base_cl_ph_rank": "",
        "acid_base_flux_nullity": "",
        "valid_wt_capacity_available": False,
        "wt_gate_evaluable": False,
        "secretion_target_read": False,
        "interpretation": (
            "ionic sensitivity only; absolute flow/time uncertified and no WT-valid denominator"
        ),
    }


def _failure_row(
    label: str,
    requested: tuple[float, float],
    failure: tuple[float, float, str],
    path: str,
) -> dict[str, Any]:
    fp, fk = requested
    fail_fp, fail_fk, reason = failure
    return {
        "record_type": "continuation_failure",
        "label": label,
        "scenario": "ae4_knockout",
        "topology_model_id": split.CationTopology(fp, fk).model_id,
        "apical_pump_fraction": fp,
        "apical_k_fraction": fk,
        "fraction_provenance": (
            "predeclared parotid/model sensitivity; uncalibrated for mouse submandibular acini"
        ),
        "continuation_path": path,
        "root_success": False,
        "failure_at_fraction": f"{fail_fp};{fail_fk}",
        "failure_reason": reason,
        "target_stage": "LOC-B",
        "valid_wt_capacity_available": False,
        "wt_gate_evaluable": False,
        "secretion_target_read": False,
        "interpretation": (
            "failed declared continuation is retained; it is not a proof that no other branch exists"
        ),
    }


def _slice_rows(
    label: str,
    fp: float,
    fk: float,
    slices: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    model = split.SplitCationChassis(
        ZeroAE4(), topology=split.CationTopology(fp, fk)
    )
    for slice_name in ("model_anchored", "wt_reference_cations"):
        state = np.asarray(slices[slice_name].state)
        raw, _ = model.evaluate(
            state, calcium=0.05, pka_activation=0.0, scenario=AE4_KNOCKOUT
        )
        rows.append(
            {
                "record_type": "conditional_target_slice",
                "label": label,
                "scenario": "ae4_knockout",
                "topology_model_id": split.CationTopology(fp, fk).model_id,
                "apical_pump_fraction": fp,
                "apical_k_fraction": fk,
                "fraction_provenance": (
                    "predeclared sensitivity; slice free coordinates are declared assumptions"
                ),
                "root_success": bool(np.max(np.abs(raw)) <= 1e-8),
                "max_abs_raw_rhs": float(np.max(np.abs(raw))),
                "target_stage": "LOC-B",
                "conditional_slice": slice_name,
                "raw_rhs_vector": ";".join(f"{value:.12g}" for value in raw),
                "valid_wt_capacity_available": False,
                "wt_gate_evaluable": False,
                "secretion_target_read": False,
                "interpretation": (
                    "conditional Cl/pH-compatible state only; unmeasured Na/K/height/lumen coordinates are not data"
                ),
            }
        )
    return rows


def _wt_capacity_attempt_rows() -> list[dict[str, Any]]:
    """Repeat the WT-Cl-only capacity attempt on four nested M1 controls.

    SR2-123/slow is the smallest conditional transporter-assay survivor.  Its
    transporter parameters are frozen before Stage B.  One capacity and seven
    steady coordinates are attempted against seven balances plus WT chloride;
    WT pH remains an acceptance gate and is not included in the objective.
    """

    with (RESULTS / "frozen_survivor_parameters.csv").open(newline="") as handle:
        survivors = list(csv.DictReader(handle))
    frozen = next(
        row
        for row in survivors
        if row["model_id"] == "SR2_SHARED_123" and row["gauge_id"] == "slow_common"
    )
    definition = CANDIDATE_DEFINITIONS[frozen["model_id"]]
    gauge = float(frozen["common_barrier_multiplier"])
    overrides = {
        edge.edge_id: gauge
        for edge in definition.edges
        if edge.edge_id.endswith("common_cl_flip")
    }

    def cycle_parameters(capacity: float) -> CycleParameters:
        return CycleParameters(
            capacity=capacity,
            na_attempt_scale=1.0,
            k_attempt_scale=float(frozen["k_attempt_scale"]),
            ec50_na_mM=float(frozen["na_loaded_energy_lump_mM"]),
            ec50_k_mM=float(frozen["k_loaded_energy_lump_mM"]),
            pka_barrier_fold=float(frozen["pka_barrier_fold"]),
            edge_barrier_overrides=overrides,
        )

    lower = np.concatenate((fixed_solvers.REDUCED_LOWER, np.asarray((-12.0,))))
    upper = np.concatenate((fixed_solvers.REDUCED_UPPER, np.asarray((16.0,))))
    rows = []
    nested = (
        ("M0_control", 0.0, 0.0),
        ("M1_pump_20pct", 0.2, 0.0),
        ("M1_k_20pct", 0.0, 0.2),
        ("M1_coupled_20pct", 0.2, 0.2),
        ("M1_coupled_40pct", 0.4, 0.4),
    )
    for label, fp, fk in nested:
        topology = split.CationTopology(fp, fk)
        reference = split.SplitCationChassis(
            StateResolvedAE4Adapter(definition, cycle_parameters(1.0)),
            topology=topology,
        )

        def residual(coordinates: np.ndarray) -> np.ndarray:
            state = fixed_solvers._reduced_to_state(coordinates[:7], reference)
            if not np.all(np.isfinite(state)) or not 2.0 < state[2] < 200.0:
                return np.full(8, 1e6)
            capacity = math.exp(float(coordinates[7]))
            model = split.SplitCationChassis(
                StateResolvedAE4Adapter(definition, cycle_parameters(capacity)),
                topology=topology,
            )
            try:
                raw, diagnostics = model.evaluate(state, calcium=0.05)
            except (ValueError, FloatingPointError, OverflowError):
                return np.full(8, 1e6)
            if float(diagnostics["h_i"]) <= 0.0:
                return np.full(8, 1e6)
            return np.concatenate(
                (
                    np.asarray(raw) / fixed_solvers.RAW_SCALES,
                    np.asarray(
                        (
                            (state[5] - fixed_solvers.WT_CL_MM)
                            / fixed_solvers.WT_CL_SEM_MM,
                        )
                    ),
                )
            )

        starts = [
            np.concatenate(
                (
                    fixed_solvers._state_to_reduced(
                        fixed_solvers.BASELINE_STATE, reference
                    ),
                    np.asarray((log_capacity,)),
                )
            )
            for log_capacity in (-8.0, -3.0, 2.0, 7.0, 12.0)
        ]
        fits = [
            least_squares(
                residual,
                np.clip(start, lower, upper),
                bounds=(lower, upper),
                x_scale="jac",
                max_nfev=1200,
                xtol=1e-11,
                ftol=1e-11,
                gtol=1e-11,
            )
            for start in starts
        ]
        best = min(fits, key=lambda fit: fit.cost)
        state = fixed_solvers._reduced_to_state(best.x[:7], reference)
        capacity = math.exp(float(best.x[7]))
        model = split.SplitCationChassis(
            StateResolvedAE4Adapter(definition, cycle_parameters(capacity)),
            topology=topology,
        )
        raw, diagnostics = model.evaluate(state, calcium=0.05)
        raw_error = float(np.max(np.abs(raw)))
        cl_z = float(
            (state[5] - fixed_solvers.WT_CL_MM) / fixed_solvers.WT_CL_SEM_MM
        )
        ph_z = float(
            (diagnostics["pH_i"] - fixed_solvers.WT_PH)
            / fixed_solvers.WT_PH_SEM
        )
        calibration_valid = bool(raw_error <= 1e-8 and abs(cl_z) <= 2.0)
        gate_pass = bool(
            calibration_valid
            and abs(ph_z) <= 2.0
            and 5.0 <= state[3] <= 80.0
            and 70.0 <= state[4] <= 170.0
            and 0.5 <= diagnostics["cell_volume_pL"] <= 5.0
        )
        rows.append(
            {
                "record_type": "wt_capacity_attempt",
                "label": label,
                "scenario": "wt",
                "topology_model_id": topology.model_id,
                "apical_pump_fraction": fp,
                "apical_k_fraction": fk,
                "fraction_provenance": (
                    "predeclared sensitivity; not a measured submandibular fraction"
                ),
                "continuation_path": "five-start bounded joint WT-Cl solve",
                "root_success": bool(raw_error <= 1e-8),
                "failure_at_fraction": "",
                "failure_reason": (
                    ""
                    if calibration_valid
                    else (
                        f"raw closure {raw_error:.9g} exceeds 1e-8; "
                        f"derived height={state[2]:.9g} um"
                    )
                ),
                "max_abs_raw_rhs": raw_error,
                "na_l_mM": state[0],
                "k_l_mM": state[1],
                "height_um": state[2],
                "na_i_mM": state[3],
                "k_i_mM": state[4],
                "cl_i_mM": state[5],
                "hco3_i_mM": state[6],
                "pH_i": diagnostics["pH_i"],
                "cl_target_residual_mM": "",
                "pH_target_residual": "",
                "normalized_ionic_distance": "",
                "wt_cl_residual_mM": state[5] - fixed_solvers.WT_CL_MM,
                "wt_cl_z": cl_z,
                "wt_pH_residual": diagnostics["pH_i"] - fixed_solvers.WT_PH,
                "wt_pH_z": ph_z,
                "optimizer_cost": best.cost,
                "attempted_capacity": capacity,
                "calibration_valid": calibration_valid,
                "wt_gate2_pass": gate_pass,
                "derived_height_boundary_hit": bool(state[2] >= 199.99),
                "ae2_test_status": "N/A_upstream_WT_calibration_failed",
                "target_stage": "CAL-WT",
                "conditional_slice": "not_applicable",
                "raw_rhs_vector": ";".join(f"{value:.12g}" for value in raw),
                "valid_wt_capacity_available": calibration_valid,
                "wt_gate_evaluable": calibration_valid,
                "secretion_target_read": False,
                "interpretation": (
                    f"SR2_SHARED_123/slow; attempted capacity={capacity:.12g}; "
                    f"WT Gate-2 pass={gate_pass}; pH was a gate, not a fit target; "
                    "AE2 dynamic test N/A because upstream WT calibration failed"
                ),
            }
        )
    return rows


def generate() -> Path:
    base = _stage_a_gate_and_root()
    rows: list[dict[str, Any]] = _wt_capacity_attempt_rows()

    pump, pump_failure = _trace(_line(0.4, "pump"), base)
    potassium, k_failure = _trace(_line(0.4, "k"), base)
    diagonal, diagonal_failure = _trace(_line(0.4, "diagonal"), base)
    if diagonal_failure is not None:
        raise RuntimeError(f"coupled continuation unexpectedly failed: {diagonal_failure}")

    configurations: dict[
        str, tuple[float, float, np.ndarray, dict[str, Any], float, str]
    ] = {}
    # M0 is present on every path; use the pump path copy.
    state, diagnostics, error = pump[(0.0, 0.0)]
    configurations["M0_control"] = (0.0, 0.0, state, diagnostics, error, "pump-line")
    for label, store, pair, path in (
        ("M1_pump_20pct", pump, (0.2, 0.0), "pump-line"),
        ("M1_k_20pct", potassium, (0.0, 0.2), "k-line"),
        ("M1_coupled_20pct", diagonal, (0.2, 0.2), "diagonal"),
        ("M1_coupled_40pct", diagonal, (0.4, 0.4), "diagonal"),
    ):
        state, diagnostics, error = store[pair]
        configurations[label] = (*pair, state, diagnostics, error, path)

    diagonal_20_state = diagonal[(0.2, 0.2)][0]
    for label, endpoint in (
        ("M1_pump40_k20", (0.4, 0.2)),
        ("M1_pump20_k40", (0.2, 0.4)),
    ):
        traced, failure = _trace(_cross((0.2, 0.2), endpoint), diagonal_20_state)
        if failure is not None:
            rows.append(_failure_row(label, endpoint, failure, "cross-from-coupled-20"))
            continue
        state, diagnostics, error = traced[endpoint]
        configurations[label] = (
            *endpoint,
            state,
            diagnostics,
            error,
            "cross-from-coupled-20",
        )

    # The historical area ratio is an initialization/sensitivity, not a
    # localization measurement.  It is reached from the exact 40/40 root.
    state_40 = diagonal[(0.4, 0.4)][0]
    area_trace, area_failure = _trace(((0.398, 0.398),), state_40)
    if area_failure is None:
        state, diagnostics, error = area_trace[(0.398, 0.398)]
        configurations["M1_coupled_area_sensitivity"] = (
            0.398,
            0.398,
            state,
            diagnostics,
            error,
            "backstep-from-coupled-40",
        )
    else:
        rows.append(
            _failure_row(
                "M1_coupled_area_sensitivity",
                (0.398, 0.398),
                area_failure,
                "backstep-from-coupled-40",
            )
        )

    if pump_failure is not None:
        rows.append(
            _failure_row("M1_pump_40pct", (0.4, 0.0), pump_failure, "pump-line")
        )
    if k_failure is not None:
        rows.append(
            _failure_row("M1_k_40pct", (0.0, 0.4), k_failure, "k-line")
        )

    base_chassis = FixedChassis(ZeroAE4())
    manifold = released_stage_b_manifold(
        base_chassis,
        stage_a_failed=True,
        cl_i_mM=CL_TARGET_MM,
        pH_i=PH_TARGET,
    )
    target_slices = dict(
        evaluated_target_slices(
            base_chassis, base, manifold, stage_a_failed=True
        )
    )
    for label, (fp, fk, state, diagnostics, error, path) in configurations.items():
        rows.append(_prediction_row(label, fp, fk, state, diagnostics, error, path))
        rows.extend(_slice_rows(label, fp, fk, target_slices))

    audit = split.acid_base_identifiability_audit()
    rows.append(
        {
            "record_type": "acid_base_identifiability",
            "label": "M2_conserved_carbon_buffer_NHE1",
            "scenario": "structural",
            "topology_model_id": "M2_NOT_CALIBRATABLE_FROM_RELEASED_PAIR",
            "fraction_provenance": "not_applicable",
            "root_success": False,
            "target_stage": "LOC-B",
            "acid_base_full_rank": audit["full_reaction_signature_rank"],
            "acid_base_cl_ph_rank": audit["released_cl_ph_direct_rank"],
            "acid_base_flux_nullity": audit["reaction_flux_nullity_under_cl_ph"],
            "valid_wt_capacity_available": False,
            "wt_gate_evaluable": False,
            "secretion_target_read": False,
            "interpretation": audit["conclusion"],
        }
    )
    rows.append(
        {
            "record_type": "wt_gate_status",
            "label": "all_round5_sensitivities",
            "scenario": "wt",
            "topology_model_id": "M1/M2",
            "fraction_provenance": (
                "no gland-matched submandibular split fraction and Stage A supplied zero valid capacities"
            ),
            "root_success": False,
            "target_stage": "CAL-WT",
            "valid_wt_capacity_available": False,
            "wt_gate_evaluable": False,
            "secretion_target_read": False,
            "interpretation": (
                "no Round-5 row can be frozen as a WT-validated predictive model; null ionic scans remain localization sensitivities"
            ),
        }
    )

    fields = (
        "record_type",
        "label",
        "scenario",
        "topology_model_id",
        "apical_pump_fraction",
        "apical_k_fraction",
        "fraction_provenance",
        "continuation_path",
        "root_success",
        "failure_at_fraction",
        "failure_reason",
        "max_abs_raw_rhs",
        "na_l_mM",
        "k_l_mM",
        "height_um",
        "na_i_mM",
        "k_i_mM",
        "cl_i_mM",
        "hco3_i_mM",
        "pH_i",
        "cl_target_residual_mM",
        "pH_target_residual",
        "normalized_ionic_distance",
        "wt_cl_residual_mM",
        "wt_cl_z",
        "wt_pH_residual",
        "wt_pH_z",
        "optimizer_cost",
        "attempted_capacity",
        "calibration_valid",
        "wt_gate2_pass",
        "derived_height_boundary_hit",
        "ae2_test_status",
        "target_stage",
        "conditional_slice",
        "raw_rhs_vector",
        "acid_base_full_rank",
        "acid_base_cl_ph_rank",
        "acid_base_flux_nullity",
        "valid_wt_capacity_available",
        "wt_gate_evaluable",
        "secretion_target_read",
        "interpretation",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return OUTPUT


if __name__ == "__main__":
    print(generate())
