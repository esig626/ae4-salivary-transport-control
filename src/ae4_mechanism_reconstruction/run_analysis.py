"""Generate all Task-12 AE4 comparison tables, trajectories, and figures."""

from __future__ import annotations

import csv
from dataclasses import fields, replace
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .chassis import (
    AE4_KNOCKOUT,
    BASELINE_STATE,
    STATE_ORDER,
    FixedChassis,
)
from .comparison import (
    BALANCE_NAMES,
    HELDOUT_TARGET_NAMES,
    HISTORICAL_WT_AE4_CONTRIBUTION,
    PATHWAY_SIGNATURES,
    ROUND2_BATH_HCO3_MILLIMOLAR,
    ROUND2_EQUILIBRIUM_CONSTANTS,
    ROUND2_NA_CAPACITY_FRACTIONS,
    ROUND2_SATURATION_STRENGTHS,
    CandidateComparison,
    RestingEquilibrium,
    StimulusComparison,
    PKA_ASSAY_FOLDS,
    PKA_IMMEDIATE_MODEL_ONSET,
    REST_MULTISTART_COUNT,
    REST_MULTISTART_SEED,
    REST_VOLUME_RELATIVE_TOLERANCE,
    calibrate_candidate,
    candidate_chassis,
    compare_calibrated_candidate,
    independently_constrained_parameters,
    pathway_localization,
    primary_comparisons,
    replay_trajectory,
    round2_grid,
    stimulus_activation_comparisons,
    transported_k_fraction_bound,
)


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/12_ae4_mechanism_reconstruction"
FIGURES = ROOT / "figures/12_ae4_mechanism_reconstruction"
HELDOUT_AE4_KO_FLOW_RATIO = 0.65

GENOTYPE_REST_TARGETS: Mapping[str, Mapping[str, Any]] = {
    "COMMON_AE4_KO_REST": {
        "target_label": "primary_AE4_KO_rest",
        "cl_mM": 36.50,
        "cl_SEM_mM": 1.60,
        "pH": 6.89,
        "pH_SEM": 0.02,
    },
    "C1_AE2_KO_REST": {
        "target_label": "primary_AE2_KO_rest_evaluated_on_C1_lineage_chassis",
        "cl_mM": 54.50,
        "cl_SEM_mM": 1.80,
        "pH": 6.95,
        "pH_SEM": 0.05,
        "own_control_cl_mM": 53.40,
        "own_control_cl_SEM_mM": 1.80,
        "own_control_pH": 6.87,
        "own_control_pH_SEM": 0.01,
    },
}


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_jsonable(item) for item in value.tolist()]
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, (bool, np.bool_)) or value is None:
        return bool(value) if value is not None else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, str):
        return value
    return str(value)


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ""
                    if value is None or (isinstance(value, float) and not math.isfinite(value))
                    else json.dumps(value, sort_keys=True)
                    if isinstance(value, (dict, list, tuple))
                    else value
                    for key, value in row.items()
                }
            )


def _comparison_columns() -> list[str]:
    return [field.name for field in fields(CandidateComparison)]


def _dataclass_columns(record_type: type[Any]) -> list[str]:
    return [field.name for field in fields(record_type)]


def _union_columns(rows: Iterable[Mapping[str, Any]]) -> list[str]:
    columns: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                columns.append(key)
                seen.add(key)
    return columns


def _genotype_rest_rows(
    rest_dicts: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Join direct genotype-specific resting targets after root prediction."""

    output: list[dict[str, Any]] = []
    for row in rest_dicts:
        mechanism_id = str(row["mechanism_id"])
        target = GENOTYPE_REST_TARGETS.get(mechanism_id)
        if target is None and mechanism_id.endswith("_AE2_KO_REST"):
            target = {
                **GENOTYPE_REST_TARGETS["C1_AE2_KO_REST"],
                "target_label": "primary_AE2_KO_rest_candidate_specific_prediction",
            }
        if target is None:
            continue
        cl_z = (float(row["cl_i_mM"]) - float(target["cl_mM"])) / float(
            target["cl_SEM_mM"]
        )
        ph_z = (float(row["pH_i"]) - float(target["pH"])) / float(
            target["pH_SEM"]
        )
        genotype_row = {
            **row,
            "target_label": target["target_label"],
            "genotype_cl_target_mM": target["cl_mM"],
            "genotype_cl_target_SEM_mM": target["cl_SEM_mM"],
            "genotype_cl_target_z": cl_z,
            "genotype_pH_target": target["pH"],
            "genotype_pH_target_SEM": target["pH_SEM"],
            "genotype_pH_target_z": ph_z,
            "genotype_rest_gate_pass": bool(
                abs(cl_z) <= 2.0
                and abs(ph_z) <= 2.0
                and float(row["max_abs_raw_rhs"]) <= 1e-8
            ),
        }
        if "own_control_cl_mM" in target:
            genotype_row.update(
                {
                    "own_control_cl_target_mM": target["own_control_cl_mM"],
                    "own_control_cl_target_SEM_mM": target[
                        "own_control_cl_SEM_mM"
                    ],
                    "own_control_cl_target_z": (
                        float(row["cl_i_mM"])
                        - float(target["own_control_cl_mM"])
                    )
                    / float(target["own_control_cl_SEM_mM"]),
                    "own_control_pH_target": target["own_control_pH"],
                    "own_control_pH_target_SEM": target["own_control_pH_SEM"],
                    "own_control_pH_target_z": (
                        float(row["pH_i"]) - float(target["own_control_pH"])
                    )
                    / float(target["own_control_pH_SEM"]),
                }
            )
        output.append(genotype_row)
    return output


def _downsample_mask(times: np.ndarray) -> np.ndarray:
    # The dense 0.1-unit grid remains in memory for integration.  Export only
    # a compact set of protocol landmarks; _trajectory_rows additionally
    # retains the first point on both sides of every calcium/PKA input change.
    landmarks = np.asarray((0.0, 50.0, 100.0, 150.0, 200.0), dtype=float)
    protocol_landmark = np.any(
        np.isclose(times[:, None], landmarks[None, :], atol=2e-9), axis=1
    )
    right_limit = np.isclose(times, 100.0 + 1e-7, atol=2e-10)
    return protocol_landmark | right_limit


def _trajectory_rows(
    all_runs: Mapping[str, Mapping[str, Any]],
    comparison_round: Mapping[str, int],
) -> Iterable[dict[str, Any]]:
    for mechanism_id, scenarios in all_runs.items():
        for scenario_name, result in scenarios.items():
            keep = _downsample_mask(result.time)
            diagnostics = np.asarray(result.diagnostics, dtype=object)
            calcium = np.asarray(
                [float(item["calcium"]) for item in diagnostics], dtype=float
            )
            pka = np.asarray(
                [float(item.get("pka_activation", 0.0)) for item in diagnostics],
                dtype=float,
            )
            # Retain the first saved point on each side of any input change.
            for values in (calcium, pka):
                changed = np.flatnonzero(~np.isclose(np.diff(values), 0.0, atol=1e-14))
                keep[changed] = True
                keep[changed + 1] = True
            for time, state, diagnostic in zip(
                result.time[keep], result.state[keep], diagnostics[keep], strict=True
            ):
                raw_max = max(abs(float(value)) for value in diagnostic["raw_rhs"].values())
                check_max = max(abs(float(value)) for value in diagnostic["checks"].values())
                row: dict[str, Any] = {
                    "mechanism_id": mechanism_id,
                    "round": comparison_round.get(mechanism_id, 3),
                    "scenario": scenario_name,
                    "time": float(time),
                    "calcium": float(diagnostic["calcium"]),
                    "pka_activation": float(diagnostic.get("pka_activation", 0.0)),
                    "q_total": float(diagnostic["q_total"]),
                    "pH_i": float(diagnostic["pH_i"]),
                    "ae4_cycle_flux": float(diagnostic["ae4_cycle_flux"]),
                    "max_abs_raw_rhs": raw_max,
                    "max_abs_conservation_check": check_max,
                    "solver_method": result.method,
                    "non_ae4_chassis_fingerprint": result.parameter_fingerprint,
                }
                row.update(
                    {
                        name: float(value)
                        for name, value in zip(STATE_ORDER, state, strict=True)
                    }
                )
                yield row


def _round2_best_bdf(
    grid: list[dict[str, Any]],
    rows: list[CandidateComparison],
    shared_bdf_ko: Any,
) -> tuple[CandidateComparison, Mapping[str, Any], dict[str, Any]]:
    finite_rows = [
        row
        for row in rows
        if row.wt_baseline_admissible
        and math.isfinite(row.ae4_ko_endpoint_flow_ratio)
    ]
    best_screen = min(finite_rows, key=lambda row: row.ae4_ko_endpoint_flow_ratio)
    definition = next(item for item in grid if item["variant_id"] == best_screen.mechanism_id)
    parameters = independently_constrained_parameters(
        na_weight=float(definition["na_capacity_fraction"]),
        k_weight=float(definition["k_capacity_fraction"]),
        equilibrium_constant=float(definition["equilibrium_constant"]),
        saturation_strength=float(definition["saturation_strength"]),
    )
    calibration = calibrate_candidate(
        str(definition["candidate_id"]),
        parameters=parameters,
        external_hco3_mM=float(definition["bath_hco3_mM"]),
    )
    chassis = candidate_chassis(calibration)
    common_ko = replay_trajectory(chassis, shared_bdf_ko, AE4_KNOCKOUT)
    row, runs = compare_calibrated_candidate(
        calibration,
        round_number=2,
        mechanism_id=best_screen.mechanism_id,
        solver_method="BDF",
        precomputed_simulations={"ae4_knockout": common_ko},
    )
    crosscheck = {
        "mechanism_id": row.mechanism_id,
        "radau_endpoint_ratio": best_screen.ae4_ko_endpoint_flow_ratio,
        "bdf_endpoint_ratio": row.ae4_ko_endpoint_flow_ratio,
        "absolute_endpoint_ratio_difference": abs(
            row.ae4_ko_endpoint_flow_ratio - best_screen.ae4_ko_endpoint_flow_ratio
        ),
        "radau_integrated_ratio": best_screen.ae4_ko_integrated_flow_ratio,
        "bdf_integrated_ratio": row.ae4_ko_integrated_flow_ratio,
        "absolute_integrated_ratio_difference": abs(
            row.ae4_ko_integrated_flow_ratio - best_screen.ae4_ko_integrated_flow_ratio
        ),
        "interpretation": "diagnostic pure-Na/energetic-bias boundary; not source-admissible",
    }
    return row, runs, crosscheck


def _c1_radau_crosscheck(bdf_row: CandidateComparison) -> dict[str, Any]:
    calibration = calibrate_candidate("C1")
    radau_row, _ = compare_calibrated_candidate(
        calibration, round_number=1, solver_method="Radau"
    )
    return {
        "mechanism_id": "C1",
        "bdf_endpoint_ratio": bdf_row.ae4_ko_endpoint_flow_ratio,
        "radau_endpoint_ratio": radau_row.ae4_ko_endpoint_flow_ratio,
        "absolute_endpoint_ratio_difference": abs(
            bdf_row.ae4_ko_endpoint_flow_ratio - radau_row.ae4_ko_endpoint_flow_ratio
        ),
        "bdf_integrated_ratio": bdf_row.ae4_ko_integrated_flow_ratio,
        "radau_integrated_ratio": radau_row.ae4_ko_integrated_flow_ratio,
        "absolute_integrated_ratio_difference": abs(
            bdf_row.ae4_ko_integrated_flow_ratio - radau_row.ae4_ko_integrated_flow_ratio
        ),
        "right_limit_and_poststep_grid": "100+1e-7 then 0.1-unit spacing",
    }


def _figure_value(record: CandidateComparison | Mapping[str, Any], name: str) -> Any:
    return record[name] if isinstance(record, Mapping) else getattr(record, name)


def _figure_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() == "true"


def _write_closure_prediction_tradeoff(
    primary: Iterable[CandidateComparison | Mapping[str, Any]],
) -> None:
    """Render the closure panel from live rows or the frozen comparison CSV."""

    rows = list(primary)
    plt.rcParams.update(
        {"font.size": 9, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    default_cluster: list[tuple[float, float]] = []
    label_offsets = {"C1": (3, 2), "C2": (-18, 10), "C3": (6, -13)}
    for row in rows:
        mechanism_id = str(_figure_value(row, "mechanism_id"))
        x = max(float(_figure_value(row, "calibration_max_abs_error")), 1e-16)
        y = float(_figure_value(row, "ae4_ko_endpoint_flow_ratio"))
        admissible = _figure_bool(_figure_value(row, "wt_baseline_admissible"))
        ax.scatter(
            x,
            y,
            marker="o" if admissible else "x",
            color="#1f77b4" if admissible else "#757575",
        )
        if mechanism_id in label_offsets:
            ax.annotate(
                mechanism_id,
                (x, y),
                xytext=label_offsets[mechanism_id],
                textcoords="offset points",
                fontsize=7,
            )
        else:
            default_cluster.append((x, y))
    if default_cluster:
        cluster_x = float(
            np.exp(np.mean(np.log([item[0] for item in default_cluster])))
        )
        cluster_y = float(np.median([item[1] for item in default_cluster]))
        ax.annotate(
            "C4-C8 defaults\n(overlapping)",
            (cluster_x, cluster_y),
            xytext=(8, -16),
            textcoords="offset points",
            fontsize=7,
        )
    ax.axvline(1e-8, color="#ef6c00", linestyle="--", label="WT vector tolerance")
    ax.axhline(
        HELDOUT_AE4_KO_FLOW_RATIO,
        color="#c62828",
        linestyle="--",
        label="held-out ratio",
    )
    ax.set_xscale("log")
    ax.set_xlabel("WT AE4 contribution max residual (chassis units)")
    ax.set_ylabel("t=200 AE4-KO / WT flow ratio")
    ax.set_title("No source-supported mixed-cation row crosses the closure gate")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "ae4_closure_prediction_tradeoff.png", dpi=180)
    plt.close(fig)


def render_closure_prediction_tradeoff_from_frozen_results() -> None:
    """Re-render the static closure panel without repeating model simulations."""

    FIGURES.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "model_comparison.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        primary = [
            row for row in csv.DictReader(handle) if int(row["round"]) == 1
        ]
    _write_closure_prediction_tradeoff(primary)


def _write_figures(
    primary: list[CandidateComparison],
    closed_round2: list[CandidateComparison],
    round2_grid_rows: list[dict[str, Any]],
    stimulus_rows: list[StimulusComparison],
    capacity_profile: list[dict[str, Any]],
) -> None:
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    # Round 2 has its own closure/K-fraction diagnostics.  Mixing its long
    # variant IDs into this panel obscures the primary C1--C8 comparison.
    selected = primary
    positions = np.arange(len(selected))
    endpoint = np.array([row.ae4_ko_endpoint_flow_ratio for row in selected])
    integrated = np.array([row.ae4_ko_integrated_flow_ratio for row in selected])
    colors = ["#1f77b4" if row.wt_baseline_admissible else "#9e9e9e" for row in selected]
    ax.scatter(positions - 0.12, endpoint, marker="o", c=colors, label="t=200 flow")
    ax.scatter(positions + 0.12, integrated, marker="s", c=colors, label="post-step integral")
    ax.axhline(HELDOUT_AE4_KO_FLOW_RATIO, color="#c62828", linestyle="--", label="held-out primary ratio 0.65")
    ax.axhline(1.0, color="black", linewidth=0.8)
    labels = [row.mechanism_id for row in selected]
    ax.set_xticks(positions, labels, rotation=75, ha="right")
    ax.set_ylabel("AE4-KO / WT flow ratio")
    ax.set_ylim(0.60, max(1.06, float(np.nanmax(integrated)) + 0.01))
    ax.set_title("Held-out phenotype predictions after WT-only calibration")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURES / "ae4_candidate_predictions.png", dpi=180)
    plt.close(fig)

    _write_closure_prediction_tradeoff(primary)

    parallel = [
        row
        for row in round2_grid_rows
        if row["candidate_id"] == "C3"
        and row["grid_axis"] == "parallel_Na_K_capacity"
    ]
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    ax.plot(
        [row["k_capacity_fraction"] for row in parallel],
        [row["calibration_max_abs_error"] for row in parallel],
        marker="o",
        color="#6a1b9a",
    )
    ax.axhline(1e-8, color="#ef6c00", linestyle="--", label="closure tolerance")
    ax.set_yscale("log")
    ax.set_xlabel("K branch capacity fraction")
    ax.set_ylabel("WT AE4-vector max residual")
    ax.set_title("Fixed historical closure drives the K branch to zero")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "ae4_k_fraction_closure.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    onset_rows = [
        row
        for row in stimulus_rows
        if row.simulation_status == "pass"
        and math.isclose(row.pka_activation_onset, PKA_IMMEDIATE_MODEL_ONSET)
        and row.pka_stimulated_fold <= max(PKA_ASSAY_FOLDS)
    ]
    families = sorted({row.family for row in onset_rows})
    palette = plt.cm.tab10(np.linspace(0.0, 1.0, max(1, len(families))))
    for color, family in zip(palette, families, strict=True):
        family_rows = sorted(
            [row for row in onset_rows if row.family == family],
            key=lambda row: (row.pka_stimulated_fold, row.activity_scale),
        )
        # Capacity-interval mechanisms have five values per fold.  Plot their
        # full envelope; singleton families reduce to a line.
        folds = sorted({row.pka_stimulated_fold for row in family_rows})
        low = []
        high = []
        middle = []
        for fold in folds:
            ratios = sorted(
                row.ae4_ko_endpoint_flow_ratio
                for row in family_rows
                if row.pka_stimulated_fold == fold
            )
            low.append(ratios[0])
            high.append(ratios[-1])
            middle.append(ratios[len(ratios) // 2])
        ax.plot(folds, middle, marker="o", color=color, label=family)
        if any(not math.isclose(a, b) for a, b in zip(low, high, strict=True)):
            ax.fill_between(folds, low, high, color=color, alpha=0.15)
    ax.axhline(HELDOUT_AE4_KO_FLOW_RATIO, color="#c62828", linestyle="--", label="held-out 0.65")
    ax.axhline(1.0, color="black", linewidth=0.8)
    ax.set_xlabel("Stimulated AE4 activity fold")
    ax.set_ylabel("t=200 AE4-KO / WT flow ratio")
    ax.set_title("KO-blind PKA-fold and rest-capacity envelopes")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURES / "ae4_pka_activation_envelope.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    profile_213 = [
        row
        for row in capacity_profile
        if row.get("core_candidate_id") in ("C6_213", "C7_213", "C8_213_SAT")
        and row.get("profile_kind") == "posthoc_KO_blind_finite_scan"
    ]
    for core_id in ("C6_213", "C7_213", "C8_213_SAT"):
        rows = sorted(
            [row for row in profile_213 if row["core_candidate_id"] == core_id],
            key=lambda row: float(row["activity_scale"]),
        )
        ax.plot(
            [float(row["cl_target_z"]) for row in rows],
            [float(row["pH_target_z"]) for row in rows],
            marker=".",
            linewidth=1.0,
            label=core_id,
        )
    ax.axvspan(-2.0, 2.0, color="#66bb6a", alpha=0.10)
    ax.axhspan(-2.0, 2.0, color="#66bb6a", alpha=0.10)
    ax.axvline(-2.0, color="#388e3c", linestyle="--", linewidth=0.8)
    ax.axvline(2.0, color="#388e3c", linestyle="--", linewidth=0.8)
    ax.axhline(-2.0, color="#388e3c", linestyle="--", linewidth=0.8)
    ax.axhline(2.0, color="#388e3c", linestyle="--", linewidth=0.8)
    ax.set_xlabel("WT resting Cl residual (SEM)")
    ax.set_ylabel("WT resting pH residual (SEM)")
    ax.set_title("2:1:3 capacity continuation reaches the resting-data box")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "ae4_rest_capacity_tradeoff.png", dpi=180)
    plt.close(fig)


def run() -> dict[str, Any]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    primary_rows, primary_runs = primary_comparisons()
    shared_bdf_ko = primary_runs["C1"]["ae4_knockout"]
    grid, round2_rows, round2_runs = round2_grid(
        shared_ae4_knockout=shared_bdf_ko
    )

    best_bdf_row, best_bdf_runs, best_crosscheck = _round2_best_bdf(
        grid, round2_rows, shared_bdf_ko
    )
    for index, row in enumerate(round2_rows):
        if row.mechanism_id == best_bdf_row.mechanism_id:
            round2_rows[index] = best_bdf_row
            break
    round2_runs[best_bdf_row.mechanism_id] = best_bdf_runs
    for record in grid:
        if record["variant_id"] == best_bdf_row.mechanism_id:
            record.update(
                {
                    "simulation_status": best_bdf_row.simulation_status,
                    "ae4_ko_endpoint_flow_ratio": best_bdf_row.ae4_ko_endpoint_flow_ratio,
                    "ae4_ko_integrated_flow_ratio": best_bdf_row.ae4_ko_integrated_flow_ratio,
                    "wt_endpoint_flow": best_bdf_row.wt_endpoint_flow,
                    "canonical_solver": "BDF",
                }
            )
            break

    c1_row = next(row for row in primary_rows if row.mechanism_id == "C1")
    c1_crosscheck = _c1_radau_crosscheck(c1_row)

    (
        frozen_stimulus_rows,
        stimulus_runs,
        rest_rows,
        thermo_profile,
        capacity_profile,
    ) = stimulus_activation_comparisons()
    # Join the held-out phenotype only after candidate definitions, resting
    # capacities, PKA folds, simulations, and capacity-interval points freeze.
    stimulus_rows = [
        replace(
            row,
            heldout_AE4_KO_ratio=HELDOUT_AE4_KO_FLOW_RATIO,
            heldout_residual_endpoint=(
                row.ae4_ko_endpoint_flow_ratio - HELDOUT_AE4_KO_FLOW_RATIO
                if math.isfinite(row.ae4_ko_endpoint_flow_ratio)
                else math.nan
            ),
            heldout_residual_integrated=(
                row.ae4_ko_integrated_flow_ratio - HELDOUT_AE4_KO_FLOW_RATIO
                if math.isfinite(row.ae4_ko_integrated_flow_ratio)
                else math.nan
            ),
        )
        for row in frozen_stimulus_rows
    ]

    all_rows = primary_rows + round2_rows
    all_runs = {**primary_runs, **round2_runs, **stimulus_runs}
    rounds = {row.mechanism_id: row.round for row in all_rows}
    rounds.update({row.mechanism_id: row.round for row in stimulus_rows})
    comparison_dicts = [row.as_dict() for row in all_rows]

    _write_csv(
        RESULTS / "model_comparison.csv",
        comparison_dicts,
        _comparison_columns(),
    )
    _write_json(
        RESULTS / "model_comparison.json",
        {
            "anti_leak_protocol": {
                "calibration": "WT balance plus per-candidate independent assay constants only",
                "heldout_names": HELDOUT_TARGET_NAMES,
                "heldout_value_read_after_freeze": HELDOUT_AE4_KO_FLOW_RATIO,
                "knockout_value_available_to_calibrator": False,
            },
            "rows": comparison_dicts,
        },
    )

    grid_columns = list(grid[0].keys())
    _write_csv(RESULTS / "round2_envelope.csv", grid, grid_columns)

    rest_dicts = [row.as_dict() for row in rest_rows]
    _write_csv(
        RESULTS / "resting_equilibria.csv",
        rest_dicts,
        _dataclass_columns(RestingEquilibrium),
    )
    _write_json(
        RESULTS / "resting_equilibria.json",
        {
            "targets": {
                "WT_Cl_mM": 50.10,
                "WT_Cl_SEM_mM": 1.50,
                "WT_pH": 6.91,
                "WT_pH_SEM": 0.07,
                "published_volume_pL": 1.3,
                "volume_sensitivity_relative": REST_VOLUME_RELATIVE_TOLERANCE,
            },
            "multistart": {
                "production_default_count": REST_MULTISTART_COUNT,
                "seed": REST_MULTISTART_SEED,
                "interpretation": "local basin robustness, not global uniqueness",
            },
            "knockout_target_read_by_solver": False,
            "rows": rest_dicts,
        },
    )
    genotype_rest = _genotype_rest_rows(rest_dicts)
    _write_csv(
        RESULTS / "genotype_resting_equilibria.csv",
        genotype_rest,
        _union_columns(genotype_rest),
    )
    _write_json(
        RESULTS / "genotype_resting_equilibria.json",
        {
            "primary_resting_targets": {
                "WT_Cl_mM": "50.10 +/- 1.50 SEM",
                "AE4_KO_Cl_mM": "36.50 +/- 1.60 SEM",
                "WT_pH": "6.91 +/- 0.07 SEM",
                "AE4_KO_pH": "6.89 +/- 0.02 SEM",
                "AE2_KO_Cl_mM": "54.50 +/- 1.80 SEM",
                "AE2_KO_pH": "6.95 +/- 0.05 SEM",
                "AE2_own_control_Cl_mM": "53.40 +/- 1.80 SEM",
                "AE2_own_control_pH": "6.87 +/- 0.01 SEM",
            },
            "interpretation": (
                "AE4-null root is candidate-independent; C1 and surviving 2:1:3 "
                "AE2-null roots are candidate/capacity-specific local predictions; "
                "all are distinct from stimulated endpoints"
            ),
            "rows": genotype_rest,
        },
    )

    _write_csv(
        RESULTS / "thermodynamic_rest_profile.csv",
        thermo_profile,
        _union_columns(thermo_profile),
    )
    _write_json(
        RESULTS / "thermodynamic_rest_profile.json",
        {"knockout_target_read_by_selector": False, "rows": thermo_profile},
    )
    _write_csv(
        RESULTS / "rest_capacity_envelope.csv",
        capacity_profile,
        _union_columns(capacity_profile),
    )
    _write_json(
        RESULTS / "rest_capacity_envelope.json",
        {
            "scope": (
                "post-hoc KO-blind WT-rest robustness continuation from 1e-8 "
                "to 1e4 times each volume-normalized capacity; axes and points "
                "frozen before held-out simulation"
            ),
            "gates": {
                "rest_cl_abs_SEM": 2.0,
                "rest_pH_abs_SEM": 2.0,
                "volume_relative": REST_VOLUME_RELATIVE_TOLERANCE,
                "volume_tolerance_provenance": (
                    "explicit new +/-10% robustness sensitivity, not source-reported uncertainty"
                ),
                "rest_raw_rhs": 1e-8,
            },
            "knockout_target_read_by_variant_selector": False,
            "rows": capacity_profile,
        },
    )
    stimulus_dicts = [row.as_dict() for row in stimulus_rows]
    _write_csv(
        RESULTS / "pka_activation_envelope.csv",
        stimulus_dicts,
        _dataclass_columns(StimulusComparison),
    )
    _write_json(
        RESULTS / "pka_activation_envelope.json",
        {
            "anti_leak_protocol": {
                "folds_and_capacity_points_frozen_before_heldout_join": True,
                "heldout_value_joined_after_simulation": HELDOUT_AE4_KO_FLOW_RATIO,
                "knockout_target_read_by_calibrator": False,
                "knockout_target_read_by_variant_selector": False,
            },
            "rows": stimulus_dicts,
        },
    )

    trajectory_columns = [
        "mechanism_id",
        "round",
        "scenario",
        "time",
        *STATE_ORDER,
        "calcium",
        "pka_activation",
        "q_total",
        "pH_i",
        "ae4_cycle_flux",
        "max_abs_raw_rhs",
        "max_abs_conservation_check",
        "solver_method",
        "non_ae4_chassis_fingerprint",
    ]
    _write_csv(
        RESULTS / "trajectories.csv",
        _trajectory_rows(all_runs, rounds),
        trajectory_columns,
    )

    base = FixedChassis()
    demand = base.ae4_source_demand()
    partition = base.baseline_partition_mismatch()
    correction = [
        -float(partition["amount_balance_residual"][name])
        for name in BALANCE_NAMES
    ]
    ideal_magnitude = float(partition["expected_na_k_magnitude"])
    ideal_correction = (-ideal_magnitude, ideal_magnitude, 0.0, 0.0)
    localization = pathway_localization(ideal_correction)
    localization_columns = list(localization[0].keys())
    _write_csv(
        RESULTS / "pathway_localization.csv",
        localization,
        localization_columns,
    )

    manifest = {
        "non_ae4_chassis_fingerprint": base.parameter_fingerprint,
        "fixed_WT_state": dict(zip(STATE_ORDER, map(float, BASELINE_STATE), strict=True)),
        "required_WT_AE4_contribution": dict(
            zip(BALANCE_NAMES, map(float, HISTORICAL_WT_AE4_CONTRIBUTION), strict=True)
        ),
        "chassis_computed_WT_AE4_demand": demand,
        "knockout_target_read_by_calibrator": False,
        "heldout_target_names_only": HELDOUT_TARGET_NAMES,
        "round2_grid_declared_without_KO": {
            "Na_capacity_fractions": ROUND2_NA_CAPACITY_FRACTIONS,
            "equilibrium_constants": ROUND2_EQUILIBRIUM_CONSTANTS,
            "saturation_strengths": ROUND2_SATURATION_STRENGTHS,
            "AE4_only_bath_HCO3_mM_sensitivities": ROUND2_BATH_HCO3_MILLIMOLAR,
            "K_eq_note": "values other than 1 are diagnostic energetic-bias requirements, not source-supported passive AE4 parameters",
            "bath_note": "40/42.9 mM alter only AE4 driving force and are not coherent whole-chassis baths",
        },
        "round3_PKA_and_moving_rest_protocol": {
            "PKA_folds": [1.0, *PKA_ASSAY_FOLDS, 10.0],
            "PKA_fold_1_6_note": "total-WT exchange sensitivity, not AE4-specific",
            "immediate_model_onset": PKA_IMMEDIATE_MODEL_ONSET,
            "time_note": "chassis code coordinate has no certified mapping to physical minutes",
            "capacity_continuation": (
                "post-hoc but KO-blind 1e-8--1e4x robustness envelope; axes, "
                "intervals, and log-space prediction points frozen before KO runs"
            ),
            "volume_tolerance": (
                "new +/-10% sensitivity around 1.3 pL; not a source-reported uncertainty"
            ),
            "heldout_value_available_to_variant_selector": False,
        },
        "solver": {
            "canonical_admissible": "SciPy BDF, rtol=atol=1e-6",
            "screening": "SciPy Radau, rtol=atol=1e-6",
            "moving_rest_PKA": "SciPy Radau, rtol=1e-7, atol=1e-8",
            "discontinuity_grid": "100+1e-7 right limit; post-step 0.1 spacing",
            "trajectory_export": (
                "fine 0.1-unit grid retained internally for integrals; CSV exports "
                "t=0,50,100/right-limit,150,200 plus both sides of any additional "
                "input discontinuity"
            ),
        },
    }
    _write_json(RESULTS / "calibration_manifest.json", manifest)
    _write_json(
        RESULTS / "solver_crosscheck.json",
        {
            "C1": c1_crosscheck,
            "round2_gentle_upper_envelope": best_crosscheck,
            "independent_C6_112_PKA_F3_BDF": {
                "wt_endpoint_flow": 0.001210263574,
                "ae4_ko_endpoint_flow": 0.001358119889,
                "ae4_ko_endpoint_ratio": 1.122168690,
                "ae4_ko_integrated_ratio": 1.130126790,
                "provenance": "independent adversarial Agent-E BDF rerun",
            },
            "independent_representative_213_BDF_Radau": {
                "cores": ["C6_213", "C7_213", "C8_213_SAT"],
                "capacity_points_per_core": 1,
                "folds": [1.0, 1.25, 1.6, 3.0, 10.0],
                "BDF_tolerances": {"rtol": 1e-7, "atol": 1e-8},
                "Radau_tolerances": {"rtol": 1e-7, "atol": 1e-8},
                "reported_upper_bound_absolute_endpoint_ratio_difference": 3e-10,
                "reported_upper_bound_absolute_integrated_ratio_difference": 3e-9,
                "reported_upper_bound_conservation_residual": 4.98e-14,
                "maximum_thermodynamic_violation": 0.0,
                "provenance": (
                    "independent numerical-reproducibility reviewer dual-solver "
                    "crosscheck at one in-band point per core (15 rows total); the "
                    "reviewer also ran BDF over all 30 positive AE4-specific "
                    "F=1.25/F=3 frozen-grid combinations, while this production run "
                    "supplies the matching full Radau grid"
                ),
            },
            "independent_source_positive_30_BDF_vs_production_Radau": {
                "cores": ["C6_213", "C7_213", "C8_213_SAT"],
                "capacity_points_per_core": 5,
                "folds": [1.25, 3.0],
                "row_count": 30,
                "BDF_tolerances": {"rtol": 1e-7, "atol": 1e-8},
                "Radau_tolerances": {"rtol": 1e-7, "atol": 1e-8},
                "all_BDF_endpoint_ratios_above_one": True,
                "all_BDF_integrated_ratios_above_one": True,
                "minimum_BDF_endpoint_ratio": 1.013777757953735,
                "minimum_BDF_integrated_ratio": 1.0209390004383792,
                "maximum_absolute_endpoint_ratio_difference": (
                    1.3592778014270834e-9
                ),
                "maximum_absolute_integrated_ratio_difference": (
                    3.370460355256455e-9
                ),
                "maximum_production_Radau_conservation_residual": (
                    4.973799150320701e-14
                ),
                "maximum_thermodynamic_violation": 0.0,
                "provenance": (
                    "independent numerical-reproducibility reviewer BDF rerun "
                    "of all five frozen capacity points for three cores and the "
                    "two positive AE4-specific folds (30 rows), compared with "
                    "the production Radau grid"
                ),
            },
        },
    )

    closed = [row for row in all_rows if row.wt_baseline_admissible]
    fully_supported = [
        record for record in grid if bool(record["fully_source_supported"])
    ]
    endpoint_best = min(closed, key=lambda row: row.ae4_ko_endpoint_flow_ratio)
    integrated_best = min(closed, key=lambda row: row.ae4_ko_integrated_flow_ratio)
    top_pathway = localization[0]
    immediate_rows = [
        row
        for row in stimulus_rows
        if row.simulation_status == "pass"
        and math.isclose(row.pka_activation_onset, PKA_IMMEDIATE_MODEL_ONSET)
    ]
    source_rest_rows = [
        row
        for row in immediate_rows
        if row.source_complete_under_new_volume_sensitivity
    ]
    strict_source_complete_rows = [row for row in immediate_rows if row.source_complete]
    source_moving_rows = [
        row
        for row in immediate_rows
        if row.assay_admissible
        and row.pka_stimulated_fold in (1.25, 3.0)
    ]
    source_rest_endpoint_best = min(
        source_rest_rows, key=lambda row: row.ae4_ko_endpoint_flow_ratio
    )
    source_rest_integrated_best = min(
        source_rest_rows, key=lambda row: row.ae4_ko_integrated_flow_ratio
    )
    source_moving_endpoint_best = min(
        source_moving_rows, key=lambda row: row.ae4_ko_endpoint_flow_ratio
    )
    no_activation_rest_rows = [
        row
        for row in immediate_rows
        if row.resting_quantified_gate_pass
        and row.assay_admissible
        and row.thermodynamic_status == "pass"
        and math.isclose(row.pka_stimulated_fold, 1.0)
    ]
    no_activation_endpoint_best = min(
        no_activation_rest_rows, key=lambda row: row.ae4_ko_endpoint_flow_ratio
    )
    total_wt_f16_rows = [
        row
        for row in immediate_rows
        if row.resting_quantified_gate_pass
        and row.assay_admissible
        and row.thermodynamic_status == "pass"
        and math.isclose(row.pka_stimulated_fold, 1.6)
    ]
    total_wt_f16_endpoint_best = min(
        total_wt_f16_rows, key=lambda row: row.ae4_ko_endpoint_flow_ratio
    )
    total_wt_f16_integrated_best = min(
        total_wt_f16_rows, key=lambda row: row.ae4_ko_integrated_flow_ratio
    )
    interval_summaries = []
    for core_id in ("C6_213", "C7_213", "C8_213_SAT"):
        core_rows = [row for row in source_rest_rows if row.core_candidate_id == core_id]
        if not core_rows:
            continue
        capacity_rows = [
            row
            for row in capacity_profile
            if row.get("core_candidate_id") == core_id
            and row.get("profile_kind") == "frozen_interval_prediction_point"
        ]
        interval_summaries.append(
            {
                "core_candidate_id": core_id,
                "capacity_interval_lower": min(
                    float(row["capacity_interval_lower"]) for row in capacity_rows
                ),
                "capacity_interval_upper": max(
                    float(row["capacity_interval_upper"]) for row in capacity_rows
                ),
                "prediction_points": len(
                    {float(row.activity_scale) for row in core_rows}
                ),
                "endpoint_ratio_min": min(
                    row.ae4_ko_endpoint_flow_ratio for row in core_rows
                ),
                "endpoint_ratio_max": max(
                    row.ae4_ko_endpoint_flow_ratio for row in core_rows
                ),
                "integrated_ratio_min": min(
                    row.ae4_ko_integrated_flow_ratio for row in core_rows
                ),
                "integrated_ratio_max": max(
                    row.ae4_ko_integrated_flow_ratio for row in core_rows
                ),
                "fold3_endpoint_ratio_min": min(
                    row.ae4_ko_endpoint_flow_ratio
                    for row in core_rows
                    if math.isclose(row.pka_stimulated_fold, 3.0)
                ),
                "fold3_endpoint_ratio_max": max(
                    row.ae4_ko_endpoint_flow_ratio
                    for row in core_rows
                    if math.isclose(row.pka_stimulated_fold, 3.0)
                ),
                "fold3_integrated_ratio_min": min(
                    row.ae4_ko_integrated_flow_ratio
                    for row in core_rows
                    if math.isclose(row.pka_stimulated_fold, 3.0)
                ),
                "fold3_integrated_ratio_max": max(
                    row.ae4_ko_integrated_flow_ratio
                    for row in core_rows
                    if math.isclose(row.pka_stimulated_fold, 3.0)
                ),
                "ae2_KO_endpoint_ratio_min": min(
                    row.ae2_ko_endpoint_flow_ratio for row in core_rows
                ),
                "ae2_KO_endpoint_ratio_max": max(
                    row.ae2_ko_endpoint_flow_ratio for row in core_rows
                ),
            }
        )
    ae4_ko_rest = next(
        row for row in rest_rows if row.mechanism_id == "COMMON_AE4_KO_REST"
    )
    experimental_upper = HELDOUT_AE4_KO_FLOW_RATIO + 0.047
    flow_compatible_rows = [
        row
        for row in source_rest_rows
        if HELDOUT_AE4_KO_FLOW_RATIO - 0.047
        <= row.ae4_ko_endpoint_flow_ratio
        <= experimental_upper
        and HELDOUT_AE4_KO_FLOW_RATIO - 0.047
        <= row.ae4_ko_integrated_flow_ratio
        <= experimental_upper
    ]
    ae4_ko_genotype_rest = next(
        row for row in genotype_rest if row["mechanism_id"] == "COMMON_AE4_KO_REST"
    )
    capacity_positions = (
        "near_lower_boundary",
        "lower_quartile",
        "midpoint",
        "upper_quartile",
        "near_upper_boundary",
    )
    ae2_rest_pass_by_core_position: dict[tuple[str, str], bool] = {}
    for genotype_row in genotype_rest:
        mechanism_id = str(genotype_row["mechanism_id"])
        for core_id in ("C6_213", "C7_213", "C8_213_SAT"):
            for position in capacity_positions:
                if mechanism_id == f"{core_id}_{position}_AE2_KO_REST":
                    ae2_rest_pass_by_core_position[(core_id, position)] = bool(
                        genotype_row["genotype_rest_gate_pass"]
                    )

    def matching_ae2_rest_pass(row: StimulusComparison) -> bool:
        position = next(
            (
                item
                for item in capacity_positions
                if f"_capacity_{item}_" in row.mechanism_id
            ),
            "",
        )
        return ae2_rest_pass_by_core_position.get(
            (row.core_candidate_id, position), False
        )

    # Flow agreement alone can never identify a mechanism.  The global AE4KO
    # rest phenotype, candidate-specific AE2KO rest/dynamic phenotype,
    # conservation, thermodynamics, and audit are additional hard gates.
    hard_gate_success_rows = [
        row
        for row in flow_compatible_rows
        if bool(ae4_ko_genotype_rest["genotype_rest_gate_pass"])
        and matching_ae2_rest_pass(row)
        and abs(row.ae2_ko_endpoint_flow_ratio - 1.0) <= 0.05
        and row.max_conservation_residual <= 1e-10
        and row.max_thermodynamic_violation <= 1e-12
    ]
    substantive_conclusion = (
        "INTERNAL GATES PASS; EXTERNAL AUDIT AND COMPLEXITY REVIEW REQUIRED"
        if hard_gate_success_rows
        else "AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED"
    )
    survivor_ae2_rest_rows = [
        row
        for row in genotype_rest
        if str(row["mechanism_id"]).startswith(("C6_213_", "C7_213_", "C8_213_SAT_"))
    ]
    summary = {
        "substantive_conclusion": substantive_conclusion,
        "scope": (
            "fixed historical seven-state non-AE4 vector field; fixed-row, moved-root, "
            "post-hoc KO-blind finite capacity continuation, and source-bounded PKA gates"
        ),
        "heldout_AE4_KO_flow_ratio": HELDOUT_AE4_KO_FLOW_RATIO,
        "heldout_AE4_KO_flow_ratio_SEM": 0.047,
        "primary_candidates_tested": len(primary_rows),
        "round2_grid_rows": len(grid),
        "WT_closed_round2_rows": len(round2_rows),
        "fully_source_supported_and_WT_closed_rows": len(fully_supported),
        "source_rest_thermo_PKA_rows_under_new_volume_sensitivity": len(source_rest_rows),
        "strict_source_complete_rows_without_new_volume_sensitivity": len(
            strict_source_complete_rows
        ),
        "source_rest_thermo_PKA_mechanism_cores_under_new_volume_sensitivity": sorted(
            {row.core_candidate_id for row in source_rest_rows}
        ),
        "flow_compatible_source_rest_rows": len(flow_compatible_rows),
        "all_hard_gate_success_rows": len(hard_gate_success_rows),
        "hard_gate_failures": {
            "candidate_independent_AE4_KO_rest_gate_pass": bool(
                ae4_ko_genotype_rest["genotype_rest_gate_pass"]
            ),
            "candidate_specific_AE2_KO_rest_passing_rows": sum(
                bool(row["genotype_rest_gate_pass"])
                for row in genotype_rest
                if str(row["mechanism_id"]).startswith(
                    ("C6_213_", "C7_213_", "C8_213_SAT_")
                )
            ),
            "note": "flow compatibility is necessary but not sufficient",
            "positive_identification_requires_external_adversarial_audit_and_complexity_review": True,
            "AE2_dynamic_tolerance_note": (
                "+/-5% is an explicit new qualitative-normality sensitivity, "
                "not an experimental confidence interval"
            ),
        },
        "C1": c1_row.as_dict(),
        "generous_closed_family_endpoint_bound": {
            "mechanism_id": endpoint_best.mechanism_id,
            "KO_over_WT": endpoint_best.ae4_ko_endpoint_flow_ratio,
            "maximum_reduction_fraction": 1.0 - endpoint_best.ae4_ko_endpoint_flow_ratio,
            "source_admissible": False,
        },
        "generous_closed_family_integrated_bound": {
            "mechanism_id": integrated_best.mechanism_id,
            "KO_over_WT": integrated_best.ae4_ko_integrated_flow_ratio,
            "maximum_reduction_fraction": max(
                0.0, 1.0 - integrated_best.ae4_ko_integrated_flow_ratio
            ),
            "source_admissible": False,
        },
        "quantified_maxima_by_domain": {
            "strict_source_plus_rest_plus_thermodynamic_gates": {
                "row_count": len(strict_source_complete_rows),
                "quantified_bound": None,
                "reason": (
                    "strict set is empty; no numerical bound is inferred from an empty set"
                ),
            },
            "source_mechanism_plus_rest_thermo_under_new_volume_sensitivity": {
                "tested_domain": (
                    "C6/C7/C8-SAT 2:1:3 connected rest-admissible capacity bands "
                    "within the finite 1e-8--1e4x continuation using a new post-hoc "
                    "+/-10% volume sensitivity; independently AE4-specific folds 1.25 and 3"
                ),
                "minimum_endpoint_KO_over_WT": source_rest_endpoint_best.ae4_ko_endpoint_flow_ratio,
                "minimum_endpoint_row": source_rest_endpoint_best.mechanism_id,
                "maximum_endpoint_reduction_fraction": max(
                    0.0, 1.0 - source_rest_endpoint_best.ae4_ko_endpoint_flow_ratio
                ),
                "minimum_integrated_KO_over_WT": source_rest_integrated_best.ae4_ko_integrated_flow_ratio,
                "minimum_integrated_row": source_rest_integrated_best.mechanism_id,
                "maximum_integrated_reduction_fraction": max(
                    0.0, 1.0 - source_rest_integrated_best.ae4_ko_integrated_flow_ratio
                ),
                "no_activation_reference_minimum_endpoint_ratio": (
                    no_activation_endpoint_best.ae4_ko_endpoint_flow_ratio
                ),
                "no_activation_reference_row": no_activation_endpoint_best.mechanism_id,
                "F1_6_total_WT_non_AE4_specific_sensitivity": {
                    "minimum_endpoint_ratio": (
                        total_wt_f16_endpoint_best.ae4_ko_endpoint_flow_ratio
                    ),
                    "minimum_endpoint_row": total_wt_f16_endpoint_best.mechanism_id,
                    "minimum_integrated_ratio": (
                        total_wt_f16_integrated_best.ae4_ko_integrated_flow_ratio
                    ),
                    "minimum_integrated_row": total_wt_f16_integrated_best.mechanism_id,
                    "interpretation": (
                        "total-WT exchanger evidence; excluded from the AE4-specific positive-fold bound"
                    ),
                },
            },
            "source_supported_moving_roots_even_if_rest_gate_fails": {
                "minimum_endpoint_KO_over_WT": source_moving_endpoint_best.ae4_ko_endpoint_flow_ratio,
                "minimum_endpoint_row": source_moving_endpoint_best.mechanism_id,
                "maximum_endpoint_reduction_fraction": max(
                    0.0, 1.0 - source_moving_endpoint_best.ae4_ko_endpoint_flow_ratio
                ),
            },
            "generous_static_diagnostic_envelope": {
                "minimum_endpoint_KO_over_WT": endpoint_best.ae4_ko_endpoint_flow_ratio,
                "minimum_endpoint_row": endpoint_best.mechanism_id,
                "maximum_endpoint_reduction_fraction": 1.0
                - endpoint_best.ae4_ko_endpoint_flow_ratio,
                "source_admissible": False,
            },
            "PKA_and_novel_K_recruitment": {
                "source_rest_capacity_band_summaries": interval_summaries,
                "interpretation": (
                    "all independently bounded positive activation folds move the "
                    "rest-admissible thermodynamic predictions in the wrong direction"
                ),
            },
        },
        "inverse_phenotype_diagnostic": {
            "required_C1_WT_endpoint_gain_factor": c1_row.ae4_ko_endpoint_flow_ratio
            / HELDOUT_AE4_KO_FLOW_RATIO,
            "required_C1_WT_integrated_gain_factor": c1_row.ae4_ko_integrated_flow_ratio
            / HELDOUT_AE4_KO_FLOW_RATIO,
            "interpretation": "diagnostic only; not an optimizer target or calibrated parameter",
        },
        "transported_K_fraction_upper_bound_at_WT_vector_tolerance": transported_k_fraction_bound(),
        "common_KO_endpoint_flow": c1_row.ae4_ko_endpoint_flow,
        "candidate_independent_AE4_KO_rest_falsification": {
            "state": list(ae4_ko_rest.state),
            "cell_volume_pL": ae4_ko_rest.cell_volume_pL,
            "cl_i_mM": ae4_ko_rest.cl_i_mM,
            "primary_AE4_KO_cl_i_mM": 36.50,
            "pH_i": ae4_ko_rest.pH_i,
            "primary_AE4_KO_pH": 6.89,
            "interpretation": (
                "AE4-null rest is candidate-independent and falsifies the direct Cl/pH "
                "phenotypes; volume/HCO3 are separate model-reference diagnostics. "
                "This localizes a missing non-AE4 acid-base/chloride module."
            ),
        },
        "candidate_specific_survivor_AE2_KO_rest_falsification": {
            "row_count": len(survivor_ae2_rest_rows),
            "direct_Cl_target_mM": 54.50,
            "direct_Cl_SEM_mM": 1.80,
            "direct_pH_target": 6.95,
            "direct_pH_SEM": 0.05,
            "all_genotype_rest_gates_pass": all(
                bool(row["genotype_rest_gate_pass"])
                for row in survivor_ae2_rest_rows
            ),
            "passing_row_count": sum(
                bool(row["genotype_rest_gate_pass"])
                for row in survivor_ae2_rest_rows
            ),
            "cl_i_mM_range": [
                min(float(row["cl_i_mM"]) for row in survivor_ae2_rest_rows),
                max(float(row["cl_i_mM"]) for row in survivor_ae2_rest_rows),
            ],
            "cl_target_z_range": [
                min(float(row["genotype_cl_target_z"]) for row in survivor_ae2_rest_rows),
                max(float(row["genotype_cl_target_z"]) for row in survivor_ae2_rest_rows),
            ],
            "pH_target_z_range": [
                min(float(row["genotype_pH_target_z"]) for row in survivor_ae2_rest_rows),
                max(float(row["genotype_pH_target_z"]) for row in survivor_ae2_rest_rows),
            ],
            "scope": "candidate/capacity-specific local roots; no genotype target fit",
        },
        "pathway_localization": {
            "baseline_correction_only": True,
            "top_one_component": top_pathway,
            "no_single_component_exact": top_pathway[
                "best_one_component_residual_fraction"
            ]
            > 0.0,
            "localized_module": (
                "acid-base/chloride homeostasis module (specific protein unresolved), with a "
                "separate two-dimensional Na/K correction needed for exact mixed-cation WT closure"
            ),
            "phenotype_caveat": "balance alignment localizes WT closure, not a fitted explanation of the held-out phenotype",
        },
        "round4": {
            "status": "analytically screened; no admissible one-component correction",
            "reason": (
                "NHE1 alone cannot cancel the K residual; NaK-ATPase alone leaves "
                "1/sqrt(26)=19.6116% of the cation correction. Exact repair needs "
                "two independent non-AE4 directions, and no source-calibrated dynamic magnitude exists."
            ),
        },
    }
    _write_json(RESULTS / "summary.json", summary)

    _write_figures(primary_rows, round2_rows, grid, stimulus_rows, capacity_profile)
    return summary


def main() -> None:
    summary = run()
    print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
