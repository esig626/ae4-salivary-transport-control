"""WT-only cross-root G5 repair profile for Task 13B.

This runner reads only the frozen WT calibration panel.  It screens every
rest-passing branch under the fixed 600-s CCh+IPR protocol and the fixed
9--10 uL/min graphical minute envelope.  The source-motivated N0/N1/N2 NKCC1
nesting is evaluated with the mandatory R1 beta/cAMP/PKA-to-AE4 reference.
No null-genotype model or phenotype target is loaded here.

The selected WT root receives the exhaustive NKCC1 profile in
``run_nkcc_stimulation_profile.py``.  This cross-root runner uses a declared
representative panel spanning the nominal P15 scale and two deliberately
out-of-scale flow-shape probes.  A repair is eligible only at the historical
0.058-to-0.55 uM calcium anchor and within the nominal NKCC1 gain profile.
The 0.35/0.75 uM calcium rows are Tier-3 diagnostics and can never select a
production model.
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
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .calibration import ModelVariant, WTCalibrationParameters, build_wt_model
from .camp_pka import R1EffectiveActivation
from .model import ModernFullModel, WT
from .nkcc_stimulation import (
    N0BaselineNkcc1,
    N1AlgebraicNkcc1,
    N2EffectiveDynamicNkcc1,
    attach_stimulated_nkcc1,
    isolated_nkcc1_fractional_cl_uptake_s,
)
from .states import CORE_STATE_NAMES
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    WTFlowEnvelope,
    fit_wt_flow_scale,
    physical_time_grid,
    sha256_file,
    simulate_wt,
    wt_flow_envelope_feasibility,
)


SELECTED_MEMBER_ID = "R1_G125_P21_PROSE_TREFERENCE"
PROFILE_FILENAME = "wt_g5_cross_root_repair_profile.csv"
SOLVER_FILENAME = "wt_g5_cross_root_repair_solver.csv"
SUMMARY_FILENAME = "wt_g5_cross_root_repair_summary.json"
HASH_FILENAME = "wt_g5_cross_root_repair_hashes.csv"


@dataclass(frozen=True)
class FrozenRootCase:
    root_id: str
    variant: ModelVariant
    calibration: WTCalibrationParameters
    core_state: tuple[float, ...]
    selected_reference: bool


@dataclass(frozen=True)
class NkccRepairSpec:
    member_id: str
    family: str
    fully_activated_multiplier: float
    tau_activation_s: float | None
    gain_evidence_status: str
    kinetic_evidence_status: str

    @property
    def nominal_gain_profile(self) -> bool:
        return self.fully_activated_multiplier <= 1.75

    def construct(self, *, stimulated_calcium_uM: float = 0.55) -> Any:
        if self.family == "N0_BASELINE":
            return N0BaselineNkcc1()
        if self.family == "N1_CCH_ALGEBRAIC":
            return N1AlgebraicNkcc1(
                fully_activated_multiplier=self.fully_activated_multiplier,
                stimulated_calcium_uM=stimulated_calcium_uM,
            )
        if self.family == "N2_CCH_EFFECTIVE_DYNAMIC":
            if self.tau_activation_s is None:
                raise ValueError("N2 requires an explicit sensitivity time")
            return N2EffectiveDynamicNkcc1(
                fully_activated_multiplier=self.fully_activated_multiplier,
                tau_activation_s=self.tau_activation_s,
                stimulated_calcium_uM=stimulated_calcium_uM,
            )
        raise ValueError(f"unknown NKCC1 repair family {self.family!r}")


@dataclass(frozen=True)
class CalciumProbe:
    label: str
    stimulated_calcium_uM: float
    evidence_status: str
    eligible_for_selection: bool


@dataclass(frozen=True)
class SimulationRequest:
    root: FrozenRootCase
    nkcc: NkccRepairSpec
    calcium: CalciumProbe
    panel_stage: str
    solver_label: str = PRODUCTION_RADAU.label


def declared_cross_root_nkcc_panel() -> tuple[NkccRepairSpec, ...]:
    """Representative predeclared N0/N1/N2 panel used on all nine roots."""

    gains = (1.75, 3.0, 8.0)
    taus_s = (10.0, 30.0, 120.0, 300.0)
    rows = [
        NkccRepairSpec(
            member_id="N0_BASELINE_M1",
            family="N0_BASELINE",
            fully_activated_multiplier=1.0,
            tau_activation_s=None,
            gain_evidence_status="EXACT_BASELINE_CONTROL",
            kinetic_evidence_status="NO_ADDED_KINETIC_STATE",
        )
    ]
    for multiplier in gains:
        gain_status = (
            "P15_DIRECTION_NOMINAL_GAIN_PROFILE"
            if multiplier <= 1.75
            else "FLOW_SHAPE_PROBE_OUTSIDE_P15_NOMINAL_SCALE"
        )
        token = str(multiplier).replace(".", "P")
        rows.append(
            NkccRepairSpec(
                member_id=f"N1_M{token}",
                family="N1_CCH_ALGEBRAIC",
                fully_activated_multiplier=multiplier,
                tau_activation_s=None,
                gain_evidence_status=gain_status,
                kinetic_evidence_status="ALGEBRAIC_NO_ADDED_TIME_CONSTANT",
            )
        )
        for tau_s in taus_s:
            rows.append(
                NkccRepairSpec(
                    member_id=f"N2_M{token}_TAU{int(tau_s)}S",
                    family="N2_CCH_EFFECTIVE_DYNAMIC",
                    fully_activated_multiplier=multiplier,
                    tau_activation_s=tau_s,
                    gain_evidence_status=gain_status,
                    kinetic_evidence_status=(
                        "UNMEASURED_PREDECLARED_WT_SENSITIVITY_NOT_FIT"
                    ),
                )
            )
    if len(rows) != 16:
        raise AssertionError("cross-root panel must contain 16 nested members")
    return tuple(rows)


def declared_calcium_probes() -> tuple[CalciumProbe, ...]:
    """Historical anchor plus two explicitly non-selecting Tier-3 probes."""

    return (
        CalciumProbe(
            "CA035_TIER3_LOW_PROBE",
            0.35,
            "TIER_3_MODEL_PROBE_NOT_HISTORICAL_MEASUREMENT",
            False,
        ),
        CalciumProbe(
            "CA055_HISTORICAL_ANCHOR",
            0.55,
            "HISTORICAL_IMPLEMENTATION_AMPLITUDE_SENSITIVITY",
            True,
        ),
        CalciumProbe(
            "CA075_TIER3_HIGH_PROBE",
            0.75,
            "TIER_3_MODEL_PROBE_NOT_HISTORICAL_MEASUREMENT",
            False,
        ),
    )


def _bool(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"expected serialized boolean, got {value!r}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_retained_wt_roots(results_directory: Path) -> tuple[FrozenRootCase, ...]:
    """Reconstruct all and only rest-passing roots from frozen WT artifacts."""

    table_path = results_directory / "wt_root_table.csv"
    state_path = results_directory / "wt_root_states.csv"
    parameter_path = results_directory / "wt_parameter_candidates.csv"
    selected_path = results_directory / "wt_selected_candidate.json"
    selected = json.loads(selected_path.read_text(encoding="utf-8"))
    selected_id = str(selected["construction"]["root_id"])

    states: dict[str, dict[str, float]] = {}
    for row in _read_rows(state_path):
        states.setdefault(row["root_id"], {})[row["state_name"]] = float(row["value"])
    parameters: dict[str, dict[str, float]] = {}
    for row in _read_rows(parameter_path):
        parameters.setdefault(row["root_id"], {})[row["parameter"]] = float(row["value"])

    retained: list[FrozenRootCase] = []
    for row in _read_rows(table_path):
        if not row["root_id"] or not _bool(row["passes_wt_gate"]):
            continue
        root_id = row["root_id"]
        state_values = states.get(root_id, {})
        if set(state_values) != set(CORE_STATE_NAMES):
            raise ValueError(f"{root_id}: incomplete or unexpected core-state ledger")
        parameter_values = parameters.get(root_id, {})
        required_parameters = set(WTCalibrationParameters.__dataclass_fields__)
        if not required_parameters.issubset(parameter_values):
            raise ValueError(f"{root_id}: incomplete calibration-parameter ledger")
        retained.append(
            FrozenRootCase(
                root_id=root_id,
                variant=ModelVariant(
                    variant_id=row["variant_id"],
                    mixed_bath_na_attempt_fraction=float(row["routing_na_fraction"]),
                    apical_pump_fraction=float(row["apical_pump_fraction"]),
                    apical_k_fraction=float(row["apical_k_fraction"]),
                ),
                calibration=WTCalibrationParameters(
                    **{
                        name: parameter_values[name]
                        for name in WTCalibrationParameters.__dataclass_fields__
                    }
                ),
                core_state=tuple(state_values[name] for name in CORE_STATE_NAMES),
                selected_reference=root_id == selected_id,
            )
        )
    retained.sort(key=lambda item: item.root_id)
    if len(retained) != 9:
        raise AssertionError(f"expected nine frozen rest-passing roots, found {len(retained)}")
    if sum(item.selected_reference for item in retained) != 1:
        raise AssertionError("selected WT reference must occur exactly once in retained panel")
    return tuple(retained)


def _build_model(request: SimulationRequest) -> Any:
    ae4_regulation = R1EffectiveActivation(tau_activation_s=30.0)
    resting = build_wt_model(
        request.root.variant,
        request.root.calibration,
        regulatory_model=ae4_regulation,
    )
    protocol = SecretagogueProtocol(
        stimulated_calcium_uM=request.calcium.stimulated_calcium_uM
    )
    template = ModernFullModel(
        parameters=resting.parameters,
        stimulus=protocol,
        regulatory_model=ae4_regulation,
        ae4_parameters=resting.ae4_parameters,
        ae4_evaluator=resting.ae4_evaluator,
    )
    return attach_stimulated_nkcc1(
        template,
        request.nkcc.construct(
            stimulated_calcium_uM=request.calcium.stimulated_calcium_uM
        ),
    )


def _empty_metrics(request: SimulationRequest, message: str) -> dict[str, Any]:
    return {
        **_request_metadata(request),
        "success": False,
        "message": message,
        "positive_core": False,
        "numerical_gate_pass": False,
        "minute_flow_shape_gate_pass": False,
        "source_eligible_for_g5_selection": False,
        "g5_repair_candidate_pass": False,
    }


def _request_metadata(request: SimulationRequest) -> dict[str, Any]:
    return {
        "root_id": request.root.root_id,
        "selected_reference": request.root.selected_reference,
        "variant_id": request.root.variant.variant_id,
        "routing_na_fraction": request.root.variant.mixed_bath_na_attempt_fraction,
        "apical_pump_fraction": request.root.variant.apical_pump_fraction,
        "apical_k_fraction": request.root.variant.apical_k_fraction,
        "panel_stage": request.panel_stage,
        "calcium_profile": request.calcium.label,
        "stimulated_calcium_uM": request.calcium.stimulated_calcium_uM,
        "calcium_evidence_status": request.calcium.evidence_status,
        "nkcc_member_id": request.nkcc.member_id,
        "nkcc_family": request.nkcc.family,
        "nkcc_fully_activated_multiplier": request.nkcc.fully_activated_multiplier,
        "nkcc_tau_activation_s": request.nkcc.tau_activation_s,
        "nkcc_gain_evidence_status": request.nkcc.gain_evidence_status,
        "nkcc_kinetic_evidence_status": request.nkcc.kinetic_evidence_status,
        "ae4_regulatory_member": SELECTED_MEMBER_ID,
        "solver_label": request.solver_label,
    }


def _simulate_request(request: SimulationRequest) -> dict[str, Any]:
    try:
        model = _build_model(request)
        solver = (
            PRODUCTION_RADAU
            if request.solver_label == PRODUCTION_RADAU.label
            else PRODUCTION_BDF
        )
        grid = physical_time_grid(duration_s=600.0, sample_step_s=5.0)
        trajectory = simulate_wt(
            model,
            request.root.core_state,
            family=f"{SELECTED_MEMBER_ID}+{request.nkcc.member_id}",
            initial_condition="baseline",
            solver=solver,
            time_s=grid,
        )
        if not trajectory.success:
            return _empty_metrics(request, trajectory.message)
        minute_times = np.arange(60.0, 601.0, 60.0)
        minute_flow = np.interp(minute_times, trajectory.time_s, trajectory.flow_pL_s)
        flow_shape = wt_flow_envelope_feasibility(trajectory, WTFlowEnvelope())
        flow_scale = fit_wt_flow_scale(trajectory, WTFlowEnvelope())
        scaled_minute = flow_scale.convert(minute_flow)
        conservation_ratios = {
            name: trajectory.max_abs_conservation_residuals[name] / tolerance
            for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
        }
        max_conservation_ratio = max(conservation_ratios.values())
        numerical_gate = bool(
            trajectory.success
            and trajectory.positive_core
            and np.all(trajectory.flow_pL_s >= 0.0)
            and max_conservation_ratio <= 1.0
        )
        source_eligible = bool(
            request.calcium.eligible_for_selection
            and request.nkcc.nominal_gain_profile
        )
        indices = {
            "onset": 1,
            "minute_1": int(np.argmin(np.abs(trajectory.time_s - 60.0))),
            "minute_10": int(np.argmin(np.abs(trajectory.time_s - 600.0))),
        }
        evaluations = {
            label: model.evaluate(
                float(trajectory.time_s[index]),
                trajectory.states[:, index],
                genotype=WT,
            )
            for label, index in indices.items()
        }
        onset_uptake = isolated_nkcc1_fractional_cl_uptake_s(
            model,
            float(trajectory.time_s[indices["onset"]]),
            trajectory.states[:, indices["onset"]],
            genotype=WT,
        )
        endpoint_uptake = isolated_nkcc1_fractional_cl_uptake_s(
            model,
            float(trajectory.time_s[indices["minute_10"]]),
            trajectory.states[:, indices["minute_10"]],
            genotype=WT,
        )
        return {
            **_request_metadata(request),
            "success": trajectory.success,
            "message": trajectory.message,
            "positive_core": trajectory.positive_core,
            "all_flow_nonnegative": bool(np.all(trajectory.flow_pL_s >= 0.0)),
            "max_dimensionless_conservation_ratio": max_conservation_ratio,
            "numerical_gate_pass": numerical_gate,
            "minute_flow_60_pL_s": float(minute_flow[0]),
            "minute_flow_600_pL_s": float(minute_flow[-1]),
            "minute_flow_max_min_ratio": flow_shape.model_max_min_ratio,
            "envelope_max_min_ratio": flow_shape.envelope_max_min_ratio,
            "feasible_multiplier_lower_uL_min_per_pL_s": (
                flow_shape.feasible_multiplier_lower_uL_min_per_pL_s
            ),
            "feasible_multiplier_upper_uL_min_per_pL_s": (
                flow_shape.feasible_multiplier_upper_uL_min_per_pL_s
            ),
            "minute_flow_shape_gate_pass": flow_shape.compatible,
            "mean_scale_multiplier_uL_min_per_pL_s": (
                flow_scale.multiplier_uL_min_per_pL_s
            ),
            "mean_scale_implied_effective_cell_count": flow_scale.effective_cell_count,
            "implied_cell_count_status": "WARNING_NO_PRIMARY_GEOMETRY_RANGE",
            "scaled_minute_min_uL_min": float(np.min(scaled_minute)),
            "scaled_minute_max_uL_min": float(np.max(scaled_minute)),
            "minutes_in_9_10_envelope": int(
                np.sum((scaled_minute >= 9.0) & (scaled_minute <= 10.0))
            ),
            "nkcc_capacity_multiplier_onset": float(
                evaluations["onset"].diagnostics.regulatory[
                    "nkcc1_capacity_multiplier"
                ]
            ),
            "nkcc_capacity_multiplier_60_s": float(
                evaluations["minute_1"].diagnostics.regulatory[
                    "nkcc1_capacity_multiplier"
                ]
            ),
            "nkcc_capacity_multiplier_600_s": float(
                evaluations["minute_10"].diagnostics.regulatory[
                    "nkcc1_capacity_multiplier"
                ]
            ),
            "nkcc_flux_60_fmol_s": float(
                evaluations["minute_1"].diagnostics.homeostasis.nkcc1_inward_fmol_s
            ),
            "nkcc_flux_600_fmol_s": float(
                evaluations["minute_10"].diagnostics.homeostasis.nkcc1_inward_fmol_s
            ),
            "isolated_nkcc_fractional_cl_rate_onset_s_inv": onset_uptake,
            "isolated_nkcc_fractional_cl_rate_600_s_inv": endpoint_uptake,
            "uptake_observation_status": (
                "MECHANISTIC_DIAGNOSTIC_NOT_EQUATED_TO_SPQ_WITHOUT_MAP"
            ),
            "endpoint_cl_mM": float(trajectory.cell_cl_mM[-1]),
            "endpoint_ph": float(trajectory.cell_ph[-1]),
            "endpoint_volume_pL": float(trajectory.cell_volume_pL[-1]),
            "endpoint_flow_pL_s": float(trajectory.flow_pL_s[-1]),
            "cumulative_flow_pL": float(trajectory.cumulative_flow_pL[-1]),
            "source_eligible_for_g5_selection": source_eligible,
            "g5_repair_candidate_pass": bool(
                numerical_gate and flow_shape.compatible and source_eligible
            ),
        }
    except Exception as exc:  # profile every declared case and record failures
        return _empty_metrics(request, f"{type(exc).__name__}: {exc}")


def _run_requests(
    requests: Sequence[SimulationRequest], workers: int
) -> list[dict[str, Any]]:
    if workers <= 1:
        return [_simulate_request(request) for request in requests]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_simulate_request, requests, chunksize=1))


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


def _best_row(rows: Iterable[Mapping[str, Any]]) -> Mapping[str, Any]:
    usable = [
        row
        for row in rows
        if bool(row.get("numerical_gate_pass"))
        and math.isfinite(float(row.get("minute_flow_max_min_ratio", math.inf)))
    ]
    if not usable:
        raise ValueError("no numerically usable row in declared subset")
    return min(usable, key=lambda row: float(row["minute_flow_max_min_ratio"]))


def run_cross_root_profile(results_directory: Path, *, workers: int = 1) -> Mapping[str, Any]:
    roots = load_retained_wt_roots(results_directory)
    nkcc_panel = declared_cross_root_nkcc_panel()
    calcium = {item.label: item for item in declared_calcium_probes()}
    historical = calcium["CA055_HISTORICAL_ANCHOR"]

    primary_requests = tuple(
        SimulationRequest(root, member, historical, "PRIMARY_CROSS_ROOT_NKCC_SCREEN")
        for root in roots
        for member in nkcc_panel
    )
    primary_rows = _run_requests(primary_requests, workers)

    # Calcium amplitude probes are deliberately post-screen diagnostics.  The
    # deterministic rule is frozen here: for each root, profile the numerically
    # best historical-anchor law under 0.35 and 0.75 uM.  These rows can never
    # select a production candidate.
    request_by_key = {
        (request.root.root_id, request.nkcc.member_id): request
        for request in primary_requests
    }
    calcium_requests: list[SimulationRequest] = []
    for root in roots:
        root_rows = [row for row in primary_rows if row["root_id"] == root.root_id]
        winner = _best_row(root_rows)
        base_request = request_by_key[(root.root_id, str(winner["nkcc_member_id"]))]
        for label in ("CA035_TIER3_LOW_PROBE", "CA075_TIER3_HIGH_PROBE"):
            calcium_requests.append(
                SimulationRequest(
                    base_request.root,
                    base_request.nkcc,
                    calcium[label],
                    "POSTSCREEN_CALCIUM_TIER3_PROBE",
                )
            )
    calcium_rows = _run_requests(tuple(calcium_requests), workers)
    rows = primary_rows + calcium_rows

    # Dual-method checks cover both the best nominal and best unrestricted
    # historical-anchor row on every retained root.  Duplicates are removed.
    solver_requests: dict[tuple[str, str], SimulationRequest] = {}
    for root in roots:
        root_primary = [row for row in primary_rows if row["root_id"] == root.root_id]
        unrestricted = _best_row(root_primary)
        nominal = _best_row(
            row for row in root_primary if bool(row["source_eligible_for_g5_selection"])
        )
        for winner in (unrestricted, nominal):
            base = request_by_key[(root.root_id, str(winner["nkcc_member_id"]))]
            solver_requests[(root.root_id, base.nkcc.member_id)] = SimulationRequest(
                base.root,
                base.nkcc,
                base.calcium,
                "BEST_PER_ROOT_SOLVER_CROSSCHECK",
                PRODUCTION_BDF.label,
            )
    bdf_rows = _run_requests(tuple(solver_requests.values()), workers)
    radau_by_key = {
        (row["root_id"], row["nkcc_member_id"]): row for row in primary_rows
    }
    solver_rows: list[dict[str, Any]] = []
    for bdf in bdf_rows:
        key = (str(bdf["root_id"]), str(bdf["nkcc_member_id"]))
        radau = radau_by_key[key]
        solver_rows.append(
            {
                "root_id": key[0],
                "nkcc_member_id": key[1],
                "radau_success": radau.get("success", False),
                "bdf_success": bdf.get("success", False),
                "radau_numerical_gate_pass": radau.get("numerical_gate_pass", False),
                "bdf_numerical_gate_pass": bdf.get("numerical_gate_pass", False),
                "radau_flow_shape_gate_pass": radau.get(
                    "minute_flow_shape_gate_pass", False
                ),
                "bdf_flow_shape_gate_pass": bdf.get(
                    "minute_flow_shape_gate_pass", False
                ),
                "radau_minute_max_min_ratio": radau.get(
                    "minute_flow_max_min_ratio"
                ),
                "bdf_minute_max_min_ratio": bdf.get("minute_flow_max_min_ratio"),
                "abs_ratio_difference": abs(
                    float(radau.get("minute_flow_max_min_ratio", math.nan))
                    - float(bdf.get("minute_flow_max_min_ratio", math.nan))
                ),
                "radau_endpoint_flow_pL_s": radau.get("endpoint_flow_pL_s"),
                "bdf_endpoint_flow_pL_s": bdf.get("endpoint_flow_pL_s"),
                "abs_endpoint_flow_difference_pL_s": abs(
                    float(radau.get("endpoint_flow_pL_s", math.nan))
                    - float(bdf.get("endpoint_flow_pL_s", math.nan))
                ),
                "scope": (
                    "TWO_SCIPY_ALGORITHMS_NOT_INDEPENDENT_IMPLEMENTATIONS"
                ),
            }
        )

    profile_path = results_directory / PROFILE_FILENAME
    solver_path = results_directory / SOLVER_FILENAME
    _write_rows(profile_path, rows)
    _write_rows(solver_path, solver_rows)

    primary_eligible = [
        row for row in primary_rows if bool(row["source_eligible_for_g5_selection"])
    ]
    passing = [row for row in primary_eligible if bool(row["g5_repair_candidate_pass"])]
    best_nominal = _best_row(primary_eligible)
    best_all = _best_row(primary_rows)
    input_artifacts = (
        "wt_root_table.csv",
        "wt_root_states.csv",
        "wt_parameter_candidates.csv",
        "wt_selected_candidate.json",
    )
    summary = {
        "analysis": "WT_ONLY_G5_CROSS_ROOT_STIMULATED_NKCC1_REPAIR_PROFILE",
        "classification": (
            "SOURCE_ELIGIBLE_CROSS_ROOT_REPAIR_EXISTS"
            if passing
            else "NO_SOURCE_ELIGIBLE_CROSS_ROOT_REPAIR_PASSES_FIXED_WT_FLOW_SHAPE"
        ),
        "retained_root_count": len(roots),
        "retained_root_ids": [root.root_id for root in roots],
        "selected_reference_root": next(
            root.root_id for root in roots if root.selected_reference
        ),
        "primary_cross_root_member_count": len(nkcc_panel),
        "primary_simulation_count": len(primary_rows),
        "tier3_calcium_probe_count": len(calcium_rows),
        "solver_crosscheck_count": len(solver_rows),
        "fixed_protocol_duration_s": 600.0,
        "historical_calcium_anchor_uM": [0.058, 0.55],
        "tier3_calcium_probe_values_uM": [0.35, 0.75],
        "fixed_wt_graphical_envelope_uL_min": [9.0, 10.0],
        "flow_scale_status": (
            "TIER_2_WT_OBSERVATION_SCALE; IMPLIED_CELL_COUNT_WARNING_WITHOUT_"
            "PRIMARY_GEOMETRY_RANGE"
        ),
        "ae4_regulatory_screen_member": SELECTED_MEMBER_ID,
        "ae4_regulatory_ensemble_followup_required_if_candidate_passes": True,
        "nkcc_dynamic_kinetics_identified": False,
        "nkcc_kinetic_policy": (
            "N2_TAUS_ARE_PREDECLARED_UNMEASURED_SENSITIVITIES_NOT_POINT_FITS"
        ),
        "source_eligible_candidate_count": len(passing),
        "source_eligible_candidate_ids": [
            {"root_id": row["root_id"], "nkcc_member_id": row["nkcc_member_id"]}
            for row in passing
        ],
        "best_source_eligible_row": dict(best_nominal),
        "best_unrestricted_probe_row": dict(best_all),
        "all_integrations_success": all(bool(row["success"]) for row in rows),
        "all_numerical_gates_pass": all(
            bool(row["numerical_gate_pass"]) for row in rows
        ),
        "input_artifact_sha256": {
            name: sha256_file(results_directory / name) for name in input_artifacts
        },
        "declared_nkcc_panel_sha256": _sha256_payload(
            [asdict(item) for item in nkcc_panel]
        ),
        "declared_calcium_probe_sha256": _sha256_payload(
            [asdict(item) for item in declared_calcium_probes()]
        ),
        "profile_sha256": sha256_file(profile_path),
        "solver_crosscheck_sha256": sha256_file(solver_path),
        "firewall": "WT_CALIBRATION_AND_WT_DYNAMIC_EVIDENCE_ONLY",
    }
    summary_path = results_directory / SUMMARY_FILENAME
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    hash_rows = [
        {"artifact": path.name, "sha256": sha256_file(path), "scope": "WT_ONLY"}
        for path in (profile_path, solver_path, summary_path)
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
    result = run_cross_root_profile(args.results_directory, workers=args.workers)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = (
    "CalciumProbe",
    "FrozenRootCase",
    "NkccRepairSpec",
    "declared_calcium_probes",
    "declared_cross_root_nkcc_panel",
    "load_retained_wt_roots",
    "run_cross_root_profile",
)
