"""Independent positivity-preserving root and continuation utilities."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from ae4_mechanism_reconstruction.chassis import BASELINE_STATE, FixedChassis, Scenario, WT

from .adapter import StateResolvedAE4Adapter
from .cycles import CANDIDATE_DEFINITIONS, CycleDefinition, CycleParameters
from .fitting import TransporterFitResult
from .protocols import historical_calcium_step, immediate_pka_step


STATE_LOWER = np.array((60.0, 0.5, 2.0, 5.0, 50.0, 1.0, 1.0))
STATE_UPPER = np.array((180.0, 30.0, 100.0, 100.0, 180.0, 100.0, 100.0))
REDUCED_LOWER = np.array((60.0, 0.5, 5.0, 50.0, 1.0, 1.0, math.log(1e-6)))
REDUCED_UPPER = np.array((180.0, 30.0, 100.0, 180.0, 100.0, 100.0, math.log(1e-2)))
RAW_SCALES = np.array((1e-4, 1e-5, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3))
WT_CL_MM = 50.10
WT_CL_SEM_MM = 1.50
WT_PH = 6.91
WT_PH_SEM = 0.07


@dataclass(frozen=True)
class RootRecord:
    model_id: str
    gauge_id: str
    scenario: str
    capacity: float
    root_id: str
    solver_success: bool
    nfev: int
    cost: float
    max_abs_raw_rhs: float
    na_l_mM: float
    k_l_mM: float
    height_um: float
    na_i_mM: float
    k_i_mM: float
    cl_i_mM: float
    hco3_i_mM: float
    pH_i: float
    volume_pL: float
    cl_z_wt: float
    ph_z_wt: float
    positivity_pass: bool
    wt_physiology_gate_pass: bool
    alternate_root_count: int
    initial_start_count: int
    calibration_target_ids: str
    heldout_target_read: bool
    failure_reason: str

    @property
    def state(self) -> NDArray[np.float64]:
        return np.array(
            (
                self.na_l_mM,
                self.k_l_mM,
                self.height_um,
                self.na_i_mM,
                self.k_i_mM,
                self.cl_i_mM,
                self.hco3_i_mM,
            )
        )

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CapacityCalibration:
    model_id: str
    gauge_id: str
    capacity: float
    state: tuple[float, ...]
    solver_success: bool
    nfev: int
    cost: float
    active_bound_count: int
    max_abs_raw_rhs: float
    cl_z_wt: float
    pH_i: float
    volume_pL: float
    calibration_valid: bool
    calibration_target_ids: str
    heldout_target_read: bool
    failure_reason: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def parameters_from_fit(
    fit: TransporterFitResult, *, capacity: float
) -> CycleParameters:
    definition = CANDIDATE_DEFINITIONS[fit.model_id]
    overrides = {
        edge.edge_id: fit.common_barrier_multiplier
        for edge in definition.edges
        if edge.edge_id.endswith("common_cl_flip")
    }
    return CycleParameters(
        capacity=capacity,
        na_attempt_scale=fit.na_attempt_scale,
        k_attempt_scale=fit.k_attempt_scale,
        ec50_na_mM=fit.na_loaded_energy_lump_mM,
        ec50_k_mM=fit.k_loaded_energy_lump_mM,
        pka_barrier_fold=fit.pka_barrier_fold,
        edge_barrier_overrides=overrides,
    )


def build_fixed_chassis(
    definition: CycleDefinition,
    parameters: CycleParameters,
) -> FixedChassis:
    return FixedChassis(
        StateResolvedAE4Adapter(definition, parameters),
        calcium_protocol=historical_calcium_step,
        pka_protocol=(immediate_pka_step if definition.pka_switch_edge_ids else (lambda _t: 0.0)),
        pka_discontinuities=((100.0,) if definition.pka_switch_edge_ids else ()),
    )


def _diagnostics_or_none(
    chassis: FixedChassis, state: NDArray[np.float64], scenario: Scenario
) -> tuple[NDArray[np.float64], dict[str, object]] | None:
    try:
        raw, diagnostics = chassis.evaluate(
            state, calcium=0.05, pka_activation=0.0, scenario=scenario
        )
    except (ValueError, FloatingPointError, OverflowError):
        return None
    h_i = float(diagnostics["h_i"])
    if h_i <= 0.0 or not np.all(np.isfinite(raw)):
        return None
    return np.asarray(raw), dict(diagnostics)


def _state_to_reduced(state: Sequence[float], chassis: FixedChassis) -> np.ndarray:
    y = np.asarray(state, dtype=float)
    p = chassis.parameters
    h_i = y[6] + y[5] + p.impermeant_charge_amount / y[2] - y[3] - y[4]
    if h_i <= 0.0 or not math.isfinite(h_i):
        h_i = 10.0 ** (3.0 - WT_PH)
    return np.array((y[0], y[1], y[3], y[4], y[5], y[6], math.log(h_i)))


def _reduced_to_state(reduced: Sequence[float], chassis: FixedChassis) -> np.ndarray:
    na_l, k_l, na_i, k_i, cl_i, hco3_i, log_h_i = map(float, reduced)
    h_i = math.exp(log_h_i)
    denominator = na_i + k_i + h_i - cl_i - hco3_i
    height = chassis.parameters.impermeant_charge_amount / denominator
    return np.array((na_l, k_l, height, na_i, k_i, cl_i, hco3_i))


def _starts(
    seed: int,
    count: int,
    chassis: FixedChassis,
    preferred: Sequence[float] | None = None,
) -> tuple[np.ndarray, ...]:
    rows = [_state_to_reduced(BASELINE_STATE, chassis)]
    if preferred is not None:
        candidate = np.asarray(preferred, dtype=float)
        if candidate.shape == (7,):
            reduced = _state_to_reduced(candidate, chassis)
            if not np.allclose(reduced, rows[0]):
                rows.append(reduced)
    generator = np.random.default_rng(seed)
    while len(rows) < count:
        candidate = generator.uniform(REDUCED_LOWER, REDUCED_UPPER)
        state = _reduced_to_state(candidate, chassis)
        if np.all(np.isfinite(state)) and 2.0 < state[2] < 200.0:
            rows.append(candidate)
    return tuple(rows)


def _root_residual(
    reduced: NDArray[np.float64], chassis: FixedChassis, scenario: Scenario
) -> NDArray[np.float64]:
    state = _reduced_to_state(reduced, chassis)
    if not np.all(np.isfinite(state)) or not 2.0 < state[2] < 200.0:
        return np.full(7, 1e6)
    evaluated = _diagnostics_or_none(chassis, state, scenario)
    if evaluated is None:
        return np.full(7, 1e6)
    raw, _ = evaluated
    return raw / RAW_SCALES


def solve_fixed_capacity_roots(
    chassis: FixedChassis,
    *,
    model_id: str,
    gauge_id: str,
    capacity: float,
    scenario: Scenario = WT,
    seed: int = 20260827,
    start_count: int = 12,
    preferred: Sequence[float] | None = None,
) -> tuple[RootRecord, ...]:
    fits = []
    for start in _starts(seed, start_count, chassis, preferred):
        fit = least_squares(
            _root_residual,
            np.clip(start, REDUCED_LOWER, REDUCED_UPPER),
            args=(chassis, scenario),
            bounds=(REDUCED_LOWER, REDUCED_UPPER),
            x_scale="jac",
            max_nfev=2000,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        state = _reduced_to_state(fit.x, chassis)
        evaluated = _diagnostics_or_none(chassis, state, scenario)
        if evaluated is None:
            continue
        raw, diagnostics = evaluated
        if fit.success and np.max(np.abs(raw)) <= 1e-8:
            fits.append((fit, state, diagnostics, float(np.max(np.abs(raw)))))
    distinct = []
    for row in sorted(fits, key=lambda item: item[0].cost):
        if not any(np.allclose(row[1], known[1], rtol=1e-6, atol=1e-5) for known in distinct):
            distinct.append(row)
    records = []
    for index, (fit, state, diagnostics, raw_error) in enumerate(distinct):
        ph = float(diagnostics["pH_i"])
        volume = float(diagnostics["cell_volume_pL"])
        cl_z = (state[5] - WT_CL_MM) / WT_CL_SEM_MM
        ph_z = (ph - WT_PH) / WT_PH_SEM
        physiological = bool(
            abs(cl_z) <= 2.0
            and abs(ph_z) <= 2.0
            and 5.0 <= state[3] <= 80.0
            and 70.0 <= state[4] <= 170.0
            and 0.5 <= volume <= 5.0
        )
        reasons = []
        if abs(cl_z) > 2.0:
            reasons.append(f"WT Cl {cl_z:.3f} SEM from E15-05")
        if abs(ph_z) > 2.0:
            reasons.append(f"WT pH {ph_z:.3f} SEM from E15-07")
        if not 0.5 <= volume <= 5.0:
            reasons.append("volume outside broad declared 0.5-5 pL domain")
        records.append(
            RootRecord(
                model_id=model_id,
                gauge_id=gauge_id,
                scenario=scenario.name,
                capacity=capacity,
                root_id=f"{model_id}:{gauge_id}:{scenario.name}:r{index}",
                solver_success=True,
                nfev=int(fit.nfev),
                cost=float(fit.cost),
                max_abs_raw_rhs=raw_error,
                na_l_mM=float(state[0]),
                k_l_mM=float(state[1]),
                height_um=float(state[2]),
                na_i_mM=float(state[3]),
                k_i_mM=float(state[4]),
                cl_i_mM=float(state[5]),
                hco3_i_mM=float(state[6]),
                pH_i=ph,
                volume_pL=volume,
                cl_z_wt=float(cl_z),
                ph_z_wt=float(ph_z),
                positivity_pass=bool(np.all(state > 0.0)),
                wt_physiology_gate_pass=physiological,
                alternate_root_count=len(distinct),
                initial_start_count=start_count,
                calibration_target_ids="E15-05;E15-07_WT",
                heldout_target_read=False,
                failure_reason="; ".join(reasons),
            )
        )
    return tuple(records)


def calibrate_wt_capacity(
    definition: CycleDefinition,
    fit_row: TransporterFitResult,
    *,
    seed: int = 20260827,
) -> CapacityCalibration:
    """Estimate one chassis capacity from WT balance + Cl/pH only."""

    reference_chassis = build_fixed_chassis(definition, parameters_from_fit(fit_row, capacity=1.0))
    lower = np.concatenate((REDUCED_LOWER, np.array((-12.0,))))
    upper = np.concatenate((REDUCED_UPPER, np.array((16.0,))))

    def residual(coordinates: NDArray[np.float64]) -> NDArray[np.float64]:
        state = _reduced_to_state(coordinates[:7], reference_chassis)
        if not np.all(np.isfinite(state)) or not 2.0 < state[2] < 200.0:
            return np.full(8, 1e6)
        capacity = math.exp(float(coordinates[7]))
        parameters = parameters_from_fit(fit_row, capacity=capacity)
        chassis = build_fixed_chassis(definition, parameters)
        evaluated = _diagnostics_or_none(chassis, state, WT)
        if evaluated is None:
            return np.full(8, 1e6)
        raw, diagnostics = evaluated
        # One fitted capacity is identified by one WT observation (E15-05).
        # WT pH remains an independent acceptance gate, not a second target
        # squeezed into an overdetermined compromise.
        return np.concatenate(
            (raw / RAW_SCALES, np.array(((state[5] - WT_CL_MM) / WT_CL_SEM_MM,)))
        )

    generator = np.random.default_rng(seed)
    starts = []
    for log_capacity in np.linspace(-8.0, 10.0, 6):
        starts.append(
            np.concatenate(
                (_state_to_reduced(BASELINE_STATE, reference_chassis), np.array((log_capacity,)))
            )
        )
    for _ in range(2):
        starts.append(
            np.concatenate(
                (
                    generator.uniform(REDUCED_LOWER, REDUCED_UPPER),
                    np.array((generator.uniform(-10.0, 12.0),)),
                )
            )
        )
    fits = [
        least_squares(
            residual,
            np.clip(start, lower, upper),
            bounds=(lower, upper),
            x_scale="jac",
            max_nfev=2000,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        for start in starts
    ]
    best = min(fits, key=lambda result: result.cost)
    capacity = math.exp(float(best.x[7]))
    state = _reduced_to_state(best.x[:7], reference_chassis)
    parameters = parameters_from_fit(fit_row, capacity=capacity)
    chassis = build_fixed_chassis(definition, parameters)
    evaluated = _diagnostics_or_none(chassis, state, WT)
    if evaluated is None:
        raw_error = ph = volume = math.nan
    else:
        raw, diagnostics = evaluated
        raw_error = float(np.max(np.abs(raw)))
        ph = float(diagnostics["pH_i"])
        volume = float(diagnostics["cell_volume_pL"])
    cl_z = float((state[5] - WT_CL_MM) / WT_CL_SEM_MM)
    active = int(np.count_nonzero(best.active_mask))
    valid = bool(
        best.success
        and active == 0
        and math.isfinite(raw_error)
        and raw_error <= 1e-8
        and abs(cl_z) <= 0.25
        and 2.0 < state[2] < 200.0
    )
    failures = []
    if not best.success:
        failures.append(str(best.message))
    if active:
        failures.append(f"{active} active optimization bounds")
    if not math.isfinite(raw_error) or raw_error > 1e-8:
        failures.append(f"raw closure {raw_error:.6g} exceeds 1e-8")
    if abs(cl_z) > 0.25:
        failures.append(f"WT Cl calibration residual {cl_z:.3f} SEM")
    if not 2.0 < state[2] < 200.0:
        failures.append("height outside open 2-200 um domain")
    return CapacityCalibration(
        model_id=definition.model_id,
        gauge_id=fit_row.gauge_id,
        capacity=capacity,
        state=tuple(float(value) for value in state),
        solver_success=bool(best.success),
        nfev=int(best.nfev),
        cost=float(best.cost),
        active_bound_count=active,
        max_abs_raw_rhs=raw_error,
        cl_z_wt=cl_z,
        pH_i=ph,
        volume_pL=volume,
        calibration_valid=valid,
        calibration_target_ids="E15-05",
        heldout_target_read=False,
        failure_reason="; ".join(failures),
    )
