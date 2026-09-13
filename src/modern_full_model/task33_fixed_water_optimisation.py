"""Task 33: one bounded R09 WT optimisation followed by two held out AE4 tests.

The immutable numerical contract is stored at
``analysis/33_fixed_water_nkcc70_ae4_optimisation/optimisation_contract.md``.
This module deliberately has two execution stages so the accepted WT payload
can be published and hashed before either held out genotype is evaluated.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, fields, is_dataclass, replace
from enum import Enum
import hashlib
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import Bounds, least_squares, minimize

from .ae4_routing_only import MODEL_ID as AE4_ROUTING_MODEL_ID
from .ae4_routing_only import routing_only_adapter
from .camp_pka import AE4Construct, R1EffectiveActivation, RegulatoryGain
from .membranes import current_to_fmol_s
from .model import Genotype, ModernFullModel, WT
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .nhe1_cha2009 import NHE1_CHA_2009, published_cha_kinetics
from .states import CORE_STATE_NAMES
from .task30_nhe1_repair import (
    AMOUNT_ROWS,
    BOUNDARY_RELATIVE_TOLERANCE,
    CELL_VOLUME_UPPER_PL,
    CHARGE_TOLERANCE_FMOL,
    COORDINATE_LOWER,
    COORDINATE_NAMES,
    COORDINATE_UPPER,
    CURRENT_TOLERANCE_A,
    INDEPENDENT_RHS_SCALES,
    INDEPENDENT_ROWS,
    JACOBIAN_RANK_RELATIVE_TOLERANCE,
    OMITTED_CHARGE_ROWS,
    OMITTED_ROW_TOLERANCE,
    ORDINARY_NA_UPPER_MM,
    REGULATORY_TOLERANCE_S_INV,
    ROOT_SCALED_TOLERANCE,
    VOLUME_ROWS,
    _conservation_audit,
    _flux_row,
    _protocol,
    _rest_observables,
    load_background,
    state_from_coordinates,
)
from .task31_nhe1_repair import _finite_capacity_audit, task31_parameters
from .validation import PRODUCTION_RADAU, sha256_file, sha256_object


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analysis/33_fixed_water_nkcc70_ae4_optimisation"
RESULTS = ROOT / "results/33_fixed_water_nkcc70_ae4_optimisation"
CONTRACT = ANALYSIS / "optimisation_contract.md"
CONTRACT_SHA256 = "aef582891dfff7f368fc36ee7768dd9456e43da29529a04e64c6759c7642526c"
PREPARED_HEAD = "3a8560d39eb6771fce76895e285c682d5d3333ac"
CONTRACT_COMMIT = "1d2c1429cd410efb9764b78afc42a2eba441b5df"
BRANCH = "codex/task-33-fixed-water-nkcc70-ae4-optimisation"

NHE1_CARRIER_AMOUNT_FMOL = 2.339370005697548e-05
REST_Q_OUT_TARGET_PL_S = 0.0010757073853493032
REST_Q_OUT_ABS_TOLERANCE_PL_S = 1.0757073853493032e-09
REST_Q_OUT_REL_TOLERANCE = 1.0e-6
STIMULATED_CALCIUM_UM = 0.25
ACTIVATION_FACTOR = 1.10

PARAMETER_NAMES = (
    "nkcc1_capacity",
    "ae4_carrier_amount",
    "ae2_capacity",
    "nak_pump_capacity",
    "k_conductance",
    "cacc_conductance",
)
REFERENCE_VALUES = np.asarray((0.32, 0.1501264265070657, 0.005, 0.08, 1.4e-08, 3.14e-08))
MULTIPLIER_LOWER = np.asarray((0.25, 0.1, 0.25, 0.5, 0.5, 0.5))
MULTIPLIER_UPPER = np.asarray((2.0, 50.0, 4.0, 2.0, 2.0, 2.0))
LOG_LOWER = np.log(MULTIPLIER_LOWER)
LOG_UPPER = np.log(MULTIPLIER_UPPER)
INITIAL_MULTIPLIERS = np.asarray((0.7134869177422251, 15.870644383966862, 1.0, 1.0, 1.0, 1.0))
INITIAL_LOG_MULTIPLIERS = np.log(INITIAL_MULTIPLIERS)

OPTIMISER_VECTOR_CEILING = 19
HARD_WT_VECTOR_CEILING = 20
STATIONARY_SOLVE_CEILING = 30
STATIONARY_RESIDUAL_CEILING = 40_000
PER_SOLVE_RESIDUAL_CEILING = 1_600
PREFREEZE_DYNAMIC_CEILING = 3
POSTFREEZE_DYNAMIC_CEILING = 2
WT_WALL_SECONDS_CEILING = 1_200.0

AE4_FIVE_PERCENT = Genotype("AE4_5_PERCENT", ae4_expression=0.05)
AE4_ZERO = Genotype("AE4_ZERO", ae4_expression=0.0)
HELD_OUT_GENOTYPES = (AE4_FIVE_PERCENT, AE4_ZERO)


class Task33Stop(RuntimeError):
    """A global Task 33 hard stop."""


class LocalStationaryStop(RuntimeError):
    """The immutable per solve residual ceiling was reached."""


def _jsonable(value: Any) -> Any:
    if is_dataclass(value):
        return {field.name: _jsonable(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    names: list[str] = []
    for row in rows:
        for key in row:
            if key not in names:
                names.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows([{key: _jsonable(value) for key, value in row.items()} for row in rows])


@dataclass
class Budget:
    distinct_wt_parameter_vectors: int = 0
    stationary_solver_calls: int = 0
    stationary_residual_evaluations: int = 0
    prefreeze_wt_dynamic_integrations: int = 0
    postfreeze_genotype_dynamic_integrations: int = 0
    numerical_wt_wall_seconds: float = 0.0
    optimiser_calls: int = 0
    deterministic_restarts: int = 0
    wt_accepted: bool = False

    def __post_init__(self) -> None:
        self._clock = time.perf_counter()
        self._timing_wt = not self.wt_accepted

    @classmethod
    def load(cls) -> "Budget":
        payload = json.loads((RESULTS / "budget.json").read_text(encoding="utf-8"))
        names = {field.name for field in fields(cls)}
        return cls(**{key: payload[key] for key in names if key in payload})

    def tick_wt(self) -> None:
        now = time.perf_counter()
        if self._timing_wt:
            self.numerical_wt_wall_seconds += now - self._clock
        self._clock = now
        if not self.wt_accepted and self.numerical_wt_wall_seconds >= WT_WALL_SECONDS_CEILING:
            raise Task33Stop("1200 s WT numerical wall time ceiling reached")

    def start_wt_vector(self) -> None:
        self.tick_wt()
        if self.distinct_wt_parameter_vectors >= HARD_WT_VECTOR_CEILING:
            raise Task33Stop("20 distinct WT parameter vector ceiling reached")
        self.distinct_wt_parameter_vectors += 1

    def start_stationary(self, *, prefreeze: bool) -> None:
        if prefreeze:
            self.tick_wt()
        if self.stationary_solver_calls >= STATIONARY_SOLVE_CEILING:
            raise Task33Stop("30 stationary solve ceiling reached")
        self.stationary_solver_calls += 1

    def count_stationary_residual(self, *, prefreeze: bool) -> None:
        if prefreeze:
            self.tick_wt()
        if self.stationary_residual_evaluations >= STATIONARY_RESIDUAL_CEILING:
            raise Task33Stop("40000 stationary residual evaluation ceiling reached")
        self.stationary_residual_evaluations += 1

    def start_dynamic(self, *, prefreeze: bool) -> None:
        if prefreeze:
            self.tick_wt()
            if self.prefreeze_wt_dynamic_integrations >= PREFREEZE_DYNAMIC_CEILING:
                raise Task33Stop("3 pre-freeze WT dynamic integration ceiling reached")
            self.prefreeze_wt_dynamic_integrations += 1
        else:
            if self.postfreeze_genotype_dynamic_integrations >= POSTFREEZE_DYNAMIC_CEILING:
                raise Task33Stop("2 post-freeze genotype dynamic integration ceiling reached")
            self.postfreeze_genotype_dynamic_integrations += 1

    def accept_wt(self) -> None:
        self.tick_wt()
        self.wt_accepted = True
        self._timing_wt = False

    def payload(self) -> dict[str, Any]:
        if self._timing_wt:
            self.tick_wt()
        payload = asdict(self)
        payload.update({
            "scientific_workers": 1,
            "blas_threads": 1,
            "optimiser_vector_ceiling": OPTIMISER_VECTOR_CEILING,
            "hard_wt_parameter_vector_ceiling": HARD_WT_VECTOR_CEILING,
            "stationary_solver_call_ceiling": STATIONARY_SOLVE_CEILING,
            "stationary_residual_evaluation_ceiling": STATIONARY_RESIDUAL_CEILING,
            "per_solve_residual_ceiling": PER_SOLVE_RESIDUAL_CEILING,
            "prefreeze_dynamic_ceiling": PREFREEZE_DYNAMIC_CEILING,
            "postfreeze_dynamic_ceiling": POSTFREEZE_DYNAMIC_CEILING,
            "wt_wall_seconds_ceiling": WT_WALL_SECONDS_CEILING,
        })
        payload["limits_respected"] = bool(
            self.distinct_wt_parameter_vectors <= HARD_WT_VECTOR_CEILING
            and self.stationary_solver_calls <= STATIONARY_SOLVE_CEILING
            and self.stationary_residual_evaluations <= STATIONARY_RESIDUAL_CEILING
            and self.prefreeze_wt_dynamic_integrations <= PREFREEZE_DYNAMIC_CEILING
            and self.postfreeze_genotype_dynamic_integrations <= POSTFREEZE_DYNAMIC_CEILING
            and self.numerical_wt_wall_seconds <= WT_WALL_SECONDS_CEILING
        )
        return payload

    def save(self) -> None:
        _write_json(RESULTS / "budget.json", self.payload())


def _check_environment() -> None:
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        if os.environ.get(name) != "1":
            raise RuntimeError(name + " must equal 1 before Python starts")
    if sha256_file(CONTRACT) != CONTRACT_SHA256:
        raise RuntimeError("immutable optimisation contract hash mismatch")


def task31_reference_wt() -> dict[str, Any]:
    payload = json.loads(
        (ROOT / "results/31_nhe1_mechanistic_repair/calibration.json").read_text(encoding="utf-8")
    )
    matches = [
        item for item in payload["evaluations"]
        if item["background"] == "R09"
        and item["genotype"] == "WT"
        and item["admissible"]
        and item["carrier_amount_fmol"] == NHE1_CARRIER_AMOUNT_FMOL
    ]
    if len(matches) != 1:
        raise RuntimeError("expected exactly one accepted Task 31 R09 WT reference")
    reference = matches[0]
    if reference["fluxes"]["water_lumen_outflow_pL_s"] != REST_Q_OUT_TARGET_PL_S:
        raise RuntimeError("Task 31 WT outflow reference differs from fixed Task 33 target")
    return reference


def multipliers_from_log(log_multipliers: ArrayLike) -> NDArray[np.float64]:
    z = np.asarray(log_multipliers, dtype=float)
    if z.shape != (len(PARAMETER_NAMES),) or np.any(~np.isfinite(z)):
        raise ValueError("log multiplier vector must be finite and length six")
    if np.any(z < LOG_LOWER - 1e-14) or np.any(z > LOG_UPPER + 1e-14):
        raise ValueError("log multiplier vector lies outside immutable bounds")
    return np.exp(z)


def absolute_free_parameters(multipliers: ArrayLike) -> dict[str, float]:
    values = REFERENCE_VALUES * np.asarray(multipliers, dtype=float)
    return {name: float(value) for name, value in zip(PARAMETER_NAMES, values)}


def build_model(multipliers: ArrayLike, *, stimulated: bool) -> tuple[Any, Any]:
    m = np.asarray(multipliers, dtype=float)
    if m.shape != (6,) or np.any(m < MULTIPLIER_LOWER) or np.any(m > MULTIPLIER_UPPER):
        raise ValueError("Task 33 parameter multipliers violate immutable bounds")
    background = load_background("R09")
    base = task31_parameters(background, carrier_amount_fmol=NHE1_CARRIER_AMOUNT_FMOL)
    parameters = replace(
        base,
        homeostasis=replace(
            base.homeostasis,
            nkcc1_capacity_fmol_s=float(REFERENCE_VALUES[0] * m[0]),
            ae2_capacity_fmol_s=float(REFERENCE_VALUES[2] * m[2]),
        ),
        membranes=replace(
            base.membranes,
            nak_capacity_fmol_s=float(REFERENCE_VALUES[3] * m[3]),
            g_k_total_S=float(REFERENCE_VALUES[4] * m[4]),
            g_cl_apical_S=float(REFERENCE_VALUES[5] * m[5]),
        ),
    )
    ae4_parameters = replace(
        background.ae4_parameters,
        carrier_amount_fmol=float(REFERENCE_VALUES[1] * m[1]),
    )
    regulation = R1EffectiveActivation(
        tau_activation_s=30.0,
        gain=RegulatoryGain(
            basal_capacity_multiplier=1.0,
            fully_activated_increment=0.25,
            coupling_scale=1.0,
        ),
        construct=AE4Construct.WT,
    )
    base_model = ModernFullModel(
        parameters,
        stimulus=_protocol(stimulated=stimulated),
        regulatory_model=regulation,
        ae4_parameters=ae4_parameters,
        ae4_evaluator=routing_only_adapter,
    )
    model = attach_stimulated_nkcc1(
        base_model,
        N1AlgebraicNkcc1(
            fully_activated_multiplier=1.75,
            resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.10,
        ),
    )
    if model.parameters.homeostasis.nhe1_model != NHE1_CHA_2009:
        raise RuntimeError("Task 31 NHE1 model changed")
    if model.parameters.homeostasis.nhe1_cha_carrier_amount_fmol != NHE1_CARRIER_AMOUNT_FMOL:
        raise RuntimeError("Task 31 NHE1 carrier amount changed")
    return model, background


def changed_parameter_paths(multipliers: ArrayLike) -> tuple[str, ...]:
    """Return differences from Task 31 to support the mutation whitelist test."""
    model, background = build_model(multipliers, stimulated=False)
    reference = task31_parameters(background, carrier_amount_fmol=NHE1_CARRIER_AMOUNT_FMOL)
    changed: list[str] = []

    def walk(prefix: str, left: Any, right: Any) -> None:
        if is_dataclass(left) and is_dataclass(right):
            for field in fields(left):
                walk(prefix + field.name + ".", getattr(left, field.name), getattr(right, field.name))
        elif left != right:
            changed.append(prefix[:-1])

    walk("parameters.", reference, model.parameters)
    if background.ae4_parameters.carrier_amount_fmol != model.ae4_parameters.carrier_amount_fmol:
        changed.append("ae4_parameters.carrier_amount_fmol")
    for name in ("common_cl_attempt_rate_s", "na_loaded_attempt_rate_s", "k_loaded_attempt_rate_s", "cooperative_gate"):
        if getattr(background.ae4_parameters, name) != getattr(model.ae4_parameters, name):
            changed.append("ae4_parameters." + name)
    return tuple(sorted(changed))


def _extended_flux_row(evaluation: Any, model: Any) -> dict[str, Any]:
    row = dict(_flux_row(evaluation, model))
    currents = evaluation.diagnostics.membranes.currents_A
    f = model.parameters.constants.faraday_C_mol
    para_cl = current_to_fmol_s(float(currents["para_cl"]), valence=-1, faraday_C_mol=f)
    lumen = evaluation.diagnostics.observables.lumen_concentrations_mM
    q_out = float(evaluation.diagnostics.water.lumen_outflow_pL_s)
    nkcc = max(float(row["nkcc1_cl_in_fmol_s"]), 0.0)
    ae4 = max(float(row["ae4_cl_cell_fmol_s"]), 0.0)
    ae2 = max(float(row["ae2_cl_in_hco3_out_fmol_s"]), 0.0)
    pool = nkcc + ae4 + ae2
    row.update({
        "paracellular_cl_lumen_to_bath_fmol_s": float(para_cl),
        "chloride_outflow_fmol_s": q_out * float(lumen["cl"]),
        "positive_chloride_loading_pool_fmol_s": pool,
        "nkcc1_positive_pool_share": nkcc / pool if pool > 0.0 else math.nan,
        "ae4_positive_pool_share": ae4 / pool if pool > 0.0 else math.nan,
        "ae2_positive_pool_share": ae2 / pool if pool > 0.0 else math.nan,
    })
    return row


def solve_stationary(
    model: Any,
    *,
    genotype: Genotype,
    start_coordinates: ArrayLike,
    multipliers: ArrayLike,
    label: str,
    budget: Budget,
    prefreeze: bool,
) -> dict[str, Any]:
    """One bounded stationary solve with exact actual residual accounting."""
    budget.start_stationary(prefreeze=prefreeze)
    started = time.perf_counter()
    start = np.asarray(start_coordinates, dtype=float)
    if start.shape != (10,) or np.any(start <= COORDINATE_LOWER) or np.any(start >= COORDINATE_UPPER):
        raise ValueError(label + ": stationary start is outside inherited bounds")
    local_count = 0
    best_x = start.copy()
    best_cost = math.inf

    def counted(candidate: NDArray[np.float64], *, audit: bool = False) -> NDArray[np.float64]:
        nonlocal local_count, best_x, best_cost
        reserved_limit = PER_SOLVE_RESIDUAL_CEILING if audit else PER_SOLVE_RESIDUAL_CEILING - 13
        if local_count >= reserved_limit:
            raise LocalStationaryStop("1600 actual residual evaluations per solve")
        budget.count_stationary_residual(prefreeze=prefreeze)
        local_count += 1
        try:
            state = state_from_coordinates(model, candidate)
            raw = model.rhs(0.0, state, genotype=genotype)
            values = np.asarray(raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES, dtype=float)
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            values = np.full(10, 1.0e6)
        cost = float(values @ values) / 2.0
        if cost < best_cost:
            best_x = np.asarray(candidate, dtype=float).copy()
            best_cost = cost
        return values

    try:
        try:
            fit = least_squares(
                counted,
                start,
                bounds=(COORDINATE_LOWER, COORDINATE_UPPER),
                max_nfev=PER_SOLVE_RESIDUAL_CEILING,
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
                x_scale="jac",
            )
        except LocalStationaryStop as exc:
            final_x = best_x.copy()
            f0 = counted(final_x, audit=True)
            jac = np.empty((10, 10))
            for column in range(10):
                dx = math.sqrt(np.finfo(float).eps) * max(abs(final_x[column]), 1.0)
                if final_x[column] + dx >= COORDINATE_UPPER[column]:
                    dx = -dx
                shifted = final_x.copy()
                shifted[column] += dx
                jac[:, column] = (counted(shifted, audit=True) - f0) / dx
            fit = SimpleNamespace(
                x=final_x,
                jac=jac,
                success=False,
                message=str(exc),
                nfev=None,
                cost=float(f0 @ f0) / 2.0,
                optimality=float(np.max(np.abs(jac.T @ f0))),
            )
        coordinates = np.asarray(fit.x, dtype=float)
        state = state_from_coordinates(model, coordinates)
        budget.count_stationary_residual(prefreeze=prefreeze)
        local_count += 1
        evaluation = model.evaluate(0.0, state, genotype=genotype)
        budget.count_stationary_residual(prefreeze=prefreeze)
        local_count += 1
        later = model.evaluate(600.0, state, genotype=genotype)
    finally:
        if prefreeze:
            budget.tick_wt()
        _ = time.perf_counter() - started

    raw = np.asarray(evaluation.rhs, dtype=float)
    scaled = raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES
    span = COORDINATE_UPPER - COORDINATE_LOWER
    jacobian = np.asarray(fit.jac, dtype=float) * span[None, :]
    singular = np.linalg.svd(jacobian, compute_uv=False)
    threshold = JACOBIAN_RANK_RELATIVE_TOLERANCE * singular[0] if singular.size else math.inf
    rank = int(np.sum(singular > threshold))
    relative_boundary_distance = np.minimum(
        coordinates - COORDINATE_LOWER, COORDINATE_UPPER - coordinates
    ) / span
    boundary_hits = [
        name for name, distance in zip(COORDINATE_NAMES, relative_boundary_distance)
        if distance <= BOUNDARY_RELATIVE_TOLERANCE
    ]
    charge_max = max(abs(float(value)) for value in evaluation.diagnostics.state_charge_fmol.values())
    current_max = max(abs(float(value)) for value in evaluation.diagnostics.membranes.current_residuals_A.values())
    regulatory_max = float(np.max(np.abs(raw[12:]))) if raw.size > 12 else 0.0
    conservation, conservation_ratio, conservation_pass = _conservation_audit(evaluation)
    capacity = _finite_capacity_audit(evaluation, model)
    observables = _rest_observables(evaluation, state)
    fluxes = _extended_flux_row(evaluation, model)
    full_rhs_pass = bool(
        np.max(np.abs(scaled)) <= ROOT_SCALED_TOLERANCE
        and np.max(np.abs(raw[OMITTED_CHARGE_ROWS])) <= OMITTED_ROW_TOLERANCE
        and regulatory_max <= REGULATORY_TOLERANCE_S_INV
    )
    numerical_pass = bool(
        fit.success
        and full_rhs_pass
        and charge_max <= CHARGE_TOLERANCE_FMOL
        and current_max <= CURRENT_TOLERANCE_A
        and conservation_pass
        and rank == 10
        and not boundary_hits
        and np.all(state[:12] > 0.0)
        and np.array_equal(raw, later.rhs)
    )
    concentration_keys = (
        "na_i_mM", "k_i_mM", "cl_i_mM", "tic_i_mM", "hco3_i_mM",
        "na_l_mM", "k_l_mM", "cl_l_mM", "tic_l_mM", "hco3_l_mM",
    )
    positive_concentrations = bool(
        all(math.isfinite(observables[key]) and observables[key] > 0.0 for key in concentration_keys)
    )
    general_physical_pass = bool(
        numerical_pass
        and capacity["pass"]
        and positive_concentrations
        and observables["na_i_mM"] <= ORDINARY_NA_UPPER_MM
        and 60.0 <= observables["k_i_mM"] <= 210.0
        and observables["tic_i_mM"] < 100.0
        and observables["hco3_i_mM"] < 100.0
        and 0.0 < observables["volume_i_pL"] < CELL_VOLUME_UPPER_PL
        and observables["volume_l_pL"] > 0.0
    )
    exact_ae4_zero = bool(
        genotype.ae4_expression != 0.0
        or all(float(fluxes[key]) == 0.0 for key in (
            "ae4_na_cell_fmol_s", "ae4_k_cell_fmol_s", "ae4_cl_cell_fmol_s", "ae4_hco3_cell_fmol_s"
        ))
    )
    return {
        "background": "R09",
        "genotype": genotype.name,
        "ae4_expression": float(genotype.ae4_expression),
        "label": label,
        "parameter_multipliers": {name: float(value) for name, value in zip(PARAMETER_NAMES, multipliers)},
        "absolute_free_parameters": absolute_free_parameters(multipliers),
        "nhe1_model": model.parameters.homeostasis.nhe1_model,
        "nhe1_carrier_amount_fmol": model.parameters.homeostasis.nhe1_cha_carrier_amount_fmol,
        "optimizer_success": bool(fit.success),
        "optimizer_message": str(fit.message),
        "optimizer_nfev": int(fit.nfev) if fit.nfev is not None else None,
        "actual_residual_evaluations": local_count,
        "cost": float(fit.cost),
        "optimality": float(fit.optimality),
        "coordinates": coordinates.tolist(),
        "core_state": state[:12].tolist(),
        "core_state_sha256": sha256_object(state[:12].tolist()),
        "raw_rhs": raw.tolist(),
        "scaled_independent_rhs": scaled.tolist(),
        "max_abs_scaled_independent_rhs": float(np.max(np.abs(scaled))),
        "max_abs_amount_rhs_fmol_s": float(np.max(np.abs(raw[AMOUNT_ROWS]))),
        "max_abs_volume_rhs_pL_s": float(np.max(np.abs(raw[VOLUME_ROWS]))),
        "max_abs_omitted_rhs_fmol_s": float(np.max(np.abs(raw[OMITTED_CHARGE_ROWS]))),
        "max_abs_regulatory_rhs_s_inv": regulatory_max,
        "max_abs_bulk_charge_fmol": charge_max,
        "max_abs_current_residual_A": current_max,
        "max_dimensionless_conservation_ratio": conservation_ratio,
        "conservation_residuals": conservation,
        "normalised_jacobian_singular_values": singular.tolist(),
        "independent_jacobian_rank": rank,
        "boundary_hits": boundary_hits,
        "rest_protocol_time_invariant_rhs": bool(np.array_equal(raw, later.rhs)),
        "positivity_pass": bool(np.all(state[:12] > 0.0)),
        "positive_concentrations_pass": positive_concentrations,
        "full_rhs_pass": full_rhs_pass,
        "conservation_pass": conservation_pass,
        "finite_capacity": capacity,
        "numerical_pass": numerical_pass,
        "general_physical_pass": general_physical_pass,
        "exact_ae4_zero": exact_ae4_zero,
        "admissible": bool(general_physical_pass and exact_ae4_zero),
        "observables": observables,
        "fluxes": fluxes,
        "active_parameter_sha256": sha256_object(model.parameters),
        "active_ae4_parameter_sha256": sha256_object(model.ae4_parameters),
    }


def _soft_objective(solution: Mapping[str, Any], multipliers: NDArray[np.float64], reference: Mapping[str, Any]) -> float:
    logs = np.log(multipliers)
    parameter_distances = np.where(
        logs >= 0.0,
        logs / np.log(MULTIPLIER_UPPER),
        logs / np.abs(np.log(MULTIPLIER_LOWER)),
    )
    coordinates = np.asarray(solution["coordinates"], dtype=float)
    reference_coordinates = np.asarray(reference["coordinates"], dtype=float)
    state_distances = (coordinates - reference_coordinates) / (COORDINATE_UPPER - COORDINATE_LOWER)
    chloride_distance = (float(solution["observables"]["cl_i_mM"]) - 60.30496692587399) / 20.0
    return float(parameter_distances @ parameter_distances + 0.1 * (state_distances @ state_distances) + 0.1 * chloride_distance**2)


def _constraint_margins(solution: Mapping[str, Any]) -> dict[str, float]:
    o = solution["observables"]
    f = solution["fluxes"]
    capacity_max = max(float(value) for value in solution["finite_capacity"]["utilisation"].values())
    q_error = abs(float(f["water_lumen_outflow_pL_s"]) - REST_Q_OUT_TARGET_PL_S)
    share = float(f["nkcc1_positive_pool_share"])
    ae4_share = float(f["ae4_positive_pool_share"])
    concentration_keys = (
        "na_i_mM", "k_i_mM", "cl_i_mM", "tic_i_mM", "hco3_i_mM",
        "na_l_mM", "k_l_mM", "cl_l_mM", "tic_l_mM", "hco3_l_mM",
    )
    return {
        "optimizer_success": 1.0 if solution["optimizer_success"] else -1.0,
        "scaled_independent_rhs": 1.0 - float(solution["max_abs_scaled_independent_rhs"]) / ROOT_SCALED_TOLERANCE,
        "omitted_charge_rhs": 1.0 - float(solution["max_abs_omitted_rhs_fmol_s"]) / OMITTED_ROW_TOLERANCE,
        "regulatory_rhs": 1.0 - float(solution["max_abs_regulatory_rhs_s_inv"]) / REGULATORY_TOLERANCE_S_INV,
        "bulk_charge": 1.0 - float(solution["max_abs_bulk_charge_fmol"]) / CHARGE_TOLERANCE_FMOL,
        "current_closure": 1.0 - float(solution["max_abs_current_residual_A"]) / CURRENT_TOLERANCE_A,
        "conservation": 1.0 - float(solution["max_dimensionless_conservation_ratio"]),
        "jacobian_rank": 1.0 if solution["independent_jacobian_rank"] == 10 else -1.0,
        "no_coordinate_boundary": 1.0 if not solution["boundary_hits"] else -1.0,
        "positive_core": 1.0 if solution["positivity_pass"] else -1.0,
        "time_invariant_rest": 1.0 if solution["rest_protocol_time_invariant_rhs"] else -1.0,
        "finite_capacity": (1.0 + 1e-12) - capacity_max,
        "q_out_target": 1.0 - q_error / REST_Q_OUT_ABS_TOLERANCE_PL_S,
        "nkcc1_share": 1.0 - abs(share - 0.70) / 0.01,
        "cl_i_lower": (float(o["cl_i_mM"]) - 45.0) / 20.0,
        "cl_i_upper": (65.0 - float(o["cl_i_mM"])) / 20.0,
        "ph_i_lower": (float(o["ph_i"]) - 6.84) / 0.14,
        "ph_i_upper": (6.98 - float(o["ph_i"])) / 0.14,
        "ae4_positive": float(f["ae4_cl_cell_fmol_s"]) / max(abs(float(f["positive_chloride_loading_pool_fmol_s"])), 1e-15),
        "ae4_non_negligible": (ae4_share - 0.05) / 0.05,
        "na_i_upper": (30.0 - float(o["na_i_mM"])) / 30.0,
        "positive_concentrations": min(float(o[key]) for key in concentration_keys) / 100.0,
        "positive_volumes": min(float(o["volume_i_pL"]), float(o["volume_l_pL"])),
        "cell_volume_upper": (3.0 - float(o["volume_i_pL"])) / 3.0,
        "k_i_lower": (float(o["k_i_mM"]) - 60.0) / 150.0,
        "k_i_upper": (210.0 - float(o["k_i_mM"])) / 150.0,
        "tic_i_upper": (100.0 - float(o["tic_i_mM"])) / 100.0,
        "hco3_i_upper": (100.0 - float(o["hco3_i_mM"])) / 100.0,
    }


def _stationary_feasible(solution: Mapping[str, Any]) -> bool:
    margins = solution["constraint_margins"]
    o = solution["observables"]
    f = solution["fluxes"]
    return bool(
        all(math.isfinite(float(value)) and float(value) >= 0.0 for value in margins.values())
        and float(f["ae4_cl_cell_fmol_s"]) > 0.0
        and 0.0 < float(o["volume_i_pL"]) < 3.0
        and float(o["volume_l_pL"]) > 0.0
    )


class WTEvaluator:
    def __init__(self, budget: Budget, reference: Mapping[str, Any]):
        self.budget = budget
        self.reference = reference
        self.cache: dict[tuple[str, ...], dict[str, Any]] = {}
        self.records: list[dict[str, Any]] = []

    @staticmethod
    def key(z: ArrayLike) -> tuple[str, ...]:
        return tuple(float(value).hex() for value in np.asarray(z, dtype=float))

    def _warm_start(self, z: NDArray[np.float64]) -> NDArray[np.float64]:
        numerical = [record for record in self.records if record["solution"]["numerical_pass"]]
        if not numerical:
            return np.asarray(self.reference["coordinates"], dtype=float)
        selected = min(
            numerical,
            key=lambda record: (
                float(np.linalg.norm(z - np.asarray(record["log_multipliers"], dtype=float))),
                int(record["evaluation_index"]),
            ),
        )
        return np.asarray(selected["solution"]["coordinates"], dtype=float)

    def evaluate(self, z_values: ArrayLike) -> dict[str, Any]:
        z = np.asarray(z_values, dtype=float)
        key = self.key(z)
        if key in self.cache:
            return self.cache[key]
        if len(self.records) >= OPTIMISER_VECTOR_CEILING:
            raise Task33Stop("declared 19 vector optimiser ceiling reached")
        multipliers = multipliers_from_log(z)
        self.budget.start_wt_vector()
        model, _ = build_model(multipliers, stimulated=False)
        solution = solve_stationary(
            model,
            genotype=WT,
            start_coordinates=self._warm_start(z),
            multipliers=multipliers,
            label=f"WT_vector_{len(self.records) + 1}",
            budget=self.budget,
            prefreeze=True,
        )
        objective = _soft_objective(solution, multipliers, self.reference)
        margins = _constraint_margins(solution)
        solution["soft_objective"] = objective
        solution["constraint_margins"] = margins
        solution["wt_stationary_feasible"] = _stationary_feasible(solution)
        violation = float(sum(max(0.0, -float(value)) ** 2 for value in margins.values()))
        record = {
            "evaluation_index": len(self.records) + 1,
            "log_multipliers": z.tolist(),
            "multipliers": multipliers.tolist(),
            "soft_objective": objective,
            "constraint_violation": violation,
            "minimum_constraint_margin": min(float(value) for value in margins.values()),
            "solution": solution,
        }
        self.records.append(record)
        self.cache[key] = record
        self._checkpoint_trace()
        self.budget.save()
        print(json.dumps({
            "wt_vector": record["evaluation_index"],
            "multipliers": record["multipliers"],
            "objective": objective,
            "violation": violation,
            "stationary_feasible": solution["wt_stationary_feasible"],
            "q_out": solution["fluxes"]["water_lumen_outflow_pL_s"],
            "nkcc1_share": solution["fluxes"]["nkcc1_positive_pool_share"],
            "cl_i": solution["observables"]["cl_i_mM"],
            "ph_i": solution["observables"]["ph_i"],
        }), flush=True)
        return record

    def _checkpoint_trace(self) -> None:
        rows = []
        for record in self.records:
            solution = record["solution"]
            row: dict[str, Any] = {
                "evaluation_index": record["evaluation_index"],
                "soft_objective": record["soft_objective"],
                "constraint_violation": record["constraint_violation"],
                "minimum_constraint_margin": record["minimum_constraint_margin"],
                "stationary_feasible": solution["wt_stationary_feasible"],
                "numerical_pass": solution["numerical_pass"],
                "actual_residual_evaluations": solution["actual_residual_evaluations"],
                "q_out_pL_s": solution["fluxes"]["water_lumen_outflow_pL_s"],
                "q_out_relative_error": (
                    solution["fluxes"]["water_lumen_outflow_pL_s"] - REST_Q_OUT_TARGET_PL_S
                ) / REST_Q_OUT_TARGET_PL_S,
                "nkcc1_positive_pool_share": solution["fluxes"]["nkcc1_positive_pool_share"],
                "ae4_positive_pool_share": solution["fluxes"]["ae4_positive_pool_share"],
                "cl_i_mM": solution["observables"]["cl_i_mM"],
                "ph_i": solution["observables"]["ph_i"],
                "na_i_mM": solution["observables"]["na_i_mM"],
            }
            row.update({"multiplier_" + name: value for name, value in solution["parameter_multipliers"].items()})
            rows.append(row)
        _write_csv(RESULTS / "optimisation_trace.csv", rows)

    def objective(self, z: ArrayLike) -> float:
        return float(self.evaluate(z)["soft_objective"])

    def constraints(self, z: ArrayLike) -> NDArray[np.float64]:
        record = self.evaluate(z)
        return np.asarray(list(record["solution"]["constraint_margins"].values()), dtype=float)

    def least_infeasible(self) -> dict[str, Any]:
        return min(
            self.records,
            key=lambda record: (
                record["constraint_violation"],
                record["soft_objective"],
                record["evaluation_index"],
            ),
        )


def run_local_optimisation(evaluator: WTEvaluator, budget: Budget) -> list[dict[str, Any]]:
    bounds = Bounds(LOG_LOWER, LOG_UPPER)
    constraint = {"type": "ineq", "fun": evaluator.constraints}
    budget.optimiser_calls += 1
    first = minimize(
        evaluator.objective,
        INITIAL_LOG_MULTIPLIERS,
        method="COBYLA",
        bounds=bounds,
        constraints=constraint,
        options={"rhobeg": 0.35, "tol": 1e-4, "catol": 1e-8, "maxiter": OPTIMISER_VECTOR_CEILING, "disp": False},
    )
    optimiser_summaries = [{
        "call": 1,
        "success": bool(first.success),
        "message": str(first.message),
        "fun": float(first.fun),
        "nfev_reported": int(first.nfev),
        "distinct_vectors_after_call": len(evaluator.records),
    }]
    remaining = OPTIMISER_VECTOR_CEILING - len(evaluator.records)
    if remaining >= 7:
        budget.optimiser_calls += 1
        budget.deterministic_restarts += 1
        start = np.asarray(evaluator.least_infeasible()["log_multipliers"], dtype=float)
        second = minimize(
            evaluator.objective,
            start,
            method="COBYLA",
            bounds=bounds,
            constraints=constraint,
            options={"rhobeg": 0.05, "tol": 1e-4, "catol": 1e-8, "maxiter": remaining, "disp": False},
        )
        optimiser_summaries.append({
            "call": 2,
            "success": bool(second.success),
            "message": str(second.message),
            "fun": float(second.fun),
            "nfev_reported": int(second.nfev),
            "distinct_vectors_after_call": len(evaluator.records),
        })
    _write_json(RESULTS / "optimiser_summary.json", optimiser_summaries)
    budget.save()
    return optimiser_summaries


def run_dynamic(
    rest_solution: Mapping[str, Any],
    *,
    genotype: Genotype,
    budget: Budget,
    prefreeze: bool,
) -> dict[str, Any]:
    budget.start_dynamic(prefreeze=prefreeze)
    multipliers = np.asarray([rest_solution["parameter_multipliers"][name] for name in PARAMETER_NAMES])
    model, _ = build_model(multipliers, stimulated=True)
    rest_model, _ = build_model(multipliers, stimulated=False)
    y0 = np.r_[np.asarray(rest_solution["core_state"], dtype=float), model.initial_state()[12:]]
    onset = model.evaluate(0.0, y0, genotype=genotype)
    rest_onset = rest_model.evaluate(0.0, y0, genotype=genotype)
    onset_match = bool(np.array_equal(onset.rhs, rest_onset.rhs))
    grid = np.r_[1e-6, np.arange(1.0, 601.0)]

    def rhs(time_s: float, state: NDArray[np.float64]) -> NDArray[np.float64]:
        if prefreeze:
            budget.tick_wt()
        return model.rhs(time_s, state, genotype=genotype)

    result = solve_ivp(
        rhs,
        (float(grid[0]), 600.0),
        y0,
        method=PRODUCTION_RADAU.method,
        rtol=PRODUCTION_RADAU.rtol,
        atol=PRODUCTION_RADAU.atol_vector(model.state_names),
        t_eval=grid,
        max_step=PRODUCTION_RADAU.max_step_s,
    )
    complete = bool(result.y.shape[1] == grid.size and result.t.size and result.t[-1] == 600.0)
    times = np.r_[0.0, np.asarray(result.t, dtype=float)]
    states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    flux_rows: list[dict[str, Any]] = []
    observables_rows: list[dict[str, float]] = []
    conservation_ratios: list[float] = []
    capacity_passes: list[bool] = []
    evaluations_valid = True
    for time_s, state in zip(times, states.T):
        if prefreeze:
            budget.tick_wt()
        try:
            evaluation = model.evaluate(float(time_s), state, genotype=genotype)
            flux_rows.append(_extended_flux_row(evaluation, model))
            observables_rows.append(_rest_observables(evaluation, state))
            conservation_ratios.append(_conservation_audit(evaluation)[1])
            capacity_passes.append(_finite_capacity_audit(evaluation, model)["pass"])
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            evaluations_valid = False
            break
    if not evaluations_valid or not flux_rows:
        summary = {
            "genotype": genotype.name,
            "ae4_expression": genotype.ae4_expression,
            "solver_success": bool(result.success),
            "solver_message": str(result.message),
            "complete": complete,
            "onset_matches_rest_rhs": onset_match,
            "valid": False,
            "activation_pass": False if prefreeze else None,
            "failure": "trajectory diagnostics were not physically evaluable",
        }
        budget.save()
        return summary
    flow = np.asarray([float(row["water_lumen_outflow_pL_s"]) for row in flux_rows])
    cumulative = float(np.trapezoid(flow, times))
    mean_mask = times >= 60.0
    mean_60_600 = float(np.trapezoid(flow[mean_mask], times[mean_mask]) / 540.0)
    concentration_keys = (
        "na_i_mM", "k_i_mM", "cl_i_mM", "tic_i_mM", "hco3_i_mM",
        "na_l_mM", "k_l_mM", "cl_l_mM", "tic_l_mM", "hco3_l_mM",
    )
    positive_concentrations = all(
        math.isfinite(float(row[key])) and float(row[key]) > 0.0
        for row in observables_rows for key in concentration_keys
    )
    positive_volumes = all(
        math.isfinite(float(row[key])) and float(row[key]) > 0.0
        for row in observables_rows for key in ("volume_i_pL", "volume_l_pL")
    )
    finite_states = bool(np.all(np.isfinite(states)))
    positive_core = bool(np.all(states[:12] > 0.0))
    nonnegative_flow = bool(np.all(np.isfinite(flow)) and np.all(flow >= 0.0))
    valid = bool(
        result.success
        and complete
        and onset_match
        and evaluations_valid
        and finite_states
        and positive_core
        and positive_concentrations
        and positive_volumes
        and nonnegative_flow
        and max(conservation_ratios) <= 1.0
        and all(capacity_passes)
    )
    activation_pass = bool(
        valid
        and mean_60_600 >= ACTIVATION_FACTOR * REST_Q_OUT_TARGET_PL_S
        and cumulative >= ACTIVATION_FACTOR * 600.0 * REST_Q_OUT_TARGET_PL_S
        and flow[-1] >= REST_Q_OUT_TARGET_PL_S
    ) if prefreeze else None
    integral_keys = (
        "nkcc1_cl_in_fmol_s",
        "ae4_cl_cell_fmol_s",
        "ae2_cl_in_hco3_out_fmol_s",
        "cacc_cl_cell_to_lumen_fmol_s",
        "paracellular_cl_lumen_to_bath_fmol_s",
        "chloride_outflow_fmol_s",
        "water_bath_to_cell_pL_s",
        "water_cell_to_lumen_pL_s",
        "water_bath_to_lumen_pL_s",
        "water_lumen_outflow_pL_s",
    )
    integrated = {
        key.replace("_fmol_s", "_fmol").replace("_pL_s", "_pL"): float(
            np.trapezoid(np.asarray([float(row[key]) for row in flux_rows]), times)
        )
        for key in integral_keys
    }
    storage = {
        name.replace("_fmol", "_storage_fmol").replace("_pL", "_storage_pL"): float(states[index, -1] - states[index, 0])
        for index, name in enumerate(CORE_STATE_NAMES)
    }
    summary = {
        "genotype": genotype.name,
        "ae4_expression": float(genotype.ae4_expression),
        "calcium_uM": STIMULATED_CALCIUM_UM,
        "interval_s": [0.0, 600.0],
        "solver": PRODUCTION_RADAU.label,
        "solver_success": bool(result.success),
        "solver_message": str(result.message),
        "complete": complete,
        "sample_count": int(times.size),
        "onset_matches_rest_rhs": onset_match,
        "finite_states": finite_states,
        "positive_core": positive_core,
        "positive_concentrations": positive_concentrations,
        "positive_volumes": positive_volumes,
        "nonnegative_flow": nonnegative_flow,
        "max_dimensionless_conservation_ratio": float(max(conservation_ratios)),
        "finite_capacity_pass": bool(all(capacity_passes)),
        "rest_q_out_reference_pL_s": REST_Q_OUT_TARGET_PL_S,
        "mean_q_out_60_600_pL_s": mean_60_600,
        "mean_q_out_activation_ratio": mean_60_600 / REST_Q_OUT_TARGET_PL_S,
        "cumulative_secretion_600_pL": cumulative,
        "cumulative_activation_ratio": cumulative / (600.0 * REST_Q_OUT_TARGET_PL_S),
        "q_out_600_pL_s": float(flow[-1]),
        "endpoint_to_rest_q_out_ratio": float(flow[-1] / REST_Q_OUT_TARGET_PL_S),
        "integrated_fluxes": integrated,
        "ionic_and_volume_storage": storage,
        "endpoint_observables": observables_rows[-1],
        "endpoint_fluxes": flux_rows[-1],
        "valid": valid,
        "activation_pass": activation_pass,
    }
    budget.save()
    return summary


def _flatten_rest(solution: Mapping[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "background": solution["background"],
        "genotype": solution["genotype"],
        "ae4_expression": solution["ae4_expression"],
        "admissible": solution["admissible"],
        "numerical_pass": solution["numerical_pass"],
        "general_physical_pass": solution["general_physical_pass"],
        "exact_ae4_zero": solution["exact_ae4_zero"],
        "actual_residual_evaluations": solution["actual_residual_evaluations"],
        "max_abs_scaled_independent_rhs": solution["max_abs_scaled_independent_rhs"],
        "max_abs_omitted_rhs_fmol_s": solution["max_abs_omitted_rhs_fmol_s"],
        "max_abs_bulk_charge_fmol": solution["max_abs_bulk_charge_fmol"],
        "max_abs_current_residual_A": solution["max_abs_current_residual_A"],
        "max_dimensionless_conservation_ratio": solution["max_dimensionless_conservation_ratio"],
    }
    row.update(solution["observables"])
    return row


def _parameter_payload(multipliers: NDArray[np.float64]) -> dict[str, Any]:
    rest_model, _ = build_model(multipliers, stimulated=False)
    stimulated_model, _ = build_model(multipliers, stimulated=True)
    core = {
        "whole_cell_parameters": _jsonable(rest_model.parameters),
        "ae4_parameters": _jsonable(rest_model.ae4_parameters),
        "routing_model": AE4_ROUTING_MODEL_ID,
        "ae4_regulatory_model": {
            "family": rest_model.regulatory_model.ae4_regulatory_model.family,
            "tau_activation_s": 30.0,
            "basal_capacity_multiplier": 1.0,
            "fully_activated_increment": 0.25,
            "coupling_scale": 1.0,
            "construct": "WT",
        },
        "nkcc1_regulatory_model": {
            "family": rest_model.regulatory_model.nkcc1_regulatory_model.family,
            "fully_activated_multiplier": 1.75,
            "resting_calcium_uM": 0.058,
            "stimulated_calcium_uM": 0.10,
        },
        "rest_protocol": _jsonable(rest_model.stimulus),
        "stimulated_protocol": _jsonable(stimulated_model.stimulus),
        "wt_genotype": _jsonable(WT),
        "free_parameter_names": PARAMETER_NAMES,
        "multipliers": {name: float(value) for name, value in zip(PARAMETER_NAMES, multipliers)},
        "absolute_free_parameters": absolute_free_parameters(multipliers),
        "nhe1_kinetics": _jsonable(published_cha_kinetics()),
    }
    return {**core, "payload_sha256": sha256_object(core)}


def _write_prefreeze_failure(
    *, reason: str, evaluator: WTEvaluator | None, budget: Budget, optimiser: Any = None
) -> dict[str, Any]:
    if evaluator is not None and evaluator.records:
        # Make the compact trace agree with the already consumed counters even
        # if the optimiser exits immediately after its final allowed callback.
        evaluator._checkpoint_trace()
    best = evaluator.least_infeasible() if evaluator is not None and evaluator.records else None
    if best is not None:
        _write_json(RESULTS / "wt_state.json", best["solution"])
        _write_csv(RESULTS / "wt_flux_ledger.csv", [best["solution"]["fluxes"]])
    else:
        _write_json(RESULTS / "wt_state.json", {"status": "NO_WT_VECTOR_COMPLETED"})
        _write_csv(RESULTS / "wt_flux_ledger.csv", [])
    _write_json(RESULTS / "wt_dynamic_activation.json", {
        "status": "NO_ACCEPTED_WT_DYNAMIC_ACTIVATION",
        "evaluations": [],
    })
    budget._timing_wt = False
    budget.save()
    verification = {
        "status": "WT_OPTIMISATION_FAILED_WITHIN_PREDECLARED_BOUNDS",
        "reason": reason,
        "contract_sha256": sha256_file(CONTRACT),
        "contract_commit": CONTRACT_COMMIT,
        "prepared_head": PREPARED_HEAD,
        "branch": BRANCH,
        "budget": budget.payload(),
        "genotypes_evaluated": [],
        "bounds_changed": False,
        "targets_changed": False,
        "alternative_mechanism_diagnosed": False,
    }
    _write_json(RESULTS / "verification.json", verification)
    text = (
        "# Task 33 final answer\n\n"
        "`WT_OPTIMISATION_FAILED_WITHIN_PREDECLARED_BOUNDS`\n\n"
        f"Reason: {reason}. No genotype was evaluated, no bound or target was changed, and no further mechanism was diagnosed.\n"
    )
    (ANALYSIS / "final_answer.md").write_text(text, encoding="utf-8")
    _write_json(RESULTS / "stage.json", {"status": "WT_OPTIMISATION_FAILED_WITHIN_PREDECLARED_BOUNDS"})
    return verification


def run_wt_stage() -> dict[str, Any]:
    _check_environment()
    if RESULTS.exists() and any(RESULTS.iterdir()):
        raise RuntimeError("Task 33 results directory is not empty; numerical rerun is forbidden")
    budget = Budget()
    reference = task31_reference_wt()
    evaluator = WTEvaluator(budget, reference)
    optimiser_summary: Any = None
    try:
        optimiser_summary = run_local_optimisation(evaluator, budget)
        feasible = sorted(
            (record for record in evaluator.records if record["solution"]["wt_stationary_feasible"]),
            key=lambda record: (record["soft_objective"], record["evaluation_index"]),
        )
        dynamic_attempts: list[dict[str, Any]] = []
        accepted: dict[str, Any] | None = None
        for record in feasible[:PREFREEZE_DYNAMIC_CEILING]:
            dynamic = run_dynamic(record["solution"], genotype=WT, budget=budget, prefreeze=True)
            dynamic["wt_vector_evaluation_index"] = record["evaluation_index"]
            dynamic["soft_objective"] = record["soft_objective"]
            dynamic_attempts.append(dynamic)
            _write_json(RESULTS / "wt_dynamic_activation.json", {
                "status": "IN_PROGRESS",
                "evaluations": dynamic_attempts,
            })
            if dynamic["activation_pass"]:
                accepted = {"record": record, "dynamic": dynamic}
                break
        if accepted is None:
            if budget.prefreeze_wt_dynamic_integrations >= PREFREEZE_DYNAMIC_CEILING:
                reason = "three pre-freeze WT dynamic integrations completed without an activation pass"
            elif not feasible:
                reason = "the single local optimisation produced no stationary feasible WT candidate"
            else:
                reason = "all stationary feasible WT candidates failed the fixed activation gate"
            return _write_prefreeze_failure(
                reason=reason, evaluator=evaluator, budget=budget, optimiser=optimiser_summary
            )
        budget.accept_wt()
        solution = accepted["record"]["solution"]
        dynamic = accepted["dynamic"]
        multipliers = np.asarray([solution["parameter_multipliers"][name] for name in PARAMETER_NAMES])
        frozen_parameters = _parameter_payload(multipliers)
        boundary_parameters = [
            name for name, value, lower, upper in zip(
                PARAMETER_NAMES, multipliers, MULTIPLIER_LOWER, MULTIPLIER_UPPER
            ) if math.isclose(float(value), float(lower), rel_tol=0.0, abs_tol=1e-8)
            or math.isclose(float(value), float(upper), rel_tol=0.0, abs_tol=1e-8)
        ]
        _write_json(RESULTS / "frozen_wt_parameters.json", frozen_parameters)
        _write_json(RESULTS / "wt_state.json", solution)
        _write_csv(RESULTS / "wt_flux_ledger.csv", [solution["fluxes"]])
        _write_json(RESULTS / "wt_dynamic_activation.json", {
            "status": "WT_ACTIVATION_ACCEPTED_AND_FROZEN",
            "accepted": dynamic,
            "evaluations": dynamic_attempts,
        })
        budget.save()
        verification = {
            "status": "WT_FROZEN_AWAITING_REMOTE_CHECKPOINT",
            "contract_sha256": sha256_file(CONTRACT),
            "contract_commit": CONTRACT_COMMIT,
            "prepared_head": PREPARED_HEAD,
            "branch": BRANCH,
            "frozen_wt_payload_sha256": frozen_parameters["payload_sha256"],
            "selected_wt_vector": accepted["record"]["evaluation_index"],
            "wt_stationary_feasible": solution["wt_stationary_feasible"],
            "wt_activation_pass": dynamic["activation_pass"],
            "boundary_dependent_parameters": boundary_parameters,
            "nhe1_model": solution["nhe1_model"],
            "nhe1_carrier_amount_fmol": solution["nhe1_carrier_amount_fmol"],
            "genotypes_evaluated": [],
            "budget": budget.payload(),
        }
        _write_json(RESULTS / "verification.json", verification)
        _write_json(RESULTS / "stage.json", {
            "status": "WT_FROZEN_AWAITING_REMOTE_CHECKPOINT",
            "frozen_wt_payload_sha256": frozen_parameters["payload_sha256"],
        })
        return verification
    except Task33Stop as exc:
        return _write_prefreeze_failure(
            reason=str(exc), evaluator=evaluator, budget=budget, optimiser=optimiser_summary
        )


def _load_frozen() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    stage = json.loads((RESULTS / "stage.json").read_text(encoding="utf-8"))
    if stage["status"] != "WT_FROZEN_AWAITING_REMOTE_CHECKPOINT":
        raise RuntimeError("WT has not reached the pre-genotype frozen stage")
    parameters = json.loads((RESULTS / "frozen_wt_parameters.json").read_text(encoding="utf-8"))
    wt_state = json.loads((RESULTS / "wt_state.json").read_text(encoding="utf-8"))
    wt_dynamic = json.loads((RESULTS / "wt_dynamic_activation.json").read_text(encoding="utf-8"))["accepted"]
    multipliers = np.asarray([parameters["multipliers"][name] for name in PARAMETER_NAMES])
    rebuilt = _parameter_payload(multipliers)
    if rebuilt != parameters or rebuilt["payload_sha256"] != stage["frozen_wt_payload_sha256"]:
        raise RuntimeError("frozen WT parameter payload does not reconstruct exactly")
    return parameters, wt_state, wt_dynamic


def _decomposition_rows(
    wt_state: Mapping[str, Any],
    genotype_states: Sequence[Mapping[str, Any]],
    wt_dynamic: Mapping[str, Any],
    genotype_dynamics: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rest_keys = (
        "nkcc1_cl_in_fmol_s", "ae4_cl_cell_fmol_s", "ae2_cl_in_hco3_out_fmol_s",
        "cacc_cl_cell_to_lumen_fmol_s", "paracellular_cl_lumen_to_bath_fmol_s",
        "chloride_outflow_fmol_s", "water_bath_to_cell_pL_s", "water_cell_to_lumen_pL_s",
        "water_bath_to_lumen_pL_s", "water_lumen_outflow_pL_s",
    )
    for state in genotype_states:
        for key in rest_keys:
            value = float(state["fluxes"][key])
            wt_value = float(wt_state["fluxes"][key])
            rows.append({
                "genotype": state["genotype"],
                "ae4_expression": state["ae4_expression"],
                "phase": "REST",
                "component": key,
                "value": value,
                "frozen_wt_value": wt_value,
                "change_from_frozen_wt": value - wt_value,
                "units": "pL/s" if "water_" in key else "fmol/s",
            })
        dynamic = genotype_dynamics.get(state["genotype"])
        if dynamic and dynamic.get("valid"):
            for key, value in dynamic["integrated_fluxes"].items():
                wt_value = float(wt_dynamic["integrated_fluxes"][key])
                rows.append({
                    "genotype": state["genotype"],
                    "ae4_expression": state["ae4_expression"],
                    "phase": "STIMULATED_0_600",
                    "component": key,
                    "value": value,
                    "frozen_wt_value": wt_value,
                    "change_from_frozen_wt": float(value) - wt_value,
                    "units": "pL" if key.endswith("_pL") else "fmol",
                })
            for key, value in dynamic["ionic_and_volume_storage"].items():
                wt_value = float(wt_dynamic["ionic_and_volume_storage"][key])
                rows.append({
                    "genotype": state["genotype"],
                    "ae4_expression": state["ae4_expression"],
                    "phase": "STIMULATED_STORAGE_0_600",
                    "component": key,
                    "value": value,
                    "frozen_wt_value": wt_value,
                    "change_from_frozen_wt": float(value) - wt_value,
                    "units": "pL" if key.endswith("_pL") else "fmol",
                })
    return rows


def run_genotype_stage(*, freeze_commit: str) -> dict[str, Any]:
    _check_environment()
    if len(freeze_commit) != 40 or any(character not in "0123456789abcdef" for character in freeze_commit):
        raise ValueError("freeze commit must be a forty character lowercase Git SHA")
    frozen_parameters, wt_state, wt_dynamic = _load_frozen()
    budget = Budget.load()
    if not budget.wt_accepted:
        raise RuntimeError("budget does not record an accepted frozen WT")
    multipliers = np.asarray([frozen_parameters["multipliers"][name] for name in PARAMETER_NAMES])
    rest_model, _ = build_model(multipliers, stimulated=False)
    genotype_states: list[dict[str, Any]] = []
    genotype_dynamics: dict[str, dict[str, Any]] = {}
    for genotype in HELD_OUT_GENOTYPES:
        state = solve_stationary(
            rest_model,
            genotype=genotype,
            start_coordinates=wt_state["coordinates"],
            multipliers=multipliers,
            label=genotype.name + "_from_frozen_WT",
            budget=budget,
            prefreeze=False,
        )
        genotype_states.append(state)
        if state["admissible"]:
            genotype_dynamics[genotype.name] = run_dynamic(
                state, genotype=genotype, budget=budget, prefreeze=False
            )
        budget.save()
    _write_csv(RESULTS / "genotype_rest_states.csv", [_flatten_rest(state) for state in genotype_states])
    _write_csv(RESULTS / "genotype_flux_ledger.csv", [
        {"genotype": state["genotype"], "ae4_expression": state["ae4_expression"], **state["fluxes"]}
        for state in genotype_states
    ])
    comparison: list[dict[str, Any]] = [{
        "genotype": "WT",
        "ae4_expression": 1.0,
        "rest_admissible": True,
        "dynamic_run": True,
        "dynamic_valid": wt_dynamic["valid"],
        "cumulative_secretion_600_pL": wt_dynamic["cumulative_secretion_600_pL"],
        "secretion_relative_to_frozen_wt": 1.0,
        "mean_q_out_60_600_pL_s": wt_dynamic["mean_q_out_60_600_pL_s"],
        "q_out_600_pL_s": wt_dynamic["q_out_600_pL_s"],
    }]
    for state in genotype_states:
        dynamic = genotype_dynamics.get(state["genotype"])
        valid = bool(dynamic and dynamic.get("valid"))
        comparison.append({
            "genotype": state["genotype"],
            "ae4_expression": state["ae4_expression"],
            "rest_admissible": state["admissible"],
            "dynamic_run": dynamic is not None,
            "dynamic_valid": valid,
            "cumulative_secretion_600_pL": dynamic["cumulative_secretion_600_pL"] if valid else "",
            "secretion_relative_to_frozen_wt": (
                dynamic["cumulative_secretion_600_pL"] / wt_dynamic["cumulative_secretion_600_pL"]
                if valid else ""
            ),
            "mean_q_out_60_600_pL_s": dynamic["mean_q_out_60_600_pL_s"] if valid else "",
            "q_out_600_pL_s": dynamic["q_out_600_pL_s"] if valid else "",
        })
    _write_csv(RESULTS / "secretion_comparison.csv", comparison)
    _write_csv(
        RESULTS / "chloride_flux_decomposition.csv",
        _decomposition_rows(wt_state, genotype_states, wt_dynamic, genotype_dynamics),
    )
    _write_json(RESULTS / "genotype_dynamics.json", genotype_dynamics)
    budget.save()
    verification = {
        "status": "TASK_33_COMPLETED",
        "contract_sha256": sha256_file(CONTRACT),
        "contract_commit": CONTRACT_COMMIT,
        "freeze_publication_commit": freeze_commit,
        "frozen_wt_payload_sha256": frozen_parameters["payload_sha256"],
        "frozen_payload_reconstructed_exactly": True,
        "nhe1_model": rest_model.parameters.homeostasis.nhe1_model,
        "nhe1_carrier_amount_fmol": rest_model.parameters.homeostasis.nhe1_cha_carrier_amount_fmol,
        "held_out_ae4_expression_levels": [state["ae4_expression"] for state in genotype_states],
        "non_ae4_refit_after_freeze": False,
        "genotype_rest_admissibility": {state["genotype"]: state["admissible"] for state in genotype_states},
        "genotype_dynamic_runs": list(genotype_dynamics),
        "budget": budget.payload(),
    }
    _write_json(RESULTS / "verification.json", verification)
    _write_json(RESULTS / "stage.json", {
        "status": "TASK_33_COMPLETED",
        "freeze_publication_commit": freeze_commit,
        "frozen_wt_payload_sha256": frozen_parameters["payload_sha256"],
    })
    lines = [
        "# Task 33 final answer",
        "",
        "The bounded WT optimisation passed the stationary and fixed calcium activation gates. The complete WT payload was frozen before the two held out AE4 expression tests.",
        "",
        "## Frozen WT",
        "",
        f"* Parameter payload SHA256: `{frozen_parameters['payload_sha256']}`.",
        f"* Resting outflow: `{wt_state['fluxes']['water_lumen_outflow_pL_s']:.12g} pL/s`.",
        f"* Intracellular chloride: `{wt_state['observables']['cl_i_mM']:.9g} mM`; pH: `{wt_state['observables']['ph_i']:.9g}`; Na: `{wt_state['observables']['na_i_mM']:.9g} mM`.",
        f"* NKCC1 positive chloride loading share: `{wt_state['fluxes']['nkcc1_positive_pool_share']:.9g}`; AE4 share: `{wt_state['fluxes']['ae4_positive_pool_share']:.9g}`.",
        f"* Stimulated mean outflow ratio: `{wt_dynamic['mean_q_out_activation_ratio']:.9g}`; cumulative ratio to the resting flow integral: `{wt_dynamic['cumulative_activation_ratio']:.9g}`; endpoint ratio: `{wt_dynamic['endpoint_to_rest_q_out_ratio']:.9g}`.",
        "",
        "## Held out AE4 tests",
        "",
        "| AE4 expression | REST admissible | Cumulative secretion (pL) | Relative to frozen WT |",
        "|---:|:---:|---:|---:|",
    ]
    for row in comparison[1:]:
        cumulative = f"{row['cumulative_secretion_600_pL']:.9g}" if row["dynamic_valid"] else "not run or invalid"
        ratio = f"{row['secretion_relative_to_frozen_wt']:.9g}" if row["dynamic_valid"] else "not available"
        lines.append(f"| {row['ae4_expression']:.2f} | {row['rest_admissible']} | {cumulative} | {ratio} |")
    lines.extend((
        "",
        "Signed resting and integrated chloride replacement terms, ionic storage, and water effects are recorded in `chloride_flux_decomposition.csv`. No intermediate AE4 level, alternative calcium input, R10 background, AE2 loss, or post-freeze refit was evaluated.",
        "",
    ))
    (ANALYSIS / "final_answer.md").write_text("\n".join(lines), encoding="utf-8")
    return verification


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="stage", required=True)
    subparsers.add_parser("wt")
    genotype_parser = subparsers.add_parser("genotypes")
    genotype_parser.add_argument("--freeze-commit", required=True)
    arguments = parser.parse_args(argv)
    result = run_wt_stage() if arguments.stage == "wt" else run_genotype_stage(
        freeze_commit=arguments.freeze_commit
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
