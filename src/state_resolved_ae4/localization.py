"""Stage-B residual localization for the state-resolved AE4 investigation.

This module deliberately separates three questions which cannot be collapsed:

1. what the two released AE4-null observations constrain about the state;
2. what source would close the *fixed* reduced chassis at a chosen member of
   that state family; and
3. which independently supported module signatures or local equilibrium
   responses can supply that source.

The 2015 null secretion magnitude and trajectory never enter this module.
Stage-B entry is explicit: callers must pass ``stage_a_failed=True`` before a
released target state or correction can be constructed.  The historical
seven-state chassis uses uncertified flux, water, and time scales, so chemical,
current, and water blocks are never combined into a dimensioned Euclidean norm
without an explicit caller-supplied scaling.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import combinations
import math
from typing import Callable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ae4_mechanism_reconstruction.chassis import (
    AE4_KNOCKOUT,
    ChassisParameters,
    FixedChassis,
    ZeroAE4,
)


# Expanded amount/current topology used by the independently audited chassis
# signatures.  ``q_*`` are positive-charge equivalents in the stated membrane
# direction, not water fluxes.  ``v_i`` is a separate cell-water/volume row.
EXPANDED_COORDINATES = (
    "na_i",
    "k_i",
    "cl_i",
    "hco3_i",
    "h_i",
    "co2_i",
    "na_l",
    "k_l",
    "cl_l",
    "hco3_l",
    "q_apical_cell_to_lumen",
    "q_basolateral_cell_to_bath",
    "q_tight_lumen_to_bath",
    "v_i",
)

CELL_CHEMICAL_ROWS = tuple(range(6))
LUMEN_CHEMICAL_ROWS = tuple(range(6, 10))
CURRENT_ROWS = tuple(range(10, 13))
WATER_ROWS = (13,)
REDUCED_MEASURED_ROWS = (0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13)


def _vector(**entries: float) -> tuple[int | float, ...]:
    values = [0] * len(EXPANDED_COORDINATES)
    lookup = {name: index for index, name in enumerate(EXPANDED_COORDINATES)}
    for name, value in entries.items():
        values[lookup[name]] = value
    return tuple(values)


# Positive directions and signs exactly match the Agent-D source audit.  These
# are per-event topology signatures; they do not assert a measured capacity.
MODULE_SIGNATURES: Mapping[str, tuple[int | float, ...]] = {
    "basolateral_nkcc1": _vector(na_i=1, k_i=1, cl_i=2),
    "basolateral_ae2": _vector(cl_i=1, hco3_i=-1),
    "basolateral_nhe1": _vector(na_i=1, h_i=-1),
    "basolateral_nak_pump": _vector(
        na_i=-3, k_i=2, q_basolateral_cell_to_bath=1
    ),
    "apical_nak_pump": _vector(
        na_i=-3,
        k_i=2,
        na_l=3,
        k_l=-2,
        q_apical_cell_to_lumen=1,
    ),
    "basolateral_k_channel": _vector(
        k_i=-1, q_basolateral_cell_to_bath=1
    ),
    "apical_k_channel": _vector(
        k_i=-1, k_l=1, q_apical_cell_to_lumen=1
    ),
    "apical_cl_channel": _vector(
        cl_i=-1, cl_l=1, q_apical_cell_to_lumen=-1
    ),
    "effective_apical_hco3_loss": _vector(
        hco3_i=-1, hco3_l=1, q_apical_cell_to_lumen=-1
    ),
    "paracellular_na": _vector(na_l=1, q_tight_lumen_to_bath=-1),
    "paracellular_k": _vector(k_l=1, q_tight_lumen_to_bath=-1),
    "carbon_hydration": _vector(hco3_i=1, h_i=1, co2_i=-1),
    "basolateral_water_entry": _vector(v_i=1),
    "apical_water_exit": _vector(v_i=-1),
}


def signature_matrix(
    names: Sequence[str] | None = None,
    *,
    rows: Sequence[int] | None = None,
) -> NDArray[np.float64]:
    """Return module signatures as columns in the declared order."""

    selected = tuple(MODULE_SIGNATURES) if names is None else tuple(names)
    row_index = tuple(range(len(EXPANDED_COORDINATES))) if rows is None else tuple(rows)
    return np.asarray(
        [[MODULE_SIGNATURES[name][row] for name in selected] for row in row_index],
        dtype=float,
    )


def exact_rank(values: Sequence[Sequence[int | float | Fraction]]) -> int:
    """Compute structural rank by exact rational Gaussian elimination.

    Integer topology matrices should not acquire a tolerance-dependent rank.
    Floats are converted through their decimal strings; callers should reserve
    this routine for exact signature matrices rather than numerical Jacobians.
    """

    matrix = [[Fraction(str(value)) for value in row] for row in values]
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("all matrix rows must have equal length")
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row_index, row in enumerate(matrix):
            if row_index == rank or not row[column]:
                continue
            factor = row[column]
            matrix[row_index] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(row, matrix[rank], strict=True)
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def exact_signature_rank(
    names: Sequence[str] | None = None,
    *,
    rows: Sequence[int] | None = None,
) -> int:
    matrix = signature_matrix(names, rows=rows)
    return exact_rank(matrix.tolist())


def topology_difference(name_apical: str, name_basolateral: str) -> NDArray[np.float64]:
    """Return the exact apical-minus-basolateral signature difference."""

    return np.asarray(MODULE_SIGNATURES[name_apical], dtype=float) - np.asarray(
        MODULE_SIGNATURES[name_basolateral], dtype=float
    )


@dataclass(frozen=True)
class StageBTargetManifold:
    """The complete fixed-chassis state family allowed by released Cl/pH.

    The free coordinates are ``(na_l, k_l, height, na_i, k_i)``.  The observed
    chloride fixes state coordinate 5 and pH fixes free proton concentration.
    Intracellular bicarbonate is then derived from exact cellular
    electroneutrality.  Thus two observations reduce seven states to a
    five-dimensional manifold; no null Na, K, or volume value is invented.

    In coordinates ``u=1/height`` the constraint is affine:

    ``-Na_i - K_i + X*u + HCO3_i = H_i - Cl_i``.
    """

    cl_i_mM: float
    pH_i: float
    impermeant_charge_amount: float

    @property
    def h_i_mM(self) -> float:
        return 1000.0 * 10.0 ** (-self.pH_i)

    @property
    def dimension(self) -> int:
        return 5

    @property
    def affine_coefficients(self) -> tuple[float, float, float, float]:
        # Order: Na_i, K_i, u=1/height, HCO3_i.
        return (-1.0, -1.0, self.impermeant_charge_amount, 1.0)

    @property
    def affine_rhs(self) -> float:
        return self.h_i_mM - self.cl_i_mM

    def hco3_i_mM(self, *, na_i: float, k_i: float, height: float) -> float:
        if not all(math.isfinite(value) for value in (na_i, k_i, height)):
            raise ValueError("free coordinates must be finite")
        if height <= 0.0:
            raise ValueError("height must be positive")
        return (
            self.h_i_mM
            - self.cl_i_mM
            - self.impermeant_charge_amount / height
            + na_i
            + k_i
        )

    def state(
        self,
        *,
        na_l: float,
        k_l: float,
        height: float,
        na_i: float,
        k_i: float,
    ) -> NDArray[np.float64]:
        hco3_i = self.hco3_i_mM(na_i=na_i, k_i=k_i, height=height)
        state = np.asarray(
            (na_l, k_l, height, na_i, k_i, self.cl_i_mM, hco3_i),
            dtype=float,
        )
        if np.any(~np.isfinite(state)) or np.any(state <= 0.0):
            raise ValueError("this free-coordinate selection is not a positive state")
        return state

    def residuals(self, state: ArrayLike) -> Mapping[str, float]:
        y = np.asarray(state, dtype=float)
        if y.shape != (7,):
            raise ValueError("expected a seven-state vector")
        h_i = y[6] + y[5] + self.impermeant_charge_amount / y[2] - y[3] - y[4]
        p_h = math.log10(1000.0 / h_i) if h_i > 0.0 else math.nan
        return {
            "cl_i_mM": float(y[5] - self.cl_i_mM),
            "pH_i": float(p_h - self.pH_i),
            "electroneutrality_h_mM": float(h_i - self.h_i_mM),
        }


def released_stage_b_manifold(
    chassis: FixedChassis,
    *,
    stage_a_failed: bool,
    cl_i_mM: float,
    pH_i: float,
) -> StageBTargetManifold:
    """Construct the Stage-B target only after an explicit Stage-A failure."""

    if not stage_a_failed:
        raise PermissionError("Stage-B ionic targets remain sealed until Stage A fails")
    if cl_i_mM <= 0.0 or not 0.0 < pH_i < 14.0:
        raise ValueError("released Cl and pH targets must be physical")
    return StageBTargetManifold(
        cl_i_mM=float(cl_i_mM),
        pH_i=float(pH_i),
        impermeant_charge_amount=chassis.parameters.impermeant_charge_amount,
    )


@dataclass(frozen=True)
class BalanceCorrection:
    """Expanded source required to make one chosen target state stationary."""

    state: tuple[float, ...]
    required: tuple[float, ...]
    determined: tuple[bool, ...]
    raw_rhs: tuple[float, ...]
    pH_i: float
    note: str

    def as_mapping(self) -> Mapping[str, float | None]:
        return {
            name: (float(value) if known else None)
            for name, value, known in zip(
                EXPANDED_COORDINATES, self.required, self.determined, strict=True
            )
        }


def required_balance_correction(
    chassis: FixedChassis,
    state: ArrayLike,
    *,
    stage_a_failed: bool,
    calcium: float = 0.05,
) -> BalanceCorrection:
    """Evaluate the AE4-null source demand at one target-manifold member.

    This does not estimate the free target coordinates.  It computes a
    conditional source vector for the caller's explicit choice.  CO2 and
    luminal HCO3 demands are undetermined because the fixed chassis eliminates
    or omits those states.  Current rows are zero source demands after the
    existing QSS voltage closure; an electrogenic added module therefore needs
    a compensating current pathway and a re-solved voltage closure.
    """

    if not stage_a_failed:
        raise PermissionError("Stage-B correction remains sealed until Stage A fails")
    y = np.asarray(state, dtype=float)
    zero = FixedChassis(
        ZeroAE4(),
        parameters=chassis.parameters,
        calcium_protocol=chassis.calcium_protocol,
        pka_protocol=chassis.pka_protocol,
        pka_discontinuities=chassis.pka_discontinuities,
    )
    raw, diagnostics = zero.evaluate(y, calcium=calcium, scenario=AE4_KNOCKOUT)
    amount_na = float(diagnostics["amount_rate_na"])
    amount_k = float(diagnostics["amount_rate_k"])
    amount_cl = float(diagnostics["amount_rate_cl"])
    amount_hco3 = float(diagnostics["amount_rate_hco3"])
    # This is exact in the reduced chassis by differentiated electroneutrality.
    amount_h = -amount_na - amount_k + amount_cl + amount_hco3
    values = np.zeros(len(EXPANDED_COORDINATES), dtype=float)
    known = np.ones(len(EXPANDED_COORDINATES), dtype=bool)
    values[0:5] = -np.asarray(
        (amount_na, amount_k, amount_cl, amount_hco3, amount_h), dtype=float
    )
    values[5] = math.nan
    known[5] = False
    values[6] = -float(raw[0])
    values[7] = -float(raw[1])
    # The reduced lumen chloride balance equals the sum of Na and K balances.
    values[8] = -float(raw[0] + raw[1])
    values[9] = math.nan
    known[9] = False
    values[10:13] = 0.0
    values[13] = -float(raw[2])
    return BalanceCorrection(
        state=tuple(map(float, y)),
        required=tuple(map(float, values)),
        determined=tuple(map(bool, known)),
        raw_rhs=tuple(map(float, raw)),
        pH_i=float(diagnostics["pH_i"]),
        note=(
            "conditional on unmeasured Na_i, K_i, height, Na_l, and K_l; "
            "chemical and water entries retain different uncertified units"
        ),
    )


def candidate_independent_ae4_null_rhs(
    chassis_builders: Sequence[Callable[[], FixedChassis]],
    state: ArrayLike,
    *,
    time: float = 0.0,
    atol: float = 1e-13,
) -> tuple[NDArray[np.float64], float]:
    """Verify that all AE4 laws vanish to the same null vector field.

    Returns the common raw RHS and maximum pairwise absolute discrepancy.
    This theorem is conditional on candidates entering only through the AE4
    source boundary and on a common non-AE4 parameter fingerprint.
    """

    if not chassis_builders:
        raise ValueError("at least one chassis builder is required")
    chassis = [builder() for builder in chassis_builders]
    fingerprints = {item.parameter_fingerprint for item in chassis}
    if len(fingerprints) != 1:
        raise ValueError("candidate-independence requires one non-AE4 fingerprint")
    fields = [item.raw_rhs(time, state, scenario=AE4_KNOCKOUT) for item in chassis]
    common = np.asarray(fields[0], dtype=float)
    maximum = max(float(np.max(np.abs(np.asarray(field) - common))) for field in fields)
    if maximum > atol:
        raise AssertionError(f"AE4-null vector fields differ by {maximum:.6g}")
    return common, maximum


def numerical_jacobian(
    function: Callable[[NDArray[np.float64]], ArrayLike],
    point: ArrayLike,
    *,
    relative_step: float = 1e-6,
) -> NDArray[np.float64]:
    """Central-difference Jacobian with coordinate-scaled steps."""

    x = np.asarray(point, dtype=float)
    base = np.asarray(function(x), dtype=float)
    jacobian = np.empty((base.size, x.size), dtype=float)
    for index in range(x.size):
        step = relative_step * max(1.0, abs(float(x[index])))
        plus = x.copy()
        minus = x.copy()
        plus[index] += step
        minus[index] -= step
        jacobian[:, index] = (
            np.asarray(function(plus), dtype=float)
            - np.asarray(function(minus), dtype=float)
        ) / (2.0 * step)
    return jacobian


def ionic_observables(chassis: FixedChassis, state: ArrayLike) -> NDArray[np.float64]:
    y = np.asarray(state, dtype=float)
    p = chassis.parameters
    h_i = y[6] + y[5] + p.impermeant_charge_amount / y[2] - y[3] - y[4]
    return np.asarray((y[5], math.log10(1000.0 / h_i)), dtype=float)


def _reduced_coordinates(chassis: FixedChassis, state: ArrayLike) -> NDArray[np.float64]:
    """Electroneutral coordinates with positive H represented by ``log(H)``."""

    y = np.asarray(state, dtype=float)
    if y.shape != (7,):
        raise ValueError("expected a seven-state vector")
    p = chassis.parameters
    h_i = y[6] + y[5] + p.impermeant_charge_amount / y[2] - y[3] - y[4]
    if h_i <= 0.0:
        raise ValueError("state has nonpositive proton concentration")
    return np.asarray((y[0], y[1], y[3], y[4], y[5], y[6], math.log(h_i)))


def _state_from_reduced(
    chassis: FixedChassis, reduced: ArrayLike
) -> NDArray[np.float64]:
    x = np.asarray(reduced, dtype=float)
    if x.shape != (7,):
        raise ValueError("expected seven reduced coordinates")
    na_l, k_l, na_i, k_i, cl_i, hco3_i, log_h_i = map(float, x)
    h_i = math.exp(log_h_i)
    denominator = na_i + k_i + h_i - cl_i - hco3_i
    if denominator <= 0.0:
        raise ValueError("reduced coordinates imply nonpositive cell height")
    height = chassis.parameters.impermeant_charge_amount / denominator
    return np.asarray((na_l, k_l, height, na_i, k_i, cl_i, hco3_i), dtype=float)


@dataclass(frozen=True)
class LocalResponse:
    parameter: str
    d_cl_d_log_parameter: float
    d_pH_d_log_parameter: float
    raw_state_response: tuple[float, ...]
    jacobian_condition: float

    @property
    def ionic_vector(self) -> NDArray[np.float64]:
        return np.asarray((self.d_cl_d_log_parameter, self.d_pH_d_log_parameter))


def local_parameter_response(
    chassis: FixedChassis,
    root: ArrayLike,
    parameter: str,
    *,
    log_step: float = 1e-5,
    state_step: float = 1e-5,
) -> LocalResponse:
    """Implicit-function response of the AE4-null root to one log parameter.

    Changing pump or K-channel parameters here re-solves the fixed chassis's
    voltage equations inside every RHS evaluation, so their current feedback
    is retained.  This is a local diagnostic, not a finite retuning and not a
    claim that a genotype changes that parameter.
    """

    if not hasattr(chassis.parameters, parameter):
        raise KeyError(parameter)
    baseline_value = float(getattr(chassis.parameters, parameter))
    if baseline_value <= 0.0:
        raise ValueError("log-parameter response requires a positive parameter")
    y = np.asarray(root, dtype=float)
    reduced_root = _reduced_coordinates(chassis, y)

    def field(state: NDArray[np.float64], parameters: ChassisParameters) -> NDArray[np.float64]:
        probe = FixedChassis(
            ZeroAE4(),
            parameters=parameters,
            calcium_protocol=chassis.calcium_protocol,
            pka_protocol=chassis.pka_protocol,
            pka_discontinuities=chassis.pka_discontinuities,
        )
        return probe.raw_rhs(0.0, state, scenario=AE4_KNOCKOUT)

    jacobian = numerical_jacobian(
        lambda reduced: field(
            _state_from_reduced(chassis, reduced), chassis.parameters
        ),
        reduced_root,
        relative_step=state_step,
    )
    plus_parameters = replace(
        chassis.parameters, **{parameter: baseline_value * math.exp(log_step)}
    )
    minus_parameters = replace(
        chassis.parameters, **{parameter: baseline_value * math.exp(-log_step)}
    )
    d_field = (
        field(y, plus_parameters) - field(y, minus_parameters)
    ) / (2.0 * log_step)
    reduced_response = -np.linalg.solve(jacobian, d_field)
    state_map_jacobian = numerical_jacobian(
        lambda reduced: _state_from_reduced(chassis, reduced),
        reduced_root,
        relative_step=state_step,
    )
    state_response = state_map_jacobian @ reduced_response
    ionic_response = np.asarray(
        (reduced_response[4], -reduced_response[6] / math.log(10.0)), dtype=float
    )
    return LocalResponse(
        parameter=parameter,
        d_cl_d_log_parameter=float(ionic_response[0]),
        d_pH_d_log_parameter=float(ionic_response[1]),
        raw_state_response=tuple(map(float, state_response)),
        jacobian_condition=float(np.linalg.cond(jacobian)),
    )


@dataclass(frozen=True)
class SparseProjection:
    names: tuple[str, ...]
    coefficients: tuple[float, ...]
    alignment: float
    residual_fraction: float
    exact_within_tolerance: bool
    condition_number: float
    positive_only: bool


def _projection_metrics(target: NDArray[np.float64], prediction: NDArray[np.float64]) -> tuple[float, float]:
    target_norm = float(np.linalg.norm(target))
    prediction_norm = float(np.linalg.norm(prediction))
    if target_norm == 0.0:
        raise ValueError("target direction must be nonzero")
    alignment = (
        float(np.dot(target, prediction) / (target_norm * prediction_norm))
        if prediction_norm > 0.0
        else 0.0
    )
    residual = float(np.linalg.norm(target - prediction) / target_norm)
    return alignment, residual


def sparse_response_projections(
    target: ArrayLike,
    responses: Mapping[str, ArrayLike],
    *,
    scales: ArrayLike,
    maximum_modules: int = 2,
    positive_only: bool = False,
    exact_tolerance: float = 1e-10,
) -> tuple[SparseProjection, ...]:
    """Enumerate sparse local fits with an explicit coordinate scaling.

    Requiring ``scales`` prevents a silent Euclidean norm across quantities
    such as millimolar chloride and dimensionless pH.  Unit scales, reported
    SEM scales, and other choices remain separate sensitivity analyses.
    """

    target_vector = np.asarray(target, dtype=float)
    if target_vector.ndim != 1:
        raise ValueError("target must be one-dimensional")
    scale_vector = np.asarray(scales, dtype=float)
    if scale_vector.shape != target_vector.shape or np.any(scale_vector <= 0.0):
        raise ValueError("scales must be positive and match target")
    weighted_target = target_vector / scale_vector
    names = tuple(responses)
    rows: list[SparseProjection] = []
    for count in range(1, maximum_modules + 1):
        for subset in combinations(names, count):
            matrix = np.column_stack(
                [np.asarray(responses[name], dtype=float) / scale_vector for name in subset]
            )
            if matrix.shape[0] != weighted_target.size:
                raise ValueError("every response must match the target dimension")
            if positive_only:
                # For one or two columns, exhaustive active-set enumeration is
                # exact and avoids adding an optimizer dependency here.
                candidates = [np.zeros(count)]
                for active_count in range(1, count + 1):
                    for active in combinations(range(count), active_count):
                        active_matrix = matrix[:, active]
                        coefficients, *_ = np.linalg.lstsq(
                            active_matrix, weighted_target, rcond=None
                        )
                        if np.all(coefficients >= 0.0):
                            candidate = np.zeros(count)
                            candidate[list(active)] = coefficients
                            candidates.append(candidate)
                coefficients = min(
                    candidates,
                    key=lambda item: float(np.linalg.norm(weighted_target - matrix @ item)),
                )
            else:
                coefficients, *_ = np.linalg.lstsq(matrix, weighted_target, rcond=None)
            prediction = matrix @ coefficients
            alignment, residual = _projection_metrics(weighted_target, prediction)
            rows.append(
                SparseProjection(
                    names=tuple(subset),
                    coefficients=tuple(map(float, coefficients)),
                    alignment=alignment,
                    residual_fraction=residual,
                    exact_within_tolerance=bool(residual <= exact_tolerance),
                    condition_number=float(np.linalg.cond(matrix)),
                    positive_only=positive_only,
                )
            )
    return tuple(sorted(rows, key=lambda row: (row.residual_fraction, len(row.names), row.names)))


def fixed_row_cation_repair(
    magnitude: float,
) -> Mapping[str, object]:
    """Exact Task-12 cation correction and its minimal positive repair.

    ``d=(-m,+m)`` is matched by ``m/2`` pump cycles plus ``m/2`` NHE cycles.
    Pump alone leaves the exact relative residual ``1/sqrt(26)``.
    """

    if magnitude < 0.0 or not math.isfinite(magnitude):
        raise ValueError("magnitude must be finite and nonnegative")
    target = np.asarray((-magnitude, magnitude), dtype=float)
    pump = np.asarray((-3.0, 2.0))
    nhe = np.asarray((1.0, 0.0))
    coefficient = float(np.dot(target, pump) / np.dot(pump, pump)) if magnitude else 0.0
    return {
        "target": tuple(map(float, target)),
        "pump_best_coefficient": coefficient,
        "pump_alignment_absolute": 5.0 / math.sqrt(26.0),
        "pump_residual_fraction": 1.0 / math.sqrt(26.0),
        "positive_exact_coefficients": {
            "basolateral_nak_pump": magnitude / 2.0,
            "basolateral_nhe1": magnitude / 2.0,
        },
        "exact_reconstruction": tuple(map(float, magnitude / 2.0 * (pump + nhe))),
    }


def evaluated_target_slices(
    chassis: FixedChassis,
    root: ArrayLike,
    manifold: StageBTargetManifold,
    *,
    stage_a_failed: bool,
) -> Mapping[str, BalanceCorrection]:
    """Two transparent conditional slices through the five-dimensional family.

    ``model_anchored`` holds every unmeasured coordinate at the fixed-chassis
    null root.  ``wt_reference_cations`` substitutes the historical WT Na/K
    and volume coordinates while retaining the root lumen coordinates.  These
    are sensitivity slices, not estimates of an unmeasured null state.
    """

    y = np.asarray(root, dtype=float)
    if y.shape != (7,):
        raise ValueError("expected a seven-state root")
    model_anchored = manifold.state(
        na_l=float(y[0]),
        k_l=float(y[1]),
        height=float(y[2]),
        na_i=float(y[3]),
        k_i=float(y[4]),
    )
    wt_reference = manifold.state(
        na_l=float(y[0]),
        k_l=float(y[1]),
        height=28.7,
        na_i=25.0,
        k_i=120.0,
    )
    return {
        "model_anchored": required_balance_correction(
            chassis, model_anchored, stage_a_failed=stage_a_failed
        ),
        "wt_reference_cations": required_balance_correction(
            chassis, wt_reference, stage_a_failed=stage_a_failed
        ),
    }


@dataclass(frozen=True)
class CorrectionFamilyJacobian:
    free_coordinates: tuple[str, ...]
    determined_coordinates: tuple[str, ...]
    jacobian: tuple[tuple[float, ...], ...]
    singular_values: tuple[float, ...]
    numerical_rank: int


def correction_family_jacobian(
    chassis: FixedChassis,
    manifold: StageBTargetManifold,
    free_point: ArrayLike,
    *,
    stage_a_failed: bool,
    relative_step: float = 1e-5,
) -> CorrectionFamilyJacobian:
    """Differentiate required source across the five free target coordinates.

    ``free_point`` is ordered ``(Na_l,K_l,height,Na_i,K_i)``.  The derivative
    is reported only for determined correction rows and is a local numerical
    diagnostic.  Its rank tests whether the source vector is constant across
    the observationally indistinguishable target family.
    """

    if not stage_a_failed:
        raise PermissionError("Stage-B correction remains sealed until Stage A fails")
    point = np.asarray(free_point, dtype=float)
    if point.shape != (5,):
        raise ValueError("free_point must have five coordinates")
    determined = tuple(
        index
        for index in range(len(EXPANDED_COORDINATES))
        if index not in (5, 9)
    )

    def correction(free: NDArray[np.float64]) -> NDArray[np.float64]:
        state = manifold.state(
            na_l=float(free[0]),
            k_l=float(free[1]),
            height=float(free[2]),
            na_i=float(free[3]),
            k_i=float(free[4]),
        )
        result = required_balance_correction(
            chassis, state, stage_a_failed=stage_a_failed
        )
        return np.asarray(result.required, dtype=float)[list(determined)]

    jacobian = numerical_jacobian(correction, point, relative_step=relative_step)
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    tolerance = max(jacobian.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    return CorrectionFamilyJacobian(
        free_coordinates=("na_l", "k_l", "height", "na_i", "k_i"),
        determined_coordinates=tuple(EXPANDED_COORDINATES[index] for index in determined),
        jacobian=tuple(tuple(map(float, row)) for row in jacobian),
        singular_values=tuple(map(float, singular_values)),
        numerical_rank=rank,
    )


__all__ = [
    "BalanceCorrection",
    "CELL_CHEMICAL_ROWS",
    "CURRENT_ROWS",
    "CorrectionFamilyJacobian",
    "EXPANDED_COORDINATES",
    "LUMEN_CHEMICAL_ROWS",
    "LocalResponse",
    "MODULE_SIGNATURES",
    "REDUCED_MEASURED_ROWS",
    "SparseProjection",
    "StageBTargetManifold",
    "WATER_ROWS",
    "candidate_independent_ae4_null_rhs",
    "correction_family_jacobian",
    "evaluated_target_slices",
    "exact_rank",
    "exact_signature_rank",
    "fixed_row_cation_repair",
    "ionic_observables",
    "local_parameter_response",
    "numerical_jacobian",
    "released_stage_b_manifold",
    "required_balance_correction",
    "signature_matrix",
    "sparse_response_projections",
    "topology_difference",
]
