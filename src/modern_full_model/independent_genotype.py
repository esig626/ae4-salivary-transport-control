"""Independent, target-free genotype continuation for native Task 13B roots.

This module is deliberately separate from :mod:`genotype_evaluation` and the
production model.  It transcribes the ten physical charge-manifold
coordinates, their broad search bounds, the ten independent steady-state
rows, and the local continuation rule a second time.  No phenotype target or
target reader is present.

All structural parameters and the fitted intracellular OTHER osmole amount
remain fixed.  Only the ten resting-state coordinates are solved as AE4 or
AE2 expression is continued from one to exact zero on a 21-point grid.
Branch loss is returned as data; it is never repaired by fitting a parameter.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .independent import CORE_SIZE, IndependentWholeCell
from .independent_native import (
    NATIVE_ROOT_ROWS,
    state_from_independent_rest_coordinates,
)


TRANSPORTERS = ("AE4", "AE2")
EXPRESSION_GRID = tuple(float(value) for value in np.linspace(1.0, 0.0, 21))

# Independent transcription of the broad WT coordinate domains.  They are
# search domains, not physiological acceptance intervals.
COORDINATE_NAMES = (
    "na_i_mM",
    "k_i_mM",
    "tic_i_mM",
    "ph_i",
    "volume_i_pL",
    "na_l_mM",
    "k_l_mM",
    "tic_l_mM",
    "ph_l",
    "volume_l_pL",
)
COORDINATE_LOWER = np.asarray(
    (4.0, 60.0, 1.0, 6.2, 0.5, 20.0, 0.05, 0.2, 5.0, 0.02),
    dtype=float,
)
COORDINATE_UPPER = np.asarray(
    (60.0, 210.0, 100.0, 8.0, 3.0, 280.0, 60.0, 100.0, 9.0, 0.80),
    dtype=float,
)
INDEPENDENT_ROW_SCALES = np.asarray(
    (0.05, 0.05, 0.05, 0.02, 1.0e-4, 0.02, 0.02, 0.02, 0.02, 1.0e-4),
    dtype=float,
)


@dataclass(frozen=True)
class IndependentGenotypeAttempt:
    """One bounded solve from one declared local start."""

    start_id: str
    optimizer_success: bool
    admissible_state: bool
    passes_numerical_gate: bool
    max_abs_scaled_residual: float
    nfev: int
    cost: float
    optimality: float
    message: str
    coordinates: tuple[float, ...]


@dataclass(frozen=True)
class IndependentGenotypeRoot:
    """One independently retained root at fixed transporter expression."""

    root_id: str
    transporter: str
    expression: float
    member_start_ids: tuple[str, ...]
    coordinates: tuple[float, ...]
    core_state: tuple[float, ...]
    max_abs_scaled_residual: float
    max_abs_omitted_charge_rhs_fmol_s: float
    max_abs_bulk_charge_residual_fmol: float
    max_abs_current_residual_A: float
    max_abs_regulatory_rhs_s_inv: float
    jacobian_singular_values: tuple[float, ...]
    jacobian_rank: int
    jacobian_nullity: int
    boundary_hits: tuple[str, ...]
    exact_deleted_transport_zero: bool
    passes_numerical_gate: bool
    gate_failures: tuple[str, ...]


@dataclass(frozen=True)
class IndependentGenotypeStep:
    """All independent attempts and retained branches at one expression."""

    expression: float
    attempts: tuple[IndependentGenotypeAttempt, ...]
    roots: tuple[IndependentGenotypeRoot, ...]
    selected_root_id: str | None
    selected_distance_from_previous_normalized: float | None
    connected_to_previous: bool
    status: str

    @property
    def selected_root(self) -> IndependentGenotypeRoot | None:
        if self.selected_root_id is None:
            return None
        return next(root for root in self.roots if root.root_id == self.selected_root_id)


@dataclass(frozen=True)
class IndependentGenotypeContinuation:
    """Complete or branch-lost WT-to-zero continuation report."""

    transporter: str
    expression_grid: tuple[float, ...]
    steps: tuple[IndependentGenotypeStep, ...]
    completed_to_exact_zero: bool
    branch_failure: str | None

    @property
    def connected_null_root(self) -> IndependentGenotypeRoot | None:
        if not self.completed_to_exact_zero or not self.steps:
            return None
        root = self.steps[-1].selected_root
        if root is None or root.expression != 0.0:
            return None
        return root


def _checked_transporter(transporter: str) -> str:
    value = str(transporter).upper()
    if value not in TRANSPORTERS:
        raise ValueError(f"transporter must be one of {TRANSPORTERS}")
    return value


def independent_coordinates_from_state(
    model: IndependentWholeCell, state: ArrayLike
) -> NDArray[np.float64]:
    """Decode the ten physiological coordinates without a production helper."""

    vector = np.asarray(state, dtype=float)
    if vector.ndim != 1 or vector.size < CORE_SIZE:
        raise ValueError("independent state must contain the twelve-state core")
    evaluation = model.evaluate(0.0, vector)
    return np.asarray(
        (
            evaluation.cell_concentrations_mM["na"],
            evaluation.cell_concentrations_mM["k"],
            evaluation.cell_concentrations_mM["tic"],
            evaluation.cell_speciation.ph,
            vector[5],
            evaluation.lumen_concentrations_mM["na"],
            evaluation.lumen_concentrations_mM["k"],
            evaluation.lumen_concentrations_mM["tic"],
            evaluation.lumen_speciation.ph,
            vector[11],
        ),
        dtype=float,
    )


def _basal_suffix(model: IndependentWholeCell) -> tuple[float, ...]:
    ae4 = () if model.regulation is None else tuple(model.regulation.basal_state)
    nkcc = () if model.nkcc_regulation is None else tuple(model.nkcc_regulation.basal_state)
    return tuple(float(value) for value in (*ae4, *nkcc))


def _state_from_coordinates(
    model: IndependentWholeCell, coordinates: ArrayLike
) -> NDArray[np.float64]:
    core = state_from_independent_rest_coordinates(model.parameters, coordinates)
    return np.r_[core, np.asarray(_basal_suffix(model), dtype=float)]


def _evaluate_expression(
    model: IndependentWholeCell,
    state: ArrayLike,
    transporter: str,
    expression: float,
) -> Any:
    if transporter == "AE4":
        return model.evaluate(0.0, state, ae4_expression=expression)
    return model.evaluate(0.0, state, ae2_expression=expression)


def _exact_deleted_transport_zero(
    evaluation: Any, transporter: str, expression: float
) -> bool:
    if expression != 0.0:
        return True
    if transporter == "AE4":
        return bool(
            float(evaluation.fluxes_fmol_s["ae4_cl"]) == 0.0
            and float(evaluation.conservation_residuals["ae4_charge_fmol_s"]) == 0.0
        )
    return bool(float(evaluation.fluxes_fmol_s["ae2"]) == 0.0)


def _one_attempt(
    model: IndependentWholeCell,
    *,
    transporter: str,
    expression: float,
    start_id: str,
    start: NDArray[np.float64],
    max_nfev: int,
) -> tuple[IndependentGenotypeAttempt, IndependentGenotypeRoot | None]:
    span = COORDINATE_UPPER - COORDINATE_LOWER

    def residual(coordinates: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            state = _state_from_coordinates(model, coordinates)
            evaluation = _evaluate_expression(
                model, state, transporter, expression
            )
            return (
                np.asarray(evaluation.rhs, dtype=float)[list(NATIVE_ROOT_ROWS)]
                / INDEPENDENT_ROW_SCALES
            )
        except (
            ValueError,
            FloatingPointError,
            OverflowError,
            ZeroDivisionError,
            np.linalg.LinAlgError,
        ):
            return np.full(10, 1.0e6)

    fit = least_squares(
        residual,
        np.asarray(start, dtype=float),
        bounds=(COORDINATE_LOWER, COORDINATE_UPPER),
        max_nfev=max_nfev,
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        x_scale="jac",
    )
    coordinates = np.asarray(fit.x, dtype=float)
    scaled_max = float(np.max(np.abs(residual(coordinates))))
    admissible = True
    root: IndependentGenotypeRoot | None = None
    try:
        state = _state_from_coordinates(model, coordinates)
        evaluation = _evaluate_expression(model, state, transporter, expression)
        raw = np.asarray(evaluation.rhs, dtype=float)
        omitted_max = float(np.max(np.abs(raw[[4, 10]])))
        fixed = model.parameters.geometry.fixed_cell_anion_equivalents_fmol
        charge_max = float(
            max(
                abs(state[0] + state[1] - state[2] - state[4] - fixed),
                abs(state[6] + state[7] - state[8] - state[10]),
            )
        )
        current_max = float(
            max(abs(float(value)) for value in evaluation.current_residuals_A.values())
        )
        regulatory_max = (
            float(np.max(np.abs(raw[CORE_SIZE:])))
            if raw.size > CORE_SIZE
            else 0.0
        )
        normalized_jacobian = np.asarray(fit.jac, dtype=float) * span[None, :]
        singular = np.linalg.svd(normalized_jacobian, compute_uv=False)
        threshold = 1.0e-8 * singular[0] if singular.size else math.inf
        rank = int(np.count_nonzero(singular > threshold))
        relative_boundary_distance = np.minimum(
            coordinates - COORDINATE_LOWER,
            COORDINATE_UPPER - coordinates,
        ) / span
        hits = tuple(
            name
            for name, distance in zip(COORDINATE_NAMES, relative_boundary_distance)
            if distance <= 1.0e-5
        )
        exact_zero = _exact_deleted_transport_zero(
            evaluation, transporter, expression
        )
        failures: list[str] = []
        if not fit.success:
            failures.append("optimizer did not converge")
        if scaled_max > 1.0e-7:
            failures.append("scaled independent residual exceeds 1e-7")
        if omitted_max > 1.0e-9:
            failures.append("omitted charge row exceeds 1e-9 fmol/s")
        if charge_max > 1.0e-9:
            failures.append("bulk charge residual exceeds 1e-9 fmol")
        if current_max > 1.0e-18:
            failures.append("membrane current residual exceeds 1e-18 A")
        if regulatory_max > 1.0e-10:
            failures.append("basal regulatory suffix is nonstationary")
        if rank != 10:
            failures.append("normalized ten-coordinate Jacobian is rank deficient")
        if hits:
            failures.append("root touches a declared coordinate boundary")
        if not exact_zero:
            failures.append("deleted transport is not exact zero")
        if np.any(np.asarray(state[:CORE_SIZE]) <= 0.0):
            failures.append("core state is not strictly positive")
        root = IndependentGenotypeRoot(
            root_id="",
            transporter=transporter,
            expression=expression,
            member_start_ids=(start_id,),
            coordinates=tuple(float(value) for value in coordinates),
            core_state=tuple(float(value) for value in state[:CORE_SIZE]),
            max_abs_scaled_residual=scaled_max,
            max_abs_omitted_charge_rhs_fmol_s=omitted_max,
            max_abs_bulk_charge_residual_fmol=charge_max,
            max_abs_current_residual_A=current_max,
            max_abs_regulatory_rhs_s_inv=regulatory_max,
            jacobian_singular_values=tuple(float(value) for value in singular),
            jacobian_rank=rank,
            jacobian_nullity=10 - rank,
            boundary_hits=hits,
            exact_deleted_transport_zero=exact_zero,
            passes_numerical_gate=not failures,
            gate_failures=tuple(failures),
        )
    except (
        ValueError,
        FloatingPointError,
        OverflowError,
        ZeroDivisionError,
        np.linalg.LinAlgError,
    ):
        admissible = False
    attempt = IndependentGenotypeAttempt(
        start_id=start_id,
        optimizer_success=bool(fit.success),
        admissible_state=admissible,
        passes_numerical_gate=bool(root and root.passes_numerical_gate),
        max_abs_scaled_residual=scaled_max,
        nfev=int(fit.nfev),
        cost=float(fit.cost),
        optimality=float(fit.optimality),
        message=str(fit.message),
        coordinates=tuple(float(value) for value in coordinates),
    )
    return attempt, root if root is not None and root.passes_numerical_gate else None


def _cluster_roots(
    roots: Sequence[IndependentGenotypeRoot],
    *,
    transporter: str,
    expression: float,
) -> tuple[IndependentGenotypeRoot, ...]:
    span = COORDINATE_UPPER - COORDINATE_LOWER
    ordered = tuple(
        sorted(
            roots,
            key=lambda item: (
                tuple(float(value) for value in item.coordinates),
                float(item.max_abs_scaled_residual),
                tuple(item.member_start_ids),
            ),
        )
    )
    points = tuple(np.asarray(root.coordinates, dtype=float) / span for root in ordered)
    if any(point.shape != span.shape or np.any(~np.isfinite(point)) for point in points):
        raise ValueError("independent root coordinates must be finite ten-vectors")
    parents = list(range(len(ordered)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)

    tolerance = 2.0e-5
    for left in range(len(points)):
        for right in range(left + 1, len(points)):
            if float(np.max(np.abs(points[left] - points[right]))) <= tolerance:
                union(left, right)
    grouped: dict[int, list[IndependentGenotypeRoot]] = {}
    for index, root in enumerate(ordered):
        grouped.setdefault(find(index), []).append(root)
    clusters = [
        grouped[key]
        for key in sorted(
            grouped,
            key=lambda component: tuple(
                float(value) for value in grouped[component][0].coordinates
            ),
        )
    ]
    retained: list[IndependentGenotypeRoot] = []
    for index, cluster in enumerate(clusters):
        best = min(
            cluster,
            key=lambda item: (
                item.max_abs_scaled_residual,
                tuple(item.coordinates),
                tuple(item.member_start_ids),
            ),
        )
        retained.append(
            IndependentGenotypeRoot(
                **{
                    **best.__dict__,
                    "root_id": f"{transporter}_E{expression:.8f}:IB{index:02d}",
                    "member_start_ids": tuple(
                        sorted(
                            start_id
                            for member in cluster
                            for start_id in member.member_start_ids
                        )
                    ),
                }
            )
        )
    return tuple(retained)


def continue_independent_genotype_expression(
    model: IndependentWholeCell,
    wt_resting_state: ArrayLike,
    *,
    transporter: str,
    expression_grid: Sequence[float] = EXPRESSION_GRID,
    local_start_fraction: float = 1.0e-3,
    maximum_connected_step_fraction: float = 0.25,
    max_nfev: int = 3000,
) -> IndependentGenotypeContinuation:
    """Continue one fixed native root without fitting or target scoring."""

    name = _checked_transporter(transporter)
    grid = tuple(float(value) for value in expression_grid)
    if len(grid) < 2 or grid[0] != 1.0 or grid[-1] != 0.0:
        raise ValueError("expression grid must begin at one and end at exact zero")
    if any(right >= left for left, right in zip(grid, grid[1:])):
        raise ValueError("expression grid must be strictly decreasing")
    if not 0.0 < local_start_fraction < 0.05:
        raise ValueError("local start fraction must lie in (0, 0.05)")
    if not 0.0 < maximum_connected_step_fraction <= 1.0:
        raise ValueError("maximum connected step fraction must lie in (0, 1]")

    core = np.asarray(wt_resting_state, dtype=float)
    if core.shape != (CORE_SIZE,):
        raise ValueError("WT resting root must contain exactly twelve core states")
    full = np.r_[core, np.asarray(_basal_suffix(model), dtype=float)]
    previous = independent_coordinates_from_state(model, full)
    if np.any(previous <= COORDINATE_LOWER) or np.any(previous >= COORDINATE_UPPER):
        raise ValueError("WT resting root lies outside independent continuation bounds")
    span = COORDINATE_UPPER - COORDINATE_LOWER
    interior = 1.0e-10 * span
    steps: list[IndependentGenotypeStep] = []
    failure: str | None = None

    for expression in grid:
        starts = (
            ("previous_connected_root", previous.copy()),
            (
                "local_minus",
                np.clip(
                    previous - local_start_fraction * span,
                    COORDINATE_LOWER + interior,
                    COORDINATE_UPPER - interior,
                ),
            ),
            (
                "local_plus",
                np.clip(
                    previous + local_start_fraction * span,
                    COORDINATE_LOWER + interior,
                    COORDINATE_UPPER - interior,
                ),
            ),
        )
        attempts: list[IndependentGenotypeAttempt] = []
        admissible_roots: list[IndependentGenotypeRoot] = []
        for start_id, start in starts:
            attempt, root = _one_attempt(
                model,
                transporter=name,
                expression=expression,
                start_id=start_id,
                start=start,
                max_nfev=max_nfev,
            )
            attempts.append(attempt)
            if root is not None:
                admissible_roots.append(root)
        roots = _cluster_roots(
            admissible_roots,
            transporter=name,
            expression=expression,
        )
        if not roots:
            failure = f"no admissible independent root at {name} expression {expression:.8f}"
            steps.append(
                IndependentGenotypeStep(
                    expression=expression,
                    attempts=tuple(attempts),
                    roots=(),
                    selected_root_id=None,
                    selected_distance_from_previous_normalized=None,
                    connected_to_previous=False,
                    status="NO_ADMISSIBLE_ROOT",
                )
            )
            break
        distances = {
            root.root_id: float(
                np.max(
                    np.abs(np.asarray(root.coordinates, dtype=float) - previous) / span
                )
            )
            for root in roots
        }
        selected = min(
            roots,
            key=lambda item: (distances[item.root_id], item.max_abs_scaled_residual),
        )
        distance = distances[selected.root_id]
        connected = bool(distance <= maximum_connected_step_fraction)
        if not connected:
            status = "ROOT_FOUND_BUT_BRANCH_JUMP_FLAGGED"
            failure = (
                f"nearest independent root jump {distance:.6g} exceeds limit "
                f"at {name} expression {expression:.8f}"
            )
        elif len(roots) > 1:
            status = "CONNECTED_WITH_ALTERNATE_ROOTS_RETAINED"
        else:
            status = "CONNECTED"
        steps.append(
            IndependentGenotypeStep(
                expression=expression,
                attempts=tuple(attempts),
                roots=roots,
                selected_root_id=selected.root_id,
                selected_distance_from_previous_normalized=distance,
                connected_to_previous=connected,
                status=status,
            )
        )
        if not connected:
            break
        previous = np.asarray(selected.coordinates, dtype=float)

    completed = bool(
        len(steps) == len(grid)
        and steps[-1].expression == 0.0
        and steps[-1].connected_to_previous
        and steps[-1].selected_root is not None
        and steps[-1].selected_root.exact_deleted_transport_zero
    )
    return IndependentGenotypeContinuation(
        transporter=name,
        expression_grid=grid,
        steps=tuple(steps),
        completed_to_exact_zero=completed,
        branch_failure=None if completed else failure or "continuation incomplete",
    )


__all__ = (
    "COORDINATE_LOWER",
    "COORDINATE_NAMES",
    "COORDINATE_UPPER",
    "EXPRESSION_GRID",
    "INDEPENDENT_ROW_SCALES",
    "IndependentGenotypeAttempt",
    "IndependentGenotypeContinuation",
    "IndependentGenotypeRoot",
    "IndependentGenotypeStep",
    "TRANSPORTERS",
    "continue_independent_genotype_expression",
    "independent_coordinates_from_state",
)
