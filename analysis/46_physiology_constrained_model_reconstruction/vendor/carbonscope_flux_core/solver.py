"""Small self-contained HiGHS FBA/FVA/VFFVA core for Task 46.

This is a deliberately reduced local port of the generic flux-analysis ideas in
`esig626/CarbonScope` commit a11e25f176a2cc70de22c3f36a1c0e623248393a.
It is not the CarbonScope package and has no network or write dependency on that
repository.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from queue import Empty, Queue
import threading
from typing import Sequence

import numpy as np
import pandas as pd

from .exceptions import AnalysisError
from .results import FBAResult, FVAResult, PrimalDiagnostics
from .schema import FluxModel
from .validation import CanonicalModelError, validate_flux_model

FEASIBILITY_TOLERANCE = 1e-7
OBJECTIVE_TOLERANCE = 1e-7
SOLVER_FEASIBILITY_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class CompiledFluxLP:
    reaction_ids: tuple[str, ...]
    balanced_metabolite_ids: tuple[str, ...]
    balance_matrix: tuple[tuple[float, ...], ...]
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]
    objective: tuple[float, ...]
    objective_direction: str


@dataclass(frozen=True, slots=True)
class RetainedObjective:
    fraction: float
    optimum: float
    bound: float
    sense: str


@dataclass(frozen=True, slots=True)
class PreparedFluxRegion:
    lp: CompiledFluxLP
    fba: FBAResult
    retention: RetainedObjective | None


def _highspy():
    try:
        import highspy
    except ImportError as error:
        raise AnalysisError(
            "Task 46 CBM requires highspy; install highspy>=1.11,<1.13"
        ) from error
    return highspy


def compile_flux_lp(model: FluxModel) -> CompiledFluxLP:
    try:
        validate_flux_model(model)
    except (CanonicalModelError, AttributeError) as error:
        raise AnalysisError(f"invalid flux model: {error}") from error
    reaction_ids = tuple(r.reaction_id for r in model.reactions)
    balanced = tuple(
        m.metabolite_id for m in model.metabolites if m.steady_state_balanced
    )
    rindex = {rid: i for i, rid in enumerate(reaction_ids)}
    mindex = {mid: i for i, mid in enumerate(balanced)}
    matrix = np.zeros((len(balanced), len(reaction_ids)), dtype=float)
    for j, reaction in enumerate(model.reactions):
        for term in reaction.stoichiometric_terms:
            row = mindex.get(term.metabolite_id)
            if row is not None:
                matrix[row, j] += float(term.coefficient)
    lower = tuple(float(r.lower_bound) for r in model.reactions)
    upper = tuple(float(r.upper_bound) for r in model.reactions)
    objective = [0.0] * len(reaction_ids)
    for term in model.objective.terms:
        objective[rindex[term.reaction_id]] += float(term.coefficient)
    direction = "max" if model.objective.direction == "maximise" else "min"
    return CompiledFluxLP(
        reaction_ids,
        balanced,
        tuple(tuple(float(v) for v in row) for row in matrix),
        lower,
        upper,
        tuple(objective),
        direction,
    )


def _csr_rows(matrix: np.ndarray) -> tuple[list[int], list[int], list[float]]:
    starts = [0]
    indices: list[int] = []
    values: list[float] = []
    for row in matrix:
        for j, value in enumerate(row):
            if value != 0.0:
                indices.append(j)
                values.append(float(value))
        starts.append(len(indices))
    return starts, indices, values


def _build_solver(
    lp: CompiledFluxLP,
    costs: Sequence[float],
    direction: str,
    retention: RetainedObjective | None = None,
):
    highspy = _highspy()
    solver = highspy.Highs()
    solver.setOptionValue("output_flag", False)
    solver.setOptionValue("threads", 1)
    solver.setOptionValue("solver", "simplex")
    solver.setOptionValue("primal_feasibility_tolerance", SOLVER_FEASIBILITY_TOLERANCE)
    solver.setOptionValue("dual_feasibility_tolerance", SOLVER_FEASIBILITY_TOLERANCE)
    n = len(lp.reaction_ids)
    solver.addCols(
        n,
        list(map(float, costs)),
        list(lp.lower_bounds),
        list(lp.upper_bounds),
        0,
        [0] * (n + 1),
        [],
        [],
    )
    matrix = np.asarray(lp.balance_matrix, dtype=float)
    starts, indices, values = _csr_rows(matrix)
    lower_rows = [0.0] * matrix.shape[0]
    upper_rows = [0.0] * matrix.shape[0]
    if retention is not None:
        row = np.asarray(lp.objective, dtype=float)
        nz = np.flatnonzero(row)
        indices.extend(int(i) for i in nz)
        values.extend(float(row[i]) for i in nz)
        starts.append(len(indices))
        if retention.sense == ">=":
            lower_rows.append(float(retention.bound))
            upper_rows.append(highspy.kHighsInf)
        else:
            lower_rows.append(-highspy.kHighsInf)
            upper_rows.append(float(retention.bound))
    solver.addRows(
        len(lower_rows),
        lower_rows,
        upper_rows,
        len(indices),
        starts,
        indices,
        values,
    )
    if direction == "max":
        solver.setMaximize()
    else:
        solver.setMinimize()
    return solver


def _objective_value(costs: Sequence[float], fluxes: Sequence[float]) -> float:
    return math.fsum(float(c) * float(v) for c, v in zip(costs, fluxes))


def _validate_solution(
    lp: CompiledFluxLP,
    fluxes: Sequence[float],
    reported_objective: float,
    costs: Sequence[float],
) -> PrimalDiagnostics:
    vector = np.asarray(fluxes, dtype=float)
    if vector.shape != (len(lp.reaction_ids),) or not np.isfinite(vector).all():
        raise AnalysisError("solver returned malformed/non-finite flux vector")
    lower = np.asarray(lp.lower_bounds, dtype=float)
    upper = np.asarray(lp.upper_bounds, dtype=float)
    lower_violation = float(np.max(np.maximum(lower - vector, 0.0), initial=0.0))
    upper_violation = float(np.max(np.maximum(vector - upper, 0.0), initial=0.0))
    matrix = np.asarray(lp.balance_matrix, dtype=float)
    residual = float(np.max(np.abs(matrix @ vector), initial=0.0)) if matrix.size else 0.0
    recalculated = _objective_value(costs, vector)
    objective_error = abs(float(reported_objective) - recalculated)
    diagnostics = PrimalDiagnostics(
        lower_violation,
        upper_violation,
        residual,
        objective_error,
    )
    if max(lower_violation, upper_violation, residual) > FEASIBILITY_TOLERANCE:
        raise AnalysisError(f"primal feasibility validation failed: {diagnostics}")
    if objective_error > OBJECTIVE_TOLERANCE * max(1.0, abs(recalculated)):
        raise AnalysisError(f"objective validation failed: {diagnostics}")
    return diagnostics


def _solve(
    lp: CompiledFluxLP,
    costs: Sequence[float],
    direction: str,
    *,
    retention: RetainedObjective | None = None,
    operation: str = "LP",
) -> tuple[list[float], float, str, PrimalDiagnostics]:
    highspy = _highspy()
    solver = _build_solver(lp, costs, direction, retention)
    solver.run()
    status = solver.getModelStatus()
    status_name = solver.modelStatusToString(status)
    if status != highspy.HighsModelStatus.kOptimal:
        raise AnalysisError(f"{operation} failed: HiGHS status {status_name}")
    fluxes = [float(v) for v in solver.getSolution().col_value]
    objective = float(solver.getObjectiveValue())
    diagnostics = _validate_solution(lp, fluxes, objective, costs)
    return fluxes, objective, status_name, diagnostics


def run_highs_fba(model: FluxModel) -> FBAResult:
    lp = compile_flux_lp(model)
    fluxes, objective, _, diagnostics = _solve(
        lp,
        lp.objective,
        lp.objective_direction,
        operation="FBA",
    )
    return FBAResult(
        objective,
        "optimal",
        lp.objective_direction,
        pd.Series(fluxes, index=lp.reaction_ids, dtype=float),
        diagnostics,
    )


def prepare_highs_flux_region(
    model: FluxModel,
    fraction_of_optimum: float | None = None,
) -> PreparedFluxRegion:
    lp = compile_flux_lp(model)
    fba = run_highs_fba(model)
    if fraction_of_optimum is None:
        return PreparedFluxRegion(lp, fba, None)
    if isinstance(fraction_of_optimum, bool):
        raise AnalysisError("fraction_of_optimum must be numeric")
    fraction = float(fraction_of_optimum)
    if not math.isfinite(fraction) or not 0.0 < fraction <= 1.0:
        raise AnalysisError("fraction_of_optimum must be in (0, 1]")
    if lp.objective_direction == "max":
        if fba.objective_value < 0 and fraction < 1:
            raise AnalysisError("fractional retention is ambiguous for negative maximum")
        retention = RetainedObjective(
            fraction,
            fba.objective_value,
            fba.objective_value * fraction,
            ">=",
        )
    else:
        if fba.objective_value > 0 and fraction < 1:
            raise AnalysisError("fractional retention is ambiguous for positive minimum")
        retention = RetainedObjective(
            fraction,
            fba.objective_value,
            fba.objective_value * fraction,
            "<=",
        )
    return PreparedFluxRegion(lp, fba, retention)


def _validate_retention(
    lp: CompiledFluxLP,
    fluxes: Sequence[float],
    retention: RetainedObjective | None,
) -> None:
    if retention is None:
        return
    value = _objective_value(lp.objective, fluxes)
    violation = (
        max(retention.bound - value, 0.0)
        if retention.sense == ">="
        else max(value - retention.bound, 0.0)
    )
    if violation > OBJECTIVE_TOLERANCE * max(1.0, abs(retention.bound)):
        raise AnalysisError(
            f"FVA endpoint violates retained objective: value={value:g}, "
            f"bound={retention.bound:g}"
        )


def run_highs_fva_reference(
    model: FluxModel,
    fraction_of_optimum: float | None = None,
) -> FVAResult:
    prepared = prepare_highs_flux_region(model, fraction_of_optimum)
    lp = prepared.lp
    minima: list[float] = []
    maxima: list[float] = []
    n = len(lp.reaction_ids)
    for j, reaction_id in enumerate(lp.reaction_ids):
        costs = [0.0] * n
        costs[j] = 1.0
        low_flux, low, _, _ = _solve(
            lp,
            costs,
            "min",
            retention=prepared.retention,
            operation=f"FVA minimum {reaction_id}",
        )
        _validate_retention(lp, low_flux, prepared.retention)
        high_flux, high, _, _ = _solve(
            lp,
            costs,
            "max",
            retention=prepared.retention,
            operation=f"FVA maximum {reaction_id}",
        )
        _validate_retention(lp, high_flux, prepared.retention)
        minima.append(float(low))
        maxima.append(float(high))
    ranges = pd.DataFrame(
        {"minimum": minima, "maximum": maxima},
        index=lp.reaction_ids,
    )
    return FVAResult(
        ranges,
        None if prepared.retention is None else prepared.retention.fraction,
        prepared.fba.objective_value,
        lp.objective_direction,
    )


class _ReusableWorker:
    def __init__(self, prepared: PreparedFluxRegion):
        self.prepared = prepared
        self.lp = prepared.lp
        self.solver = _build_solver(
            self.lp,
            [0.0] * len(self.lp.reaction_ids),
            "max",
            prepared.retention,
        )

    def begin_pass(self, direction: str) -> None:
        if direction == "max":
            self.solver.setMaximize()
        else:
            self.solver.setMinimize()

    def solve_endpoint(self, j: int, direction: str) -> float:
        highspy = _highspy()
        self.solver.changeColCost(j, 1.0)
        try:
            self.solver.run()
            status = self.solver.getModelStatus()
            if status != highspy.HighsModelStatus.kOptimal:
                raise AnalysisError(
                    f"VFFVA {direction} endpoint {self.lp.reaction_ids[j]} failed: "
                    f"{self.solver.modelStatusToString(status)}"
                )
            value = float(self.solver.getObjectiveValue())
            fluxes = [float(v) for v in self.solver.getSolution().col_value]
            costs = [0.0] * len(self.lp.reaction_ids)
            costs[j] = 1.0
            _validate_solution(self.lp, fluxes, value, costs)
            _validate_retention(self.lp, fluxes, self.prepared.retention)
            return value
        finally:
            self.solver.changeColCost(j, 0.0)


def run_prepared_highs_vffva(
    prepared: PreparedFluxRegion,
    *,
    workers: int | None = None,
    chunk_size: int = 16,
) -> FVAResult:
    """Shared-memory FVA with one reusable single-thread HiGHS model per worker."""

    n = len(prepared.lp.reaction_ids)
    if workers is None:
        workers = min(4, n)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise AnalysisError("workers must be a positive integer")
    workers = min(workers, n)
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or chunk_size < 1:
        raise AnalysisError("chunk_size must be a positive integer")

    minima = [math.nan] * n
    maxima = [math.nan] * n
    failure: list[BaseException] = []
    failure_lock = threading.Lock()

    def execute_pass(direction: str, destination: list[float]) -> None:
        queue: Queue[tuple[int, ...]] = Queue()
        for start in range(0, n, chunk_size):
            queue.put(tuple(range(start, min(start + chunk_size, n))))

        def target() -> None:
            worker = _ReusableWorker(prepared)
            worker.begin_pass(direction)
            while True:
                try:
                    chunk = queue.get_nowait()
                except Empty:
                    return
                try:
                    if failure:
                        continue
                    for j in chunk:
                        if failure:
                            break
                        destination[j] = worker.solve_endpoint(j, direction)
                except BaseException as error:
                    with failure_lock:
                        if not failure:
                            failure.append(error)
                finally:
                    queue.task_done()

        threads = [threading.Thread(target=target) for _ in range(workers)]
        for thread in threads:
            thread.start()
        queue.join()
        for thread in threads:
            thread.join()
        if failure:
            raise AnalysisError(str(failure[0])) from None

    execute_pass("max", maxima)
    execute_pass("min", minima)
    if not all(math.isfinite(v) for v in minima + maxima):
        raise AnalysisError("VFFVA failed to return all endpoints")
    for j, (low, high) in enumerate(zip(minima, maxima)):
        if low > high + 1e-9:
            raise AnalysisError(
                f"VFFVA endpoint inversion for {prepared.lp.reaction_ids[j]}: "
                f"{low:g}>{high:g}"
            )
    ranges = pd.DataFrame(
        {"minimum": minima, "maximum": maxima},
        index=prepared.lp.reaction_ids,
    )
    return FVAResult(
        ranges,
        None if prepared.retention is None else prepared.retention.fraction,
        prepared.fba.objective_value,
        prepared.lp.objective_direction,
    )


def run_highs_vffva(
    model: FluxModel,
    fraction_of_optimum: float | None = None,
    *,
    workers: int | None = None,
    chunk_size: int = 16,
) -> FVAResult:
    return run_prepared_highs_vffva(
        prepare_highs_flux_region(model, fraction_of_optimum),
        workers=workers,
        chunk_size=chunk_size,
    )


def assert_fva_equivalent(
    reference: FVAResult,
    fast: FVAResult,
    *,
    atol: float = 1e-8,
    rtol: float = 1e-8,
) -> None:
    if tuple(reference.ranges.index) != tuple(fast.ranges.index):
        raise AnalysisError("FVA result reaction ordering differs")
    if tuple(reference.ranges.columns) != tuple(fast.ranges.columns):
        raise AnalysisError("FVA result columns differ")
    if not np.allclose(
        reference.ranges.to_numpy(dtype=float),
        fast.ranges.to_numpy(dtype=float),
        atol=atol,
        rtol=rtol,
    ):
        delta = np.max(
            np.abs(
                reference.ranges.to_numpy(dtype=float)
                - fast.ranges.to_numpy(dtype=float)
            )
        )
        raise AnalysisError(f"reference FVA and VFFVA differ; max abs delta={delta:g}")
