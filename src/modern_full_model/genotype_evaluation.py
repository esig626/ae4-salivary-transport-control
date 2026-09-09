"""Target-firewalled genotype continuation and dynamic evaluation.

This module is the genotype-side counterpart of the WT-only calibration and
dynamic-contract modules.  It deliberately contains no experimental
knockout-secretion values, target intervals, or target-ledger readers.  Its
only comparisons are model-to-model ratios formed *after* a WT model and its
observation scale have been frozen elsewhere.

Resting deletion is not approximated by evaluating a null vector field at the
WT state.  Instead, transporter expression is continued from one to zero on
the same exact cell/lumen bulk-charge manifolds used for WT calibration.  At
each expression value ten independent steady-state equations are solved for
the ten physical charge-manifold coordinates.  All model parameters,
stimulus inputs, and the finite intracellular OTHER osmole amount remain
unchanged.  Local multistarts are clustered and retained; the branch nearest
the preceding root is selected without reference to a phenotype.

The frozen production dynamic helper instantiates the predeclared
``R3_G125_P21_PROSE_TREFERENCE`` member, the 0.10-uM effective calcium
sensitivity, and the algebraic 1.75-fold NKCC1 capacity response.  Those
choices are inherited WT-only construction decisions.  They are not inferred
from a genotype trajectory.  AE4 deletion leaves the cAMP/PKA regulatory
states present but multiplies every AE4 conserved source by exact zero.

Provenance
----------
``DERIVED_CONSTRAINT``
    Charge-manifold elimination, ten-row independence, exact-zero deletion,
    branch-distance normalization, and model-to-model ratio algebra.
``PRIMARY_MEASUREMENT``
    The 600-s CCh+IPR protocol duration and separate muscarinic/beta arms.
``NEW_MODELING_DECISION``
    Expression grid, local continuation starts, numerical tolerances, the
    effective calcium sensitivity, and the unmeasured regulatory time scales.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

from .calibration import (
    CORE_INDEPENDENT_ROWS,
    CORE_OMITTED_CHARGE_ROWS,
    REST_COORDINATE_NAMES,
    WTCalibrationSpec,
    decode_rest_coordinates,
    encode_electroneutral_coordinates,
)
from .model import Genotype, ModernFullModel, WT
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .states import CORE_STATE_NAMES
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    SolverSpecification,
    StimulusArm,
    attach_basal_regulation,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
)


FROZEN_REGULATORY_MEMBER_ID = "R3_G125_P21_PROSE_TREFERENCE"
FROZEN_EFFECTIVE_CALCIUM_UM = 0.10
FROZEN_NKCC1_FULLY_ACTIVATED_MULTIPLIER = 1.75
GENOTYPE_SOLVER_RELATIVE_TOLERANCE = 1.0e-4
DEFAULT_EXPRESSION_GRID = tuple(float(value) for value in np.linspace(1.0, 0.0, 21))
SUPPORTED_TRANSPORTERS = ("AE4", "AE2")


@dataclass(frozen=True)
class ContinuationAttempt:
    """One bounded ten-coordinate root attempt at fixed expression."""

    start_id: str
    optimizer_success: bool
    admissible_state: bool
    converged_root: bool
    max_abs_scaled_independent_rhs: float
    cost: float
    optimality: float
    nfev: int
    message: str
    coordinates: tuple[float, ...]


@dataclass(frozen=True)
class GenotypeRoot:
    """One retained mathematical root at a fixed expression value."""

    root_id: str
    transporter: str
    expression: float
    member_start_ids: tuple[str, ...]
    coordinates: tuple[float, ...]
    core_state: tuple[float, ...]
    max_abs_scaled_independent_rhs: float
    max_abs_amount_rhs_fmol_s: float
    max_abs_volume_rhs_pL_s: float
    max_abs_omitted_charge_rhs_fmol_equivalent_s: float
    max_abs_regulatory_rhs_s_inv: float
    max_abs_charge_residual_fmol: float
    max_abs_current_residual_A: float
    normalized_jacobian_singular_values: tuple[float, ...]
    normalized_jacobian_rank: int
    normalized_jacobian_nullity: int
    boundary_hits: tuple[str, ...]
    exact_deleted_transport_zero: bool
    passes_numerical_gate: bool
    gate_failures: tuple[str, ...]


@dataclass(frozen=True)
class ContinuationStep:
    """All attempts and retained roots at one expression value."""

    expression: float
    attempts: tuple[ContinuationAttempt, ...]
    roots: tuple[GenotypeRoot, ...]
    selected_root_id: str | None
    selected_distance_from_previous_normalized: float | None
    connected_to_previous: bool
    alternate_root_count: int
    status: str

    @property
    def selected_root(self) -> GenotypeRoot | None:
        if self.selected_root_id is None:
            return None
        return next(root for root in self.roots if root.root_id == self.selected_root_id)


@dataclass(frozen=True)
class GenotypeContinuation:
    """Connected expression branch with immutable model-construction metadata."""

    transporter: str
    expression_grid: tuple[float, ...]
    fixed_other_impermeant_osmoles_fmol: float
    independent_equation_count: int
    fitted_coordinate_count: int
    fitted_parameter_count: int
    steps: tuple[ContinuationStep, ...]
    completed_to_exact_zero: bool
    branch_failure: str | None

    @property
    def connected_null_root(self) -> GenotypeRoot | None:
        if not self.completed_to_exact_zero or not self.steps:
            return None
        root = self.steps[-1].selected_root
        if root is None or root.expression != 0.0:
            return None
        return root


@dataclass(frozen=True)
class GenotypeTrajectory:
    """One target-free WT or genotype dynamic trajectory."""

    genotype_name: str
    transporter: str
    solver: SolverSpecification
    success: bool
    message: str
    time_s: NDArray[np.float64]
    state_names: tuple[str, ...]
    states: NDArray[np.float64]
    flow_pL_s: NDArray[np.float64]
    cumulative_flow_pL: NDArray[np.float64]
    cell_na_mM: NDArray[np.float64]
    cell_k_mM: NDArray[np.float64]
    cell_cl_mM: NDArray[np.float64]
    cell_ph: NDArray[np.float64]
    cell_volume_pL: NDArray[np.float64]
    lumen_na_mM: NDArray[np.float64]
    lumen_k_mM: NDArray[np.float64]
    lumen_cl_mM: NDArray[np.float64]
    lumen_ph: NDArray[np.float64]
    lumen_volume_pL: NDArray[np.float64]
    ae4_capacity_multiplier: NDArray[np.float64]
    max_abs_deleted_transport_flux_fmol_s: float
    max_abs_conservation_residuals: Mapping[str, float]
    max_dimensionless_conservation_ratio: float
    positive_core: bool
    all_flow_nonnegative: bool
    numerical_gate_pass: bool


@dataclass(frozen=True)
class GenotypePairSummary:
    """Generic genotype/WT ratios with no experimental acceptance criteria."""

    transporter: str
    solver_label: str
    wt_success: bool
    genotype_success: bool
    paired_numerical_gate_pass: bool
    cumulative_ratio_at_landmarks: Mapping[str, float]
    flow_ratio_at_landmarks: Mapping[str, float]
    wt_endpoint: Mapping[str, float]
    genotype_endpoint: Mapping[str, float]
    genotype_to_wt_endpoint_ratios: Mapping[str, float]
    wt_max_dimensionless_conservation_ratio: float
    genotype_max_dimensionless_conservation_ratio: float


@dataclass(frozen=True)
class DualSolverGenotypeEvaluation:
    """Radau/BDF paired evaluations and numerical crosschecks."""

    transporter: str
    evaluations: Mapping[str, GenotypePairSummary]
    wt_trajectories: Mapping[str, GenotypeTrajectory]
    genotype_trajectories: Mapping[str, GenotypeTrajectory]
    cross_solver_relative_differences: Mapping[str, float]
    cross_solver_relative_tolerance: float
    cross_solver_gate_pass: bool
    both_solver_pairs_pass: bool


def _checked_transporter(transporter: str) -> str:
    name = str(transporter).upper()
    if name not in SUPPORTED_TRANSPORTERS:
        raise ValueError(f"transporter must be one of {SUPPORTED_TRANSPORTERS}")
    return name


def genotype_with_expression(transporter: str, expression: float) -> Genotype:
    """Return a single-transporter expression perturbation on a WT background."""

    name = _checked_transporter(transporter)
    value = float(expression)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("expression must be finite and lie in [0, 1]")
    label = f"{name}_EXPRESSION_{value:.8f}"
    if name == "AE4":
        return Genotype(label, ae4_expression=value)
    return Genotype(label, ae2_expression=value)


def frozen_regulatory_member() -> Any:
    """Return the exact predeclared R3/G1.25/reference-timing member."""

    matches = tuple(
        member
        for member in pre_reveal_regulatory_ensemble()
        if member.member_id == FROZEN_REGULATORY_MEMBER_ID
    )
    if len(matches) != 1:
        raise AssertionError("frozen regulatory member must occur exactly once")
    member = matches[0]
    if member.family != "R3" or member.fully_activated_multiplier != 1.25:
        raise AssertionError("frozen regulatory member metadata changed")
    return member


def build_frozen_genotype_dynamic_model(template: Any) -> Any:
    """Clone one frozen structural row into the genotype dynamic contract.

    ``template`` may be a plain :class:`ModernFullModel` or a wrapped model;
    only its already-frozen parameters and AE4 transport module are reused.
    A new wrapper is built so the selected regulatory and NKCC1 contracts are
    explicit and cannot be inherited accidentally from an older generation.
    """

    member = frozen_regulatory_member()
    protocol = SecretagogueProtocol(
        arm=StimulusArm.CCH_IPR,
        stimulated_calcium_uM=FROZEN_EFFECTIVE_CALCIUM_UM,
    )
    base = ModernFullModel(
        parameters=template.parameters,
        stimulus=protocol,
        regulatory_model=member.regulatory_model,
        ae4_parameters=template.ae4_parameters,
        ae4_evaluator=template.ae4_evaluator,
    )
    nkcc1 = N1AlgebraicNkcc1(
        fully_activated_multiplier=FROZEN_NKCC1_FULLY_ACTIVATED_MULTIPLIER,
        resting_calcium_uM=protocol.resting_calcium_uM,
        stimulated_calcium_uM=protocol.stimulated_calcium_uM,
    )
    model = attach_stimulated_nkcc1(base, nkcc1)
    before = float(template.parameters.geometry.cell_impermeant_osmoles_fmol)
    after = float(model.parameters.geometry.cell_impermeant_osmoles_fmol)
    if before != after:
        raise AssertionError("dynamic cloning changed frozen OTHER osmoles")
    return model


def frozen_dynamic_contract_metadata() -> Mapping[str, Any]:
    """Return target-free, machine-readable genotype protocol metadata."""

    member = frozen_regulatory_member()
    return {
        "regulatory_member_id": member.member_id,
        "regulatory_family": member.family,
        "regulatory_gain_label": member.gain_label,
        "regulatory_fully_activated_multiplier": member.fully_activated_multiplier,
        "regulatory_timing_label": member.timing_label,
        "regulatory_timing_parameters": dict(member.timing_parameters),
        "effective_calcium_uM": FROZEN_EFFECTIVE_CALCIUM_UM,
        "nkcc1_family": "N1_CCH_ALGEBRAIC",
        "nkcc1_fully_activated_multiplier": (
            FROZEN_NKCC1_FULLY_ACTIVATED_MULTIPLIER
        ),
        "stimulus_arm": StimulusArm.CCH_IPR.value,
        "duration_s": 600.0,
        "calibration_firewall": (
            "all values frozen from WT/source construction before genotype evaluation"
        ),
    }


def _expression_grid(values: Sequence[float]) -> tuple[float, ...]:
    grid = tuple(float(value) for value in values)
    if len(grid) < 2 or grid[0] != 1.0 or grid[-1] != 0.0:
        raise ValueError("expression grid must begin at one and end at exact zero")
    if any(not math.isfinite(value) or not 0.0 <= value <= 1.0 for value in grid):
        raise ValueError("expression grid values must be finite and lie in [0, 1]")
    if any(right >= left for left, right in zip(grid, grid[1:])):
        raise ValueError("expression grid must be strictly decreasing")
    return grid


def _root_state_with_basal_regulation(model: Any, core_state: ArrayLike) -> NDArray[np.float64]:
    candidate = np.asarray(core_state, dtype=float)
    if candidate.shape not in ((len(CORE_STATE_NAMES),), (model.layout.size,)):
        raise ValueError("WT root must contain the twelve-state core or full model state")
    return attach_basal_regulation(model, candidate)


def _core_coordinates(model: Any, resting_state: ArrayLike) -> NDArray[np.float64]:
    full = _root_state_with_basal_regulation(model, resting_state)
    coordinates = decode_rest_coordinates(model, full)
    return np.asarray(coordinates[: len(REST_COORDINATE_NAMES)], dtype=float)


def _regulatory_suffix(model: Any) -> NDArray[np.float64]:
    initial = np.asarray(model.initial_state(), dtype=float)
    return initial[len(CORE_STATE_NAMES) :].copy()


def _state_from_core_coordinates(
    model: Any,
    coordinates: ArrayLike,
    regulatory_suffix: ArrayLike,
) -> NDArray[np.float64]:
    core = np.asarray(coordinates, dtype=float)
    suffix = np.asarray(regulatory_suffix, dtype=float)
    return encode_electroneutral_coordinates(model, np.r_[core, suffix])


def _scaled_residual(
    model: Any,
    coordinates: NDArray[np.float64],
    regulatory_suffix: NDArray[np.float64],
    genotype: Genotype,
    spec: WTCalibrationSpec,
) -> NDArray[np.float64]:
    state = _state_from_core_coordinates(model, coordinates, regulatory_suffix)
    raw = model.rhs(0.0, state, genotype=genotype)
    return raw[np.asarray(CORE_INDEPENDENT_ROWS, dtype=int)] / np.asarray(
        spec.independent_rhs_scales, dtype=float
    )


def _unit_maxima(raw_rhs: ArrayLike) -> tuple[float, float, float, float]:
    raw = np.asarray(raw_rhs, dtype=float)
    amount_rows = np.asarray((0, 1, 2, 3, 6, 7, 8, 9), dtype=int)
    volume_rows = np.asarray((5, 11), dtype=int)
    amount = float(np.max(np.abs(raw[amount_rows])))
    volume = float(np.max(np.abs(raw[volume_rows])))
    regulatory = float(np.max(np.abs(raw[12:]))) if len(raw) > 12 else 0.0
    omitted = float(
        np.max(np.abs(raw[np.asarray(CORE_OMITTED_CHARGE_ROWS, dtype=int)]))
    )
    return amount, volume, regulatory, omitted


def _boundary_hits(
    coordinates: NDArray[np.float64],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    relative_tolerance: float,
) -> tuple[str, ...]:
    relative_distance = np.minimum(coordinates - lower, upper - coordinates) / (
        upper - lower
    )
    return tuple(
        name
        for name, distance in zip(REST_COORDINATE_NAMES, relative_distance)
        if distance <= relative_tolerance
    )


def _deleted_transport_is_exact_zero(
    evaluation: Any, transporter: str, expression: float
) -> bool:
    if expression != 0.0:
        return True
    if transporter == "AE4":
        contribution = evaluation.diagnostics.ae4
        values = (
            contribution.na_cell_fmol_s,
            contribution.k_cell_fmol_s,
            contribution.cl_cell_fmol_s,
            contribution.hco3_cell_fmol_s,
            contribution.transported_charge_fmol_s,
        )
        return all(float(value) == 0.0 for value in values)
    return float(evaluation.diagnostics.homeostasis.ae2_inward_fmol_s) == 0.0


def _attempt_and_root(
    *,
    model: Any,
    transporter: str,
    expression: float,
    genotype: Genotype,
    regulatory_suffix: NDArray[np.float64],
    spec: WTCalibrationSpec,
    start_id: str,
    start: NDArray[np.float64],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    max_nfev: int,
) -> tuple[ContinuationAttempt, GenotypeRoot | None]:
    span = upper - lower

    def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            return _scaled_residual(
                model, candidate, regulatory_suffix, genotype, spec
            )
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(len(CORE_INDEPENDENT_ROWS), 1.0e6)

    fit = least_squares(
        residual,
        start,
        bounds=(lower, upper),
        max_nfev=max_nfev,
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        x_scale="jac",
    )
    coordinates = np.asarray(fit.x, dtype=float)
    scaled_max = float(np.max(np.abs(residual(coordinates))))
    admissible = True
    root: GenotypeRoot | None = None
    try:
        state = _state_from_core_coordinates(model, coordinates, regulatory_suffix)
        evaluation = model.evaluate(0.0, state, genotype=genotype)
        amount_max, volume_max, regulatory_max, omitted_max = _unit_maxima(
            evaluation.rhs
        )
        charge_max = max(
            abs(float(value))
            for value in evaluation.diagnostics.state_charge_fmol.values()
        )
        current_max = max(
            abs(float(value))
            for value in evaluation.diagnostics.membranes.current_residuals_A.values()
        )
        normalized_jacobian = np.asarray(fit.jac, dtype=float) * span[None, :]
        singular = np.linalg.svd(normalized_jacobian, compute_uv=False)
        threshold = (
            spec.jacobian_rank_relative_tolerance * singular[0]
            if singular.size
            else math.inf
        )
        rank = int(np.sum(singular > threshold))
        hits = _boundary_hits(
            coordinates, lower, upper, spec.boundary_relative_tolerance
        )
        exact_zero = _deleted_transport_is_exact_zero(
            evaluation, transporter, expression
        )
        failures: list[str] = []
        if not fit.success:
            failures.append("optimizer did not report convergence")
        if scaled_max > spec.root_scaled_tolerance:
            failures.append("independent steady-state residual exceeds tolerance")
        if omitted_max > spec.omitted_row_raw_tolerance:
            failures.append("omitted charge-dependent row exceeds tolerance")
        if charge_max > spec.charge_tolerance_fmol:
            failures.append("bulk-charge manifold residual exceeds tolerance")
        if current_max > spec.current_tolerance_A:
            failures.append("membrane current closure exceeds tolerance")
        if regulatory_max > 1.0e-10:
            failures.append("basal regulatory state is not stationary")
        if rank != len(REST_COORDINATE_NAMES):
            failures.append("normalized ten-by-ten Jacobian is rank deficient")
        if hits:
            failures.append("root touches a declared search boundary")
        if not exact_zero:
            failures.append("deleted transporter contribution is not exact zero")
        if not np.all(np.asarray(state[: len(CORE_STATE_NAMES)]) > 0.0):
            failures.append("core state is not strictly positive")
        root = GenotypeRoot(
            root_id="",
            transporter=transporter,
            expression=expression,
            member_start_ids=(start_id,),
            coordinates=tuple(float(value) for value in coordinates),
            core_state=tuple(float(value) for value in state[: len(CORE_STATE_NAMES)]),
            max_abs_scaled_independent_rhs=scaled_max,
            max_abs_amount_rhs_fmol_s=amount_max,
            max_abs_volume_rhs_pL_s=volume_max,
            max_abs_omitted_charge_rhs_fmol_equivalent_s=omitted_max,
            max_abs_regulatory_rhs_s_inv=regulatory_max,
            max_abs_charge_residual_fmol=charge_max,
            max_abs_current_residual_A=current_max,
            normalized_jacobian_singular_values=tuple(
                float(value) for value in singular
            ),
            normalized_jacobian_rank=rank,
            normalized_jacobian_nullity=len(REST_COORDINATE_NAMES) - rank,
            boundary_hits=hits,
            exact_deleted_transport_zero=exact_zero,
            passes_numerical_gate=not failures,
            gate_failures=tuple(failures),
        )
    except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
        admissible = False
    attempt = ContinuationAttempt(
        start_id=start_id,
        optimizer_success=bool(fit.success),
        admissible_state=admissible,
        converged_root=bool(
            admissible and scaled_max <= spec.root_scaled_tolerance
        ),
        max_abs_scaled_independent_rhs=scaled_max,
        cost=float(fit.cost),
        optimality=float(fit.optimality),
        nfev=int(fit.nfev),
        message=str(fit.message),
        coordinates=tuple(float(value) for value in coordinates),
    )
    return attempt, root if root is not None and root.passes_numerical_gate else None


def _cluster_roots(
    roots: Sequence[GenotypeRoot],
    *,
    expression: float,
    transporter: str,
    span: NDArray[np.float64],
    relative_tolerance: float,
) -> tuple[GenotypeRoot, ...]:
    if not math.isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise ValueError("root clustering tolerance must be finite and positive")
    scale = np.asarray(span, dtype=float)
    if scale.ndim != 1 or np.any(~np.isfinite(scale)) or np.any(scale <= 0.0):
        raise ValueError("root clustering span must be a finite positive vector")
    ordered = tuple(
        sorted(
            roots,
            key=lambda item: (
                tuple(float(value) for value in item.coordinates),
                float(item.max_abs_scaled_independent_rhs),
                tuple(item.member_start_ids),
            ),
        )
    )
    points = tuple(np.asarray(root.coordinates, dtype=float) / scale for root in ordered)
    if any(point.shape != scale.shape or np.any(~np.isfinite(point)) for point in points):
        raise ValueError("root coordinates must match the finite clustering span")
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

    for left in range(len(points)):
        for right in range(left + 1, len(points)):
            if float(np.max(np.abs(points[left] - points[right]))) <= relative_tolerance:
                union(left, right)
    grouped: dict[int, list[GenotypeRoot]] = {}
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
    retained: list[GenotypeRoot] = []
    for index, cluster in enumerate(clusters):
        best = min(
            cluster,
            key=lambda item: (
                item.max_abs_scaled_independent_rhs,
                tuple(item.coordinates),
                tuple(item.member_start_ids),
            ),
        )
        retained.append(
            GenotypeRoot(
                **{
                    **asdict(best),
                    "root_id": (
                        f"{transporter}_E{expression:.8f}:B{index:02d}"
                    ),
                    "member_start_ids": tuple(
                        sorted(
                            start_id
                            for item in cluster
                            for start_id in item.member_start_ids
                        )
                    ),
                }
            )
        )
    return tuple(retained)


def continue_genotype_expression(
    model: Any,
    wt_resting_state: ArrayLike,
    *,
    transporter: str,
    expression_grid: Sequence[float] = DEFAULT_EXPRESSION_GRID,
    spec: WTCalibrationSpec | None = None,
    local_start_count: int = 3,
    local_start_fraction: float = 1.0e-3,
    maximum_connected_step_fraction: float = 0.25,
    max_nfev: int = 3000,
) -> GenotypeContinuation:
    """Continue one frozen WT root to exact transporter deletion.

    No parameter is optimized.  The ten fitted variables at every step are
    exactly ``REST_COORDINATE_NAMES``; OTHER is read once from the frozen
    parameter object and reported unchanged.  Every numerically admissible
    local branch is retained.  Selection uses only normalized distance from
    the preceding selected root.
    """

    name = _checked_transporter(transporter)
    grid = _expression_grid(expression_grid)
    declared = spec or WTCalibrationSpec()
    if local_start_count not in (1, 3):
        raise ValueError("local_start_count must be one or three")
    if not 0.0 < local_start_fraction < 0.05:
        raise ValueError("local_start_fraction must lie in (0, 0.05)")
    if not 0.0 < maximum_connected_step_fraction <= 1.0:
        raise ValueError("maximum connected step fraction must lie in (0, 1]")
    if max_nfev < 1:
        raise ValueError("max_nfev must be positive")

    lower = declared.coordinate_lower.copy()
    upper = declared.coordinate_upper.copy()
    span = upper - lower
    previous = _core_coordinates(model, wt_resting_state)
    if np.any(previous <= lower) or np.any(previous >= upper):
        raise ValueError("frozen WT root must lie strictly inside continuation bounds")
    regulatory_suffix = _regulatory_suffix(model)
    fixed_other = float(model.parameters.geometry.cell_impermeant_osmoles_fmol)
    steps: list[ContinuationStep] = []
    failure: str | None = None

    for expression in grid:
        genotype = genotype_with_expression(name, expression)
        starts: list[tuple[str, NDArray[np.float64]]] = [
            ("previous_connected_root", previous.copy())
        ]
        if local_start_count == 3:
            interior = 1.0e-10 * span
            starts.extend(
                (
                    (
                        "local_minus",
                        np.clip(
                            previous - local_start_fraction * span,
                            lower + interior,
                            upper - interior,
                        ),
                    ),
                    (
                        "local_plus",
                        np.clip(
                            previous + local_start_fraction * span,
                            lower + interior,
                            upper - interior,
                        ),
                    ),
                )
            )
        attempts: list[ContinuationAttempt] = []
        admissible_roots: list[GenotypeRoot] = []
        for start_id, start in starts:
            attempt, root = _attempt_and_root(
                model=model,
                transporter=name,
                expression=expression,
                genotype=genotype,
                regulatory_suffix=regulatory_suffix,
                spec=declared,
                start_id=start_id,
                start=start,
                lower=lower,
                upper=upper,
                max_nfev=max_nfev,
            )
            attempts.append(attempt)
            if root is not None:
                admissible_roots.append(root)
        roots = _cluster_roots(
            admissible_roots,
            expression=expression,
            transporter=name,
            span=span,
            relative_tolerance=declared.cluster_relative_tolerance,
        )
        if not roots:
            failure = f"no admissible root at {name} expression {expression:.8f}"
            steps.append(
                ContinuationStep(
                    expression=expression,
                    attempts=tuple(attempts),
                    roots=(),
                    selected_root_id=None,
                    selected_distance_from_previous_normalized=None,
                    connected_to_previous=False,
                    alternate_root_count=0,
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
            key=lambda item: (
                distances[item.root_id], item.max_abs_scaled_independent_rhs
            ),
        )
        distance = distances[selected.root_id]
        connected = bool(distance <= maximum_connected_step_fraction)
        status = "CONNECTED"
        if not connected:
            status = "ROOT_FOUND_BUT_BRANCH_JUMP_FLAGGED"
            failure = (
                f"nearest root jump {distance:.6g} exceeds declared normalized limit "
                f"at {name} expression {expression:.8f}"
            )
        elif len(roots) > 1:
            status = "CONNECTED_WITH_ALTERNATE_ROOTS_RETAINED"
        steps.append(
            ContinuationStep(
                expression=expression,
                attempts=tuple(attempts),
                roots=roots,
                selected_root_id=selected.root_id,
                selected_distance_from_previous_normalized=distance,
                connected_to_previous=connected,
                alternate_root_count=max(0, len(roots) - 1),
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
    return GenotypeContinuation(
        transporter=name,
        expression_grid=grid,
        fixed_other_impermeant_osmoles_fmol=fixed_other,
        independent_equation_count=len(CORE_INDEPENDENT_ROWS),
        fitted_coordinate_count=len(REST_COORDINATE_NAMES),
        fitted_parameter_count=0,
        steps=tuple(steps),
        completed_to_exact_zero=completed,
        branch_failure=None if completed else failure or "continuation incomplete",
    )


def continue_ae4_expression(
    model: Any, wt_resting_state: ArrayLike, **kwargs: Any
) -> GenotypeContinuation:
    """Convenience wrapper for WT-to-AE4-null resting continuation."""

    return continue_genotype_expression(
        model, wt_resting_state, transporter="AE4", **kwargs
    )


def continue_ae2_expression(
    model: Any, wt_resting_state: ArrayLike, **kwargs: Any
) -> GenotypeContinuation:
    """Convenience wrapper for independent WT-to-AE2-null validation."""

    return continue_genotype_expression(
        model, wt_resting_state, transporter="AE2", **kwargs
    )


def _deleted_transport_magnitude(evaluation: Any, transporter: str) -> float:
    if transporter == "AE4":
        contribution = evaluation.diagnostics.ae4
        return float(
            max(
                abs(contribution.na_cell_fmol_s),
                abs(contribution.k_cell_fmol_s),
                abs(contribution.cl_cell_fmol_s),
                abs(contribution.hco3_cell_fmol_s),
                abs(contribution.transported_charge_fmol_s),
            )
        )
    return float(abs(evaluation.diagnostics.homeostasis.ae2_inward_fmol_s))


def simulate_genotype(
    model: Any,
    resting_state: ArrayLike,
    *,
    transporter: str,
    genotype: Genotype,
    solver: SolverSpecification,
    time_s: ArrayLike | None = None,
) -> GenotypeTrajectory:
    """Integrate one frozen WT or deletion trajectory without target scoring."""

    name = _checked_transporter(transporter)
    grid = (
        physical_time_grid(duration_s=600.0, sample_step_s=1.0)
        if time_s is None
        else np.asarray(time_s, dtype=float)
    )
    if (
        grid.ndim != 1
        or grid.size < 3
        or grid[0] != 0.0
        or np.any(~np.isfinite(grid))
        or np.any(np.diff(grid) <= 0.0)
    ):
        raise ValueError("time grid must start at zero and increase through an endpoint")
    y0 = _root_state_with_basal_regulation(model, resting_state)
    integration_grid = grid[1:]
    result = model.solve_dynamics(
        (float(integration_grid[0]), float(integration_grid[-1])),
        initial_state=y0,
        genotype=genotype,
        method=solver.method,
        rtol=solver.rtol,
        atol=solver.atol_vector(model.state_names),
        t_eval=integration_grid,
        max_step_s=solver.max_step_s,
    )
    complete = bool(result.y.shape[1] == integration_grid.size)
    if complete:
        states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    else:
        grid = np.r_[0.0, np.asarray(result.t, dtype=float)]
        states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    sample_count = grid.size
    flow = np.full(sample_count, np.nan)
    cell_na = np.full(sample_count, np.nan)
    cell_k = np.full(sample_count, np.nan)
    cell_cl = np.full(sample_count, np.nan)
    cell_ph = np.full(sample_count, np.nan)
    cell_volume = np.full(sample_count, np.nan)
    lumen_na = np.full(sample_count, np.nan)
    lumen_k = np.full(sample_count, np.nan)
    lumen_cl = np.full(sample_count, np.nan)
    lumen_ph = np.full(sample_count, np.nan)
    lumen_volume = np.full(sample_count, np.nan)
    multiplier = np.full(sample_count, np.nan)
    conservation_max = {
        key: 0.0 for key in CONSERVATION_RESIDUAL_TOLERANCES
    }
    deleted_max = 0.0
    for index, (time_value, state) in enumerate(zip(grid, states.T)):
        evaluation = model.evaluate(float(time_value), state, genotype=genotype)
        diagnostics = evaluation.diagnostics
        decoded = model.layout.decode(state)
        ci = diagnostics.observables.cell_concentrations_mM
        li = diagnostics.observables.lumen_concentrations_mM
        flow[index] = diagnostics.water.lumen_outflow_pL_s
        cell_na[index] = ci["na"]
        cell_k[index] = ci["k"]
        cell_cl[index] = ci["cl"]
        cell_ph[index] = diagnostics.observables.cell_acid_base.ph
        cell_volume[index] = decoded.cell.volume_pL
        lumen_na[index] = li["na"]
        lumen_k[index] = li["k"]
        lumen_cl[index] = li["cl"]
        lumen_ph[index] = diagnostics.observables.lumen_acid_base.ph
        lumen_volume[index] = decoded.lumen.volume_pL
        multiplier[index] = float(
            diagnostics.regulatory.get(
                "capacity_multiplier", diagnostics.regulatory.get("activity", 1.0)
            )
        )
        if (
            (name == "AE4" and genotype.ae4_expression == 0.0)
            or (name == "AE2" and genotype.ae2_expression == 0.0)
        ):
            deleted_max = max(
                deleted_max, _deleted_transport_magnitude(evaluation, name)
            )
        unknown = set(diagnostics.conservation_residuals) - set(conservation_max)
        if unknown:
            raise KeyError(
                f"conservation diagnostics lack declared tolerances: {sorted(unknown)}"
            )
        for key, value in diagnostics.conservation_residuals.items():
            conservation_max[key] = max(conservation_max[key], abs(float(value)))
    cumulative = np.r_[0.0, cumulative_trapezoid(flow, grid)]
    conservation_ratio = max(
        conservation_max[key] / tolerance
        for key, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
    )
    positive = bool(np.all(states[: len(CORE_STATE_NAMES)] > 0.0))
    nonnegative_flow = bool(np.all(np.isfinite(flow)) and np.all(flow >= 0.0))
    exact_deletion = bool(deleted_max == 0.0)
    deletion_required = bool(
        (name == "AE4" and genotype.ae4_expression == 0.0)
        or (name == "AE2" and genotype.ae2_expression == 0.0)
    )
    numerical_gate = bool(
        result.success
        and complete
        and positive
        and nonnegative_flow
        and conservation_ratio <= 1.0
        and (not deletion_required or exact_deletion)
    )
    return GenotypeTrajectory(
        genotype_name=genotype.name,
        transporter=name,
        solver=solver,
        success=bool(result.success and complete),
        message=str(result.message),
        time_s=grid,
        state_names=tuple(model.state_names),
        states=states,
        flow_pL_s=flow,
        cumulative_flow_pL=cumulative,
        cell_na_mM=cell_na,
        cell_k_mM=cell_k,
        cell_cl_mM=cell_cl,
        cell_ph=cell_ph,
        cell_volume_pL=cell_volume,
        lumen_na_mM=lumen_na,
        lumen_k_mM=lumen_k,
        lumen_cl_mM=lumen_cl,
        lumen_ph=lumen_ph,
        lumen_volume_pL=lumen_volume,
        ae4_capacity_multiplier=multiplier,
        max_abs_deleted_transport_flux_fmol_s=deleted_max,
        max_abs_conservation_residuals=conservation_max,
        max_dimensionless_conservation_ratio=float(conservation_ratio),
        positive_core=positive,
        all_flow_nonnegative=nonnegative_flow,
        numerical_gate_pass=numerical_gate,
    )


def _positive_ratio(numerator: float, denominator: float, label: str) -> float:
    values = (float(numerator), float(denominator))
    if not all(math.isfinite(value) for value in values) or denominator <= 0.0:
        raise ValueError(f"{label} ratio requires a finite positive WT denominator")
    return float(numerator / denominator)


def summarize_genotype_pair(
    wt: GenotypeTrajectory,
    genotype: GenotypeTrajectory,
    *,
    landmarks_s: Sequence[float] = (60.0, 120.0, 180.0, 600.0),
) -> GenotypePairSummary:
    """Form generic ratios against the frozen WT trajectory.

    Landmark times are reporting coordinates only.  This function declares no
    expected direction, interval, or pass/fail phenotype criterion.
    """

    if wt.transporter != genotype.transporter:
        raise ValueError("paired trajectories must use the same transporter label")
    if wt.solver != genotype.solver:
        raise ValueError("paired trajectories must use the same solver specification")
    if wt.time_s.shape != genotype.time_s.shape or not np.array_equal(
        wt.time_s, genotype.time_s
    ):
        raise ValueError("paired trajectories must share the exact time grid")
    landmarks = tuple(float(value) for value in landmarks_s)
    if (
        not landmarks
        or any(not math.isfinite(value) or value <= 0.0 for value in landmarks)
        or any(right <= left for left, right in zip(landmarks, landmarks[1:]))
        or landmarks[-1] > wt.time_s[-1]
    ):
        raise ValueError("landmarks must be positive, increasing, and inside the trajectory")
    cumulative_ratios: dict[str, float] = {}
    flow_ratios: dict[str, float] = {}
    for landmark in landmarks:
        key = f"{landmark:g}_s"
        cumulative_ratios[key] = _positive_ratio(
            float(np.interp(landmark, genotype.time_s, genotype.cumulative_flow_pL)),
            float(np.interp(landmark, wt.time_s, wt.cumulative_flow_pL)),
            f"cumulative at {key}",
        )
        flow_ratios[key] = _positive_ratio(
            float(np.interp(landmark, genotype.time_s, genotype.flow_pL_s)),
            float(np.interp(landmark, wt.time_s, wt.flow_pL_s)),
            f"flow at {key}",
        )
    wt_endpoint = {
        "flow_pL_s": float(wt.flow_pL_s[-1]),
        "cumulative_flow_pL": float(wt.cumulative_flow_pL[-1]),
        "cell_na_mM": float(wt.cell_na_mM[-1]),
        "cell_k_mM": float(wt.cell_k_mM[-1]),
        "cell_cl_mM": float(wt.cell_cl_mM[-1]),
        "cell_ph": float(wt.cell_ph[-1]),
        "cell_volume_pL": float(wt.cell_volume_pL[-1]),
        "lumen_na_mM": float(wt.lumen_na_mM[-1]),
        "lumen_k_mM": float(wt.lumen_k_mM[-1]),
        "lumen_cl_mM": float(wt.lumen_cl_mM[-1]),
        "lumen_ph": float(wt.lumen_ph[-1]),
        "lumen_volume_pL": float(wt.lumen_volume_pL[-1]),
    }
    genotype_endpoint = {
        "flow_pL_s": float(genotype.flow_pL_s[-1]),
        "cumulative_flow_pL": float(genotype.cumulative_flow_pL[-1]),
        "cell_na_mM": float(genotype.cell_na_mM[-1]),
        "cell_k_mM": float(genotype.cell_k_mM[-1]),
        "cell_cl_mM": float(genotype.cell_cl_mM[-1]),
        "cell_ph": float(genotype.cell_ph[-1]),
        "cell_volume_pL": float(genotype.cell_volume_pL[-1]),
        "lumen_na_mM": float(genotype.lumen_na_mM[-1]),
        "lumen_k_mM": float(genotype.lumen_k_mM[-1]),
        "lumen_cl_mM": float(genotype.lumen_cl_mM[-1]),
        "lumen_ph": float(genotype.lumen_ph[-1]),
        "lumen_volume_pL": float(genotype.lumen_volume_pL[-1]),
    }
    endpoint_ratios = {
        key: _positive_ratio(genotype_endpoint[key], wt_endpoint[key], key)
        for key in wt_endpoint
        if key != "cell_ph" and key != "lumen_ph"
    }
    # pH is logarithmic, so report differences rather than invalid ratios.
    endpoint_ratios["cell_ph_difference"] = (
        genotype_endpoint["cell_ph"] - wt_endpoint["cell_ph"]
    )
    endpoint_ratios["lumen_ph_difference"] = (
        genotype_endpoint["lumen_ph"] - wt_endpoint["lumen_ph"]
    )
    return GenotypePairSummary(
        transporter=wt.transporter,
        solver_label=wt.solver.label,
        wt_success=wt.success,
        genotype_success=genotype.success,
        paired_numerical_gate_pass=bool(
            wt.numerical_gate_pass and genotype.numerical_gate_pass
        ),
        cumulative_ratio_at_landmarks=cumulative_ratios,
        flow_ratio_at_landmarks=flow_ratios,
        wt_endpoint=wt_endpoint,
        genotype_endpoint=genotype_endpoint,
        genotype_to_wt_endpoint_ratios=endpoint_ratios,
        wt_max_dimensionless_conservation_ratio=(
            wt.max_dimensionless_conservation_ratio
        ),
        genotype_max_dimensionless_conservation_ratio=(
            genotype.max_dimensionless_conservation_ratio
        ),
    )


def _relative_difference(left: ArrayLike, right: ArrayLike, floor: float) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("cross-solver arrays must be finite and aligned")
    scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), float(floor))
    return float(np.max(np.abs(a - b) / scale))


def evaluate_dual_solver_genotype_pair(
    model: Any,
    wt_resting_state: ArrayLike,
    genotype_resting_state: ArrayLike,
    *,
    transporter: str,
    time_s: ArrayLike | None = None,
    landmarks_s: Sequence[float] = (60.0, 120.0, 180.0, 600.0),
) -> DualSolverGenotypeEvaluation:
    """Run frozen WT and deletion roots with both production stiff solvers."""

    name = _checked_transporter(transporter)
    genotype = genotype_with_expression(name, 0.0)
    grid = (
        physical_time_grid(duration_s=600.0, sample_step_s=1.0)
        if time_s is None
        else np.asarray(time_s, dtype=float)
    )
    summaries: dict[str, GenotypePairSummary] = {}
    wt_trajectories: dict[str, GenotypeTrajectory] = {}
    genotype_trajectories: dict[str, GenotypeTrajectory] = {}
    for solver in (PRODUCTION_RADAU, PRODUCTION_BDF):
        wt_trajectory = simulate_genotype(
            model,
            wt_resting_state,
            transporter=name,
            genotype=WT,
            solver=solver,
            time_s=grid,
        )
        genotype_trajectory = simulate_genotype(
            model,
            genotype_resting_state,
            transporter=name,
            genotype=genotype,
            solver=solver,
            time_s=grid,
        )
        wt_trajectories[solver.label] = wt_trajectory
        genotype_trajectories[solver.label] = genotype_trajectory
        summaries[solver.label] = summarize_genotype_pair(
            wt_trajectory,
            genotype_trajectory,
            landmarks_s=landmarks_s,
        )
    wt_radau = wt_trajectories[PRODUCTION_RADAU.label]
    wt_bdf = wt_trajectories[PRODUCTION_BDF.label]
    genotype_radau = genotype_trajectories[PRODUCTION_RADAU.label]
    genotype_bdf = genotype_trajectories[PRODUCTION_BDF.label]
    cross = {
        "wt_state": _relative_difference(wt_radau.states, wt_bdf.states, 1.0e-10),
        "wt_flow": _relative_difference(
            wt_radau.flow_pL_s, wt_bdf.flow_pL_s, 1.0e-12
        ),
        "wt_cumulative": _relative_difference(
            wt_radau.cumulative_flow_pL,
            wt_bdf.cumulative_flow_pL,
            1.0e-12,
        ),
        "genotype_state": _relative_difference(
            genotype_radau.states, genotype_bdf.states, 1.0e-10
        ),
        "genotype_flow": _relative_difference(
            genotype_radau.flow_pL_s, genotype_bdf.flow_pL_s, 1.0e-12
        ),
        "genotype_cumulative": _relative_difference(
            genotype_radau.cumulative_flow_pL,
            genotype_bdf.cumulative_flow_pL,
            1.0e-12,
        ),
    }
    cross_solver_gate = bool(
        cross
        and all(
            math.isfinite(float(value))
            and float(value) <= GENOTYPE_SOLVER_RELATIVE_TOLERANCE
            for value in cross.values()
        )
    )
    return DualSolverGenotypeEvaluation(
        transporter=name,
        evaluations=summaries,
        wt_trajectories=wt_trajectories,
        genotype_trajectories=genotype_trajectories,
        cross_solver_relative_differences=cross,
        cross_solver_relative_tolerance=GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        cross_solver_gate_pass=cross_solver_gate,
        both_solver_pairs_pass=bool(
            cross_solver_gate
            and all(
                summary.paired_numerical_gate_pass
                for summary in summaries.values()
            )
        ),
    )


def continuation_to_dict(report: GenotypeContinuation) -> Mapping[str, Any]:
    """Convert a continuation report to a stable JSON-compatible mapping."""

    return asdict(report)


__all__ = (
    "DEFAULT_EXPRESSION_GRID",
    "FROZEN_EFFECTIVE_CALCIUM_UM",
    "FROZEN_NKCC1_FULLY_ACTIVATED_MULTIPLIER",
    "FROZEN_REGULATORY_MEMBER_ID",
    "GENOTYPE_SOLVER_RELATIVE_TOLERANCE",
    "ContinuationAttempt",
    "ContinuationStep",
    "DualSolverGenotypeEvaluation",
    "GenotypeContinuation",
    "GenotypePairSummary",
    "GenotypeRoot",
    "GenotypeTrajectory",
    "build_frozen_genotype_dynamic_model",
    "continuation_to_dict",
    "continue_ae2_expression",
    "continue_ae4_expression",
    "continue_genotype_expression",
    "evaluate_dual_solver_genotype_pair",
    "frozen_dynamic_contract_metadata",
    "frozen_regulatory_member",
    "genotype_with_expression",
    "simulate_genotype",
    "summarize_genotype_pair",
)
