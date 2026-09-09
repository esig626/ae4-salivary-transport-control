"""Independent numerical checks for the Task 13 reconstruction.

This module intentionally does not import the production cycle or solver
implementation.  It provides a second, small numerical route for the final
audit:

* CTMC stationary probabilities from an SVD nullspace;
* the same probabilities from generator cofactors (matrix-tree theorem);
* explicit edge and branch currents for ring and two-state multigraph cycles;
* affinity calculations directly from stoichiometric log activities; and
* positive-root solving in logarithmic coordinates with ``scipy.root``.

The generator convention is row-oriented: ``Q[i, j]`` is the rate from state
``i`` to state ``j`` and each row sums to zero.  Stationarity is therefore
``Q.T @ p == 0``.  Keeping this convention and implementation local is part of
the independent-reproduction firewall; no fitted or held-out target is read
here.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Callable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import root
from scipy.optimize import nnls


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class StationaryCheck:
    """Two independent stationary solutions and their numerical diagnostics."""

    svd: FloatArray
    matrix_tree: FloatArray
    max_probability_difference: float
    svd_residual_inf: float
    matrix_tree_residual_inf: float


@dataclass(frozen=True)
class RootCheck:
    """Result of a positive root solved in log coordinates."""

    state: FloatArray
    log_state: FloatArray
    residual: FloatArray
    residual_inf: float
    success: bool
    message: str
    function_evaluations: int


@dataclass(frozen=True)
class SharedBranchCheck:
    """Stationary audit of the reduced SR2 two-conformation multigraph."""

    occupancy_o: float
    occupancy_i: float
    j_cl: float
    j_na: float
    j_k: float
    affinity_na: float
    affinity_k: float
    carrier_residual: float
    entropy_production: float


@dataclass(frozen=True)
class FixedNullDiagnostics:
    """Derived observables for the independently coded fixed-chassis null."""

    pH_i: float
    h_i_mM: float
    cell_volume_pL: float
    q_total: float
    j_nak: float
    j_nkcc1: float
    j_ae2: float
    j_nhe1: float
    j_buffer: float
    j_cl_signed: float
    j_k: float
    current_residual_apical: float
    current_residual_basolateral: float
    amount_rate_na: float
    amount_rate_k: float
    amount_rate_cl: float
    amount_rate_hco3: float
    lumen_na_rhs: float
    lumen_k_rhs: float
    volume_rate: float


@dataclass(frozen=True)
class SplitNullDiagnostics:
    """Independent diagnostics for the Round-5 split-cation chassis."""

    pH_i: float
    h_i_mM: float
    cell_volume_pL: float
    q_total: float
    j_nak_apical: float
    j_nak_basolateral: float
    j_k_apical: float
    j_k_basolateral: float
    j_t_na: float
    j_t_k: float
    j_cl_signed: float
    current_residual_apical: float
    current_residual_basolateral: float
    amount_rate_na: float
    amount_rate_k: float
    amount_rate_cl: float
    amount_rate_hco3: float
    ae4_na_source: float
    ae4_k_source: float
    ae4_cl_source: float
    ae4_hco3_source: float


@dataclass(frozen=True)
class FrozenSR2Sources:
    """Independent SR2 1:2:3 source calculation at one whole-cell state."""

    na_i: float
    k_i: float
    cl_i: float
    hco3_i: float
    na_branch_current: float
    k_branch_current: float
    stationary_residual_inf: float


@dataclass(frozen=True)
class ProjectionCheck:
    """Independent linear and sign-constrained source projection."""

    coefficients_unconstrained: FloatArray
    residual_unconstrained: FloatArray
    residual_norm_unconstrained: float
    coefficients_nonnegative: FloatArray
    residual_nonnegative: FloatArray
    residual_norm_nonnegative: float
    numerical_rank: int


def _validated_generator(generator: ArrayLike) -> FloatArray:
    q = np.asarray(generator, dtype=float)
    if q.ndim != 2 or q.shape[0] != q.shape[1] or q.shape[0] < 2:
        raise ValueError("generator must be a square matrix of dimension >= 2")
    if not np.all(np.isfinite(q)):
        raise ValueError("generator contains a non-finite value")
    off_diagonal = q.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if np.min(off_diagonal) < -1e-14:
        raise ValueError("off-diagonal generator entries must be nonnegative")
    scale = max(1.0, float(np.max(np.abs(q))))
    if float(np.max(np.abs(np.sum(q, axis=1)))) > 1e-11 * scale:
        raise ValueError("generator rows must sum to zero")
    if np.max(np.diag(q)) > 1e-14:
        raise ValueError("generator diagonal entries must be nonpositive")
    return q


def generator_from_edges(
    n_states: int,
    edges: Sequence[tuple[int, int, float]],
) -> FloatArray:
    """Build a row-oriented CTMC generator from directed rates.

    Parallel directed edges are allowed and are summed.  ``edges`` entries are
    ``(source_index, target_index, rate)``.
    """

    if n_states < 2:
        raise ValueError("n_states must be at least two")
    q = np.zeros((n_states, n_states), dtype=float)
    for source, target, rate in edges:
        if source == target:
            raise ValueError("self transitions are not explicit generator edges")
        if not 0 <= source < n_states or not 0 <= target < n_states:
            raise ValueError("edge index lies outside the state domain")
        if not math.isfinite(rate) or rate <= 0.0:
            raise ValueError("every transition rate must be finite and positive")
        q[source, target] += float(rate)
    q[np.diag_indices(n_states)] = -np.sum(q, axis=1)
    return _validated_generator(q)


def ring_generator(forward: ArrayLike, reverse: ArrayLike) -> FloatArray:
    """Build the generator of a reversible ring.

    ``forward[i]`` is the rate ``i -> (i+1) mod n`` and ``reverse[i]`` is the
    rate ``(i+1) mod n -> i`` on the same edge.
    """

    u = np.asarray(forward, dtype=float)
    w = np.asarray(reverse, dtype=float)
    if u.ndim != 1 or w.shape != u.shape or len(u) < 3:
        raise ValueError("forward and reverse must be equal-length vectors (n >= 3)")
    directed: list[tuple[int, int, float]] = []
    for i, (ui, wi) in enumerate(zip(u, w, strict=True)):
        j = (i + 1) % len(u)
        directed.append((i, j, float(ui)))
        directed.append((j, i, float(wi)))
    return generator_from_edges(len(u), directed)


def stationary_distribution_svd(generator: ArrayLike) -> FloatArray:
    """Solve ``Q.T p = 0`` from the smallest right singular vector."""

    q = _validated_generator(generator)
    _u, _s, vh = np.linalg.svd(q.T, full_matrices=True)
    p = np.asarray(vh[-1, :], dtype=float)
    if float(np.sum(p)) < 0.0:
        p = -p
    total = float(np.sum(p))
    if not math.isfinite(total) or total <= 0.0:
        raise RuntimeError("SVD null vector could not be normalized")
    p /= total
    tolerance = 5e-12 * max(1.0, float(np.max(np.abs(q))))
    if float(np.min(p)) < -tolerance:
        raise RuntimeError("stationary null vector has a material negative component")
    p[p < 0.0] = 0.0
    p /= float(np.sum(p))
    return p


def stationary_distribution_matrix_tree(generator: ArrayLike) -> FloatArray:
    """Solve stationarity from principal cofactors of ``-Q``.

    For a row-oriented irreducible generator, the cofactor obtained by deleting
    row and column ``i`` is proportional to the stationary probability of
    state ``i``.  This route shares no linear solve or normalization equation
    with the production implementation.
    """

    q = _validated_generator(generator)
    laplacian = -q
    weights = np.empty(q.shape[0], dtype=float)
    for i in range(q.shape[0]):
        minor = np.delete(np.delete(laplacian, i, axis=0), i, axis=1)
        sign, logabs = np.linalg.slogdet(minor)
        if sign <= 0.0 or not math.isfinite(float(logabs)):
            raise RuntimeError("generator is reducible or a matrix-tree cofactor is nonpositive")
        weights[i] = float(logabs)
    weights -= float(np.max(weights))
    probabilities = np.exp(weights)
    probabilities /= float(np.sum(probabilities))
    return probabilities


def stationary_crosscheck(generator: ArrayLike) -> StationaryCheck:
    """Return SVD and matrix-tree solutions with residual diagnostics."""

    q = _validated_generator(generator)
    p_svd = stationary_distribution_svd(q)
    p_tree = stationary_distribution_matrix_tree(q)
    return StationaryCheck(
        svd=p_svd,
        matrix_tree=p_tree,
        max_probability_difference=float(np.max(np.abs(p_svd - p_tree))),
        svd_residual_inf=float(np.linalg.norm(q.T @ p_svd, ord=np.inf)),
        matrix_tree_residual_inf=float(np.linalg.norm(q.T @ p_tree, ord=np.inf)),
    )


def ring_edge_currents(
    stationary: ArrayLike,
    forward: ArrayLike,
    reverse: ArrayLike,
) -> FloatArray:
    """Return oriented stationary currents on every reversible ring edge."""

    p = np.asarray(stationary, dtype=float)
    u = np.asarray(forward, dtype=float)
    w = np.asarray(reverse, dtype=float)
    if p.ndim != 1 or u.shape != p.shape or w.shape != p.shape:
        raise ValueError("stationary, forward, and reverse vectors must have equal shape")
    return np.asarray(
        [p[i] * u[i] - p[(i + 1) % len(p)] * w[i] for i in range(len(p))],
        dtype=float,
    )


def log_rate_affinity(forward: ArrayLike, reverse: ArrayLike) -> float:
    """Compute ``log(prod(forward)/prod(reverse))`` without overflow."""

    u = np.asarray(forward, dtype=float)
    w = np.asarray(reverse, dtype=float)
    if u.shape != w.shape or u.ndim != 1 or len(u) == 0:
        raise ValueError("forward and reverse must be nonempty equal-length vectors")
    if np.any(~np.isfinite(u)) or np.any(~np.isfinite(w)) or np.any(u <= 0.0) or np.any(w <= 0.0):
        raise ValueError("all rates must be finite and positive")
    return float(np.sum(np.log(u)) - np.sum(np.log(w)))


def chemical_affinity(
    activities: Mapping[str, float],
    reactant_coefficients: Mapping[str, float],
    product_coefficients: Mapping[str, float],
) -> float:
    """Return ``log(Q_reactants/Q_products)`` for a passive forward event."""

    answer = 0.0
    species = set(reactant_coefficients) | set(product_coefficients)
    for name in species:
        activity = float(activities[name])
        if not math.isfinite(activity) or activity <= 0.0:
            raise ValueError(f"activity {name!r} must be finite and positive")
        coefficient = float(reactant_coefficients.get(name, 0.0)) - float(
            product_coefficients.get(name, 0.0)
        )
        answer += coefficient * math.log(activity)
    return float(answer)


def exact_rational_rank(matrix: Sequence[Sequence[int | float | Fraction]]) -> int:
    """Compute exact row rank by Fraction Gaussian elimination.

    This is deliberately separate from all floating-point projection code and
    is intended for the integer module-signature matrices in Stage B.
    """

    rows = [[Fraction(value) for value in row] for row in matrix]
    if not rows:
        return 0
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows must have equal length")
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_entry = rows[pivot_row][column]
        rows[pivot_row] = [value / pivot_entry for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_component
                for value, pivot_component in zip(rows[row], rows[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def source_projection(target: ArrayLike, signatures: ArrayLike) -> ProjectionCheck:
    """Project a target source onto module columns by two solver routes.

    The unconstrained calculation uses an SVD-backed pseudoinverse through
    ``numpy.linalg.lstsq``.  The sign-restricted calculation uses the Lawson-
    Hanson nonnegative least-squares algorithm in ``scipy.optimize.nnls``.
    ``signatures`` has balance coordinates in rows and modules in columns.
    """

    vector = np.asarray(target, dtype=float)
    matrix = np.asarray(signatures, dtype=float)
    if vector.ndim != 1 or matrix.ndim != 2 or matrix.shape[0] != vector.size:
        raise ValueError("signature rows must match the one-dimensional target")
    if np.any(~np.isfinite(vector)) or np.any(~np.isfinite(matrix)):
        raise ValueError("projection inputs must be finite")
    coefficients, _sum_squares, rank, _singular_values = np.linalg.lstsq(
        matrix, vector, rcond=None
    )
    residual = vector - matrix @ coefficients
    coefficients_nnls, _reported_norm = nnls(matrix, vector)
    residual_nnls = vector - matrix @ coefficients_nnls
    return ProjectionCheck(
        coefficients_unconstrained=np.asarray(coefficients, dtype=float),
        residual_unconstrained=np.asarray(residual, dtype=float),
        residual_norm_unconstrained=float(np.linalg.norm(residual)),
        coefficients_nonnegative=np.asarray(coefficients_nnls, dtype=float),
        residual_nonnegative=np.asarray(residual_nnls, dtype=float),
        residual_norm_nonnegative=float(np.linalg.norm(residual_nnls)),
        numerical_rank=int(rank),
    )


def shared_two_state_check(
    *,
    a: float,
    b: float,
    c_na: float,
    d_na: float,
    c_k: float,
    d_k: float,
) -> SharedBranchCheck:
    """Audit the exact two-state SR2 multigraph independently.

    ``a`` is the Cl-loading ``O -> I`` rate and ``b`` its reverse.  Each
    ``c_s`` is the cation-loaded return rate ``I -> O`` and ``d_s`` its
    reverse.  Parallel channels are retained when calculating branch currents,
    although the stationary generator contains only their summed rates.
    """

    rates = (a, b, c_na, d_na, c_k, d_k)
    if any(not math.isfinite(x) or x <= 0.0 for x in rates):
        raise ValueError("all two-state multigraph rates must be finite and positive")
    q = generator_from_edges(
        2,
        (
            (0, 1, a),
            (1, 0, b),
            (1, 0, c_na),
            (0, 1, d_na),
            (1, 0, c_k),
            (0, 1, d_k),
        ),
    )
    p_o, p_i = stationary_distribution_svd(q)
    j_cl = float(a * p_o - b * p_i)
    j_na = float(c_na * p_i - d_na * p_o)
    j_k = float(c_k * p_i - d_k * p_o)
    affinity_na = float(math.log(a * c_na / (b * d_na)))
    affinity_k = float(math.log(a * c_k / (b * d_k)))
    entropy = float(j_na * affinity_na + j_k * affinity_k)
    return SharedBranchCheck(
        occupancy_o=float(p_o),
        occupancy_i=float(p_i),
        j_cl=j_cl,
        j_na=j_na,
        j_k=j_k,
        affinity_na=affinity_na,
        affinity_k=affinity_k,
        carrier_residual=float(j_cl - j_na - j_k),
        entropy_production=entropy,
    )


def solve_positive_root(
    residual: Callable[[FloatArray], ArrayLike],
    initial_state: ArrayLike,
    *,
    method: str = "hybr",
    tolerance: float = 1e-10,
) -> RootCheck:
    """Solve a positive-state residual with ``scipy.root`` in log coordinates.

    The callback receives the physical positive state.  Solving for ``log(x)``
    independently enforces positivity without the bounds/least-squares route
    used by the main analysis.
    """

    x0 = np.asarray(initial_state, dtype=float)
    if x0.ndim != 1 or len(x0) == 0 or np.any(~np.isfinite(x0)) or np.any(x0 <= 0.0):
        raise ValueError("initial_state must be a finite, positive vector")

    def log_residual(log_state: FloatArray) -> FloatArray:
        physical = np.exp(log_state)
        value = np.asarray(residual(physical), dtype=float)
        if value.shape != x0.shape:
            raise ValueError("residual shape must match initial_state shape")
        if np.any(~np.isfinite(value)):
            raise FloatingPointError("residual returned a non-finite value")
        return value

    solution = root(log_residual, np.log(x0), method=method, tol=tolerance)
    physical = np.exp(np.asarray(solution.x, dtype=float))
    final_residual = np.asarray(residual(physical), dtype=float)
    return RootCheck(
        state=physical,
        log_state=np.asarray(solution.x, dtype=float),
        residual=final_residual,
        residual_inf=float(np.linalg.norm(final_residual, ord=np.inf)),
        success=bool(solution.success),
        message=str(solution.message),
        function_evaluations=int(solution.nfev),
    )


# Frozen Task-12 comparative-chassis constants.  They are repeated here rather
# than imported from the production chassis so this audit cannot inherit its
# evaluator or root solver by accident.  Their provenance and equations are in
# analysis/12_ae4_mechanism_reconstruction/chassis.md.
_FIXED_NULL = {
    "buffer_scale": 0.199652558943136,
    "co2_transport_scale": 1970.0,
    "buffer_forward": 0.132,
    "buffer_reverse": 312.0,
    "impermeant_charge_amount": 2439.50353087137,
    "co2_elimination_denominator": 3939.9145458622197,
    "co2_external_lumen_sum": 13.2,
    "nak_capacity": 0.00138020583534274,
    "nak_rate": 1_305_000.0,
    "k_e": 5.3,
    "nak_alpha": 0.641,
    "thermal_voltage_mv": 26.7137302360965,
    "faraday": 96485.3365,
    "cacc_half_activation": 100.26,
    "kcc_half_activation": 100.3,
    "cacc_hill": 1.49,
    "kcc_hill": 1.7,
    "tight_na_conductance": 357.863312079523,
    "tight_k_conductance": 25.893178360928,
    "cacc_conductance": 20_468_119.8347071,
    "kcc_conductance": 87_495_178.2340436,
    "nkcc_capacity": 0.0063812,
    "nkcc_forward": 157.55,
    "nkcc_reverse": 2.0096e-5,
    "nkcc_denominator_forward": 1.0306,
    "nkcc_denominator_reverse": 1.3852e-6,
    "nhe1_capacity": 0.674737308094733,
    "ae2_capacity": 0.399084872489615,
    "water_apical": 0.000579346121973713,
    "water_basolateral": 5.45083716910752e-5,
    "water_paracellular": 7.31969562708724e-6,
    "lumen_impermeant_osmolyte": 48.8001357236744,
    "external_osmotic_sum": 288.1,
    "area_to_litre_per_height": 4.52961672473868e-14,
    "cl_e": 102.6,
    "na_e": 140.2,
    "hco3_e_transport": 21.0,
    "ae2_cl_half": 5.6,
    "ae2_hco3_half": 0.0001,
    "nhe1_na_half": 15.0,
    "nhe1_h_half": 0.00045,
    "h_e": 3.8904514499428e-5,
}


def _fixed_null_evaluation(
    state: ArrayLike,
    *,
    calcium: float,
    ae2_activity_scale: float,
) -> tuple[FloatArray, FixedNullDiagnostics]:
    """Evaluate the seven-state AE4-null chassis from documented equations.

    The membrane-potential closure is evaluated by the explicit inverse of its
    two-by-two matrix rather than the main evaluator's linear-system call.
    """

    y = np.asarray(state, dtype=float)
    if y.shape != (7,) or np.any(~np.isfinite(y)) or np.any(y <= 0.0):
        raise ValueError("fixed-null state must be a finite positive seven-vector")
    if not math.isfinite(calcium) or calcium <= 0.0:
        raise ValueError("calcium must be finite and positive")
    if not math.isfinite(ae2_activity_scale) or ae2_activity_scale < 0.0:
        raise ValueError("ae2_activity_scale must be finite and nonnegative")
    p = _FIXED_NULL
    na_l, k_l, height, na_i, k_i, cl_i, hco3_i = map(float, y)

    h_i = hco3_i + cl_i + p["impermeant_charge_amount"] / height - na_i - k_i
    if h_i <= 0.0:
        raise ValueError("eliminated proton coordinate is nonpositive")
    co2_i = (
        p["co2_transport_scale"] * p["co2_external_lumen_sum"]
        - p["buffer_scale"] * p["buffer_reverse"] * hco3_i * h_i
    ) / p["co2_elimination_denominator"]
    cl_l = na_l + k_l

    na_m = 1e-3 * na_i
    k_e_m = 1e-3 * p["k_e"]
    j_nak = p["nak_capacity"] * p["nak_rate"] * k_e_m**2 * na_m**3 / (
        k_e_m**2 + p["nak_alpha"] * na_m**3
    )

    gate_cl = 1.0 / (1.0 + (p["cacc_half_activation"] / calcium) ** p["cacc_hill"])
    gate_k = 1.0 / (1.0 + (p["kcc_half_activation"] / calcium) ** p["kcc_hill"])
    g_cl = p["cacc_conductance"] * gate_cl
    g_k = p["kcc_conductance"] * gate_k
    v_cl = p["thermal_voltage_mv"] * math.log(cl_l / cl_i)
    v_k = p["thermal_voltage_mv"] * math.log(p["k_e"] / k_i)
    v_t_na = p["thermal_voltage_mv"] * math.log(na_l / p["na_e"])
    v_t_k = p["thermal_voltage_mv"] * math.log(k_l / p["k_e"])
    g_t = p["tight_na_conductance"] + p["tight_k_conductance"]
    m00 = g_cl + g_t
    m11 = g_k + g_t
    determinant = m00 * m11 - g_t**2
    rhs0 = (
        -g_cl * v_cl
        + p["tight_na_conductance"] * v_t_na
        + p["tight_k_conductance"] * v_t_k
    )
    rhs1 = (
        -p["faraday"] * j_nak
        + g_k * v_k
        - p["tight_na_conductance"] * v_t_na
        - p["tight_k_conductance"] * v_t_k
    )
    v_a = (m11 * rhs0 + g_t * rhs1) / determinant
    v_b = (g_t * rhs0 + m00 * rhs1) / determinant
    v_t = v_a - v_b
    j_cl_signed = g_cl * (v_a + v_cl) / p["faraday"]
    j_k = g_k * (v_b - v_k) / p["faraday"]
    j_t_na = p["tight_na_conductance"] * (v_t - v_t_na) / p["faraday"]
    j_t_k = p["tight_k_conductance"] * (v_t - v_t_k) / p["faraday"]

    q_a = p["water_apical"] * (
        2.0 * (na_l + k_l - na_i - k_i - h_i)
        - co2_i
        + p["lumen_impermeant_osmolyte"]
    )
    q_b = p["water_basolateral"] * (
        2.0 * (na_i + k_i + h_i) + co2_i - p["external_osmotic_sum"]
    )
    q_t = p["water_paracellular"] * (
        2.0 * (na_l + k_l)
        + p["lumen_impermeant_osmolyte"]
        - p["external_osmotic_sum"]
    )
    q_total = q_a + q_t
    volume_rate = q_b - q_a

    fourth_order = na_i * k_i * cl_i**2
    j_nkcc = p["nkcc_capacity"] * (
        (p["nkcc_forward"] - p["nkcc_reverse"] * fourth_order)
        / (p["nkcc_denominator_forward"] + p["nkcc_denominator_reverse"] * fourth_order)
    )
    j_ae2 = ae2_activity_scale * p["ae2_capacity"] * (
        (p["cl_e"] / (p["cl_e"] + p["ae2_cl_half"]))
        * (hco3_i / (hco3_i + p["ae2_hco3_half"]))
        - (cl_i / (cl_i + p["ae2_cl_half"]))
        * (p["hco3_e_transport"] / (p["hco3_e_transport"] + p["ae2_hco3_half"]))
    )
    j_nhe1 = p["nhe1_capacity"] * (
        (p["na_e"] / (p["na_e"] + p["nhe1_na_half"]))
        * (h_i / (p["nhe1_h_half"] + h_i))
        - (na_i / (na_i + p["nhe1_na_half"]))
        * (p["h_e"] / (p["nhe1_h_half"] + p["h_e"]))
    )
    j_buffer = p["buffer_scale"] * (
        p["buffer_forward"] * co2_i - p["buffer_reverse"] * hco3_i * h_i
    )

    amount_na = j_nkcc - 3.0 * j_nak + j_nhe1
    amount_k = j_nkcc + 2.0 * j_nak - j_k
    amount_cl = 2.0 * j_nkcc + j_ae2 + j_cl_signed
    amount_hco3 = j_buffer - j_ae2
    residual = np.asarray(
        (
            j_t_na - q_total * na_l,
            j_t_k - q_total * k_l,
            volume_rate,
            (amount_na - volume_rate * na_i) / height,
            (amount_k - volume_rate * k_i) / height,
            (amount_cl - volume_rate * cl_i) / height,
            (amount_hco3 - volume_rate * hco3_i) / height,
        ),
        dtype=float,
    )
    diagnostics = FixedNullDiagnostics(
        pH_i=math.log10(1000.0 / h_i),
        h_i_mM=h_i,
        cell_volume_pL=p["area_to_litre_per_height"] * height * 1e12,
        q_total=q_total,
        j_nak=j_nak,
        j_nkcc1=j_nkcc,
        j_ae2=j_ae2,
        j_nhe1=j_nhe1,
        j_buffer=j_buffer,
        j_cl_signed=j_cl_signed,
        j_k=j_k,
        current_residual_apical=j_cl_signed + j_t_na + j_t_k,
        current_residual_basolateral=j_nak + j_k - j_t_na - j_t_k,
        amount_rate_na=amount_na,
        amount_rate_k=amount_k,
        amount_rate_cl=amount_cl,
        amount_rate_hco3=amount_hco3,
        lumen_na_rhs=float(residual[0]),
        lumen_k_rhs=float(residual[1]),
        volume_rate=volume_rate,
    )
    return residual, diagnostics


def fixed_chassis_ae4_null_residual(
    state: ArrayLike,
    *,
    calcium: float = 0.05,
    ae2_activity_scale: float = 1.0,
) -> FloatArray:
    """Independent fixed-chassis AE4-null seven-equation residual."""

    residual, _diagnostics = _fixed_null_evaluation(
        state, calcium=calcium, ae2_activity_scale=ae2_activity_scale
    )
    return residual


def fixed_chassis_ae4_null_diagnostics(
    state: ArrayLike,
    *,
    calcium: float = 0.05,
    ae2_activity_scale: float = 1.0,
) -> FixedNullDiagnostics:
    """Return independent observables for an AE4-null fixed-chassis state."""

    _residual, diagnostics = _fixed_null_evaluation(
        state, calcium=calcium, ae2_activity_scale=ae2_activity_scale
    )
    return diagnostics


def frozen_sr2_shared_123_sources(
    state: ArrayLike,
    capacity: float,
) -> FrozenSR2Sources:
    """Evaluate the frozen SR2_SHARED_123/slow source independently.

    The eight-state graph, local-detailed-balance rate construction, and
    frozen pre-Stage-B assay lumps are restated here.  In particular, this
    function does not call the production cycle registry, adapter, stationary
    solver, or Round-5 optimizer.
    """

    y = np.asarray(state, dtype=float)
    if y.shape != (7,) or np.any(~np.isfinite(y)) or np.any(y <= 0.0):
        raise ValueError("state must be a finite positive seven-vector")
    if not math.isfinite(capacity) or capacity < 0.0:
        raise ValueError("capacity must be finite and nonnegative")
    _na_l, _k_l, _height, na_i, k_i, cl_i, hco3_i = map(float, y)

    names = (
        "Eo",
        "EoCl1",
        "EiCl1",
        "Ei",
        "EiNa2B3",
        "EoNa2B3",
        "EiK2B3",
        "EoK2B3",
    )
    index = {name: i for i, name in enumerate(names)}
    ec50_na = 58.06159101325025
    ec50_k = 74.17835667688657
    k_attempt = 1.1026527100857708
    loaded_bicarbonate_energy = 3.0 * math.log(25.0)
    energies = {
        "Eo": 0.0,
        "EoCl1": 0.0,
        "EiCl1": 0.0,
        "Ei": 0.0,
        "EiNa2B3": loaded_bicarbonate_energy + 2.0 * math.log(ec50_na),
        "EoNa2B3": loaded_bicarbonate_energy + 2.0 * math.log(ec50_na),
        "EiK2B3": loaded_bicarbonate_energy + 2.0 * math.log(ec50_k),
        "EoK2B3": loaded_bicarbonate_energy + 2.0 * math.log(ec50_k),
    }
    activities = {
        "na_i": na_i,
        "k_i": k_i,
        "cl_i": cl_i,
        "hco3_i": hco3_i,
        "na_e": _FIXED_NULL["na_e"],
        "k_e": _FIXED_NULL["k_e"],
        "cl_e": _FIXED_NULL["cl_e"],
        "hco3_e": _FIXED_NULL["hco3_e_transport"],
    }
    # (id, source, target, attempt multiplier, reservoir stoichiometry).
    edges = (
        ("common_cl_bind", "Eo", "EoCl1", 1.0, {"cl_e": 1.0}),
        ("common_cl_flip", "EoCl1", "EiCl1", 0.1, {}),
        ("common_cl_release", "EiCl1", "Ei", 1.0, {"cl_i": -1.0}),
        (
            "na_loaded_bind",
            "Ei",
            "EiNa2B3",
            1.0,
            {"na_i": 2.0, "hco3_i": 3.0},
        ),
        ("na_loaded_flip", "EiNa2B3", "EoNa2B3", 1.0, {}),
        (
            "na_loaded_release",
            "EoNa2B3",
            "Eo",
            1.0,
            {"na_e": -2.0, "hco3_e": -3.0},
        ),
        (
            "k_loaded_bind",
            "Ei",
            "EiK2B3",
            k_attempt,
            {"k_i": 2.0, "hco3_i": 3.0},
        ),
        ("k_loaded_flip", "EiK2B3", "EoK2B3", k_attempt, {}),
        (
            "k_loaded_release",
            "EoK2B3",
            "Eo",
            k_attempt,
            {"k_e": -2.0, "hco3_e": -3.0},
        ),
    )
    directed: list[tuple[int, int, float]] = []
    rate_rows: dict[str, tuple[int, int, float, float]] = {}
    for edge_id, source, target, attempt, stoichiometry in edges:
        delta = energies[target] - energies[source]
        log_forward = math.log(attempt) - 0.5 * delta
        log_reverse = math.log(attempt) + 0.5 * delta
        for species, coefficient in stoichiometry.items():
            if coefficient > 0.0:
                log_forward += coefficient * math.log(activities[species])
            elif coefficient < 0.0:
                log_reverse += -coefficient * math.log(activities[species])
        forward = math.exp(log_forward)
        reverse = math.exp(log_reverse)
        source_index = index[source]
        target_index = index[target]
        directed.extend(
            (
                (source_index, target_index, forward),
                (target_index, source_index, reverse),
            )
        )
        rate_rows[edge_id] = (source_index, target_index, forward, reverse)
    generator = generator_from_edges(len(names), directed)
    occupancy = stationary_distribution_svd(generator)

    def current(edge_id: str) -> float:
        source, target, forward, reverse = rate_rows[edge_id]
        return capacity * (
            occupancy[source] * forward - occupancy[target] * reverse
        )

    j_na = float(current("na_loaded_flip"))
    j_k = float(current("k_loaded_flip"))
    j_total = j_na + j_k
    return FrozenSR2Sources(
        na_i=-2.0 * j_na,
        k_i=-2.0 * j_k,
        cl_i=j_total,
        hco3_i=-3.0 * j_total,
        na_branch_current=j_na,
        k_branch_current=j_k,
        stationary_residual_inf=float(
            np.linalg.norm(generator.T @ occupancy, ord=np.inf)
        ),
    )


def _split_cation_evaluation(
    state: ArrayLike,
    *,
    apical_pump_fraction: float,
    apical_k_fraction: float,
    calcium: float,
    ae2_activity_scale: float,
    ae4_sources: Sequence[float] = (0.0, 0.0, 0.0, 0.0),
) -> tuple[FloatArray, SplitNullDiagnostics]:
    """Evaluate the Round-5 split equations by an explicit second route."""

    y = np.asarray(state, dtype=float)
    sources = np.asarray(ae4_sources, dtype=float)
    if y.shape != (7,) or np.any(~np.isfinite(y)) or np.any(y <= 0.0):
        raise ValueError("split-cation state must be a finite positive seven-vector")
    if sources.shape != (4,) or np.any(~np.isfinite(sources)):
        raise ValueError("AE4 sources must be a finite four-vector")
    if not math.isfinite(calcium) or calcium <= 0.0:
        raise ValueError("calcium must be finite and positive")
    if not math.isfinite(ae2_activity_scale) or ae2_activity_scale < 0.0:
        raise ValueError("ae2_activity_scale must be finite and nonnegative")
    fp = float(apical_pump_fraction)
    fk = float(apical_k_fraction)
    if not 0.0 <= fp <= 1.0 or not 0.0 <= fk <= 1.0:
        raise ValueError("split fractions must lie in [0, 1]")
    p = _FIXED_NULL
    na_l, k_l, height, na_i, k_i, cl_i, hco3_i = map(float, y)

    h_i = hco3_i + cl_i + p["impermeant_charge_amount"] / height - na_i - k_i
    if h_i <= 0.0:
        raise ValueError("eliminated proton coordinate is nonpositive")
    co2_i = (
        p["co2_transport_scale"] * p["co2_external_lumen_sum"]
        - p["buffer_scale"] * p["buffer_reverse"] * hco3_i * h_i
    ) / p["co2_elimination_denominator"]
    cl_l = na_l + k_l

    def pump_turnover(local_k_o: float) -> float:
        na_m = 1e-3 * na_i
        ko_m = 1e-3 * local_k_o
        return p["nak_rate"] * ko_m**2 * na_m**3 / (
            ko_m**2 + p["nak_alpha"] * na_m**3
        )

    j_pump_a = p["nak_capacity"] * fp * pump_turnover(k_l)
    j_pump_b = p["nak_capacity"] * (1.0 - fp) * pump_turnover(p["k_e"])
    gate_cl = 1.0 / (
        1.0 + (p["cacc_half_activation"] / calcium) ** p["cacc_hill"]
    )
    gate_k = 1.0 / (
        1.0 + (p["kcc_half_activation"] / calcium) ** p["kcc_hill"]
    )
    g_cl = p["cacc_conductance"] * gate_cl
    g_k_total = p["kcc_conductance"] * gate_k
    g_k_a = fk * g_k_total
    g_k_b = (1.0 - fk) * g_k_total
    v_cl = p["thermal_voltage_mv"] * math.log(cl_l / cl_i)
    v_k_a = p["thermal_voltage_mv"] * math.log(k_l / k_i)
    v_k_b = p["thermal_voltage_mv"] * math.log(p["k_e"] / k_i)
    v_t_na = p["thermal_voltage_mv"] * math.log(na_l / p["na_e"])
    v_t_k = p["thermal_voltage_mv"] * math.log(k_l / p["k_e"])
    g_t_na = p["tight_na_conductance"]
    g_t_k = p["tight_k_conductance"]
    g_t = g_t_na + g_t_k
    m00 = g_cl + g_k_a + g_t
    m11 = g_k_b + g_t
    determinant = m00 * m11 - g_t**2
    rhs0 = (
        -g_cl * v_cl
        + g_k_a * v_k_a
        + g_t_na * v_t_na
        + g_t_k * v_t_k
        - p["faraday"] * j_pump_a
    )
    rhs1 = (
        -p["faraday"] * j_pump_b
        + g_k_b * v_k_b
        - g_t_na * v_t_na
        - g_t_k * v_t_k
    )
    v_a = (m11 * rhs0 + g_t * rhs1) / determinant
    v_b = (g_t * rhs0 + m00 * rhs1) / determinant
    v_t = v_a - v_b
    j_cl_signed = g_cl * (v_a + v_cl) / p["faraday"]
    j_k_a = g_k_a * (v_a - v_k_a) / p["faraday"]
    j_k_b = g_k_b * (v_b - v_k_b) / p["faraday"]
    j_t_na = g_t_na * (v_t - v_t_na) / p["faraday"]
    j_t_k = g_t_k * (v_t - v_t_k) / p["faraday"]

    q_a = p["water_apical"] * (
        2.0 * (na_l + k_l - na_i - k_i - h_i)
        - co2_i
        + p["lumen_impermeant_osmolyte"]
    )
    q_b = p["water_basolateral"] * (
        2.0 * (na_i + k_i + h_i) + co2_i - p["external_osmotic_sum"]
    )
    q_t = p["water_paracellular"] * (
        2.0 * (na_l + k_l)
        + p["lumen_impermeant_osmolyte"]
        - p["external_osmotic_sum"]
    )
    q_total = q_a + q_t
    volume_rate = q_b - q_a
    fourth_order = na_i * k_i * cl_i**2
    j_nkcc = p["nkcc_capacity"] * (
        (p["nkcc_forward"] - p["nkcc_reverse"] * fourth_order)
        / (
            p["nkcc_denominator_forward"]
            + p["nkcc_denominator_reverse"] * fourth_order
        )
    )
    j_ae2 = ae2_activity_scale * p["ae2_capacity"] * (
        (p["cl_e"] / (p["cl_e"] + p["ae2_cl_half"]))
        * (hco3_i / (hco3_i + p["ae2_hco3_half"]))
        - (cl_i / (cl_i + p["ae2_cl_half"]))
        * (
            p["hco3_e_transport"]
            / (p["hco3_e_transport"] + p["ae2_hco3_half"])
        )
    )
    j_nhe1 = p["nhe1_capacity"] * (
        (p["na_e"] / (p["na_e"] + p["nhe1_na_half"]))
        * (h_i / (p["nhe1_h_half"] + h_i))
        - (na_i / (na_i + p["nhe1_na_half"]))
        * (p["h_e"] / (p["nhe1_h_half"] + p["h_e"]))
    )
    j_buffer = p["buffer_scale"] * (
        p["buffer_forward"] * co2_i
        - p["buffer_reverse"] * hco3_i * h_i
    )

    ae4_na, ae4_k, ae4_cl, ae4_hco3 = map(float, sources)
    j_pump_total = j_pump_a + j_pump_b
    amount_na = j_nkcc - 3.0 * j_pump_total + j_nhe1 + ae4_na
    amount_k = (
        j_nkcc + 2.0 * j_pump_total - j_k_a - j_k_b + ae4_k
    )
    amount_cl = 2.0 * j_nkcc + j_ae2 + j_cl_signed + ae4_cl
    amount_hco3 = j_buffer - j_ae2 + ae4_hco3
    lumen_na_source = j_t_na + 3.0 * j_pump_a
    lumen_k_source = j_t_k + j_k_a - 2.0 * j_pump_a
    residual = np.asarray(
        (
            lumen_na_source - q_total * na_l,
            lumen_k_source - q_total * k_l,
            volume_rate,
            (amount_na - volume_rate * na_i) / height,
            (amount_k - volume_rate * k_i) / height,
            (amount_cl - volume_rate * cl_i) / height,
            (amount_hco3 - volume_rate * hco3_i) / height,
        ),
        dtype=float,
    )
    diagnostics = SplitNullDiagnostics(
        pH_i=math.log10(1000.0 / h_i),
        h_i_mM=h_i,
        cell_volume_pL=p["area_to_litre_per_height"] * height * 1e12,
        q_total=q_total,
        j_nak_apical=j_pump_a,
        j_nak_basolateral=j_pump_b,
        j_k_apical=j_k_a,
        j_k_basolateral=j_k_b,
        j_t_na=j_t_na,
        j_t_k=j_t_k,
        j_cl_signed=j_cl_signed,
        current_residual_apical=(
            j_cl_signed + j_k_a + j_pump_a + j_t_na + j_t_k
        ),
        current_residual_basolateral=(
            j_pump_b + j_k_b - j_t_na - j_t_k
        ),
        amount_rate_na=amount_na,
        amount_rate_k=amount_k,
        amount_rate_cl=amount_cl,
        amount_rate_hco3=amount_hco3,
        ae4_na_source=ae4_na,
        ae4_k_source=ae4_k,
        ae4_cl_source=ae4_cl,
        ae4_hco3_source=ae4_hco3,
    )
    return residual, diagnostics


def split_cation_ae4_null_residual(
    state: ArrayLike,
    *,
    apical_pump_fraction: float,
    apical_k_fraction: float,
    calcium: float = 0.05,
    ae2_activity_scale: float = 1.0,
) -> FloatArray:
    """Independent seven-equation AE4-null residual for a split topology."""

    residual, _ = _split_cation_evaluation(
        state,
        apical_pump_fraction=apical_pump_fraction,
        apical_k_fraction=apical_k_fraction,
        calcium=calcium,
        ae2_activity_scale=ae2_activity_scale,
    )
    return residual


def split_cation_ae4_null_diagnostics(
    state: ArrayLike,
    *,
    apical_pump_fraction: float,
    apical_k_fraction: float,
    calcium: float = 0.05,
    ae2_activity_scale: float = 1.0,
) -> SplitNullDiagnostics:
    """Independent diagnostics for a split-topology AE4-null state."""

    _, diagnostics = _split_cation_evaluation(
        state,
        apical_pump_fraction=apical_pump_fraction,
        apical_k_fraction=apical_k_fraction,
        calcium=calcium,
        ae2_activity_scale=ae2_activity_scale,
    )
    return diagnostics


def split_cation_sr2_wt_residual(
    state: ArrayLike,
    *,
    capacity: float,
    apical_pump_fraction: float,
    apical_k_fraction: float,
    calcium: float = 0.05,
) -> FloatArray:
    """Independent Round-5 WT residual using frozen SR2_SHARED_123/slow."""

    source = frozen_sr2_shared_123_sources(state, capacity)
    residual, _ = _split_cation_evaluation(
        state,
        apical_pump_fraction=apical_pump_fraction,
        apical_k_fraction=apical_k_fraction,
        calcium=calcium,
        ae2_activity_scale=1.0,
        ae4_sources=(source.na_i, source.k_i, source.cl_i, source.hco3_i),
    )
    return residual


def split_cation_null_coordinate_residual(
    coordinates: ArrayLike,
    *,
    apical_pump_fraction: float,
    apical_k_fraction: float,
    calcium: float = 0.05,
) -> FloatArray:
    """Split-topology null residual on the positive proton coordinate."""

    state = fixed_null_coordinates_to_state(coordinates)
    return split_cation_ae4_null_residual(
        state,
        apical_pump_fraction=apical_pump_fraction,
        apical_k_fraction=apical_k_fraction,
        calcium=calcium,
    )


def fixed_null_state_to_coordinates(state: ArrayLike) -> FloatArray:
    """Replace ``HCO3_i`` by the positive eliminated-proton coordinate.

    This coordinate system preserves the historical seven equations but makes
    the tiny positive ``h_i`` domain explicit during an independent log-root
    solve.  It avoids relying on the production solver's bounds.
    """

    y = np.asarray(state, dtype=float)
    if y.shape != (7,) or np.any(~np.isfinite(y)) or np.any(y <= 0.0):
        raise ValueError("state must be a finite positive seven-vector")
    height, na_i, k_i, cl_i, hco3_i = map(float, (y[2], y[3], y[4], y[5], y[6]))
    h_i = hco3_i + cl_i + _FIXED_NULL["impermeant_charge_amount"] / height - na_i - k_i
    if h_i <= 0.0:
        raise ValueError("state implies a nonpositive eliminated-proton coordinate")
    return np.asarray([y[0], y[1], height, na_i, k_i, cl_i, h_i], dtype=float)


def fixed_null_coordinates_to_state(coordinates: ArrayLike) -> FloatArray:
    """Map ``(Na_l,K_l,H,Na_i,K_i,Cl_i,h_i)`` to the seven-state vector."""

    z = np.asarray(coordinates, dtype=float)
    if z.shape != (7,) or np.any(~np.isfinite(z)) or np.any(z <= 0.0):
        raise ValueError("coordinates must be a finite positive seven-vector")
    na_l, k_l, height, na_i, k_i, cl_i, h_i = map(float, z)
    hco3_i = h_i - cl_i - _FIXED_NULL["impermeant_charge_amount"] / height + na_i + k_i
    if hco3_i <= 0.0:
        raise ValueError("coordinates imply nonpositive intracellular bicarbonate")
    return np.asarray([na_l, k_l, height, na_i, k_i, cl_i, hco3_i], dtype=float)


def fixed_chassis_ae4_null_coordinate_residual(
    coordinates: ArrayLike,
    *,
    calcium: float = 0.05,
    ae2_activity_scale: float = 1.0,
) -> FloatArray:
    """Fixed-null residual on a positivity-preserving proton coordinate."""

    state = fixed_null_coordinates_to_state(coordinates)
    return fixed_chassis_ae4_null_residual(
        state, calcium=calcium, ae2_activity_scale=ae2_activity_scale
    )


def stage_b_target_state(
    free_coordinates: ArrayLike,
    *,
    cl_i_mM: float = 36.5,
    pH_i: float = 6.89,
) -> FloatArray:
    """Construct one member of the released Cl/pH target manifold.

    ``free_coordinates`` is ``(Na_l,K_l,height,Na_i,K_i)``.  No unmeasured
    knockout Na, K, lumen, or volume value is selected by this function.
    """

    free = np.asarray(free_coordinates, dtype=float)
    if free.shape != (5,) or np.any(~np.isfinite(free)) or np.any(free <= 0.0):
        raise ValueError("free_coordinates must be a finite positive five-vector")
    if cl_i_mM <= 0.0 or not 0.0 < pH_i < 14.0:
        raise ValueError("Cl and pH targets must be physical")
    na_l, k_l, height, na_i, k_i = map(float, free)
    h_i = 1000.0 * 10.0 ** (-pH_i)
    hco3_i = (
        h_i
        - cl_i_mM
        - _FIXED_NULL["impermeant_charge_amount"] / height
        + na_i
        + k_i
    )
    if hco3_i <= 0.0:
        raise ValueError("free coordinates imply nonpositive bicarbonate")
    return np.asarray(
        [na_l, k_l, height, na_i, k_i, cl_i_mM, hco3_i], dtype=float
    )


def fixed_null_required_balance(state: ArrayLike) -> FloatArray:
    """Return the 14-coordinate conditional source demand at one state.

    Coordinate order is ``Na_i,K_i,Cl_i,HCO3_i,H_i,CO2_i,Na_l,K_l,Cl_l,
    HCO3_l,q_apical,q_basolateral,q_tight,V_i``.  The eliminated ``CO2_i``
    and absent luminal bicarbonate demands are returned as NaN.  Charge rows
    are zero because the pre-existing voltage closure is re-solved.
    """

    raw, diagnostics = _fixed_null_evaluation(
        state, calcium=0.05, ae2_activity_scale=1.0
    )
    amount_h = (
        -diagnostics.amount_rate_na
        - diagnostics.amount_rate_k
        + diagnostics.amount_rate_cl
        + diagnostics.amount_rate_hco3
    )
    demand = np.zeros(14, dtype=float)
    demand[:5] = -np.asarray(
        [
            diagnostics.amount_rate_na,
            diagnostics.amount_rate_k,
            diagnostics.amount_rate_cl,
            diagnostics.amount_rate_hco3,
            amount_h,
        ]
    )
    demand[5] = math.nan
    demand[6] = -raw[0]
    demand[7] = -raw[1]
    demand[8] = -(raw[0] + raw[1])
    demand[9] = math.nan
    demand[10:13] = 0.0
    demand[13] = -raw[2]
    return demand


def central_difference_jacobian(
    function: Callable[[FloatArray], ArrayLike],
    point: ArrayLike,
    *,
    relative_step: float = 1e-5,
) -> FloatArray:
    """Coordinate-scaled central-difference Jacobian for audit calculations."""

    x = np.asarray(point, dtype=float)
    if x.ndim != 1 or np.any(~np.isfinite(x)):
        raise ValueError("point must be a finite vector")
    baseline = np.asarray(function(x), dtype=float)
    jacobian = np.empty((baseline.size, x.size), dtype=float)
    for i in range(x.size):
        step = relative_step * max(1.0, abs(float(x[i])))
        plus = x.copy()
        minus = x.copy()
        plus[i] += step
        minus[i] -= step
        jacobian[:, i] = (
            np.asarray(function(plus), dtype=float)
            - np.asarray(function(minus), dtype=float)
        ) / (2.0 * step)
    return jacobian


def stage_b_correction_jacobian(
    free_coordinates: ArrayLike,
    *,
    cl_i_mM: float = 36.5,
    pH_i: float = 6.89,
    relative_step: float = 1e-5,
) -> FloatArray:
    """Differentiate the 12 determined source rows across the target family."""

    free = np.asarray(free_coordinates, dtype=float)
    determined = np.asarray([0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13])

    def demand(point: FloatArray) -> FloatArray:
        state = stage_b_target_state(point, cl_i_mM=cl_i_mM, pH_i=pH_i)
        return fixed_null_required_balance(state)[determined]

    return central_difference_jacobian(demand, free, relative_step=relative_step)
