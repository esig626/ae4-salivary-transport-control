"""Common, knockout-blind calibration and comparison utilities.

The fixed-chassis evaluator is deliberately kept separate from the AE4
mechanism registry.  This module owns the comparison protocol: fit AE4-only
quantities to the wild-type resting balance, freeze them, and only then expose
the held-out knockout outputs.  No knockout observation is accepted by any
calibration function in this file.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.optimize import brentq, least_squares

from .chassis import (
    AE2_KNOCKOUT,
    AE4_KNOCKOUT,
    BASELINE_STATE,
    WT,
    AE4Environment,
    FixedChassis,
    SimulationResult,
)
from .mechanisms import (
    AE4Context,
    AE4Mechanism,
    AE4Parameters,
    AE4Contribution,
    CANDIDATE_SPECS,
    CandidateSpec,
    candidate_specs,
    evaluate_ae4,
)


BALANCE_NAMES = ("na_i", "k_i", "cl_i", "hco3_i")

# Historical implementation evidence: the internally closed seven-state WT
# row requires this AE4 contribution.  It is a WT closure constraint, not an
# independent assay and not a knockout calibration target.
HISTORICAL_WT_AE4_CYCLE_FLUX = 0.03888333867647999
HISTORICAL_WT_AE4_CONTRIBUTION = np.array(
    [
        -HISTORICAL_WT_AE4_CYCLE_FLUX,
        0.0,
        HISTORICAL_WT_AE4_CYCLE_FLUX,
        -2.0 * HISTORICAL_WT_AE4_CYCLE_FLUX,
    ],
    dtype=float,
)

PRIMARY_CALIBRATION_TARGET_NAMES = (
    "historical_WT_AE4_balance_contribution",
    "primary_Na_EC50_49_mM",
    "primary_Na_Hill_2.0",
    "primary_K_EC50_62_mM",
    "primary_K_Hill_1.8",
)

# Names only.  Values are joined after all calibration and simulation steps.
# Keeping them out of the calibration API makes leakage mechanically testable.
HELDOUT_TARGET_NAMES = (
    "primary_AE4_KO_stimulated_saliva_ratio",
    "primary_AE2_KO_stimulated_saliva_ratio",
)

# Independent resting measurements used to audit candidate-specific WT roots.
# The uncertainty values are reported SEMs from Peña-Münzenmayer et al. 2015
# Table 1.  They constrain resting-state acceptability, never the knockout-flow
# calibration or any AE4 kinetic parameter.
WT_REST_CL_MILLIMOLAR = 50.10
WT_REST_CL_SEM_MILLIMOLAR = 1.50
WT_REST_PH = 6.91
WT_REST_PH_SEM = 0.07
PUBLISHED_REST_NA_MILLIMOLAR = 20.0
PUBLISHED_REST_K_MILLIMOLAR = 120.0
PUBLISHED_REST_HCO3_MILLIMOLAR = 12.0
PUBLISHED_REST_VOLUME_PL = 1.3
REST_VOLUME_RELATIVE_TOLERANCE = 0.10

REST_RHS_SCALES = np.array(
    (1e-4, 1e-5, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3), dtype=float
)
REST_REDUCED_LOWER = np.array(
    (60.0, 0.5, 5.0, 50.0, 10.0, 1.0, math.log(1e-6)), dtype=float
)
REST_REDUCED_UPPER = np.array(
    (180.0, 30.0, 80.0, 180.0, 100.0, 100.0, math.log(1e-2)),
    dtype=float,
)

# The 1.25 recombinant forskolin result and approximate 1.6 total-WT and 3.0
# native AE4-only Fig. 1 responses are independent assay constraints.  F=1 is
# included only as the no-activation reference; F=10 is an explicit inverse-
# diagnostic stress test outside the source-supported envelope.
PKA_ASSAY_FOLDS = (1.25, 1.6, 3.0)
PKA_REFERENCE_FOLD = 1.0
PKA_INVERSE_DIAGNOSTIC_FOLD = 10.0
# Model-input conventions in the chassis-native (uncertified physical-time)
# coordinate.  They are not transcriptions of the 2--3 minute observation.
PKA_IMMEDIATE_MODEL_ONSET = 100.0
PKA_TIMING_SENSITIVITY_ONSETS = (120.0, 130.0)
REST_MULTISTART_SEED = 20260827
REST_MULTISTART_COUNT = 20


def comparison_time_grid() -> np.ndarray:
    """Output grid resolving the right side of the calcium discontinuity.

    Pre-step values are saved every chassis time unit.  The stimulated segment
    includes an explicit right limit and 0.1-unit spacing so its integral does
    not interpolate across the discontinuous calcium input.
    """

    pre = np.arange(0.0, 100.0 + 0.5, 1.0)
    post = np.arange(100.1, 200.0 + 0.05, 0.1)
    return np.concatenate((pre, np.array([100.0 + 1e-7]), post))


@dataclass(frozen=True)
class ScaleFit:
    """One-parameter non-negative least-squares fit to the WT balance vector."""

    activity_scale: float
    fitted_contribution: tuple[float, float, float, float]
    residual: tuple[float, float, float, float]
    l2_error: float
    max_abs_error: float


@dataclass(frozen=True)
class CandidateComparison:
    """Flat record used for both JSON and CSV comparison tables."""

    mechanism_id: str
    round: int
    family: str
    n_free_ae4_parameters: int
    n_fitted_ae4_parameters: int
    complexity_penalty: int
    calibration_target_names: str
    calibration_target_values: str
    independently_fixed_ae4_parameters: str
    new_modelling_assumptions: str
    mechanism_equivalence_alias: str
    heldout_target_names: str
    knockout_target_read_by_calibrator: bool
    activity_scale: float
    ae4_density_or_turnover_source: str
    non_ae4_chassis_fingerprint: str
    units_convention: str
    calibration_l2_error: float
    calibration_max_abs_error: float
    wt_baseline_closure_error: float
    wt_baseline_admissible: bool
    assay_admissible: bool
    thermodynamic_status: str
    electroneutrality_status: str
    simulation_status: str
    solver_method: str
    wt_endpoint_flow: float
    ae4_ko_endpoint_flow: float
    ae2_ko_endpoint_flow: float
    ae4_ko_endpoint_flow_ratio: float
    ae2_ko_endpoint_flow_ratio: float
    ae4_ko_integrated_flow_ratio: float
    ae2_ko_integrated_flow_ratio: float
    ae4_ko_cl_i_change_mM: float
    ae2_ko_cl_i_change_mM: float
    ae4_ko_pH_change: float
    ae2_ko_pH_change: float
    wt_endpoint_max_abs_raw_rhs: float
    ae4_ko_endpoint_max_abs_raw_rhs: float
    ae2_ko_endpoint_max_abs_raw_rhs: float
    endpoint_interpretation: str
    max_conservation_residual: float
    max_thermodynamic_violation: float
    failure_reason: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CalibratedCandidate:
    """Candidate, independently fixed assay values, and WT-only scale fit."""

    candidate_id: str
    parameters: AE4Parameters
    external_hco3_mM: float
    scale_fit: ScaleFit
    unit_contribution: tuple[float, float, float, float]
    baseline_diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class RestingEquilibrium:
    """Knockout-blind candidate-specific resting-root audit row."""

    mechanism_id: str
    scenario: str
    activity_scale: float
    calibration_mode: str
    solver_success: bool
    solver_message: str
    solver_nfev: int
    solver_cost: float
    solver_optimality: float
    active_bound_count: int
    initial_start_count: int
    initial_start_seed: int | None
    converged_root_count: int
    distinct_converged_root_count: int
    state: tuple[float, float, float, float, float, float, float]
    max_abs_raw_rhs: float
    pH_i: float
    cell_volume_pL: float
    cl_i_mM: float
    cl_target_z: float
    pH_target_z: float
    na_i_minus_published_mM: float
    k_i_minus_published_mM: float
    hco3_i_minus_published_mM: float
    volume_minus_published_pL: float
    quantified_rest_gate_pass: bool
    knockout_flow_target_read_by_solver: bool
    non_ae4_chassis_fingerprint: str
    failure_reason: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class StimulusComparison:
    """Source-bounded beta/cAMP/PKA activation prediction row."""

    mechanism_id: str
    core_candidate_id: str
    round: int
    family: str
    pka_stimulated_fold: float
    pka_activation_onset: float
    pka_fold_status: str
    timing_status: str
    resting_state_mode: str
    resting_state: str
    resting_max_abs_raw_rhs: float
    resting_cl_target_z: float
    resting_pH_target_z: float
    resting_quantified_gate_pass: bool
    assay_admissible: bool
    thermodynamic_status: str
    source_complete: bool
    source_complete_under_new_volume_sensitivity: bool
    n_fitted_ae4_parameters: int
    n_free_ae4_parameters: int
    calibration_target_names: str
    calibration_target_values: str
    independently_fixed_ae4_parameters: str
    new_modelling_assumptions: str
    knockout_target_read_by_calibrator: bool
    knockout_target_read_by_variant_selector: bool
    activity_scale: float
    non_ae4_chassis_fingerprint: str
    solver_method: str
    simulation_status: str
    wt_endpoint_flow: float
    ae4_ko_endpoint_flow: float
    ae2_ko_endpoint_flow: float
    ae4_ko_endpoint_flow_ratio: float
    ae2_ko_endpoint_flow_ratio: float
    ae4_ko_integrated_flow_ratio: float
    ae2_ko_integrated_flow_ratio: float
    ae4_ko_flow_ratio_at_t120_code_units: float
    ae4_ko_flow_ratio_at_t130_code_units: float
    wt_endpoint_cl_i_mM: float
    ae4_ko_endpoint_cl_i_mM: float
    ae2_ko_endpoint_cl_i_mM: float
    wt_endpoint_pH_i: float
    ae4_ko_endpoint_pH_i: float
    ae2_ko_endpoint_pH_i: float
    wt_endpoint_max_abs_raw_rhs: float
    ae4_ko_endpoint_max_abs_raw_rhs: float
    ae2_ko_endpoint_max_abs_raw_rhs: float
    max_conservation_residual: float
    max_thermodynamic_violation: float
    heldout_AE4_KO_ratio: float
    heldout_residual_endpoint: float
    heldout_residual_integrated: float
    failure_reason: str

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class BathOverrideMechanism:
    """Inject a common physical bath into AE4 without changing the chassis.

    The hard fixed-chassis comparison uses its documented 21 mM exchanger
    bath for every candidate.  The same adapter is used to carry 40 and
    42.9 mM as explicitly separate acid/base sensitivity branches; changing
    this AE4 context never mutates a non-AE4 equation.
    """

    mechanism: AE4Mechanism
    external_hco3_mM: float

    @property
    def spec(self) -> CandidateSpec:
        return self.mechanism.spec

    def evaluate(self, environment: AE4Environment, activity_scale: float = 1.0) -> Any:
        context = AE4Context(
            na_i=environment.na_i,
            k_i=environment.k_i,
            cl_i=environment.cl_i,
            hco3_i=environment.hco3_i,
            na_e=environment.na_e,
            k_e=environment.k_e,
            cl_e=environment.cl_e,
            hco3_e=self.external_hco3_mM,
        )
        return self.mechanism.evaluate(
            context,
            activity_scale=activity_scale,
            stimulus_activation=float(getattr(environment, "pka_activation", 0.0)),
        )


@dataclass(frozen=True)
class CommonPkaGateMechanism:
    """Apply the source-supported macroscopic PKA gate to any core law.

    The registered C1/C3 stimulus candidates implement this internally.  The
    adapter generalizes the common positive gate to reversible C6/C7 and
    saturating C8 cores not represented by a dedicated registry wrapper.  It
    preserves core affinity/sign and species signature while testing the
    independently observed positive activity fold.
    """

    mechanism: AE4Mechanism
    stimulated_fold: float

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.stimulated_fold)) or self.stimulated_fold < 1.0:
            raise ValueError("stimulated_fold must be finite and at least one")

    @property
    def spec(self) -> CandidateSpec:
        return self.mechanism.spec

    def evaluate(
        self,
        context: Any,
        activity_scale: float = 1.0,
        stimulus_activation: float | None = None,
    ) -> AE4Contribution:
        basal = self.mechanism.evaluate(
            context,
            activity_scale=activity_scale,
            stimulus_activation=0.0,
        )
        activation = (
            float(getattr(context, "pka_activation", 0.0))
            if stimulus_activation is None
            else float(stimulus_activation)
        )
        if not 0.0 <= activation <= 1.0:
            raise ValueError("stimulus_activation must lie in [0, 1]")
        multiplier = 1.0 + activation * (self.stimulated_fold - 1.0)
        diagnostics = dict(basal.diagnostics)
        for branch_key in ("na_flux", "k_flux"):
            if branch_key in diagnostics:
                diagnostics[branch_key] = multiplier * float(diagnostics[branch_key])
        diagnostics.update(
            {
                "pka_activation": activation,
                "pka_stimulated_fold": self.stimulated_fold,
                "pka_activity_multiplier": multiplier,
                "pka_interpretation": "common macroscopic activity multiplier adapter",
            }
        )
        return replace(
            basal,
            cycle_flux=multiplier * basal.cycle_flux,
            cl_i=multiplier * basal.cl_i,
            na_i=multiplier * basal.na_i,
            k_i=multiplier * basal.k_i,
            hco3_i=multiplier * basal.hco3_i,
            transported_charge=multiplier * basal.transported_charge,
            diagnostics=diagnostics,
        )


def contribution_vector(contribution: Any) -> np.ndarray:
    """Return an AE4 amount-balance vector in the canonical row order."""

    if isinstance(contribution, Mapping):
        return np.asarray([contribution[name] for name in BALANCE_NAMES], dtype=float)
    return np.asarray(
        [getattr(contribution, name) for name in BALANCE_NAMES], dtype=float
    )


def fit_nonnegative_activity_scale(
    unit_contribution: Sequence[float],
    target: Sequence[float] = HISTORICAL_WT_AE4_CONTRIBUTION,
) -> ScaleFit:
    """Fit one AE4 multiplier using only the fixed WT balance requirement.

    This is the closed-form solution of ``min_{a >= 0} ||a*s-t||_2``.  The
    target is a four-species WT balance vector; this function has no knockout
    argument and reads no phenotype ledger.
    """

    source = np.asarray(unit_contribution, dtype=float)
    wanted = np.asarray(target, dtype=float)
    if source.shape != (4,) or wanted.shape != (4,):
        raise ValueError("AE4 contribution and WT target must each have four rows")
    if not np.all(np.isfinite(source)) or not np.all(np.isfinite(wanted)):
        raise ValueError("AE4 calibration vectors must be finite")
    norm_sq = float(source @ source)
    if norm_sq <= 0.0:
        scale = 0.0
    else:
        scale = max(0.0, float(source @ wanted) / norm_sq)
    fitted = scale * source
    residual = fitted - wanted
    return ScaleFit(
        activity_scale=scale,
        fitted_contribution=tuple(float(value) for value in fitted),
        residual=tuple(float(value) for value in residual),
        l2_error=float(np.linalg.norm(residual)),
        max_abs_error=float(np.max(np.abs(residual))),
    )


def ratio(numerator: float, denominator: float) -> float:
    """A finite ratio helper that reports undefined comparisons as NaN."""

    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return math.nan
    return numerator / denominator


def integrated_flow(times: Sequence[float], flows: Sequence[float], start: float = 100.0) -> float:
    """Integrate stimulated code-native flow after the archived calcium step."""

    t = np.asarray(times, dtype=float)
    q = np.asarray(flows, dtype=float)
    mask = t > start
    if t.shape != q.shape or np.count_nonzero(mask) < 2:
        return math.nan
    return float(np.trapezoid(q[mask], t[mask]))


def cosine_alignment(target: Sequence[float], signature: Sequence[float]) -> float:
    """Absolute cosine used to rank balance pathways independent of sign."""

    lhs = np.asarray(target, dtype=float)
    rhs = np.asarray(signature, dtype=float)
    denom = float(np.linalg.norm(lhs) * np.linalg.norm(rhs))
    if denom == 0.0:
        return 0.0
    return abs(float(lhs @ rhs) / denom)


def projected_residual_fraction(target: Sequence[float], signature: Sequence[float]) -> float:
    """Residual norm after the best signed one-dimensional pathway correction."""

    lhs = np.asarray(target, dtype=float)
    rhs = np.asarray(signature, dtype=float)
    lhs_norm = float(np.linalg.norm(lhs))
    rhs_norm_sq = float(rhs @ rhs)
    if lhs_norm == 0.0:
        return 0.0
    if rhs_norm_sq == 0.0:
        return 1.0
    coefficient = float(lhs @ rhs) / rhs_norm_sq
    return float(np.linalg.norm(lhs - coefficient * rhs) / lhs_norm)


def stable_fingerprint(payload: Any) -> str:
    """Hash a JSON-serializable non-AE4 chassis description."""

    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def finite_max(values: Iterable[float], default: float = math.nan) -> float:
    finite = [abs(float(value)) for value in values if math.isfinite(float(value))]
    return max(finite) if finite else default


def baseline_context(external_hco3_mM: float = 40.0) -> AE4Context:
    """Source-labelled resting concentrations used for AE4-only calibration."""

    return AE4Context(
        na_i=float(BASELINE_STATE[3]),
        k_i=float(BASELINE_STATE[4]),
        cl_i=float(BASELINE_STATE[5]),
        hco3_i=float(BASELINE_STATE[6]),
        na_e=140.2,
        k_e=5.3,
        cl_e=102.6,
        hco3_e=float(external_hco3_mM),
    )


def independently_constrained_parameters(**overrides: float) -> AE4Parameters:
    """Return AE4 parameters fixed from the primary cation dose responses.

    The Na/K capacity weights use the measured maximal responses 1.5 and 1.6.
    EC50, Hill, and gate floors are the primary-paper fits.  The only fitted
    value in :func:`calibrate_candidate` is ``activity_scale``.
    """

    values: dict[str, float] = {
        "activity_scale": 1.0,
        "na_weight": 1.5,
        "k_weight": 1.6,
        "ec50_na": 49.0,
        "ec50_k": 62.0,
        "hill_na": 2.0,
        "hill_k": 1.8,
        "gate_floor_na": 0.3 / 1.5,
        "gate_floor_k": 0.4 / 1.6,
    }
    values.update({name: float(value) for name, value in overrides.items()})
    return AE4Parameters(**values)


def calibrate_candidate(
    candidate_id: str,
    *,
    parameters: AE4Parameters | None = None,
    external_hco3_mM: float | None = None,
) -> CalibratedCandidate:
    """Fit the AE4 capacity to WT closure without reading knockout data."""

    spec = CANDIDATE_SPECS[candidate_id]
    if not spec.active:
        raise ValueError(f"cannot calibrate analytically rejected {candidate_id}")
    bath_hco3 = 21.0
    if external_hco3_mM is not None:
        bath_hco3 = float(external_hco3_mM)
    unit_parameters = parameters or independently_constrained_parameters()
    unit_parameters = replace(unit_parameters, activity_scale=1.0)
    unit = evaluate_ae4(candidate_id, baseline_context(bath_hco3), unit_parameters)
    source = contribution_vector(unit)
    fit = fit_nonnegative_activity_scale(source)
    calibrated_parameters = replace(unit_parameters, activity_scale=fit.activity_scale)
    calibrated = evaluate_ae4(
        candidate_id, baseline_context(bath_hco3), calibrated_parameters
    )
    return CalibratedCandidate(
        candidate_id=candidate_id,
        parameters=calibrated_parameters,
        external_hco3_mM=bath_hco3,
        scale_fit=fit,
        unit_contribution=tuple(float(value) for value in source),
        baseline_diagnostics={
            "cycle_flux": float(calibrated.cycle_flux),
            "affinity": calibrated.affinity,
            "valid": bool(calibrated.valid),
            "transported_charge": float(calibrated.transported_charge),
            **dict(calibrated.diagnostics),
        },
    )


def candidate_chassis(calibration: CalibratedCandidate) -> FixedChassis:
    mechanism = AE4Mechanism(calibration.candidate_id, calibration.parameters)
    return FixedChassis(
        BathOverrideMechanism(mechanism, calibration.external_hco3_mM)
    )


def _reduced_rest_coordinates(state: Sequence[float], chassis: FixedChassis) -> np.ndarray:
    y = np.asarray(state, dtype=float)
    p = chassis.parameters
    h_i = y[6] + y[5] + p.impermeant_charge_amount / y[2] - y[3] - y[4]
    if not math.isfinite(h_i) or h_i <= 0.0:
        h_i = 10.0 ** (3.0 - WT_REST_PH)
    return np.array((y[0], y[1], y[3], y[4], y[5], y[6], math.log(h_i)))


def _state_from_reduced_rest(
    reduced: Sequence[float], chassis: FixedChassis
) -> np.ndarray:
    na_l, k_l, na_i, k_i, cl_i, hco3_i, log_h_i = map(float, reduced)
    h_i = math.exp(log_h_i)
    denominator = na_i + k_i + h_i - cl_i - hco3_i
    height = chassis.parameters.impermeant_charge_amount / denominator
    return np.array((na_l, k_l, height, na_i, k_i, cl_i, hco3_i))


def _rest_multistart_states(
    chassis: FixedChassis,
    *,
    preferred: Sequence[float] | None = None,
    seed: int = REST_MULTISTART_SEED,
    count: int = REST_MULTISTART_COUNT,
) -> tuple[np.ndarray, ...]:
    """Create a deterministic, knockout-blind broad resting-state multistart.

    The draw is made in the reduced electroneutral coordinate system.  It is
    a basin-robustness check, not a probability distribution over physiology
    and not a claim of global root uniqueness.
    """

    if count < 1:
        raise ValueError("rest multistart count must be positive")
    starts: list[np.ndarray] = [np.asarray(BASELINE_STATE, dtype=float)]
    if preferred is not None and not np.allclose(preferred, BASELINE_STATE):
        starts.append(np.asarray(preferred, dtype=float))
    generator = np.random.default_rng(seed)
    while len(starts) < count:
        p_h = float(generator.uniform(6.75, 7.40))
        reduced = np.array(
            (
                generator.uniform(90.0, 160.0),
                generator.uniform(1.0, 15.0),
                generator.uniform(10.0, 60.0),
                generator.uniform(80.0, 150.0),
                generator.uniform(25.0, 75.0),
                generator.uniform(5.0, 45.0),
                math.log(10.0 ** (3.0 - p_h)),
            ),
            dtype=float,
        )
        state = _state_from_reduced_rest(reduced, chassis)
        if np.all(np.isfinite(state)) and np.all(state > 0.0) and 2.0 < state[2] < 200.0:
            starts.append(state)
    return tuple(starts)


def solve_resting_equilibrium(
    chassis: FixedChassis,
    *,
    mechanism_id: str,
    activity_scale: float,
    scenario: Any = WT,
    initial_states: Sequence[Sequence[float]] = (BASELINE_STATE,),
    initial_start_seed: int | None = None,
    calibration_mode: str = "WT-only projected capacity; seven resting RHS equations",
) -> RestingEquilibrium:
    """Solve the seven resting equations in positivity-preserving coordinates.

    Cell height is eliminated through exact electroneutrality using a positive
    log-H+ coordinate.  Capacity is frozen before the solve; no knockout flow
    target, stimulated trajectory, or held-out ratio is available here.
    """

    best: Any | None = None
    fits: list[Any] = []

    def residual(reduced: np.ndarray) -> np.ndarray:
        state = _state_from_reduced_rest(reduced, chassis)
        if not np.all(np.isfinite(state)) or not 2.0 < state[2] < 200.0:
            return np.full(7, 1e3)
        try:
            raw = chassis.raw_rhs(0.0, state, scenario=scenario)
        except (ValueError, FloatingPointError, OverflowError):
            return np.full(7, 1e3)
        return np.asarray(raw, dtype=float) / REST_RHS_SCALES

    starts = tuple(initial_states) or (BASELINE_STATE,)
    for initial in starts:
        x0 = np.clip(
            _reduced_rest_coordinates(initial, chassis),
            REST_REDUCED_LOWER + 1e-12,
            REST_REDUCED_UPPER - 1e-12,
        )
        fit = least_squares(
            residual,
            x0,
            bounds=(REST_REDUCED_LOWER, REST_REDUCED_UPPER),
            x_scale="jac",
            diff_step=1e-7,
            max_nfev=10_000,
            xtol=1e-14,
            ftol=1e-14,
            gtol=1e-14,
        )
        fits.append(fit)
        if best is None or fit.cost < best.cost:
            best = fit
    assert best is not None
    converged_states: list[np.ndarray] = []
    for fit in fits:
        candidate_state = _state_from_reduced_rest(fit.x, chassis)
        if not fit.success or np.count_nonzero(fit.active_mask):
            continue
        try:
            candidate_raw = chassis.raw_rhs(0.0, candidate_state, scenario=scenario)
        except (ValueError, FloatingPointError, OverflowError):
            continue
        if finite_max(candidate_raw) <= 1e-8:
            converged_states.append(candidate_state)
    distinct_states: list[np.ndarray] = []
    for candidate_state in converged_states:
        if not any(
            np.allclose(candidate_state, known, rtol=1e-7, atol=1e-5)
            for known in distinct_states
        ):
            distinct_states.append(candidate_state)
    state = _state_from_reduced_rest(best.x, chassis)
    failure = ""
    try:
        raw, diagnostics = chassis.evaluate(
            state,
            calcium=0.05,
            pka_activation=0.0,
            scenario=scenario,
        )
        raw_error = finite_max(raw)
        ph_i = float(diagnostics["pH_i"])
        volume = float(diagnostics["cell_volume_pL"])
    except (ValueError, FloatingPointError, OverflowError) as exc:
        raw_error = ph_i = volume = math.nan
        diagnostics = {}
        failure = f"{type(exc).__name__}: {exc}"
    cl_z = (float(state[5]) - WT_REST_CL_MILLIMOLAR) / WT_REST_CL_SEM_MILLIMOLAR
    ph_z = (ph_i - WT_REST_PH) / WT_REST_PH_SEM
    quantified_gate = bool(
        best.success
        and int(np.count_nonzero(best.active_mask)) == 0
        and math.isfinite(raw_error)
        and raw_error <= 1e-8
        and abs(cl_z) <= 2.0
        and abs(ph_z) <= 2.0
        and abs(volume / PUBLISHED_REST_VOLUME_PL - 1.0)
        <= REST_VOLUME_RELATIVE_TOLERANCE
    )
    if not quantified_gate and not failure:
        reasons = []
        if raw_error > 1e-8:
            reasons.append(f"rest closure residual {raw_error:.6g} exceeds 1e-8")
        if abs(cl_z) > 2.0:
            reasons.append(f"WT resting Cl is {cl_z:.3f} SEM from target")
        if abs(ph_z) > 2.0:
            reasons.append(f"WT resting pH is {ph_z:.3f} SEM from target")
        if abs(volume / PUBLISHED_REST_VOLUME_PL - 1.0) > REST_VOLUME_RELATIVE_TOLERANCE:
            reasons.append(
                "WT resting volume lies outside the explicit +/-10% sensitivity "
                "around the published 1.3 pL convention"
            )
        if np.count_nonzero(best.active_mask):
            reasons.append("rest root touches a declared state bound")
        failure = "; ".join(reasons)
    return RestingEquilibrium(
        mechanism_id=mechanism_id,
        scenario=str(scenario.name),
        activity_scale=float(activity_scale),
        calibration_mode=calibration_mode,
        solver_success=bool(best.success),
        solver_message=str(best.message),
        solver_nfev=int(best.nfev),
        solver_cost=float(best.cost),
        solver_optimality=float(best.optimality),
        active_bound_count=int(np.count_nonzero(best.active_mask)),
        initial_start_count=len(starts),
        initial_start_seed=initial_start_seed,
        converged_root_count=len(converged_states),
        distinct_converged_root_count=len(distinct_states),
        state=tuple(float(value) for value in state),
        max_abs_raw_rhs=raw_error,
        pH_i=ph_i,
        cell_volume_pL=volume,
        cl_i_mM=float(state[5]),
        cl_target_z=cl_z,
        pH_target_z=ph_z,
        na_i_minus_published_mM=float(state[3] - PUBLISHED_REST_NA_MILLIMOLAR),
        k_i_minus_published_mM=float(state[4] - PUBLISHED_REST_K_MILLIMOLAR),
        hco3_i_minus_published_mM=float(
            state[6] - PUBLISHED_REST_HCO3_MILLIMOLAR
        ),
        volume_minus_published_pL=volume - PUBLISHED_REST_VOLUME_PL,
        quantified_rest_gate_pass=quantified_gate,
        knockout_flow_target_read_by_solver=False,
        non_ae4_chassis_fingerprint=chassis.parameter_fingerprint,
        failure_reason=failure,
    )


C6_VOLUME_CAPACITY_BRACKETS: Mapping[str, tuple[float, float]] = {
    "C6_112": (3.0e-8, 1.0e-7),
    "C6_123": (7.94328234724e-12, 1.0e-11),
    "C6_213": (1.99526231497e-11, 2.51188643151e-11),
}

C6_VOLUME_ROOT_SEEDS: Mapping[str, tuple[float, ...]] = {
    "C6_112": (
        118.76431077024066,
        5.619198735897633,
        28.700000000023827,
        33.409182950989795,
        111.66717301586701,
        43.47645503518593,
        16.599876242763532,
    ),
    "C6_123": (
        118.78620189797043,
        5.62570369237288,
        28.700000000011226,
        37.82898649650798,
        107.2733321441732,
        40.67477276443607,
        19.427513726276604,
    ),
    "C6_213": (
        118.78254683574043,
        5.624618651576397,
        28.700000000000415,
        31.727496347856736,
        113.37047307594486,
        45.7388045997583,
        14.359147449546358,
    ),
}


def volume_calibrated_c6_root(candidate_id: str) -> RestingEquilibrium:
    """Derive a C6 capacity from a KO-blind 1.3-pL resting-volume target.

    An outer Brent solve acts on log10 capacity.  Every objective evaluation
    solves the seven resting equations in the reduced electroneutral state.
    The capacity brackets were declared from a coarse continuation before any
    knockout simulation.  The final root is then challenged with the common
    deterministic 20-start basin screen.
    """

    lower, upper = C6_VOLUME_CAPACITY_BRACKETS[candidate_id]
    preferred = C6_VOLUME_ROOT_SEEDS[candidate_id]

    def root_at(capacity: float, starts: Sequence[Sequence[float]]) -> RestingEquilibrium:
        parameters = independently_constrained_parameters(activity_scale=capacity)
        chassis = FixedChassis(
            BathOverrideMechanism(AE4Mechanism(candidate_id, parameters), 21.0)
        )
        return solve_resting_equilibrium(
            chassis,
            mechanism_id=f"{candidate_id}_VOLUME_PROFILE",
            activity_scale=capacity,
            initial_states=starts,
            calibration_mode="KO-blind log-capacity continuation",
        )

    def objective(log10_capacity: float) -> float:
        root = root_at(10.0**log10_capacity, (preferred, BASELINE_STATE))
        if not root.solver_success or root.max_abs_raw_rhs > 1e-8:
            raise RuntimeError(f"{candidate_id} volume-profile rest solve failed")
        return root.cell_volume_pL - PUBLISHED_REST_VOLUME_PL

    log10_capacity = brentq(
        objective,
        math.log10(lower),
        math.log10(upper),
        xtol=1e-12,
        rtol=1e-14,
    )
    capacity = 10.0**log10_capacity
    parameters = independently_constrained_parameters(activity_scale=capacity)
    chassis = FixedChassis(
        BathOverrideMechanism(AE4Mechanism(candidate_id, parameters), 21.0)
    )
    starts = _rest_multistart_states(chassis, preferred=preferred)
    return solve_resting_equilibrium(
        chassis,
        mechanism_id=f"{candidate_id}_VOLUME_CALIBRATED",
        activity_scale=capacity,
        initial_states=starts,
        initial_start_seed=REST_MULTISTART_SEED,
        calibration_mode=(
            "seven WT resting RHS equations; AE4 capacity derived by KO-blind "
            f"Brent continuation over [{lower:.12g}, {upper:.12g}] to the "
            "published 1.3 pL cell volume; no KO flow datum"
        ),
    )


GENERIC_VOLUME_CAPACITY_INITIALS: Mapping[str, tuple[float, ...]] = {
    "C2": (1e-8, 1e-5, 1e-2),
    "C3": (1e-8, 1e-5, 1e-2),
    "C4": (1e-8, 1e-5, 1e-2),
    "C5": (1e-10, 1e-7, 1e-4),
    "C7_112": (1e-4, 1e-2, 1.0),
    "C7_123": (1e-4, 1e-2, 1.0),
    "C7_213": (1e-4, 1e-2, 1.0),
    "C8_112_SAT": (1e-3, 1e-1, 10.0),
    "C8_123_SAT": (1e-3, 1e-1, 10.0),
    "C8_213_SAT": (1e-3, 1e-1, 10.0),
}


def volume_calibrated_candidate_root(candidate_id: str) -> RestingEquilibrium:
    """Jointly solve seven WT rest equations and the 1.3-pL volume target.

    This is a KO-blind eight-variable solve over the seven reduced resting
    coordinates and log10 AE4 capacity.  It supplies the moved-root screen for
    candidates without a convenient monotone outer-capacity continuation.
    The final capacity is frozen and its state root is independently repeated
    from the common 20-start basin set.
    """

    base = FixedChassis()
    lower = np.concatenate((REST_REDUCED_LOWER, np.array((-20.0,))))
    upper = np.concatenate((REST_REDUCED_UPPER, np.array((5.0,))))
    target_height = PUBLISHED_REST_VOLUME_PL / (
        base.parameters.area_to_litre_per_height * 1e12
    )

    def residual(coordinates: np.ndarray) -> np.ndarray:
        capacity = 10.0 ** float(coordinates[-1])
        parameters = independently_constrained_parameters(activity_scale=capacity)
        chassis = FixedChassis(
            BathOverrideMechanism(AE4Mechanism(candidate_id, parameters), 21.0)
        )
        state = _state_from_reduced_rest(coordinates[:-1], chassis)
        if not np.all(np.isfinite(state)) or not 2.0 < state[2] < 200.0:
            return np.full(8, 1e3)
        try:
            raw = chassis.raw_rhs(0.0, state, scenario=WT)
        except (ValueError, FloatingPointError, OverflowError):
            return np.full(8, 1e3)
        return np.concatenate(
            (
                np.asarray(raw, dtype=float) / REST_RHS_SCALES,
                np.array(((state[2] - target_height) / 0.1,)),
            )
        )

    reduced_baseline = _reduced_rest_coordinates(BASELINE_STATE, base)
    fits: list[Any] = []
    for capacity in GENERIC_VOLUME_CAPACITY_INITIALS[candidate_id]:
        start = np.concatenate((reduced_baseline, np.array((math.log10(capacity),))))
        fits.append(
            least_squares(
                residual,
                start,
                bounds=(lower, upper),
                x_scale="jac",
                diff_step=1e-7,
                max_nfev=10_000,
                xtol=1e-14,
                ftol=1e-14,
                gtol=1e-14,
            )
        )
    fit = min(fits, key=lambda value: value.cost)
    capacity = 10.0 ** float(fit.x[-1])
    parameters = independently_constrained_parameters(activity_scale=capacity)
    chassis = FixedChassis(
        BathOverrideMechanism(AE4Mechanism(candidate_id, parameters), 21.0)
    )
    preferred = _state_from_reduced_rest(fit.x[:-1], chassis)
    starts = _rest_multistart_states(chassis, preferred=preferred)
    return solve_resting_equilibrium(
        chassis,
        mechanism_id=f"{candidate_id}_VOLUME_CALIBRATED",
        activity_scale=capacity,
        initial_states=starts,
        initial_start_seed=REST_MULTISTART_SEED,
        calibration_mode=(
            "KO-blind joint solve of seven WT resting RHS equations plus "
            "published 1.3 pL volume for AE4 capacity; initial capacity grid="
            f"{list(GENERIC_VOLUME_CAPACITY_INITIALS[candidate_id])}; final "
            "capacity frozen for deterministic broad 20-start basin screen"
        ),
    )


def resting_equilibrium_screens() -> tuple[
    list[RestingEquilibrium], dict[str, RestingEquilibrium]
]:
    """Screen fixed-capacity and KO-blind volume-calibrated family roots.

    C8 mass-action rows are exact evaluator aliases of their corresponding C6
    rows and therefore are not counted as independent moving-root evidence.
    """

    rows: list[RestingEquilibrium] = []
    roots: dict[str, RestingEquilibrium] = {}
    for candidate_id in ("C1", "C2", "C3"):
        calibration = calibrate_candidate(candidate_id)
        chassis = candidate_chassis(calibration)
        starts = _rest_multistart_states(chassis)
        root = solve_resting_equilibrium(
            chassis,
            mechanism_id=candidate_id,
            activity_scale=calibration.parameters.activity_scale,
            initial_states=starts,
            initial_start_seed=REST_MULTISTART_SEED,
            calibration_mode=(
                "WT-only projected capacity; seven resting RHS equations; "
                "deterministic broad 20-start basin screen"
            ),
        )
        rows.append(root)
        roots[candidate_id] = root

    c1_chassis = candidate_chassis(calibrate_candidate("C1"))
    knockout_starts = _rest_multistart_states(c1_chassis)
    knockout_root = solve_resting_equilibrium(
        c1_chassis,
        mechanism_id="COMMON_AE4_KO_REST",
        activity_scale=0.0,
        scenario=AE4_KNOCKOUT,
        initial_states=knockout_starts,
        initial_start_seed=REST_MULTISTART_SEED,
        calibration_mode=(
            "genotype-specific AE4-null seven-equation rest audit; no phenotype "
            "target and no parameter fit; deterministic broad 20-start screen"
        ),
    )
    rows.append(knockout_root)
    roots["AE4_KO"] = knockout_root
    ae2_knockout_root = solve_resting_equilibrium(
        c1_chassis,
        mechanism_id="C1_AE2_KO_REST",
        activity_scale=calibrate_candidate("C1").parameters.activity_scale,
        scenario=AE2_KNOCKOUT,
        initial_states=_rest_multistart_states(c1_chassis),
        initial_start_seed=REST_MULTISTART_SEED,
        calibration_mode=(
            "genotype-specific AE2-null seven-equation rest audit; no phenotype "
            "target and no parameter fit; deterministic broad 20-start screen"
        ),
    )
    rows.append(ae2_knockout_root)
    roots["C1_AE2_KO"] = ae2_knockout_root

    for candidate_id in ("C6_112", "C6_123", "C6_213"):
        root = volume_calibrated_c6_root(candidate_id)
        rows.append(root)
        roots[candidate_id] = root

    for candidate_id in (
        "C2",
        "C3",
        "C4",
        "C5",
        "C7_112",
        "C7_123",
        "C7_213",
        "C8_112_SAT",
        "C8_123_SAT",
        "C8_213_SAT",
    ):
        root = volume_calibrated_candidate_root(candidate_id)
        rows.append(root)
        roots[f"{candidate_id}_VOLUME"] = root
    return rows, roots


CAPACITY_PROFILE_LOG10_OFFSETS = tuple(
    sorted(
        {
            -8.0,
            -6.0,
            -4.0,
            -3.0,
            -2.0,
            -1.5,
            -1.25,
            *(
                round(float(value), 6)
                for value in np.arange(-1.0, 0.0 + 0.0125, 0.025)
            ),
            0.125,
            0.25,
            0.5,
            1.0,
            2.0,
            3.0,
            4.0,
        }
    )
)
CAPACITY_INTERVAL_SELECTION_FRACTIONS = (0.001, 0.25, 0.5, 0.75, 0.999)


def _capacity_gate_margin(root: RestingEquilibrium) -> float:
    """Normalized margin for KO-blind Cl, pH, and +/-10% volume gates."""

    return min(
        1.0 - abs(root.cl_target_z) / 2.0,
        1.0 - abs(root.pH_target_z) / 2.0,
        1.0
        - abs(root.cell_volume_pL / PUBLISHED_REST_VOLUME_PL - 1.0)
        / REST_VOLUME_RELATIVE_TOLERANCE,
    )


def _rest_at_core_capacity(
    candidate_id: str,
    capacity: float,
    initial_states: Sequence[Sequence[float]],
    *,
    mechanism_id: str,
) -> RestingEquilibrium:
    parameters = independently_constrained_parameters(activity_scale=capacity)
    chassis = FixedChassis(
        BathOverrideMechanism(AE4Mechanism(candidate_id, parameters), 21.0)
    )
    return solve_resting_equilibrium(
        chassis,
        mechanism_id=mechanism_id,
        activity_scale=capacity,
        initial_states=initial_states,
        calibration_mode=(
            "KO-blind capacity continuation using WT rest closure, primary Cl/pH "
            "constraints, and explicit +/-10% published-volume sensitivity"
        ),
    )


def rest_admissible_capacity_envelopes(
    roots: Mapping[str, RestingEquilibrium],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Locate rest-admissible capacity intervals without reading KO outcomes.

    A post-hoc but knockout-blind finite continuation spans 1e-8--1e4 times each independently
    volume-normalized capacity, with 0.025-decade resolution near the observed
    transition.  Sign changes in the joint resting gate are refined with Brent
    solves.  Five log-spaced interior points are frozen for later prediction;
    no knockout output enters interval discovery or point selection.  This is
    a tested finite envelope, not proof about unsampled remote branches.
    """

    candidate_ids = (
        "C2",
        "C3",
        "C4",
        "C5",
        "C6_112",
        "C6_123",
        "C6_213",
        "C7_112",
        "C7_123",
        "C7_213",
        "C8_112_SAT",
        "C8_123_SAT",
        "C8_213_SAT",
    )
    profile: list[dict[str, Any]] = []
    selections: list[dict[str, Any]] = []
    for candidate_id in candidate_ids:
        volume_key = (
            candidate_id
            if candidate_id.startswith("C6_")
            else f"{candidate_id}_VOLUME"
        )
        reference = roots[volume_key]
        reference_capacity = reference.activity_scale
        evaluated: list[tuple[float, RestingEquilibrium, float]] = []
        previous: Sequence[float] = reference.state
        for offset in CAPACITY_PROFILE_LOG10_OFFSETS:
            capacity = reference_capacity * 10.0**offset
            root = _rest_at_core_capacity(
                candidate_id,
                capacity,
                (previous, reference.state, BASELINE_STATE),
                mechanism_id=f"{candidate_id}_CAPACITY_PROFILE_{offset:+.3f}",
            )
            previous = root.state
            margin = _capacity_gate_margin(root)
            evaluated.append((offset, root, margin))
            record = root.as_dict()
            record.update(
                {
                    "core_candidate_id": candidate_id,
                    "reference_volume_capacity": reference_capacity,
                    "log10_capacity_offset": offset,
                    "capacity_gate_margin": margin,
                    "capacity_gate_pass": bool(
                        root.max_abs_raw_rhs <= 1e-8 and margin >= 0.0
                    ),
                    "profile_kind": "posthoc_KO_blind_finite_scan",
                    "knockout_target_read_by_variant_selector": False,
                }
            )
            profile.append(record)

        passing = [
            index
            for index, (_, root, margin) in enumerate(evaluated)
            if root.max_abs_raw_rhs <= 1e-8 and margin >= 0.0
        ]
        if not passing:
            continue
        # All observed pass sets are a single connected interval.  Refuse to
        # conceal disconnected regions if a later mechanism change creates one.
        if passing != list(range(min(passing), max(passing) + 1)):
            raise RuntimeError(f"disconnected rest-admissible interval for {candidate_id}")
        first = min(passing)
        last = max(passing)
        if first == 0 or last == len(evaluated) - 1:
            raise RuntimeError(
                f"rest-admissible interval touches finite scan boundary for {candidate_id}"
            )

        def margin_at(offset: float) -> float:
            capacity = reference_capacity * 10.0**offset
            root = _rest_at_core_capacity(
                candidate_id,
                capacity,
                (reference.state, BASELINE_STATE),
                mechanism_id=f"{candidate_id}_CAPACITY_BOUNDARY",
            )
            return _capacity_gate_margin(root)

        lower_offset = evaluated[first][0]
        if first > 0:
            lower_offset = brentq(
                margin_at,
                evaluated[first - 1][0],
                evaluated[first][0],
                xtol=1e-10,
                rtol=1e-12,
            )
        upper_offset = evaluated[last][0]
        if last + 1 < len(evaluated):
            upper_offset = brentq(
                margin_at,
                evaluated[last][0],
                evaluated[last + 1][0],
                xtol=1e-10,
                rtol=1e-12,
            )

        for fraction in CAPACITY_INTERVAL_SELECTION_FRACTIONS:
            offset = lower_offset + fraction * (upper_offset - lower_offset)
            capacity = reference_capacity * 10.0**offset
            preliminary = _rest_at_core_capacity(
                candidate_id,
                capacity,
                (reference.state, BASELINE_STATE),
                mechanism_id=f"{candidate_id}_REST_ADMISSIBLE_{fraction:.2f}",
            )
            parameters = independently_constrained_parameters(activity_scale=capacity)
            chassis = FixedChassis(
                BathOverrideMechanism(AE4Mechanism(candidate_id, parameters), 21.0)
            )
            starts = _rest_multistart_states(
                chassis,
                preferred=preliminary.state,
                count=10,
            )
            selected_root = solve_resting_equilibrium(
                chassis,
                mechanism_id=f"{candidate_id}_REST_ADMISSIBLE_{fraction:.2f}",
                activity_scale=capacity,
                initial_states=starts,
                initial_start_seed=REST_MULTISTART_SEED,
                calibration_mode=(
                    "KO-blind interior point of the connected rest-admissible interval "
                    "within the post-hoc KO-blind finite 1e-8--1e4x envelope "
                    "whose axes were frozen before held-out simulation; "
                    f"log-interval fraction={fraction:.3f}; "
                    "10-start basin screen; no KO phenotype available"
                ),
            )
            selected_margin = _capacity_gate_margin(selected_root)
            if (
                not selected_root.quantified_rest_gate_pass
                or selected_margin < -1e-10
            ):
                raise RuntimeError(
                    f"selected rest-admissible point failed its gate: {candidate_id}"
                )
            selection = {
                "core_candidate_id": candidate_id,
                "capacity_position": {
                    0.001: "near_lower_boundary",
                    0.25: "lower_quartile",
                    0.5: "midpoint",
                    0.75: "upper_quartile",
                    0.999: "near_upper_boundary",
                }[fraction],
                "capacity_interval_lower": reference_capacity * 10.0**lower_offset,
                "capacity_interval_upper": reference_capacity * 10.0**upper_offset,
                "capacity_log_interval_fraction": fraction,
                "reference_volume_capacity": reference_capacity,
                "root": selected_root,
                "knockout_target_read_by_variant_selector": False,
            }
            selections.append(selection)
            selected_record = selected_root.as_dict()
            selected_record.update(
                {
                    key: value for key, value in selection.items() if key != "root"
                }
            )
            selected_record.update(
                {
                    "log10_capacity_offset": offset,
                    "capacity_gate_margin": selected_margin,
                    "capacity_gate_pass": bool(selected_margin >= -1e-10),
                    "profile_kind": "frozen_interval_prediction_point",
                }
            )
            profile.append(selected_record)
    return profile, selections


def thermodynamic_rest_profile(
    volume_calibrated_capacity: float = 3.2298432383311024e-8,
) -> list[dict[str, Any]]:
    """Predeclared C6_112 capacity continuation exposing rest tradeoffs."""

    capacities = (
        1e-12,
        1e-10,
        1e-9,
        1e-8,
        3e-8,
        float(volume_calibrated_capacity),
        1e-7,
    )
    records: list[dict[str, Any]] = []
    previous: Sequence[float] = BASELINE_STATE
    for capacity in capacities:
        parameters = independently_constrained_parameters(activity_scale=capacity)
        chassis = FixedChassis(
            BathOverrideMechanism(AE4Mechanism("C6_112", parameters), 21.0)
        )
        root = solve_resting_equilibrium(
            chassis,
            mechanism_id=f"C6_112_G{capacity:.12g}",
            activity_scale=capacity,
            initial_states=(previous, BASELINE_STATE),
            calibration_mode="predeclared log-capacity thermodynamic rest continuation",
        )
        previous = root.state
        _, diagnostics = chassis.evaluate(
            root.state, calcium=0.05, pka_activation=0.0, scenario=WT
        )
        ae4_diag = diagnostics["ae4_diagnostics"]
        record = root.as_dict()
        record.update(
            {
                "na_branch_flux": ae4_diag.get("na_flux"),
                "k_branch_flux": ae4_diag.get("k_flux"),
                "na_affinity": ae4_diag.get("affinity_na"),
                "k_affinity": ae4_diag.get("affinity_k"),
                "net_cycle_flux": diagnostics["ae4_cycle_flux"],
                "thermodynamic_status": (
                    "pass"
                    if bool(ae4_diag.get("thermodynamic_valid", False))
                    else "fail"
                ),
                "volume_calibration_point": math.isclose(
                    capacity,
                    volume_calibrated_capacity,
                    rel_tol=1e-10,
                    abs_tol=1e-20,
                ),
            }
        )
        records.append(record)
    return records


def _trajectory_values(result: SimulationResult, name: str) -> np.ndarray:
    return np.asarray([float(row[name]) for row in result.diagnostics], dtype=float)


def _maximum_check_residual(results: Iterable[SimulationResult]) -> float:
    values: list[float] = []
    for result in results:
        for row in result.diagnostics:
            values.extend(float(value) for value in row["checks"].values())
    return finite_max(values)


def _maximum_thermodynamic_violation(
    results: Iterable[SimulationResult], spec: CandidateSpec
) -> float:
    if not spec.thermodynamic_sign_controlled:
        return math.nan
    violations: list[float] = []
    for result in results:
        for row in result.diagnostics:
            diagnostics = row["ae4_diagnostics"]
            for branch in ("na", "k"):
                flux = diagnostics.get(f"{branch}_flux")
                affinity = diagnostics.get(f"affinity_{branch}")
                if flux is not None and affinity is not None:
                    violations.append(max(0.0, -float(flux) * float(affinity)))
    return max(violations, default=0.0)


def _assay_admissible(spec: CandidateSpec, parameters: AE4Parameters) -> bool:
    # The direct primary flux experiments exclude the Na-only negative control.
    if spec.cation_mode == "Na only" or spec.stimulus_mode == "pka_k_recruit":
        return False
    # A zero-capacity branch is useful for Round-2 localization but contradicts
    # direct movement of both Na and K and is not fully assay-admissible.
    return parameters.na_weight > 0.0 and parameters.k_weight > 0.0


def _candidate_calibration_metadata(
    spec: CandidateSpec, calibration: CalibratedCandidate
) -> tuple[tuple[str, ...], dict[str, Any], dict[str, Any], list[str]]:
    targets = ["fixed_chassis_WT_AE4_balance_vector"]
    values: dict[str, Any] = {
        "WT_AE4_vector": list(HISTORICAL_WT_AE4_CONTRIBUTION),
        "AE4_bath_HCO3_mM": calibration.external_hco3_mM,
    }
    fixed: dict[str, Any] = {}
    additional_assumptions: list[str] = []
    core_id = spec.core_candidate_id or spec.candidate_id
    if core_id == "C1":
        fixed.update(
            {
                "k_forward": "historical 0.0192",
                "k_reverse": "historical folded 0.00159305881398359",
                "reverse_HCO3_mM": 21.0,
            }
        )
    elif core_id == "C2":
        fixed.update(
            {
                "k_forward": "2018 printed 0.0192",
                "k_reverse": "2018 printed 1.3e-5",
            }
        )
    else:
        assay_ratio = 1.5 / (1.5 + 1.6)
        actual_ratio = calibration.parameters.na_weight / (
            calibration.parameters.na_weight + calibration.parameters.k_weight
        )
        if math.isclose(actual_ratio, assay_ratio, rel_tol=0.0, abs_tol=1e-14):
            targets.append("primary_Na_K_maximal_response_capacity_ratio")
            values.update({"Na_Rmax": 1.5, "K_Rmax": 1.6})
            fixed.update({"na_weight": 1.5, "k_weight": 1.6})
        else:
            targets.append("predeclared_Round2_Na_K_capacity_interpolation")
            values.update(
                {
                    "na_capacity_fraction": actual_ratio,
                    "k_capacity_fraction": 1.0 - actual_ratio,
                }
            )
            additional_assumptions.append(
                "Na/K capacity fraction is a Round-2 localization axis, not an independent assay fit"
            )
    if core_id in ("C4", "C5"):
        targets.append("primary_cation_dose_response_shape")
        values.update(
            {
                "Na_EC50_mM": 49.0,
                "Na_Hill": 2.0,
                "K_EC50_mM": 62.0,
                "K_Hill": 1.8,
            }
        )
        fixed.update(
            {
                "ec50_na": 49.0,
                "hill_na": 2.0,
                "ec50_k": 62.0,
                "hill_k": 1.8,
            }
        )
    if core_id == "C5":
        values.update({"Na_Rmin": 0.3, "K_Rmin": 0.4})
        fixed.update({"gate_floor_na": 0.2, "gate_floor_k": 0.25})
    if spec.candidate_id.endswith("_SAT"):
        additional_assumptions.append(
            f"saturation_strength={calibration.parameters.saturation_strength:g} is an unconstrained modelling choice"
        )
    if spec.candidate_id.startswith(("C5", "C6", "C7", "C8")) and not math.isclose(
        calibration.parameters.equilibrium_constant, 1.0
    ):
        additional_assumptions.append(
            f"K_eq={calibration.parameters.equilibrium_constant:g} is an unsupported energetic-bias diagnostic"
        )
    if not math.isclose(calibration.external_hco3_mM, 21.0):
        additional_assumptions.append(
            f"AE4-only HCO3 bath={calibration.external_hco3_mM:g} mM is not a coherent whole-chassis bath"
        )
    if spec.stimulus_mode != "none":
        targets.append("primary_beta_cAMP_PKA_activity_fold")
        values.update(
            {
                "PKA_stimulated_fold": calibration.parameters.pka_stimulated_fold,
                "source_assay_envelope": [1.25, 1.6, 3.0],
            }
        )
        fixed["pka_stimulated_fold"] = calibration.parameters.pka_stimulated_fold
    return tuple(targets), values, fixed, additional_assumptions


def _effective_complexity(
    spec: CandidateSpec, calibration: CalibratedCandidate
) -> int:
    complexity = 1  # WT-fitted activity scale.
    core_id = spec.core_candidate_id or spec.candidate_id
    if core_id not in ("C1", "C2"):
        assay_ratio = 1.5 / (1.5 + 1.6)
        actual_ratio = calibration.parameters.na_weight / (
            calibration.parameters.na_weight + calibration.parameters.k_weight
        )
        if not math.isclose(actual_ratio, assay_ratio, rel_tol=0.0, abs_tol=1e-14):
            complexity += 1
    if spec.candidate_id.startswith(("C5", "C6", "C7", "C8")) and not math.isclose(
        calibration.parameters.equilibrium_constant, 1.0
    ):
        complexity += 1
    if spec.candidate_id.endswith("_SAT"):
        complexity += 1
    return complexity


def _density_source(spec: CandidateSpec, calibration: CalibratedCandidate) -> str:
    core_id = spec.core_candidate_id or spec.candidate_id
    if core_id == "C1":
        stimulus = (
            "; stimulated fold independently fixed from 2021 beta/cAMP assay"
            if spec.stimulus_mode != "none"
            else ""
        )
        return "WT balance activity scale; historical empirical turnover constants" + stimulus
    if core_id == "C2":
        return "WT balance activity scale; 2018 printed empirical turnover constants"
    assay_ratio = 1.5 / (1.5 + 1.6)
    actual_ratio = calibration.parameters.na_weight / (
        calibration.parameters.na_weight + calibration.parameters.k_weight
    )
    ratio_source = (
        "primary Rmax Na/K ratio"
        if math.isclose(actual_ratio, assay_ratio, rel_tol=0.0, abs_tol=1e-14)
        else "Round-2 Na/K interpolation"
    )
    return f"WT balance activity scale; branch ratio from {ratio_source}"


def _mechanism_alias(candidate_id: str) -> str:
    if candidate_id.startswith("C8_") and candidate_id.endswith("_MA"):
        return "C6_" + candidate_id.split("_")[1]
    return ""


def compare_calibrated_candidate(
    calibration: CalibratedCandidate,
    *,
    round_number: int = 1,
    mechanism_id: str | None = None,
    closure_tolerance: float = 1e-8,
    contribution_tolerance: float = 1e-8,
    solver_method: str = "Radau",
    precomputed_simulations: Mapping[str, SimulationResult] | None = None,
) -> tuple[CandidateComparison, Mapping[str, SimulationResult]]:
    """Freeze a calibration and simulate WT, AE4-KO, and AE2-KO."""

    spec = CANDIDATE_SPECS[calibration.candidate_id]
    chassis = candidate_chassis(calibration)
    baseline_raw, baseline = chassis.evaluate(BASELINE_STATE, calcium=0.05, scenario=WT)
    closure_error = finite_max(baseline_raw)
    closure_ok = (
        closure_error <= closure_tolerance
        and calibration.scale_fit.max_abs_error <= contribution_tolerance
    )
    simulations: dict[str, SimulationResult] = dict(precomputed_simulations or {})
    status = "pass"
    reason = ""
    try:
        if "wt" not in simulations:
            simulations["wt"] = chassis.simulate(
                WT, method=solver_method, t_eval=comparison_time_grid()
            )
        if "ae4_knockout" not in simulations:
            simulations["ae4_knockout"] = chassis.simulate(
                AE4_KNOCKOUT,
                method=solver_method,
                t_eval=comparison_time_grid(),
            )
        if "ae2_knockout" not in simulations:
            simulations["ae2_knockout"] = chassis.simulate(
                AE2_KNOCKOUT,
                method=solver_method,
                t_eval=comparison_time_grid(),
            )
    except (RuntimeError, ValueError, FloatingPointError, OverflowError) as exc:
        status = "failed"
        reason = f"{type(exc).__name__}: {exc}"
        simulations = {}

    nan = math.nan
    if simulations:
        wt = simulations["wt"]
        ae4_ko = simulations["ae4_knockout"]
        ae2_ko = simulations["ae2_knockout"]
        wt_endpoint = float(wt.endpoint["q_total"])
        ae4_endpoint = float(ae4_ko.endpoint["q_total"])
        ae2_endpoint = float(ae2_ko.endpoint["q_total"])
        wt_integral = integrated_flow(wt.time, _trajectory_values(wt, "q_total"))
        ae4_integral = integrated_flow(
            ae4_ko.time, _trajectory_values(ae4_ko, "q_total")
        )
        ae2_integral = integrated_flow(
            ae2_ko.time, _trajectory_values(ae2_ko, "q_total")
        )
        conservation = _maximum_check_residual(simulations.values())
        thermo_violation = _maximum_thermodynamic_violation(
            simulations.values(), spec
        )
        cl_ae4 = float(ae4_ko.state[-1, 5] - wt.state[-1, 5])
        cl_ae2 = float(ae2_ko.state[-1, 5] - wt.state[-1, 5])
        ph_ae4 = float(ae4_ko.endpoint["pH_i"] - wt.endpoint["pH_i"])
        ph_ae2 = float(ae2_ko.endpoint["pH_i"] - wt.endpoint["pH_i"])
        wt_endpoint_residual = finite_max(wt.endpoint["raw_rhs"].values())
        ae4_endpoint_residual = finite_max(
            ae4_ko.endpoint["raw_rhs"].values()
        )
        ae2_endpoint_residual = finite_max(
            ae2_ko.endpoint["raw_rhs"].values()
        )
    else:
        wt_endpoint = ae4_endpoint = ae2_endpoint = nan
        wt_integral = ae4_integral = ae2_integral = nan
        conservation = thermo_violation = cl_ae4 = cl_ae2 = ph_ae4 = ph_ae2 = nan
        wt_endpoint_residual = ae4_endpoint_residual = ae2_endpoint_residual = nan

    if not closure_ok:
        detail = (
            f"WT closure failed: raw={closure_error:.6g}, "
            f"AE4-vector={calibration.scale_fit.max_abs_error:.6g}"
        )
        reason = f"{reason}; {detail}".lstrip("; ")
    assay_ok = _assay_admissible(spec, calibration.parameters)
    if not assay_ok:
        detail = "direct primary Na/K transport evidence not satisfied"
        reason = f"{reason}; {detail}".lstrip("; ")
    thermo_status = (
        "pass"
        if spec.thermodynamic_sign_controlled and thermo_violation <= 1e-12
        else "fail"
        if spec.thermodynamic_sign_controlled
        else "empirical_not_sign_controlled"
    )

    (
        calibration_names,
        calibration_values,
        independently_fixed,
        additional_assumptions,
    ) = (
        _candidate_calibration_metadata(spec, calibration)
    )
    complexity = _effective_complexity(spec, calibration)
    row = CandidateComparison(
        mechanism_id=mechanism_id or calibration.candidate_id,
        round=round_number,
        family=spec.family,
        n_free_ae4_parameters=complexity,
        n_fitted_ae4_parameters=1,
        complexity_penalty=complexity,
        calibration_target_names=";".join(calibration_names),
        calibration_target_values=json.dumps(calibration_values, sort_keys=True),
        independently_fixed_ae4_parameters=json.dumps(
            independently_fixed, sort_keys=True
        ),
        new_modelling_assumptions="; ".join(
            (*spec.new_assumptions, *additional_assumptions)
        ),
        mechanism_equivalence_alias=_mechanism_alias(spec.candidate_id),
        heldout_target_names=";".join(HELDOUT_TARGET_NAMES),
        knockout_target_read_by_calibrator=False,
        activity_scale=calibration.parameters.activity_scale,
        ae4_density_or_turnover_source=_density_source(spec, calibration),
        non_ae4_chassis_fingerprint=chassis.parameter_fingerprint,
        units_convention="fixed historical chassis-native ratios; absolute flow/time uncertified",
        calibration_l2_error=calibration.scale_fit.l2_error,
        calibration_max_abs_error=calibration.scale_fit.max_abs_error,
        wt_baseline_closure_error=closure_error,
        wt_baseline_admissible=closure_ok,
        assay_admissible=assay_ok,
        thermodynamic_status=thermo_status,
        electroneutrality_status="pass" if bool(baseline["ae4_valid"]) else "fail",
        simulation_status=status,
        solver_method="+".join(
            sorted({result.method for result in simulations.values()})
        )
        if simulations
        else solver_method,
        wt_endpoint_flow=wt_endpoint,
        ae4_ko_endpoint_flow=ae4_endpoint,
        ae2_ko_endpoint_flow=ae2_endpoint,
        ae4_ko_endpoint_flow_ratio=ratio(ae4_endpoint, wt_endpoint),
        ae2_ko_endpoint_flow_ratio=ratio(ae2_endpoint, wt_endpoint),
        ae4_ko_integrated_flow_ratio=ratio(ae4_integral, wt_integral),
        ae2_ko_integrated_flow_ratio=ratio(ae2_integral, wt_integral),
        ae4_ko_cl_i_change_mM=cl_ae4,
        ae2_ko_cl_i_change_mM=cl_ae2,
        ae4_ko_pH_change=ph_ae4,
        ae2_ko_pH_change=ph_ae2,
        wt_endpoint_max_abs_raw_rhs=wt_endpoint_residual,
        ae4_ko_endpoint_max_abs_raw_rhs=ae4_endpoint_residual,
        ae2_ko_endpoint_max_abs_raw_rhs=ae2_endpoint_residual,
        endpoint_interpretation="t=200 protocol endpoint; not assumed to be an exact steady state",
        max_conservation_residual=conservation,
        max_thermodynamic_violation=thermo_violation,
        failure_reason=reason,
    )
    return row, simulations


def replay_trajectory(
    chassis: FixedChassis,
    state_source: SimulationResult,
    scenario: Any,
) -> SimulationResult:
    """Re-evaluate diagnostics on a trajectory known to have the same RHS.

    This is used only for exact-zero AE4 trajectories.  At AE4 activity zero,
    the state RHS is independent of the candidate law, so reintegration would
    repeat identical work.  Candidate-specific affinities and diagnostics are
    still recomputed at every saved state.
    """

    rows = []
    for time, state in zip(state_source.time, state_source.state, strict=True):
        _, diagnostic = chassis.evaluate(
            state,
            calcium=float(chassis.calcium_protocol(float(time))),
            pka_activation=float(chassis.pka_protocol(float(time))),
            scenario=scenario,
        )
        rows.append(diagnostic)
    return SimulationResult(
        time=state_source.time.copy(),
        state=state_source.state.copy(),
        diagnostics=tuple(rows),
        scenario=scenario,
        parameter_fingerprint=chassis.parameter_fingerprint,
        method=state_source.method,
        rtol=state_source.rtol,
        atol=state_source.atol,
    )


def primary_comparisons() -> tuple[list[CandidateComparison], dict[str, Mapping[str, SimulationResult]]]:
    """Run every active source-justified Round-1 registry candidate."""

    # Stimulus candidates have a separate independently bounded Round-3
    # protocol below.  The Round-1 static table contains only ungated cores.
    calibrations = [
        calibrate_candidate(spec.candidate_id)
        for spec in candidate_specs(include_rejected=False)
        if spec.stimulus_mode == "none"
    ]
    # The AE4-KO state trajectory is exactly mechanism-independent.  Likewise,
    # all candidates whose nonnegative WT fit is zero share the same WT and
    # AE2-KO state trajectories.  Integrate each such RHS once, then replay
    # candidate-specific diagnostics.
    reference = candidate_chassis(calibrations[0])
    # The source-matched BDF solver is canonical for the only WT-closed row and
    # for the common held-out KO trajectory.  Radau is a positivity-safe screen
    # for closure-failing candidates and is cross-checked against BDF below.
    shared_ae4_ko = reference.simulate(
        AE4_KNOCKOUT, method="BDF", t_eval=comparison_time_grid()
    )
    zero_calibration = next(
        (item for item in calibrations if item.parameters.activity_scale == 0.0),
        None,
    )
    shared_zero_wt: SimulationResult | None = None
    shared_zero_ae2: SimulationResult | None = None
    if zero_calibration is not None:
        zero_chassis = candidate_chassis(zero_calibration)
        # WT with zero calibrated capacity has the same RHS as AE4 KO.
        shared_zero_wt = replay_trajectory(zero_chassis, shared_ae4_ko, WT)
        shared_zero_ae2 = zero_chassis.simulate(
            AE2_KNOCKOUT, method="Radau", t_eval=comparison_time_grid()
        )

    rows: list[CandidateComparison] = []
    trajectories: dict[str, Mapping[str, SimulationResult]] = {}
    for calibration in calibrations:
        chassis = candidate_chassis(calibration)
        precomputed: dict[str, SimulationResult] = {
            "ae4_knockout": replay_trajectory(
                chassis, shared_ae4_ko, AE4_KNOCKOUT
            )
        }
        if calibration.parameters.activity_scale == 0.0:
            assert shared_zero_wt is not None and shared_zero_ae2 is not None
            precomputed["wt"] = replay_trajectory(chassis, shared_zero_wt, WT)
            precomputed["ae2_knockout"] = replay_trajectory(
                chassis, shared_zero_ae2, AE2_KNOCKOUT
            )
        solver = "BDF" if calibration.candidate_id == "C1" else "Radau"
        row, runs = compare_calibrated_candidate(
            calibration,
            round_number=1,
            solver_method=solver,
            precomputed_simulations=precomputed,
        )
        rows.append(row)
        trajectories[row.mechanism_id] = runs
    return rows, trajectories


def _step_activation(onset: float) -> Any:
    def protocol(time: float) -> float:
        return 0.0 if float(time) <= onset else 1.0

    return protocol


def _flow_at(result: SimulationResult, time: float) -> float:
    index = int(np.argmin(np.abs(np.asarray(result.time) - time)))
    if abs(float(result.time[index]) - time) > 1e-6:
        return math.nan
    return float(result.diagnostics[index]["q_total"])


def _stimulus_law_and_metadata(
    mechanism_id: str,
    fold: float,
    roots: Mapping[str, RestingEquilibrium],
    *,
    rest_override: RestingEquilibrium | None = None,
) -> tuple[Any, RestingEquilibrium, CandidateSpec, float, dict[str, Any]]:
    if mechanism_id in ("C1_PKA_COMMON", "C1_PKA_K_RECRUIT"):
        calibration = calibrate_candidate(
            mechanism_id,
            parameters=replace(
                independently_constrained_parameters(),
                pka_stimulated_fold=fold,
            ),
        )
        law: Any = AE4Mechanism(mechanism_id, calibration.parameters)
        rest = roots["C1"]
        core_id = "C1"
        assumptions = {
            "C1_PKA_COMMON": "Na-only common macroscopic activity gate",
            "C1_PKA_K_RECRUIT": (
                "novel selectively recruited K branch; direction-locked and "
                "thermodynamically uncertified"
            ),
        }[mechanism_id]
        return law, rest, CANDIDATE_SPECS[mechanism_id], calibration.parameters.activity_scale, {
            "core_id": core_id,
            "assumptions": assumptions,
            "calibration": "exact historical WT resting vector plus independent PKA fold",
        }
    if mechanism_id == "C3_PKA_COMMON":
        calibration = calibrate_candidate(
            mechanism_id,
            parameters=replace(
                independently_constrained_parameters(),
                pka_stimulated_fold=fold,
            ),
        )
        return (
            AE4Mechanism(mechanism_id, calibration.parameters),
            roots["C3"],
            CANDIDATE_SPECS[mechanism_id],
            calibration.parameters.activity_scale,
            {
                "core_id": "C3",
                "assumptions": "common macroscopic gate on both empirical Na/K branches",
                "calibration": "WT-only projected capacity, moved WT rest root, independent PKA fold",
            },
        )
    if mechanism_id in (
        "C6_112_PKA_COMMON",
        "C6_123_PKA_COMMON",
        "C6_213_PKA_COMMON",
    ):
        core_id = mechanism_id.removesuffix("_PKA_COMMON")
        rest = rest_override or roots[core_id]
        capacity = rest.activity_scale
        parameters = replace(
            independently_constrained_parameters(activity_scale=capacity),
            pka_stimulated_fold=fold,
        )
        if mechanism_id in CANDIDATE_SPECS:
            law = AE4Mechanism(mechanism_id, parameters)
            spec = CANDIDATE_SPECS[mechanism_id]
        else:
            base = AE4Mechanism(core_id, parameters)
            law = CommonPkaGateMechanism(base, fold)
            spec = CANDIDATE_SPECS[core_id]
        return (
            law,
            rest,
            spec,
            capacity,
            {
                "core_id": core_id,
                "assumptions": "common macroscopic gate preserving the reversible Keq=1 core",
                "calibration": (
                    "KO-blind rest-admissible capacity interval plus independent PKA fold"
                    if rest_override is not None
                    else "seven WT resting equations plus independent 1.3 pL volume and PKA fold"
                ),
                "skip_dynamic": rest_override is None
                and core_id in ("C6_123", "C6_213"),
            },
        )
    if mechanism_id in (
        "C7_213_PKA_COMMON",
        "C8_213_SAT_PKA_COMMON",
    ):
        if rest_override is None:
            raise ValueError(f"{mechanism_id} requires a frozen rest-interval point")
        core_id = mechanism_id.removesuffix("_PKA_COMMON")
        capacity = rest_override.activity_scale
        parameters = replace(
            independently_constrained_parameters(activity_scale=capacity),
            pka_stimulated_fold=fold,
        )
        base = AE4Mechanism(core_id, parameters)
        return (
            CommonPkaGateMechanism(base, fold),
            rest_override,
            CANDIDATE_SPECS[core_id],
            capacity,
            {
                "core_id": core_id,
                "assumptions": (
                    "common macroscopic PKA gate preserving the reversible Keq=1 "
                    "core; new mapping of the native activity assay onto this law"
                ),
                "calibration": (
                    "KO-blind rest-admissible capacity interval plus independent PKA fold"
                ),
                "skip_dynamic": False,
            },
        )
    raise KeyError(mechanism_id)


def stimulus_activation_comparisons() -> tuple[
    list[StimulusComparison],
    dict[str, Mapping[str, SimulationResult]],
    list[RestingEquilibrium],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    """Run the source-bounded PKA envelope after every calibration is frozen."""

    rest_rows, roots = resting_equilibrium_screens()
    thermo_profile = thermodynamic_rest_profile(roots["C6_112"].activity_scale)
    capacity_profile, interval_selections = rest_admissible_capacity_envelopes(roots)
    mechanisms = (
        "C1_PKA_COMMON",
        "C1_PKA_K_RECRUIT",
        "C3_PKA_COMMON",
        "C6_112_PKA_COMMON",
    )
    definitions: list[
        tuple[str, float, float, str | None, RestingEquilibrium | None]
    ] = [
        (mechanism_id, fold, PKA_IMMEDIATE_MODEL_ONSET, None, None)
        for mechanism_id in mechanisms
        for fold in (
            PKA_REFERENCE_FOLD,
            *PKA_ASSAY_FOLDS,
            PKA_INVERSE_DIAGNOSTIC_FOLD,
        )
    ]
    definitions.extend(
        (mechanism_id, 3.0, onset, None, None)
        for mechanism_id in mechanisms
        for onset in PKA_TIMING_SENSITIVITY_ONSETS
    )
    # These source-left-viable stoichiometries first receive the same KO-blind
    # volume-normalized moving-rest screen.  Their independent resting Cl gate
    # fails, so decision criterion 1 prohibits a stimulated phenotype run.
    definitions.extend(
        (
            ("C6_123_PKA_COMMON", 3.0, PKA_IMMEDIATE_MODEL_ONSET, None, None),
            ("C6_213_PKA_COMMON", 3.0, PKA_IMMEDIATE_MODEL_ONSET, None, None),
        )
    )
    stimulus_id_by_core = {
        "C6_213": "C6_213_PKA_COMMON",
        "C7_213": "C7_213_PKA_COMMON",
        "C8_213_SAT": "C8_213_SAT_PKA_COMMON",
    }
    for selection in interval_selections:
        core_id = str(selection["core_candidate_id"])
        if core_id not in stimulus_id_by_core:
            raise RuntimeError(
                f"rest-admissible source core lacks a PKA wrapper: {core_id}"
            )
        rest = selection["root"]
        assert isinstance(rest, RestingEquilibrium)
        rest_rows.append(rest)
        ae2_parameters = independently_constrained_parameters(
            activity_scale=rest.activity_scale
        )
        ae2_rest_chassis = FixedChassis(
            BathOverrideMechanism(AE4Mechanism(core_id, ae2_parameters), 21.0)
        )
        ae2_rest = solve_resting_equilibrium(
            ae2_rest_chassis,
            mechanism_id=(
                f"{core_id}_{selection['capacity_position']}_AE2_KO_REST"
            ),
            activity_scale=rest.activity_scale,
            scenario=AE2_KNOCKOUT,
            initial_states=(rest.state, roots["C1_AE2_KO"].state, BASELINE_STATE),
            calibration_mode=(
                "candidate-specific AE2-null rest prediction at a frozen KO-blind "
                "WT-rest capacity point; three local starts; no genotype target fit"
            ),
        )
        rest_rows.append(ae2_rest)
        label = f"capacity_{selection['capacity_position']}"
        definitions.extend(
            (
                stimulus_id_by_core[core_id],
                fold,
                PKA_IMMEDIATE_MODEL_ONSET,
                label,
                rest,
            )
            for fold in (
                PKA_REFERENCE_FOLD,
                *PKA_ASSAY_FOLDS,
                PKA_INVERSE_DIAGNOSTIC_FOLD,
            )
        )

    common_ko: dict[str, SimulationResult] = {}
    rows: list[StimulusComparison] = []
    trajectories: dict[str, Mapping[str, SimulationResult]] = {}
    for mechanism_id, fold, onset, capacity_label, rest_override in definitions:
        law, rest, spec, activity_scale, metadata = _stimulus_law_and_metadata(
            mechanism_id,
            fold,
            roots,
            rest_override=rest_override,
        )
        if rest_override is not None:
            selection = next(
                item
                for item in interval_selections
                if item["root"] is rest_override
            )
            metadata = {
                **metadata,
                "capacity_interval_lower": selection["capacity_interval_lower"],
                "capacity_interval_upper": selection["capacity_interval_upper"],
                "capacity_position": selection["capacity_position"],
                "capacity_log_interval_fraction": selection[
                    "capacity_log_interval_fraction"
                ],
                "assumptions": (
                    f"{metadata['assumptions']}; +/-10% volume is a new post-hoc "
                    "KO-blind robustness sensitivity, not source-reported uncertainty"
                ),
            }
        protocol = _step_activation(onset)
        chassis = FixedChassis(
            BathOverrideMechanism(law, 21.0),
            pka_protocol=protocol,
            pka_discontinuities=(onset,),
        )
        root_key = (
            f"{metadata['core_id']}:{rest.activity_scale:.17g}:"
            f"{rest.mechanism_id}"
        )
        capacity_token = "" if capacity_label is None else f"_{capacity_label}"
        variant_id = f"{mechanism_id}{capacity_token}_F{fold:g}_onset{onset:g}"
        skip_dynamic = bool(metadata.get("skip_dynamic", False))
        status = "not_run_failed_rest_gate" if skip_dynamic else "pass"
        reason = ""
        runs: dict[str, SimulationResult] = {}
        if skip_dynamic:
            reason = (
                "stimulated run withheld under decision criterion 1 after the "
                "KO-blind volume-calibrated root failed independent resting data"
            )
        else:
            if root_key not in common_ko:
                common_ko[root_key] = chassis.simulate(
                    AE4_KNOCKOUT,
                    method="Radau",
                    rtol=1e-7,
                    atol=1e-8,
                    t_eval=comparison_time_grid(),
                    y0=rest.state,
                )
            runs["ae4_knockout"] = replay_trajectory(
                chassis, common_ko[root_key], AE4_KNOCKOUT
            )
            try:
                runs["wt"] = chassis.simulate(
                    WT,
                    method="Radau",
                    rtol=1e-7,
                    atol=1e-8,
                    t_eval=comparison_time_grid(),
                    y0=rest.state,
                )
                if onset == PKA_IMMEDIATE_MODEL_ONSET and fold in (
                    PKA_REFERENCE_FOLD,
                    *PKA_ASSAY_FOLDS,
                ):
                    runs["ae2_knockout"] = chassis.simulate(
                        AE2_KNOCKOUT,
                        method="Radau",
                        rtol=1e-7,
                        atol=1e-8,
                        t_eval=comparison_time_grid(),
                        y0=rest.state,
                    )
            except (RuntimeError, ValueError, FloatingPointError, OverflowError) as exc:
                status = "failed"
                reason = f"{type(exc).__name__}: {exc}"

        nan = math.nan
        if status == "pass":
            wt = runs["wt"]
            ae4_ko = runs["ae4_knockout"]
            ae2_ko = runs.get("ae2_knockout")
            wt_endpoint = float(wt.endpoint["q_total"])
            ae4_endpoint = float(ae4_ko.endpoint["q_total"])
            ae2_endpoint = (
                math.nan if ae2_ko is None else float(ae2_ko.endpoint["q_total"])
            )
            wt_integral = integrated_flow(wt.time, _trajectory_values(wt, "q_total"))
            ae4_integral = integrated_flow(
                ae4_ko.time, _trajectory_values(ae4_ko, "q_total")
            )
            ae2_integral = (
                math.nan
                if ae2_ko is None
                else integrated_flow(
                    ae2_ko.time, _trajectory_values(ae2_ko, "q_total")
                )
            )
            endpoint_ratio = ratio(ae4_endpoint, wt_endpoint)
            integrated_ratio = ratio(ae4_integral, wt_integral)
            conservation = _maximum_check_residual(runs.values())
            thermo_violation = _maximum_thermodynamic_violation(
                runs.values(), spec
            )
            early20 = ratio(_flow_at(ae4_ko, 120.0), _flow_at(wt, 120.0))
            early30 = ratio(_flow_at(ae4_ko, 130.0), _flow_at(wt, 130.0))
            wt_raw = finite_max(wt.endpoint["raw_rhs"].values())
            ae4_raw = finite_max(ae4_ko.endpoint["raw_rhs"].values())
            ae2_raw = (
                math.nan
                if ae2_ko is None
                else finite_max(ae2_ko.endpoint["raw_rhs"].values())
            )
        else:
            wt_endpoint = ae4_endpoint = ae2_endpoint = nan
            wt_integral = ae4_integral = ae2_integral = nan
            endpoint_ratio = integrated_ratio = conservation = thermo_violation = nan
            early20 = early30 = wt_raw = ae4_raw = ae2_raw = nan
            wt = ae4_ko = ae2_ko = None

        fold_status = {
            1.0: "no_activation_reference",
            1.25: "2021_heterologous_forskolin_direct_measurement",
            1.6: "2021_total_WT_exchange_intermediate_sensitivity_not_AE4_specific",
            3.0: "2021_native_AE4_only_approximate_digitized_fold",
            10.0: "out_of_assay_inverse_diagnostic",
        }[fold]
        timing_status = (
            "model_input_immediate_step_at_t100_code_units"
            if onset == PKA_IMMEDIATE_MODEL_ONSET
            else "post_step_code_coordinate_onset_sensitivity_physical_time_uncertified"
        )
        assay_ok = bool(
            str(metadata["core_id"]) != "C1"
            and mechanism_id != "C1_PKA_K_RECRUIT"
        )
        thermo_status = (
            "pass"
            if spec.thermodynamic_sign_controlled
            and (
                (math.isfinite(thermo_violation) and thermo_violation <= 1e-12)
                or skip_dynamic
            )
            else "empirical_not_sign_controlled"
            if mechanism_id != "C1_PKA_K_RECRUIT"
            else "fail_direction_locked_K_recruitment"
        )
        source_complete_under_new_volume_sensitivity = bool(
            assay_ok
            and rest.quantified_rest_gate_pass
            and thermo_status == "pass"
            and fold in (1.25, 3.0)
        )
        source_complete = bool(
            source_complete_under_new_volume_sensitivity
            and rest_override is None
        )
        if not rest.quantified_rest_gate_pass:
            reason = f"{reason}; {rest.failure_reason}".strip("; ")
        if not assay_ok:
            reason = (
                f"{reason}; candidate does not satisfy the complete direct Na/K assay gate"
            ).strip("; ")
        row = StimulusComparison(
            mechanism_id=variant_id,
            core_candidate_id=str(metadata["core_id"]),
            round=3,
            family=mechanism_id,
            pka_stimulated_fold=fold,
            pka_activation_onset=onset,
            pka_fold_status=fold_status,
            timing_status=timing_status,
            resting_state_mode=(
                "common candidate-specific WT root at t=0; each genotype evolves "
                "independently for 100 code units before stimulation"
            ),
            resting_state=json.dumps(list(rest.state)),
            resting_max_abs_raw_rhs=rest.max_abs_raw_rhs,
            resting_cl_target_z=rest.cl_target_z,
            resting_pH_target_z=rest.pH_target_z,
            resting_quantified_gate_pass=rest.quantified_rest_gate_pass,
            assay_admissible=assay_ok,
            thermodynamic_status=thermo_status,
            source_complete=source_complete,
            source_complete_under_new_volume_sensitivity=(
                source_complete_under_new_volume_sensitivity
            ),
            n_fitted_ae4_parameters=1,
            n_free_ae4_parameters=1,
            calibration_target_names=(
                "WT_rest_closure;independent_rest_observables;primary_2021_PKA_activity_fold"
            ),
            calibration_target_values=json.dumps(
                {
                    "fold": fold,
                    "source_fold_envelope": list(PKA_ASSAY_FOLDS),
                    "rest_Cl_mM": WT_REST_CL_MILLIMOLAR,
                    "rest_pH": WT_REST_PH,
                    "rest_volume_pL": PUBLISHED_REST_VOLUME_PL,
                    "rest_volume_relative_tolerance": REST_VOLUME_RELATIVE_TOLERANCE,
                    "capacity_interval_lower": metadata.get(
                        "capacity_interval_lower"
                    ),
                    "capacity_interval_upper": metadata.get(
                        "capacity_interval_upper"
                    ),
                    "capacity_position": metadata.get("capacity_position"),
                },
                sort_keys=True,
            ),
            independently_fixed_ae4_parameters=str(metadata["calibration"]),
            new_modelling_assumptions=str(metadata["assumptions"]),
            knockout_target_read_by_calibrator=False,
            knockout_target_read_by_variant_selector=False,
            activity_scale=activity_scale,
            non_ae4_chassis_fingerprint=chassis.parameter_fingerprint,
            solver_method="Radau rtol=1e-7 atol=1e-8",
            simulation_status=status,
            wt_endpoint_flow=wt_endpoint,
            ae4_ko_endpoint_flow=ae4_endpoint,
            ae2_ko_endpoint_flow=ae2_endpoint,
            ae4_ko_endpoint_flow_ratio=endpoint_ratio,
            ae2_ko_endpoint_flow_ratio=ratio(ae2_endpoint, wt_endpoint),
            ae4_ko_integrated_flow_ratio=integrated_ratio,
            ae2_ko_integrated_flow_ratio=ratio(ae2_integral, wt_integral),
            ae4_ko_flow_ratio_at_t120_code_units=early20,
            ae4_ko_flow_ratio_at_t130_code_units=early30,
            wt_endpoint_cl_i_mM=nan if wt is None else float(wt.state[-1, 5]),
            ae4_ko_endpoint_cl_i_mM=nan if ae4_ko is None else float(ae4_ko.state[-1, 5]),
            ae2_ko_endpoint_cl_i_mM=nan if ae2_ko is None else float(ae2_ko.state[-1, 5]),
            wt_endpoint_pH_i=nan if wt is None else float(wt.endpoint["pH_i"]),
            ae4_ko_endpoint_pH_i=nan if ae4_ko is None else float(ae4_ko.endpoint["pH_i"]),
            ae2_ko_endpoint_pH_i=nan if ae2_ko is None else float(ae2_ko.endpoint["pH_i"]),
            wt_endpoint_max_abs_raw_rhs=wt_raw,
            ae4_ko_endpoint_max_abs_raw_rhs=ae4_raw,
            ae2_ko_endpoint_max_abs_raw_rhs=ae2_raw,
            max_conservation_residual=conservation,
            max_thermodynamic_violation=thermo_violation,
            heldout_AE4_KO_ratio=math.nan,
            heldout_residual_endpoint=math.nan,
            heldout_residual_integrated=math.nan,
            failure_reason=reason,
        )
        rows.append(row)
        if runs:
            trajectories[variant_id] = runs
    return rows, trajectories, rest_rows, thermo_profile, capacity_profile


PATHWAY_SIGNATURES: Mapping[str, tuple[float, float, float, float]] = {
    "NaK_ATPase": (-3.0, 2.0, 0.0, 0.0),
    "NHE1": (1.0, 0.0, 0.0, 0.0),
    "CaKC_K_efflux": (0.0, -1.0, 0.0, 0.0),
    "NKCC1": (1.0, 1.0, 2.0, 0.0),
    "AE2": (0.0, 0.0, 1.0, -1.0),
    "buffer_acid_base": (0.0, 0.0, 0.0, 1.0),
    "volume_dilution_at_WT": (-25.0, -120.0, -50.0, -10.0),
}


PATHWAY_MISMATCH_EVIDENCE: Mapping[str, str] = {
    "NaK_ATPase": "historical M basis recovered but published density/whole-cell scale unresolved",
    "NHE1": "historical first-power proton law differs from published squared law",
    "CaKC_K_efflux": "historical gates/conductance differ materially from the published values",
    "NKCC1": "coefficient convention corrected, but density/whole-cell scale remains unresolved",
    "AE2": "historical 1e-4 mM constant conflicts with the printed 1e4 mM value",
    "buffer_acid_base": "historical CO2 QSS/sign and fitted buffer scaling differ from published dynamics",
    "volume_dilution_at_WT": "historical water coefficients are fitted and not the published permeability set",
}


def pathway_localization(
    correction_vector: Sequence[float],
) -> list[dict[str, Any]]:
    """Rank one-component balance directions without changing the chassis."""

    target = np.asarray(correction_vector, dtype=float)
    rows = []
    for pathway, signature in PATHWAY_SIGNATURES.items():
        rows.append(
            {
                "pathway": pathway,
                "correction_vector": [float(value) for value in target],
                "signature": list(signature),
                "absolute_cosine_alignment": cosine_alignment(target, signature),
                "best_one_component_residual_fraction": projected_residual_fraction(
                    target, signature
                ),
                "historical_vs_published_mismatch": PATHWAY_MISMATCH_EVIDENCE[pathway],
            }
        )
    rows.sort(
        key=lambda row: (
            -float(row["absolute_cosine_alignment"]),
            float(row["best_one_component_residual_fraction"]),
            str(row["pathway"]),
        )
    )
    for index, row in enumerate(rows, start=1):
        row["rank"] = index
    return rows


ROUND2_NA_CAPACITY_FRACTIONS = (
    0.0,
    0.25,
    1.5 / (1.5 + 1.6),
    0.5,
    0.75,
    1.0,
)
ROUND2_EQUILIBRIUM_CONSTANTS = (0.01, 0.1, 1.0, 10.0, 100.0)
ROUND2_SATURATION_STRENGTHS = (0.01, 0.1, 1.0, 10.0, 100.0)
ROUND2_BATH_HCO3_MILLIMOLAR = (21.0, 40.0, 42.9)


def _round2_variant_definitions() -> list[dict[str, Any]]:
    """Predeclared structural grid; it contains no phenotype value."""

    variants: list[dict[str, Any]] = []

    def add(
        candidate_id: str,
        *,
        f_na: float,
        equilibrium_constant: float = 1.0,
        saturation_strength: float = 1.0,
        bath_hco3_mM: float = 21.0,
        axis: str,
    ) -> None:
        variants.append(
            {
                "candidate_id": candidate_id,
                "na_capacity_fraction": float(f_na),
                "k_capacity_fraction": float(1.0 - f_na),
                "equilibrium_constant": float(equilibrium_constant),
                "saturation_strength": float(saturation_strength),
                "bath_hco3_mM": float(bath_hco3_mM),
                "grid_axis": axis,
            }
        )

    for f_na in ROUND2_NA_CAPACITY_FRACTIONS:
        add("C3", f_na=f_na, axis="parallel_Na_K_capacity")
        add("C4", f_na=f_na, axis="cooperative_occupancy_capacity")

    thermo_ids = (
        "C5",
        "C6_112",
        "C6_123",
        "C6_213",
        "C7_112",
        "C7_123",
        "C7_213",
        "C8_112_MA",
        "C8_123_MA",
        "C8_213_MA",
    )
    for candidate_id in thermo_ids:
        for f_na in ROUND2_NA_CAPACITY_FRACTIONS:
            for equilibrium_constant in ROUND2_EQUILIBRIUM_CONSTANTS:
                add(
                    candidate_id,
                    f_na=f_na,
                    equilibrium_constant=equilibrium_constant,
                    axis="Na_K_capacity_x_equilibrium",
                )

    for candidate_id in ("C8_112_SAT", "C8_123_SAT", "C8_213_SAT"):
        for f_na in ROUND2_NA_CAPACITY_FRACTIONS:
            for equilibrium_constant in ROUND2_EQUILIBRIUM_CONSTANTS:
                for saturation_strength in ROUND2_SATURATION_STRENGTHS:
                    add(
                        candidate_id,
                        f_na=f_na,
                        equilibrium_constant=equilibrium_constant,
                        saturation_strength=saturation_strength,
                        axis="Na_K_capacity_x_equilibrium_x_saturation",
                    )

    # Acid/base bath conventions are a separate sensitivity and never replace
    # the 21 mM fixed-chassis primary context.
    assay_fraction = 1.5 / (1.5 + 1.6)
    for candidate_id in ("C2", "C3"):
        for bath in ROUND2_BATH_HCO3_MILLIMOLAR:
            add(
                candidate_id,
                f_na=assay_fraction,
                bath_hco3_mM=bath,
                axis="AE4_bath_HCO3_sensitivity",
            )
    return variants


def round2_grid(
    *,
    shared_ae4_knockout: SimulationResult | None = None,
) -> tuple[
    list[dict[str, Any]],
    list[CandidateComparison],
    dict[str, Mapping[str, SimulationResult]],
]:
    """Evaluate the predeclared interpolation grid and simulate closed rows.

    Closure-failing rows are not promoted to phenotype predictions: their WT
    state is not the calibrated state.  Round 1 already executes a full
    diagnostic trajectory for every source-family representative.  Round 2
    uses the full grid to localize the closure boundary and integrates every
    newly WT-closed interpolation with Radau.
    """

    records: list[dict[str, Any]] = []
    comparisons: list[CandidateComparison] = []
    trajectories: dict[str, Mapping[str, SimulationResult]] = {}
    shared = shared_ae4_knockout

    for index, definition in enumerate(_round2_variant_definitions(), start=1):
        f_na = float(definition["na_capacity_fraction"])
        parameters = independently_constrained_parameters(
            na_weight=f_na,
            k_weight=1.0 - f_na,
            equilibrium_constant=float(definition["equilibrium_constant"]),
            saturation_strength=float(definition["saturation_strength"]),
        )
        calibration = calibrate_candidate(
            str(definition["candidate_id"]),
            parameters=parameters,
            external_hco3_mM=float(definition["bath_hco3_mM"]),
        )
        chassis = candidate_chassis(calibration)
        raw, _ = chassis.evaluate(BASELINE_STATE, calcium=0.05, scenario=WT)
        raw_error = finite_max(raw)
        closure_ok = (
            calibration.scale_fit.max_abs_error <= 1e-8 and raw_error <= 1e-8
        )
        uses_thermodynamic_keq = str(definition["candidate_id"]).startswith(
            ("C5", "C6", "C7", "C8")
        )
        energetic_bias_supported = (
            not uses_thermodynamic_keq
            or float(definition["equilibrium_constant"]) == 1.0
        )
        coherent_fixed_bath = float(definition["bath_hco3_mM"]) == 21.0
        variant_id = (
            f"R2_{index:04d}_{definition['candidate_id']}_"
            f"fNa{f_na:.8g}_Keq{definition['equilibrium_constant']:.8g}_"
            f"sat{definition['saturation_strength']:.8g}_"
            f"B{definition['bath_hco3_mM']:.8g}"
        )
        record = {
            "variant_id": variant_id,
            **definition,
            "activity_scale": calibration.parameters.activity_scale,
            "unit_contribution": list(calibration.unit_contribution),
            "fitted_contribution": list(
                calibration.scale_fit.fitted_contribution
            ),
            "calibration_residual": list(calibration.scale_fit.residual),
            "calibration_l2_error": calibration.scale_fit.l2_error,
            "calibration_max_abs_error": calibration.scale_fit.max_abs_error,
            "wt_baseline_closure_error": raw_error,
            "wt_baseline_admissible": closure_ok,
            "assay_admissible": _assay_admissible(
                CANDIDATE_SPECS[calibration.candidate_id], parameters
            ),
            "energetic_bias_source_supported": energetic_bias_supported,
            "coherent_fixed_chassis_bath": coherent_fixed_bath,
            "fully_source_supported": bool(
                closure_ok
                and _assay_admissible(
                    CANDIDATE_SPECS[calibration.candidate_id], parameters
                )
                and energetic_bias_supported
                and coherent_fixed_bath
            ),
            "knockout_target_read_by_calibrator": False,
            "simulation_status": "not_run_round2_closure_localization",
            "ae4_ko_endpoint_flow_ratio": math.nan,
            "ae4_ko_integrated_flow_ratio": math.nan,
            "wt_endpoint_flow": math.nan,
        }
        if closure_ok:
            if shared is None:
                shared = chassis.simulate(
                    AE4_KNOCKOUT,
                    method="Radau",
                    t_eval=comparison_time_grid(),
                )
            precomputed = {
                "ae4_knockout": replay_trajectory(
                    chassis, shared, AE4_KNOCKOUT
                )
            }
            comparison, runs = compare_calibrated_candidate(
                calibration,
                round_number=2,
                mechanism_id=variant_id,
                solver_method="Radau",
                precomputed_simulations=precomputed,
            )
            comparisons.append(comparison)
            trajectories[variant_id] = runs
            record.update(
                {
                    "simulation_status": comparison.simulation_status,
                    "ae4_ko_endpoint_flow_ratio": comparison.ae4_ko_endpoint_flow_ratio,
                    "ae4_ko_integrated_flow_ratio": comparison.ae4_ko_integrated_flow_ratio,
                    "wt_endpoint_flow": comparison.wt_endpoint_flow,
                }
            )
        records.append(record)
    return records, comparisons, trajectories


def transported_k_fraction_bound(
    *, contribution_tolerance: float = 1e-8,
    wt_cycle_flux: float = HISTORICAL_WT_AE4_CYCLE_FLUX,
) -> float:
    """Exact upper bound from the fixed WT K-balance source demand."""

    return contribution_tolerance / abs(wt_cycle_flux)
