"""Nested whole-cell chassis extensions for Task 13 Round 5.

The :class:`SplitCationChassis` is deliberately a *nested* extension of the
fixed seven-state comparison chassis.  It redistributes, but does not increase,
the total Na/K-ATPase capacity and total Ca-activated K conductance between the
apical and basolateral membranes.  At zero apical fractions it reproduces the
fixed chassis algebraically.  Nonzero fractions use local extracellular K and
the local membrane voltage, add the corresponding lumen Na/K amount terms, and
re-solve both membrane-current equations.

This module does not contain or accept an AE4-null saliva target.  The released
Stage-B chloride and pH observations may be passed only through functions that
require ``stage_a_failed=True``.  The historical balance, water, and time scales
remain uncertified; no output here is converted to a physical secretion rate.

The carbon/buffer routines below are structural conservation tests, not a
calibrated salivary carbon model.  They make explicit why adding carbon and
NHE1 dynamics requires measurements of total inorganic carbon, buffer pools,
and boundary fluxes that the AE4-null experiment did not report.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares, root

from ae4_mechanism_reconstruction.chassis import (
    AE4_KNOCKOUT,
    BASELINE_STATE,
    STATE_ORDER,
    WT,
    AE4Law,
    ChassisParameters,
    FixedChassis,
    Scenario,
)


@dataclass(frozen=True)
class CationTopology:
    """Capacity fractions assigned to the apical membrane.

    The fractions are not measured for mouse submandibular acini.  A nonzero
    value is therefore a parotid-informed topology sensitivity, not a fitted
    physiological parameter.  Total capacity/conductance is conserved.
    """

    apical_pump_fraction: float = 0.0
    apical_k_fraction: float = 0.0

    def __post_init__(self) -> None:
        for name, value in asdict(self).items():
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be finite and in [0, 1]")

    @property
    def model_id(self) -> str:
        fp = self.apical_pump_fraction
        fk = self.apical_k_fraction
        if fp == 0.0 and fk == 0.0:
            return "M0_FIXED_BASOLATERAL"
        if fp > 0.0 and fk == 0.0:
            return "M1_APICAL_PUMP_ONLY"
        if fp == 0.0 and fk > 0.0:
            return "M1_APICAL_K_ONLY"
        return "M1_COUPLED_APICAL_PUMP_K"


class SplitCationChassis(FixedChassis):
    """Seven-state chassis with local apical/basal pump and K pathways.

    The retained state is ``(Na_l,K_l,height,Na_i,K_i,Cl_i,HCO3_i)``.  Thus
    this class is suitable for the first, smallest topology test (M1), not for
    claiming that the acid/base and lumen reductions have been modernized.
    """

    def __init__(
        self,
        ae4_law: AE4Law | Any | None = None,
        *,
        topology: CationTopology | None = None,
        parameters: ChassisParameters | None = None,
        **kwargs: Any,
    ) -> None:
        self.topology = topology or CationTopology()
        super().__init__(ae4_law, parameters=parameters, **kwargs)
        payload = json.dumps(
            {
                "base_fingerprint": self._parameter_fingerprint,
                "topology": asdict(self.topology),
                "equation_version": "split-cation-v1",
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        self._split_parameter_fingerprint = hashlib.sha256(
            payload.encode("utf-8")
        ).hexdigest()

    @property
    def parameter_fingerprint(self) -> str:
        return self._split_parameter_fingerprint

    def _pump_turnover(self, na_i: float, local_k_o: float) -> float:
        """Historical pump kinetics evaluated with the local external K."""

        p = self.parameters
        na_m = p.mm_to_molar * na_i
        ko_m = p.mm_to_molar * local_k_o
        return p.nak_rate * ko_m**2 * na_m**3 / (
            ko_m**2 + p.nak_alpha * na_m**3
        )

    def _evaluate(
        self,
        y: ArrayLike,
        *,
        calcium: float,
        pka_activation: float = 0.0,
        scenario: Scenario = WT,
        include_diagnostics: bool,
    ) -> tuple[NDArray[np.float64], Mapping[str, Any]]:
        # Reuse only the retained non-cation laws and algebraic reductions.  All
        # cation currents, voltages, lumen terms, and dependent checks are
        # recomputed below.  Calling the parent with diagnostics also gives one
        # consistent AE4 evaluation rather than evaluating a CTMC twice.
        _, retained = super()._evaluate(
            y,
            calcium=calcium,
            pka_activation=pka_activation,
            scenario=scenario,
            include_diagnostics=True,
        )
        state = self._state(y, require_all_positive=True)
        p = self.parameters
        na_l, k_l, height, na_i, k_i, cl_i, hco3_i = map(float, state)
        h_i = float(retained["h_i"])
        cl_l = float(retained["cl_l"])

        fp = self.topology.apical_pump_fraction
        fk = self.topology.apical_k_fraction
        j_pump_a = p.nak_capacity * fp * self._pump_turnover(na_i, k_l)
        j_pump_b = p.nak_capacity * (1.0 - fp) * self._pump_turnover(
            na_i, p.k_e
        )

        p_cl = 1.0 / (1.0 + (p.cacc_half_activation / calcium) ** p.cacc_hill)
        p_k = 1.0 / (1.0 + (p.kcc_half_activation / calcium) ** p.kcc_hill)
        g_cl = p.cacc_conductance * p_cl
        g_k_total = p.kcc_conductance * p_k
        g_k_a = fk * g_k_total
        g_k_b = (1.0 - fk) * g_k_total

        v_cl_mag = p.thermal_voltage_mv * math.log(cl_l / cl_i)
        v_k_a = p.thermal_voltage_mv * math.log(k_l / k_i)
        v_k_b = p.thermal_voltage_mv * math.log(p.k_e / k_i)
        v_t_na = p.thermal_voltage_mv * math.log(na_l / p.na_e)
        v_t_k = p.thermal_voltage_mv * math.log(k_l / p.k_e)
        g_t_na = p.tight_na_conductance
        g_t_k = p.tight_k_conductance
        g_t = g_t_na + g_t_k

        # I_a-I_t=0 and I_b+I_t=0 in particle-rate/Faraday units.  Pump
        # current is one outward positive-charge equivalent per pump cycle.
        matrix = np.asarray(
            (
                (g_cl + g_k_a + g_t, -g_t),
                (-g_t, g_k_b + g_t),
            ),
            dtype=float,
        )
        rhs = np.asarray(
            (
                -g_cl * v_cl_mag
                + g_k_a * v_k_a
                + g_t_na * v_t_na
                + g_t_k * v_t_k
                - p.faraday * j_pump_a,
                -p.faraday * j_pump_b
                + g_k_b * v_k_b
                - g_t_na * v_t_na
                - g_t_k * v_t_k,
            ),
            dtype=float,
        )
        try:
            v_a, v_b = np.linalg.solve(matrix, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError("singular split membrane-current closure") from exc
        v_t = float(v_a - v_b)
        j_cl_signed = g_cl * (v_a + v_cl_mag) / p.faraday
        j_k_a = g_k_a * (v_a - v_k_a) / p.faraday
        j_k_b = g_k_b * (v_b - v_k_b) / p.faraday
        j_t_na = g_t_na * (v_t - v_t_na) / p.faraday
        j_t_k = g_t_k * (v_t - v_t_k) / p.faraday

        q_total = float(retained["q_total"])
        volume_rate = float(retained["volume_rate"])
        j_nkcc = float(retained["j_nkcc1"])
        j_ae2 = float(retained["j_ae2"])
        j_nhe1 = float(retained["j_nhe1"])
        j_buffer = float(retained["j_buffer"])
        ae4_na = float(retained["ae4_na_source"])
        ae4_k = float(retained["ae4_k_source"])
        ae4_cl = float(retained["ae4_cl_source"])
        ae4_hco3 = float(retained["ae4_hco3_source"])

        j_pump_total = j_pump_a + j_pump_b
        amount_na = j_nkcc - 3.0 * j_pump_total + j_nhe1 + ae4_na
        amount_k = (
            j_nkcc + 2.0 * j_pump_total - j_k_a - j_k_b + ae4_k
        )
        amount_cl = 2.0 * j_nkcc + j_ae2 + j_cl_signed + ae4_cl
        amount_hco3 = j_buffer - j_ae2 + ae4_hco3

        # Complete M1 lumen Na/K balances.  Fixed lumen volume implies the
        # existing q_total outflow; Cl remains the redundant electroneutral row.
        lumen_na_source = j_t_na + 3.0 * j_pump_a
        lumen_k_source = j_t_k + j_k_a - 2.0 * j_pump_a
        raw_rhs = np.asarray(
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
        if not include_diagnostics:
            return raw_rhs, {}

        apical_current = j_cl_signed + j_k_a + j_pump_a + j_t_na + j_t_k
        basolateral_current = j_pump_b + j_k_b - j_t_na - j_t_k
        inferred_h_amount = -amount_na - amount_k + amount_cl + amount_hco3
        lumen_cl_redundancy = raw_rhs[0] + raw_rhs[1] - (
            -j_cl_signed - q_total * cl_l
        )
        checks = dict(retained["checks"])
        checks.update(
            {
                "apical_current_residual": apical_current,
                "basolateral_current_residual": basolateral_current,
                "proton_balance_identity_residual": inferred_h_amount
                - (j_buffer - j_nhe1),
                "lumen_chloride_redundancy_residual": lumen_cl_redundancy,
                "pump_capacity_partition_residual": (
                    fp + (1.0 - fp) - 1.0
                ),
                "k_conductance_partition_residual": (
                    fk + (1.0 - fk) - 1.0
                ),
            }
        )
        diagnostics = dict(retained)
        diagnostics.update(
            {
                "v_a": float(v_a),
                "v_b": float(v_b),
                "v_t": v_t,
                "j_nak": j_pump_total,
                "j_nak_apical": j_pump_a,
                "j_nak_basolateral": j_pump_b,
                "j_k": j_k_a + j_k_b,
                "j_k_apical": j_k_a,
                "j_k_basolateral": j_k_b,
                "j_cl_signed": j_cl_signed,
                "j_t_na": j_t_na,
                "j_t_k": j_t_k,
                "lumen_na_source": lumen_na_source,
                "lumen_k_source": lumen_k_source,
                "amount_rate_na": amount_na,
                "amount_rate_k": amount_k,
                "amount_rate_cl": amount_cl,
                "amount_rate_hco3": amount_hco3,
                "checks": checks,
                "cation_topology": asdict(self.topology),
                "topology_model_id": self.topology.model_id,
                "parameter_fingerprint": self.parameter_fingerprint,
                "raw_rhs": dict(zip(STATE_ORDER, map(float, raw_rhs), strict=True)),
            }
        )
        return raw_rhs, diagnostics


# The same broad positive domain and row scales used for the fixed-chassis
# comparison.  They are numerical declarations, not physiological priors.
SPLIT_STATE_LOWER = np.asarray((60.0, 0.5, 2.0, 5.0, 50.0, 1.0, 1.0))
SPLIT_STATE_UPPER = np.asarray((180.0, 30.0, 100.0, 100.0, 180.0, 100.0, 100.0))
SPLIT_RAW_SCALES = np.asarray((1e-4, 1e-5, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3))


@dataclass(frozen=True)
class SplitRoot:
    """One positive resting root of a fixed split topology."""

    success: bool
    state: tuple[float, ...]
    max_abs_raw_rhs: float
    pH_i: float
    nfev: int
    distinct_root_count: int
    failure_reason: str = ""


def _split_root_residual(
    log_state: NDArray[np.float64],
    chassis: SplitCationChassis,
    scenario: Scenario,
    calcium: float,
) -> NDArray[np.float64]:
    state = np.exp(log_state)
    try:
        raw, diagnostics = chassis.evaluate(
            state, calcium=calcium, pka_activation=0.0, scenario=scenario
        )
    except (ValueError, FloatingPointError, OverflowError):
        return np.full(7, 1e6)
    if float(diagnostics["h_i"]) <= 0.0 or not np.all(np.isfinite(raw)):
        return np.full(7, 1e6)
    return np.asarray(raw, dtype=float) / SPLIT_RAW_SCALES


def solve_split_roots(
    chassis: SplitCationChassis,
    *,
    scenario: Scenario = AE4_KNOCKOUT,
    calcium: float = 0.05,
    preferred: ArrayLike | None = None,
    start_count: int = 6,
    seed: int = 20260827,
) -> tuple[SplitRoot, ...]:
    """Find distinct positive roots with bounded log-state multistart."""

    if start_count < 1:
        raise ValueError("start_count must be positive")
    starts = [np.asarray(BASELINE_STATE, dtype=float)]
    if preferred is not None:
        row = np.asarray(preferred, dtype=float)
        if row.shape != (7,):
            raise ValueError("preferred root must contain seven states")
        if not np.allclose(row, starts[0]):
            starts.append(row)
    generator = np.random.default_rng(seed)
    while len(starts) < start_count:
        starts.append(
            np.exp(
                generator.uniform(
                    np.log(SPLIT_STATE_LOWER), np.log(SPLIT_STATE_UPPER)
                )
            )
        )
    fits: list[tuple[Any, NDArray[np.float64], Mapping[str, Any], float]] = []
    for start in starts:
        # A square hybrid solve follows the continuous root branch much more
        # accurately than a bounded least-squares solve in this ill-scaled
        # historical coordinate system.  Bounds are checked after convergence;
        # the bounded route below remains a fallback and supplies a second
        # numerical path.
        hybrid = root(
            _split_root_residual,
            np.log(np.clip(start, SPLIT_STATE_LOWER, SPLIT_STATE_UPPER)),
            args=(chassis, scenario, calcium),
            method="hybr",
            options={"xtol": 1e-10, "maxfev": 5000},
        )
        hybrid_state = np.exp(hybrid.x)
        hybrid_in_domain = bool(
            np.all(hybrid_state >= SPLIT_STATE_LOWER)
            and np.all(hybrid_state <= SPLIT_STATE_UPPER)
        )
        if hybrid_in_domain:
            try:
                hybrid_raw, hybrid_diagnostics = chassis.evaluate(
                    hybrid_state,
                    calcium=calcium,
                    pka_activation=0.0,
                    scenario=scenario,
                )
                hybrid_error = float(np.max(np.abs(hybrid_raw)))
            except (ValueError, FloatingPointError, OverflowError):
                hybrid_error = math.inf
                hybrid_diagnostics = {}
            if hybrid_error <= 1e-8 and float(
                hybrid_diagnostics.get("h_i", -1.0)
            ) > 0.0:
                fits.append(
                    (hybrid, hybrid_state, hybrid_diagnostics, hybrid_error)
                )
                continue
        fit = least_squares(
            _split_root_residual,
            np.log(np.clip(start, SPLIT_STATE_LOWER, SPLIT_STATE_UPPER)),
            args=(chassis, scenario, calcium),
            bounds=(np.log(SPLIT_STATE_LOWER), np.log(SPLIT_STATE_UPPER)),
            x_scale="jac",
            max_nfev=4000,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        state = np.exp(fit.x)
        # Refine the bounded minimum as a square root.  This prevents a small
        # scaled least-squares cost from being mislabeled as stationarity.
        refined = root(
            _split_root_residual,
            fit.x,
            args=(chassis, scenario, calcium),
            method="hybr",
            options={"xtol": 1e-10, "maxfev": 5000},
        )
        refined_state = np.exp(refined.x)
        if (
            np.all(refined_state >= SPLIT_STATE_LOWER)
            and np.all(refined_state <= SPLIT_STATE_UPPER)
        ):
            state = refined_state
            accepted_fit = refined
        else:
            accepted_fit = fit
        try:
            raw, diagnostics = chassis.evaluate(
                state, calcium=calcium, pka_activation=0.0, scenario=scenario
            )
        except (ValueError, FloatingPointError, OverflowError):
            continue
        raw_error = float(np.max(np.abs(raw)))
        if raw_error <= 1e-8 and float(diagnostics["h_i"]) > 0.0:
            fits.append((accepted_fit, state, diagnostics, raw_error))
    distinct: list[tuple[Any, NDArray[np.float64], Mapping[str, Any], float]] = []
    for row in sorted(fits, key=lambda item: item[3]):
        if not any(
            np.allclose(row[1], known[1], rtol=1e-6, atol=1e-5)
            for known in distinct
        ):
            distinct.append(row)
    return tuple(
        SplitRoot(
            success=True,
            state=tuple(map(float, state)),
            max_abs_raw_rhs=raw_error,
            pH_i=float(diagnostics["pH_i"]),
            nfev=int(getattr(fit, "nfev", 0)),
            distinct_root_count=len(distinct),
        )
        for fit, state, diagnostics, raw_error in distinct
    )


@dataclass(frozen=True)
class SplitScanRow:
    """Stage-B ionic sensitivity row; no saliva observation is present."""

    topology_model_id: str
    apical_pump_fraction: float
    apical_k_fraction: float
    fraction_provenance: str
    root_success: bool
    root_id: str
    max_abs_raw_rhs: float
    na_l_mM: float
    k_l_mM: float
    height_um: float
    na_i_mM: float
    k_i_mM: float
    cl_i_mM: float
    hco3_i_mM: float
    pH_i: float
    cl_target_residual_mM: float
    pH_target_residual: float
    normalized_ionic_distance: float
    distinct_root_count: int
    target_stage: str
    secretion_target_read: bool
    interpretation: str


def scan_split_topologies(
    ae4_law: Any,
    *,
    stage_a_failed: bool,
    cl_target_mM: float,
    pH_target: float,
    pump_fractions: Sequence[float],
    k_fractions: Sequence[float],
    cl_scale_mM: float = 1.6,
    pH_scale: float = 0.02,
    parameters: ChassisParameters | None = None,
    preferred: ArrayLike | None = None,
    start_count: int = 4,
) -> tuple[SplitScanRow, ...]:
    """Re-solve the AE4-null root across predeclared topology sensitivities.

    Fractions are never optimized here.  The caller must predeclare the grid,
    and every row is retained.  ``cl_scale_mM`` and ``pH_scale`` only normalize
    the released ionic residuals; they are not likelihood assumptions.
    """

    if not stage_a_failed:
        raise PermissionError("Stage-B ionic targets remain sealed until Stage A fails")
    if cl_target_mM <= 0.0 or not 0.0 < pH_target < 14.0:
        raise ValueError("released ionic targets must be physical")
    if cl_scale_mM <= 0.0 or pH_scale <= 0.0:
        raise ValueError("target scales must be positive")
    rows: list[SplitScanRow] = []
    continuation = None if preferred is None else np.asarray(preferred, dtype=float)
    for fp in pump_fractions:
        for fk in k_fractions:
            topology = CationTopology(float(fp), float(fk))
            chassis = SplitCationChassis(
                ae4_law,
                topology=topology,
                parameters=parameters,
            )
            roots = solve_split_roots(
                chassis,
                scenario=AE4_KNOCKOUT,
                preferred=continuation,
                start_count=start_count,
            )
            if not roots:
                rows.append(
                    SplitScanRow(
                        topology_model_id=topology.model_id,
                        apical_pump_fraction=float(fp),
                        apical_k_fraction=float(fk),
                        fraction_provenance=(
                            "predeclared sensitivity; uncalibrated for mouse "
                            "submandibular acini"
                        ),
                        root_success=False,
                        root_id="",
                        max_abs_raw_rhs=math.nan,
                        na_l_mM=math.nan,
                        k_l_mM=math.nan,
                        height_um=math.nan,
                        na_i_mM=math.nan,
                        k_i_mM=math.nan,
                        cl_i_mM=math.nan,
                        hco3_i_mM=math.nan,
                        pH_i=math.nan,
                        cl_target_residual_mM=math.nan,
                        pH_target_residual=math.nan,
                        normalized_ionic_distance=math.inf,
                        distinct_root_count=0,
                        target_stage="LOC-B",
                        secretion_target_read=False,
                        interpretation="no positive root found in declared domain",
                    )
                )
                continue
            for index, root in enumerate(roots):
                state = np.asarray(root.state)
                d_cl = float(state[5] - cl_target_mM)
                d_ph = float(root.pH_i - pH_target)
                distance = math.hypot(d_cl / cl_scale_mM, d_ph / pH_scale)
                rows.append(
                    SplitScanRow(
                        topology_model_id=topology.model_id,
                        apical_pump_fraction=float(fp),
                        apical_k_fraction=float(fk),
                        fraction_provenance=(
                            "predeclared sensitivity; uncalibrated for mouse "
                            "submandibular acini"
                        ),
                        root_success=True,
                        root_id=f"fp{fp:.6g}_fk{fk:.6g}_r{index}",
                        max_abs_raw_rhs=root.max_abs_raw_rhs,
                        na_l_mM=float(state[0]),
                        k_l_mM=float(state[1]),
                        height_um=float(state[2]),
                        na_i_mM=float(state[3]),
                        k_i_mM=float(state[4]),
                        cl_i_mM=float(state[5]),
                        hco3_i_mM=float(state[6]),
                        pH_i=root.pH_i,
                        cl_target_residual_mM=d_cl,
                        pH_target_residual=d_ph,
                        normalized_ionic_distance=distance,
                        distinct_root_count=root.distinct_root_count,
                        target_stage="LOC-B",
                        secretion_target_read=False,
                        interpretation=(
                            "ionic localization sensitivity only; neither slice "
                            "fraction is a gland-matched measurement"
                        ),
                    )
                )
            continuation = np.asarray(roots[0].state)
    return tuple(rows)


CARBON_COORDINATES = (
    "na_i",
    "hco3_i",
    "h_i",
    "co2_i",
    "buffer_h",
    "buffer_minus",
)

# Exact per-event signatures.  NHE1 is positive for Na entry/H exit; hydration
# and buffer dissociation have their chemical forward directions.
ACID_BASE_SIGNATURES: Mapping[str, tuple[int, ...]] = {
    "nhe1": (1, 0, -1, 0, 0, 0),
    "carbon_hydration": (0, 1, 1, -1, 0, 0),
    "buffer_dissociation": (0, 0, 1, 0, -1, 1),
}


@dataclass(frozen=True)
class CarbonBufferRates:
    """Rates and exact conservation residuals for a closed reaction block."""

    hydration_rate: float
    buffer_dissociation_rate: float
    sources: Mapping[str, float]
    inorganic_carbon_residual: float
    buffer_site_residual: float
    charge_residual: float


def conserved_carbon_buffer_rates(
    *,
    co2_i: float,
    hco3_i: float,
    h_i: float,
    buffer_h: float,
    buffer_minus: float,
    hydration_forward: float,
    hydration_reverse: float,
    buffer_dissociation: float,
    buffer_association: float,
) -> CarbonBufferRates:
    """Return mass-action sources for conserved carbon and buffer reactions."""

    values = (
        co2_i,
        hco3_i,
        h_i,
        buffer_h,
        buffer_minus,
        hydration_forward,
        hydration_reverse,
        buffer_dissociation,
        buffer_association,
    )
    if any(not math.isfinite(value) or value < 0.0 for value in values):
        raise ValueError("states and rate constants must be finite and nonnegative")
    hydration = hydration_forward * co2_i - hydration_reverse * h_i * hco3_i
    dissociation = (
        buffer_dissociation * buffer_h
        - buffer_association * h_i * buffer_minus
    )
    sources = {
        "na_i": 0.0,
        "hco3_i": hydration,
        "h_i": hydration + dissociation,
        "co2_i": -hydration,
        "buffer_h": -dissociation,
        "buffer_minus": dissociation,
    }
    carbon_residual = sources["hco3_i"] + sources["co2_i"]
    buffer_residual = sources["buffer_h"] + sources["buffer_minus"]
    charge_residual = (
        sources["h_i"]
        - sources["hco3_i"]
        - sources["buffer_minus"]
    )
    return CarbonBufferRates(
        hydration_rate=hydration,
        buffer_dissociation_rate=dissociation,
        sources=sources,
        inorganic_carbon_residual=carbon_residual,
        buffer_site_residual=buffer_residual,
        charge_residual=charge_residual,
    )


def _exact_rank(matrix: Sequence[Sequence[int | float]]) -> int:
    rows = [[Fraction(str(value)) for value in row] for row in matrix]
    if not rows:
        return 0
    rank = 0
    width = len(rows[0])
    for column in range(width):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[row], rows[rank], strict=True)
            ]
        rank += 1
    return rank


def acid_base_identifiability_audit() -> Mapping[str, int | str]:
    """Exact ranks showing what released Cl/pH cannot distinguish.

    Cl is absent from this reaction block and pH supplies one scalar H row.
    Consequently, the three independent reaction signatures have rank one in
    the released-observation projection and leave a two-dimensional flux
    nullspace even before unknown carbon/buffer totals and boundary fluxes are
    counted.
    """

    matrix = np.column_stack(tuple(ACID_BASE_SIGNATURES.values())).tolist()
    full_rank = _exact_rank(matrix)
    # Observation rows are (Cl_i,H_i).  Every module has zero direct Cl source.
    observed = [
        [0, 0, 0],
        [
            ACID_BASE_SIGNATURES[name][2]
            for name in ACID_BASE_SIGNATURES
        ],
    ]
    observation_rank = _exact_rank(observed)
    return {
        "full_reaction_signature_rank": full_rank,
        "released_cl_ph_direct_rank": observation_rank,
        "reaction_flux_nullity_under_cl_ph": len(ACID_BASE_SIGNATURES)
        - observation_rank,
        "split_topology_direct_acid_base_rank": 0,
        "conclusion": (
            "released Cl/pH cannot distinguish indirect cation-current feedback "
            "from NHE1/carbon/buffer flux; total carbon, buffer pool, or a "
            "matched flux measurement is required"
        ),
    }


__all__ = [
    "ACID_BASE_SIGNATURES",
    "CARBON_COORDINATES",
    "CarbonBufferRates",
    "CationTopology",
    "SplitCationChassis",
    "SplitRoot",
    "SplitScanRow",
    "acid_base_identifiability_audit",
    "conserved_carbon_buffer_rates",
    "scan_split_topologies",
    "solve_split_roots",
]
