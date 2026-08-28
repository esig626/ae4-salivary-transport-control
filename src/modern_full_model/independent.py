"""Independent Task 13B equation transcription and numerical checks.

This module intentionally does **not** call the production whole-cell object
or any of its ``evaluate``/``rhs`` methods.  The conserved-coordinate, acid/base,
homeostasis, AE4, membrane-current, water, and outflow equations are
transcribed a second time so that agreement is capable of detecting an error
in the production assembly.  Only immutable parameter containers and a
stimulus callable are shared.

The transcription uses the same declared physical units as the production
model: amounts are fmol, volumes are pL, concentrations are mM, currents are
A, and time is seconds.  Its charge-manifold root refinement removes the two
dependent chloride coordinates and the two corresponding redundant source
rows instead of asking an unconstrained least-squares problem to discover a
conserved charge leaf.

No held-out target, genotype ratio, or phenotype value is present here.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any, Callable, Mapping, Protocol, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

from .parameters import FullModelParameters


FMOL_PER_MOL = 1.0e15
CORE_SIZE = 12
FREE_CORE_INDICES: tuple[int, ...] = (0, 1, 3, 4, 5, 6, 7, 9, 10, 11)
INDEPENDENT_CORE_ROWS: tuple[int, ...] = FREE_CORE_INDICES


class IndependentStimulus(Protocol):
    def __call__(self, time_s: float) -> Any: ...


@dataclass(frozen=True)
class IndependentConstantStimulus:
    calcium_uM: float
    beta_input: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.calcium_uM) or self.calcium_uM < 0.0:
            raise ValueError("independent calcium input must be nonnegative")
        if not math.isfinite(self.beta_input) or not 0.0 <= self.beta_input <= 1.0:
            raise ValueError("independent beta input must lie in [0, 1]")

    def __call__(self, time_s: float) -> Mapping[str, float]:
        if not math.isfinite(time_s):
            raise ValueError("independent stimulus time must be finite")
        return {"calcium_uM": self.calcium_uM, "beta_input": self.beta_input}


@dataclass(frozen=True)
class IndependentSecretagogueProtocol:
    """Independent transcription of the frozen combined CCh/IPR input."""

    onset_s: float = 0.0
    duration_s: float = 600.0
    resting_calcium_uM: float = 0.05
    stimulated_calcium_uM: float = 0.55
    beta_occupancy_on: float = 1.0

    def __post_init__(self) -> None:
        values = (
            self.onset_s,
            self.duration_s,
            self.resting_calcium_uM,
            self.stimulated_calcium_uM,
            self.beta_occupancy_on,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("independent protocol values must be finite")
        if self.onset_s < 0.0 or self.duration_s <= 0.0:
            raise ValueError("independent protocol onset/duration are invalid")
        if min(self.resting_calcium_uM, self.stimulated_calcium_uM) < 0.0:
            raise ValueError("independent calcium inputs must be nonnegative")
        if not 0.0 <= self.beta_occupancy_on <= 1.0:
            raise ValueError("independent beta occupancy must lie in [0, 1]")

    @property
    def end_s(self) -> float:
        return self.onset_s + self.duration_s

    def __call__(self, time_s: float) -> Mapping[str, float]:
        if not math.isfinite(time_s):
            raise ValueError("independent protocol time must be finite")
        active = self.onset_s < time_s <= self.end_s
        return {
            "calcium_uM": (
                self.stimulated_calcium_uM if active else self.resting_calcium_uM
            ),
            "beta_input": self.beta_occupancy_on if active else 0.0,
        }


def independent_physical_time_grid(
    duration_s: float = 600.0, sample_step_s: float = 5.0
) -> NDArray[np.float64]:
    """Seconds grid with an exact basal row and positive onset right-limit."""

    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise ValueError("independent duration must be positive")
    if not math.isfinite(sample_step_s) or sample_step_s <= 0.0:
        raise ValueError("independent sample step must be positive")
    regular = np.arange(sample_step_s, duration_s, sample_step_s, dtype=float)
    return np.concatenate((np.asarray((0.0, 1.0e-6)), regular, np.asarray((duration_s,))))


@dataclass(frozen=True)
class IndependentAE4Kinetics:
    """Plain AE4 parameter record, independent of the production evaluator."""

    carrier_amount_fmol: float
    common_cl_attempt_rate_s: float
    na_loaded_attempt_rate_s: float
    k_loaded_attempt_rate_s: float
    cooperative_gate: Any | None = None


@dataclass(frozen=True)
class IndependentCalibrationDefinition:
    """Inputs needed to reconstruct one WT calibration outside its builder."""

    mixed_bath_na_attempt_fraction: float
    apical_pump_fraction: float
    apical_k_fraction: float
    ae4_carrier_amount_fmol_at_unit_rate_gauge: float
    nhe1_capacity_fmol_s: float
    common_membrane_conductance_scale: float
    cell_other_impermeant_osmoles_fmol: float
    resting_calcium_uM: float = 0.058

    def __post_init__(self) -> None:
        fractions = (
            self.mixed_bath_na_attempt_fraction,
            self.apical_pump_fraction,
            self.apical_k_fraction,
        )
        if not all(math.isfinite(value) and 0.0 < value < 1.0 for value in fractions):
            raise ValueError("independent calibration fractions must lie in (0, 1)")
        positives = (
            self.ae4_carrier_amount_fmol_at_unit_rate_gauge,
            self.nhe1_capacity_fmol_s,
            self.common_membrane_conductance_scale,
            self.cell_other_impermeant_osmoles_fmol,
            self.resting_calcium_uM,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in positives):
            raise ValueError("independent calibration values must be positive")

    def reconstruct_parameters(
        self, base: FullModelParameters | None = None
    ) -> tuple[FullModelParameters, IndependentAE4Kinetics]:
        """Apply the frozen calibration without importing its build helper."""

        source = FullModelParameters() if base is None else base
        scale = self.common_membrane_conductance_scale
        parameters = replace(
            source,
            homeostasis=replace(
                source.homeostasis,
                nhe1_capacity_fmol_s=self.nhe1_capacity_fmol_s,
            ),
            membranes=replace(
                source.membranes,
                apical_pump_fraction=self.apical_pump_fraction,
                apical_k_fraction=self.apical_k_fraction,
                g_k_total_S=source.membranes.g_k_total_S * scale,
                g_cl_apical_S=source.membranes.g_cl_apical_S * scale,
                g_basolateral_background_S=(
                    source.membranes.g_basolateral_background_S * scale
                ),
            ),
            geometry=replace(
                source.geometry,
                cell_impermeant_osmoles_fmol=(
                    self.cell_other_impermeant_osmoles_fmol
                ),
            ),
        )
        fraction = self.mixed_bath_na_attempt_fraction
        ae4 = IndependentAE4Kinetics(
            carrier_amount_fmol=(
                self.ae4_carrier_amount_fmol_at_unit_rate_gauge
            ),
            common_cl_attempt_rate_s=1.0,
            na_loaded_attempt_rate_s=2.0 * fraction,
            k_loaded_attempt_rate_s=2.0 * (1.0 - fraction),
        )
        return parameters, ae4


@dataclass(frozen=True)
class IndependentRegulation:
    """Numerically explicit regulatory system used by the transcription.

    ``rhs_function`` and ``output_function`` must be independently written
    formulas.  They must not delegate to the production regulatory object.
    """

    family: str
    state_names: tuple[str, ...]
    basal_state: tuple[float, ...]
    rhs_function: Callable[[float, Sequence[float], float], Sequence[float]]
    output_function: Callable[[float, Sequence[float], float], float]

    def __post_init__(self) -> None:
        if len(self.state_names) != len(self.basal_state):
            raise ValueError("regulatory names and basal state must have equal length")
        if len(set(self.state_names)) != len(self.state_names):
            raise ValueError("regulatory state names must be unique")
        if not all(math.isfinite(value) for value in self.basal_state):
            raise ValueError("regulatory basal state must be finite")

    def rhs(self, time_s: float, state: Sequence[float], beta_input: float) -> tuple[float, ...]:
        values = tuple(float(value) for value in self.rhs_function(time_s, state, beta_input))
        if len(values) != len(self.state_names) or not all(math.isfinite(value) for value in values):
            raise ValueError("independent regulatory RHS has invalid shape or value")
        return values

    def capacity_multiplier(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> float:
        value = float(self.output_function(time_s, state, beta_input))
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("independent regulatory capacity must be finite and nonnegative")
        return value


@dataclass(frozen=True)
class IndependentN2NkccRegulation:
    """Independent one-state CCh/Ca-to-NKCC1 capacity hypothesis.

    The calcium arm is protocol-normalized rather than interpreted as a
    measured dose response.  Only the capacity of the already reversible
    1 Na : 1 K : 2 Cl law is multiplied; its affinity and stoichiometry are
    unchanged.
    """

    fully_activated_multiplier: float
    tau_activation_s: float
    resting_calcium_uM: float = 0.058
    stimulated_calcium_uM: float = 0.55
    family: str = "N2_CCH_EFFECTIVE_DYNAMIC"
    state_names: tuple[str, ...] = ("nkcc1_effective_activation_fraction",)
    basal_state: tuple[float, ...] = (0.0,)

    def __post_init__(self) -> None:
        if (
            not math.isfinite(self.fully_activated_multiplier)
            or self.fully_activated_multiplier < 1.0
        ):
            raise ValueError("N2 fully activated multiplier must be finite and >= 1")
        if not math.isfinite(self.tau_activation_s) or self.tau_activation_s <= 0.0:
            raise ValueError("N2 activation time must be finite and positive")
        if (
            not math.isfinite(self.resting_calcium_uM)
            or not math.isfinite(self.stimulated_calcium_uM)
            or self.resting_calcium_uM < 0.0
            or self.stimulated_calcium_uM <= self.resting_calcium_uM
        ):
            raise ValueError("N2 calcium endpoints must be finite and ordered")

    def normalized_calcium_arm(self, calcium_uM: float) -> float:
        calcium = float(calcium_uM)
        if not math.isfinite(calcium) or calcium < 0.0:
            raise ValueError("N2 calcium input must be finite and nonnegative")
        raw = (calcium - self.resting_calcium_uM) / (
            self.stimulated_calcium_uM - self.resting_calcium_uM
        )
        return min(1.0, max(0.0, raw))

    @staticmethod
    def _checked_activation(state: Sequence[float]) -> float:
        if len(tuple(state)) != 1:
            raise ValueError("N2 requires one activation fraction")
        activation = float(tuple(state)[0])
        if not math.isfinite(activation) or not -1.0e-4 <= activation <= 1.0 + 1.0e-4:
            raise ValueError("N2 activation fraction left its integration trial band")
        return min(1.0, max(0.0, activation))

    def rhs(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> tuple[float, ...]:
        if not math.isfinite(time_s) or time_s < 0.0:
            raise ValueError("N2 time must be finite and nonnegative")
        activation = self._checked_activation(state)
        target = self.normalized_calcium_arm(calcium_uM)
        return ((target - activation) / self.tau_activation_s,)

    def activation_fraction(self, state: Sequence[float]) -> float:
        return self._checked_activation(state)

    def capacity_multiplier(self, state: Sequence[float]) -> float:
        activation = self._checked_activation(state)
        return 1.0 + (self.fully_activated_multiplier - 1.0) * activation


@dataclass(frozen=True)
class IndependentN1NkccRegulation:
    """Independent algebraic CCh/Ca-to-NKCC1 capacity hypothesis."""

    fully_activated_multiplier: float
    resting_calcium_uM: float = 0.058
    stimulated_calcium_uM: float = 0.55
    family: str = "N1_CCH_ALGEBRAIC"
    state_names: tuple[str, ...] = ()
    basal_state: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        if (
            not math.isfinite(self.fully_activated_multiplier)
            or self.fully_activated_multiplier < 1.0
        ):
            raise ValueError("N1 fully activated multiplier must be finite and >= 1")
        if (
            not math.isfinite(self.resting_calcium_uM)
            or not math.isfinite(self.stimulated_calcium_uM)
            or self.resting_calcium_uM < 0.0
            or self.stimulated_calcium_uM <= self.resting_calcium_uM
        ):
            raise ValueError("N1 calcium endpoints must be finite and ordered")

    def normalized_calcium_arm(self, calcium_uM: float) -> float:
        calcium = float(calcium_uM)
        if not math.isfinite(calcium) or calcium < 0.0:
            raise ValueError("N1 calcium input must be finite and nonnegative")
        raw = (calcium - self.resting_calcium_uM) / (
            self.stimulated_calcium_uM - self.resting_calcium_uM
        )
        return min(1.0, max(0.0, raw))

    @staticmethod
    def _check_empty_state(state: Sequence[float]) -> None:
        if tuple(state):
            raise ValueError("N1 has no dynamic state")

    def rhs(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> tuple[float, ...]:
        if not math.isfinite(time_s) or time_s < 0.0:
            raise ValueError("N1 time must be finite and nonnegative")
        self._check_empty_state(state)
        self.normalized_calcium_arm(calcium_uM)
        return ()

    def activation_fraction(self, state: Sequence[float]) -> float:
        self._check_empty_state(state)
        # The actual algebraic activation is evaluated from calcium in the
        # capacity method; this state-only accessor is unused by N1.
        return 0.0

    def capacity_multiplier(
        self, state: Sequence[float], calcium_uM: float | None = None
    ) -> float:
        self._check_empty_state(state)
        if calcium_uM is None:
            raise ValueError("N1 capacity evaluation requires calcium")
        activation = self.normalized_calcium_arm(calcium_uM)
        return 1.0 + (self.fully_activated_multiplier - 1.0) * activation


def independent_r0_regulation(
    *,
    basal_multiplier: float,
    fully_activated_increment: float,
    coupling_scale: float = 1.0,
    construct_response_scale: float = 1.0,
) -> IndependentRegulation:
    """Zero-state instantaneous comparison map (not a dynamic model)."""

    if min(
        basal_multiplier,
        fully_activated_increment,
        coupling_scale,
        construct_response_scale,
    ) < 0.0:
        raise ValueError("R0 capacity multipliers must be nonnegative")
    return IndependentRegulation(
        family="R0",
        state_names=(),
        basal_state=(),
        rhs_function=lambda _t, _y, _beta: (),
        output_function=lambda _t, _y, beta: basal_multiplier
        + fully_activated_increment
        * coupling_scale
        * construct_response_scale
        * float(beta),
    )


def independent_r1_regulation(
    *,
    tau_activation_s: float,
    basal_multiplier: float,
    fully_activated_increment: float,
    coupling_scale: float = 1.0,
    construct_response_scale: float = 1.0,
) -> IndependentRegulation:
    """One effective beta/PKA/AE4 activation state.

    This constructor is an independent formula utility, not a claim that R1
    is the retained Task 13B regulatory family.
    """

    if tau_activation_s <= 0.0:
        raise ValueError("R1 activation time must be positive")
    if min(
        basal_multiplier,
        fully_activated_increment,
        coupling_scale,
        construct_response_scale,
    ) < 0.0:
        raise ValueError("R1 capacity multipliers must be nonnegative")

    return IndependentRegulation(
        family="R1",
        state_names=("ae4_effective_activation_fraction",),
        basal_state=(0.0,),
        rhs_function=lambda _t, y, beta: (
            (float(beta) - float(y[0])) / tau_activation_s,
        ),
        output_function=lambda _t, y, _beta: basal_multiplier
        + fully_activated_increment
        * coupling_scale
        * construct_response_scale
        * float(y[0]),
    )


def independent_r2_regulation(
    *,
    tau_camp_s: float,
    tau_activation_s: float,
    basal_multiplier: float,
    fully_activated_increment: float,
    coupling_scale: float = 1.0,
    construct_response_scale: float = 1.0,
) -> IndependentRegulation:
    """Effective cAMP state followed by an AE4 activation state."""

    values = (tau_camp_s, tau_activation_s)
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("R2 time and activation scales must be positive")
    if min(
        basal_multiplier,
        fully_activated_increment,
        coupling_scale,
        construct_response_scale,
    ) < 0.0:
        raise ValueError("R2 capacity multipliers must be nonnegative")

    def rhs(_t: float, y: Sequence[float], beta: float) -> tuple[float, float]:
        camp, active = (float(value) for value in y)
        return (
            (float(beta) - camp) / tau_camp_s,
            (camp - active) / tau_activation_s,
        )

    return IndependentRegulation(
        family="R2",
        state_names=("camp_effective_fraction", "ae4_activation_fraction"),
        basal_state=(0.0, 0.0),
        rhs_function=rhs,
        output_function=lambda _t, y, _beta: basal_multiplier
        + fully_activated_increment
        * coupling_scale
        * construct_response_scale
        * float(y[1]),
    )


def independent_r3_regulation(
    *,
    tau_pka_s: float,
    forward_regulation_rate_s: float,
    reverse_regulation_rate_s: float,
    basal_multiplier: float,
    fully_activated_increment: float,
    coupling_scale: float = 1.0,
    construct_response_scale: float = 1.0,
) -> IndependentRegulation:
    """Effective PKA state plus a phosphorylation-like AE4 fraction.

    The second state is not labelled a demonstrated phosphorylated fraction:
    the primary evidence establishes pathway and S173 dependence but does not
    measure an elementary phosphorylation rate.
    """

    rates = (tau_pka_s, forward_regulation_rate_s, reverse_regulation_rate_s)
    if not all(math.isfinite(value) and value > 0.0 for value in rates):
        raise ValueError("R3 time and forward/reverse rates must be positive")
    if min(
        basal_multiplier,
        fully_activated_increment,
        coupling_scale,
        construct_response_scale,
    ) < 0.0:
        raise ValueError("R3 capacity multipliers must be nonnegative")

    def rhs(_t: float, y: Sequence[float], beta: float) -> tuple[float, float]:
        pka, regulated = (float(value) for value in y)
        return (
            (float(beta) - pka) / tau_pka_s,
            forward_regulation_rate_s * pka * (1.0 - regulated)
            - reverse_regulation_rate_s * (1.0 - pka) * regulated,
        )

    return IndependentRegulation(
        family="R3",
        state_names=(
            "pka_effective_activity_fraction",
            "ae4_regulated_fraction",
        ),
        basal_state=(0.0, 0.0),
        rhs_function=rhs,
        output_function=lambda _t, y, _beta: basal_multiplier
        + fully_activated_increment
        * coupling_scale
        * construct_response_scale
        * float(y[1]),
    )


@dataclass(frozen=True)
class IndependentSpeciation:
    ph: float
    h_mM: float
    co2_mM: float
    hco3_mM: float
    co3_mM: float
    total_alkalinity_mM: float
    alkalinity_residual_mM: float


@dataclass(frozen=True)
class IndependentEvaluation:
    rhs: NDArray[np.float64]
    cell_concentrations_mM: Mapping[str, float]
    lumen_concentrations_mM: Mapping[str, float]
    cell_speciation: IndependentSpeciation
    lumen_speciation: IndependentSpeciation
    capacity_multiplier: float
    nkcc_capacity_multiplier: float
    nkcc_activation_fraction: float
    flow_pL_s: float
    voltages_V: Mapping[str, float]
    current_residuals_A: Mapping[str, float]
    fluxes_fmol_s: Mapping[str, float]
    conservation_residuals: Mapping[str, float]


@dataclass(frozen=True)
class IndependentRootResult:
    state: NDArray[np.float64]
    success: bool
    max_abs_scaled_residual: float
    max_abs_amount_source_fmol_s: float
    max_abs_volume_source_pL_s: float
    max_abs_regulatory_rhs_s_inv: float
    nfev: int
    cost: float
    scaled_jacobian_singular_values: tuple[float, ...]
    jacobian_rank: int
    jacobian_nullity: int
    message: str


@dataclass(frozen=True)
class IndependentTrajectory:
    method: str
    onset_handling: str
    success: bool
    message: str
    time_s: NDArray[np.float64]
    states: NDArray[np.float64]
    flow_pL_s: NDArray[np.float64]
    cumulative_flow_pL: NDArray[np.float64]
    cell_cl_mM: NDArray[np.float64]
    cell_ph: NDArray[np.float64]
    capacity_multiplier: NDArray[np.float64]
    nkcc_capacity_multiplier: NDArray[np.float64]
    nkcc_activation_fraction: NDArray[np.float64]
    max_abs_conservation_fmol_s: float
    max_abs_water_accounting_pL_s: float


def _positive(*values: float) -> None:
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("independent thermodynamic inputs must be finite and positive")


def _carbon_fractions(ph: float, p: FullModelParameters) -> tuple[float, float, float]:
    h_M = 10.0 ** (-ph)
    k1 = 10.0 ** (-p.acid_base.carbon_pka1)
    k2 = 10.0 ** (-p.acid_base.carbon_pka2)
    denominator = h_M * h_M + k1 * h_M + k1 * k2
    return h_M * h_M / denominator, k1 * h_M / denominator, k1 * k2 / denominator


def independent_total_alkalinity_mM(
    ph: float,
    total_carbon_mM: float,
    buffer_total_mM: float,
    buffer_pka: float,
    parameters: FullModelParameters,
) -> float:
    _, alpha_hco3, alpha_co3 = _carbon_fractions(ph, parameters)
    h_mM = 1.0e3 * 10.0 ** (-ph)
    oh_mM = 1.0e3 * 10.0 ** (ph - parameters.acid_base.water_pkw)
    buffer_minus = buffer_total_mM / (1.0 + 10.0 ** (buffer_pka - ph))
    return (
        total_carbon_mM * alpha_hco3
        + 2.0 * total_carbon_mM * alpha_co3
        + buffer_minus
        + oh_mM
        - h_mM
    )


def independent_speciate(
    total_carbon_mM: float,
    total_alkalinity_mM: float,
    buffer_total_mM: float,
    buffer_pka: float,
    parameters: FullModelParameters,
) -> IndependentSpeciation:
    _positive(total_carbon_mM)
    if buffer_total_mM < 0.0:
        raise ValueError("buffer total cannot be negative")

    def residual(ph: float) -> float:
        return independent_total_alkalinity_mM(
            ph, total_carbon_mM, buffer_total_mM, buffer_pka, parameters
        ) - total_alkalinity_mM

    lower = parameters.acid_base.ph_lower
    upper = parameters.acid_base.ph_upper
    if residual(lower) * residual(upper) > 0.0:
        raise ValueError("alkalinity is outside the independent speciation bracket")
    ph = float(brentq(residual, lower, upper, xtol=5.0e-13, rtol=5.0e-15))
    alpha_co2, alpha_hco3, alpha_co3 = _carbon_fractions(ph, parameters)
    h_mM = 1.0e3 * 10.0 ** (-ph)
    co2 = total_carbon_mM * alpha_co2
    hco3 = total_carbon_mM * alpha_hco3
    co3 = total_carbon_mM * alpha_co3
    return IndependentSpeciation(
        ph=ph,
        h_mM=h_mM,
        co2_mM=co2,
        hco3_mM=hco3,
        co3_mM=co3,
        total_alkalinity_mM=total_alkalinity_mM,
        alkalinity_residual_mM=residual(ph),
    )


def _nernst(c_first: float, c_second: float, z: int, thermal_voltage: float) -> float:
    _positive(c_first, c_second, thermal_voltage)
    return thermal_voltage / z * math.log(c_second / c_first)


def _current_to_flux(current_A: float, z: int, faraday: float) -> float:
    return current_A / (z * faraday) * FMOL_PER_MOL


def _pump(capacity: float, na_i: float, k_out: float, na_half: float, k_half: float) -> float:
    na_gate = na_i**3 / (na_i**3 + na_half**3)
    k_gate = k_out**2 / (k_out**2 + k_half**2)
    return capacity * na_gate * k_gate


def _ae4_sources(
    ci: Mapping[str, float],
    bath: Mapping[str, float],
    ae4_parameters: Any | None,
    capacity_multiplier: float,
    expression_scale: float,
) -> Mapping[str, float]:
    zero = {
        "na": 0.0,
        "k": 0.0,
        "cl": 0.0,
        "tic": 0.0,
        "alk": 0.0,
        "affinity_na": 0.0,
        "affinity_k": 0.0,
        "charge": 0.0,
        "entropy_over_r": 0.0,
    }
    if ae4_parameters is None or expression_scale == 0.0:
        return zero
    p = ae4_parameters
    gate_na = 1.0
    gate_k = 1.0
    gate = getattr(p, "cooperative_gate", None)
    if gate is not None:
        def empirical_factor(cation: float, ec50: float, hill: float, floor: float) -> float:
            occupancy = 1.0 / (1.0 + math.exp(hill * math.log(ec50 / cation)))
            return floor + (1.0 - floor) * occupancy

        gate_na = empirical_factor(bath["na"], gate.ec50_na_mM, gate.hill_na, gate.floor_na)
        gate_k = empirical_factor(bath["k"], gate.ec50_k_mM, gate.hill_k, gate.floor_k)

    log_cl = math.log(bath["cl"] / ci["cl"])
    log_bic = 2.0 * math.log(ci["hco3"] / bath["hco3"])
    log_na = math.log(ci["na"] / bath["na"]) + log_bic
    log_k = math.log(ci["k"] / bath["k"]) + log_bic

    def pair(attempt: float, log_ratio: float) -> tuple[float, float]:
        return attempt * math.exp(0.5 * log_ratio), attempt * math.exp(-0.5 * log_ratio)

    cl_f, cl_r = pair(p.common_cl_attempt_rate_s, log_cl)
    na_f, na_r = pair(p.na_loaded_attempt_rate_s * gate_na, log_na)
    k_f, k_r = pair(p.k_loaded_attempt_rate_s * gate_k, log_k)
    outward_to_inward = cl_f + na_r + k_r
    inward_to_outward = cl_r + na_f + k_f
    sigma = outward_to_inward + inward_to_outward
    occ_out = inward_to_outward / sigma
    occ_in = outward_to_inward / sigma
    carrier = p.carrier_amount_fmol * capacity_multiplier * expression_scale
    j_na = carrier * (na_f * occ_in - na_r * occ_out)
    j_k = carrier * (k_f * occ_in - k_r * occ_out)
    total = j_na + j_k
    affinity_na = log_cl + log_na
    affinity_k = log_cl + log_k
    sources = {
        "na": -j_na,
        "k": -j_k,
        "cl": total,
        "tic": -2.0 * total,
        "alk": -2.0 * total,
        "affinity_na": affinity_na,
        "affinity_k": affinity_k,
        "charge": -j_na - j_k - total + 2.0 * total,
        "entropy_over_r": j_na * affinity_na + j_k * affinity_k,
    }
    return sources


def state_from_free_charge_coordinates(
    free: ArrayLike,
    *,
    fixed_cell_anion_equivalents_fmol: float,
    regulation: Sequence[float] = (),
) -> NDArray[np.float64]:
    """Recover both chloride amounts from exact bulk electroneutrality."""

    x = np.asarray(free, dtype=float)
    if x.shape != (len(FREE_CORE_INDICES),):
        raise ValueError("free core vector must have ten coordinates")
    y = np.empty(CORE_SIZE + len(regulation), dtype=float)
    y[list(FREE_CORE_INDICES)] = x
    y[2] = y[0] + y[1] - y[4] - fixed_cell_anion_equivalents_fmol
    y[8] = y[6] + y[7] - y[10]
    if regulation:
        y[CORE_SIZE:] = np.asarray(regulation, dtype=float)
    return y


def free_charge_coordinates(state: ArrayLike) -> NDArray[np.float64]:
    y = np.asarray(state, dtype=float)
    if y.ndim != 1 or y.size < CORE_SIZE:
        raise ValueError("state is smaller than the independent core")
    return np.asarray(y[list(FREE_CORE_INDICES)], dtype=float)


def independent_nearby_charge_preserving_states(
    state: ArrayLike, *, fraction: float = 1.0e-3
) -> Mapping[str, NDArray[np.float64]]:
    """Construct local amount/volume perturbations without adding bulk charge."""

    baseline = np.asarray(state, dtype=float)
    if (
        baseline.ndim != 1
        or baseline.size < CORE_SIZE
        or not np.all(np.isfinite(baseline))
        or np.any(baseline[:CORE_SIZE] <= 0.0)
    ):
        raise ValueError("nearby-state baseline must contain a positive finite core")
    if not math.isfinite(fraction) or not 0.0 < fraction < 0.05:
        raise ValueError("nearby-state fraction must lie in (0, 0.05)")
    variants: dict[str, NDArray[np.float64]] = {"baseline": baseline.copy()}
    for label, index, factor in (
        ("cell_volume_down", 5, 1.0 - fraction),
        ("cell_volume_up", 5, 1.0 + fraction),
        ("lumen_volume_down", 11, 1.0 - fraction),
        ("lumen_volume_up", 11, 1.0 + fraction),
    ):
        candidate = baseline.copy()
        candidate[index] *= factor
        variants[label] = candidate
    cell_delta = fraction * min(float(baseline[0]), float(baseline[2]))
    lumen_delta = fraction * min(float(baseline[6]), float(baseline[8]))
    for label, cation_index, anion_index, delta, sign in (
        ("cell_nacl_down", 0, 2, cell_delta, -1.0),
        ("cell_nacl_up", 0, 2, cell_delta, 1.0),
        ("lumen_nacl_down", 6, 8, lumen_delta, -1.0),
        ("lumen_nacl_up", 6, 8, lumen_delta, 1.0),
    ):
        candidate = baseline.copy()
        candidate[cation_index] += sign * delta
        candidate[anion_index] += sign * delta
        if np.any(candidate[:CORE_SIZE] <= 0.0):
            raise ValueError(f"{label} produces a nonpositive core state")
        variants[label] = candidate
    baseline_cell_charge = baseline[0] + baseline[1] - baseline[2] - baseline[4]
    baseline_lumen_charge = baseline[6] + baseline[7] - baseline[8] - baseline[10]
    for label, candidate in variants.items():
        cell_charge = candidate[0] + candidate[1] - candidate[2] - candidate[4]
        lumen_charge = candidate[6] + candidate[7] - candidate[8] - candidate[10]
        if abs(cell_charge - baseline_cell_charge) > 1.0e-12:
            raise AssertionError(f"{label} changes intracellular bulk charge")
        if abs(lumen_charge - baseline_lumen_charge) > 1.0e-12:
            raise AssertionError(f"{label} changes luminal bulk charge")
    return variants


class IndependentWholeCell:
    """Second implementation of the modern whole-cell equations."""

    def __init__(
        self,
        parameters: FullModelParameters,
        *,
        ae4_parameters: Any | None = None,
        regulation: IndependentRegulation | None = None,
        nkcc_regulation: (
            IndependentN1NkccRegulation | IndependentN2NkccRegulation | None
        ) = None,
        stimulus: IndependentStimulus | None = None,
    ) -> None:
        self.parameters = parameters
        self.ae4_parameters = ae4_parameters
        self.regulation = regulation
        self.nkcc_regulation = nkcc_regulation
        self.stimulus = stimulus
        if (
            regulation is not None
            and nkcc_regulation is not None
            and set(regulation.state_names) & set(nkcc_regulation.state_names)
        ):
            raise ValueError("AE4 and NKCC1 regulatory state names must be disjoint")
        p = parameters
        bath_ta = independent_total_alkalinity_mM(
            p.bath.ph,
            p.bath.tic_mM,
            p.bath.buffer_total_mM,
            p.acid_base.lumen_buffer_pka,
            p,
        )
        self.bath_speciation = independent_speciate(
            p.bath.tic_mM,
            bath_ta,
            p.bath.buffer_total_mM,
            p.acid_base.lumen_buffer_pka,
            p,
        )

    @property
    def state_names(self) -> tuple[str, ...]:
        core = (
            "na_i_fmol", "k_i_fmol", "cl_i_fmol", "tic_i_fmol",
            "alk_i_fmol", "volume_i_pL", "na_l_fmol", "k_l_fmol",
            "cl_l_fmol", "tic_l_fmol", "alk_l_fmol", "volume_l_pL",
        )
        ae4 = () if self.regulation is None else self.regulation.state_names
        nkcc = () if self.nkcc_regulation is None else self.nkcc_regulation.state_names
        return core + ae4 + nkcc

    def _stimulus(self, time_s: float) -> tuple[float, float]:
        if self.stimulus is None:
            return 0.10, 0.0
        value = self.stimulus(time_s)
        if isinstance(value, Mapping):
            return float(value["calcium_uM"]), float(value["beta_input"])
        return float(value.calcium_uM), float(value.beta_input)

    def evaluate(
        self,
        time_s: float,
        state: ArrayLike,
        *,
        ae4_expression: float = 1.0,
        ae2_expression: float = 1.0,
    ) -> IndependentEvaluation:
        y = np.asarray(state, dtype=float)
        ae4_state_count = (
            0 if self.regulation is None else len(self.regulation.state_names)
        )
        nkcc_state_count = (
            0 if self.nkcc_regulation is None else len(self.nkcc_regulation.state_names)
        )
        expected = CORE_SIZE + ae4_state_count + nkcc_state_count
        if y.shape != (expected,) or not np.all(np.isfinite(y)) or np.any(y[:CORE_SIZE] <= 0.0):
            raise ValueError("independent state is not finite, positive, or correctly sized")
        p = self.parameters
        ci = {
            "na": y[0] / y[5], "k": y[1] / y[5], "cl": y[2] / y[5],
            "tic": y[3] / y[5], "alk": y[4] / y[5],
        }
        li = {
            "na": y[6] / y[11], "k": y[7] / y[11], "cl": y[8] / y[11],
            "tic": y[9] / y[11], "alk": y[10] / y[11],
        }
        cell_ab = independent_speciate(
            ci["tic"], ci["alk"], p.geometry.cell_buffer_total_fmol / y[5],
            p.acid_base.cell_buffer_pka, p,
        )
        lumen_ab = independent_speciate(
            li["tic"], li["alk"], p.geometry.lumen_buffer_total_fmol / y[11],
            p.acid_base.lumen_buffer_pka, p,
        )
        ci["hco3"] = cell_ab.hco3_mM
        ci["co2"] = cell_ab.co2_mM
        li["hco3"] = lumen_ab.hco3_mM
        li["co2"] = lumen_ab.co2_mM
        bath = {
            "na": p.bath.na_mM, "k": p.bath.k_mM, "cl": p.bath.cl_mM,
            "tic": p.bath.tic_mM, "hco3": self.bath_speciation.hco3_mM,
            "co2": self.bath_speciation.co2_mM, "h": self.bath_speciation.h_mM,
        }

        calcium, beta = self._stimulus(time_s)
        if self.regulation is None:
            capacity = 1.0
            regulatory_rhs: tuple[float, ...] = ()
        else:
            regulatory_state = tuple(
                float(value) for value in y[CORE_SIZE : CORE_SIZE + ae4_state_count]
            )
            capacity = self.regulation.capacity_multiplier(time_s, regulatory_state, beta)
            regulatory_rhs = self.regulation.rhs(time_s, regulatory_state, beta)
        if self.nkcc_regulation is None:
            nkcc_capacity = 1.0
            nkcc_activation = 0.0
            nkcc_regulatory_rhs: tuple[float, ...] = ()
        else:
            nkcc_state = tuple(
                float(value) for value in y[CORE_SIZE + ae4_state_count :]
            )
            if isinstance(self.nkcc_regulation, IndependentN1NkccRegulation):
                nkcc_activation = self.nkcc_regulation.normalized_calcium_arm(calcium)
                nkcc_capacity = self.nkcc_regulation.capacity_multiplier(
                    nkcc_state, calcium
                )
            else:
                nkcc_capacity = self.nkcc_regulation.capacity_multiplier(nkcc_state)
                nkcc_activation = self.nkcc_regulation.activation_fraction(nkcc_state)
            nkcc_regulatory_rhs = self.nkcc_regulation.rhs(
                time_s, nkcc_state, calcium
            )

        # Reversible homeostasis.
        width = p.homeostasis.thermodynamic_saturation_log_width
        a_nkcc = math.log(
            bath["na"] * bath["k"] * bath["cl"]**2
            / (ci["na"] * ci["k"] * ci["cl"]**2)
        )
        a_nhe = math.log(bath["na"] * cell_ab.h_mM / (ci["na"] * bath["h"]))
        a_ae2 = math.log(bath["cl"] * ci["hco3"] / (ci["cl"] * bath["hco3"]))
        j_nkcc = (
            p.homeostasis.nkcc1_capacity_fmol_s
            * nkcc_capacity
            * math.tanh(a_nkcc / width)
        )
        j_nhe = p.homeostasis.nhe1_capacity_fmol_s * math.tanh(a_nhe / width)
        j_ae2 = (
            p.homeostasis.ae2_capacity_fmol_s
            * ae2_expression
            * math.tanh(a_ae2 / width)
        )
        home = {
            "na": j_nkcc + j_nhe, "k": j_nkcc, "cl": 2.0 * j_nkcc + j_ae2,
            "tic": -j_ae2, "alk": j_nhe - j_ae2,
        }
        ae4 = _ae4_sources(ci, bath, self.ae4_parameters, capacity, ae4_expression)

        # Two-membrane electrical closure, solved analytically as a 2x2 system.
        mp = p.membranes
        jpa = _pump(
            mp.nak_capacity_fmol_s * mp.apical_pump_fraction,
            ci["na"], li["k"], mp.nak_na_half_mM, mp.nak_k_half_mM,
        )
        jpb = _pump(
            mp.nak_capacity_fmol_s * (1.0 - mp.apical_pump_fraction),
            ci["na"], bath["k"], mp.nak_na_half_mM, mp.nak_k_half_mM,
        )
        if calcium == 0.0:
            ca_gate = 0.0
        else:
            ca_power = calcium**mp.calcium_hill
            ca_gate = ca_power / (ca_power + mp.calcium_half_uM**mp.calcium_hill)
        gka = mp.g_k_total_S * ca_gate * mp.apical_k_fraction + mp.g_apical_background_S
        gkb = mp.g_k_total_S * ca_gate * (1.0 - mp.apical_k_fraction) + mp.g_basolateral_background_S
        gcla = mp.g_cl_apical_S * ca_gate
        vt = p.constants.thermal_voltage_V
        eka = _nernst(ci["k"], li["k"], +1, vt)
        ekb = _nernst(ci["k"], bath["k"], +1, vt)
        ecla = _nernst(ci["cl"], li["cl"], -1, vt)
        para = {
            "na": (mp.g_para_na_S, +1, _nernst(li["na"], bath["na"], +1, vt)),
            "k": (mp.g_para_k_S, +1, _nernst(li["k"], bath["k"], +1, vt)),
            "cl": (mp.g_para_cl_S, -1, _nernst(li["cl"], bath["cl"], -1, vt)),
            "hco3": (mp.g_para_hco3_S, -1, _nernst(li["hco3"], bath["hco3"], -1, vt)),
        }
        gp = sum(item[0] for item in para.values())
        ipa = p.constants.faraday_C_mol * jpa / FMOL_PER_MOL
        ipb = p.constants.faraday_C_mol * jpb / FMOL_PER_MOL
        ga = gka + gcla
        ca = -gka * eka - gcla * ecla + ipa
        cb = -gkb * ekb + ipb
        cp = -sum(g * reversal for g, _z, reversal in para.values())
        a11, a12, a21, a22 = ga + gp, -gp, -gp, gkb + gp
        b1, b2 = -ca + cp, -cb - cp
        determinant = a11 * a22 - a12 * a21
        if determinant == 0.0:
            raise ValueError("independent current closure is singular")
        va = (b1 * a22 - a12 * b2) / determinant
        vb = (a11 * b2 - b1 * a21) / determinant
        vtrans = vb - va
        ika = gka * (va - eka)
        ikb = gkb * (vb - ekb)
        icla = gcla * (va - ecla)
        ipara = {name: g * (vtrans - reversal) for name, (g, _z, reversal) in para.items()}
        ipara_total = math.fsum(ipara.values())
        ia_total = ika + icla + ipa
        ib_total = ikb + ipb
        jka = _current_to_flux(ika, +1, p.constants.faraday_C_mol)
        jkb = _current_to_flux(ikb, +1, p.constants.faraday_C_mol)
        jcla = _current_to_flux(icla, -1, p.constants.faraday_C_mol)
        jpara = {
            name: _current_to_flux(ipara[name], para[name][1], p.constants.faraday_C_mol)
            for name in para
        }
        elec_cell = {
            "na": -3.0 * (jpa + jpb),
            "k": 2.0 * (jpa + jpb) - jka - jkb,
            "cl": -jcla, "tic": 0.0, "alk": 0.0,
        }
        elec_lumen = {
            "na": 3.0 * jpa - jpara["na"],
            "k": -2.0 * jpa + jka - jpara["k"],
            "cl": jcla - jpara["cl"],
            "tic": -jpara["hco3"], "alk": -jpara["hco3"],
        }

        # Osmotic water balance.  The finite buffer is one osmole per molecule.
        osm_cell = (
            ci["na"] + ci["k"] + ci["cl"] + ci["tic"]
            + (p.geometry.cell_impermeant_osmoles_fmol + p.geometry.cell_buffer_total_fmol) / y[5]
        )
        osm_lumen = li["na"] + li["k"] + li["cl"] + li["tic"]
        osm_bath = (
            bath["na"] + bath["k"] + bath["cl"] + bath["tic"]
            + p.bath.untracked_osmolyte_mM
        )
        wp = p.water
        qb = wp.basolateral_hydraulic_pL_s_mOsm * (osm_cell - osm_bath)
        qa = wp.apical_hydraulic_pL_s_mOsm * (osm_lumen - osm_cell)
        qt = wp.paracellular_hydraulic_pL_s_mOsm * (osm_lumen - osm_bath)
        qout = wp.outflow_rate_s * max(y[11] - wp.lumen_dead_volume_pL, 0.0)
        dvi = qb - qa
        dvl = qa + qt - qout
        outflow = {name: qout * li[name] for name in ("na", "k", "cl", "tic", "alk")}
        jco2b = p.homeostasis.co2_basolateral_permeability_fmol_s_mM * (bath["co2"] - ci["co2"])
        jco2a = p.homeostasis.co2_apical_permeability_fmol_s_mM * (li["co2"] - ci["co2"])

        rhs = np.zeros(expected, dtype=float)
        rhs[0] = home["na"] + ae4["na"] + elec_cell["na"]
        rhs[1] = home["k"] + ae4["k"] + elec_cell["k"]
        rhs[2] = home["cl"] + ae4["cl"] + elec_cell["cl"]
        rhs[3] = home["tic"] + ae4["tic"] + jco2b + jco2a
        rhs[4] = home["alk"] + ae4["alk"]
        rhs[5] = dvi
        rhs[6] = elec_lumen["na"] - outflow["na"]
        rhs[7] = elec_lumen["k"] - outflow["k"]
        rhs[8] = elec_lumen["cl"] - outflow["cl"]
        rhs[9] = elec_lumen["tic"] - jco2a - outflow["tic"]
        rhs[10] = elec_lumen["alk"] - outflow["alk"]
        rhs[11] = dvl
        if regulatory_rhs:
            rhs[CORE_SIZE : CORE_SIZE + ae4_state_count] = regulatory_rhs
        if nkcc_regulatory_rhs:
            rhs[CORE_SIZE + ae4_state_count :] = nkcc_regulatory_rhs

        cell_charge_rate = rhs[0] + rhs[1] - rhs[2] - rhs[4]
        lumen_charge_rate = rhs[6] + rhs[7] - rhs[8] - rhs[10]
        outflow_charge = outflow["na"] + outflow["k"] - outflow["cl"] - outflow["alk"]
        external_carbon = home["tic"] + ae4["tic"] + jco2b + elec_lumen["tic"] - outflow["tic"]
        conservation = {
            "cell_charge_rate_fmol_s": cell_charge_rate,
            "lumen_charge_minus_outflow_fmol_s": lumen_charge_rate + outflow_charge,
            "carbon_accounting_fmol_s": rhs[3] + rhs[9] - external_carbon,
            "water_accounting_pL_s": dvi + dvl - (qb + qt - qout),
            "ae4_charge_fmol_s": ae4["charge"],
        }
        fluxes = {
            "nkcc1": j_nkcc, "nhe1": j_nhe, "ae2": j_ae2,
            "ae4_cl": ae4["cl"], "ae4_entropy_over_r": ae4["entropy_over_r"],
            "pump_apical": jpa, "pump_basolateral": jpb,
            "co2_bath_to_cell": jco2b, "co2_lumen_to_cell": jco2a,
        }
        return IndependentEvaluation(
            rhs=rhs,
            cell_concentrations_mM=ci,
            lumen_concentrations_mM=li,
            cell_speciation=cell_ab,
            lumen_speciation=lumen_ab,
            capacity_multiplier=capacity,
            nkcc_capacity_multiplier=nkcc_capacity,
            nkcc_activation_fraction=nkcc_activation,
            flow_pL_s=qout,
            voltages_V={"apical": va, "basolateral": vb, "transepithelial": vtrans},
            current_residuals_A={"apical": ia_total - ipara_total, "basolateral": ib_total + ipara_total},
            fluxes_fmol_s=fluxes,
            conservation_residuals=conservation,
        )

    def rhs(self, time_s: float, state: ArrayLike, *, ae4_expression: float = 1.0, ae2_expression: float = 1.0) -> NDArray[np.float64]:
        return self.evaluate(
            time_s, state, ae4_expression=ae4_expression, ae2_expression=ae2_expression
        ).rhs

    def refine_resting_root(
        self,
        candidate_state: ArrayLike,
        *,
        max_nfev: int = 5000,
        relative_bound: float = 0.20,
        row_scales: ArrayLike | None = None,
        ae4_expression: float = 1.0,
        ae2_expression: float = 1.0,
    ) -> IndependentRootResult:
        """Independently refine a WT root on its electroneutral charge leaf."""

        candidate = np.asarray(candidate_state, dtype=float)
        if candidate.shape != (len(self.state_names),):
            raise ValueError("candidate state does not match independent model")
        if not 0.0 < relative_bound < 0.90:
            raise ValueError("relative root bound must lie in (0, 0.9)")
        free0 = free_charge_coordinates(candidate)
        fixed = self.parameters.geometry.fixed_cell_anion_equivalents_fmol
        regulation = tuple(candidate[CORE_SIZE:])
        lower = np.maximum(free0 * (1.0 - relative_bound), 1.0e-12)
        upper = free0 * (1.0 + relative_bound)
        scales = (
            np.maximum(np.abs(candidate[list(INDEPENDENT_CORE_ROWS)]) / 100.0, 1.0e-9)
            if row_scales is None else np.asarray(row_scales, dtype=float)
        )
        if scales.shape != (len(INDEPENDENT_CORE_ROWS),) or np.any(scales <= 0.0):
            raise ValueError("root row scales must be ten positive values")

        def residual(free: NDArray[np.float64]) -> NDArray[np.float64]:
            state = state_from_free_charge_coordinates(
                free, fixed_cell_anion_equivalents_fmol=fixed, regulation=regulation
            )
            if np.any(state[:CORE_SIZE] <= 0.0):
                return np.full(len(INDEPENDENT_CORE_ROWS), 1.0e9)
            try:
                raw = self.rhs(
                    0.0,
                    state,
                    ae4_expression=ae4_expression,
                    ae2_expression=ae2_expression,
                )
            except (ValueError, OverflowError, FloatingPointError):
                return np.full(len(INDEPENDENT_CORE_ROWS), 1.0e9)
            return raw[list(INDEPENDENT_CORE_ROWS)] / scales

        fit = least_squares(
            residual,
            np.clip(free0, lower, upper),
            bounds=(lower, upper),
            max_nfev=max_nfev,
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
        )
        state = state_from_free_charge_coordinates(
            fit.x, fixed_cell_anion_equivalents_fmol=fixed, regulation=regulation
        )
        raw = self.rhs(
            0.0,
            state,
            ae4_expression=ae4_expression,
            ae2_expression=ae2_expression,
        )
        independent_raw = raw[list(INDEPENDENT_CORE_ROWS)]
        scaled = independent_raw / scales
        # ``fit.jac`` already has dimensionless scaled-residual rows.  Multiply
        # each column by its declared free-coordinate scale before computing
        # singular values, so amount (fmol) and volume (pL) columns are not
        # silently compared in their native numerical magnitudes.
        coordinate_scales = np.maximum(np.abs(free0), 1.0e-12)
        scaled_jacobian = np.asarray(fit.jac, dtype=float) @ np.diag(coordinate_scales)
        singular = np.linalg.svd(scaled_jacobian, compute_uv=False)
        threshold = max(fit.jac.shape) * np.finfo(float).eps * max(singular, default=0.0)
        rank = int(np.count_nonzero(singular > threshold))
        regulatory_max = (
            float(np.max(np.abs(raw[CORE_SIZE:]))) if raw.size > CORE_SIZE else 0.0
        )
        return IndependentRootResult(
            state=state,
            success=bool(
                fit.success
                and np.max(np.abs(scaled)) <= 1.0e-6
                and regulatory_max <= 1.0e-9
            ),
            max_abs_scaled_residual=float(np.max(np.abs(scaled))),
            max_abs_amount_source_fmol_s=float(
                np.max(np.abs(raw[[0, 1, 3, 4, 6, 7, 9, 10]]))
            ),
            max_abs_volume_source_pL_s=float(np.max(np.abs(raw[[5, 11]]))),
            max_abs_regulatory_rhs_s_inv=regulatory_max,
            nfev=int(fit.nfev),
            cost=float(fit.cost),
            scaled_jacobian_singular_values=tuple(float(value) for value in singular),
            jacobian_rank=rank,
            jacobian_nullity=int(fit.jac.shape[1] - rank),
            message=str(fit.message),
        )

    def integrate(
        self,
        initial_state: ArrayLike,
        time_s: ArrayLike,
        *,
        method: str,
        rtol: float = 1.0e-7,
        atol: float | ArrayLike = 1.0e-10,
        max_step_s: float = 2.0,
        ae4_expression: float = 1.0,
        ae2_expression: float = 1.0,
        preserve_basal_zero: bool = False,
    ) -> IndependentTrajectory:
        """Integrate the independent equations with Radau or BDF.

        A scalar ``atol`` is interpreted by SciPy in each coordinate's native
        unit (fmol, pL, or dimensionless regulatory fraction).  Production
        crosschecks should therefore also include a vector tolerance derived
        from the frozen state scales.
        """

        if method not in {"Radau", "BDF"}:
            raise ValueError("independent trajectory requires Radau or BDF")
        grid = np.asarray(time_s, dtype=float)
        if grid.ndim != 1 or grid.size < 2 or np.any(np.diff(grid) <= 0.0):
            raise ValueError("time grid must be strictly increasing")
        y0 = np.asarray(initial_state, dtype=float)
        if preserve_basal_zero:
            if grid.size < 3 or grid[0] != 0.0 or grid[1] <= 0.0:
                raise ValueError(
                    "basal-zero handling requires zero, a positive right-limit, and an endpoint"
                )
            integration_grid = grid[1:]
            integration_span = (float(integration_grid[0]), float(integration_grid[-1]))
        else:
            integration_grid = grid
            integration_span = (float(grid[0]), float(grid[-1]))
        result = solve_ivp(
            lambda t, y: self.rhs(
                t, y, ae4_expression=ae4_expression, ae2_expression=ae2_expression
            ),
            integration_span,
            y0,
            method=method,
            t_eval=integration_grid,
            rtol=rtol,
            atol=atol,
            max_step=max_step_s,
        )
        complete = result.y.shape[1] == integration_grid.size
        if preserve_basal_zero:
            times = np.concatenate((np.asarray((0.0,)), np.asarray(result.t, dtype=float)))
            states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
        else:
            times = np.asarray(result.t, dtype=float)
            states = np.asarray(result.y, dtype=float)
        flow = np.empty(times.size, dtype=float)
        chloride = np.empty(times.size, dtype=float)
        ph = np.empty(times.size, dtype=float)
        capacity = np.empty(times.size, dtype=float)
        nkcc_capacity = np.empty(times.size, dtype=float)
        nkcc_activation = np.empty(times.size, dtype=float)
        conservation_fmol_max = 0.0
        water_residual_max = 0.0
        for index, (time_value, vector) in enumerate(zip(times, states.T)):
            evaluated = self.evaluate(
                float(time_value), vector,
                ae4_expression=ae4_expression, ae2_expression=ae2_expression,
            )
            flow[index] = evaluated.flow_pL_s
            chloride[index] = evaluated.cell_concentrations_mM["cl"]
            ph[index] = evaluated.cell_speciation.ph
            capacity[index] = evaluated.capacity_multiplier
            nkcc_capacity[index] = evaluated.nkcc_capacity_multiplier
            nkcc_activation[index] = evaluated.nkcc_activation_fraction
            conservation_fmol_max = max(
                conservation_fmol_max,
                max(
                    abs(float(value))
                    for name, value in evaluated.conservation_residuals.items()
                    if name != "water_accounting_pL_s"
                ),
            )
            water_residual_max = max(
                water_residual_max,
                abs(float(evaluated.conservation_residuals["water_accounting_pL_s"])),
            )
        cumulative = np.zeros(times.size, dtype=float)
        if times.size > 1:
            cumulative[1:] = np.cumsum(
                0.5 * (flow[:-1] + flow[1:]) * np.diff(times)
            )
        return IndependentTrajectory(
            method=method,
            onset_handling=(
                "ANALYTIC_BASAL_ZERO_THEN_POSITIVE_RIGHT_LIMIT"
                if preserve_basal_zero
                else "CONTINUOUS_GRID"
            ),
            success=bool(result.success and complete and times.size == grid.size),
            message=str(result.message),
            time_s=times,
            states=states,
            flow_pL_s=flow,
            cumulative_flow_pL=cumulative,
            cell_cl_mM=chloride,
            cell_ph=ph,
            capacity_multiplier=capacity,
            nkcc_capacity_multiplier=nkcc_capacity,
            nkcc_activation_fraction=nkcc_activation,
            max_abs_conservation_fmol_s=float(conservation_fmol_max),
            max_abs_water_accounting_pL_s=float(water_residual_max),
        )


__all__ = (
    "CORE_SIZE",
    "FREE_CORE_INDICES",
    "INDEPENDENT_CORE_ROWS",
    "IndependentEvaluation",
    "IndependentAE4Kinetics",
    "IndependentCalibrationDefinition",
    "IndependentConstantStimulus",
    "IndependentN1NkccRegulation",
    "IndependentN2NkccRegulation",
    "IndependentSecretagogueProtocol",
    "IndependentRegulation",
    "IndependentRootResult",
    "IndependentSpeciation",
    "IndependentTrajectory",
    "IndependentWholeCell",
    "free_charge_coordinates",
    "independent_nearby_charge_preserving_states",
    "independent_r0_regulation",
    "independent_r1_regulation",
    "independent_r2_regulation",
    "independent_r3_regulation",
    "independent_physical_time_grid",
    "independent_speciate",
    "independent_total_alkalinity_mM",
    "state_from_free_charge_coordinates",
)
