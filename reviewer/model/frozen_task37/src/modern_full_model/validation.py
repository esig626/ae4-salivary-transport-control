"""Pre-reveal WT dynamic validation for the Task 13B full model.

This module is intentionally target-firewalled.  It contains no AE4-null
secretion observation and it never evaluates a knockout trajectory.  It
provides the physical stimulation protocol, deterministic nearby-state
perturbations, dual-method solver comparisons, WT flow scaling, and immutable
manifest helpers needed before the held-out data can be opened.

Scientific scope
----------------
The 2015 gland protocol used 0.3 uM carbachol (CCh) plus 5 uM
isoproterenol (IPR) for 10 minutes.  CCh is represented only through the
calcium input to the membrane module; IPR is represented only through the
beta input to the regulatory module.  They are never substituted for one
another.  The historical-lineage approximately 0.05-to-0.55 uM calcium step
is retained as an explicitly labelled input sensitivity because the primary
gland study did not measure the intracellular calcium waveform.  The 0.058 uM
resting value matches the frozen WT-root convention; its old code-time onset
is discarded.  This module's time coordinate is seconds and the source
protocol duration is 600 s.

The 2021 study establishes beta/adenylate-cyclase/PKA-dependent AE4
activation and S173 dependence, but it does not report an acute cAMP, PKA,
AE4-regulatory-state, or phosphorylation time course.  Consequently this
module reports regulatory time constants but does not optimize them or score
them against a fabricated onset target.

Provenance categories: ``PRIMARY_MEASUREMENT`` (agonist doses and duration,
WT flow envelope), ``HISTORICAL_IMPLEMENTATION`` (calcium amplitudes),
``DERIVED_CONSTRAINT`` (charge-preserving perturbations and unit conversion),
and ``NEW_MODELING_DECISION`` (solver/tolerance and sensitivity choices).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass, replace
from enum import Enum
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid

from .model import ModernFullModel, Stimulus, WT
from .states import CORE_STATE_NAMES, amount_charge_equivalents_fmol


class StimulusArm(str, Enum):
    """Experimentally distinct secretagogue arms."""

    REST = "REST"
    CCH_ONLY = "CCH_ONLY"
    IPR_ONLY = "IPR_ONLY"
    CCH_IPR = "CCH_IPR"


@dataclass(frozen=True)
class SecretagogueProtocol:
    """Physically timed CCh/Ca and IPR/beta protocol.

    At the exact onset coordinate the basal input is returned so a resting
    initial condition is well defined.  The stimulated right limit is used at
    every positive time after onset.  This convention is recorded in the
    frozen manifest and is tested explicitly.
    """

    arm: StimulusArm = StimulusArm.CCH_IPR
    onset_s: float = 0.0
    duration_s: float = 600.0
    cch_dose_uM: float = 0.3
    ipr_dose_uM: float = 5.0
    resting_calcium_uM: float = 0.058
    stimulated_calcium_uM: float = 0.55
    beta_occupancy_on: float = 1.0

    def __post_init__(self) -> None:
        values = (
            self.onset_s,
            self.duration_s,
            self.cch_dose_uM,
            self.ipr_dose_uM,
            self.resting_calcium_uM,
            self.stimulated_calcium_uM,
            self.beta_occupancy_on,
        )
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("stimulus values must be finite")
        if self.onset_s < 0.0 or self.duration_s <= 0.0:
            raise ValueError("onset must be nonnegative and duration positive")
        if self.cch_dose_uM < 0.0 or self.ipr_dose_uM < 0.0:
            raise ValueError("agonist doses must be nonnegative")
        if self.resting_calcium_uM < 0.0 or self.stimulated_calcium_uM < 0.0:
            raise ValueError("calcium inputs must be nonnegative")
        if not 0.0 <= self.beta_occupancy_on <= 1.0:
            raise ValueError("beta occupancy must lie in [0, 1]")

    @property
    def end_s(self) -> float:
        return self.onset_s + self.duration_s

    @property
    def contains_cch(self) -> bool:
        return self.arm in (StimulusArm.CCH_ONLY, StimulusArm.CCH_IPR)

    @property
    def contains_ipr(self) -> bool:
        return self.arm in (StimulusArm.IPR_ONLY, StimulusArm.CCH_IPR)

    def is_active(self, time_s: float) -> bool:
        if not math.isfinite(time_s):
            raise ValueError("time must be finite")
        return bool(self.onset_s < time_s <= self.end_s)

    def cch_input(self, time_s: float) -> float:
        return float(self.cch_dose_uM if self.contains_cch and self.is_active(time_s) else 0.0)

    def ipr_input(self, time_s: float) -> float:
        return float(self.ipr_dose_uM if self.contains_ipr and self.is_active(time_s) else 0.0)

    def calcium_input(self, time_s: float) -> float:
        if self.contains_cch and self.is_active(time_s):
            return float(self.stimulated_calcium_uM)
        return float(self.resting_calcium_uM)

    def beta_input(self, time_s: float) -> float:
        if self.contains_ipr and self.is_active(time_s):
            return float(self.beta_occupancy_on)
        return 0.0

    def __call__(self, time_s: float) -> Stimulus:
        return Stimulus(
            calcium_uM=self.calcium_input(time_s),
            beta_input=self.beta_input(time_s),
        )

    def provenance(self) -> Mapping[str, str]:
        return {
            "agonist_doses_and_duration": "PRIMARY_MEASUREMENT:P15",
            "cch_to_calcium_arm": "SOURCE_SUPPORTED_MECHANISTIC_CONSTRAINT",
            "ipr_to_beta_arm": "PRIMARY_MEASUREMENT:P15,P21",
            "calcium_amplitudes": "HISTORICAL_IMPLEMENTATION_SENSITIVITY",
            "calcium_onset_shape": "NEW_MODELING_DECISION:instantaneous_step",
            "beta_occupancy_scale": "NEW_MODELING_DECISION:normalized_input",
            "physical_time": "PRIMARY_PROTOCOL_SECONDS; no inherited rhs_scale",
        }


@dataclass(frozen=True)
class SolverSpecification:
    """Numerical choices frozen before any held-out evaluation."""

    label: str
    method: str
    rtol: float
    amount_atol_fmol: float
    volume_atol_pL: float
    regulatory_atol_fraction: float
    max_step_s: float

    def __post_init__(self) -> None:
        if self.method not in {"Radau", "BDF"}:
            raise ValueError("only the declared stiff methods Radau and BDF are allowed")
        if (
            self.rtol <= 0.0
            or self.amount_atol_fmol <= 0.0
            or self.volume_atol_pL <= 0.0
            or self.regulatory_atol_fraction <= 0.0
            or self.max_step_s <= 0.0
        ):
            raise ValueError("solver tolerances and max step must be positive")

    def atol_vector(self, state_names: Sequence[str]) -> NDArray[np.float64]:
        """Return absolute tolerances in each coordinate's native unit."""

        values = []
        for name in state_names:
            if name.endswith("_fmol"):
                values.append(self.amount_atol_fmol)
            elif name.endswith("_pL"):
                values.append(self.volume_atol_pL)
            else:
                values.append(self.regulatory_atol_fraction)
        return np.asarray(values, dtype=float)


PRODUCTION_RADAU = SolverSpecification(
    "production_radau", "Radau", 1.0e-7, 1.0e-10, 1.0e-12, 1.0e-10, 2.0
)
PRODUCTION_BDF = SolverSpecification(
    "production_bdf", "BDF", 1.0e-7, 1.0e-10, 1.0e-12, 1.0e-10, 2.0
)
LOOSE_RADAU = SolverSpecification(
    "loose_radau", "Radau", 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-8, 5.0
)
TIGHT_RADAU = SolverSpecification(
    "tight_radau", "Radau", 1.0e-8, 1.0e-11, 1.0e-13, 1.0e-11, 1.0
)


# Each conservation diagnostic retains its physical unit.  These numerical
# acceptance choices are frozen before held-out evaluation; they are not
# biological uncertainty bands.
CONSERVATION_RESIDUAL_UNITS: Mapping[str, str] = {
    "cell_bulk_charge_rate_fmol_s": "fmol charge equivalents s^-1",
    "lumen_bulk_charge_rate_minus_outflow_fmol_s": "fmol charge equivalents s^-1",
    "carbon_accounting_fmol_s": "fmol s^-1",
    "buffer_site_accounting_fmol_s": "fmol sites s^-1",
    "water_volume_accounting_pL_s": "pL s^-1",
    "homeostasis_charge_fmol_s": "fmol charge equivalents s^-1",
    "ae4_charge_fmol_s": "fmol charge equivalents s^-1",
    "apical_current_A": "A",
    "basolateral_current_A": "A",
    "cell_speciation_alkalinity_mM": "mM",
    "lumen_speciation_alkalinity_mM": "mM",
}

CONSERVATION_RESIDUAL_TOLERANCES: Mapping[str, float] = {
    "cell_bulk_charge_rate_fmol_s": 1.0e-10,
    "lumen_bulk_charge_rate_minus_outflow_fmol_s": 1.0e-10,
    "carbon_accounting_fmol_s": 1.0e-10,
    "buffer_site_accounting_fmol_s": 1.0e-10,
    "water_volume_accounting_pL_s": 1.0e-12,
    "homeostasis_charge_fmol_s": 1.0e-10,
    "ae4_charge_fmol_s": 1.0e-10,
    "apical_current_A": 1.0e-20,
    "basolateral_current_A": 1.0e-20,
    "cell_speciation_alkalinity_mM": 1.0e-9,
    "lumen_speciation_alkalinity_mM": 1.0e-9,
}


@dataclass(frozen=True)
class WTFlowEnvelope:
    """Permitted WT-only whole-gland flow scale from P15 Fig. 1A."""

    lower_uL_min: float = 9.0
    upper_uL_min: float = 10.0
    scale_target_uL_min: float = 9.5
    source: str = "P15 Fig. 1A graphical approximation; n=6 glands"

    def __post_init__(self) -> None:
        if not 0.0 < self.lower_uL_min <= self.scale_target_uL_min <= self.upper_uL_min:
            raise ValueError("WT flow envelope must contain its positive scale target")


@dataclass(frozen=True)
class FlowScale:
    """One WT-only map from single-cell pL/s to whole-gland uL/min."""

    multiplier_uL_min_per_pL_s: float
    effective_cell_count: float
    calibration_statistic: str
    target_uL_min: float
    model_statistic_pL_s: float
    source: str

    def convert(self, flow_pL_s: ArrayLike) -> NDArray[np.float64]:
        return np.asarray(flow_pL_s, dtype=float) * self.multiplier_uL_min_per_pL_s


@dataclass(frozen=True)
class FlowEnvelopeFeasibility:
    """Scale-free compatibility of minute landmarks with the WT envelope."""

    compatible: bool
    sampled_min_flow_pL_s: float
    sampled_max_flow_pL_s: float
    model_max_min_ratio: float
    envelope_max_min_ratio: float
    feasible_multiplier_lower_uL_min_per_pL_s: float
    feasible_multiplier_upper_uL_min_per_pL_s: float
    sampled_times_s: tuple[float, ...]


@dataclass(frozen=True)
class Trajectory:
    """One solver result plus decoded WT observables."""

    family: str
    initial_condition: str
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
    capacity_multiplier: NDArray[np.float64]
    max_abs_conservation_residuals: Mapping[str, float]

    @property
    def positive_core(self) -> bool:
        return bool(np.all(self.states[: len(CORE_STATE_NAMES), :] > 0.0))


@dataclass(frozen=True)
class ArmResponse:
    """WT onset diagnostic for one experimentally distinct agonist arm."""

    family: str
    arm: str
    regulatory_initialization: str
    calcium_uM: float
    beta_input: float
    capacity_multiplier: float
    cell_cl_mM: float
    d_cell_cl_mM_s: float
    fractional_cell_cl_rate_s: float
    flow_pL_s: float


@dataclass(frozen=True)
class RegulatoryEnsembleMember:
    """One predeclared, non-fitted R0--R3/gain/timing sensitivity."""

    member_id: str
    family: str
    gain_label: str
    fully_activated_multiplier: float
    timing_label: str
    timing_parameters: Mapping[str, float]
    kinetic_evidence_status: str
    regulatory_model: Any

    def metadata(self) -> Mapping[str, Any]:
        return {
            "member_id": self.member_id,
            "family": self.family,
            "gain_label": self.gain_label,
            "fully_activated_multiplier": self.fully_activated_multiplier,
            "timing_label": self.timing_label,
            "timing_parameters": dict(self.timing_parameters),
            "kinetic_evidence_status": self.kinetic_evidence_status,
            "regulatory_model_class": type(self.regulatory_model).__name__,
        }


def pre_reveal_regulatory_ensemble() -> tuple[RegulatoryEnsembleMember, ...]:
    """Return the frozen R0--R3, gain, and timing sensitivity ensemble.

    Neither the two gain readings nor the fast/reference/slow kinetic choices
    are point estimates.  They are declared before held-out evaluation to
    distinguish model-class eligibility from kinetic identifiability.
    """

    from .camp_pka import (
        R0Static,
        R1EffectiveActivation,
        R2CampActivation,
        R3PkaRegulatedFraction,
        RegulatoryGain,
    )

    gains = (
        ("G125_P21_PROSE", 0.25),
        ("G170_P21_FIGURE_VISUAL", 0.70),
    )
    timings = (
        ("TFAST", 3.0, 10.0),
        ("TREFERENCE", 10.0, 30.0),
        ("TSLOW", 30.0, 90.0),
    )
    members: list[RegulatoryEnsembleMember] = []
    for gain_label, increment in gains:
        gain = RegulatoryGain(fully_activated_increment=increment)
        multiplier = 1.0 + increment
        members.append(
            RegulatoryEnsembleMember(
                member_id=f"R0_{gain_label}_STATIC",
                family="R0",
                gain_label=gain_label,
                fully_activated_multiplier=multiplier,
                timing_label="STATIC",
                timing_parameters={},
                kinetic_evidence_status="STATIC_NEGATIVE_CONTROL",
                regulatory_model=R0Static(gain=gain),
            )
        )
        for timing_label, upstream_s, downstream_s in timings:
            common_status = "UNMEASURED_PREDECLARED_SENSITIVITY_NOT_FIT"
            members.extend(
                (
                    RegulatoryEnsembleMember(
                        member_id=f"R1_{gain_label}_{timing_label}",
                        family="R1",
                        gain_label=gain_label,
                        fully_activated_multiplier=multiplier,
                        timing_label=timing_label,
                        timing_parameters={"tau_activation_s": downstream_s},
                        kinetic_evidence_status=common_status,
                        regulatory_model=R1EffectiveActivation(
                            tau_activation_s=downstream_s, gain=gain
                        ),
                    ),
                    RegulatoryEnsembleMember(
                        member_id=f"R2_{gain_label}_{timing_label}",
                        family="R2",
                        gain_label=gain_label,
                        fully_activated_multiplier=multiplier,
                        timing_label=timing_label,
                        timing_parameters={
                            "tau_camp_s": upstream_s,
                            "tau_activation_s": downstream_s,
                        },
                        kinetic_evidence_status=common_status,
                        regulatory_model=R2CampActivation(
                            tau_camp_s=upstream_s,
                            tau_activation_s=downstream_s,
                            gain=gain,
                        ),
                    ),
                    RegulatoryEnsembleMember(
                        member_id=f"R3_{gain_label}_{timing_label}",
                        family="R3",
                        gain_label=gain_label,
                        fully_activated_multiplier=multiplier,
                        timing_label=timing_label,
                        timing_parameters={
                            "tau_pka_s": upstream_s,
                            "forward_regulation_rate_s": 1.0 / downstream_s,
                            "reverse_regulation_rate_s": 1.0 / downstream_s,
                        },
                        kinetic_evidence_status=common_status,
                        regulatory_model=R3PkaRegulatedFraction(
                            tau_pka_s=upstream_s,
                            forward_regulation_rate_s=1.0 / downstream_s,
                            reverse_regulation_rate_s=1.0 / downstream_s,
                            gain=gain,
                        ),
                    ),
                )
            )
    if len({member.member_id for member in members}) != len(members):
        raise AssertionError("regulatory ensemble member IDs must be unique")
    return tuple(members)


def physical_time_grid(duration_s: float = 600.0, sample_step_s: float = 1.0) -> NDArray[np.float64]:
    """Return a seconds grid that resolves the right side of the onset step."""

    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise ValueError("duration must be finite and positive")
    if not math.isfinite(sample_step_s) or sample_step_s <= 0.0:
        raise ValueError("sample step must be finite and positive")
    regular = np.arange(sample_step_s, duration_s, sample_step_s, dtype=float)
    return np.concatenate((np.asarray((0.0, 1.0e-6)), regular, np.asarray((duration_s,))))


def _regulatory_initial_state(model: ModernFullModel, *, beta_input: float = 0.0) -> tuple[float, ...]:
    regulatory = model.regulatory_model
    if regulatory is None:
        return ()
    value = regulatory.initial_state(beta_input=beta_input)
    if hasattr(value, "as_vector"):
        value = value.as_vector()
    elif hasattr(value, "values") and not isinstance(value, Mapping):
        value = value.values
    if isinstance(value, Mapping):
        return tuple(float(value[name]) for name in model.layout.regulatory_names)
    return tuple(float(item) for item in value)


def attach_basal_regulation(model: ModernFullModel, state: ArrayLike) -> NDArray[np.float64]:
    """Attach or reset the regulatory suffix to its beta-free basal state."""

    candidate = np.asarray(state, dtype=float)
    core_size = len(CORE_STATE_NAMES)
    if candidate.shape == (core_size,):
        full = np.concatenate((candidate, np.asarray(_regulatory_initial_state(model))))
    elif candidate.shape == (model.layout.size,):
        full = candidate.copy()
        if model.layout.size > core_size:
            full[core_size:] = _regulatory_initial_state(model)
    else:
        raise ValueError(
            f"root state must have {core_size} core or {model.layout.size} full coordinates"
        )
    return model.layout.validate(full)


def nearby_charge_preserving_states(
    model: ModernFullModel,
    state: ArrayLike,
    *,
    fraction: float = 1.0e-3,
) -> Mapping[str, NDArray[np.float64]]:
    """Construct deterministic nearby WT states without injecting bulk charge.

    Volume perturbations change concentrations without changing conserved
    amount charge.  NaCl-pair perturbations add/remove equal cation and anion
    amounts.  Regulatory coordinates remain at their frozen basal values.
    """

    if not math.isfinite(fraction) or not 0.0 < fraction < 0.05:
        raise ValueError("nearby-state fraction must lie in (0, 0.05)")
    baseline = attach_basal_regulation(model, state)
    variants: dict[str, NDArray[np.float64]] = {"baseline": baseline.copy()}
    for label, factor in (("cell_volume_down", 1.0 - fraction), ("cell_volume_up", 1.0 + fraction)):
        candidate = baseline.copy()
        candidate[5] *= factor
        variants[label] = model.layout.validate(candidate)
    for label, factor in (
        ("lumen_volume_down", 1.0 - fraction),
        ("lumen_volume_up", 1.0 + fraction),
    ):
        candidate = baseline.copy()
        candidate[11] *= factor
        variants[label] = model.layout.validate(candidate)
    delta = fraction * min(float(baseline[0]), float(baseline[2]))
    for label, sign in (("cell_nacl_down", -1.0), ("cell_nacl_up", 1.0)):
        candidate = baseline.copy()
        candidate[0] += sign * delta
        candidate[2] += sign * delta
        variants[label] = model.layout.validate(candidate)
    lumen_delta = fraction * min(float(baseline[6]), float(baseline[8]))
    for label, sign in (("lumen_nacl_down", -1.0), ("lumen_nacl_up", 1.0)):
        candidate = baseline.copy()
        candidate[6] += sign * lumen_delta
        candidate[8] += sign * lumen_delta
        variants[label] = model.layout.validate(candidate)

    decoded_base = model.layout.decode(baseline)
    base_cell_charge = amount_charge_equivalents_fmol(
        decoded_base.cell,
        fixed_anion_equivalents_fmol=(
            model.parameters.geometry.fixed_cell_anion_equivalents_fmol
        ),
    )
    base_lumen_charge = amount_charge_equivalents_fmol(decoded_base.lumen)
    for label, candidate in variants.items():
        decoded = model.layout.decode(candidate)
        cell_charge = amount_charge_equivalents_fmol(
            decoded.cell,
            fixed_anion_equivalents_fmol=(
                model.parameters.geometry.fixed_cell_anion_equivalents_fmol
            ),
        )
        lumen_charge = amount_charge_equivalents_fmol(decoded.lumen)
        if abs(cell_charge - base_cell_charge) > 1.0e-12:
            raise AssertionError(f"{label} changes intracellular bulk charge")
        if abs(lumen_charge - base_lumen_charge) > 1.0e-12:
            raise AssertionError(f"{label} changes luminal bulk charge")
    return variants


def wt_arm_response(
    template: ModernFullModel,
    state: ArrayLike,
    *,
    family: str,
    arm: StimulusArm,
    regulatory_initialization: str = "basal_at_onset",
) -> ArmResponse:
    """Evaluate a WT-only CCh, IPR, or combined onset response.

    The reported fractional chloride rate is an exact model concentration
    derivative.  It is not silently equated with an SPQ fluorescence slope;
    such a comparison requires an explicit assay observation map.  The
    optional ``beta_preconditioned`` state is useful for the endpoint-style
    2015/2021 uptake contrasts, whose pretreatment timing does not identify an
    acute regulatory time constant.
    """

    if arm is StimulusArm.REST:
        time_s = 0.0
    else:
        time_s = 1.0e-6
    if isinstance(template.stimulus, SecretagogueProtocol):
        protocol = replace(template.stimulus, arm=arm)
    else:
        protocol = SecretagogueProtocol(arm=arm)
    model = ModernFullModel(
        parameters=template.parameters,
        stimulus=protocol,
        regulatory_model=template.regulatory_model,
        ae4_parameters=template.ae4_parameters,
        ae4_evaluator=template.ae4_evaluator,
    )
    vector = attach_basal_regulation(model, state)
    if regulatory_initialization == "beta_preconditioned":
        if model.layout.size > len(CORE_STATE_NAMES):
            vector[len(CORE_STATE_NAMES) :] = _regulatory_initial_state(
                model, beta_input=protocol.beta_input(time_s)
            )
    elif regulatory_initialization != "basal_at_onset":
        raise ValueError(
            "regulatory initialization must be basal_at_onset or beta_preconditioned"
        )
    evaluation = model.evaluate(time_s, vector, genotype=WT)
    decoded = model.layout.decode(vector)
    chloride_mM = decoded.cell.cl_fmol / decoded.cell.volume_pL
    d_chloride_mM_s = (
        evaluation.rhs[2] - chloride_mM * evaluation.rhs[5]
    ) / decoded.cell.volume_pL
    regulatory = evaluation.diagnostics.regulatory
    capacity_multiplier = float(
        regulatory.get("capacity_multiplier", regulatory.get("activity", 1.0))
    )
    stimulus = evaluation.diagnostics.stimulus
    return ArmResponse(
        family=str(family),
        arm=arm.value,
        regulatory_initialization=regulatory_initialization,
        calcium_uM=float(stimulus.calcium_uM),
        beta_input=float(stimulus.beta_input),
        capacity_multiplier=capacity_multiplier,
        cell_cl_mM=float(chloride_mM),
        d_cell_cl_mM_s=float(d_chloride_mM_s),
        fractional_cell_cl_rate_s=float(d_chloride_mM_s / chloride_mM),
        flow_pL_s=float(evaluation.diagnostics.water.lumen_outflow_pL_s),
    )


def simulate_wt(
    model: ModernFullModel,
    initial_state: ArrayLike,
    *,
    family: str,
    initial_condition: str,
    solver: SolverSpecification,
    time_s: ArrayLike | None = None,
) -> Trajectory:
    """Integrate and decode one strictly WT trajectory."""

    grid = physical_time_grid() if time_s is None else np.asarray(time_s, dtype=float)
    if grid.ndim != 1 or grid.size < 3 or grid[0] != 0.0 or np.any(np.diff(grid) <= 0.0):
        raise ValueError(
            "time grid must contain basal zero, a positive right-limit, and an endpoint"
        )
    y0 = attach_basal_regulation(model, initial_state)
    # The source protocol is discontinuous at stimulation onset.  Preserve the
    # exact basal row at t=0 analytically, then start the numerical IVP at the
    # positive right-limit where both inputs are smooth.  This prevents an
    # implicit solver stage from straddling the discontinuity and spuriously
    # stepping a bounded regulatory fraction outside [0, 1].
    integration_grid = grid[1:]
    result = model.solve_dynamics(
        (float(integration_grid[0]), float(integration_grid[-1])),
        initial_state=y0,
        genotype=WT,
        method=solver.method,
        rtol=solver.rtol,
        atol=solver.atol_vector(model.state_names),
        t_eval=integration_grid,
        max_step_s=solver.max_step_s,
    )
    complete = result.y.shape[1] == integration_grid.size
    if complete:
        states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    else:
        # Preserve the basal row and partial result with exact alignment.
        grid = np.concatenate((np.asarray((0.0,)), np.asarray(result.t, dtype=float)))
        states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    count = grid.size
    flow = np.full(count, np.nan, dtype=float)
    na = np.full(count, np.nan, dtype=float)
    k = np.full(count, np.nan, dtype=float)
    cl = np.full(count, np.nan, dtype=float)
    ph = np.full(count, np.nan, dtype=float)
    volume = np.full(count, np.nan, dtype=float)
    multiplier = np.full(count, np.nan, dtype=float)
    conservation_max = {name: 0.0 for name in CONSERVATION_RESIDUAL_UNITS}
    for index, (time_value, vector) in enumerate(zip(grid, states.T)):
        evaluation = model.evaluate(float(time_value), vector, genotype=WT)
        diagnostics = evaluation.diagnostics
        ci = diagnostics.observables.cell_concentrations_mM
        flow[index] = diagnostics.water.lumen_outflow_pL_s
        na[index] = ci["na"]
        k[index] = ci["k"]
        cl[index] = ci["cl"]
        ph[index] = diagnostics.observables.cell_acid_base.ph
        volume[index] = vector[5]
        multiplier[index] = float(
            diagnostics.regulatory.get(
                "capacity_multiplier", diagnostics.regulatory.get("activity", 1.0)
            )
        )
        unknown = set(diagnostics.conservation_residuals) - set(conservation_max)
        if unknown:
            raise KeyError(f"conservation diagnostics lack declared units: {sorted(unknown)}")
        for name, value in diagnostics.conservation_residuals.items():
            conservation_max[name] = max(conservation_max[name], abs(float(value)))
    cumulative = np.concatenate(
        (np.asarray((0.0,)), cumulative_trapezoid(flow, grid))
    )
    return Trajectory(
        family=str(family),
        initial_condition=str(initial_condition),
        solver=solver,
        success=bool(result.success and complete),
        message=str(result.message),
        time_s=grid,
        state_names=tuple(model.state_names),
        states=states,
        flow_pL_s=flow,
        cumulative_flow_pL=cumulative,
        cell_na_mM=na,
        cell_k_mM=k,
        cell_cl_mM=cl,
        cell_ph=ph,
        cell_volume_pL=volume,
        capacity_multiplier=multiplier,
        max_abs_conservation_residuals=conservation_max,
    )


def fit_wt_flow_scale(
    trajectory: Trajectory,
    envelope: WTFlowEnvelope = WTFlowEnvelope(),
) -> FlowScale:
    """Fit exactly one scale to the permitted WT envelope.

    The statistic is the physical 0--600 s time average.  The implied cell
    count is reported because the scale is otherwise an opaque nuisance
    factor: 1 pL/s per cell equals ``6e-5`` uL/min per cell.
    """

    duration = float(trajectory.time_s[-1] - trajectory.time_s[0])
    if duration <= 0.0:
        raise ValueError("trajectory has no positive duration")
    statistic = float(np.trapezoid(trajectory.flow_pL_s, trajectory.time_s) / duration)
    if not math.isfinite(statistic) or statistic <= 0.0:
        raise ValueError("WT mean model flow must be finite and positive before scaling")
    multiplier = envelope.scale_target_uL_min / statistic
    return FlowScale(
        multiplier_uL_min_per_pL_s=float(multiplier),
        effective_cell_count=float(multiplier / 6.0e-5),
        calibration_statistic="time_mean_0_to_600_s",
        target_uL_min=float(envelope.scale_target_uL_min),
        model_statistic_pL_s=statistic,
        source=envelope.source,
    )


def wt_flow_envelope_feasibility(
    trajectory: Trajectory,
    envelope: WTFlowEnvelope = WTFlowEnvelope(),
    *,
    sampled_times_s: Sequence[float] = tuple(float(60 * minute) for minute in range(1, 11)),
) -> FlowEnvelopeFeasibility:
    """Test whether *any* one scale can fit all WT minute landmarks.

    This is independent of the fitted multiplier and implied cell count.  The
    9--10 uL/min bounds are a graphical approximation, so the result is a
    declared graphical-envelope gate rather than a fabricated pointwise error
    model.
    """

    times = np.asarray(sampled_times_s, dtype=float)
    if times.ndim != 1 or times.size == 0 or np.any(~np.isfinite(times)):
        raise ValueError("flow-envelope sample times must be a finite nonempty vector")
    if np.any(times < trajectory.time_s[0]) or np.any(times > trajectory.time_s[-1]):
        raise ValueError("flow-envelope sample times lie outside the trajectory")
    sampled = np.interp(times, trajectory.time_s, trajectory.flow_pL_s)
    if np.any(~np.isfinite(sampled)) or np.any(sampled <= 0.0):
        raise ValueError("minute-wise WT model flow must be finite and positive")
    scale_lower = float(np.max(envelope.lower_uL_min / sampled))
    scale_upper = float(np.min(envelope.upper_uL_min / sampled))
    minimum = float(np.min(sampled))
    maximum = float(np.max(sampled))
    return FlowEnvelopeFeasibility(
        compatible=bool(scale_lower <= scale_upper),
        sampled_min_flow_pL_s=minimum,
        sampled_max_flow_pL_s=maximum,
        model_max_min_ratio=float(maximum / minimum),
        envelope_max_min_ratio=float(envelope.upper_uL_min / envelope.lower_uL_min),
        feasible_multiplier_lower_uL_min_per_pL_s=scale_lower,
        feasible_multiplier_upper_uL_min_per_pL_s=scale_upper,
        sampled_times_s=tuple(float(value) for value in times),
    )


def solver_comparison(reference: Trajectory, comparison: Trajectory) -> Mapping[str, float | bool]:
    """Return coordinate-wise and dimensionless differences on a shared grid.

    Radau and BDF are distinct algorithms exposed by the same SciPy
    implementation.  This is a dual-method robustness check, not the
    independent reimplementation required of Task 13B Agent I.
    """

    if not np.array_equal(reference.time_s, comparison.time_s):
        raise ValueError("solver comparison requires the identical time grid")
    if reference.state_names != comparison.state_names:
        raise ValueError("solver comparison requires the identical state layout")
    # The coordinate-specific atol/rtol crossover is a native-unit floor:
    # below it, absolute rather than relative solver accuracy controls.  This
    # avoids a unitless numeric floor being reused across fmol, pL, and
    # regulatory fractions.
    reference_floor = reference.solver.atol_vector(reference.state_names) / reference.solver.rtol
    comparison_floor = (
        comparison.solver.atol_vector(comparison.state_names) / comparison.solver.rtol
    )
    coordinate_floor = np.maximum(reference_floor, comparison_floor)[:, None]
    state_scale = np.maximum(
        np.maximum(np.abs(reference.states), np.abs(comparison.states)),
        coordinate_floor,
    )
    state_relative = np.abs(reference.states - comparison.states) / state_scale
    flow_scale = np.maximum(
        np.maximum(np.abs(reference.flow_pL_s), np.abs(comparison.flow_pL_s)), 1.0e-12
    )
    flow_relative = np.abs(reference.flow_pL_s - comparison.flow_pL_s) / flow_scale
    result: dict[str, float | bool] = {
        "both_success": bool(reference.success and comparison.success),
        "max_relative_state_difference": float(np.max(state_relative)),
        "max_abs_flow_difference_pL_s": float(
            np.max(np.abs(reference.flow_pL_s - comparison.flow_pL_s))
        ),
        "max_relative_flow_difference": float(np.max(flow_relative)),
        "endpoint_cumulative_flow_relative_difference": float(
            abs(reference.cumulative_flow_pL[-1] - comparison.cumulative_flow_pL[-1])
            / max(
                abs(reference.cumulative_flow_pL[-1]),
                abs(comparison.cumulative_flow_pL[-1]),
                1.0e-12,
            )
        ),
    }
    for index, name in enumerate(reference.state_names):
        result[f"max_abs_state_difference__{name}"] = float(
            np.max(np.abs(reference.states[index] - comparison.states[index]))
        )
    return result


def regulatory_landmarks(trajectory: Trajectory) -> Mapping[str, float | str | None]:
    """Report, but never fit, regulatory activation landmarks."""

    values = trajectory.capacity_multiplier
    basal = float(values[0])
    final = float(values[-1])
    span = final - basal

    def crossing(fraction: float) -> float | None:
        if span <= 0.0:
            return None
        threshold = basal + fraction * span
        indices = np.flatnonzero(values >= threshold)
        return None if indices.size == 0 else float(trajectory.time_s[int(indices[0])])

    return {
        "family": trajectory.family,
        "basal_capacity_multiplier": basal,
        "endpoint_capacity_multiplier": final,
        "t10_s": crossing(0.10),
        "t50_s": crossing(0.50),
        "t90_s": crossing(0.90),
        "kinetic_fit_status": "UNIDENTIFIED_SENSITIVITY_NOT_FIT",
        "evidence_limit": "P21 contains no acute beta-to-cAMP/PKA/AE4 onset trace",
    }


def _jsonable(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _jsonable(item) for key, item in asdict(value).items()}
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


def canonical_json_bytes(value: Any) -> bytes:
    """Canonical encoding used for frozen pre-reveal object hashes."""

    return json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def sha256_object(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def sha256_file(path: str | Path) -> str:
    """Hash a declared non-held-out artifact."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_trajectory_csv(path: str | Path, trajectories: Iterable[Trajectory]) -> None:
    """Write tidy WT-only trajectories."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        "family",
        "initial_condition",
        "solver_label",
        "solver_method",
        "time_s",
        "flow_pL_s",
        "cumulative_flow_pL",
        "cell_na_mM",
        "cell_k_mM",
        "cell_cl_mM",
        "cell_ph",
        "cell_volume_pL",
        "capacity_multiplier",
    )
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for trajectory in trajectories:
            for index, time_value in enumerate(trajectory.time_s):
                writer.writerow(
                    {
                        "family": trajectory.family,
                        "initial_condition": trajectory.initial_condition,
                        "solver_label": trajectory.solver.label,
                        "solver_method": trajectory.solver.method,
                        "time_s": f"{float(time_value):.12g}",
                        "flow_pL_s": f"{trajectory.flow_pL_s[index]:.12g}",
                        "cumulative_flow_pL": f"{trajectory.cumulative_flow_pL[index]:.12g}",
                        "cell_na_mM": f"{trajectory.cell_na_mM[index]:.12g}",
                        "cell_k_mM": f"{trajectory.cell_k_mM[index]:.12g}",
                        "cell_cl_mM": f"{trajectory.cell_cl_mM[index]:.12g}",
                        "cell_ph": f"{trajectory.cell_ph[index]:.12g}",
                        "cell_volume_pL": f"{trajectory.cell_volume_pL[index]:.12g}",
                        "capacity_multiplier": f"{trajectory.capacity_multiplier[index]:.12g}",
                    }
                )


def assert_pre_reveal_payload(payload: Mapping[str, Any]) -> None:
    """Reject a manifest that names or embeds forbidden validation material."""

    encoded = canonical_json_bytes(payload).decode("utf-8").lower()
    forbidden_tokens = (
        "heldout_targets.csv",
        "ae4_null_secretion",
        "knockout_secretion",
        "ko_secretion",
    )
    found = [token for token in forbidden_tokens if token in encoded]
    if found:
        raise ValueError(f"pre-reveal payload contains forbidden material: {found}")


def write_json(path: str | Path, payload: Mapping[str, Any], *, pre_reveal: bool = False) -> None:
    if pre_reveal:
        assert_pre_reveal_payload(payload)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


__all__ = (
    "CONSERVATION_RESIDUAL_TOLERANCES",
    "CONSERVATION_RESIDUAL_UNITS",
    "LOOSE_RADAU",
    "PRODUCTION_BDF",
    "PRODUCTION_RADAU",
    "TIGHT_RADAU",
    "FlowScale",
    "FlowEnvelopeFeasibility",
    "ArmResponse",
    "SecretagogueProtocol",
    "RegulatoryEnsembleMember",
    "SolverSpecification",
    "StimulusArm",
    "Trajectory",
    "WTFlowEnvelope",
    "assert_pre_reveal_payload",
    "attach_basal_regulation",
    "canonical_json_bytes",
    "fit_wt_flow_scale",
    "nearby_charge_preserving_states",
    "physical_time_grid",
    "pre_reveal_regulatory_ensemble",
    "regulatory_landmarks",
    "sha256_file",
    "sha256_object",
    "simulate_wt",
    "solver_comparison",
    "write_json",
    "write_trajectory_csv",
    "wt_arm_response",
    "wt_flow_envelope_feasibility",
)
