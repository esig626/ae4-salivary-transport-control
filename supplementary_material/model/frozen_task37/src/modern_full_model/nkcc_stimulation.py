"""WT-only, source-disciplined stimulation of basolateral NKCC1.

The conservation-explicit core deliberately started with a reversible NKCC1
law whose capacity was independent of secretagogue input.  Peña-Münzenmayer
et al. (2015; DOI 10.1074/jbc.M114.612895) instead show that the native WT
muscarinic arm has a large chloride-reuptake response and that 50 uM
bumetanide removes 95.6 +/- 0.9% of the experimentally isolated NKCC1 signal.
Those observations justify testing a *capacity* response while leaving the
measured 1 Na : 1 K : 2 Cl stoichiometry and thermodynamic reversal intact.

This module provides a strict nested sequence:

``N0_BASELINE``
    The original stimulus-independent capacity (exact multiplier one).
``N1_CCH_ALGEBRAIC``
    One algebraic, calcium-arm capacity gain.  It adds one gain but no time
    constant.
``N2_CCH_EFFECTIVE_DYNAMIC``
    One effective activation fraction, ``dn/dt=(u_Ca-n)/tau_N``, followed by
    the same capacity map.  Its time constant is unmeasured and must be
    profiled as a sensitivity, never fitted to an AE4-null secretion trace.

Both extensions recover the frozen resting model exactly because the calcium
activation coordinate is zero at the declared 0.058 uM rest input.  Setting
their fully activated multiplier to one recovers N0 for every input.  The
wrapper changes only the NKCC1 capacity passed to the existing reversible
homeostasis evaluator; it does not change affinity, reversal, stoichiometry,
other transporters, water, or cAMP/PKA regulation of AE4.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any, ClassVar, Mapping, Protocol, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from .membranes import (
    HomeostasisEnvironment,
    evaluate_homeostasis,
)
from .model import Genotype, ModelEvaluation, ModernFullModel, StimulusProtocol, WT


RESTING_CALCIUM_UM = 0.058
STIMULATED_CALCIUM_UM = 0.55
_INTEGRATION_TRIAL_TOLERANCE = 1.0e-4


def _finite_time(time_s: float) -> float:
    value = float(time_s)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("time must be finite and nonnegative")
    return value


def _finite_calcium(calcium_uM: float) -> float:
    value = float(calcium_uM)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("calcium input must be finite and nonnegative")
    return value


def _state_fraction(value: float, name: str) -> float:
    checked = float(value)
    if not math.isfinite(checked):
        raise ValueError(f"{name} must be finite")
    if checked < -_INTEGRATION_TRIAL_TOLERANCE or checked > 1.0 + _INTEGRATION_TRIAL_TOLERANCE:
        raise ValueError(f"{name} must remain in [0, 1]")
    # Radau and BDF can evaluate an implicit trial stage just beyond a bounded
    # state even when the accepted solution remains physical.
    return min(1.0, max(0.0, checked))


def normalized_calcium_arm(
    calcium_uM: float,
    *,
    resting_calcium_uM: float = RESTING_CALCIUM_UM,
    stimulated_calcium_uM: float = STIMULATED_CALCIUM_UM,
) -> float:
    """Map the declared CCh/Ca sensitivity endpoints to an exact [0, 1] arm.

    This is a protocol-normalized input, not a measured NKCC1 dose-response.
    It avoids adding a half-saturation and Hill coefficient that P15 does not
    identify.  Values outside the two declared calcium endpoints saturate.
    """

    calcium = _finite_calcium(calcium_uM)
    resting = _finite_calcium(resting_calcium_uM)
    stimulated = _finite_calcium(stimulated_calcium_uM)
    if stimulated <= resting:
        raise ValueError("stimulated calcium must exceed resting calcium")
    return min(1.0, max(0.0, (calcium - resting) / (stimulated - resting)))


@dataclass(frozen=True)
class Nkcc1RegulatoryOutput:
    family: str
    activation_fraction: float
    capacity_multiplier: float
    evidence_status: str


class Nkcc1RegulatoryModel(Protocol):
    state_names: tuple[str, ...]

    def initial_state(self, calcium_uM: float) -> tuple[float, ...]: ...

    def rhs(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> tuple[float, ...]: ...

    def evaluate(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> Nkcc1RegulatoryOutput: ...


@dataclass(frozen=True)
class N0BaselineNkcc1:
    """Exact stimulus-independent control already present in the core."""

    family: ClassVar[str] = "N0_BASELINE"
    state_names: ClassVar[tuple[str, ...]] = ()

    def initial_state(self, calcium_uM: float) -> tuple[float, ...]:
        _finite_calcium(calcium_uM)
        return ()

    def rhs(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> tuple[float, ...]:
        _finite_time(time_s)
        _finite_calcium(calcium_uM)
        if tuple(state):
            raise ValueError("N0 has no state")
        return ()

    def evaluate(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> Nkcc1RegulatoryOutput:
        self.rhs(time_s, state, calcium_uM)
        return Nkcc1RegulatoryOutput(
            family=self.family,
            activation_fraction=0.0,
            capacity_multiplier=1.0,
            evidence_status="EXACT_BASELINE_CONTROL",
        )


@dataclass(frozen=True)
class N1AlgebraicNkcc1:
    """Smallest one-gain CCh/Ca capacity hypothesis."""

    fully_activated_multiplier: float = 1.75
    resting_calcium_uM: float = RESTING_CALCIUM_UM
    stimulated_calcium_uM: float = STIMULATED_CALCIUM_UM

    family: ClassVar[str] = "N1_CCH_ALGEBRAIC"
    state_names: ClassVar[tuple[str, ...]] = ()

    def __post_init__(self) -> None:
        multiplier = float(self.fully_activated_multiplier)
        if not math.isfinite(multiplier) or multiplier < 1.0:
            raise ValueError("fully activated NKCC1 multiplier must be finite and >= 1")
        normalized_calcium_arm(
            self.resting_calcium_uM,
            resting_calcium_uM=self.resting_calcium_uM,
            stimulated_calcium_uM=self.stimulated_calcium_uM,
        )

    def _activation(self, calcium_uM: float) -> float:
        return normalized_calcium_arm(
            calcium_uM,
            resting_calcium_uM=self.resting_calcium_uM,
            stimulated_calcium_uM=self.stimulated_calcium_uM,
        )

    def initial_state(self, calcium_uM: float) -> tuple[float, ...]:
        self._activation(calcium_uM)
        return ()

    def rhs(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> tuple[float, ...]:
        _finite_time(time_s)
        self._activation(calcium_uM)
        if tuple(state):
            raise ValueError("N1 has no state")
        return ()

    def evaluate(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> Nkcc1RegulatoryOutput:
        self.rhs(time_s, state, calcium_uM)
        activation = self._activation(calcium_uM)
        multiplier = 1.0 + (self.fully_activated_multiplier - 1.0) * activation
        return Nkcc1RegulatoryOutput(
            family=self.family,
            activation_fraction=activation,
            capacity_multiplier=multiplier,
            evidence_status=(
                "P15_MUSCARINIC_AND_BUMETANIDE_DIRECTION; GAIN_PROFILED_WT_ONLY"
            ),
        )


@dataclass(frozen=True)
class N2EffectiveDynamicNkcc1:
    """One effective activation state for the CCh/Ca-to-NKCC1 response."""

    fully_activated_multiplier: float = 1.75
    tau_activation_s: float = 60.0
    resting_calcium_uM: float = RESTING_CALCIUM_UM
    stimulated_calcium_uM: float = STIMULATED_CALCIUM_UM

    family: ClassVar[str] = "N2_CCH_EFFECTIVE_DYNAMIC"
    state_names: ClassVar[tuple[str, ...]] = ("nkcc1_effective_activation_fraction",)

    def __post_init__(self) -> None:
        multiplier = float(self.fully_activated_multiplier)
        tau = float(self.tau_activation_s)
        if not math.isfinite(multiplier) or multiplier < 1.0:
            raise ValueError("fully activated NKCC1 multiplier must be finite and >= 1")
        if not math.isfinite(tau) or tau <= 0.0:
            raise ValueError("NKCC1 activation time must be finite and positive")
        normalized_calcium_arm(
            self.resting_calcium_uM,
            resting_calcium_uM=self.resting_calcium_uM,
            stimulated_calcium_uM=self.stimulated_calcium_uM,
        )

    def _target(self, calcium_uM: float) -> float:
        return normalized_calcium_arm(
            calcium_uM,
            resting_calcium_uM=self.resting_calcium_uM,
            stimulated_calcium_uM=self.stimulated_calcium_uM,
        )

    def initial_state(self, calcium_uM: float) -> tuple[float, ...]:
        return (self._target(calcium_uM),)

    def rhs(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> tuple[float, ...]:
        _finite_time(time_s)
        if len(tuple(state)) != 1:
            raise ValueError("N2 requires one activation fraction")
        activation = _state_fraction(tuple(state)[0], self.state_names[0])
        return ((self._target(calcium_uM) - activation) / self.tau_activation_s,)

    def evaluate(
        self, time_s: float, state: Sequence[float], calcium_uM: float
    ) -> Nkcc1RegulatoryOutput:
        self.rhs(time_s, state, calcium_uM)
        activation = _state_fraction(tuple(state)[0], self.state_names[0])
        multiplier = 1.0 + (self.fully_activated_multiplier - 1.0) * activation
        return Nkcc1RegulatoryOutput(
            family=self.family,
            activation_fraction=activation,
            capacity_multiplier=multiplier,
            evidence_status=(
                "P15_MUSCARINIC_AND_BUMETANIDE_DIRECTION; TAU_UNMEASURED_WT_SENSITIVITY"
            ),
        )


def _object_mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return dict(value)
    if hasattr(value, "__dataclass_fields__"):
        return {name: getattr(value, name) for name in value.__dataclass_fields__}
    return {name: getattr(value, name) for name in dir(value) if not name.startswith("_")}


@dataclass(frozen=True)
class CompositeAe4NkccRegulation:
    """Concatenate mandatory AE4 cAMP/PKA states with an NKCC1 hypothesis."""

    ae4_regulatory_model: Any
    nkcc1_regulatory_model: Nkcc1RegulatoryModel
    stimulus: StimulusProtocol

    @property
    def family(self) -> str:
        ae4_family = getattr(self.ae4_regulatory_model, "family", "AE4_REGULATION")
        nkcc_family = getattr(self.nkcc1_regulatory_model, "family", "NKCC1_REGULATION")
        return f"{ae4_family}+{nkcc_family}"

    @property
    def state_names(self) -> tuple[str, ...]:
        ae4_names = tuple(str(name) for name in self.ae4_regulatory_model.state_names)
        nkcc_names = tuple(str(name) for name in self.nkcc1_regulatory_model.state_names)
        if set(ae4_names) & set(nkcc_names):
            raise ValueError("AE4 and NKCC1 regulatory state names must be disjoint")
        return ae4_names + nkcc_names

    @property
    def _ae4_size(self) -> int:
        return len(tuple(self.ae4_regulatory_model.state_names))

    def _split(self, state: Sequence[float]) -> tuple[tuple[float, ...], tuple[float, ...]]:
        values = tuple(float(value) for value in state)
        if len(values) != len(self.state_names):
            raise ValueError("composite regulatory state has the wrong length")
        return values[: self._ae4_size], values[self._ae4_size :]

    def initial_state(self, beta_input: float) -> tuple[float, ...]:
        ae4 = tuple(
            float(value)
            for value in self.ae4_regulatory_model.initial_state(beta_input=beta_input)
        )
        calcium = self.stimulus(0.0).calcium_uM
        nkcc = self.nkcc1_regulatory_model.initial_state(calcium)
        return ae4 + tuple(float(value) for value in nkcc)

    def equilibrium(self, beta_input: float) -> tuple[float, ...]:
        equilibrium = getattr(
            self.ae4_regulatory_model,
            "equilibrium",
            self.ae4_regulatory_model.initial_state,
        )
        ae4 = tuple(float(value) for value in equilibrium(beta_input=beta_input))
        calcium = self.stimulus(0.0).calcium_uM
        nkcc = self.nkcc1_regulatory_model.initial_state(calcium)
        return ae4 + tuple(float(value) for value in nkcc)

    def rhs(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> tuple[float, ...]:
        ae4_state, nkcc_state = self._split(state)
        ae4_rhs = self.ae4_regulatory_model.rhs(
            time_s, ae4_state, beta_input=beta_input
        )
        calcium = self.stimulus(time_s).calcium_uM
        nkcc_rhs = self.nkcc1_regulatory_model.rhs(time_s, nkcc_state, calcium)
        return tuple(float(value) for value in ae4_rhs) + tuple(
            float(value) for value in nkcc_rhs
        )

    def evaluate(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> Mapping[str, Any]:
        ae4_state, nkcc_state = self._split(state)
        ae4 = _object_mapping(
            self.ae4_regulatory_model.evaluate(
                time_s, ae4_state, beta_input=beta_input
            )
        )
        calcium = self.stimulus(time_s).calcium_uM
        nkcc = self.nkcc1_regulatory_model.evaluate(time_s, nkcc_state, calcium)
        return {
            **ae4,
            "family": self.family,
            "ae4_regulatory_family": ae4.get("family"),
            "nkcc1_regulatory_family": nkcc.family,
            "nkcc1_activation_fraction": nkcc.activation_fraction,
            "nkcc1_capacity_multiplier": nkcc.capacity_multiplier,
            "nkcc1_evidence_status": nkcc.evidence_status,
        }


class StimulatedNkcc1Model:
    """Conserved wrapper that applies a composite regulator to NKCC1 capacity."""

    def __init__(self, base_model: ModernFullModel) -> None:
        self.base_model = base_model
        initial = self.base_model.evaluate(0.0, self.base_model.initial_state(), genotype=WT)
        if "nkcc1_capacity_multiplier" not in initial.diagnostics.regulatory:
            raise ValueError(
                "base model must use CompositeAe4NkccRegulation before NKCC1 wrapping"
            )

    @property
    def parameters(self) -> Any:
        return self.base_model.parameters

    @property
    def stimulus(self) -> Any:
        return self.base_model.stimulus

    @property
    def regulatory_model(self) -> Any:
        return self.base_model.regulatory_model

    @property
    def ae4_parameters(self) -> Any:
        return self.base_model.ae4_parameters

    @property
    def ae4_evaluator(self) -> Any:
        return self.base_model.ae4_evaluator

    @property
    def layout(self) -> Any:
        return self.base_model.layout

    @property
    def state_names(self) -> tuple[str, ...]:
        return self.base_model.state_names

    def initial_state(self) -> NDArray[np.float64]:
        return self.base_model.initial_state()

    def _regulated_homeostasis(
        self, evaluation: ModelEvaluation, genotype: Genotype
    ) -> Any:
        observables = evaluation.diagnostics.observables
        ci = observables.cell_concentrations_mM
        bath = self.parameters.bath
        bath_ab = observables.bath_acid_base
        multiplier = float(
            evaluation.diagnostics.regulatory["nkcc1_capacity_multiplier"]
        )
        if not math.isfinite(multiplier) or multiplier < 0.0:
            raise ValueError("NKCC1 capacity multiplier must be finite and nonnegative")
        return evaluate_homeostasis(
            HomeostasisEnvironment(
                na_i_mM=ci["na"],
                k_i_mM=ci["k"],
                cl_i_mM=ci["cl"],
                h_i_mM=observables.cell_acid_base.h_mM,
                hco3_i_mM=ci["hco3"],
                na_e_mM=bath.na_mM,
                k_e_mM=bath.k_mM,
                cl_e_mM=bath.cl_mM,
                h_e_mM=bath_ab.h_mM,
                hco3_e_mM=bath_ab.hco3_mM,
            ),
            self.parameters,
            nkcc1_scale=genotype.nkcc1_expression * multiplier,
            nhe1_scale=genotype.nhe1_expression,
            ae2_scale=genotype.ae2_expression,
        )

    def evaluate(
        self, time_s: float, vector: ArrayLike, *, genotype: Genotype = WT
    ) -> ModelEvaluation:
        baseline = self.base_model.evaluate(time_s, vector, genotype=genotype)
        regulated = self._regulated_homeostasis(baseline, genotype)
        previous = baseline.diagnostics.homeostasis
        raw = np.asarray(baseline.rhs, dtype=float).copy()
        for index, field_name in (
            (0, "na_cell_fmol_s"),
            (1, "k_cell_fmol_s"),
            (2, "cl_cell_fmol_s"),
            (3, "tic_cell_fmol_s"),
            (4, "alkalinity_cell_fmol_s"),
        ):
            raw[index] += float(getattr(regulated, field_name)) - float(
                getattr(previous, field_name)
            )

        # NKCC1 adds 1 Na + 1 K + 2 Cl and therefore changes no bulk-charge,
        # carbon, buffer, water, or current residual.  Still recompute the cell
        # bulk-charge diagnostic from the returned source vector rather than
        # assuming the cancellation in an audit-facing artifact.
        conservation = dict(baseline.diagnostics.conservation_residuals)
        conservation["cell_bulk_charge_rate_fmol_s"] = float(
            raw[0] + raw[1] - raw[2] - raw[4]
        )
        conservation["homeostasis_charge_fmol_s"] = float(
            regulated.charge_source_residual_fmol_s
        )
        diagnostics = replace(
            baseline.diagnostics,
            homeostasis=regulated,
            conservation_residuals=conservation,
        )
        return replace(baseline, rhs=raw, diagnostics=diagnostics)

    def rhs(
        self, time_s: float, vector: ArrayLike, *, genotype: Genotype = WT
    ) -> NDArray[np.float64]:
        return self.evaluate(time_s, vector, genotype=genotype).rhs

    def solve_dynamics(
        self,
        time_span_s: tuple[float, float],
        *,
        initial_state: ArrayLike | None = None,
        genotype: Genotype = WT,
        method: str = "Radau",
        rtol: float = 1.0e-7,
        atol: float | ArrayLike = 1.0e-10,
        t_eval: ArrayLike | None = None,
        max_step_s: float = math.inf,
    ) -> Any:
        if method not in {"Radau", "BDF"}:
            raise ValueError("production dynamics must use method='Radau' or 'BDF'")
        y0 = self.initial_state() if initial_state is None else self.layout.validate(initial_state)
        return solve_ivp(
            lambda t, y: self.rhs(t, y, genotype=genotype),
            time_span_s,
            y0,
            method=method,
            rtol=rtol,
            atol=atol,
            t_eval=None if t_eval is None else np.asarray(t_eval, dtype=float),
            max_step=max_step_s,
        )


def attach_stimulated_nkcc1(
    template: ModernFullModel,
    nkcc1_regulatory_model: Nkcc1RegulatoryModel,
) -> StimulatedNkcc1Model:
    """Attach an NKCC1 hypothesis without changing frozen transport parameters."""

    if template.regulatory_model is None:
        raise ValueError("mandatory explicit AE4 beta/cAMP/PKA regulation is required")
    composite = CompositeAe4NkccRegulation(
        ae4_regulatory_model=template.regulatory_model,
        nkcc1_regulatory_model=nkcc1_regulatory_model,
        stimulus=template.stimulus,
    )
    base = ModernFullModel(
        parameters=template.parameters,
        stimulus=template.stimulus,
        regulatory_model=composite,
        ae4_parameters=template.ae4_parameters,
        ae4_evaluator=template.ae4_evaluator,
    )
    return StimulatedNkcc1Model(base)


def isolated_nkcc1_fractional_cl_uptake_s(
    model: StimulatedNkcc1Model,
    time_s: float,
    state: ArrayLike,
    *,
    genotype: Genotype = WT,
) -> float:
    """Return ``2 J_NKCC1 / intracellular Cl amount`` as an audit diagnostic.

    This is *not* silently equated with the P15 SPQ whole-cell uptake slope.
    The assay used a low/high-chloride depletion protocol and reports a net
    fluorescence initial slope.  The returned mechanistic rate is therefore
    useful for scale/profile rejection, but a point fit requires an explicit
    depletion-state and fluorescence observation model.
    """

    vector = model.layout.validate(state)
    evaluation = model.evaluate(time_s, vector, genotype=genotype)
    cell = model.layout.decode(vector).cell
    if cell.cl_fmol <= 0.0:
        raise ValueError("intracellular chloride amount must be positive")
    return float(2.0 * evaluation.diagnostics.homeostasis.nkcc1_inward_fmol_s / cell.cl_fmol)


__all__ = (
    "RESTING_CALCIUM_UM",
    "STIMULATED_CALCIUM_UM",
    "CompositeAe4NkccRegulation",
    "N0BaselineNkcc1",
    "N1AlgebraicNkcc1",
    "N2EffectiveDynamicNkcc1",
    "Nkcc1RegulatoryOutput",
    "StimulatedNkcc1Model",
    "attach_stimulated_nkcc1",
    "isolated_nkcc1_fractional_cl_uptake_s",
    "normalized_calcium_arm",
)
