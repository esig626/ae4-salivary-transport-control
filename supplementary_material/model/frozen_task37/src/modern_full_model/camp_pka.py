"""Nested beta/cAMP/PKA/AE4 regulatory families for Task 13B.

The 2021 primary study establishes the direction of the biological chain

``beta adrenergic input -> adenylate cyclase/cAMP -> PKA -> AE4 regulation``

but does not measure an acute cAMP trace, a PKA-activity trace, an AE4
regulatory-state trace, or an activation/deactivation time constant.  The
classes below therefore use normalized fractions and seconds-native *model
assumptions*.  Their kinetic defaults are sensitivity values, not estimates.

All retained families couple regulation to AE4 through one common capacity
multiplier.  No family modifies an individual carrier transition.  This is a
deliberately falsifiable working hypothesis: the available evidence does not
identify a transition-specific (R4) coupling.

Stable duck-typed API
---------------------
Every family exposes ``state_names``, ``initial_state(beta_input)``,
``rhs(t, state, beta_input)``, and ``evaluate(t, state, beta_input)``.
``evaluate`` returns a :class:`RegulatoryOutput` containing at least
``family``, ``activity``, and ``capacity_multiplier``.  Fractions and beta
input lie in ``[0, 1]``; time and every rate are in physical seconds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import math
from typing import ClassVar, Mapping, Sequence


# Inputs and analytic equilibria are validated on the exact physical domain.
# Implicit stiff solvers, however, evaluate their RHS at internal trial states
# that can sit just outside an invariant interval before the accepted step is
# projected back by the numerical method.  A separate, explicit trial band
# avoids turning a harmless Radau stage overshoot into solver failure.  It is
# not used for beta input, parameters, or analytic equilibrium construction.
INTEGRATION_TRIAL_TOLERANCE = 1.0e-4


class AE4Construct(str, Enum):
    """Construct-level qualitative discriminator from the 2021 study.

    ``S173A`` removes the detected PKAc-dependent increment in that construct
    assay.  This encodes dependence of the response on S173; it does **not**
    assert that PKA directly phosphorylates S173.  ``S273A`` retains the
    response.  Approximate basal mutant effects are not made kinetic facts.
    """

    WT = "WT"
    S173A = "S173A"
    S273A = "S273A"


CONSTRUCT_RESPONSE_SCALE: Mapping[AE4Construct, float] = {
    AE4Construct.WT: 1.0,
    AE4Construct.S173A: 0.0,
    AE4Construct.S273A: 1.0,
}


@dataclass(frozen=True)
class RegulatoryGain:
    """Common-capacity observation map shared by R0--R3.

    The default fully activated increment of ``0.25`` records the approximate
    25 percent forskolin statement as one *sensitivity*, not an exact fit and
    not a universal biochemical gain.  The same 2021 figure can be visually
    read as roughly a 0.70 increment; :data:`P21_GAIN_SENSITIVITIES` preserves
    both interpretations.  ``coupling_scale=0`` is useful for the qualitative
    PKA-dependence discriminator, but must not be interpreted as a specific
    quantitative H89 model because H89 is nonselective.
    """

    basal_capacity_multiplier: float = 1.0
    fully_activated_increment: float = 0.25
    coupling_scale: float = 1.0

    def __post_init__(self) -> None:
        values = (
            self.basal_capacity_multiplier,
            self.fully_activated_increment,
            self.coupling_scale,
        )
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            raise ValueError("regulatory capacity parameters must be finite and nonnegative")

    def capacity_multiplier(
        self, activity: float, construct: AE4Construct = AE4Construct.WT
    ) -> float:
        value = _fraction(activity, "regulatory activity")
        response_scale = CONSTRUCT_RESPONSE_SCALE[AE4Construct(construct)]
        multiplier = float(
            self.basal_capacity_multiplier
            + self.fully_activated_increment
            * self.coupling_scale
            * response_scale
            * value
        )
        if not math.isfinite(multiplier) or multiplier < 0.0:
            raise ValueError("computed AE4 capacity multiplier must be finite and nonnegative")
        return multiplier


@dataclass(frozen=True)
class GainSensitivity:
    """One non-fitted reading of the internally discrepant 2021 gain report."""

    label: str
    fully_activated_increment: float
    evidential_status: str


P21_GAIN_SENSITIVITIES: tuple[GainSensitivity, ...] = (
    GainSensitivity(
        "P21_PROSE_APPROXIMATE",
        0.25,
        "SENSITIVITY_ONLY_NOT_POINT_TARGET",
    ),
    GainSensitivity(
        "P21_FIGURE_VISUAL_APPROXIMATE",
        0.70,
        "SENSITIVITY_ONLY_NOT_POINT_TARGET",
    ),
)


@dataclass(frozen=True)
class RegulatoryOutput:
    """Common regulatory output consumed by :class:`ModernFullModel`."""

    family: str
    activity: float
    capacity_multiplier: float
    construct: str
    coupling_hypothesis: str = "COMMON_AE4_CAPACITY_ONLY"
    kinetic_evidence_status: str = "UNMEASURED_ACUTE_KINETICS"


def _time(time_s: float) -> float:
    value = float(time_s)
    if not math.isfinite(value):
        raise ValueError("regulatory time must be finite")
    return value


def _positive_seconds(value: float, name: str) -> float:
    checked = float(value)
    if not math.isfinite(checked) or checked <= 0.0:
        raise ValueError(f"{name} must be finite and positive in seconds")
    return checked


def _beta(beta_input: float) -> float:
    return _fraction(beta_input, "normalized beta input")


def _fraction(value: float, name: str) -> float:
    checked = float(value)
    if not math.isfinite(checked):
        raise ValueError(f"{name} must be finite")
    if checked < 0.0 or checked > 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return checked


def _integration_trial_fraction(value: float, name: str) -> float:
    checked = float(value)
    if not math.isfinite(checked):
        raise ValueError(f"{name} must be finite")
    if (
        checked < -INTEGRATION_TRIAL_TOLERANCE
        or checked > 1.0 + INTEGRATION_TRIAL_TOLERANCE
    ):
        raise ValueError(
            f"{name} lies outside [0, 1] by more than the declared "
            f"implicit-solver trial tolerance {INTEGRATION_TRIAL_TOLERANCE:g}"
        )
    return min(1.0, max(0.0, checked))


def _state(state: Sequence[float], names: Sequence[str]) -> tuple[float, ...]:
    try:
        raw = tuple(float(value) for value in state)
    except TypeError as exc:
        raise ValueError("regulatory state must be a finite sequence") from exc
    if len(raw) != len(names):
        raise ValueError(
            f"regulatory state has length {len(raw)}; expected {len(names)} for {tuple(names)}"
        )
    return tuple(
        _integration_trial_fraction(value, name)
        for value, name in zip(raw, names)
    )


def _output(
    family: str,
    activity: float,
    gain: RegulatoryGain,
    construct: AE4Construct,
) -> RegulatoryOutput:
    bounded_activity = _fraction(activity, "regulatory activity")
    return RegulatoryOutput(
        family=family,
        activity=bounded_activity,
        capacity_multiplier=gain.capacity_multiplier(bounded_activity, construct),
        construct=AE4Construct(construct).value,
    )


@dataclass(frozen=True)
class R0Static:
    """R0: zero-state instantaneous capacity control.

    R0 is the static comparison, not the mandatory dynamic biological model.
    At constant beta input it is exactly the equilibrium output map of R1 and
    of the default symmetric R2/R3 families.
    """

    gain: RegulatoryGain = field(default_factory=RegulatoryGain)
    construct: AE4Construct = AE4Construct.WT

    family: ClassVar[str] = "R0"
    state_names: ClassVar[tuple[str, ...]] = ()

    def __post_init__(self) -> None:
        AE4Construct(self.construct)

    def equilibrium(self, beta_input: float) -> tuple[float, ...]:
        _beta(beta_input)
        return ()

    def initial_state(self, beta_input: float) -> tuple[float, ...]:
        return self.equilibrium(beta_input)

    def rhs(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> tuple[float, ...]:
        _time(time_s)
        _state(state, self.state_names)
        _beta(beta_input)
        return ()

    def evaluate(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> RegulatoryOutput:
        _time(time_s)
        _state(state, self.state_names)
        return _output(self.family, _beta(beta_input), self.gain, self.construct)


@dataclass(frozen=True)
class R1EffectiveActivation:
    """R1: one effective beta/PKA/AE4 activation fraction.

    ``dx/dt = (beta - x) / tau_activation``.
    The state deliberately aggregates upstream cAMP, PKA, and the AE4
    regulatory response because the 2021 endpoint assays do not separate
    their acute kinetics.
    """

    tau_activation_s: float = 30.0
    gain: RegulatoryGain = field(default_factory=RegulatoryGain)
    construct: AE4Construct = AE4Construct.WT

    family: ClassVar[str] = "R1"
    state_names: ClassVar[tuple[str, ...]] = ("ae4_effective_activation_fraction",)

    def __post_init__(self) -> None:
        _positive_seconds(self.tau_activation_s, "tau_activation_s")
        AE4Construct(self.construct)

    def equilibrium(self, beta_input: float) -> tuple[float, ...]:
        return (_beta(beta_input),)

    def initial_state(self, beta_input: float) -> tuple[float, ...]:
        return self.equilibrium(beta_input)

    def rhs(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> tuple[float, ...]:
        _time(time_s)
        (activation,) = _state(state, self.state_names)
        target = _beta(beta_input)
        return ((target - activation) / self.tau_activation_s,)

    def evaluate(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> RegulatoryOutput:
        _time(time_s)
        (activation,) = _state(state, self.state_names)
        _beta(beta_input)
        return _output(self.family, activation, self.gain, self.construct)


@dataclass(frozen=True)
class R2CampActivation:
    """R2: normalized effective cAMP followed by AE4 activation.

    ``dc/dt = (beta - c) / tau_camp`` and
    ``dx/dt = (c - x) / tau_activation``.

    ``c`` is not a measured cAMP concentration.  It is a bounded effective
    adenylate-cyclase/cAMP coordinate necessitated only by the pathway-order
    hypothesis.  Under constant beta, the manifold ``c=beta`` is invariant;
    on that manifold R2 collapses exactly to R1 with the same downstream time
    constant.
    """

    tau_camp_s: float = 10.0
    tau_activation_s: float = 30.0
    gain: RegulatoryGain = field(default_factory=RegulatoryGain)
    construct: AE4Construct = AE4Construct.WT

    family: ClassVar[str] = "R2"
    state_names: ClassVar[tuple[str, ...]] = (
        "camp_effective_fraction",
        "ae4_activation_fraction",
    )

    def __post_init__(self) -> None:
        _positive_seconds(self.tau_camp_s, "tau_camp_s")
        _positive_seconds(self.tau_activation_s, "tau_activation_s")
        AE4Construct(self.construct)

    def equilibrium(self, beta_input: float) -> tuple[float, ...]:
        beta = _beta(beta_input)
        return beta, beta

    def initial_state(self, beta_input: float) -> tuple[float, ...]:
        return self.equilibrium(beta_input)

    def rhs(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> tuple[float, ...]:
        _time(time_s)
        camp, activation = _state(state, self.state_names)
        target = _beta(beta_input)
        return (
            (target - camp) / self.tau_camp_s,
            (camp - activation) / self.tau_activation_s,
        )

    def evaluate(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> RegulatoryOutput:
        _time(time_s)
        _, activation = _state(state, self.state_names)
        _beta(beta_input)
        return _output(self.family, activation, self.gain, self.construct)


@dataclass(frozen=True)
class R3PkaRegulatedFraction:
    """R3: upstream PKA activity plus an AE4 regulated fraction.

    The equations are

    ``dp/dt = (beta - p) / tau_pka``

    ``dq/dt = k_forward*p*(1-q) - k_reverse*(1-p)*q``.

    ``q`` is a phosphorylation/dephosphorylation-*like* regulated fraction.
    The labels ``forward`` and ``reverse`` are intentional: direct S173
    phosphorylation was not demonstrated.  When
    ``k_forward == k_reverse == 1/tau_activation``, the second equation is
    exactly ``(p-q)/tau_activation``.  R3 is then dynamically isomorphic to R2
    under ``p <-> c`` and ``q <-> x`` for any beta waveform.
    """

    tau_pka_s: float = 10.0
    forward_regulation_rate_s: float = 1.0 / 30.0
    reverse_regulation_rate_s: float = 1.0 / 30.0
    gain: RegulatoryGain = field(default_factory=RegulatoryGain)
    construct: AE4Construct = AE4Construct.WT

    family: ClassVar[str] = "R3"
    state_names: ClassVar[tuple[str, ...]] = (
        "pka_effective_activity_fraction",
        "ae4_regulated_fraction",
    )

    def __post_init__(self) -> None:
        _positive_seconds(self.tau_pka_s, "tau_pka_s")
        for value, name in (
            (self.forward_regulation_rate_s, "forward_regulation_rate_s"),
            (self.reverse_regulation_rate_s, "reverse_regulation_rate_s"),
        ):
            checked = float(value)
            if not math.isfinite(checked) or checked <= 0.0:
                raise ValueError(f"{name} must be finite and positive in s^-1")
        AE4Construct(self.construct)

    def equilibrium(self, beta_input: float) -> tuple[float, ...]:
        pka = _beta(beta_input)
        forward = self.forward_regulation_rate_s * pka
        reverse = self.reverse_regulation_rate_s * (1.0 - pka)
        denominator = forward + reverse
        if not math.isfinite(denominator) or denominator <= 0.0:
            raise ValueError("R3 equilibrium rates do not define a finite positive denominator")
        regulated = forward / denominator
        return pka, regulated

    def initial_state(self, beta_input: float) -> tuple[float, ...]:
        return self.equilibrium(beta_input)

    def rhs(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> tuple[float, ...]:
        _time(time_s)
        pka, regulated = _state(state, self.state_names)
        target = _beta(beta_input)
        return (
            (target - pka) / self.tau_pka_s,
            self.forward_regulation_rate_s * pka * (1.0 - regulated)
            - self.reverse_regulation_rate_s * (1.0 - pka) * regulated,
        )

    def evaluate(
        self, time_s: float, state: Sequence[float], beta_input: float
    ) -> RegulatoryOutput:
        _time(time_s)
        _, regulated = _state(state, self.state_names)
        _beta(beta_input)
        return _output(self.family, regulated, self.gain, self.construct)


R4_REJECTION: Mapping[str, str] = {
    "family": "R4",
    "status": "REJECTED_NOT_RETAINED",
    "reason": (
        "The 2021 experiments do not distinguish common AE4 capacity modulation "
        "from modulation of any selected carrier transition, so a state-specific "
        "transition model would add unsupported freedom."
    ),
}


def default_regulatory_families(
    *,
    gain: RegulatoryGain | None = None,
    construct: AE4Construct = AE4Construct.WT,
) -> Mapping[str, object]:
    """Return a parameter-nested R0--R3 registry for WT sensitivity runs.

    Defaults make R2 and R3 exactly isomorphic: ``tau_camp=tau_pka=10 s``
    and both downstream rates equal ``1/30 s^-1``.  These values are declared
    kinetic sensitivities only; no 2021 acute trace identifies them.
    """

    common_gain = RegulatoryGain() if gain is None else gain
    selected_construct = AE4Construct(construct)
    return {
        "R0": R0Static(gain=common_gain, construct=selected_construct),
        "R1": R1EffectiveActivation(
            tau_activation_s=30.0,
            gain=common_gain,
            construct=selected_construct,
        ),
        "R2": R2CampActivation(
            tau_camp_s=10.0,
            tau_activation_s=30.0,
            gain=common_gain,
            construct=selected_construct,
        ),
        "R3": R3PkaRegulatedFraction(
            tau_pka_s=10.0,
            forward_regulation_rate_s=1.0 / 30.0,
            reverse_regulation_rate_s=1.0 / 30.0,
            gain=common_gain,
            construct=selected_construct,
        ),
    }


# Short aliases are convenient in generation registries while the descriptive
# class names remain explicit in frozen manifests.
R0 = R0Static
R1 = R1EffectiveActivation
R2 = R2CampActivation
R3 = R3PkaRegulatedFraction


__all__ = (
    "AE4Construct",
    "CONSTRUCT_RESPONSE_SCALE",
    "GainSensitivity",
    "INTEGRATION_TRIAL_TOLERANCE",
    "P21_GAIN_SENSITIVITIES",
    "R0",
    "R0Static",
    "R1",
    "R1EffectiveActivation",
    "R2",
    "R2CampActivation",
    "R3",
    "R3PkaRegulatedFraction",
    "R4_REJECTION",
    "RegulatoryGain",
    "RegulatoryOutput",
    "default_regulatory_families",
)
