"""Pre-reveal diagnostic for the beta-sensitive, volume-gated Cl pathway.

Catalan et al. (2015, PNAS, doi:10.1073/pnas.1415739112) showed that
isoproterenol can support secretion when acinar TMEM16A is deleted and
reported slow acinar-cell swelling together with sensitivity to anion-channel
blockers.  Those experiments require a beta-associated apical anion-exit
topology distinct from the Ca-activated TMEM16A current.  They do *not*
provide a whole-cell conductance that can be mapped onto this model.

This module therefore implements only the predeclared nested diagnostic:

``V0``
    Exact zero-current control.

``V1``
    ``g_V * beta * max(V_i / V_rest - 1, 0)`` is added to the apical chloride
    conductance.  It has no additional dynamic state.  ``g_V`` is an
    unlicensed diagnostic capacity, so V1 is never production eligible.

The one-sided swelling gate creates an important anti-deadlock requirement.
Starting from an exact resting volume, some *pre-V1* process must first make
the cell swell.  V1 cannot be credited with creating the swelling that is
needed to activate V1.  :func:`assess_pre_vbeta_trajectory` makes that test
explicit and rejects a no-swelling trajectory regardless of the chosen
diagnostic conductance.

Sign convention
---------------
The voltage convention matches :mod:`modern_full_model.membranes`:
``V_a = phi_cell - phi_lumen`` and positive current is conventional current
from cell to lumen.  Chloride molar flux is ``I / (-F)``; cell and lumen
chloride sources are equal and opposite.

Provenance
----------
``PRIMARY_MEASUREMENT``
    A beta-associated, blocker-sensitive anion-exit pathway and IPR-associated
    swelling in adult mouse salivary acini (Catalan et al. 2015).
``DERIVED_CONSTRAINT``
    Nernst reversal, current-to-flux conversion, exact source/charge closure,
    and the exact-rest anti-deadlock test.
``NEW_MODELING_DECISION``
    The rectified volume gate and every numerical value of ``g_V``.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike

from .membranes import current_to_fmol_s, nernst_voltage_V


class VbetaTopology(str, Enum):
    """Nested beta/volume-sensitive apical chloride alternatives."""

    V0_EXACT_CONTROL = "V0_EXACT_CONTROL"
    V1_RECTIFIED_SWELLING_GATE = "V1_RECTIFIED_SWELLING_GATE"


CAPACITY_EVIDENCE_STATUS = "UNLICENSED_NO_WHOLE_CELL_CURRENT_MAP"
V1_PRODUCTION_ELIGIBLE = False


@dataclass(frozen=True)
class VbetaParameters:
    """Parameters for a diagnostic topology, never a fitted production path."""

    topology: VbetaTopology = VbetaTopology.V0_EXACT_CONTROL
    maximum_conductance_S: float = 0.0
    reference_cell_volume_pL: float = 1.3

    def __post_init__(self) -> None:
        object.__setattr__(self, "topology", VbetaTopology(self.topology))
        if (
            not math.isfinite(self.maximum_conductance_S)
            or self.maximum_conductance_S < 0.0
        ):
            raise ValueError("Vbeta diagnostic conductance must be finite and nonnegative")
        if (
            not math.isfinite(self.reference_cell_volume_pL)
            or self.reference_cell_volume_pL <= 0.0
        ):
            raise ValueError("Vbeta reference cell volume must be finite and positive")

    @property
    def production_eligible(self) -> bool:
        """V0 is a control and V1 lacks an absolute-capacity observation map."""

        return False


@dataclass(frozen=True)
class VbetaFlux:
    """One fixed-voltage evaluation of the diagnostic apical Cl pathway."""

    topology: str
    beta_input: float
    relative_swelling: float
    rectified_swelling: float
    effective_conductance_S: float
    chloride_reversal_V: float
    current_A: float
    chloride_cell_source_fmol_s: float
    chloride_lumen_source_fmol_s: float
    cell_charge_source_fmol_s: float
    cell_charge_from_current_fmol_s: float
    capacity_evidence_status: str
    production_eligible: bool

    @property
    def chloride_conservation_residual_fmol_s(self) -> float:
        return self.chloride_cell_source_fmol_s + self.chloride_lumen_source_fmol_s

    @property
    def charge_conversion_residual_fmol_s(self) -> float:
        return self.cell_charge_source_fmol_s - self.cell_charge_from_current_fmol_s


@dataclass(frozen=True)
class VbetaActivationAssessment:
    """Anti-deadlock result for a pre-Vbeta IPR-only volume trajectory."""

    status: str
    exact_rest_initialization: bool
    sample_count: int
    maximum_relative_swelling: float
    minimum_relative_swelling: float
    maximum_rectified_gate: float
    activation_detected: bool
    usable_rescue: bool
    production_eligible: bool
    rejection_reasons: tuple[str, ...]


def _beta(value: float) -> float:
    checked = float(value)
    if not math.isfinite(checked) or not 0.0 <= checked <= 1.0:
        raise ValueError("Vbeta beta input must be finite and lie in [0, 1]")
    return checked


def relative_swelling(cell_volume_pL: float, reference_cell_volume_pL: float) -> float:
    """Return ``V_i / V_rest - 1`` on a positive physical-volume domain."""

    values = (float(cell_volume_pL), float(reference_cell_volume_pL))
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("cell and reference volumes must be finite and positive")
    return values[0] / values[1] - 1.0


def vbeta_effective_conductance_S(
    parameters: VbetaParameters,
    *,
    beta_input: float,
    cell_volume_pL: float,
) -> float:
    """Return the extra apical Cl conductance for V0 or V1.

    V0 returns the floating-point value ``0.0`` exactly, even if a nonzero
    probe conductance was supplied.  This guarantees exact nesting rather
    than merely an asymptotically small perturbation.
    """

    beta = _beta(beta_input)
    swelling = relative_swelling(cell_volume_pL, parameters.reference_cell_volume_pL)
    if parameters.topology is VbetaTopology.V0_EXACT_CONTROL:
        return 0.0
    return parameters.maximum_conductance_S * beta * max(swelling, 0.0)


def evaluate_vbeta_current(
    parameters: VbetaParameters,
    *,
    beta_input: float,
    cell_volume_pL: float,
    v_apical_V: float,
    chloride_cell_mM: float,
    chloride_lumen_mM: float,
    thermal_voltage_V: float,
    faraday_C_mol: float,
) -> VbetaFlux:
    """Evaluate Vbeta at a supplied apical voltage.

    This fixed-voltage function is an isolated diagnostic surface, not a
    production membrane solver.  If a future source supplies ``g_V``, its
    effective conductance must be included *inside* the two-membrane algebraic
    current closure, rather than adding this current after closure.
    """

    beta = _beta(beta_input)
    values = (
        float(v_apical_V),
        float(chloride_cell_mM),
        float(chloride_lumen_mM),
        float(thermal_voltage_V),
        float(faraday_C_mol),
    )
    if not math.isfinite(values[0]):
        raise ValueError("apical voltage must be finite")
    if not all(math.isfinite(value) and value > 0.0 for value in values[1:]):
        raise ValueError("chloride, thermal voltage, and Faraday constant must be positive")

    swelling = relative_swelling(cell_volume_pL, parameters.reference_cell_volume_pL)
    conductance = vbeta_effective_conductance_S(
        parameters,
        beta_input=beta,
        cell_volume_pL=cell_volume_pL,
    )
    reversal = nernst_voltage_V(
        chloride_cell_mM,
        chloride_lumen_mM,
        valence=-1,
        thermal_voltage_V=thermal_voltage_V,
    )
    current = 0.0 if conductance == 0.0 else conductance * (v_apical_V - reversal)
    chloride_cell_to_lumen = current_to_fmol_s(
        current,
        valence=-1,
        faraday_C_mol=faraday_C_mol,
    )
    chloride_cell_source = -chloride_cell_to_lumen
    # Cell conventional-charge source equals minus outward conventional
    # current divided by F.  Units are fmol charge-equivalents/s.
    charge_from_current = -current / faraday_C_mol * 1.0e15
    return VbetaFlux(
        topology=parameters.topology.value,
        beta_input=beta,
        relative_swelling=swelling,
        rectified_swelling=max(swelling, 0.0),
        effective_conductance_S=conductance,
        chloride_reversal_V=reversal,
        current_A=current,
        chloride_cell_source_fmol_s=chloride_cell_source,
        chloride_lumen_source_fmol_s=-chloride_cell_source,
        cell_charge_source_fmol_s=-chloride_cell_source,
        cell_charge_from_current_fmol_s=charge_from_current,
        capacity_evidence_status=CAPACITY_EVIDENCE_STATUS,
        production_eligible=False,
    )


def assess_pre_vbeta_trajectory(
    cell_volume_pL: ArrayLike,
    beta_input: ArrayLike,
    *,
    reference_cell_volume_pL: float,
    exact_rest_tolerance: float = 1.0e-12,
    activation_tolerance: float = 1.0e-12,
) -> VbetaActivationAssessment:
    """Test whether pre-V1 dynamics can open the one-sided swelling gate.

    The inputs must come from an IPR-only trajectory evaluated *without* V1.
    This prevents circular credit: V1 cannot be used to generate its own
    activating perturbation.  A trajectory initialized away from the exact
    resting volume is rejected rather than treated as evidence of activation.
    """

    volumes = np.asarray(cell_volume_pL, dtype=float).reshape(-1)
    beta = np.asarray(beta_input, dtype=float).reshape(-1)
    if volumes.size == 0 or volumes.shape != beta.shape:
        raise ValueError("Vbeta volume and beta trajectories must be nonempty and aligned")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("Vbeta volume trajectory must be finite and positive")
    if np.any(~np.isfinite(beta)) or np.any((beta < 0.0) | (beta > 1.0)):
        raise ValueError("Vbeta beta trajectory must lie in [0, 1]")
    for value, name in (
        (reference_cell_volume_pL, "reference_cell_volume_pL"),
        (exact_rest_tolerance, "exact_rest_tolerance"),
        (activation_tolerance, "activation_tolerance"),
    ):
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ValueError(f"{name} must be finite and positive")

    relative = volumes / float(reference_cell_volume_pL) - 1.0
    rectified = beta * np.maximum(relative, 0.0)
    exact_rest = bool(abs(relative[0]) <= exact_rest_tolerance)
    maximum_gate = float(np.max(rectified))
    activated = bool(maximum_gate > activation_tolerance)
    reasons: list[str] = []
    if not exact_rest:
        reasons.append("TRAJECTORY_NOT_INITIALIZED_AT_EXACT_REST")
    if not activated:
        reasons.append("NO_PRE_VBETA_SWELLING_TO_OPEN_ONE_SIDED_GATE")
    reasons.append(CAPACITY_EVIDENCE_STATUS)

    if not exact_rest:
        status = "REJECTED_INVALID_ANTI_DEADLOCK_INITIALIZATION"
    elif not activated:
        status = "REJECTED_EXACT_REST_ANTI_DEADLOCK"
    else:
        status = "DIAGNOSTIC_GATE_OPENS_BUT_CAPACITY_UNLICENSED"
    return VbetaActivationAssessment(
        status=status,
        exact_rest_initialization=exact_rest,
        sample_count=int(volumes.size),
        maximum_relative_swelling=float(np.max(relative)),
        minimum_relative_swelling=float(np.min(relative)),
        maximum_rectified_gate=maximum_gate,
        activation_detected=activated,
        # Even an activated gate is not a usable production rescue without a
        # source-supported whole-cell conductance/observation map.
        usable_rescue=False,
        production_eligible=V1_PRODUCTION_ELIGIBLE,
        rejection_reasons=tuple(reasons),
    )


__all__: Sequence[str] = (
    "CAPACITY_EVIDENCE_STATUS",
    "V1_PRODUCTION_ELIGIBLE",
    "VbetaActivationAssessment",
    "VbetaFlux",
    "VbetaParameters",
    "VbetaTopology",
    "assess_pre_vbeta_trajectory",
    "evaluate_vbeta_current",
    "relative_swelling",
    "vbeta_effective_conductance_S",
)
