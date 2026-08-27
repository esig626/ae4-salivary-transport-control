"""Fixed-seven-state-chassis adapter for explicit AE4 CTMC cycles."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

from .chemistry import CarbonateChemistry, DEFAULT_CARBONATE_CHEMISTRY
from .cycles import (
    CycleDefinition,
    CycleEvaluation,
    CycleParameters,
    ReservoirActivities,
    evaluate_cycle,
)


class UnsupportedChassisChemistry(RuntimeError):
    """Raised when a carbonate cycle is inserted into a bicarbonate-only chassis."""


@dataclass(frozen=True)
class StateResolvedBalance:
    """Duck-typed balance consumed by ``FixedChassis``."""

    cycle_flux: float
    na_i: float
    k_i: float
    cl_i: float
    hco3_i: float
    transported_charge: float
    affinity: float
    valid: bool = True
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class StateResolvedAE4Adapter:
    """Quasi-steady CTMC source boundary for the Task-12 fixed chassis.

    The transporter states are solved at stationary occupancy for every
    whole-cell RHS call.  This preserves state competition and edge-specific
    regulation but assumes transporter relaxation is fast.  The absolute
    capacity stays in the chassis-native, physically uncertified flux unit.
    """

    definition: CycleDefinition
    parameters: CycleParameters

    def evaluate(self, environment: Any, *, activity_scale: float) -> StateResolvedBalance:
        if activity_scale < 0.0 or not math.isfinite(activity_scale):
            raise ValueError("activity_scale must be finite and nonnegative")
        # Exact short circuit is important: all explicit AE4 state families
        # must reduce to the same non-AE4 vector field after gene deletion.
        if activity_scale == 0.0 or self.parameters.capacity == 0.0:
            return StateResolvedBalance(
                cycle_flux=0.0,
                na_i=0.0,
                k_i=0.0,
                cl_i=0.0,
                hco3_i=0.0,
                transported_charge=0.0,
                affinity=math.nan,
                diagnostics={
                    "mechanism_id": self.definition.model_id,
                    "state_resolved": True,
                    "knockout_short_circuit": True,
                },
            )
        if "co3_i" in self.definition.required_reservoir_species:
            raise UnsupportedChassisChemistry(
                f"{self.definition.model_id} carries explicit carbonate, but the fixed "
                "seven-state chassis has no carbonate/proton state"
            )
        reservoirs = ReservoirActivities(
            na_i=float(environment.na_i),
            k_i=float(environment.k_i),
            cl_i=float(environment.cl_i),
            hco3_i=float(environment.hco3_i),
            na_e=float(environment.na_e),
            k_e=float(environment.k_e),
            cl_e=float(environment.cl_e),
            hco3_e=float(environment.hco3_e),
        )
        evaluation = evaluate_cycle(
            self.definition,
            reservoirs,
            self.parameters,
            pka_activation=float(getattr(environment, "pka_activation", 0.0)),
        )
        scale = float(activity_scale)
        sources = evaluation.intracellular_sources
        branch_affinities = tuple(evaluation.branch_affinities.values())
        # A multicycle shared network has no meaningful arithmetic-mean
        # affinity.  Keep the scalar only for a true single-cycle adapter;
        # branch affinities and total entropy production remain in diagnostics.
        affinity = float(branch_affinities[0]) if len(branch_affinities) == 1 else math.nan
        return StateResolvedBalance(
            cycle_flux=scale * evaluation.cycle_flux,
            na_i=scale * sources["na_i"],
            k_i=scale * sources["k_i"],
            cl_i=scale * sources["cl_i"],
            hco3_i=scale * sources["hco3_i"],
            transported_charge=scale * evaluation.transported_charge_rate,
            affinity=affinity,
            valid=(
                evaluation.local_detailed_balance_residual <= 1e-10
                and abs(evaluation.transported_charge_rate) <= 1e-10
            ),
            diagnostics=_diagnostics(evaluation),
        )


def _diagnostics(evaluation: CycleEvaluation) -> Mapping[str, Any]:
    return {
        "mechanism_id": evaluation.model_id,
        "state_resolved": True,
        "occupancy": dict(evaluation.occupancy),
        "edge_currents": dict(evaluation.edge_currents),
        "branch_currents": dict(evaluation.branch_currents),
        "branch_affinities": dict(evaluation.branch_affinities),
        "entropy_production": evaluation.entropy_production,
        "generator_residual": evaluation.generator_residual,
        "local_detailed_balance_residual": evaluation.local_detailed_balance_residual,
        "pka_activation": evaluation.pka_activation,
    }


def transporter_reservoirs_with_carbonate(
    *,
    na_i: float,
    k_i: float,
    cl_i: float,
    hco3_i: float,
    ph_i: float,
    na_e: float,
    k_e: float,
    cl_e: float,
    hco3_e: float,
    ph_e: float,
    chemistry: CarbonateChemistry = DEFAULT_CARBONATE_CHEMISTRY,
) -> ReservoirActivities:
    """Build an explicit SR4 transporter context without hiding carbonate."""

    return ReservoirActivities(
        na_i=na_i,
        k_i=k_i,
        cl_i=cl_i,
        hco3_i=hco3_i,
        na_e=na_e,
        k_e=k_e,
        cl_e=cl_e,
        hco3_e=hco3_e,
        co3_i=chemistry.carbonate_from_bicarbonate(hco3_i, ph_i),
        co3_e=chemistry.carbonate_from_bicarbonate(hco3_e, ph_e),
    )
