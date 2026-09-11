"""Task 20: one reversible, equal-weight Na+K AE4 cycle, without slip.

The 1:1:2 working stoichiometry and donor-side bookkeeping are declared
coarse-grained assumptions, not experimentally identified microscopic laws.
The legacy SR2 evaluator is not modified or used to compute this affinity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from .model import ModernFullModel
from .nkcc_stimulation import StimulatedNkcc1Model
from .transporters import AE4Parameters


MODEL_ID = "POOLED_CATION_112_NO_SLIP"
LEGACY_MODEL_ID = "SR2_SHARED_112_QSS"


@dataclass(frozen=True)
class PooledAE4Environment:
    """mM activities; one cation may vanish, but each pooled side is positive."""

    na_i_mM: float
    k_i_mM: float
    cl_i_mM: float
    hco3_i_mM: float
    na_o_mM: float
    k_o_mM: float
    cl_o_mM: float
    hco3_o_mM: float

    def __post_init__(self):
        for name, value in vars(self).items():
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
            if name.startswith(("cl_", "hco3_")) and value == 0:
                raise ValueError(f"{name} must be positive")
        for side in ("i", "o"):
            pool = getattr(self, f"na_{side}_mM") + getattr(self, f"k_{side}_mM")
            if not math.isfinite(pool) or pool <= 0:
                raise ValueError("Both pooled cation activities must be positive")


@dataclass(frozen=True)
class PooledAE4Parameters:
    """Unmeasured kinetic lump: carrier (fmol) times attempt rate (s^-1).

    Task 20 fixes attempt_rate_s=1, independently of legacy routing. The
    carrier reference is the inherited AE4 carrier for log-fold accounting;
    its cross-architecture magnitude is a gauge convention, not a measurement.
    """

    carrier_amount_fmol: float
    attempt_rate_s: float = 1.0

    def __post_init__(self):
        if not math.isfinite(self.carrier_amount_fmol) or self.carrier_amount_fmol < 0:
            raise ValueError("carrier_amount_fmol must be finite and nonnegative")
        if not math.isfinite(self.attempt_rate_s) or self.attempt_rate_s <= 0:
            raise ValueError("attempt_rate_s must be finite and positive")


@dataclass(frozen=True)
class PooledAE4Flux:
    model_id: str
    units: str
    j_ae4_fmol_s: float
    affinity: float
    k_forward_s: float
    k_reverse_s: float
    donor_na_fraction: float
    donor_side: str
    source_na_i_fmol_s: float
    source_k_i_fmol_s: float
    source_cl_i_fmol_s: float
    source_hco3_i_fmol_s: float
    source_tic_i_fmol_s: float
    source_alkalinity_i_fmol_s: float
    entropy_production_over_r_fmol_s: float
    transported_charge_equivalents_fmol_s: float
    local_detailed_balance_residual: float
    regulation_gain: float

    @property
    def intracellular_conserved_sources_fmol_s(self):
        return {"na_i": self.source_na_i_fmol_s, "k_i": self.source_k_i_fmol_s,
                "cl_i": self.source_cl_i_fmol_s, "tic_i": self.source_tic_i_fmol_s,
                "alkalinity_i": self.source_alkalinity_i_fmol_s}


def pooled_affinity(environment: PooledAE4Environment) -> float:
    e = environment
    return math.fsum((
        math.log(e.cl_o_mM) - math.log(e.cl_i_mM),
        math.log(e.na_i_mM + e.k_i_mM) - math.log(e.na_o_mM + e.k_o_mM),
        2 * (math.log(e.hco3_i_mM) - math.log(e.hco3_o_mM)),
    ))


def evaluate_pooled_ae4(environment: PooledAE4Environment,
                        parameters: PooledAE4Parameters, *,
                        regulation_gain: float = 1.0) -> PooledAE4Flux:
    if not math.isfinite(regulation_gain) or regulation_gain < 0:
        raise ValueError("regulation_gain must be finite and nonnegative")
    a = pooled_affinity(environment)
    forward = parameters.attempt_rate_s * math.exp(a / 2)
    reverse = parameters.attempt_rate_s * math.exp(-a / 2)
    # sinh avoids subtractive cancellation near equilibrium. No affinity clip,
    # branch affinity, selectivity fit, or donor interpolation is introduced.
    j = parameters.carrier_amount_fmol * regulation_gain * (
        2 * parameters.attempt_rate_s * math.sinh(a / 2)
    )
    side = "i" if j >= 0 else "o"
    na = getattr(environment, f"na_{side}_mM")
    k = getattr(environment, f"k_{side}_mM")
    p = na / (na + k)
    source_na = -p * j
    source_k = -(1 - p) * j
    source_hco3 = -2 * j
    charge = math.fsum((source_na, source_k, -j, -source_hco3))
    if not all(math.isfinite(x) for x in (j, forward, reverse, source_hco3)):
        raise FloatingPointError("Nonfinite pooled AE4 flux or effective rate")
    return PooledAE4Flux(
        model_id=MODEL_ID, units="fmol/s", j_ae4_fmol_s=j, affinity=a,
        k_forward_s=forward, k_reverse_s=reverse,
        donor_na_fraction=p, donor_side=("intracellular" if side == "i" else "extracellular"),
        source_na_i_fmol_s=source_na, source_k_i_fmol_s=source_k,
        source_cl_i_fmol_s=j, source_hco3_i_fmol_s=source_hco3,
        source_tic_i_fmol_s=source_hco3, source_alkalinity_i_fmol_s=source_hco3,
        entropy_production_over_r_fmol_s=j * a,
        transported_charge_equivalents_fmol_s=charge,
        local_detailed_balance_residual=abs(math.log(forward) - math.log(reverse) - a),
        regulation_gain=regulation_gain,
    )


def pooled_ae4_adapter(environment: Mapping[str, float], parameters: PooledAE4Parameters,
                       regulation_gain: float) -> PooledAE4Flux:
    """Use the existing production evaluator interface without changing it."""
    values = {f"{ion}_{side}_mM": environment[f"{ion}_{'e' if side == 'o' else side}_mM"]
              for ion in ("na", "k", "cl", "hco3") for side in ("i", "o")}
    return evaluate_pooled_ae4(PooledAE4Environment(**values), parameters,
                               regulation_gain=regulation_gain)


def select_ae4_mechanism(model, mechanism: str, *, parameters=None):
    """Select a law while preserving every non-AE4 parameter and regulator.

    Returning from pooled to legacy requires the original AE4Parameters; no
    legacy attempt rates are inferred from the pooled kinetic lump.
    """
    if mechanism == MODEL_ID:
        selected = parameters or PooledAE4Parameters(model.ae4_parameters.carrier_amount_fmol)
        if not isinstance(selected, PooledAE4Parameters):
            raise TypeError("Pooled mechanism requires PooledAE4Parameters")
        evaluator = pooled_ae4_adapter
    elif mechanism == LEGACY_MODEL_ID:
        selected = parameters or model.ae4_parameters
        if not isinstance(selected, AE4Parameters) or selected.cooperative_gate is not None:
            raise TypeError("SR2 legacy selection requires its original ungated AE4Parameters")
        evaluator = None
    else:
        raise ValueError(f"Unknown AE4 mechanism: {mechanism}")
    base = ModernFullModel(parameters=model.parameters, stimulus=model.stimulus,
                           regulatory_model=model.regulatory_model,
                           ae4_parameters=selected, ae4_evaluator=evaluator)
    return StimulatedNkcc1Model(base) if isinstance(model, StimulatedNkcc1Model) else base
