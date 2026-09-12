"""Change AE4 cation routing while preserving the inherited net chloride law.

This is a controlled routing intervention, not a new fitted AE4 rate law.
The original shared carrier evaluator supplies one net cycle J. Its cation
source is divided by donor concentrations with equal Na/K weighting. No
parameter, anion source, capacity or regulator is changed.

The inherited microscopic branch affinities describe the original model;
they are not presented as thermodynamic validation of this routing change.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from .model import ModernFullModel
from .nkcc_stimulation import StimulatedNkcc1Model
from .transporters import AE4Environment, AE4Parameters, evaluate_ae4_qss


MODEL_ID = "SR2_NET_CHLORIDE_DONOR_CATION_ROUTING"


@dataclass(frozen=True)
class RoutedAE4Flux:
    model_id: str
    units: str
    j_ae4_fmol_s: float
    j_na_fmol_s: float
    j_k_fmol_s: float
    donor_na_fraction: float
    donor_side: str
    source_na_i_fmol_s: float
    source_k_i_fmol_s: float
    source_cl_i_fmol_s: float
    source_hco3_i_fmol_s: float
    source_tic_i_fmol_s: float
    source_alkalinity_i_fmol_s: float
    transported_charge_equivalents_fmol_s: float
    regulation_gain: float
    reference_j_na_fmol_s: float
    reference_j_k_fmol_s: float

    @property
    def intracellular_conserved_sources_fmol_s(self):
        return {"na_i": self.source_na_i_fmol_s,
                "k_i": self.source_k_i_fmol_s,
                "cl_i": self.source_cl_i_fmol_s,
                "tic_i": self.source_tic_i_fmol_s,
                "alkalinity_i": self.source_alkalinity_i_fmol_s}


def evaluate_routed_ae4(environment: AE4Environment,
                        parameters: AE4Parameters, *,
                        regulation_gain: float = 1.0) -> RoutedAE4Flux:
    reference = evaluate_ae4_qss(environment, parameters,
                                 regulation_gain=regulation_gain)
    # Preserve exactly the existing net chloride driving law, including all
    # attempt rates, carrier occupancy, capacity and regulation.
    j = reference.source_cl_i_fmol_s
    side = "i" if j >= 0.0 else "o"
    na = getattr(environment, f"na_{side}_mM")
    k = getattr(environment, f"k_{side}_mM")
    fraction = na / (na + k)
    j_na = fraction * j
    j_k = (1.0 - fraction) * j
    source_na, source_k = -j_na, -j_k
    charge = math.fsum((source_na, source_k, -j,
                        -reference.source_alkalinity_i_fmol_s))
    return RoutedAE4Flux(
        model_id=MODEL_ID, units="fmol/s", j_ae4_fmol_s=j,
        j_na_fmol_s=j_na, j_k_fmol_s=j_k,
        donor_na_fraction=fraction,
        donor_side="intracellular" if side == "i" else "extracellular",
        source_na_i_fmol_s=source_na, source_k_i_fmol_s=source_k,
        source_cl_i_fmol_s=reference.source_cl_i_fmol_s,
        source_hco3_i_fmol_s=reference.source_hco3_i_fmol_s,
        source_tic_i_fmol_s=reference.source_tic_i_fmol_s,
        source_alkalinity_i_fmol_s=reference.source_alkalinity_i_fmol_s,
        transported_charge_equivalents_fmol_s=charge,
        regulation_gain=regulation_gain,
        reference_j_na_fmol_s=reference.j_na_fmol_s,
        reference_j_k_fmol_s=reference.j_k_fmol_s,
    )


def routing_only_adapter(environment: Mapping[str, float],
                         parameters: AE4Parameters,
                         regulation_gain: float) -> RoutedAE4Flux:
    values = {f"{ion}_{side}_mM": environment[
        f"{ion}_{'e' if side == 'o' else side}_mM"]
        for ion in ("na", "k", "cl", "hco3") for side in ("i", "o")}
    return evaluate_routed_ae4(AE4Environment(**values), parameters,
                               regulation_gain=regulation_gain)


def with_routing_only(model):
    """Retain the complete existing model and replace only its AE4 evaluator."""
    if model.ae4_evaluator is not None or not isinstance(model.ae4_parameters, AE4Parameters):
        raise ValueError("Routing intervention requires the inherited SR2 AE4 law")
    base = ModernFullModel(parameters=model.parameters, stimulus=model.stimulus,
                           regulatory_model=model.regulatory_model,
                           ae4_parameters=model.ae4_parameters,
                           ae4_evaluator=routing_only_adapter)
    return StimulatedNkcc1Model(base) if isinstance(model, StimulatedNkcc1Model) else base
