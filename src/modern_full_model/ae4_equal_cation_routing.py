"""Task 40: fixed equal Na/K routing of the unchanged inherited AE4 cycle.

Only the ensemble cation partition changes. The shared-carrier QSS evaluator
still determines J4, including its concentration dependence and regulation.
No donor concentration is read to choose a fraction; no fraction is tunable.
This structural hypothesis is not a new microscopic thermodynamic claim.
"""
from dataclasses import dataclass
import math
from typing import Mapping

from .transporters import AE4Environment, AE4Parameters, evaluate_ae4_qss

MODEL_ID = "SR2_NET_CHLORIDE_EQUAL_NA_K_CATION_ROUTING"


@dataclass(frozen=True)
class EqualRoutedAE4Flux:
    j_ae4_fmol_s: float
    regulation_gain: float
    model_id: str = MODEL_ID
    units: str = "fmol/s"

    @property
    def j_na_fmol_s(self):
        return 0.5 * self.j_ae4_fmol_s

    @property
    def j_k_fmol_s(self):
        return 0.5 * self.j_ae4_fmol_s

    @property
    def source_na_i_fmol_s(self):
        return -0.5 * self.j_ae4_fmol_s

    @property
    def source_k_i_fmol_s(self):
        return -0.5 * self.j_ae4_fmol_s

    @property
    def source_cl_i_fmol_s(self):
        return self.j_ae4_fmol_s

    @property
    def source_hco3_i_fmol_s(self):
        return -2.0 * self.j_ae4_fmol_s

    source_tic_i_fmol_s = source_hco3_i_fmol_s
    source_alkalinity_i_fmol_s = source_hco3_i_fmol_s

    @property
    def transported_charge_equivalents_fmol_s(self):
        return math.fsum((self.source_na_i_fmol_s, self.source_k_i_fmol_s,
                          -self.source_cl_i_fmol_s,
                          -self.source_alkalinity_i_fmol_s))

    @property
    def intracellular_conserved_sources_fmol_s(self):
        return {"na_i": self.source_na_i_fmol_s,
                "k_i": self.source_k_i_fmol_s,
                "cl_i": self.source_cl_i_fmol_s,
                "tic_i": self.source_tic_i_fmol_s,
                "alkalinity_i": self.source_alkalinity_i_fmol_s}


def evaluate_equal_routed_ae4(environment: AE4Environment,
                              parameters: AE4Parameters, *,
                              regulation_gain: float = 1.0) -> EqualRoutedAE4Flux:
    reference = evaluate_ae4_qss(environment, parameters,
                                 regulation_gain=regulation_gain)
    return EqualRoutedAE4Flux(reference.source_cl_i_fmol_s, regulation_gain)


def equal_routing_adapter(environment: Mapping[str, float],
                          parameters: AE4Parameters,
                          regulation_gain: float) -> EqualRoutedAE4Flux:
    values = {f"{ion}_{side}_mM": environment[
        f"{ion}_{'e' if side == 'o' else side}_mM"]
        for ion in ("na", "k", "cl", "hco3") for side in ("i", "o")}
    return evaluate_equal_routed_ae4(AE4Environment(**values), parameters,
                                      regulation_gain=regulation_gain)
