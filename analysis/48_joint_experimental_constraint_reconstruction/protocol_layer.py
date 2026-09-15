"""Task 48 protocol plumbing; inherited intracellular equations are unchanged.

The positive B+ bath is an explicitly qualified external reservoir projection.
Unsupported source assays raise a descriptive error instead of supplying a
surrogate observation. No parameter inference occurs in this module.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from modern_full_model import parameters as par
from modern_full_model.acid_base import total_alkalinity_mM
from modern_full_model.model import ModernFullModel, Genotype
from modern_full_model.transporters import AE4Parameters
from modern_full_model.ae4_equal_cation_routing import equal_routing_adapter
from modern_full_model.nkcc1_palk2010 import Nkcc1Kinetics
from modern_full_model.nkcc_stimulation import (
    N1AlgebraicNkcc1, CompositeAe4NkccRegulation, StimulatedNkcc1Model,
)
from modern_full_model.camp_pka import R1EffectiveActivation, RegulatoryGain, AE4Construct
from modern_full_model.nbc_minimal import MinimalNbcModel, MinimalNbcParameters
from modern_full_model.validation import SecretagogueProtocol, StimulusArm

FROZEN_INPUT = ROOT / "analysis/47_dynamic_potassium_recycling_reconstruction/input/active_parameters.json"
GENOTYPES = {
    "WT": Genotype("WT"),
    "AE4_KO": Genotype("AE4_KO", ae4_expression=0),
    "AE2_KO": Genotype("AE2_KO", ae2_expression=0),
    "DOUBLE_KO": Genotype("DOUBLE_KO", ae4_expression=0, ae2_expression=0),
    "NKCC1_KO": Genotype("NKCC1_KO", nkcc1_expression=0),
    "NHE1_KO": Genotype("NHE1_KO", nhe1_expression=0),
}


class UnrepresentableProtocol(ValueError):
    pass


@dataclass(frozen=True)
class Interventions:
    bumetanide: bool = False
    t16ainh_a01: bool = False
    eipa: bool = False
    ethoxyzolamide: bool = False

    def effective_genotype(self, genotype: Genotype) -> Genotype:
        """Acute masks affect their target, never the pre-assay resting state."""
        return replace(
            genotype,
            nkcc1_expression=0.0 if self.bumetanide else genotype.nkcc1_expression,
            nhe1_expression=0.0 if self.eipa else genotype.nhe1_expression,
        )


BLOCKED_PROTOCOLS = {
    "NKCC_ISOLATION": ["Bicarbonate absent outside inherited positive domain",
                       "Carbonic anhydrase kinetics absent",
                       "Selected NKCC law has no extracellular substrate dependence",
                       "Source conditioning and slope windows not numerically specified"],
    "NKCC_ISOLATION_BUMETANIDE": ["Same missing protocol as NKCC_ISOLATION"],
    "EXCHANGER_ISOLATION": ["Low chloride predepletion uses a bath independent NKCC law",
                           "Conditioning durations and source regression windows unavailable",
                           "Isolated acinar apical bath access not specified"],
    "EXCHANGER_IPR": ["EXCHANGER_ISOLATION limitations and IPR exposure time unavailable"],
    "EXCHANGER_NA_FREE": ["Extracellular Na=0 lies outside inherited rate and Nernst domains",
                          "EXCHANGER_ISOLATION limitations"],
}


def require_protocol(name: str) -> None:
    if name in BLOCKED_PROTOCOLS:
        raise UnrepresentableProtocol("; ".join(BLOCKED_PROTOCOLS[name]))
    if name not in {"REST", "CCH_ONLY", "IPR_ONLY", "CCH_IPR", "NHE", "NHE_EIPA"}:
        raise KeyError(name)


def shared_payload() -> dict:
    return json.loads(FROZEN_INPUT.read_text())


def payload_hash(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


class SourceBathCore(ModernFullModel):
    """Retain cellular validation and account explicitly for bath spectators.

    Only the positive high chloride B+ projection is admitted. Calcium and
    magnesium supply +4 mEq/L; the small additional carbonate/titrant term
    follows the inherited speciation convention. Neither is an inferred
    cellular parameter. The original intracellular equations are inherited.
    """
    def _validate_reference_electroneutrality(self) -> None:
        p = self.parameters
        b = p.bath
        if not np.allclose([b.na_mM, b.k_mM, b.cl_mM, b.ph], [145., 4.3, 128.3, 7.4], rtol=0, atol=1e-12):
            raise ValueError("SourceBathCore requires the declared B+ recipe")
        ta_bath = self._bath_speciation.total_alkalinity_mM
        spectator_charge = 4.0 + ta_bath - 25.0
        bath_charge = b.na_mM + b.k_mM - b.cl_mM - ta_bath + spectator_charge
        initial = p.initial
        cell_ta = total_alkalinity_mM(
            ph=initial.cell_ph, total_carbon_mM=initial.cell_tic_mM,
            buffer_total_mM=p.geometry.cell_buffer_total_fmol / initial.cell_volume_pL,
            buffer_pka=p.acid_base.cell_buffer_pka, parameters=p.acid_base,
        )
        lumen_ta = total_alkalinity_mM(
            ph=initial.lumen_ph, total_carbon_mM=initial.lumen_tic_mM,
            buffer_total_mM=p.geometry.lumen_buffer_total_fmol / initial.lumen_volume_pL,
            buffer_pka=p.acid_base.lumen_buffer_pka, parameters=p.acid_base,
        )
        residuals = {
            "bath_charge_mM": bath_charge,
            "cell_charge_fmol": (initial.cell_na_mM + initial.cell_k_mM - initial.cell_cl_mM - cell_ta) * initial.cell_volume_pL - p.geometry.fixed_cell_anion_equivalents_fmol,
            "lumen_charge_fmol": (initial.lumen_na_mM + initial.lumen_k_mM - initial.lumen_cl_mM - lumen_ta) * initial.lumen_volume_pL,
        }
        if max(map(abs, residuals.values())) > 1e-9:
            raise ValueError(f"Charge accounting failed: {residuals}")
        if abs(self._bath_speciation.hco3_mM - 25.) > 1e-8:
            raise ValueError("Source bicarbonate was not preserved")
        self.bath_projection = {
            "qualification": "effective positive B+ reservoir; unreported HEPES/titrant chemistry omitted",
            "tracked_ions_mM": {"Na": b.na_mM, "K": b.k_mM, "Cl": b.cl_mM, "HCO3": self._bath_speciation.hco3_mM},
            "TIC_mM": b.tic_mM, "divalent_charge_mEq_L": 4.0,
            "carbonate_titrant_accounting_mEq_L": ta_bath - 25.,
            "total_spectator_charge_mEq_L": spectator_charge,
            "spectator_osmoles_mM": b.untracked_osmolyte_mM,
            "charge_checks": residuals,
        }


class ProtocolModel:
    """Ensure every public RHS and integration entry applies acute masks."""
    def __init__(self, model, interventions):
        self.model = model
        self.task48_interventions = interventions

    def __getattr__(self, name):
        return getattr(self.model, name)

    def evaluate(self, t, state, *, genotype=GENOTYPES["WT"]):
        return self.model.evaluate(t, state, genotype=self.task48_interventions.effective_genotype(genotype))

    def rhs(self, t, state, *, genotype=GENOTYPES["WT"]):
        return self.evaluate(t, state, genotype=genotype).rhs

    def solve_dynamics(self, time_span_s, *, initial_state, genotype=GENOTYPES["WT"],
                       method="BDF", rtol=1e-8, atol=1e-10, t_eval=None, max_step_s=2.):
        return solve_ivp(lambda t, y: self.rhs(t, y, genotype=genotype),
                         time_span_s, initial_state, method=method, rtol=rtol,
                         atol=atol, t_eval=t_eval, max_step=max_step_s)


def build_model(arm="REST", *, interventions=Interventions(), payload=None, source_bath=True):
    require_protocol(arm)
    if arm in {"NHE", "NHE_EIPA"}:
        if arm == "NHE_EIPA": interventions = replace(interventions, eipa=True)
        arm = "CCH_IPR"
    if interventions.ethoxyzolamide:
        raise UnrepresentableProtocol("No carbonic anhydrase kinetic parameter in equilibrium carbonate chemistry")
    data = shared_payload() if payload is None else payload
    types = dict(constants=par.PhysicalConstants, acid_base=par.AcidBaseParameters,
                 bath=par.BathParameters, geometry=par.GeometryParameters,
                 homeostasis=par.HomeostasisParameters, membranes=par.MembraneParameters,
                 water=par.WaterParameters, initial=par.InitialConditions)
    parameters = par.FullModelParameters(**{k: cls(**data["parameters"][k]) for k, cls in types.items()})
    if source_bath:
        acid = parameters.acid_base
        tic = 25. * (1. + 10. ** (acid.carbon_pka1 - 7.4) + 10. ** (7.4 - acid.carbon_pka2))
        parameters = replace(parameters, bath=replace(parameters.bath,
            na_mM=145., k_mM=4.3, cl_mM=128.3, tic_mM=tic, ph=7.4,
            buffer_total_mM=0., untracked_osmolyte_mM=17.))
    if interventions.t16ainh_a01:
        parameters = replace(parameters, membranes=replace(parameters.membranes, g_cl_apical_S=0.))
    stimulus_data = {**data["stimulus"], "arm": StimulusArm(arm)}
    stimulus = SecretagogueProtocol(**stimulus_data)
    rp = dict(data["regulatory_model"]["ae4_regulatory_model"])
    rp["construct"] = AE4Construct(rp["construct"])
    rp["gain"] = RegulatoryGain(**rp["gain"])
    regulation = CompositeAe4NkccRegulation(
        R1EffectiveActivation(**rp),
        N1AlgebraicNkcc1(**data["regulatory_model"]["nkcc1_regulatory_model"]), stimulus)
    core_type = SourceBathCore if source_bath else ModernFullModel
    core = core_type(parameters, stimulus=stimulus, regulatory_model=regulation,
        ae4_parameters=AE4Parameters(**data["ae4_parameters"]),
        ae4_evaluator=equal_routing_adapter,
        nkcc1_kinetics=Nkcc1Kinetics(**data["nkcc1_kinetics"]))
    model = MinimalNbcModel(StimulatedNkcc1Model(core), MinimalNbcParameters(**data["nbc_parameters"]))
    model.task48_interventions = interventions
    model.task48_payload_hash = payload_hash(data)
    model.task48_bath_projection = getattr(core, "bath_projection", {"qualification": "inherited bath for nesting tests only"})
    return ProtocolModel(model, interventions)


def evaluate(model, t, state, genotype):
    return model.evaluate(t, state, genotype=genotype)


def original_seed():
    return np.asarray(json.loads((ROOT / "results/40_ae4_equal_cation_routing/wt_rest.json").read_text())["state_vector"])


def concentration_derivative(state, rhs, amount_index=2, volume_index=5):
    volume = float(state[volume_index])
    c = float(state[amount_index]) / volume
    return (float(rhs[amount_index]) - c * float(rhs[volume_index])) / volume


def spq_inverse_ratio(chloride, initial_chloride, *, quenching_per_mM):
    if quenching_per_mM is None:
        raise ValueError("Source optical calibration unavailable: specify a conditional quenching coefficient")
    if not math.isfinite(quenching_per_mM) or quenching_per_mM <= 0:
        raise ValueError("Quenching coefficient must be positive")
    return (1 + quenching_per_mM * np.asarray(chloride)) / (1 + quenching_per_mM * initial_chloride)


def uptake_slope(times, observations, *, window):
    if window is None:
        raise ValueError("A documented recovery window is required; initial exit is not uptake")
    t, y = np.asarray(times), np.asarray(observations)
    selected = (t >= window[0]) & (t <= window[1])
    if selected.sum() < 3 or window[1] <= window[0]:
        raise ValueError("At least three points in a positive duration window required")
    t = t[selected] - t[selected].mean()
    return float(t @ y[selected] / (t @ t))


def gland_integrals(times, flow, bicarbonate):
    t, q, b = map(np.asarray, (times, flow, bicarbonate))
    return {"water_pL": float(np.trapezoid(q, t)), "bicarbonate_fmol": float(np.trapezoid(q * b, t))}


def gaussian_ratio_contrast(predicted_ratio, control_mean, control_se, ko_mean, ko_se):
    return (ko_mean - predicted_ratio * control_mean) / math.hypot(ko_se, predicted_ratio * control_se)


def rest_cache_key(model, genotype):
    return payload_hash({"shared_parameters": model.task48_payload_hash,
                         "genotype": asdict(genotype), "bath": asdict(model.parameters.bath),
                         "interventions": asdict(model.task48_interventions),
                         "stimulus": asdict(model.stimulus)})
