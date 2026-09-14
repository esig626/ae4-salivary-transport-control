"""Integrated conservation-explicit modern salivary acinar-cell model.

This module assembles the independently testable acid-base, homeostasis,
electrical, water, AE4, and beta/cAMP/PKA modules.  It is deliberately not a
calibration script: no AE4-null secretion target appears in this source.

Stable public API
-----------------
``ModernFullModel.state_names``
    Full named state vector, including regulatory suffix states.
``ModernFullModel.initial_state()``
    Auditable positive reference state on a physical amount/volume scale.
``ModernFullModel.evaluate(t, y, genotype=...)``
    RHS plus pH, currents, source terms, and conservation diagnostics.
``ModernFullModel.rhs(t, y, genotype=...)``
    Solver-compatible RHS in fmol/s, pL/s, and regulatory-state/s units.
``ModernFullModel.solve_resting_root(...)``
    Diagnostic bounded residual attempt only; never production eligible.
``ModernFullModel.solve_dynamics(..., method='Radau'|'BDF')``
    Stiff dynamic integration on an explicit seconds axis.

AE4 and regulatory objects use narrow duck-typed interfaces so their dedicated
Task 13B modules remain separable.  AE4 deletion is an exact source-level zero
without changing any unrelated parameter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Callable, Mapping, Protocol, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .acid_base import AcidBaseSpeciation, speciate, total_alkalinity_mM
from .membranes import (
    ElectricalEnvironment,
    HomeostasisEnvironment,
    HomeostasisFluxes,
    MembraneFluxes,
    evaluate_homeostasis,
    evaluate_membrane_closure,
)
from .parameters import FullModelParameters
from .states import (
    CORE_STATE_NAMES,
    CompartmentAmounts,
    StateLayout,
    WholeCellState,
    amount_charge_equivalents_fmol,
    compartment_from_concentrations,
)
from .water import (
    OsmoticEnvironment,
    WaterFluxes,
    cell_impermeant_osmoles_fmol,
    compartment_osmolarity_mOsm,
    evaluate_water_fluxes,
)


@dataclass(frozen=True)
class Stimulus:
    calcium_uM: float
    beta_input: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.calcium_uM) or self.calcium_uM < 0.0:
            raise ValueError("calcium stimulus must be finite and nonnegative")
        if not math.isfinite(self.beta_input) or not 0.0 <= self.beta_input <= 1.0:
            raise ValueError("beta input must be finite and lie in [0, 1]")


class StimulusProtocol(Protocol):
    def __call__(self, time_s: float) -> Stimulus: ...


@dataclass(frozen=True)
class ConstantStimulus:
    calcium_uM: float = 0.10
    beta_input: float = 0.0

    def __call__(self, time_s: float) -> Stimulus:
        if not math.isfinite(time_s):
            raise ValueError("time must be finite")
        return Stimulus(self.calcium_uM, self.beta_input)


@dataclass(frozen=True)
class Genotype:
    name: str
    ae4_expression: float = 1.0
    ae2_expression: float = 1.0
    nkcc1_expression: float = 1.0
    nhe1_expression: float = 1.0

    def __post_init__(self) -> None:
        for field_name in (
            "ae4_expression",
            "ae2_expression",
            "nkcc1_expression",
            "nhe1_expression",
        ):
            value = getattr(self, field_name)
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{field_name} must be finite and nonnegative")


WT = Genotype("WT")
AE4_NULL = Genotype("AE4_NULL", ae4_expression=0.0)
AE2_NULL = Genotype("AE2_NULL", ae2_expression=0.0)


@dataclass(frozen=True)
class CellObservables:
    cell_concentrations_mM: Mapping[str, float]
    lumen_concentrations_mM: Mapping[str, float]
    cell_acid_base: AcidBaseSpeciation
    lumen_acid_base: AcidBaseSpeciation
    bath_acid_base: AcidBaseSpeciation
    osmolarities_mOsm: Mapping[str, float]


@dataclass(frozen=True)
class AE4Contribution:
    na_cell_fmol_s: float = 0.0
    k_cell_fmol_s: float = 0.0
    cl_cell_fmol_s: float = 0.0
    hco3_cell_fmol_s: float = 0.0
    transported_charge_fmol_s: float = 0.0
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    @property
    def tic_cell_fmol_s(self) -> float:
        return self.hco3_cell_fmol_s

    @property
    def alkalinity_cell_fmol_s(self) -> float:
        return self.hco3_cell_fmol_s

    @property
    def charge_source_residual_fmol_s(self) -> float:
        return (
            self.na_cell_fmol_s
            + self.k_cell_fmol_s
            - self.cl_cell_fmol_s
            - self.hco3_cell_fmol_s
        )


@dataclass(frozen=True)
class ModelDiagnostics:
    stimulus: Stimulus
    regulatory: Mapping[str, Any]
    observables: CellObservables
    homeostasis: HomeostasisFluxes
    ae4: AE4Contribution
    membranes: MembraneFluxes
    water: WaterFluxes
    co2_fluxes_fmol_s: Mapping[str, float]
    outflow_sources_fmol_s: Mapping[str, float]
    state_charge_fmol: Mapping[str, float]
    conservation_residuals: Mapping[str, float]


@dataclass(frozen=True)
class ModelEvaluation:
    rhs: NDArray[np.float64]
    diagnostics: ModelDiagnostics


@dataclass(frozen=True)
class RestingRoot:
    """Bounded root-attempt result.

    ``max_abs_raw_rhs`` mixes fmol/s, pL/s, and (when present) regulatory
    state/s, so it is a debug summary only and is never a physical norm or an
    acceptance metric.  ``max_abs_scaled_rhs`` uses declared row scales.
    """

    state: NDArray[np.float64]
    success: bool
    max_abs_scaled_rhs: float
    max_abs_raw_rhs: float
    cost: float
    nfev: int
    message: str
    production_eligible: bool
    limitations: tuple[str, ...]


def _object_mapping(value: Any) -> Mapping[str, Any]:
    if value is None:
        return {}
    if isinstance(value, Mapping):
        return dict(value)
    if hasattr(value, "__dataclass_fields__"):
        return {name: getattr(value, name) for name in value.__dataclass_fields__}
    return {name: getattr(value, name) for name in dir(value) if not name.startswith("_")}


from .nkcc1_palk2010 import Nkcc1Kinetics


class ModernFullModel:
    """Twelve-state conserved whole-cell core plus optional regulatory states."""

    def __init__(
        self,
        parameters: FullModelParameters | None = None,
        *,
        stimulus: StimulusProtocol | None = None,
        regulatory_model: Any | None = None,
        ae4_parameters: Any | None = None,
        ae4_evaluator: Callable[..., Any] | None = None,
        nkcc1_kinetics: Nkcc1Kinetics = Nkcc1Kinetics(),
    ) -> None:
        self.parameters = parameters or FullModelParameters()
        self.stimulus = stimulus or ConstantStimulus()
        self.regulatory_model = regulatory_model
        self.ae4_parameters = ae4_parameters
        self.ae4_evaluator = ae4_evaluator
        self.nkcc1_kinetics = nkcc1_kinetics
        regulatory_names = (
            tuple(str(name) for name in regulatory_model.state_names)
            if regulatory_model is not None
            else ()
        )
        self.layout = StateLayout(regulatory_names)
        self._bath_speciation = speciate(
            total_carbon_mM=self.parameters.bath.tic_mM,
            total_alkalinity_mM_value=total_alkalinity_mM(
                ph=self.parameters.bath.ph,
                total_carbon_mM=self.parameters.bath.tic_mM,
                buffer_total_mM=self.parameters.bath.buffer_total_mM,
                buffer_pka=self.parameters.acid_base.lumen_buffer_pka,
                parameters=self.parameters.acid_base,
            ),
            buffer_total_mM=self.parameters.bath.buffer_total_mM,
            buffer_pka=self.parameters.acid_base.lumen_buffer_pka,
            parameters=self.parameters.acid_base,
        )
        self._validate_reference_electroneutrality()

    def _validate_reference_electroneutrality(self) -> None:
        """Reject stale derived charge constants after parameter changes.

        Bath chloride, the fixed-cell anion amount, and reference luminal
        chloride are derived for the active bath/reference tuple.  They are
        stored explicitly for provenance and freezing, so changing a parent
        concentration requires changing the corresponding derived value too.
        Silent charge drift is not allowed.
        """

        p = self.parameters
        tolerance = 1.0e-9
        bath_charge_mM = (
            p.bath.na_mM
            + p.bath.k_mM
            - p.bath.cl_mM
            - self._bath_speciation.total_alkalinity_mM
        )
        initial = p.initial
        cell_ta_mM = total_alkalinity_mM(
            ph=initial.cell_ph,
            total_carbon_mM=initial.cell_tic_mM,
            buffer_total_mM=(
                p.geometry.cell_buffer_total_fmol / initial.cell_volume_pL
            ),
            buffer_pka=p.acid_base.cell_buffer_pka,
            parameters=p.acid_base,
        )
        cell_charge_fmol = (
            (
                initial.cell_na_mM
                + initial.cell_k_mM
                - initial.cell_cl_mM
                - cell_ta_mM
            )
            * initial.cell_volume_pL
            - p.geometry.fixed_cell_anion_equivalents_fmol
        )
        lumen_ta_mM = total_alkalinity_mM(
            ph=initial.lumen_ph,
            total_carbon_mM=initial.lumen_tic_mM,
            buffer_total_mM=(
                p.geometry.lumen_buffer_total_fmol / initial.lumen_volume_pL
            ),
            buffer_pka=p.acid_base.lumen_buffer_pka,
            parameters=p.acid_base,
        )
        lumen_charge_fmol = (
            initial.lumen_na_mM
            + initial.lumen_k_mM
            - initial.lumen_cl_mM
            - lumen_ta_mM
        ) * initial.lumen_volume_pL
        residuals = {
            "bath_charge_mM": bath_charge_mM,
            "cell_reference_charge_fmol": cell_charge_fmol,
            "lumen_reference_charge_fmol": lumen_charge_fmol,
        }
        stale = {name: value for name, value in residuals.items() if abs(value) > tolerance}
        if stale:
            detail = ", ".join(f"{name}={value:.12g}" for name, value in stale.items())
            raise ValueError(
                "derived electroneutrality parameters are inconsistent with the active "
                f"bath/reference tuple ({detail}); update bath Cl, fixed-cell anion, "
                "or luminal Cl from the same tuple"
            )

    @property
    def state_names(self) -> tuple[str, ...]:
        return self.layout.names

    def initial_state(self) -> NDArray[np.float64]:
        """Return the positive, named reference state without fitting any KO data."""

        p = self.parameters
        initial = p.initial
        cell_buffer_mM = p.geometry.cell_buffer_total_fmol / initial.cell_volume_pL
        lumen_buffer_mM = p.geometry.lumen_buffer_total_fmol / initial.lumen_volume_pL
        cell_ta = total_alkalinity_mM(
            ph=initial.cell_ph,
            total_carbon_mM=initial.cell_tic_mM,
            buffer_total_mM=cell_buffer_mM,
            buffer_pka=p.acid_base.cell_buffer_pka,
            parameters=p.acid_base,
        )
        lumen_ta = total_alkalinity_mM(
            ph=initial.lumen_ph,
            total_carbon_mM=initial.lumen_tic_mM,
            buffer_total_mM=lumen_buffer_mM,
            buffer_pka=p.acid_base.lumen_buffer_pka,
            parameters=p.acid_base,
        )
        core = WholeCellState(
            cell=compartment_from_concentrations(
                na_mM=initial.cell_na_mM,
                k_mM=initial.cell_k_mM,
                cl_mM=initial.cell_cl_mM,
                tic_mM=initial.cell_tic_mM,
                alkalinity_mM=cell_ta,
                volume_pL=initial.cell_volume_pL,
            ),
            lumen=compartment_from_concentrations(
                na_mM=initial.lumen_na_mM,
                k_mM=initial.lumen_k_mM,
                cl_mM=initial.lumen_cl_mM,
                tic_mM=initial.lumen_tic_mM,
                alkalinity_mM=lumen_ta,
                volume_pL=initial.lumen_volume_pL,
            ),
            regulation=self._initial_regulatory_state(),
        )
        return self.layout.encode(core)

    def _initial_regulatory_state(self) -> tuple[float, ...]:
        if self.regulatory_model is None:
            return ()
        beta = self.stimulus(0.0).beta_input
        value = self.regulatory_model.initial_state(beta_input=beta)
        if hasattr(value, "as_vector"):
            value = value.as_vector()
        elif hasattr(value, "values") and not isinstance(value, Mapping):
            value = value.values
        if isinstance(value, Mapping):
            return tuple(float(value[name]) for name in self.layout.regulatory_names)
        return tuple(float(item) for item in value)

    def _decode_observables(self, state: WholeCellState) -> CellObservables:
        p = self.parameters
        ci = state.cell.concentrations()
        cl = state.lumen.concentrations()
        cell_ab = speciate(
            total_carbon_mM=ci.tic_mM,
            total_alkalinity_mM_value=ci.alkalinity_mM,
            buffer_total_mM=p.geometry.cell_buffer_total_fmol / state.cell.volume_pL,
            buffer_pka=p.acid_base.cell_buffer_pka,
            parameters=p.acid_base,
        )
        lumen_ab = speciate(
            total_carbon_mM=cl.tic_mM,
            total_alkalinity_mM_value=cl.alkalinity_mM,
            buffer_total_mM=p.geometry.lumen_buffer_total_fmol / state.lumen.volume_pL,
            buffer_pka=p.acid_base.lumen_buffer_pka,
            parameters=p.acid_base,
        )
        osm_cell = compartment_osmolarity_mOsm(
            na_mM=ci.na_mM,
            k_mM=ci.k_mM,
            cl_mM=ci.cl_mM,
            total_carbon_mM=ci.tic_mM,
            impermeant_osmoles_fmol=cell_impermeant_osmoles_fmol(
                other_impermeant_fmol=p.geometry.cell_impermeant_osmoles_fmol,
                finite_buffer_fmol=p.geometry.cell_buffer_total_fmol,
            ),
            volume_pL=state.cell.volume_pL,
        )
        osm_lumen = compartment_osmolarity_mOsm(
            na_mM=cl.na_mM,
            k_mM=cl.k_mM,
            cl_mM=cl.cl_mM,
            total_carbon_mM=cl.tic_mM,
        )
        bath = p.bath
        osm_bath = compartment_osmolarity_mOsm(
            na_mM=bath.na_mM,
            k_mM=bath.k_mM,
            cl_mM=bath.cl_mM,
            total_carbon_mM=bath.tic_mM,
            untracked_osmolyte_mM=bath.untracked_osmolyte_mM,
        )
        return CellObservables(
            cell_concentrations_mM={
                "na": ci.na_mM,
                "k": ci.k_mM,
                "cl": ci.cl_mM,
                "tic": ci.tic_mM,
                "alkalinity": ci.alkalinity_mM,
                "hco3": cell_ab.hco3_mM,
                "co2": cell_ab.co2_mM,
                "co3": cell_ab.co3_mM,
            },
            lumen_concentrations_mM={
                "na": cl.na_mM,
                "k": cl.k_mM,
                "cl": cl.cl_mM,
                "tic": cl.tic_mM,
                "alkalinity": cl.alkalinity_mM,
                "hco3": lumen_ab.hco3_mM,
                "co2": lumen_ab.co2_mM,
                "co3": lumen_ab.co3_mM,
            },
            cell_acid_base=cell_ab,
            lumen_acid_base=lumen_ab,
            bath_acid_base=self._bath_speciation,
            osmolarities_mOsm={"cell": osm_cell, "lumen": osm_lumen, "bath": osm_bath},
        )

    def _regulatory_evaluation(
        self, time_s: float, regulatory_state: Sequence[float], beta_input: float
    ) -> tuple[float, tuple[float, ...], Mapping[str, Any]]:
        if self.regulatory_model is None:
            return 1.0, (), {
                "family": "R0_UNCONFIGURED_STATIC_CONTROL",
                "activity": 0.0,
                "capacity_multiplier": 1.0,
            }
        output = self.regulatory_model.evaluate(time_s, regulatory_state, beta_input)
        mapping = _object_mapping(output)
        multiplier = float(mapping.get("capacity_multiplier", mapping.get("activity", 1.0)))
        if not math.isfinite(multiplier) or multiplier < 0.0:
            raise ValueError("AE4 regulatory multiplier must be finite and nonnegative")
        rates = self.regulatory_model.rhs(time_s, regulatory_state, beta_input)
        return multiplier, tuple(float(value) for value in rates), mapping

    def _evaluate_ae4(
        self,
        observables: CellObservables,
        *,
        regulatory_multiplier: float,
        expression_scale: float,
    ) -> AE4Contribution:
        if expression_scale == 0.0:
            return AE4Contribution(
                diagnostics={"status": "GENOTYPE_EXACT_ZERO", "expression_scale": 0.0}
            )
        if self.ae4_parameters is None and self.ae4_evaluator is None:
            return AE4Contribution(
                diagnostics={
                    "status": "AE4_MODULE_NOT_CONFIGURED",
                    "expression_scale": expression_scale,
                }
            )

        evaluator = self.ae4_evaluator
        environment: Any
        if evaluator is None:
            from .transporters import AE4Environment, evaluate_ae4_qss

            evaluator = evaluate_ae4_qss
            ci = observables.cell_concentrations_mM
            bath = self.parameters.bath
            environment = AE4Environment(
                na_i_mM=ci["na"],
                k_i_mM=ci["k"],
                cl_i_mM=ci["cl"],
                hco3_i_mM=ci["hco3"],
                na_o_mM=bath.na_mM,
                k_o_mM=bath.k_mM,
                cl_o_mM=bath.cl_mM,
                hco3_o_mM=self._bath_speciation.hco3_mM,
            )
        else:
            ci = observables.cell_concentrations_mM
            environment = {
                "na_i_mM": ci["na"],
                "k_i_mM": ci["k"],
                "cl_i_mM": ci["cl"],
                "hco3_i_mM": ci["hco3"],
                "na_e_mM": self.parameters.bath.na_mM,
                "k_e_mM": self.parameters.bath.k_mM,
                "cl_e_mM": self.parameters.bath.cl_mM,
                "hco3_e_mM": self._bath_speciation.hco3_mM,
            }
        if self.ae4_evaluator is None:
            result = evaluator(
                environment,
                self.ae4_parameters,
                regulation_gain=regulatory_multiplier,
            )
        else:
            result = evaluator(environment, self.ae4_parameters, regulatory_multiplier)
        result_map = _object_mapping(result)
        sources = getattr(result, "intracellular_conserved_sources_fmol_s", None)
        if sources is None:
            sources = result_map.get(
                "intracellular_sources_fmol_s", result_map.get("intracellular_sources", {})
            )
        if not isinstance(sources, Mapping):
            sources = _object_mapping(sources)

        def source(*keys: str) -> float:
            for key in keys:
                if key in sources:
                    return float(sources[key]) * expression_scale
                if key in result_map:
                    return float(result_map[key]) * expression_scale
            return 0.0

        charge = float(
            result_map.get(
                "transported_charge_fmol_s",
                result_map.get(
                    "transported_charge_equivalents_fmol_s",
                    result_map.get("transported_charge", result_map.get("charge_fmol_s", 0.0)),
                ),
            )
        ) * expression_scale
        return AE4Contribution(
            na_cell_fmol_s=source("na", "na_i", "na_cell_fmol_s", "na_i_fmol_s"),
            k_cell_fmol_s=source("k", "k_i", "k_cell_fmol_s", "k_i_fmol_s"),
            cl_cell_fmol_s=source("cl", "cl_i", "cl_cell_fmol_s", "cl_i_fmol_s"),
            hco3_cell_fmol_s=source(
                "alkalinity_i",
                "hco3",
                "hco3_i",
                "hco3_cell_fmol_s",
                "hco3_i_fmol_s",
            ),
            transported_charge_fmol_s=charge,
            diagnostics={**result_map, "expression_scale": expression_scale},
        )

    def evaluate(
        self, time_s: float, vector: ArrayLike, *, genotype: Genotype = WT
    ) -> ModelEvaluation:
        state = self.layout.decode(vector)
        observables = self._decode_observables(state)
        stimulus = self.stimulus(time_s)
        regulation, regulatory_rhs, regulatory_diagnostics = self._regulatory_evaluation(
            time_s, state.regulation, stimulus.beta_input
        )
        ci = observables.cell_concentrations_mM
        li = observables.lumen_concentrations_mM
        bath = self.parameters.bath
        bath_ab = observables.bath_acid_base
        homeostasis = evaluate_homeostasis(
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
            nkcc1_scale=genotype.nkcc1_expression,
            nkcc1_kinetics=self.nkcc1_kinetics,
            nhe1_scale=genotype.nhe1_expression,
            ae2_scale=genotype.ae2_expression,
        )
        ae4 = self._evaluate_ae4(
            observables,
            regulatory_multiplier=regulation,
            expression_scale=genotype.ae4_expression,
        )
        membranes = evaluate_membrane_closure(
            ElectricalEnvironment(
                na_i_mM=ci["na"],
                k_i_mM=ci["k"],
                cl_i_mM=ci["cl"],
                hco3_i_mM=ci["hco3"],
                na_l_mM=li["na"],
                k_l_mM=li["k"],
                cl_l_mM=li["cl"],
                hco3_l_mM=li["hco3"],
                na_e_mM=bath.na_mM,
                k_e_mM=bath.k_mM,
                cl_e_mM=bath.cl_mM,
                hco3_e_mM=bath_ab.hco3_mM,
                calcium_uM=stimulus.calcium_uM,
            ),
            self.parameters,
        )
        water = evaluate_water_fluxes(
            OsmoticEnvironment(
                osm_cell_mOsm=observables.osmolarities_mOsm["cell"],
                osm_lumen_mOsm=observables.osmolarities_mOsm["lumen"],
                osm_bath_mOsm=observables.osmolarities_mOsm["bath"],
                lumen_volume_pL=state.lumen.volume_pL,
            ),
            self.parameters,
        )

        hp = self.parameters.homeostasis
        j_co2_b = hp.co2_basolateral_permeability_fmol_s_mM * (
            bath_ab.co2_mM - ci["co2"]
        )
        # Positive apical value is lumen -> cell, so it cancels between the two
        # finite compartments exactly.
        j_co2_a = hp.co2_apical_permeability_fmol_s_mM * (
            li["co2"] - ci["co2"]
        )

        cell_electrical = membranes.cell_sources_fmol_s
        lumen_electrical = membranes.lumen_sources_fmol_s
        q_out = water.lumen_outflow_pL_s
        outflow = {
            "na": q_out * li["na"],
            "k": q_out * li["k"],
            "cl": q_out * li["cl"],
            "tic": q_out * li["tic"],
            "alkalinity": q_out * li["alkalinity"],
        }

        raw = np.zeros(self.layout.size, dtype=float)
        raw[0] = homeostasis.na_cell_fmol_s + ae4.na_cell_fmol_s + cell_electrical["na"]
        raw[1] = homeostasis.k_cell_fmol_s + ae4.k_cell_fmol_s + cell_electrical["k"]
        raw[2] = homeostasis.cl_cell_fmol_s + ae4.cl_cell_fmol_s + cell_electrical["cl"]
        raw[3] = (
            homeostasis.tic_cell_fmol_s
            + ae4.tic_cell_fmol_s
            + cell_electrical["tic"]
            + j_co2_b
            + j_co2_a
        )
        raw[4] = (
            homeostasis.alkalinity_cell_fmol_s
            + ae4.alkalinity_cell_fmol_s
            + cell_electrical["alkalinity"]
        )
        raw[5] = water.cell_volume_source_pL_s
        raw[6] = lumen_electrical["na"] - outflow["na"]
        raw[7] = lumen_electrical["k"] - outflow["k"]
        raw[8] = lumen_electrical["cl"] - outflow["cl"]
        raw[9] = lumen_electrical["tic"] - j_co2_a - outflow["tic"]
        raw[10] = lumen_electrical["alkalinity"] - outflow["alkalinity"]
        raw[11] = water.lumen_volume_source_pL_s
        if regulatory_rhs:
            raw[len(CORE_STATE_NAMES) :] = regulatory_rhs

        q_cell = amount_charge_equivalents_fmol(
            state.cell,
            fixed_anion_equivalents_fmol=(
                self.parameters.geometry.fixed_cell_anion_equivalents_fmol
            ),
        )
        q_lumen = amount_charge_equivalents_fmol(state.lumen)
        qdot_cell = raw[0] + raw[1] - raw[2] - raw[4]
        qdot_lumen = raw[6] + raw[7] - raw[8] - raw[10]
        outflow_charge = (
            outflow["na"] + outflow["k"] - outflow["cl"] - outflow["alkalinity"]
        )
        external_carbon_source = (
            homeostasis.tic_cell_fmol_s
            + ae4.tic_cell_fmol_s
            + j_co2_b
            + lumen_electrical["tic"]
            - outflow["tic"]
        )
        carbon_residual = raw[3] + raw[9] - external_carbon_source
        # qdot_lumen includes removal of any pre-existing nonzero bulk charge.
        expected_lumen_qdot = -outflow_charge
        conservation = {
            "cell_bulk_charge_rate_fmol_s": qdot_cell,
            "lumen_bulk_charge_rate_minus_outflow_fmol_s": (
                qdot_lumen - expected_lumen_qdot
            ),
            "carbon_accounting_fmol_s": carbon_residual,
            "buffer_site_accounting_fmol_s": 0.0,
            "water_volume_accounting_pL_s": water.volume_conservation_residual_pL_s,
            "homeostasis_charge_fmol_s": homeostasis.charge_source_residual_fmol_s,
            "ae4_charge_fmol_s": ae4.charge_source_residual_fmol_s,
            "apical_current_A": membranes.current_residuals_A["apical"],
            "basolateral_current_A": membranes.current_residuals_A["basolateral"],
            "cell_speciation_alkalinity_mM": observables.cell_acid_base.alkalinity_residual_mM,
            "lumen_speciation_alkalinity_mM": observables.lumen_acid_base.alkalinity_residual_mM,
        }
        return ModelEvaluation(
            rhs=raw,
            diagnostics=ModelDiagnostics(
                stimulus=stimulus,
                regulatory=regulatory_diagnostics,
                observables=observables,
                homeostasis=homeostasis,
                ae4=ae4,
                membranes=membranes,
                water=water,
                co2_fluxes_fmol_s={"bath_to_cell": j_co2_b, "lumen_to_cell": j_co2_a},
                outflow_sources_fmol_s=outflow,
                state_charge_fmol={"cell": q_cell, "lumen": q_lumen},
                conservation_residuals=conservation,
            ),
        )

    def rhs(self, time_s: float, vector: ArrayLike, *, genotype: Genotype = WT) -> NDArray[np.float64]:
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
        """Integrate with either independent SciPy stiff method."""

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

    def solve_resting_root(
        self,
        *,
        preferred: ArrayLike | None = None,
        genotype: Genotype = WT,
        start_count: int = 5,
        max_nfev: int = 3000,
        core_lower_factor: float = 0.10,
        core_upper_factor: float = 10.0,
    ) -> RestingRoot:
        """Run a diagnostic bounded residual minimization, not a production root.

        Bounds are relative to the preferred state so an optimizer cannot
        convert a failed physiological-root search into an unreported remote
        branch.  Calibration code should tighten these bounds from declared WT
        ranges rather than widening them until a root appears.

        This helper does not certify steady-state row rank, explicitly
        parameterize the charge manifold, track branches, or provide a
        sufficiently rich multistart design.  Its return is therefore always
        marked ``production_eligible=False``.  Task 13B calibration must use a
        separate root-geometry implementation before accepting a WT model.
        """

        if start_count < 1:
            raise ValueError("start_count must be positive")
        if not 0.0 < core_lower_factor < 1.0 < core_upper_factor:
            raise ValueError("root factors must satisfy 0 < lower < 1 < upper")
        reference = (
            self.initial_state()
            if preferred is None
            else self.layout.validate(preferred).copy()
        )
        lower = np.maximum(reference * core_lower_factor, 1.0e-12)
        upper = reference * core_upper_factor
        if self.layout.regulatory_names:
            lower[len(CORE_STATE_NAMES) :] = 0.0
            upper[len(CORE_STATE_NAMES) :] = 1.0
        # A 100 s source-to-state scale makes amount, volume, and regulatory
        # rows comparable without claiming a calibrated physical relaxation.
        scales = np.maximum(np.abs(reference) / 100.0, 1.0e-8)
        starts = [reference]
        factors = (0.85, 1.15, 0.70, 1.30)
        for factor in factors[: max(0, start_count - 1)]:
            candidate = reference.copy()
            candidate[: len(CORE_STATE_NAMES)] *= factor
            if self.layout.regulatory_names:
                candidate[len(CORE_STATE_NAMES) :] = np.clip(
                    candidate[len(CORE_STATE_NAMES) :], 1.0e-8, 1.0 - 1.0e-8
                )
            starts.append(candidate)

        results = []
        for start in starts:
            def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
                try:
                    return self.rhs(0.0, candidate, genotype=genotype) / scales
                except (ValueError, FloatingPointError, OverflowError):
                    return np.full(self.layout.size, 1.0e6)

            fit = least_squares(
                residual,
                np.clip(start, lower + 1.0e-14, upper),
                bounds=(lower, upper),
                max_nfev=max_nfev,
                xtol=1.0e-11,
                ftol=1.0e-11,
                gtol=1.0e-11,
            )
            raw = self.rhs(0.0, fit.x, genotype=genotype)
            results.append((fit, raw))
        fit, raw = min(results, key=lambda item: np.max(np.abs(item[1] / scales)))
        scaled_max = float(np.max(np.abs(raw / scales)))
        return RestingRoot(
            state=np.asarray(fit.x, dtype=float),
            success=bool(fit.success and scaled_max <= 1.0e-6),
            max_abs_scaled_rhs=scaled_max,
            max_abs_raw_rhs=float(np.max(np.abs(raw))),
            cost=float(fit.cost),
            nfev=int(fit.nfev),
            message=str(fit.message),
            production_eligible=False,
            limitations=(
                "steady-state row rank not certified",
                "charge manifold monitored but not parameterized",
                "multistart design is diagnostic and not space filling",
                "branch uniqueness not established",
            ),
        )


__all__ = (
    "AE2_NULL",
    "AE4_NULL",
    "WT",
    "AE4Contribution",
    "CellObservables",
    "ConstantStimulus",
    "Genotype",
    "ModelDiagnostics",
    "ModelEvaluation",
    "ModernFullModel",
    "RestingRoot",
    "Stimulus",
    "StimulusProtocol",
)
