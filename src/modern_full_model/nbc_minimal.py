"""Stimulus-recruited electrogenic 1 Na : 2 HCO3 NBC surrogate.

This module replaces the earlier electroneutral 1:1 draft. It implements the
minimal architecture required by the conserved alkalinity, carbon, cation and
charge balances:

    Na_o + 2 HCO3_o <-> Na_i + 2 HCO3_i

Positive flux is bath -> cell. Each inward cycle therefore adds one Na amount,
two TIC amounts and two alkalinity equivalents to the cell, carries net -1
charge into the cell, and has no direct chloride source.

The NBC is an increment recruited by the existing muscarinic/Ca secretory arm.
It is exactly zero at the Task 31 resting calcium (0.058 uM), so the accepted
Task 31 R09 resting state is an exact nesting limit. At the standard central
Ca=0.25 uM condition the increment is fully recruited. The same normalized
secretory coordinate raises the already-built mechanistic NHE1 activity from
1.0 at rest to the source-fixed 2.3-fold stimulated value. No new kinetic
state or time constant is introduced.

The default NBC capacity is a derived WT stimulated feasibility scale, not a
measured transporter abundance and not a knockout-phenotype fit. It is chosen
from the algebraic 70/30 chloride-loading constraint on the Task 31 R09 flux
scale, then corrected for self-consistent electrogenic current closure at the
Task 31 reference chemical state.

This is a structural NBCe-like test only. It does not assign a specific SLC4
isoform to mouse submandibular acinar cells.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .membranes import (
    ElectricalEnvironment,
    FMOL_PER_MOL,
    HomeostasisEnvironment,
    MembraneFluxes,
    current_to_fmol_s,
    evaluate_homeostasis,
    evaluate_membrane_closure,
    hill_activation,
    nernst_voltage_V,
)
from .model import Genotype, ModelEvaluation, WT


RESTING_CALCIUM_UM = 0.058
FULLY_RECRUITED_CALCIUM_UM = 0.25
SOURCE_FIXED_NHE1_STIMULATED_MULTIPLIER = 2.3

TASK31_R09_NKCC1_CL_FMOL_S = 0.2562444047732706
TASK31_R09_NHE1_FMOL_S = 0.007890580489555009
TARGET_NKCC1_POSITIVE_CL_SHARE = 0.70
DERIVED_AE4_CL_FMOL_S = (
    (1.0 - TARGET_NKCC1_POSITIVE_CL_SHARE)
    / TARGET_NKCC1_POSITIVE_CL_SHARE
    * TASK31_R09_NKCC1_CL_FMOL_S
)
DERIVED_STIMULATED_NHE1_FMOL_S = (
    SOURCE_FIXED_NHE1_STIMULATED_MULTIPLIER * TASK31_R09_NHE1_FMOL_S
)
DERIVED_REQUIRED_NBC_CYCLE_FMOL_S = (
    DERIVED_AE4_CL_FMOL_S - 0.5 * DERIVED_STIMULATED_NHE1_FMOL_S
)

# Full electrogenic current closure at the Task 31 R09 reference chemical state
# gives this capacity for the derived cycle flux above. It is a structural
# feasibility scale, not a measurement.
DERIVED_REFERENCE_NBC_CAPACITY_FMOL_S = 0.11570913197464398


@dataclass(frozen=True)
class MinimalNbcParameters:
    """Fixed parameters for the stimulus-recruited 1 Na : 2 HCO3 surrogate."""

    capacity_fmol_s: float = DERIVED_REFERENCE_NBC_CAPACITY_FMOL_S
    log_width: float = 2.0
    resting_calcium_uM: float = RESTING_CALCIUM_UM
    fully_recruited_calcium_uM: float = FULLY_RECRUITED_CALCIUM_UM
    nhe1_stimulated_multiplier: float = SOURCE_FIXED_NHE1_STIMULATED_MULTIPLIER

    def __post_init__(self) -> None:
        if not math.isfinite(self.capacity_fmol_s) or self.capacity_fmol_s < 0.0:
            raise ValueError("NBC capacity must be finite and nonnegative")
        if not math.isfinite(self.log_width) or self.log_width <= 0.0:
            raise ValueError("NBC log width must be finite and positive")
        if (
            not math.isfinite(self.resting_calcium_uM)
            or not math.isfinite(self.fully_recruited_calcium_uM)
            or self.resting_calcium_uM < 0.0
            or self.fully_recruited_calcium_uM <= self.resting_calcium_uM
        ):
            raise ValueError("NBC calcium endpoints must define a positive interval")
        if (
            not math.isfinite(self.nhe1_stimulated_multiplier)
            or self.nhe1_stimulated_multiplier < 1.0
        ):
            raise ValueError("stimulated NHE1 multiplier must be finite and >= 1")


@dataclass(frozen=True)
class MinimalNbcFlux:
    """NBC cycle flux, species sources and electrogenic current."""

    cycle_inward_fmol_s: float
    activation_fraction: float
    affinity_log: float
    current_A: float
    na_cell_fmol_s: float
    cl_cell_fmol_s: float
    tic_cell_fmol_s: float
    alkalinity_cell_fmol_s: float
    transported_charge_equivalents_fmol_s: float
    charge_source_residual_fmol_s: float

    @property
    def hco3_cell_fmol_s(self) -> float:
        return self.tic_cell_fmol_s


def normalized_secretory_activation(
    calcium_uM: float,
    parameters: MinimalNbcParameters | None = None,
) -> float:
    """Protocol-normalized Ca arm, exactly zero at REST and one at 0.25 uM."""

    p = parameters or MinimalNbcParameters()
    calcium = float(calcium_uM)
    if not math.isfinite(calcium) or calcium < 0.0:
        raise ValueError("calcium must be finite and nonnegative")
    return min(
        1.0,
        max(
            0.0,
            (calcium - p.resting_calcium_uM)
            / (p.fully_recruited_calcium_uM - p.resting_calcium_uM),
        ),
    )


def nhe1_stimulation_multiplier(
    activation_fraction: float,
    parameters: MinimalNbcParameters | None = None,
) -> float:
    """Interpolate only between the frozen REST gain 1 and stimulated gain 2.3."""

    p = parameters or MinimalNbcParameters()
    u = float(activation_fraction)
    if not math.isfinite(u) or not 0.0 <= u <= 1.0:
        raise ValueError("secretory activation fraction must lie in [0, 1]")
    return 1.0 + (p.nhe1_stimulated_multiplier - 1.0) * u


def nbc_affinity_log(
    *,
    na_i_mM: float,
    hco3_i_mM: float,
    na_o_mM: float,
    hco3_o_mM: float,
    v_basolateral_V: float,
    thermal_voltage_V: float,
) -> float:
    """Dimensionless inward affinity for Na_o + 2HCO3_o -> Na_i + 2HCO3_i."""

    values = (
        na_i_mM,
        hco3_i_mM,
        na_o_mM,
        hco3_o_mM,
        thermal_voltage_V,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("NBC concentrations and thermal voltage must be finite and positive")
    if not math.isfinite(v_basolateral_V):
        raise ValueError("basolateral voltage must be finite")
    return (
        math.log((na_o_mM * hco3_o_mM**2) / (na_i_mM * hco3_i_mM**2))
        + v_basolateral_V / thermal_voltage_V
    )


def evaluate_minimal_nbc(
    *,
    na_i_mM: float,
    hco3_i_mM: float,
    na_o_mM: float,
    hco3_o_mM: float,
    v_basolateral_V: float,
    thermal_voltage_V: float,
    faraday_C_mol: float,
    activation_fraction: float,
    parameters: MinimalNbcParameters | None = None,
) -> MinimalNbcFlux:
    """Evaluate the reversible electrogenic 1 Na : 2 HCO3 cotransporter."""

    p = parameters or MinimalNbcParameters()
    if not math.isfinite(faraday_C_mol) or faraday_C_mol <= 0.0:
        raise ValueError("Faraday constant must be finite and positive")
    u = float(activation_fraction)
    if not math.isfinite(u) or not 0.0 <= u <= 1.0:
        raise ValueError("secretory activation fraction must lie in [0, 1]")
    affinity = nbc_affinity_log(
        na_i_mM=na_i_mM,
        hco3_i_mM=hco3_i_mM,
        na_o_mM=na_o_mM,
        hco3_o_mM=hco3_o_mM,
        v_basolateral_V=v_basolateral_V,
        thermal_voltage_V=thermal_voltage_V,
    )
    cycle = p.capacity_fmol_s * u * math.tanh(affinity / p.log_width)

    # One inward cycle imports +1 Na and +2 bicarbonates. The net cell charge
    # source is therefore -1 equivalent per cycle. Conventional current is
    # positive cell -> bath for inward transport of one net negative charge.
    na_source = cycle
    tic_source = 2.0 * cycle
    alkalinity_source = 2.0 * cycle
    charge_source = na_source - alkalinity_source
    current_A = -charge_source / FMOL_PER_MOL * faraday_C_mol

    return MinimalNbcFlux(
        cycle_inward_fmol_s=cycle,
        activation_fraction=u,
        affinity_log=affinity,
        current_A=current_A,
        na_cell_fmol_s=na_source,
        cl_cell_fmol_s=0.0,
        tic_cell_fmol_s=tic_source,
        alkalinity_cell_fmol_s=alkalinity_source,
        transported_charge_equivalents_fmol_s=charge_source,
        charge_source_residual_fmol_s=charge_source,
    )


def _positive(*values: float) -> None:
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("electrical concentrations must be finite and positive")


def evaluate_membrane_closure_with_nbc(
    environment: ElectricalEnvironment,
    model_parameters: Any,
    *,
    nbc_parameters: MinimalNbcParameters | None = None,
    activation_fraction: float,
) -> tuple[MembraneFluxes, MinimalNbcFlux]:
    """Solve apical/basolateral QSS current closure including electrogenic NBC.

    The first current equation is linear and is eliminated analytically,
    leaving one strictly monotone scalar equation for V_b. No iterative
    parameter search is involved.
    """

    p_nbc = nbc_parameters or MinimalNbcParameters()
    u = float(activation_fraction)
    if not math.isfinite(u) or not 0.0 <= u <= 1.0:
        raise ValueError("secretory activation fraction must lie in [0, 1]")

    # Exact resting nesting: use the existing closure literally.
    if u == 0.0:
        membranes = evaluate_membrane_closure(environment, model_parameters)
        zero = evaluate_minimal_nbc(
            na_i_mM=environment.na_i_mM,
            hco3_i_mM=environment.hco3_i_mM,
            na_o_mM=environment.na_e_mM,
            hco3_o_mM=environment.hco3_e_mM,
            v_basolateral_V=membranes.v_basolateral_V,
            thermal_voltage_V=model_parameters.constants.thermal_voltage_V,
            faraday_C_mol=model_parameters.constants.faraday_C_mol,
            activation_fraction=0.0,
            parameters=p_nbc,
        )
        return membranes, zero

    e = environment
    _positive(
        e.na_i_mM,
        e.k_i_mM,
        e.cl_i_mM,
        e.hco3_i_mM,
        e.na_l_mM,
        e.k_l_mM,
        e.cl_l_mM,
        e.hco3_l_mM,
        e.na_e_mM,
        e.k_e_mM,
        e.cl_e_mM,
        e.hco3_e_mM,
    )
    if e.calcium_uM < 0.0:
        raise ValueError("calcium input must be nonnegative")

    mp = model_parameters.membranes
    constants = model_parameters.constants
    fp = mp.apical_pump_fraction
    fk = mp.apical_k_fraction

    # Pump rates are voltage-independent in the existing model. Reuse the
    # canonical evaluator to avoid creating a second pump implementation.
    canonical = evaluate_membrane_closure(e, model_parameters)
    j_pump_a = canonical.pump_apical_fmol_s
    j_pump_b = canonical.pump_basolateral_fmol_s

    gate = hill_activation(e.calcium_uM, mp.calcium_half_uM, mp.calcium_hill)
    g_k_a = mp.g_k_total_S * gate * fk + mp.g_apical_background_S
    g_k_b = mp.g_k_total_S * gate * (1.0 - fk) + mp.g_basolateral_background_S
    g_cl_a = mp.g_cl_apical_S * gate

    vt = constants.thermal_voltage_V
    e_k_a = nernst_voltage_V(
        e.k_i_mM, e.k_l_mM, valence=+1, thermal_voltage_V=vt
    )
    e_k_b = nernst_voltage_V(
        e.k_i_mM, e.k_e_mM, valence=+1, thermal_voltage_V=vt
    )
    e_cl_a = nernst_voltage_V(
        e.cl_i_mM, e.cl_l_mM, valence=-1, thermal_voltage_V=vt
    )
    para_paths = {
        "para_na": (mp.g_para_na_S, +1, e.na_l_mM, e.na_e_mM),
        "para_k": (mp.g_para_k_S, +1, e.k_l_mM, e.k_e_mM),
        "para_cl": (mp.g_para_cl_S, -1, e.cl_l_mM, e.cl_e_mM),
        "para_hco3": (mp.g_para_hco3_S, -1, e.hco3_l_mM, e.hco3_e_mM),
    }
    para_reversals = {
        name: nernst_voltage_V(c_l, c_e, valence=z, thermal_voltage_V=vt)
        for name, (_, z, c_l, c_e) in para_paths.items()
    }
    g_para = sum(item[0] for item in para_paths.values())

    pump_current_a = constants.faraday_C_mol * j_pump_a / FMOL_PER_MOL
    pump_current_b = constants.faraday_C_mol * j_pump_b / FMOL_PER_MOL
    g_apical = g_k_a + g_cl_a
    c_apical = -g_k_a * e_k_a - g_cl_a * e_cl_a + pump_current_a
    c_para = -sum(
        conductance * para_reversals[name]
        for name, (conductance, _, _, _) in para_paths.items()
    )

    denominator = g_apical + g_para
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("apical/parallel conductance closure is singular")

    def va_from_vb(v_b: float) -> float:
        return (g_para * v_b - c_apical + c_para) / denominator

    def balance(v_b: float) -> float:
        v_a = va_from_vb(v_b)
        v_trans = v_b - v_a
        i_para = sum(
            conductance * (v_trans - para_reversals[name])
            for name, (conductance, _, _, _) in para_paths.items()
        )
        nbc = evaluate_minimal_nbc(
            na_i_mM=e.na_i_mM,
            hco3_i_mM=e.hco3_i_mM,
            na_o_mM=e.na_e_mM,
            hco3_o_mM=e.hco3_e_mM,
            v_basolateral_V=v_b,
            thermal_voltage_V=vt,
            faraday_C_mol=constants.faraday_C_mol,
            activation_fraction=u,
            parameters=p_nbc,
        )
        i_basolateral_old = g_k_b * (v_b - e_k_b) + pump_current_b
        return i_basolateral_old + nbc.current_A + i_para

    lower, upper = -0.5, 0.5
    f_lower, f_upper = balance(lower), balance(upper)
    if not (math.isfinite(f_lower) and math.isfinite(f_upper)):
        raise ValueError("NBC current closure produced a non-finite voltage residual")
    if f_lower == 0.0:
        v_b = lower
    elif f_upper == 0.0:
        v_b = upper
    elif f_lower * f_upper > 0.0:
        raise ValueError("NBC current closure is not bracketed on [-0.5, 0.5] V")
    else:
        v_b = float(brentq(balance, lower, upper, xtol=1e-14, rtol=1e-13))

    v_a = va_from_vb(v_b)
    v_trans = v_b - v_a
    nbc = evaluate_minimal_nbc(
        na_i_mM=e.na_i_mM,
        hco3_i_mM=e.hco3_i_mM,
        na_o_mM=e.na_e_mM,
        hco3_o_mM=e.hco3_e_mM,
        v_basolateral_V=v_b,
        thermal_voltage_V=vt,
        faraday_C_mol=constants.faraday_C_mol,
        activation_fraction=u,
        parameters=p_nbc,
    )

    i_k_a = g_k_a * (v_a - e_k_a)
    i_cl_a = g_cl_a * (v_a - e_cl_a)
    i_k_b = g_k_b * (v_b - e_k_b)
    para_currents = {
        name: conductance * (v_trans - para_reversals[name])
        for name, (conductance, _, _, _) in para_paths.items()
    }
    i_para = sum(para_currents.values())
    i_apical = i_k_a + i_cl_a + pump_current_a
    i_basolateral_old = i_k_b + pump_current_b
    i_basolateral_total = i_basolateral_old + nbc.current_A

    j_k_a = current_to_fmol_s(
        i_k_a, valence=+1, faraday_C_mol=constants.faraday_C_mol
    )
    j_cl_a = current_to_fmol_s(
        i_cl_a, valence=-1, faraday_C_mol=constants.faraday_C_mol
    )
    j_k_b = current_to_fmol_s(
        i_k_b, valence=+1, faraday_C_mol=constants.faraday_C_mol
    )
    para_molar = {
        name: current_to_fmol_s(
            para_currents[name], valence=z, faraday_C_mol=constants.faraday_C_mol
        )
        for name, (_, z, _, _) in para_paths.items()
    }

    # Existing channel/pump sources plus the NBC conserved-coordinate source.
    cell_sources = {
        "na": -3.0 * (j_pump_a + j_pump_b) + nbc.na_cell_fmol_s,
        "k": 2.0 * (j_pump_a + j_pump_b) - j_k_a - j_k_b,
        "cl": -j_cl_a,
        "tic": nbc.tic_cell_fmol_s,
        "alkalinity": nbc.alkalinity_cell_fmol_s,
    }
    lumen_sources = {
        "na": 3.0 * j_pump_a - para_molar["para_na"],
        "k": -2.0 * j_pump_a + j_k_a - para_molar["para_k"],
        "cl": j_cl_a - para_molar["para_cl"],
        "tic": -para_molar["para_hco3"],
        "alkalinity": -para_molar["para_hco3"],
    }

    cell_charge = (
        cell_sources["na"]
        + cell_sources["k"]
        - cell_sources["cl"]
        - cell_sources["alkalinity"]
    )
    lumen_charge = (
        lumen_sources["na"]
        + lumen_sources["k"]
        - lumen_sources["cl"]
        - lumen_sources["alkalinity"]
    )
    electrical_cell_charge = (
        -(i_apical + i_basolateral_total)
        / constants.faraday_C_mol
        * FMOL_PER_MOL
    )
    electrical_lumen_charge = (
        (i_apical - i_para) / constants.faraday_C_mol * FMOL_PER_MOL
    )

    currents = {
        "k_apical": i_k_a,
        "cl_apical": i_cl_a,
        "pump_apical": pump_current_a,
        "k_basolateral": i_k_b,
        "pump_basolateral": pump_current_b,
        "nbc_basolateral": nbc.current_A,
        **para_currents,
        "apical_total": i_apical,
        "basolateral_old": i_basolateral_old,
        "basolateral_total": i_basolateral_total,
        "paracellular_total": i_para,
    }
    return (
        MembraneFluxes(
            v_apical_V=v_a,
            v_basolateral_V=v_b,
            v_transepithelial_V=v_trans,
            pump_apical_fmol_s=j_pump_a,
            pump_basolateral_fmol_s=j_pump_b,
            currents_A=currents,
            cell_sources_fmol_s=cell_sources,
            lumen_sources_fmol_s=lumen_sources,
            current_residuals_A={
                "apical": i_apical - i_para,
                "basolateral": i_basolateral_total + i_para,
            },
            charge_residuals_fmol_s={
                "cell_source_minus_current": cell_charge - electrical_cell_charge,
                "lumen_source_minus_current": lumen_charge - electrical_lumen_charge,
            },
            capacity_partition_residuals={
                "pump_fraction": fp + (1.0 - fp) - 1.0,
                "k_fraction": fk + (1.0 - fk) - 1.0,
            },
        ),
        nbc,
    )


class MinimalNbcModel:
    """Whole-cell wrapper for the coupled stimulated NHE1 + NBC architecture."""

    def __init__(
        self,
        base_model: Any,
        parameters: MinimalNbcParameters | None = None,
    ) -> None:
        self.base_model = base_model
        self.nbc_parameters = parameters or MinimalNbcParameters()

    @property
    def nkcc1_kinetics(self) -> Any:
        return self.base_model.nkcc1_kinetics

    @property
    def parameters(self) -> Any:
        return self.base_model.parameters

    @property
    def stimulus(self) -> Any:
        return self.base_model.stimulus

    @property
    def regulatory_model(self) -> Any:
        return getattr(self.base_model, "regulatory_model", None)

    @property
    def ae4_parameters(self) -> Any:
        return getattr(self.base_model, "ae4_parameters", None)

    @property
    def ae4_evaluator(self) -> Any:
        return getattr(self.base_model, "ae4_evaluator", None)

    @property
    def layout(self) -> Any:
        return self.base_model.layout

    @property
    def state_names(self) -> tuple[str, ...]:
        return self.base_model.state_names

    def initial_state(self) -> NDArray[np.float64]:
        return self.base_model.initial_state()

    def _homeostasis_with_stimulated_nhe1(
        self,
        baseline: ModelEvaluation,
        genotype: Genotype,
        *,
        nhe1_multiplier: float,
    ) -> Any:
        obs = baseline.diagnostics.observables
        ci = obs.cell_concentrations_mM
        bath = self.parameters.bath
        bath_ab = obs.bath_acid_base
        regulatory = baseline.diagnostics.regulatory
        nkcc_multiplier = float(regulatory.get("nkcc1_capacity_multiplier", 1.0))
        if not math.isfinite(nkcc_multiplier) or nkcc_multiplier < 0.0:
            raise ValueError("NKCC1 regulatory multiplier must be finite and nonnegative")
        return evaluate_homeostasis(
            HomeostasisEnvironment(
                na_i_mM=ci["na"],
                k_i_mM=ci["k"],
                cl_i_mM=ci["cl"],
                h_i_mM=obs.cell_acid_base.h_mM,
                hco3_i_mM=ci["hco3"],
                na_e_mM=bath.na_mM,
                k_e_mM=bath.k_mM,
                cl_e_mM=bath.cl_mM,
                h_e_mM=bath_ab.h_mM,
                hco3_e_mM=bath_ab.hco3_mM,
            ),
            self.parameters,
            nkcc1_scale=genotype.nkcc1_expression * nkcc_multiplier,
            nkcc1_kinetics=self.nkcc1_kinetics,
            nhe1_scale=genotype.nhe1_expression * nhe1_multiplier,
            ae2_scale=genotype.ae2_expression,
        )

    def evaluate(
        self,
        time_s: float,
        vector: ArrayLike,
        *,
        genotype: Genotype = WT,
    ) -> ModelEvaluation:
        baseline = self.base_model.evaluate(time_s, vector, genotype=genotype)
        stimulus = baseline.diagnostics.stimulus
        activation = normalized_secretory_activation(
            stimulus.calcium_uM, self.nbc_parameters
        )
        nhe_multiplier = nhe1_stimulation_multiplier(
            activation, self.nbc_parameters
        )

        # Exact Task 31 REST nesting. Add diagnostics only; do not recalculate
        # any baseline flux or voltage when the recruitment coordinate is zero.
        if activation == 0.0:
            mem = baseline.diagnostics.membranes
            obs = baseline.diagnostics.observables
            zero = evaluate_minimal_nbc(
                na_i_mM=obs.cell_concentrations_mM["na"],
                hco3_i_mM=obs.cell_concentrations_mM["hco3"],
                na_o_mM=self.parameters.bath.na_mM,
                hco3_o_mM=obs.bath_acid_base.hco3_mM,
                v_basolateral_V=mem.v_basolateral_V,
                thermal_voltage_V=self.parameters.constants.thermal_voltage_V,
                faraday_C_mol=self.parameters.constants.faraday_C_mol,
                activation_fraction=0.0,
                parameters=self.nbc_parameters,
            )
            regulatory = dict(baseline.diagnostics.regulatory)
            regulatory.update(
                {
                    "minimal_nbc_model": "STIMULUS_RECRUITED_ELECTROGENIC_1NA_2HCO3",
                    "minimal_nbc_activation_fraction": 0.0,
                    "minimal_nbc_cycle_inward_fmol_s": 0.0,
                    "minimal_nbc_current_A": 0.0,
                    "minimal_nbc_affinity_log": zero.affinity_log,
                    "minimal_nbc_capacity_fmol_s": self.nbc_parameters.capacity_fmol_s,
                    "nhe1_stimulation_multiplier": 1.0,
                }
            )
            conservation = dict(baseline.diagnostics.conservation_residuals)
            conservation.update(
                {
                    "minimal_nbc_charge_source_fmol_s": 0.0,
                    "minimal_nbc_carbon_source_fmol_s": 0.0,
                }
            )
            diagnostics = replace(
                baseline.diagnostics,
                regulatory=regulatory,
                conservation_residuals=conservation,
            )
            return replace(baseline, diagnostics=diagnostics)

        obs = baseline.diagnostics.observables
        ci = obs.cell_concentrations_mM
        li = obs.lumen_concentrations_mM
        bath = self.parameters.bath
        bath_ab = obs.bath_acid_base

        homeostasis = self._homeostasis_with_stimulated_nhe1(
            baseline, genotype, nhe1_multiplier=nhe_multiplier
        )
        membranes, nbc = evaluate_membrane_closure_with_nbc(
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
            nbc_parameters=self.nbc_parameters,
            activation_fraction=activation,
        )

        previous_homeostasis = baseline.diagnostics.homeostasis
        previous_membranes = baseline.diagnostics.membranes
        raw = np.asarray(baseline.rhs, dtype=float).copy()

        cell_fields = ("na", "k", "cl", "tic", "alkalinity")
        homeostasis_fields = (
            "na_cell_fmol_s",
            "k_cell_fmol_s",
            "cl_cell_fmol_s",
            "tic_cell_fmol_s",
            "alkalinity_cell_fmol_s",
        )
        for row, (source_name, homeostasis_name) in enumerate(
            zip(cell_fields, homeostasis_fields)
        ):
            raw[row] += float(getattr(homeostasis, homeostasis_name)) - float(
                getattr(previous_homeostasis, homeostasis_name)
            )
            raw[row] += float(membranes.cell_sources_fmol_s[source_name]) - float(
                previous_membranes.cell_sources_fmol_s[source_name]
            )

        for row, source_name in zip(range(6, 11), cell_fields):
            raw[row] += float(membranes.lumen_sources_fmol_s[source_name]) - float(
                previous_membranes.lumen_sources_fmol_s[source_name]
            )

        outflow = baseline.diagnostics.outflow_sources_fmol_s
        qdot_cell = raw[0] + raw[1] - raw[2] - raw[4]
        qdot_lumen = raw[6] + raw[7] - raw[8] - raw[10]
        outflow_charge = (
            outflow["na"] + outflow["k"] - outflow["cl"] - outflow["alkalinity"]
        )
        expected_lumen_qdot = -outflow_charge

        j_co2_b = baseline.diagnostics.co2_fluxes_fmol_s["bath_to_cell"]
        external_carbon_source = (
            homeostasis.tic_cell_fmol_s
            + baseline.diagnostics.ae4.tic_cell_fmol_s
            + membranes.cell_sources_fmol_s["tic"]
            + j_co2_b
            + membranes.lumen_sources_fmol_s["tic"]
            - outflow["tic"]
        )
        carbon_residual = raw[3] + raw[9] - external_carbon_source

        conservation = dict(baseline.diagnostics.conservation_residuals)
        conservation.update(
            {
                "cell_bulk_charge_rate_fmol_s": float(qdot_cell),
                "lumen_bulk_charge_rate_minus_outflow_fmol_s": float(
                    qdot_lumen - expected_lumen_qdot
                ),
                "carbon_accounting_fmol_s": float(carbon_residual),
                "homeostasis_charge_fmol_s": float(
                    homeostasis.charge_source_residual_fmol_s
                ),
                "minimal_nbc_charge_source_fmol_s": float(
                    nbc.charge_source_residual_fmol_s
                ),
                "minimal_nbc_carbon_source_fmol_s": float(nbc.tic_cell_fmol_s),
                "apical_current_A": float(membranes.current_residuals_A["apical"]),
                "basolateral_current_A": float(
                    membranes.current_residuals_A["basolateral"]
                ),
            }
        )

        regulatory = dict(baseline.diagnostics.regulatory)
        regulatory.update(
            {
                "minimal_nbc_model": "STIMULUS_RECRUITED_ELECTROGENIC_1NA_2HCO3",
                "minimal_nbc_activation_fraction": activation,
                "minimal_nbc_cycle_inward_fmol_s": float(nbc.cycle_inward_fmol_s),
                "minimal_nbc_current_A": float(nbc.current_A),
                "minimal_nbc_affinity_log": float(nbc.affinity_log),
                "minimal_nbc_capacity_fmol_s": self.nbc_parameters.capacity_fmol_s,
                "nhe1_stimulation_multiplier": nhe_multiplier,
                "minimal_nbc_evidence_status": (
                    "STRUCTURAL_NBCE_LIKE_TEST; 1NA_2HCO3; "
                    "CAPACITY_DERIVED_FROM_WT_STIMULATED_FEASIBILITY"
                ),
            }
        )

        diagnostics = replace(
            baseline.diagnostics,
            regulatory=regulatory,
            homeostasis=homeostasis,
            membranes=membranes,
            conservation_residuals=conservation,
        )
        return replace(baseline, rhs=raw, diagnostics=diagnostics)

    def rhs(
        self,
        time_s: float,
        vector: ArrayLike,
        *,
        genotype: Genotype = WT,
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
        y0 = (
            self.initial_state()
            if initial_state is None
            else self.layout.validate(initial_state)
        )
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


__all__ = (
    "DERIVED_AE4_CL_FMOL_S",
    "DERIVED_REFERENCE_NBC_CAPACITY_FMOL_S",
    "DERIVED_REQUIRED_NBC_CYCLE_FMOL_S",
    "DERIVED_STIMULATED_NHE1_FMOL_S",
    "FULLY_RECRUITED_CALCIUM_UM",
    "MinimalNbcFlux",
    "MinimalNbcModel",
    "MinimalNbcParameters",
    "RESTING_CALCIUM_UM",
    "SOURCE_FIXED_NHE1_STIMULATED_MULTIPLIER",
    "TARGET_NKCC1_POSITIVE_CL_SHARE",
    "evaluate_membrane_closure_with_nbc",
    "evaluate_minimal_nbc",
    "nbc_affinity_log",
    "nhe1_stimulation_multiplier",
    "normalized_secretory_activation",
)
