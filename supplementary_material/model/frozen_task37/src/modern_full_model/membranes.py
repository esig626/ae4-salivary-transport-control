"""Thermodynamic homeostasis laws and two-membrane current closure.

Sign conventions
----------------
``I > 0`` is conventional positive-charge current from the first named
compartment to the second.  For a species of valence ``z``, the corresponding
molar flux is ``J = I/(z F)``.  Thus a positive conventional chloride current
represents chloride movement opposite to the named positive direction.

``V_a = phi_cell - phi_lumen`` and ``V_b = phi_cell - phi_bath``.  Therefore
``V_t = phi_lumen - phi_bath = V_b - V_a``.  The algebraic closure is

``I_apical - I_paracellular(lumen->bath) = 0`` and
``I_basolateral + I_paracellular(lumen->bath) = 0``.

The apical/basolateral pump and K fractions partition total *capacity* and
total Ca-activated conductance.  Their values are modeling parameters, not
claimed localization measurements.  Setting both fractions to zero recovers
the historical basolateral-only cation topology as an exact nesting limit.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np

from .parameters import FullModelParameters
from .nhe1_cha2009 import (
    NHE1_CHA_2009,
    cha_nhe1_flux_fmol_s,
    published_cha_kinetics,
)


FMOL_PER_MOL = 1.0e15


def _positive(*values: float) -> None:
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("thermodynamic concentrations must be finite and positive")


def nernst_voltage_V(
    concentration_first_mM: float,
    concentration_second_mM: float,
    *,
    valence: int,
    thermal_voltage_V: float,
) -> float:
    """Equilibrium ``phi_first - phi_second`` for the named ion."""

    _positive(concentration_first_mM, concentration_second_mM, thermal_voltage_V)
    if valence == 0:
        raise ValueError("Nernst voltage is undefined for a neutral species")
    return (thermal_voltage_V / valence) * math.log(
        concentration_second_mM / concentration_first_mM
    )


def current_to_fmol_s(current_A: float, *, valence: int, faraday_C_mol: float) -> float:
    if valence == 0:
        raise ValueError("current conversion requires a charged species")
    return current_A / (valence * faraday_C_mol) * FMOL_PER_MOL


def hill_activation(input_value: float, half_value: float, exponent: float) -> float:
    if input_value < 0.0 or half_value <= 0.0 or exponent <= 0.0:
        raise ValueError("Hill-law arguments are outside their physical domain")
    if input_value == 0.0:
        return 0.0
    numerator = input_value**exponent
    return numerator / (numerator + half_value**exponent)


def _bounded_reversible_flux(
    capacity_fmol_s: float, log_activity_ratio: float, width: float
) -> float:
    """Odd saturating law that preserves the exact thermodynamic reversal."""

    if capacity_fmol_s < 0.0 or width <= 0.0:
        raise ValueError("flux capacity must be nonnegative and width positive")
    return capacity_fmol_s * math.tanh(log_activity_ratio / width)


NHE1_LEGACY_TANH = "legacy_tanh"
NHE1_VERA_SIGUENZA_2018 = "vera_siguenza_2018_eq26"


def evaluate_nhe1_legacy_tanh(
    environment: "HomeostasisEnvironment",
    parameters: FullModelParameters,
    *,
    activity_scale: float = 1.0,
) -> float:
    """Evaluate the retained generic reversible NHE1 comparator in fmol/s."""

    e = environment
    _positive(e.na_i_mM, e.h_i_mM, e.na_e_mM, e.h_e_mM)
    if not math.isfinite(activity_scale) or activity_scale < 0.0:
        raise ValueError("NHE1 activity scale must be finite and nonnegative")
    affinity = math.log((e.na_e_mM * e.h_i_mM) / (e.na_i_mM * e.h_e_mM))
    return _bounded_reversible_flux(
        parameters.homeostasis.nhe1_capacity_fmol_s * activity_scale,
        affinity,
        parameters.homeostasis.thermodynamic_saturation_log_width,
    )


def evaluate_nhe1_vera_siguenza_2018(
    environment: "HomeostasisEnvironment",
    parameters: FullModelParameters,
    *,
    activity_scale: float = 1.0,
) -> float:
    """Evaluate Vera-Sigüenza et al. (2018), Eq. 26, in fmol/s.

    Proton concentrations and ``K_H`` are in mM.  Positive flux is one Na
    entering the cell while one H leaves it.  The squared proton saturation
    factors are retained exactly as printed in the paper.
    """

    e = environment
    p = parameters.homeostasis
    _positive(
        e.na_i_mM,
        e.h_i_mM,
        e.na_e_mM,
        e.h_e_mM,
        p.nhe1_published_k_h_mM,
        p.nhe1_published_k_na_mM,
    )
    if not math.isfinite(activity_scale) or activity_scale < 0.0:
        raise ValueError("NHE1 activity scale must be finite and nonnegative")
    if not math.isfinite(p.nhe1_published_g_fmol_s) or p.nhe1_published_g_fmol_s < 0.0:
        raise ValueError("published NHE1 activity must be finite and nonnegative")

    k_h = p.nhe1_published_k_h_mM
    k_na = p.nhe1_published_k_na_mM
    forward = (e.h_i_mM / (e.h_i_mM + k_h)) ** 2 * (
        e.na_e_mM / (e.na_e_mM + k_na)
    )
    reverse = (e.na_i_mM / (e.na_i_mM + k_na)) * (
        e.h_e_mM / (e.h_e_mM + k_h)
    ) ** 2
    return p.nhe1_published_g_fmol_s * activity_scale * (forward - reverse)


def evaluate_nhe1_cha_2009(
    environment: "HomeostasisEnvironment",
    parameters: FullModelParameters,
    *,
    activity_scale: float = 1.0,
) -> float:
    """Select verified Cha Table S1/Eq. 3/Mod2 and require a salivary amount."""
    kinetics = published_cha_kinetics()
    amount = parameters.homeostasis.nhe1_cha_carrier_amount_fmol
    if amount is None:
        raise ValueError("Task 31 salivary NHE1 amount has not been calibrated")
    e = environment
    return cha_nhe1_flux_fmol_s(
        na_i_mM=e.na_i_mM, h_i_mM=e.h_i_mM,
        na_o_mM=e.na_e_mM, h_o_mM=e.h_e_mM,
        kinetics=kinetics, carrier_amount_fmol=amount,
        activity_scale=activity_scale,
    )


@dataclass(frozen=True)
class HomeostasisEnvironment:
    na_i_mM: float
    k_i_mM: float
    cl_i_mM: float
    h_i_mM: float
    hco3_i_mM: float
    na_e_mM: float
    k_e_mM: float
    cl_e_mM: float
    h_e_mM: float
    hco3_e_mM: float


@dataclass(frozen=True)
class HomeostasisFluxes:
    """Basolateral electroneutral transporter rates and cell sources."""

    nkcc1_inward_fmol_s: float
    nhe1_inward_fmol_s: float
    nhe1_model: str
    ae2_inward_fmol_s: float
    na_cell_fmol_s: float
    k_cell_fmol_s: float
    cl_cell_fmol_s: float
    tic_cell_fmol_s: float
    alkalinity_cell_fmol_s: float
    affinities: Mapping[str, float]
    charge_source_residual_fmol_s: float


def evaluate_homeostasis(
    environment: HomeostasisEnvironment,
    parameters: FullModelParameters,
    *,
    nkcc1_scale: float = 1.0,
    nhe1_scale: float = 1.0,
    ae2_scale: float = 1.0,
) -> HomeostasisFluxes:
    """Evaluate NKCC1, selected NHE1, and AE2 laws.

    Positive directions are NKCC1 influx (Na + K + 2Cl), NHE1 Na influx/H
    extrusion, and AE2 Cl influx/HCO3 extrusion.  The log-activity affinities
    are exact ideal-solution deductions.  NKCC1 and AE2 retain the bounded
    ``tanh`` law.  NHE1 is selected explicitly in ``HomeostasisParameters``.
    """

    e = environment
    _positive(
        e.na_i_mM,
        e.k_i_mM,
        e.cl_i_mM,
        e.h_i_mM,
        e.hco3_i_mM,
        e.na_e_mM,
        e.k_e_mM,
        e.cl_e_mM,
        e.h_e_mM,
        e.hco3_e_mM,
    )
    if min(nkcc1_scale, nhe1_scale, ae2_scale) < 0.0:
        raise ValueError("transporter expression scales must be nonnegative")

    affinity_nkcc = math.log(
        (e.na_e_mM * e.k_e_mM * e.cl_e_mM**2)
        / (e.na_i_mM * e.k_i_mM * e.cl_i_mM**2)
    )
    affinity_nhe = math.log(
        (e.na_e_mM * e.h_i_mM) / (e.na_i_mM * e.h_e_mM)
    )
    affinity_ae2 = math.log(
        (e.cl_e_mM * e.hco3_i_mM) / (e.cl_i_mM * e.hco3_e_mM)
    )
    width = parameters.homeostasis.thermodynamic_saturation_log_width
    j_nkcc = _bounded_reversible_flux(
        parameters.homeostasis.nkcc1_capacity_fmol_s * nkcc1_scale,
        affinity_nkcc,
        width,
    )
    nhe1_model = parameters.homeostasis.nhe1_model
    if nhe1_model == NHE1_LEGACY_TANH:
        j_nhe = evaluate_nhe1_legacy_tanh(
            environment, parameters, activity_scale=nhe1_scale
        )
    elif nhe1_model == NHE1_VERA_SIGUENZA_2018:
        j_nhe = evaluate_nhe1_vera_siguenza_2018(
            environment, parameters, activity_scale=nhe1_scale
        )
    elif nhe1_model == NHE1_CHA_2009:
        j_nhe = evaluate_nhe1_cha_2009(
            environment, parameters, activity_scale=nhe1_scale
        )
    else:  # FullModelParameters validates this, but retain a local guard.
        raise ValueError(f"unsupported NHE1 model {nhe1_model!r}")
    j_ae2 = _bounded_reversible_flux(
        parameters.homeostasis.ae2_capacity_fmol_s * ae2_scale,
        affinity_ae2,
        width,
    )

    # NHE proton extrusion raises TA; AE2 HCO3 extrusion lowers TIC and TA.
    na_source = j_nkcc + j_nhe
    k_source = j_nkcc
    cl_source = 2.0 * j_nkcc + j_ae2
    tic_source = -j_ae2
    alkalinity_source = j_nhe - j_ae2
    charge_residual = na_source + k_source - cl_source - alkalinity_source
    return HomeostasisFluxes(
        nkcc1_inward_fmol_s=j_nkcc,
        nhe1_inward_fmol_s=j_nhe,
        nhe1_model=nhe1_model,
        ae2_inward_fmol_s=j_ae2,
        na_cell_fmol_s=na_source,
        k_cell_fmol_s=k_source,
        cl_cell_fmol_s=cl_source,
        tic_cell_fmol_s=tic_source,
        alkalinity_cell_fmol_s=alkalinity_source,
        affinities={"NKCC1": affinity_nkcc, "NHE1": affinity_nhe, "AE2": affinity_ae2},
        charge_source_residual_fmol_s=charge_residual,
    )


@dataclass(frozen=True)
class ElectricalEnvironment:
    na_i_mM: float
    k_i_mM: float
    cl_i_mM: float
    hco3_i_mM: float
    na_l_mM: float
    k_l_mM: float
    cl_l_mM: float
    hco3_l_mM: float
    na_e_mM: float
    k_e_mM: float
    cl_e_mM: float
    hco3_e_mM: float
    calcium_uM: float


@dataclass(frozen=True)
class MembraneFluxes:
    v_apical_V: float
    v_basolateral_V: float
    v_transepithelial_V: float
    pump_apical_fmol_s: float
    pump_basolateral_fmol_s: float
    currents_A: Mapping[str, float]
    cell_sources_fmol_s: Mapping[str, float]
    lumen_sources_fmol_s: Mapping[str, float]
    current_residuals_A: Mapping[str, float]
    charge_residuals_fmol_s: Mapping[str, float]
    capacity_partition_residuals: Mapping[str, float]


def _pump_cycle_rate_fmol_s(
    *,
    capacity_fmol_s: float,
    na_i_mM: float,
    k_out_mM: float,
    na_half_mM: float,
    k_half_mM: float,
) -> float:
    _positive(na_i_mM, k_out_mM, na_half_mM, k_half_mM)
    if capacity_fmol_s < 0.0:
        raise ValueError("pump capacity must be nonnegative")
    na_gate = na_i_mM**3 / (na_i_mM**3 + na_half_mM**3)
    k_gate = k_out_mM**2 / (k_out_mM**2 + k_half_mM**2)
    return capacity_fmol_s * na_gate * k_gate


def evaluate_membrane_closure(
    environment: ElectricalEnvironment, parameters: FullModelParameters
) -> MembraneFluxes:
    """Solve the coupled apical/basolateral current equations exactly."""

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
    p = parameters.membranes
    constants = parameters.constants
    fp = p.apical_pump_fraction
    fk = p.apical_k_fraction

    j_pump_a = _pump_cycle_rate_fmol_s(
        capacity_fmol_s=p.nak_capacity_fmol_s * fp,
        na_i_mM=e.na_i_mM,
        k_out_mM=e.k_l_mM,
        na_half_mM=p.nak_na_half_mM,
        k_half_mM=p.nak_k_half_mM,
    )
    j_pump_b = _pump_cycle_rate_fmol_s(
        capacity_fmol_s=p.nak_capacity_fmol_s * (1.0 - fp),
        na_i_mM=e.na_i_mM,
        k_out_mM=e.k_e_mM,
        na_half_mM=p.nak_na_half_mM,
        k_half_mM=p.nak_k_half_mM,
    )
    gate = hill_activation(e.calcium_uM, p.calcium_half_uM, p.calcium_hill)
    g_k_a = p.g_k_total_S * gate * fk + p.g_apical_background_S
    g_k_b = p.g_k_total_S * gate * (1.0 - fk) + p.g_basolateral_background_S
    g_cl_a = p.g_cl_apical_S * gate

    vt = constants.thermal_voltage_V
    e_k_a = nernst_voltage_V(e.k_i_mM, e.k_l_mM, valence=+1, thermal_voltage_V=vt)
    e_k_b = nernst_voltage_V(e.k_i_mM, e.k_e_mM, valence=+1, thermal_voltage_V=vt)
    e_cl_a = nernst_voltage_V(e.cl_i_mM, e.cl_l_mM, valence=-1, thermal_voltage_V=vt)
    para_paths = {
        "para_na": (p.g_para_na_S, +1, e.na_l_mM, e.na_e_mM),
        "para_k": (p.g_para_k_S, +1, e.k_l_mM, e.k_e_mM),
        "para_cl": (p.g_para_cl_S, -1, e.cl_l_mM, e.cl_e_mM),
        "para_hco3": (p.g_para_hco3_S, -1, e.hco3_l_mM, e.hco3_e_mM),
    }
    para_reversals = {
        name: nernst_voltage_V(c_l, c_e, valence=z, thermal_voltage_V=vt)
        for name, (_, z, c_l, c_e) in para_paths.items()
    }
    g_para = sum(item[0] for item in para_paths.values())

    # I = G V + C for each membrane/current group.
    pump_current_a = constants.faraday_C_mol * j_pump_a / FMOL_PER_MOL
    pump_current_b = constants.faraday_C_mol * j_pump_b / FMOL_PER_MOL
    g_apical = g_k_a + g_cl_a
    c_apical = -g_k_a * e_k_a - g_cl_a * e_cl_a + pump_current_a
    g_basolateral = g_k_b
    c_basolateral = -g_k_b * e_k_b + pump_current_b
    c_para = -sum(
        conductance * para_reversals[name]
        for name, (conductance, _, _, _) in para_paths.items()
    )

    matrix = np.asarray(
        ((g_apical + g_para, -g_para), (-g_para, g_basolateral + g_para)),
        dtype=float,
    )
    right = np.asarray((-c_apical + c_para, -c_basolateral - c_para), dtype=float)
    try:
        v_a, v_b = (float(value) for value in np.linalg.solve(matrix, right))
    except np.linalg.LinAlgError as exc:
        raise ValueError("membrane-current closure is singular") from exc
    v_trans = v_b - v_a

    i_k_a = g_k_a * (v_a - e_k_a)
    i_cl_a = g_cl_a * (v_a - e_cl_a)
    i_k_b = g_k_b * (v_b - e_k_b)
    para_currents = {
        name: conductance * (v_trans - para_reversals[name])
        for name, (conductance, _, _, _) in para_paths.items()
    }
    i_para = sum(para_currents.values())
    i_apical = i_k_a + i_cl_a + pump_current_a
    i_basolateral = i_k_b + pump_current_b

    j_k_a = current_to_fmol_s(i_k_a, valence=+1, faraday_C_mol=constants.faraday_C_mol)
    j_cl_a = current_to_fmol_s(i_cl_a, valence=-1, faraday_C_mol=constants.faraday_C_mol)
    j_k_b = current_to_fmol_s(i_k_b, valence=+1, faraday_C_mol=constants.faraday_C_mol)
    para_molar = {
        name: current_to_fmol_s(
            para_currents[name], valence=z, faraday_C_mol=constants.faraday_C_mol
        )
        for name, (_, z, _, _) in para_paths.items()
    }

    # Pump cycle: 3 Na cell->side, 2 K side->cell.  Channel fluxes are named
    # cell->side and paracellular fluxes lumen->bath.
    cell_sources = {
        "na": -3.0 * (j_pump_a + j_pump_b),
        "k": 2.0 * (j_pump_a + j_pump_b) - j_k_a - j_k_b,
        "cl": -j_cl_a,
        "tic": 0.0,
        "alkalinity": 0.0,
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
    electrical_cell_charge = -(i_apical + i_basolateral) / constants.faraday_C_mol * FMOL_PER_MOL
    electrical_lumen_charge = (i_apical - i_para) / constants.faraday_C_mol * FMOL_PER_MOL

    currents = {
        "k_apical": i_k_a,
        "cl_apical": i_cl_a,
        "pump_apical": pump_current_a,
        "k_basolateral": i_k_b,
        "pump_basolateral": pump_current_b,
        **para_currents,
        "apical_total": i_apical,
        "basolateral_total": i_basolateral,
        "paracellular_total": i_para,
    }
    return MembraneFluxes(
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
            "basolateral": i_basolateral + i_para,
        },
        charge_residuals_fmol_s={
            "cell_source_minus_current": cell_charge - electrical_cell_charge,
            "lumen_source_minus_current": lumen_charge - electrical_lumen_charge,
        },
        capacity_partition_residuals={
            "pump_fraction": fp + (1.0 - fp) - 1.0,
            "k_fraction": fk + (1.0 - fk) - 1.0,
        },
    )


def cation_topology_source_directions() -> np.ndarray:
    """Source signatures of one apical pump cycle and one apical K efflux.

    Rows are ``(Na_i, K_i, Na_l, K_l)``.  These two columns are linearly
    independent; redistribution can alter only this cation/lumen subspace
    directly and has no direct acid-base row.
    """

    return np.asarray(((-3.0, 0.0), (+2.0, -1.0), (+3.0, 0.0), (-2.0, +1.0)))


__all__ = (
    "NHE1_LEGACY_TANH",
    "NHE1_VERA_SIGUENZA_2018",
    "NHE1_CHA_2009",
    "ElectricalEnvironment",
    "FMOL_PER_MOL",
    "HomeostasisEnvironment",
    "HomeostasisFluxes",
    "MembraneFluxes",
    "cation_topology_source_directions",
    "current_to_fmol_s",
    "evaluate_homeostasis",
    "evaluate_nhe1_legacy_tanh",
    "evaluate_nhe1_vera_siguenza_2018",
    "evaluate_nhe1_cha_2009",
    "evaluate_membrane_closure",
    "hill_activation",
    "nernst_voltage_V",
)
