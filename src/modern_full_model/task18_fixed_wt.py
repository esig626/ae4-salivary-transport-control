"""Task 18 fixed WT inverse transporter rebalancing.

This module is deliberately WT only.  It never imports a genotype runner and
never reads a genotype result.  The ten inherited conserved states are copied
without an encode/decode cycle.  At each state the chloride allocation fixes
the complete AE4 and NKCC1 transporter scales, after which the remaining
steady state balances are solved in pathway flux coordinates.  The resulting
fluxes are mapped back to existing positive capacities and the production
model is used for final verification.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, fields, replace
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import differential_evolution, minimize

from .calibration import CORE_INDEPENDENT_ROWS, WTCalibrationSpec
from .membranes import hill_activation, nernst_voltage_V
from .model import ModernFullModel, WT
from .nkcc_stimulation import StimulatedNkcc1Model
from .parameters import parameter_records
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .states import CORE_STATE_NAMES
from .task14_blind import build_model, routing, slug
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    attach_basal_regulation,
    sha256_file,
    sha256_object,
)


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/18_fixed_wt_inverse_rebalance"
ANALYSIS = REPO / "analysis/18_fixed_wt_inverse_rebalance"
BRANCH = "codex/task-18-fixed-wt-inverse-rebalance"
INHERITANCE = "fb2afcf3e23af1d19b297280ad67953370f7c9a9"
SCIENTIFIC_BASELINE = "4d4403f279fc36aec940c9e4759486d02f07e3f6"
INTAKE = "ef29cdeb80334ac44d7b1c2d10196bcac921d63f"
CONDITIONS = ("inherited_baseline", "share_0.10", "share_0.30")
TARGETS = {"share_0.10": 0.10, "share_0.30": 0.30}
CALCIUM = (0.10, 0.25, 0.50)
SPEC = WTCalibrationSpec()

# These are independent positive capacities.  In particular, pump and K
# membrane capacities are not represented as totals plus penalised fractions.
FREE_CAPACITIES = (
    "NHE1",
    "pump_apical",
    "pump_basolateral",
    "K_apical",
    "K_basolateral",
    "CaCC",
    "CO2_basolateral",
    "CO2_apical",
    "paracellular_Na",
    "paracellular_K",
    "paracellular_Cl",
    "paracellular_HCO3",
)
FORCED_CAPACITIES = ("AE4", "NKCC1")
OBJECTIVE_CAPACITIES = FORCED_CAPACITIES + FREE_CAPACITIES
SOURCE_ROWS = (
    *CORE_STATE_NAMES,
    "apical_current_fmol_charge_s",
    "basolateral_current_fmol_charge_s",
)
PREFREEZE_INPUTS = (
    "results/13B_modern_full_model/native_dynamic_contract_manifest.json",
    "results/13B_modern_full_model/native_dynamic_contract_profile.csv",
    "results/13B_modern_full_model/native_dynamic_contract_group_gate.csv",
    "results/13B_modern_full_model/native_dynamic_contract_solver_crosschecks.csv",
    "results/13B_modern_full_model/native_source_wt_roots.csv",
    "results/13B_modern_full_model/native_source_wt_summary.json",
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read(path: Path) -> Any:
    return json.loads(path.read_text())


def relative(path: Path) -> str:
    return str(path.relative_to(REPO))


def state_bytes_sha256(state: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(state, dtype=np.float64).tobytes()).hexdigest()


def baseline(manifest: Mapping[str, Any], root_id: str):
    model = build_model(manifest, root_id, 0.10)
    state = attach_basal_regulation(model, manifest["roots"][root_id]["core_state"])
    inherited_path = (
        REPO
        / "results/17_fixed_wt_flux_repartition/fixed_wt_payloads"
        / f"{slug(root_id)}__inherited_baseline.json"
    )
    inherited = read(inherited_path)
    if inherited["core_state"] != manifest["roots"][root_id]["core_state"]:
        raise AssertionError("Task 17 core state does not match the frozen WT manifest")
    if inherited["complete_state"] != state.tolist():
        raise AssertionError("Task 17 complete WT state was not copied exactly")
    if inherited["complete_state_bytes_sha256"] != state_bytes_sha256(state):
        raise AssertionError("Task 17 complete WT state byte hash changed")
    state.flags.writeable = False
    return model, state, model.evaluate(0.0, state, genotype=WT)


def _ae4_sources(evaluation) -> dict[str, float]:
    ae4 = evaluation.diagnostics.ae4
    return {
        "na": float(ae4.na_cell_fmol_s),
        "k": float(ae4.k_cell_fmol_s),
        "cl": float(ae4.cl_cell_fmol_s),
        "tic": float(ae4.tic_cell_fmol_s),
        "alkalinity": float(ae4.alkalinity_cell_fmol_s),
    }


def reference_capacities(model) -> dict[str, float]:
    hp = model.parameters.homeostasis
    mp = model.parameters.membranes
    values = {
        "AE4": model.ae4_parameters.carrier_amount_fmol,
        "NKCC1": hp.nkcc1_capacity_fmol_s,
        "NHE1": hp.nhe1_capacity_fmol_s,
        "pump_apical": mp.nak_capacity_fmol_s * mp.apical_pump_fraction,
        "pump_basolateral": mp.nak_capacity_fmol_s * (1.0 - mp.apical_pump_fraction),
        "K_apical": mp.g_k_total_S * mp.apical_k_fraction,
        "K_basolateral": mp.g_k_total_S * (1.0 - mp.apical_k_fraction),
        "CaCC": mp.g_cl_apical_S,
        "CO2_basolateral": hp.co2_basolateral_permeability_fmol_s_mM,
        "CO2_apical": hp.co2_apical_permeability_fmol_s_mM,
        "paracellular_Na": mp.g_para_na_S,
        "paracellular_K": mp.g_para_k_S,
        "paracellular_Cl": mp.g_para_cl_S,
        "paracellular_HCO3": mp.g_para_hco3_S,
    }
    if not all(math.isfinite(v) and v > 0.0 for v in values.values()):
        raise ValueError("Task 18 objective requires positive inherited capacities")
    return values


def allocation(evaluation, target: float | None) -> dict[str, Any]:
    ae4 = evaluation.diagnostics.ae4
    home = evaluation.diagnostics.homeostasis
    signed = {
        "AE4": float(ae4.cl_cell_fmol_s),
        "NKCC1": float(2.0 * home.nkcc1_inward_fmol_s),
        "AE2": float(home.ae2_inward_fmol_s),
    }
    positive = {name: value for name, value in signed.items() if value > 0.0}
    if set(positive) != {"AE4", "NKCC1"}:
        raise ValueError(f"Unexpected positive basolateral Cl loaders: {positive}")
    lplus = sum(positive.values())
    if target is None:
        lam = rho = 1.0
    else:
        if target not in TARGETS.values():
            raise ValueError("Only the predeclared 0.10 and 0.30 shares are allowed")
        lam = target * lplus / positive["AE4"]
        rho = (1.0 - target) * lplus / positive["NKCC1"]
    return {
        "target_share": target,
        "signed_basolateral_cl_fmol_s": signed,
        "positive_loader_names": tuple(sorted(positive)),
        "inherited_positive_loading_pool_fmol_s": lplus,
        "ae4_capacity_multiplier": lam,
        "other_positive_loader_common_multiplier": rho,
        "nkcc1_capacity_multiplier": rho,
    }


def fixed_requirements(model, evaluation, allocation_values: Mapping[str, Any]) -> dict[str, Any]:
    d = evaluation.diagnostics
    ae4 = _ae4_sources(evaluation)
    lam = float(allocation_values["ae4_capacity_multiplier"])
    rho = float(allocation_values["nkcc1_capacity_multiplier"])
    n = float(d.homeostasis.nkcc1_inward_fmol_s) * rho
    e = float(d.homeostasis.ae2_inward_fmol_s)
    out = {key: float(value) for key, value in d.outflow_sources_fmol_s.items()}

    nhe = e - lam * ae4["alkalinity"]
    co2_apical = out["alkalinity"] - out["tic"]
    co2_basolateral = e - lam * ae4["tic"] - co2_apical
    pump_total = (n + nhe + lam * ae4["na"]) / 3.0
    k_total = n + lam * ae4["k"] + 2.0 * pump_total
    cacc = 2.0 * n + e + lam * ae4["cl"]
    para_hco3 = -out["alkalinity"]
    para_cl = cacc - out["cl"]

    required = {
        "nkcc1_inward_fmol_s": n,
        "ae2_inward_fmol_s": e,
        "nhe1_inward_fmol_s": nhe,
        "co2_apical_lumen_to_cell_fmol_s": co2_apical,
        "co2_basolateral_bath_to_cell_fmol_s": co2_basolateral,
        "pump_total_cycles_fmol_s": pump_total,
        "k_total_cell_outward_fmol_s": k_total,
        "cacc_cl_cell_to_lumen_fmol_s": cacc,
        "paracellular_hco3_lumen_to_bath_fmol_s": para_hco3,
        "paracellular_cl_lumen_to_bath_fmol_s": para_cl,
        "outflow_fmol_s": out,
        "ae4_sources_fmol_s": {key: lam * value for key, value in ae4.items()},
    }
    finite = [value for value in required.values() if isinstance(value, float)]
    if not all(math.isfinite(value) for value in finite):
        raise ValueError("Nonfinite fixed state flux requirement")
    if min(nhe, co2_basolateral, pump_total, cacc) <= 0.0:
        raise ValueError("Existing topology has a sign incompatible neutral or pump requirement")
    return required


def fixed_environment(model, evaluation, requirements: Mapping[str, Any]) -> dict[str, float]:
    d = evaluation.diagnostics
    obs = d.observables
    ci, li = obs.cell_concentrations_mM, obs.lumen_concentrations_mM
    bath, mp = model.parameters.bath, model.parameters.membranes
    thermal = model.parameters.constants.thermal_voltage_V
    bath_hco3 = obs.bath_acid_base.hco3_mM

    def reversal(first: float, second: float, z: int) -> float:
        return nernst_voltage_V(first, second, valence=z, thermal_voltage_V=thermal)

    ref = reference_capacities(model)
    result = {
        "faraday_C_mol": model.parameters.constants.faraday_C_mol,
        "calcium_gate": hill_activation(
            d.stimulus.calcium_uM, mp.calcium_half_uM, mp.calcium_hill
        ),
        "E_K_apical_V": reversal(ci["k"], li["k"], +1),
        "E_K_basolateral_V": reversal(ci["k"], bath.k_mM, +1),
        "E_Cl_apical_V": reversal(ci["cl"], li["cl"], -1),
        "E_para_Na_V": reversal(li["na"], bath.na_mM, +1),
        "E_para_K_V": reversal(li["k"], bath.k_mM, +1),
        "E_para_Cl_V": reversal(li["cl"], bath.cl_mM, -1),
        "E_para_HCO3_V": reversal(li["hco3"], bath_hco3, -1),
        "pump_apical_flux_per_capacity": (
            d.membranes.pump_apical_fmol_s / ref["pump_apical"]
        ),
        "pump_basolateral_flux_per_capacity": (
            d.membranes.pump_basolateral_fmol_s / ref["pump_basolateral"]
        ),
        "nhe1_flux_per_capacity": d.homeostasis.nhe1_inward_fmol_s / ref["NHE1"],
        "co2_basolateral_drive_mM": (
            obs.bath_acid_base.co2_mM - ci["co2"]
        ),
        "co2_apical_drive_mM": li["co2"] - ci["co2"],
        "g_apical_background_S": mp.g_apical_background_S,
        "g_basolateral_background_S": mp.g_basolateral_background_S,
        "reference_v_apical_V": d.membranes.v_apical_V,
        "reference_v_basolateral_V": d.membranes.v_basolateral_V,
        "reference_v_transepithelial_V": d.membranes.v_transepithelial_V,
        "reference_pump_apical_fmol_s": d.membranes.pump_apical_fmol_s,
        "reference_k_apical_fmol_s": (
            d.membranes.currents_A["k_apical"]
            * 1.0e15 / model.parameters.constants.faraday_C_mol
        ),
    }
    if result["calcium_gate"] <= 0.0:
        raise ValueError("Positive basal calcium gate required")
    return result


def capacities_from_free(
    free: Sequence[float],
    model,
    requirements: Mapping[str, Any],
    environment: Mapping[str, float],
    allocation_values: Mapping[str, Any],
) -> tuple[dict[str, float], dict[str, float]]:
    """Map ``(pump fraction, scaled K flux, Va mV, Vt mV)`` to capacities."""

    q, k_scaled, va_mV, vt_mV = (float(value) for value in free)
    ref = reference_capacities(model)
    p_total = float(requirements["pump_total_cycles_fmol_s"])
    k_total = float(requirements["k_total_cell_outward_fmol_s"])
    k_scale = max(abs(k_total), abs(environment["reference_k_apical_fmol_s"]), 0.05)
    pa = q * p_total
    pb = (1.0 - q) * p_total
    ka = k_scaled * k_scale
    kb = k_total - ka
    out = requirements["outflow_fmol_s"]
    para_na = 3.0 * pa - out["na"]
    para_k = ka - 2.0 * pa - out["k"]
    para_cl = float(requirements["paracellular_cl_lumen_to_bath_fmol_s"])
    para_hco3 = float(requirements["paracellular_hco3_lumen_to_bath_fmol_s"])

    va = va_mV / 1000.0
    vt = vt_mV / 1000.0
    vb = va + vt
    faraday_scale = environment["faraday_C_mol"] / 1.0e15
    gate = environment["calcium_gate"]

    def ratio(numerator: float, denominator: float) -> float:
        if denominator == 0.0:
            return math.copysign(math.inf, numerator)
        return numerator / denominator

    gka_effective = ratio(faraday_scale * ka, va - environment["E_K_apical_V"])
    gkb_effective = ratio(faraday_scale * kb, vb - environment["E_K_basolateral_V"])
    capacities = {
        "AE4": ref["AE4"] * float(allocation_values["ae4_capacity_multiplier"]),
        "NKCC1": ref["NKCC1"] * float(allocation_values["nkcc1_capacity_multiplier"]),
        "NHE1": ratio(
            float(requirements["nhe1_inward_fmol_s"]),
            environment["nhe1_flux_per_capacity"],
        ),
        "pump_apical": ratio(pa, environment["pump_apical_flux_per_capacity"]),
        "pump_basolateral": ratio(pb, environment["pump_basolateral_flux_per_capacity"]),
        "K_apical": ratio(
            gka_effective - environment["g_apical_background_S"], gate
        ),
        "K_basolateral": ratio(
            gkb_effective - environment["g_basolateral_background_S"], gate
        ),
        "CaCC": ratio(
            -faraday_scale * float(requirements["cacc_cl_cell_to_lumen_fmol_s"]),
            gate * (va - environment["E_Cl_apical_V"]),
        ),
        "CO2_basolateral": ratio(
            float(requirements["co2_basolateral_bath_to_cell_fmol_s"]),
            environment["co2_basolateral_drive_mM"],
        ),
        "CO2_apical": ratio(
            float(requirements["co2_apical_lumen_to_cell_fmol_s"]),
            environment["co2_apical_drive_mM"],
        ),
        "paracellular_Na": ratio(
            faraday_scale * para_na, vt - environment["E_para_Na_V"]
        ),
        "paracellular_K": ratio(
            faraday_scale * para_k, vt - environment["E_para_K_V"]
        ),
        "paracellular_Cl": ratio(
            -faraday_scale * para_cl, vt - environment["E_para_Cl_V"]
        ),
        "paracellular_HCO3": ratio(
            -faraday_scale * para_hco3, vt - environment["E_para_HCO3_V"]
        ),
    }
    fluxes = {
        "pump_apical_fmol_s": pa,
        "pump_basolateral_fmol_s": pb,
        "K_apical_cell_to_lumen_fmol_s": ka,
        "K_basolateral_cell_to_bath_fmol_s": kb,
        "paracellular_Na_lumen_to_bath_fmol_s": para_na,
        "paracellular_K_lumen_to_bath_fmol_s": para_k,
        "paracellular_Cl_lumen_to_bath_fmol_s": para_cl,
        "paracellular_HCO3_lumen_to_bath_fmol_s": para_hco3,
        "v_apical_V": va,
        "v_basolateral_V": vb,
        "v_transepithelial_V": vt,
        "k_flux_scale_fmol_s": k_scale,
    }
    return capacities, fluxes


def model_with_capacities(model, capacities: Mapping[str, float]):
    if not all(
        name in capacities and math.isfinite(capacities[name]) and capacities[name] > 0.0
        for name in OBJECTIVE_CAPACITIES
    ):
        raise ValueError("All Task 18 capacity coordinates must be finite and positive")
    hp0, mp0 = model.parameters.homeostasis, model.parameters.membranes
    pump_total = capacities["pump_apical"] + capacities["pump_basolateral"]
    k_total = capacities["K_apical"] + capacities["K_basolateral"]
    hp = replace(
        hp0,
        nkcc1_capacity_fmol_s=capacities["NKCC1"],
        nhe1_capacity_fmol_s=capacities["NHE1"],
        co2_basolateral_permeability_fmol_s_mM=capacities["CO2_basolateral"],
        co2_apical_permeability_fmol_s_mM=capacities["CO2_apical"],
    )
    mp = replace(
        mp0,
        nak_capacity_fmol_s=pump_total,
        apical_pump_fraction=capacities["pump_apical"] / pump_total,
        g_k_total_S=k_total,
        apical_k_fraction=capacities["K_apical"] / k_total,
        g_cl_apical_S=capacities["CaCC"],
        g_para_na_S=capacities["paracellular_Na"],
        g_para_k_S=capacities["paracellular_K"],
        g_para_cl_S=capacities["paracellular_Cl"],
        g_para_hco3_S=capacities["paracellular_HCO3"],
    )
    ae4 = replace(model.ae4_parameters, carrier_amount_fmol=capacities["AE4"])
    base = ModernFullModel(
        parameters=replace(model.parameters, homeostasis=hp, membranes=mp),
        stimulus=model.stimulus,
        regulatory_model=model.regulatory_model,
        ae4_parameters=ae4,
        ae4_evaluator=model.ae4_evaluator,
    )
    return StimulatedNkcc1Model(base)


def log_changes(capacities: Mapping[str, float], reference: Mapping[str, float]) -> np.ndarray:
    return np.asarray(
        [math.log(capacities[name] / reference[name]) for name in OBJECTIVE_CAPACITIES],
        dtype=float,
    )


def _bounds(environment: Mapping[str, float]) -> tuple[tuple[float, float], ...]:
    # Fixed Cl and HCO3 paracellular fluxes impose this exact voltage interval.
    vt_low = 1000.0 * environment["E_para_HCO3_V"] + 1.0e-5
    vt_high = 1000.0 * environment["E_para_Cl_V"] - 1.0e-5
    if vt_low >= vt_high:
        raise ValueError("Paracellular Cl and HCO3 signs have no voltage overlap")
    va_high = 1000.0 * environment["E_Cl_apical_V"] - 1.0e-5
    return ((1.0e-8, 1.0 - 1.0e-8), (-20.0, 20.0), (-1000.0, va_high), (vt_low, vt_high))


def _valid_logs(
    free: Sequence[float], model, requirements, environment, allocation_values
) -> tuple[np.ndarray | None, dict[str, float] | None, dict[str, float] | None]:
    capacities, fluxes = capacities_from_free(
        free, model, requirements, environment, allocation_values
    )
    if not all(
        math.isfinite(capacities[name]) and capacities[name] > 0.0
        for name in OBJECTIVE_CAPACITIES
    ):
        return None, capacities, fluxes
    return log_changes(capacities, reference_capacities(model)), capacities, fluxes


def _random_feasible_seeds(
    model,
    requirements,
    environment,
    allocation_values,
    *,
    seed: int,
    draws: int = 120000,
) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    bounds = _bounds(environment)
    p_total = requirements["pump_total_cycles_fmol_s"]
    q_na = requirements["outflow_fmol_s"]["na"] / (3.0 * p_total)
    candidates: list[tuple[float, float, np.ndarray]] = []
    for index in range(draws):
        if index % 3 == 0:
            q = rng.uniform(bounds[0][0], min(bounds[0][1], max(bounds[0][0], q_na)))
        elif index % 3 == 1:
            q = min(bounds[0][1], max(bounds[0][0], q_na * rng.uniform(0.02, 2.0)))
        else:
            q = rng.beta(0.45, 2.2)
        if index % 2:
            k_scaled = rng.uniform(-3.0, 3.0)
        else:
            target_fraction = rng.uniform(-1.0, 2.0)
            k_scaled = (
                target_fraction * requirements["k_total_cell_outward_fmol_s"]
                / max(
                    abs(requirements["k_total_cell_outward_fmol_s"]),
                    abs(environment["reference_k_apical_fmol_s"]),
                    0.05,
                )
            )
        if requirements["k_total_cell_outward_fmol_s"] < 0.0 and index % 4 < 3:
            vt_hi = min(bounds[3][1], 1000.0 * environment["E_para_K_V"] - 1.0e-5)
            vt = rng.uniform(bounds[3][0], vt_hi) if vt_hi > bounds[3][0] else rng.uniform(*bounds[3])
            va_limit = min(
                bounds[2][1],
                1000.0 * environment["E_K_apical_V"] - 1.0e-5,
                1000.0 * (environment["E_K_basolateral_V"] - vt / 1000.0) - 1.0e-5,
            )
            va = rng.uniform(bounds[2][0], va_limit) if va_limit > bounds[2][0] else rng.uniform(*bounds[2])
        else:
            vt = rng.uniform(*bounds[3])
            va = rng.uniform(*bounds[2])
        free = np.asarray((q, k_scaled, va, vt), dtype=float)
        logs, _, _ = _valid_logs(
            free, model, requirements, environment, allocation_values
        )
        if logs is not None:
            candidates.append((float(np.max(np.abs(logs))), float(logs @ logs), free))
    candidates.sort(key=lambda item: (item[0], item[1]))
    return [item[2] for item in candidates[:32]]


def solve_canonical(model, evaluation, allocation_values: Mapping[str, Any], *, seed: int):
    reference = reference_capacities(model)
    requirements = fixed_requirements(model, evaluation, allocation_values)
    environment = fixed_environment(model, evaluation, requirements)
    target = allocation_values["target_share"]
    if target is None:
        membrane = evaluation.diagnostics.membranes
        convert = 1.0e15 / model.parameters.constants.faraday_C_mol
        k_scale = max(
            abs(requirements["k_total_cell_outward_fmol_s"]),
            abs(environment["reference_k_apical_fmol_s"]),
            0.05,
        )
        free = np.asarray(
            (
                membrane.pump_apical_fmol_s / requirements["pump_total_cycles_fmol_s"],
                environment["reference_k_apical_fmol_s"] / k_scale,
                1000.0 * membrane.v_apical_V,
                1000.0 * membrane.v_transepithelial_V,
            )
        )
        fluxes = {
            "pump_apical_fmol_s": membrane.pump_apical_fmol_s,
            "pump_basolateral_fmol_s": membrane.pump_basolateral_fmol_s,
            "K_apical_cell_to_lumen_fmol_s": membrane.currents_A["k_apical"] * convert,
            "K_basolateral_cell_to_bath_fmol_s": membrane.currents_A["k_basolateral"] * convert,
            "paracellular_Na_lumen_to_bath_fmol_s": membrane.currents_A["para_na"] * convert,
            "paracellular_K_lumen_to_bath_fmol_s": membrane.currents_A["para_k"] * convert,
            "paracellular_Cl_lumen_to_bath_fmol_s": -membrane.currents_A["para_cl"] * convert,
            "paracellular_HCO3_lumen_to_bath_fmol_s": -membrane.currents_A["para_hco3"] * convert,
            "v_apical_V": membrane.v_apical_V,
            "v_basolateral_V": membrane.v_basolateral_V,
            "v_transepithelial_V": membrane.v_transepithelial_V,
            "k_flux_scale_fmol_s": k_scale,
        }
        return {
            "free": free,
            "capacities": reference,
            "fluxes": fluxes,
            "requirements": requirements,
            "environment": environment,
            "stage1": {"success": True, "method": "inherited_exact", "max_abs_log_fold": 0.0},
            "stage2": {"success": True, "method": "inherited_exact", "sum_squared_log_fold": 0.0},
            "search": {"random_draws": 0, "feasible_seed_count": 1},
        }

    bounds = _bounds(environment)

    def de_objective(free):
        logs, capacities, _ = _valid_logs(
            free, model, requirements, environment, allocation_values
        )
        if logs is not None:
            return float(np.max(np.abs(logs)))
        penalty = 100.0
        assert capacities is not None
        for name in OBJECTIVE_CAPACITIES:
            ratio = capacities[name] / reference[name]
            if not math.isfinite(ratio):
                penalty += 100.0
            elif ratio <= 0.0:
                penalty += min(100.0, 1.0 + abs(ratio))
        return penalty

    seeds = _random_feasible_seeds(
        model, requirements, environment, allocation_values, seed=seed
    )
    de = differential_evolution(
        de_objective,
        bounds,
        seed=seed + 1009,
        popsize=20,
        maxiter=350,
        tol=1.0e-10,
        atol=1.0e-10,
        polish=False,
        updating="immediate",
        workers=1,
    )
    de_logs, _, _ = _valid_logs(
        de.x, model, requirements, environment, allocation_values
    )
    if de_logs is not None:
        seeds.append(np.asarray(de.x, dtype=float))
    if not seeds:
        raise RuntimeError("No physically positive fixed state parameterisation found")

    def epigraph_constraints(z):
        logs, _, _ = _valid_logs(
            z[:4], model, requirements, environment, allocation_values
        )
        if logs is None:
            return np.full(2 * len(OBJECTIVE_CAPACITIES), -1.0e6)
        return np.r_[z[4] - logs, z[4] + logs]

    stage1_results = []
    for free0 in seeds[:20]:
        logs0, _, _ = _valid_logs(
            free0, model, requirements, environment, allocation_values
        )
        if logs0 is None:
            continue
        z0 = np.r_[free0, np.max(np.abs(logs0)) + 1.0e-8]
        result = minimize(
            lambda z: z[4],
            z0,
            method="SLSQP",
            bounds=(*bounds, (0.0, 50.0)),
            constraints={"type": "ineq", "fun": epigraph_constraints},
            options={"maxiter": 4000, "ftol": 2.0e-13, "disp": False},
        )
        logs, _, _ = _valid_logs(
            result.x[:4], model, requirements, environment, allocation_values
        )
        violation = (
            math.inf if logs is None else max(0.0, float(np.max(np.abs(logs))) - result.x[4])
        )
        if logs is not None and violation <= 2.0e-7:
            stage1_results.append((float(np.max(np.abs(logs))), result, logs))
    if not stage1_results:
        raise RuntimeError("Minimax capacity optimisation did not converge")
    stage1_results.sort(key=lambda item: (item[0], float(item[2] @ item[2])))
    t_star, first, first_logs = stage1_results[0]

    t_limit = t_star + 2.0e-8

    def stage2_constraint(free):
        logs, _, _ = _valid_logs(
            free, model, requirements, environment, allocation_values
        )
        if logs is None:
            return np.full(2 * len(OBJECTIVE_CAPACITIES), -1.0e6)
        return np.r_[t_limit - logs, t_limit + logs]

    stage2_results = []
    stage2_starts = [item[1].x[:4] for item in stage1_results[:12]]
    for free0 in stage2_starts:
        result = minimize(
            lambda free: (
                1.0e6
                if _valid_logs(free, model, requirements, environment, allocation_values)[0] is None
                else float(
                    _valid_logs(free, model, requirements, environment, allocation_values)[0]
                    @ _valid_logs(free, model, requirements, environment, allocation_values)[0]
                )
            ),
            free0,
            method="SLSQP",
            bounds=bounds,
            constraints={"type": "ineq", "fun": stage2_constraint},
            options={"maxiter": 5000, "ftol": 2.0e-13, "disp": False},
        )
        logs, capacities, fluxes = _valid_logs(
            result.x, model, requirements, environment, allocation_values
        )
        if logs is not None and np.max(np.abs(logs)) <= t_limit + 3.0e-7:
            stage2_results.append(
                (not bool(result.success), float(logs @ logs), result, logs, capacities, fluxes)
            )
    if not stage2_results:
        raise RuntimeError("Lexicographic squared log fold optimisation did not converge")
    stage2_results.sort(key=lambda item: (item[0], item[1], float(np.max(np.abs(item[3])))))
    failed_flag, sumsq, second, logs, capacities, fluxes = stage2_results[0]
    assert capacities is not None and fluxes is not None
    return {
        "free": np.asarray(second.x, dtype=float),
        "capacities": capacities,
        "fluxes": fluxes,
        "requirements": requirements,
        "environment": environment,
        "stage1": {
            "success": not failed_flag,
            "accepted_feasible_limit_point": True,
            "method": "differential_evolution plus multistart SLSQP epigraph",
            "max_abs_log_fold": t_star,
            "reported_solution_max_abs_log_fold": float(np.max(np.abs(logs))),
            "independent_successful_refinements": len(stage1_results),
            "message": str(first.message),
        },
        "stage2": {
            "success": True,
            "method": "multistart SLSQP within minimax face",
            "sum_squared_log_fold": sumsq,
            "successful_refinements": len(stage2_results),
            "message": str(second.message),
        },
        "search": {
            "random_draws": 120000,
            "feasible_seed_count": len(seeds),
            "differential_evolution_success": bool(de.success),
            "differential_evolution_message": str(de.message),
            "bounds_role": "numerical search box only; no inherited capacity or voltage value is a feasibility wall",
            "solution_distance_to_search_bounds": [
                min(second.x[i] - bounds[i][0], bounds[i][1] - second.x[i])
                for i in range(4)
            ],
        },
    }


def residual_summary(evaluation) -> dict[str, Any]:
    independent = evaluation.rhs[list(CORE_INDEPENDENT_ROWS)] / SPEC.independent_rhs_scales
    conservation = dict(evaluation.diagnostics.conservation_residuals)
    return {
        "max_abs_complete_rhs": float(np.max(np.abs(evaluation.rhs))),
        "max_abs_scaled_independent_rhs": float(np.max(np.abs(independent))),
        "max_abs_omitted_rhs": float(np.max(np.abs(evaluation.rhs[[4, 10]]))),
        "conservation_residuals": conservation,
        "max_conservation_tolerance_ratio": max(
            abs(value) / CONSERVATION_RESIDUAL_TOLERANCES[name]
            for name, value in conservation.items()
        ),
        "max_abs_state_charge_fmol": max(
            abs(value) for value in evaluation.diagnostics.state_charge_fmol.values()
        ),
        "max_abs_current_A": max(
            abs(value)
            for value in evaluation.diagnostics.membranes.current_residuals_A.values()
        ),
    }


def chloride_ledger(evaluation) -> dict[str, float]:
    d = evaluation.diagnostics
    faraday = d.membranes.currents_A
    convert = 1.0e15 / 96485.33212
    return {
        "AE4": float(d.ae4.cl_cell_fmol_s),
        "NKCC1": float(2.0 * d.homeostasis.nkcc1_inward_fmol_s),
        "AE2": float(d.homeostasis.ae2_inward_fmol_s),
        "CaCC_cell_source": float(faraday["cl_apical"] * convert),
        "CaCC_cell_to_lumen": float(-faraday["cl_apical"] * convert),
        "paracellular_Cl_lumen_to_bath": float(-faraday["para_cl"] * convert),
        "lumen_outflow_Cl": float(d.outflow_sources_fmol_s["cl"]),
    }


def verify_solution(
    model,
    state: np.ndarray,
    baseline_evaluation,
    solution: Mapping[str, Any],
    allocation_values: Mapping[str, Any],
) -> tuple[Any, dict[str, Any]]:
    modified = model_with_capacities(model, solution["capacities"])
    state_before = state.tobytes()
    evaluation = modified.evaluate(0.0, state, genotype=WT)
    if state.tobytes() != state_before:
        raise AssertionError("Fixed state bytes changed during production evaluation")
    summary = residual_summary(evaluation)
    if summary["max_abs_scaled_independent_rhs"] > SPEC.root_scaled_tolerance:
        raise AssertionError(summary)
    if summary["max_abs_omitted_rhs"] > SPEC.omitted_row_raw_tolerance:
        raise AssertionError(summary)
    if summary["max_conservation_tolerance_ratio"] > 1.0:
        raise AssertionError(summary)
    ledger = chloride_ledger(evaluation)
    positive = sum(value for name, value in ledger.items() if name in {"AE4", "NKCC1"} and value > 0)
    inherited_pool = allocation_values["inherited_positive_loading_pool_fmol_s"]
    realized = ledger["AE4"] / positive
    expected = allocation_values["target_share"]
    if expected is None:
        expected = chloride_ledger(baseline_evaluation)["AE4"] / inherited_pool
    if not math.isclose(positive, inherited_pool, rel_tol=2e-11, abs_tol=2e-13):
        raise AssertionError("Positive chloride loading pool changed")
    if not math.isclose(realized, expected, rel_tol=2e-11, abs_tol=2e-13):
        raise AssertionError("AE4 chloride share was not enforced")
    voltage = evaluation.diagnostics.membranes
    target_voltage = solution["fluxes"]
    for actual, key in (
        (voltage.v_apical_V, "v_apical_V"),
        (voltage.v_basolateral_V, "v_basolateral_V"),
        (voltage.v_transepithelial_V, "v_transepithelial_V"),
    ):
        if not math.isclose(actual, target_voltage[key], rel_tol=1e-9, abs_tol=2e-11):
            raise AssertionError(f"Electrical reclosure disagrees for {key}")
    return modified, {
        **summary,
        "realized_ae4_positive_cl_share": realized,
        "realized_positive_loading_pool_fmol_s": positive,
        "chloride_ledger": ledger,
        "voltages_V": {
            "apical": voltage.v_apical_V,
            "basolateral": voltage.v_basolateral_V,
            "transepithelial": voltage.v_transepithelial_V,
        },
    }


def _rank_payload(model, state, evaluation) -> dict[str, Any]:
    """Closed and explicit voltage fixed state sensitivity audits."""

    reference = reference_capacities(model)
    names = FREE_CAPACITIES
    base = evaluation.rhs[:12]
    columns = []
    step = 1.0e-5
    for name in names:
        plus = dict(reference)
        minus = dict(reference)
        plus[name] *= math.exp(step)
        minus[name] *= math.exp(-step)
        ep = model_with_capacities(model, plus).evaluate(0.0, state, genotype=WT)
        em = model_with_capacities(model, minus).evaluate(0.0, state, genotype=WT)
        columns.append((ep.rhs[:12] - em.rhs[:12]) / (2.0 * step))
    matrix = np.column_stack(columns)
    row_scale = np.maximum(np.max(np.abs(matrix), axis=1), 1.0e-14)
    normal = matrix / row_scale[:, None]
    col_scale = np.maximum(np.linalg.norm(normal, axis=0), 1.0e-14)
    normal = normal / col_scale
    u, singular, vh = np.linalg.svd(normal, full_matrices=True)
    threshold = SPEC.jacobian_rank_relative_tolerance * singular[0]
    rank = int(np.sum(singular > threshold))
    dependent = []
    for index in range(rank, u.shape[1]):
        vector = u[:, index]
        terms = sorted(
            ((CORE_STATE_NAMES[i], float(vector[i])) for i in range(len(vector))),
            key=lambda item: abs(item[1]),
            reverse=True,
        )[:5]
        dependent.append(terms)
    closed_payload = {
        "capacity_names": names,
        "complete_rhs_rows": CORE_STATE_NAMES,
        "matrix": matrix.tolist(),
        "row_scale": row_scale.tolist(),
        "column_scale": col_scale.tolist(),
        "singular_values": singular.tolist(),
        "rank": rank,
        "nullity": len(names) - rank,
        "left_nullity": len(CORE_STATE_NAMES) - rank,
        "structurally_dependent_row_combinations": dependent,
        "baseline_max_abs_rhs": float(np.max(np.abs(base))),
    }

    environment = fixed_environment(
        model, evaluation, fixed_requirements(model, evaluation, allocation(evaluation, None))
    )
    x0 = np.r_[
        np.zeros(len(FREE_CAPACITIES)),
        1000.0 * evaluation.diagnostics.membranes.v_apical_V,
        1000.0 * evaluation.diagnostics.membranes.v_basolateral_V,
    ]
    voltage_step_mV = 1.0e-3
    augmented_columns = []
    for index in range(len(x0)):
        local_step = step if index < len(FREE_CAPACITIES) else voltage_step_mV
        plus = x0.copy()
        minus = x0.copy()
        plus[index] += local_step
        minus[index] -= local_step
        augmented_columns.append(
            (
                _open_voltage_residual(model, evaluation, environment, plus)
                - _open_voltage_residual(model, evaluation, environment, minus)
            )
            / (2.0 * local_step)
        )
    augmented = np.column_stack(augmented_columns)
    aug_row_scale = np.maximum(np.max(np.abs(augmented), axis=1), 1.0e-14)
    aug_normal = augmented / aug_row_scale[:, None]
    aug_col_scale = np.maximum(np.linalg.norm(aug_normal, axis=0), 1.0e-14)
    aug_normal = aug_normal / aug_col_scale
    aug_u, aug_singular, aug_vh = np.linalg.svd(aug_normal, full_matrices=True)
    aug_threshold = SPEC.jacobian_rank_relative_tolerance * aug_singular[0]
    aug_rank = int(np.sum(aug_singular > aug_threshold))
    aug_dependent = []
    for index in range(aug_rank, aug_u.shape[1]):
        vector = aug_u[:, index]
        terms = sorted(
            ((SOURCE_ROWS[i], float(vector[i])) for i in range(len(vector))),
            key=lambda item: abs(item[1]),
            reverse=True,
        )[:6]
        aug_dependent.append(terms)
    baseline_open = _open_voltage_residual(model, evaluation, environment, x0)
    return {
        "rank": aug_rank,
        "nullity": len(x0) - aug_rank,
        "left_nullity": len(SOURCE_ROWS) - aug_rank,
        "explicit_voltage_augmented": {
            "variable_names": (*FREE_CAPACITIES, "v_apical_mV", "v_basolateral_mV"),
            "equation_names": SOURCE_ROWS,
            "matrix": augmented.tolist(),
            "row_scale": aug_row_scale.tolist(),
            "column_scale": aug_col_scale.tolist(),
            "singular_values": aug_singular.tolist(),
            "rank": aug_rank,
            "nullity": len(x0) - aug_rank,
            "left_nullity": len(SOURCE_ROWS) - aug_rank,
            "structurally_dependent_row_combinations": aug_dependent,
            "baseline_max_abs_residual": float(np.max(np.abs(baseline_open))),
        },
        "closed_voltage_capacity_sensitivity": closed_payload,
        "electrical_closure": {
            "voltage_variables": ["v_apical_V", "v_basolateral_V"],
            "closure_rows": ["apical_current_fmol_charge_s", "basolateral_current_fmol_charge_s"],
            "treatment": "The augmented audit retains both voltages and both closure equations explicitly. The production sensitivity separately solves the same two by two closure exactly.",
            "baseline_current_residuals_A": dict(
                evaluation.diagnostics.membranes.current_residuals_A
            ),
            "baseline_voltages_V": {
                "apical": evaluation.diagnostics.membranes.v_apical_V,
                "basolateral": evaluation.diagnostics.membranes.v_basolateral_V,
                "transepithelial": evaluation.diagnostics.membranes.v_transepithelial_V,
            },
        },
        "baseline_max_abs_rhs": float(np.max(np.abs(base))),
    }


def _open_voltage_residual(model, evaluation, environment, variables) -> np.ndarray:
    """Material and current residuals before eliminating the two voltages."""

    logs = dict(zip(FREE_CAPACITIES, np.asarray(variables[: len(FREE_CAPACITIES)])))
    multiplier = {name: math.exp(float(logs[name])) for name in FREE_CAPACITIES}
    va = float(variables[-2]) / 1000.0
    vb = float(variables[-1]) / 1000.0
    vt = vb - va
    d = evaluation.diagnostics
    home, ae4, water = d.homeostasis, d.ae4, d.water
    out = d.outflow_sources_fmol_s
    faraday = environment["faraday_C_mol"]
    to_fmol = 1.0e15 / faraday
    ref = reference_capacities(model)

    n = home.nkcc1_inward_fmol_s
    e = home.ae2_inward_fmol_s
    h = home.nhe1_inward_fmol_s * multiplier["NHE1"]
    pa = d.membranes.pump_apical_fmol_s * multiplier["pump_apical"]
    pb = d.membranes.pump_basolateral_fmol_s * multiplier["pump_basolateral"]
    gate = environment["calcium_gate"]
    gka = gate * ref["K_apical"] * multiplier["K_apical"] + environment[
        "g_apical_background_S"
    ]
    gkb = gate * ref["K_basolateral"] * multiplier["K_basolateral"] + environment[
        "g_basolateral_background_S"
    ]
    gcl = gate * ref["CaCC"] * multiplier["CaCC"]
    ika = gka * (va - environment["E_K_apical_V"])
    ikb = gkb * (vb - environment["E_K_basolateral_V"])
    icl = gcl * (va - environment["E_Cl_apical_V"])
    jka = ika * to_fmol
    jkb = ikb * to_fmol
    jcl = -icl * to_fmol

    para = {}
    para_current = {}
    for name, capacity_name, reversal, valence in (
        ("na", "paracellular_Na", "E_para_Na_V", +1),
        ("k", "paracellular_K", "E_para_K_V", +1),
        ("cl", "paracellular_Cl", "E_para_Cl_V", -1),
        ("hco3", "paracellular_HCO3", "E_para_HCO3_V", -1),
    ):
        current = ref[capacity_name] * multiplier[capacity_name] * (
            vt - environment[reversal]
        )
        para_current[name] = current
        para[name] = current * to_fmol / valence
    co2b = d.co2_fluxes_fmol_s["bath_to_cell"] * multiplier["CO2_basolateral"]
    co2a = d.co2_fluxes_fmol_s["lumen_to_cell"] * multiplier["CO2_apical"]

    residual = np.zeros(14)
    residual[0] = n + h + ae4.na_cell_fmol_s - 3.0 * (pa + pb)
    residual[1] = n + ae4.k_cell_fmol_s + 2.0 * (pa + pb) - jka - jkb
    residual[2] = 2.0 * n + e + ae4.cl_cell_fmol_s - jcl
    residual[3] = -e + ae4.tic_cell_fmol_s + co2b + co2a
    residual[4] = h - e + ae4.alkalinity_cell_fmol_s
    residual[5] = water.cell_volume_source_pL_s
    residual[6] = 3.0 * pa - para["na"] - out["na"]
    residual[7] = -2.0 * pa + jka - para["k"] - out["k"]
    residual[8] = jcl - para["cl"] - out["cl"]
    residual[9] = -para["hco3"] - co2a - out["tic"]
    residual[10] = -para["hco3"] - out["alkalinity"]
    residual[11] = water.lumen_volume_source_pL_s
    pump_current_a = faraday * pa / 1.0e15
    pump_current_b = faraday * pb / 1.0e15
    i_apical = ika + icl + pump_current_a
    i_basolateral = ikb + pump_current_b
    i_para = sum(para_current.values())
    residual[12] = (i_apical - i_para) * to_fmol
    residual[13] = (i_basolateral + i_para) * to_fmol
    return residual


def prepare() -> None:
    if (OUT / "contract.json").exists():
        raise FileExistsError("Task 18 contract is immutable")
    manifest, verification = load_freeze()
    input_names = list(PREFREEZE_INPUTS) + [
        str(
            Path("results/17_fixed_wt_flux_repartition/fixed_wt_payloads")
            / f"{slug(root_id)}__inherited_baseline.json"
        )
        for root_id in sorted(manifest["roots"])
    ]
    input_hashes = {name: sha256_file(REPO / name) for name in input_names}
    if any("ae4_5pct" in name.lower() or "ae2_results" in name.lower() for name in input_hashes):
        raise AssertionError("Phenotype firewall input violation")

    provenance_rows = []
    ledger_rows = []
    ranks = {}
    for root_id in sorted(manifest["roots"]):
        model, state, evaluation = baseline(manifest, root_id)
        ref = reference_capacities(model)
        ranks[root_id] = _rank_payload(model, state, evaluation)
        allocation_values = allocation(evaluation, None)
        cl = chloride_ledger(evaluation)
        for name, value in cl.items():
            ledger_rows.append(
                {
                    "root_id": root_id,
                    "routing_family": routing(manifest, root_id),
                    "pathway": name,
                    "signed_cl_fmol_s": value,
                    "positive_basolateral_loader": name in allocation_values["positive_loader_names"],
                    "positive_loading_denominator_fmol_s": allocation_values[
                        "inherited_positive_loading_pool_fmol_s"
                    ],
                }
            )
        records = {record.name: record for record in parameter_records(model.parameters)}
        for name, value in ref.items():
            source_name = {
                "NKCC1": "homeostasis.nkcc1_capacity_fmol_s",
                "NHE1": "homeostasis.nhe1_capacity_fmol_s",
                "CaCC": "membranes.g_cl_apical_S",
                "CO2_basolateral": "homeostasis.co2_basolateral_permeability_fmol_s_mM",
                "CO2_apical": "homeostasis.co2_apical_permeability_fmol_s_mM",
                "paracellular_Na": "membranes.g_para_na_S",
                "paracellular_K": "membranes.g_para_k_S",
                "paracellular_Cl": "membranes.g_para_cl_S",
                "paracellular_HCO3": "membranes.g_para_hco3_S",
            }.get(name)
            record = records.get(source_name) if source_name else None
            provenance_rows.append(
                {
                    "root_id": root_id,
                    "capacity_coordinate": name,
                    "inherited_value": value,
                    "units": record.units if record else (
                        "fmol" if name == "AE4" else "independent membrane capacity"
                    ),
                    "provenance": record.provenance if record else "DERIVED_INDEPENDENT_COORDINATE",
                    "source_field": source_name or (
                        "ae4_parameters.carrier_amount_fmol" if name == "AE4" else "total plus membrane fraction"
                    ),
                    "task18_role": "TARGET_DETERMINED" if name in FORCED_CAPACITIES else "ADJUSTABLE",
                    "hard_inherited_bound": False,
                    "direct_same_system_measurement_interval": False,
                    "note": record.note if record else "Independent positive coordinate used by the canonical objective.",
                }
            )
        extra = (
            ("homeostasis.ae2_capacity_fmol_s", "FIXED_NONLOADING_TRANSPORTER"),
            ("membranes.g_apical_background_S", "BASELINE_ZERO_NOT_ACTIVATED"),
            ("membranes.g_basolateral_background_S", "BASELINE_ZERO_NOT_ACTIVATED"),
            ("water.apical_hydraulic_pL_s_mOsm", "AUDITED_FIXED_ZERO_VOLUME_PERTURBATION"),
            ("water.basolateral_hydraulic_pL_s_mOsm", "AUDITED_FIXED_ZERO_VOLUME_PERTURBATION"),
            ("water.paracellular_hydraulic_pL_s_mOsm", "AUDITED_FIXED_ZERO_VOLUME_PERTURBATION"),
            ("water.outflow_rate_s", "AUDITED_FIXED_ZERO_VOLUME_PERTURBATION"),
        )
        for source_name, role in extra:
            record = records[source_name]
            provenance_rows.append(
                {
                    "root_id": root_id,
                    "capacity_coordinate": source_name,
                    "inherited_value": record.value,
                    "units": record.units,
                    "provenance": record.provenance,
                    "source_field": source_name,
                    "task18_role": role,
                    "hard_inherited_bound": False,
                    "direct_same_system_measurement_interval": False,
                    "note": record.note,
                }
            )

    write_rows(OUT / "parameter_provenance.csv", provenance_rows)
    write_rows(OUT / "baseline_flux_ledger.csv", ledger_rows)
    write_json(OUT / "source_rank.json", ranks)
    contract = {
        "task": "18_FIXED_WT_INVERSE_REBALANCE",
        "created_utc": now(),
        "branch": BRANCH,
        "inheritance_commit": INHERITANCE,
        "scientific_baseline_commit": SCIENTIFIC_BASELINE,
        "intake_commit": INTAKE,
        "conditions": CONDITIONS,
        "target_shares": [0.10, 0.30],
        "calcium_uM": CALCIUM,
        "root_ids": sorted(manifest["roots"]),
        "fixed_state_rule": "Copy every inherited conserved amount, both volumes and basal regulatory coordinate exactly; no state solve, fitting, continuation, rounding or reconstruction.",
        "chloride_rule": "Scale the complete AE4 carrier to the target share and every other positive basolateral Cl loader by one proportional factor while preserving inherited Lplus.",
        "positive_loader_audit": "Signed production ledger identifies AE4 and NKCC1 only; AE2 is a negative counterflux at all ten roots.",
        "hard_constraints": [
            "fixed conserved WT state",
            "complete production steady state equations",
            "target AE4 positive Cl loading share",
            "preserved total positive basolateral Cl loading pool",
            "existing stoichiometry, topology, thermodynamics, bath, chemistry and geometry",
            "nonnegative physical capacities, membrane fractions in [0,1], finite quantities",
            "charge, current, carbon, alkalinity, water and lumen accounting",
        ],
        "reference_not_bounds": [
            "14 nS total K conductance",
            "31.4 nS CaCC conductance",
            "old NHE1 and AE4 boxes",
            "inherited pump, paracellular and CO2 capacities",
            "historical voltage targets",
        ],
        "objective_capacity_coordinates": OBJECTIVE_CAPACITIES,
        "canonical_objective": [
            "minimise largest absolute natural log fold change",
            "subject to that optimum minimise sum of squared natural log fold changes",
        ],
        "independent_membrane_coordinates": True,
        "genotype_information_in_objective": False,
        "genotype_information_read_before_freeze": False,
        "prefreeze_input_hashes": input_hashes,
        "intake_verification": verification,
        "source_rank_sha256": sha256_file(OUT / "source_rank.json"),
        "parameter_provenance_sha256": sha256_file(OUT / "parameter_provenance.csv"),
        "baseline_flux_ledger_sha256": sha256_file(OUT / "baseline_flux_ledger.csv"),
    }
    write_json(OUT / "contract.json", contract)
    print(json.dumps({"roots": len(ranks), "target_solves": 0, "genotype_reads": 0}))


def solve() -> None:
    if (OUT / "frozen_wt_manifest.json").exists():
        raise FileExistsError("Frozen WT manifest already exists")
    contract = read(OUT / "contract.json")
    if tuple(contract["conditions"]) != CONDITIONS:
        raise AssertionError("Condition contract changed")
    manifest, _ = load_freeze()
    rows = []
    for root_index, root_id in enumerate(sorted(manifest["roots"])):
        model, state, baseline_evaluation = baseline(manifest, root_id)
        state_digest = state_bytes_sha256(state)
        for condition_index, condition in enumerate(CONDITIONS):
            allocation_values = allocation(baseline_evaluation, TARGETS.get(condition))
            solution = solve_canonical(
                model,
                baseline_evaluation,
                allocation_values,
                seed=180000 + 101 * root_index + condition_index,
            )
            modified, verification = verify_solution(
                model, state, baseline_evaluation, solution, allocation_values
            )
            capacities = solution["capacities"]
            reference = reference_capacities(model)
            logs = log_changes(capacities, reference)
            candidate_id = slug(root_id) + "__" + condition
            membrane = modified.parameters.membranes
            ae4 = modified.evaluate(0.0, state).diagnostics.ae4
            turnover = abs(ae4.diagnostics["j_na_fmol_s"]) + abs(
                ae4.diagnostics["j_k_fmol_s"]
            )
            row = {
                "candidate_id": candidate_id,
                "root_id": root_id,
                "routing_family": routing(manifest, root_id),
                "condition": condition,
                "target_share": allocation_values["target_share"],
                "realized_ae4_share": verification["realized_ae4_positive_cl_share"],
                "inherited_positive_loading_pool_fmol_s": allocation_values[
                    "inherited_positive_loading_pool_fmol_s"
                ],
                "positive_loading_pool_change_fmol_s": verification[
                    "realized_positive_loading_pool_fmol_s"
                ] - allocation_values["inherited_positive_loading_pool_fmol_s"],
                "ae4_capacity_multiplier": allocation_values["ae4_capacity_multiplier"],
                "other_positive_loader_common_multiplier": allocation_values[
                    "other_positive_loader_common_multiplier"
                ],
                "fixed_wt_feasible": True,
                "status": "FIXED_WT_EXACT_REBALANCE_FEASIBLE",
                "state_coordinates_changed": 0,
                "exact_inherited_state_bytes_retained": True,
                "complete_state_bytes_sha256": state_digest,
                "largest_abs_log_fold": float(np.max(np.abs(logs))),
                "rms_log_fold": float(np.sqrt(np.mean(logs**2))),
                "sum_squared_log_fold": float(logs @ logs),
                "largest_fold_multiplier": float(math.exp(np.max(np.abs(logs)))),
                "g_k_total_S": membrane.g_k_total_S,
                "g_k_total_nS": membrane.g_k_total_S * 1.0e9,
                "g_cl_apical_S": membrane.g_cl_apical_S,
                "g_cl_apical_nS": membrane.g_cl_apical_S * 1.0e9,
                "apical_k_fraction": membrane.apical_k_fraction,
                "nak_capacity_fmol_s": membrane.nak_capacity_fmol_s,
                "apical_pump_fraction": membrane.apical_pump_fraction,
                "v_apical_mV": 1000.0 * verification["voltages_V"]["apical"],
                "v_basolateral_mV": 1000.0 * verification["voltages_V"]["basolateral"],
                "v_transepithelial_mV": 1000.0 * verification["voltages_V"]["transepithelial"],
                "max_abs_complete_rhs": verification["max_abs_complete_rhs"],
                "max_abs_scaled_independent_rhs": verification[
                    "max_abs_scaled_independent_rhs"
                ],
                "max_abs_omitted_rhs": verification["max_abs_omitted_rhs"],
                "max_conservation_tolerance_ratio": verification[
                    "max_conservation_tolerance_ratio"
                ],
                "max_abs_current_A": verification["max_abs_current_A"],
                "ae4_productive_fraction": abs(ae4.cl_cell_fmol_s) / turnover,
                "ae4_cancellation_fraction": 1.0 - abs(ae4.cl_cell_fmol_s) / turnover,
                **{
                    f"{name}_multiplier": capacities[name] / reference[name]
                    for name in OBJECTIVE_CAPACITIES
                },
            }
            payload = {
                "row": row,
                "contract_sha256": sha256_file(OUT / "contract.json"),
                "core_state": manifest["roots"][root_id]["core_state"],
                "core_state_sha256": manifest["roots"][root_id]["core_state_sha256"],
                "complete_state": state.tolist(),
                "complete_state_bytes_sha256": state_digest,
                "baseline_whole_cell_parameters": asdict(model.parameters),
                "baseline_ae4_parameters": asdict(model.ae4_parameters),
                "accepted_whole_cell_parameters": asdict(modified.parameters),
                "accepted_ae4_parameters": asdict(modified.ae4_parameters),
                "regulatory_model": asdict(modified.regulatory_model),
                "stimulus": asdict(modified.stimulus),
                "allocation": allocation_values,
                "requirements": solution["requirements"],
                "accepted_capacities": capacities,
                "reference_capacities": reference,
                "capacity_multipliers": {
                    name: capacities[name] / reference[name]
                    for name in OBJECTIVE_CAPACITIES
                },
                "free_solution": solution["free"].tolist(),
                "flux_solution": solution["fluxes"],
                "environment": solution["environment"],
                "stage1_optimisation": solution["stage1"],
                "stage2_optimisation": solution["stage2"],
                "search_diagnostics": solution["search"],
                "baseline_residuals": residual_summary(baseline_evaluation),
                "accepted_residuals": verification,
                "phenotype_information_used": False,
            }
            write_json(OUT / "wt_parameter_payloads" / f"{candidate_id}.json", payload)
            rows.append(row)
            print(json.dumps({"candidate": candidate_id, "max_log": row["largest_abs_log_fold"]}))
    if len(rows) != 30:
        raise AssertionError("Incomplete Task 18 ensemble")
    write_rows(OUT / "wt_inverse_solutions.csv", rows)


def freeze() -> None:
    path = OUT / "frozen_wt_manifest.json"
    if path.exists():
        raise FileExistsError("Frozen WT manifest is immutable")
    contract = read(OUT / "contract.json")
    payload_paths = sorted((OUT / "wt_parameter_payloads").glob("*.json"))
    if len(payload_paths) != 30:
        raise AssertionError("Expected all ten roots and three conditions")
    payloads = [read(item) for item in payload_paths]
    expected = {(root, condition) for root in contract["root_ids"] for condition in CONDITIONS}
    actual = {(item["row"]["root_id"], item["row"]["condition"]) for item in payloads}
    if actual != expected:
        raise AssertionError("Incomplete root by condition ensemble")
    if not all(item["row"]["fixed_wt_feasible"] for item in payloads):
        raise AssertionError("Freeze must preserve every feasibility decision")
    hashes = {
        relative(item): sha256_file(item)
        for item in OUT.rglob("*")
        if item.is_file() and item != path
    }
    for item in (
        Path(__file__),
        REPO / "tests/test_task18_fixed_wt.py",
        ANALYSIS / "method.md",
        ANALYSIS / "parameter_changes.md",
    ):
        hashes[relative(item)] = sha256_file(item)
    manifest = {
        "manifest_id": "TASK18_FIXED_WT_PRE_GENOTYPE_V1",
        "created_utc": now(),
        "branch": BRANCH,
        "intake_commit": INTAKE,
        "inheritance_commit": INHERITANCE,
        "scientific_baseline_commit": SCIENTIFIC_BASELINE,
        "contract_sha256": sha256_file(OUT / "contract.json"),
        "hashes": hashes,
        "all_decisions": [item["row"] for item in payloads],
        "feasible_parameterisations": [
            {
                "candidate_id": item["row"]["candidate_id"],
                "payload_path": relative(path_item),
                "payload_sha256": sha256_file(path_item),
                "whole_cell_parameters_sha256": sha256_object(
                    item["accepted_whole_cell_parameters"]
                ),
                "ae4_parameters_sha256": sha256_object(item["accepted_ae4_parameters"]),
                "core_state_sha256": item["core_state_sha256"],
                "complete_state_bytes_sha256": item["complete_state_bytes_sha256"],
            }
            for path_item, item in zip(payload_paths, payloads)
        ],
        "all_ten_roots_retained": True,
        "both_routing_families_retained": True,
        "state_solves_performed": 0,
        "genotype_outcomes_used": False,
        "genotype_result_files_read": 0,
        "pushed_checkpoint_required_before_genotype": True,
        "canonical_selection_complete": True,
    }
    write_json(path, manifest)
    print("Task 18 WT inverse solution set frozen; push this checkpoint before genotype evaluation.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "solve", "freeze"))
    args = parser.parse_args()
    {"prepare": prepare, "solve": solve, "freeze": freeze}[args.stage]()


if __name__ == "__main__":
    main()
