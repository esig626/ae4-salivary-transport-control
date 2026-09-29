"""Task52's two prescribed charge/osmotic onset projections.

Only closed-form state construction is performed. The frozen Task37 package
must first be loaded under ``task52_pre_palk`` by the task loader. Importing
this module alone imports no scientific model and performs no evaluation.
"""
from __future__ import annotations

import math


def project_onset(parameters, original_state, cl_mM, ph, mode="central"):
    """Return a measured-Cl/pH onset state and its analytical accounting.

    Central retains Task37 Na concentration and volume, solving K and TIC.
    Alternate retains Task37 Na and TIC concentrations, solving K and volume.
    Both leave lumen, regulatory coordinates and parameter objects unchanged.
    Physical acceptance and independent pH/speciation recovery belong to the
    calling Task52 verification process; no transporter or RHS is called here.
    """
    if mode not in ("central", "alternate"):
        raise ValueError("projection mode must be 'central' or 'alternate'")
    original = [float(value) for value in original_state]
    if len(original) != 13:
        raise ValueError("Task37 onset state must have exactly 13 coordinates")
    if not all(math.isfinite(value) for value in original):
        raise ValueError("original onset state must be finite")
    if any(value <= 0.0 for value in original[:12]):
        raise ValueError("original conserved core must be positive")
    chloride = float(cl_mM)
    imposed_ph = float(ph)
    if not math.isfinite(chloride) or chloride <= 0.0:
        raise ValueError("imposed chloride must be finite and positive")
    if not math.isfinite(imposed_ph):
        raise ValueError("imposed pH must be finite")
    acid = parameters.acid_base
    if not acid.ph_lower <= imposed_ph <= acid.ph_upper:
        raise ValueError("imposed pH lies outside the inherited speciation domain")

    # These imports resolve exclusively to the byte-exact Task37 snapshot.
    from task52_pre_palk.acid_base import _carbon_fractions, total_alkalinity_mM
    from task52_pre_palk.water import (
        cell_impermeant_osmoles_fmol,
        compartment_osmolarity_mOsm,
    )

    geometry = parameters.geometry
    bath = parameters.bath
    original_volume = original[5]
    sodium = original[0] / original_volume
    original_tic = original[3] / original_volume
    buffer_amount = float(geometry.cell_buffer_total_fmol)
    other_osmoles = float(geometry.cell_impermeant_osmoles_fmol)
    fixed_charge = float(geometry.fixed_cell_anion_equivalents_fmol)
    if not math.isfinite(fixed_charge):
        raise ValueError("fixed cell charge must be finite")
    impermeant_total = cell_impermeant_osmoles_fmol(
        other_impermeant_fmol=other_osmoles,
        finite_buffer_fmol=buffer_amount,
    )
    bath_osm = compartment_osmolarity_mOsm(
        na_mM=bath.na_mM,
        k_mM=bath.k_mM,
        cl_mM=bath.cl_mM,
        total_carbon_mM=bath.tic_mM,
        untracked_osmolyte_mM=bath.untracked_osmolyte_mM,
    )
    alpha_co2, alpha_hco3, alpha_co3 = _carbon_fractions(imposed_ph, acid)
    carbon_alkalinity_per_tic = alpha_hco3 + 2.0 * alpha_co3
    buffer_minus_amount = buffer_amount / (
        1.0 + 10.0 ** (acid.cell_buffer_pka - imposed_ph)
    )
    hydrogen = 1.0e3 * 10.0 ** (-imposed_ph)
    hydroxide = 1.0e3 * 10.0 ** (imposed_ph - acid.water_pkw)
    water_alkalinity = hydroxide - hydrogen
    fixed_sum = fixed_charge + buffer_minus_amount + impermeant_total

    if mode == "central":
        volume = original_volume
        numerator = bath_osm - 2.0 * chloride - water_alkalinity - fixed_sum / volume
        denominator = 1.0 + carbon_alkalinity_per_tic
        solved_coordinate = "tic_i_mM"
        tic = numerator / denominator
    else:
        tic = original_tic
        numerator = fixed_sum
        denominator = (
            bath_osm - 2.0 * chloride - water_alkalinity
            - (1.0 + carbon_alkalinity_per_tic) * tic
        )
        solved_coordinate = "volume_i_pL"
        if not math.isfinite(denominator) or denominator == 0.0:
            raise ValueError("alternate projection has a singular volume equation")
        volume = numerator / denominator

    if not all(math.isfinite(value) and value > 0.0 for value in (tic, volume)):
        raise ValueError(
            f"{mode} projection has no positive finite TIC/volume solution: "
            f"TIC={tic!r}, volume={volume!r}"
        )
    potassium = (
        chloride + carbon_alkalinity_per_tic * tic + water_alkalinity
        + (fixed_charge + buffer_minus_amount) / volume - sodium
    )
    if not math.isfinite(potassium) or potassium <= 0.0:
        raise ValueError(f"{mode} projection has no positive finite K solution")
    alkalinity = total_alkalinity_mM(
        ph=imposed_ph,
        total_carbon_mM=tic,
        buffer_total_mM=buffer_amount / volume,
        buffer_pka=acid.cell_buffer_pka,
        parameters=acid,
    )
    if not math.isfinite(alkalinity) or alkalinity <= 0.0:
        raise ValueError("projected alkalinity violates the positive core-state domain")

    projected = original.copy()
    projected[:6] = [
        original[0] if mode == "central" else sodium * volume,
        potassium * volume,
        chloride * volume,
        tic * volume,
        alkalinity * volume,
        volume,
    ]
    if not all(math.isfinite(value) and value > 0.0 for value in projected[:12]):
        raise ValueError("projected core contains a nonpositive or nonfinite coordinate")
    actual_na, actual_k, actual_cl, actual_tic, actual_ta = (
        projected[index] / volume for index in range(5)
    )
    cell_osm = compartment_osmolarity_mOsm(
        na_mM=actual_na,
        k_mM=actual_k,
        cl_mM=actual_cl,
        total_carbon_mM=actual_tic,
        impermeant_osmoles_fmol=impermeant_total,
        volume_pL=volume,
    )
    charge_residual = math.fsum((
        projected[0], projected[1], -projected[2], -projected[4], -fixed_charge,
    ))
    analytic_ta = (
        carbon_alkalinity_per_tic * actual_tic
        + buffer_minus_amount / volume + water_alkalinity
    )
    return {
        "mode": mode,
        "classification": "MEASUREMENT-CONSTRAINED ONSET PROJECTION; NOT AN EQUILIBRIUM",
        "state": projected,
        "inputs": {
            "imposed_cl_i_mM": chloride,
            "imposed_ph_i": imposed_ph,
            "original_na_i_mM": sodium,
            "original_tic_i_mM": original_tic,
            "original_volume_i_pL": original_volume,
            "bath_osmolarity_mOsm": bath_osm,
            "fixed_cell_anion_equivalents_fmol": fixed_charge,
            "cell_buffer_total_fmol": buffer_amount,
            "other_impermeant_osmoles_fmol": other_osmoles,
        },
        "linear_solution": {
            "solved_coordinate": solved_coordinate,
            "numerator": numerator,
            "denominator": denominator,
            "carbon_alkalinity_per_tic": carbon_alkalinity_per_tic,
            "buffer_minus_amount_fmol": buffer_minus_amount,
            "water_alkalinity_mM": water_alkalinity,
            "total_impermeant_osmoles_fmol": impermeant_total,
        },
        "outputs": {
            "na_i_mM": actual_na,
            "k_i_mM": actual_k,
            "cl_i_mM": actual_cl,
            "tic_i_mM": actual_tic,
            "ta_i_mM": actual_ta,
            "imposed_ph_i": imposed_ph,
            "volume_i_pL": volume,
            "cell_osmolarity_mOsm": cell_osm,
            "co2_at_imposed_ph_mM": alpha_co2 * actual_tic,
            "hco3_at_imposed_ph_mM": alpha_hco3 * actual_tic,
            "co3_at_imposed_ph_mM": alpha_co3 * actual_tic,
            "buffer_minus_at_imposed_ph_mM": buffer_minus_amount / volume,
            "buffer_h_at_imposed_ph_mM": (buffer_amount - buffer_minus_amount) / volume,
        },
        "residuals": {
            "cell_bulk_charge_fmol": charge_residual,
            "cell_minus_bath_osmolarity_mOsm": cell_osm - bath_osm,
            "imposed_chloride_mM": actual_cl - chloride,
            "retained_sodium_mM": actual_na - sodium,
            "retained_volume_pL": volume - original_volume if mode == "central" else None,
            "retained_tic_mM": actual_tic - original_tic if mode == "alternate" else None,
            "alkalinity_formula_mM": actual_ta - analytic_ta,
        },
        "unchanged_lumen_and_regulatory_coordinates": projected[6:] == original[6:],
        "normal_speciation_and_physical_gates": "PENDING_CALLER_VERIFICATION",
    }
