"""Finite-buffer carbonate chemistry in conserved coordinates.

The dynamic coordinates are total inorganic carbon (TIC) and total
alkalinity (TA), both stored as amounts by :mod:`states`.  Fast hydration,
carbonate dissociation, water ionization, and a finite monoprotic buffer are
reconstructed algebraically.  Consequently internal acid-base reactions
cannot create carbon or alkalinity, while external HCO3, CO2, H, and CO3
fluxes have explicit, charge-correct source signatures.

This is an equilibrium closure, not an infinite-buffer approximation: the
buffer amount is finite and dilution changes its concentration.  If future
WT kinetics require finite carbonic-anhydrase relaxation, the same TIC/TA
coordinates can be retained and one independent species coordinate added.

Provenance: equilibrium identities are ``DERIVED_CONSTRAINT``; the
single-effective-buffer closure is a ``NEW_MODELING_DECISION``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from .parameters import AcidBaseParameters


@dataclass(frozen=True)
class AcidBaseSpeciation:
    ph: float
    h_mM: float
    oh_mM: float
    co2_mM: float
    hco3_mM: float
    co3_mM: float
    buffer_h_mM: float
    buffer_minus_mM: float
    total_carbon_mM: float
    total_alkalinity_mM: float
    alkalinity_residual_mM: float

    @property
    def carbon_residual_mM(self) -> float:
        return self.co2_mM + self.hco3_mM + self.co3_mM - self.total_carbon_mM

    @property
    def buffer_residual_mM(self) -> float:
        return self.buffer_h_mM + self.buffer_minus_mM


def _carbon_fractions(ph: float, parameters: AcidBaseParameters) -> tuple[float, float, float]:
    h_M = 10.0 ** (-ph)
    ka1 = 10.0 ** (-parameters.carbon_pka1)
    ka2 = 10.0 ** (-parameters.carbon_pka2)
    denominator = h_M * h_M + ka1 * h_M + ka1 * ka2
    return (
        h_M * h_M / denominator,
        ka1 * h_M / denominator,
        ka1 * ka2 / denominator,
    )


def total_alkalinity_mM(
    *,
    ph: float,
    total_carbon_mM: float,
    buffer_total_mM: float,
    buffer_pka: float,
    parameters: AcidBaseParameters,
) -> float:
    """Return TA = HCO3 + 2 CO3 + B- + OH - H in mM equivalents."""

    if total_carbon_mM < 0.0 or buffer_total_mM < 0.0:
        raise ValueError("carbon and buffer totals must be nonnegative")
    alpha_co2, alpha_hco3, alpha_co3 = _carbon_fractions(ph, parameters)
    del alpha_co2
    h_mM = 1.0e3 * 10.0 ** (-ph)
    oh_mM = 1.0e3 * 10.0 ** (ph - parameters.water_pkw)
    buffer_minus = buffer_total_mM / (1.0 + 10.0 ** (buffer_pka - ph))
    return (
        total_carbon_mM * alpha_hco3
        + 2.0 * total_carbon_mM * alpha_co3
        + buffer_minus
        + oh_mM
        - h_mM
    )


def speciate(
    *,
    total_carbon_mM: float,
    total_alkalinity_mM_value: float,
    buffer_total_mM: float,
    buffer_pka: float,
    parameters: AcidBaseParameters,
) -> AcidBaseSpeciation:
    """Recover pH and explicit species from TIC, TA, and a finite buffer pool."""

    values = (total_carbon_mM, total_alkalinity_mM_value, buffer_total_mM, buffer_pka)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("acid-base inputs must be finite")
    if total_carbon_mM <= 0.0:
        raise ValueError("total inorganic carbon must be positive")
    if buffer_total_mM < 0.0:
        raise ValueError("buffer total must be nonnegative")

    def residual(ph: float) -> float:
        return total_alkalinity_mM(
            ph=ph,
            total_carbon_mM=total_carbon_mM,
            buffer_total_mM=buffer_total_mM,
            buffer_pka=buffer_pka,
            parameters=parameters,
        ) - total_alkalinity_mM_value

    lower = parameters.ph_lower
    upper = parameters.ph_upper
    r_lower = residual(lower)
    r_upper = residual(upper)
    if r_lower == 0.0:
        ph = lower
    elif r_upper == 0.0:
        ph = upper
    elif r_lower * r_upper > 0.0:
        feasible_lower = total_alkalinity_mM(
            ph=lower,
            total_carbon_mM=total_carbon_mM,
            buffer_total_mM=buffer_total_mM,
            buffer_pka=buffer_pka,
            parameters=parameters,
        )
        feasible_upper = total_alkalinity_mM(
            ph=upper,
            total_carbon_mM=total_carbon_mM,
            buffer_total_mM=buffer_total_mM,
            buffer_pka=buffer_pka,
            parameters=parameters,
        )
        raise ValueError(
            "alkalinity lies outside the speciation bracket: "
            f"target={total_alkalinity_mM_value}, feasible=[{feasible_lower}, {feasible_upper}]"
        )
    else:
        ph = float(brentq(residual, lower, upper, xtol=1.0e-12, rtol=1.0e-13))

    alpha_co2, alpha_hco3, alpha_co3 = _carbon_fractions(ph, parameters)
    co2 = total_carbon_mM * alpha_co2
    hco3 = total_carbon_mM * alpha_hco3
    co3 = total_carbon_mM * alpha_co3
    h_mM = 1.0e3 * 10.0 ** (-ph)
    oh_mM = 1.0e3 * 10.0 ** (ph - parameters.water_pkw)
    buffer_minus = buffer_total_mM / (1.0 + 10.0 ** (buffer_pka - ph))
    buffer_h = buffer_total_mM - buffer_minus
    reconstructed_ta = hco3 + 2.0 * co3 + buffer_minus + oh_mM - h_mM
    return AcidBaseSpeciation(
        ph=ph,
        h_mM=h_mM,
        oh_mM=oh_mM,
        co2_mM=co2,
        hco3_mM=hco3,
        co3_mM=co3,
        buffer_h_mM=buffer_h,
        buffer_minus_mM=buffer_minus,
        total_carbon_mM=total_carbon_mM,
        total_alkalinity_mM=total_alkalinity_mM_value,
        alkalinity_residual_mM=reconstructed_ta - total_alkalinity_mM_value,
    )


@dataclass(frozen=True)
class CarbonAlkalinitySource:
    """External source in the conserved TIC/TA coordinates (fmol/s)."""

    total_carbon_fmol_s: float
    total_alkalinity_fmol_s: float

    @property
    def charge_equivalent_source_fmol_s(self) -> float:
        """Charge source represented by acid/base species: ``-dTA/dt``."""

        return -self.total_alkalinity_fmol_s


def carbon_alkalinity_source(
    *,
    hco3_fmol_s: float = 0.0,
    co3_fmol_s: float = 0.0,
    co2_fmol_s: float = 0.0,
    h_fmol_s: float = 0.0,
    oh_fmol_s: float = 0.0,
) -> CarbonAlkalinitySource:
    """Map explicit species fluxes to conserved TIC and TA sources.

    Positive values add the named species to the compartment.  In particular,
    removing one proton (negative ``h_fmol_s``) raises alkalinity by one
    equivalent, as required for NHE1-mediated proton extrusion.
    """

    values = (hco3_fmol_s, co3_fmol_s, co2_fmol_s, h_fmol_s, oh_fmol_s)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("species fluxes must be finite")
    return CarbonAlkalinitySource(
        total_carbon_fmol_s=hco3_fmol_s + co3_fmol_s + co2_fmol_s,
        total_alkalinity_fmol_s=(
            hco3_fmol_s + 2.0 * co3_fmol_s - h_fmol_s + oh_fmol_s
        ),
    )


def closed_reaction_source_matrix() -> NDArray[np.float64]:
    """Stoichiometric matrix for hydration and buffer dissociation.

    Rows are ``(CO2, HCO3, H, BH, B-)``.  Columns are
    ``CO2 -> HCO3 + H`` and ``BH -> B- + H``.
    """

    return np.asarray(
        (
            (-1.0, 0.0),
            (+1.0, 0.0),
            (+1.0, +1.0),
            (0.0, -1.0),
            (0.0, +1.0),
        ),
        dtype=float,
    )


def closed_reaction_invariant_residuals() -> Mapping[str, NDArray[np.float64]]:
    """Return exact left-null residuals for carbon, buffer, and charge."""

    source = closed_reaction_source_matrix()
    carbon = np.asarray((1.0, 1.0, 0.0, 0.0, 0.0)) @ source
    buffer_sites = np.asarray((0.0, 0.0, 0.0, 1.0, 1.0)) @ source
    charge = np.asarray((0.0, -1.0, +1.0, 0.0, -1.0)) @ source
    return {"carbon": carbon, "buffer_sites": buffer_sites, "charge": charge}


__all__ = (
    "AcidBaseSpeciation",
    "CarbonAlkalinitySource",
    "carbon_alkalinity_source",
    "closed_reaction_invariant_residuals",
    "closed_reaction_source_matrix",
    "speciate",
    "total_alkalinity_mM",
)
