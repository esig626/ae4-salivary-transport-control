"""Osmotic water balance and finite-lumen outflow closure.

Water-flux signs are explicit:

* ``q_b``: bath -> cell;
* ``q_a``: cell -> lumen;
* ``q_t``: bath -> lumen (paracellular);
* ``q_out``: lumen -> ductal sink.

Hence ``dV_i/dt = q_b - q_a`` and
``dV_l/dt = q_a + q_t - q_out``.  The local outflow law is a transparent
compliance closure and must not be interpreted as a calibrated secretion law.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from .parameters import FullModelParameters


@dataclass(frozen=True)
class OsmoticEnvironment:
    """Osmolarities are numerical values in mOsm/L despite compact field names."""

    osm_cell_mOsm: float
    osm_lumen_mOsm: float
    osm_bath_mOsm: float
    lumen_volume_pL: float


@dataclass(frozen=True)
class WaterFluxes:
    bath_to_cell_pL_s: float
    cell_to_lumen_pL_s: float
    bath_to_lumen_pL_s: float
    lumen_outflow_pL_s: float
    cell_volume_source_pL_s: float
    lumen_volume_source_pL_s: float
    volume_conservation_residual_pL_s: float


def cell_impermeant_osmoles_fmol(
    *, other_impermeant_fmol: float, finite_buffer_fmol: float
) -> float:
    """Return all impermeant intracellular osmotic particles.

    ``other_impermeant_fmol`` explicitly excludes the finite acid-base buffer.
    Each buffer molecule remains one osmotic particle whether protonated or
    deprotonated, so the fixed buffer amount is added exactly once.
    """

    if not all(
        math.isfinite(value) and value >= 0.0
        for value in (other_impermeant_fmol, finite_buffer_fmol)
    ):
        raise ValueError("impermeant and finite-buffer osmoles must be nonnegative")
    return other_impermeant_fmol + finite_buffer_fmol


def compartment_osmolarity_mOsm(
    *,
    na_mM: float,
    k_mM: float,
    cl_mM: float,
    total_carbon_mM: float,
    impermeant_osmoles_fmol: float = 0.0,
    volume_pL: float | None = None,
    untracked_osmolyte_mM: float = 0.0,
) -> float:
    """Ideal dilute osmolarity in mOsm/L of explicitly represented solutes.

    TIC contributes one carbon-bearing particle per carbon.  Free H/OH make a
    negligible but recoverable correction and can be added later if a gate is
    sensitive to it.  Impermeant osmoles require a volume for conversion.
    """

    concentrations = (na_mM, k_mM, cl_mM, total_carbon_mM, untracked_osmolyte_mM)
    if not all(math.isfinite(value) and value >= 0.0 for value in concentrations):
        raise ValueError("osmotic concentrations must be finite and nonnegative")
    if impermeant_osmoles_fmol < 0.0:
        raise ValueError("impermeant osmoles must be nonnegative")
    if impermeant_osmoles_fmol > 0.0:
        if volume_pL is None or not math.isfinite(volume_pL) or volume_pL <= 0.0:
            raise ValueError("a positive volume is required for impermeant osmoles")
        impermeant_mM = impermeant_osmoles_fmol / volume_pL
    else:
        impermeant_mM = 0.0
    return sum(concentrations) + impermeant_mM


def evaluate_water_fluxes(
    environment: OsmoticEnvironment, parameters: FullModelParameters
) -> WaterFluxes:
    """Evaluate osmotic inflows and lumen outflow in physical seconds."""

    e = environment
    if not all(
        math.isfinite(value)
        for value in (e.osm_cell_mOsm, e.osm_lumen_mOsm, e.osm_bath_mOsm)
    ):
        raise ValueError("osmolarities must be finite")
    if e.lumen_volume_pL <= 0.0:
        raise ValueError("lumen volume must be positive")
    p = parameters.water
    q_b = p.basolateral_hydraulic_pL_s_mOsm * (
        e.osm_cell_mOsm - e.osm_bath_mOsm
    )
    q_a = p.apical_hydraulic_pL_s_mOsm * (
        e.osm_lumen_mOsm - e.osm_cell_mOsm
    )
    q_t = p.paracellular_hydraulic_pL_s_mOsm * (
        e.osm_lumen_mOsm - e.osm_bath_mOsm
    )
    q_out = p.outflow_rate_s * max(
        e.lumen_volume_pL - p.lumen_dead_volume_pL, 0.0
    )
    dvi = q_b - q_a
    dvl = q_a + q_t - q_out
    # Internal q_a cancels; bath inflows and ductal outflow are external.
    residual = dvi + dvl - (q_b + q_t - q_out)
    return WaterFluxes(
        bath_to_cell_pL_s=q_b,
        cell_to_lumen_pL_s=q_a,
        bath_to_lumen_pL_s=q_t,
        lumen_outflow_pL_s=q_out,
        cell_volume_source_pL_s=dvi,
        lumen_volume_source_pL_s=dvl,
        volume_conservation_residual_pL_s=residual,
    )


__all__ = (
    "OsmoticEnvironment",
    "WaterFluxes",
    "cell_impermeant_osmoles_fmol",
    "compartment_osmolarity_mOsm",
    "evaluate_water_fluxes",
)
