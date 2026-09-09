"""Certified steady-state identities and published-baseline audits.

This module is intentionally narrower than a full AE4 cell simulator.  It uses
only equations and numbers printed in Vera-Siguenza et al. (2018).  The paper's
unit and closure conflicts prevent an unambiguous executable reconstruction;
the functions below keep exact stoichiometric results separate from diagnostic
calculations made under explicitly stated, direct unit conversions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np


PAPER_DOI = "10.1007/s11538-017-0370-6"
PAPER_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/"
PALK_DOI = "10.1016/j.jtbi.2010.06.027"
PALK_CORRIGENDUM_URL = (
    "https://www.sciencedirect.com/science/article/pii/S0022519312005668"
)


@dataclass(frozen=True)
class BaselineState:
    """Nominal resting values printed in Table 1 of the 2018 paper."""

    na_i: float = 25.0
    k_i: float = 120.0
    cl_i: float = 50.0
    hco3_i: float = 12.1
    ph_i: float = 6.91
    co2_i: float = 6.6
    na_l: float = 118.7
    k_l: float = 5.6
    cl_l: float = 124.3
    hco3_l: float = 1.5e-4
    ph_l: float = 6.81
    na_e: float = 140.2
    k_e: float = 5.3
    cl_e: float = 102.6
    hco3_e: float = 42.9
    ph_e: float = 7.4
    co2_e: float = 1.9
    cell_volume_pl: float = 1.3
    impermeant_i_mM: float = 82.8
    lumen_osmolyte_mM: float = 48.8
    va_mV: float = -50.24
    vb_mV: float = -62.8
    ca_i_uM: float = 0.058


def proton_mM(ph: float) -> float:
    """Convert pH to millimolar proton concentration."""

    return 10.0 ** (3.0 - ph)


def ae2_turnover_factor(state: BaselineState) -> float:
    """The bracketed, dimensionless factor in published Eq. (25)."""

    k_cl = 5.6
    k_b = 1.0e4
    forward = (state.cl_e / (state.cl_e + k_cl)) * (
        state.hco3_i / (state.hco3_i + k_b)
    )
    reverse = (state.cl_i / (state.cl_i + k_cl)) * (
        state.hco3_e / (state.hco3_e + k_b)
    )
    return forward - reverse


def ae4_turnover_per_second(state: BaselineState) -> float:
    """The bracketed rate in printed Eq. (31), using its stated mM units."""

    k_plus = 1.92e-2  # mM^-4 s^-1, Table 7
    k_minus = 1.3e-5  # mM^-4 s^-1, Table 7
    forward = (
        k_plus
        * state.cl_e
        * state.hco3_i**2
        * (state.na_i + state.k_i)
    )
    reverse = (
        k_minus
        * state.cl_i
        * state.hco3_e**2
        * (state.na_e + state.k_e)
    )
    return forward - reverse


def nhe1_flux_fmol_per_s(state: BaselineState) -> float:
    """Published Eq. (26), whose activity is tabulated in fmol/s."""

    g_nhe1 = 0.0305
    k_h = 4.5e-4
    k_na = 15.0
    h_i = proton_mM(state.ph_i)
    h_e = proton_mM(state.ph_e)
    forward = (h_i / (h_i + k_h)) ** 2 * (state.na_e / (state.na_e + k_na))
    reverse = (state.na_i / (state.na_i + k_na)) * (h_e / (h_e + k_h)) ** 2
    return g_nhe1 * (forward - reverse)


def buffer_rate_mM_per_s(state: BaselineState) -> float:
    """Printed mass-action buffer expression in Appendix 8."""

    k1 = 11.0
    k_minus1 = 2.6e4
    return k1 * state.co2_i - k_minus1 * proton_mM(state.ph_i) * state.hco3_i


def recomputed_osmolarities_mM(state: BaselineState) -> dict[str, float]:
    """Recompute the three printed ionic sums from Table 1 values."""

    lumen = (
        state.k_l
        + state.na_l
        + proton_mM(state.ph_l)
        + state.cl_l
        + state.hco3_l
        + state.lumen_osmolyte_mM
    )
    cell = (
        state.k_i
        + state.na_i
        + proton_mM(state.ph_i)
        + state.cl_i
        + state.hco3_i
        + state.co2_i
        + state.impermeant_i_mM
    )
    bath = (
        state.k_e
        + state.na_e
        + proton_mM(state.ph_e)
        + state.cl_e
        + state.hco3_e
        + state.co2_e
    )
    return {"lumen": lumen, "cell": cell, "bath": bath}


def water_flux_audit() -> dict[str, float]:
    """Evaluate Eqs. (10)-(13) at the paper's reported osmolarities.

    A millimolar concentration is 1e-3 mol/L.  Therefore a permeability in
    L^2 mol^-1 s^-1 times an osmotic difference in mM gives L/s after the
    factor 1e-3.  The calculation is a diagnostic of the printed record, not
    a replacement calibration.
    """

    osm_l, osm_i, osm_e = 297.4, 296.6, 292.6
    p_a, p_b, p_t = 4.32e-12, 5.15e-11, 2.6e-13
    q_a = p_a * (osm_l - osm_i) * 1e-3
    q_b = p_b * (osm_i - osm_e) * 1e-3
    q_t = p_t * (osm_l - osm_e) * 1e-3
    q_total = q_a + q_t
    q_total_uL_min = q_total * 1e6 * 60.0
    reported_rest_uL_min = 0.24e-8
    return {
        "q_a_L_per_s": q_a,
        "q_b_L_per_s": q_b,
        "q_t_L_per_s": q_t,
        "q_total_L_per_s": q_total,
        "q_b_over_q_a": q_b / q_a,
        "q_total_uL_per_min": q_total_uL_min,
        "reported_rest_q_uL_per_min": reported_rest_uL_min,
        "q_total_over_reported": q_total_uL_min / reported_rest_uL_min,
    }


def apical_qss_current_audit(state: BaselineState) -> dict[str, float]:
    """Evaluate the printed apical QSS constraint at the Table 1 baseline.

    Conductance in nS times potential in mV is current in pA.  The signs here
    follow Eqs. (23), (36), (37), and the QSS equation printed in section 2.9.
    """

    r = 8.3144621
    temperature = 310.0
    faraday = 96485.3365
    rt_over_f_mV = r * temperature / faraday * 1e3
    p_cacc = (state.ca_i_uM / (state.ca_i_uM + 0.26)) ** 1.46
    v_cacc_mV = -rt_over_f_mV * np.log(state.cl_l / state.cl_i)
    i_cacc_pA = 71.3 * p_cacc * (state.va_mV - v_cacc_mV)

    v_t_mV = state.va_mV - state.vb_mV
    v_na_t_mV = rt_over_f_mV * np.log(state.na_l / state.na_e)
    v_k_t_mV = rt_over_f_mV * np.log(state.k_l / state.k_e)
    i_na_t_pA = 12.46 * (v_t_mV - v_na_t_mV)
    i_k_t_pA = 0.9 * (v_t_mV - v_k_t_mV)
    residual_pA = -i_cacc_pA - i_na_t_pA - i_k_t_pA
    tight_total = abs(i_na_t_pA) + abs(i_k_t_pA)
    return {
        "p_cacc": float(p_cacc),
        "v_cacc_mV": float(v_cacc_mV),
        "i_cacc_pA": float(i_cacc_pA),
        "v_t_mV": float(v_t_mV),
        "i_na_t_pA": float(i_na_t_pA),
        "i_k_t_pA": float(i_k_t_pA),
        "qss_residual_pA": float(residual_pA),
        "abs_residual_over_abs_tight_current": float(abs(residual_pA) / tight_total),
    }


def kinetic_scale_audit(state: BaselineState) -> dict[str, float | str]:
    """Compare direct printed-unit use with an inconsistent M-state reuse.

    For Ae4, Table 7's density (amol/um^3) is multiplied by the 1.3 pL
    (1300 um^3) cell volume and by Eq. (31)'s printed s^-1 bracket.  For the
    buffer, mM/s times pL equals fmol/s.  These are not asserted to be the
    authors' hidden implementation choices.  Because the rate constants are
    printed in mM^-4 s^-1, the mM calculation is the literal convention.  The
    M calculation shows the 10^12 error caused by reusing the same numerical
    constants with M-valued states; a proper unit conversion would remove it.
    """

    phi2 = ae2_turnover_factor(state)
    phi4 = ae4_turnover_per_second(state)
    phi4_if_concentrations_are_molar = phi4 * 1e-12
    j2 = 0.01807 * phi2
    cell_volume_um3 = state.cell_volume_pl * 1000.0
    ae4_amount_fmol = 0.66 * cell_volume_um3 / 1000.0
    j4 = ae4_amount_fmol * phi4
    j_nhe1 = nhe1_flux_fmol_per_s(state)
    j_buffer = buffer_rate_mM_per_s(state) * state.cell_volume_pl
    return {
        "ae4_unit_interpretation": (
            "The mM calculation is literal for the printed mM^-4 rate constants. "
            "The M-state value shows the error from reusing those same numerical "
            "constants without conversion; it is not a second physical convention."
        ),
        "ae2_turnover_factor": phi2,
        "ae2_flux_fmol_per_s": j2,
        "ae4_printed_turnover_per_s": phi4,
        "ae4_turnover_if_numeric_concentrations_are_molar_per_s": (
            phi4_if_concentrations_are_molar
        ),
        "ae4_mM_to_M_fourth_power_ratio": 1e12,
        "ae4_direct_conversion_flux_fmol_per_s": j4,
        "nhe1_flux_fmol_per_s": j_nhe1,
        "buffer_direct_conversion_flux_fmol_per_s": j_buffer,
        "buffer_over_nhe1": j_buffer / j_nhe1,
        "ae4_over_nhe1": j4 / j_nhe1,
    }


def nkcc_turnover_audit(state: BaselineState) -> dict[str, float | str]:
    """Compare the 2018 NKCC table with the upstream Palk corrigendum.

    The corrigendum gives a2 and a4 in M^-4.  Their numerical values must be
    multiplied by 1e-12 when the state is represented in mM.  This correction
    is source-backed; it does not resolve the separately fitted whole-cell
    density scale in the 2018 model.
    """

    concentration_product_mM4 = state.na_i * state.k_i * state.cl_i**2
    a1 = 157.55
    a3 = 1.0306
    a2_printed_as_mM = 2.0096e7
    a4_printed_as_mM = 1.3852e6
    a2_corrected_mM = a2_printed_as_mM * 1e-12
    a4_corrected_mM = a4_printed_as_mM * 1e-12

    literal_2018 = (
        a1 - a2_printed_as_mM * concentration_product_mM4
    ) / (a3 + a4_printed_as_mM * concentration_product_mM4)
    corrected = (
        a1 - a2_corrected_mM * concentration_product_mM4
    ) / (a3 + a4_corrected_mM * concentration_product_mM4)
    return {
        "palk_source_doi": PALK_DOI,
        "corrigendum_url": PALK_CORRIGENDUM_URL,
        "concentration_product_mM4": concentration_product_mM4,
        "literal_2018_table_turnover_per_s": literal_2018,
        "corrigendum_converted_turnover_per_s": corrected,
        "a2_corrigendum_converted_mM_minus4_s_minus2": a2_corrected_mM,
        "a4_corrigendum_converted_mM_minus4_s_minus1": a4_corrected_mM,
    }


def acid_base_equilibrium_audit(state: BaselineState) -> dict[str, float]:
    """Compare Appendix 8's printed equilibrium with the Table 1 CO2 value."""

    k1 = 11.0
    k_minus1 = 2.6e4
    predicted_co2 = (k_minus1 / k1) * proton_mM(state.ph_i) * state.hco3_i
    return {
        "equilibrium_predicted_co2_i_mM": predicted_co2,
        "reported_co2_i_mM": state.co2_i,
        "reported_over_equilibrium_prediction": state.co2_i / predicted_co2,
    }


def activity_signature_matrix(na_i: float, k_i: float) -> np.ndarray:
    """AE2/AE4 cycle signatures in (Cl, Na, K, HCO3, H, CO2) balances."""

    total_cation = na_i + k_i
    if total_cation <= 0:
        raise ValueError("Na_i + K_i must be positive")
    alpha = na_i / total_cation
    beta = k_i / total_cation
    ae2 = np.array([1.0, 0.0, 0.0, -1.0, 0.0, 0.0])
    ae4 = np.array([1.0, -alpha, -beta, -2.0, 0.0, 0.0])
    return np.column_stack((ae2, ae4))


def flux_aggregate_matrix() -> np.ndarray:
    """Map (J_Ae2, J_Ae4) to A=J2+J4 and B=J2+2J4."""

    return np.array([[1.0, 1.0], [1.0, 2.0]])


def flux_aggregates(j2: float, j4: float) -> tuple[float, float]:
    """Return shared chloride-cycle aggregate A and bicarbonate aggregate B."""

    return j2 + j4, j2 + 2.0 * j4


def recover_exchanger_fluxes(a: float, b: float) -> tuple[float, float]:
    """Exact inverse: J2=2A-B and J4=B-A."""

    return 2.0 * a - b, b - a


def forward_flux_wedge(a: float, b: float, atol: float = 1e-12) -> bool:
    """Whether aggregates are compatible with nonnegative AE2/AE4 cycle rates."""

    return a >= -atol and b >= a - atol and b <= 2.0 * a + atol


def implicit_sensitivity(f_u: np.ndarray, f_theta: np.ndarray) -> np.ndarray:
    """Compute D_theta u* = -F_u^{-1} F_theta."""

    return -np.linalg.solve(np.asarray(f_u, dtype=float), np.asarray(f_theta, dtype=float))


def finite_difference_jacobian(
    function: Callable[[np.ndarray], np.ndarray],
    point: Iterable[float],
    step: float = 1e-6,
) -> np.ndarray:
    """Central finite-difference Jacobian used for independent derivative checks."""

    x = np.asarray(tuple(point), dtype=float)
    value = np.asarray(function(x), dtype=float)
    jacobian = np.empty((value.size, x.size), dtype=float)
    for column in range(x.size):
        direction = np.zeros_like(x)
        direction[column] = step
        jacobian[:, column] = (
            np.asarray(function(x + direction), dtype=float)
            - np.asarray(function(x - direction), dtype=float)
        ) / (2.0 * step)
    return jacobian


def complete_baseline_audit() -> dict[str, object]:
    """Return all deterministic diagnostics in a JSON-serializable dictionary."""

    state = BaselineState()
    signatures = activity_signature_matrix(state.na_i, state.k_i)
    aggregate = flux_aggregate_matrix()
    return {
        "source": {"doi": PAPER_DOI, "url": PAPER_URL},
        "scope": "printed-equation audit; not a reproduced full model",
        "steady_volume_sign_result": (
            "Both printed signs reduce to q_a = q_b at steady state; the sign conflict "
            "does not affect the steady-state identities."
        ),
        "recomputed_osmolarities_mM": recomputed_osmolarities_mM(state),
        "reported_osmolarity_water_flux_audit": water_flux_audit(),
        "apical_qss_current_audit": apical_qss_current_audit(state),
        "kinetic_scale_audit": kinetic_scale_audit(state),
        "nkcc_upstream_corrigendum_audit": nkcc_turnover_audit(state),
        "acid_base_equilibrium_audit": acid_base_equilibrium_audit(state),
        "co2_sign_result": (
            "The individually printed apical and basolateral influxes sum to "
            "P_CO2*(CO2_l + CO2_e - 2*CO2_i), the negative of the aggregate "
            "J_CO2 printed in Appendix 8."
        ),
        "activity_signature_rank": int(np.linalg.matrix_rank(signatures)),
        "activity_signature_cl_hco3_minor_determinant": float(
            np.linalg.det(signatures[[0, 3], :])
        ),
        "flux_aggregate_rank": int(np.linalg.matrix_rank(aggregate)),
        "flux_aggregate_determinant": float(np.linalg.det(aggregate)),
    }
