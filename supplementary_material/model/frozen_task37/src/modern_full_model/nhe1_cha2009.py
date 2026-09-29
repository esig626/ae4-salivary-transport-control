"""Cha et al. (2009) Eqs. 1, 3, 5 and 6: eight-state NHE with Mod2.

Published final Table S1 constants: eight-state exchange, Mod2, n_H=3,
m_H=1. Source PDF provenance is recorded in Task 31 source_table_s1.json.
Printed constants are retained exactly; their rounding leaves a small
microscopic-reversibility discrepancy, reported separately from Eq. 3.
DOI: https://doi.org/10.1016/j.bpj.2009.08.053
"""
from __future__ import annotations

from dataclasses import dataclass
import math


NHE1_CHA_2009 = "cha_2009_eight_state_mod2"
AVOGADRO_PER_MOL = 6.02214076e23
INTRACELLULAR_MODIFIER_HILL = 3
EXTRACELLULAR_MODIFIER_HILL = 1


@dataclass(frozen=True)
class ChaNHE1Kinetics:
    """Independent kinetic constants; concentrations in mM, rates in ms^-1.

    The fourth rate can be supplied as printed in Table S1. When omitted,
    it follows the source's microscopic reversibility constraint (Eq. 6),
    which is useful for exact algebraic fixtures. Production uses 183.
    No cardiac density, cell geometry or background transport enters here.
    """

    k1_plus_per_ms: float
    k1_minus_per_ms: float
    k2_plus_per_ms: float
    k_na_o_mM: float
    k_na_i_mM: float
    k_h_o_mM: float
    k_h_i_mM: float
    modifier_k_i_mM: float
    modifier_k_o_mM: float
    k2_minus_table_per_ms: float | None = None

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            if value is None and name == "k2_minus_table_per_ms":
                continue
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(self.k2_minus_per_ms) or self.k2_minus_per_ms <= 0.0:
            raise ValueError("Eq. 6 must give a representable positive k2_minus")

    @property
    def k2_minus_per_ms(self) -> float:
        if self.k2_minus_table_per_ms is not None:
            return self.k2_minus_table_per_ms
        return self.k2_minus_eq6_per_ms

    @property
    def k2_minus_eq6_per_ms(self) -> float:
        # K_H^i K_Na^o / (K_H^o K_Na^i) = k1+ k2+ / (k1- k2-).
        log_rate = (
            math.log(self.k1_plus_per_ms) + math.log(self.k2_plus_per_ms)
            + math.log(self.k_h_o_mM) + math.log(self.k_na_i_mM)
            - math.log(self.k1_minus_per_ms) - math.log(self.k_h_i_mM)
            - math.log(self.k_na_o_mM)
        )
        try:
            return math.exp(log_rate)
        except OverflowError as exc:
            raise ValueError("Eq. 6 rate overflows") from exc

    @property
    def turnover_bound_per_s(self) -> float:
        """Conservative bound from Eq. 3 and 0 <= Mod2 <= 1."""
        return 1000.0 * max(
            min(self.k1_plus_per_ms, self.k2_plus_per_ms),
            min(self.k1_minus_per_ms, self.k2_minus_per_ms),
        )


def published_cha_kinetics() -> ChaNHE1Kinetics:
    """Exact printed final parameters, Document S1, Table S1, PDF page 5.

    The source N=489900 is intentionally excluded: salivary carrier amount
    is the one independent WT-only calibration parameter.
    """
    return ChaNHE1Kinetics(
        k1_plus_per_ms=10.5, k1_minus_per_ms=0.201,
        k2_plus_per_ms=15.8, k2_minus_table_per_ms=183.0,
        k_na_o_mM=195.0, k_na_i_mM=16.2,
        k_h_o_mM=1.62e-3, k_h_i_mM=6.05e-4,
        modifier_k_i_mM=3.07e-5, modifier_k_o_mM=4.8e-7,
    )


def _softplus(value: float) -> float:
    return max(value, 0.0) + math.log1p(math.exp(-abs(value)))


def _log_binding(concentration: float, dissociation: float) -> tuple[float, float]:
    """Log occupied/unoccupied fractions without concentration overflow."""
    z = math.log(concentration) - math.log(dissociation)
    return -_softplus(-z), -_softplus(z)


@dataclass(frozen=True)
class ChaNHE1Rate:
    exchange_per_s: float
    modifier: float
    regulated_per_s: float
    forward_affinity: float


def cha_nhe1_rate(
    *, na_i_mM: float, h_i_mM: float, na_o_mM: float, h_o_mM: float,
    kinetics: ChaNHE1Kinetics,
) -> ChaNHE1Rate:
    """Return per-carrier turnover; positive means Na inward/H outward.

    Eq. 3: (a*b-c*d)/(a+b+c+d), where a=k1+*P_o(Na only),
    b=k2+*P_i(H only), c=k1-*P_i(Na only), d=k2-*P_o(H only).
    A log-domain product difference avoids overflow/underflow. For the
    printed parameter set, retain its small Eq. 6 rounding discrepancy;
    do not silently replace the reported fourth transition rate.
    """
    for value in (na_i_mM, h_i_mM, na_o_mM, h_o_mM):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("NHE Na and H concentrations must be finite and positive")
    k = kinetics
    no, no_free = _log_binding(na_o_mM, k.k_na_o_mM)
    ni, ni_free = _log_binding(na_i_mM, k.k_na_i_mM)
    ho, ho_free = _log_binding(h_o_mM, k.k_h_o_mM)
    hi, hi_free = _log_binding(h_i_mM, k.k_h_i_mM)
    log_a = math.log(k.k1_plus_per_ms) + no + ho_free
    log_b = math.log(k.k2_plus_per_ms) + hi + ni_free
    log_c = math.log(k.k1_minus_per_ms) + ni + hi_free
    log_d = math.log(k.k2_minus_per_ms) + ho + no_free
    rates = (log_a, log_b, log_c, log_d)
    maximum = max(rates)
    log_sum = maximum + math.log(math.fsum(math.exp(v - maximum) for v in rates))
    affinity = math.fsum((
        math.log(na_o_mM), math.log(h_i_mM),
        -math.log(na_i_mM), -math.log(h_o_mM),
    ))
    cycle_affinity = affinity
    if k.k2_minus_table_per_ms is not None:
        cycle_affinity += math.log(k.k2_minus_eq6_per_ms / k.k2_minus_per_ms)
    if cycle_affinity >= 0.0:
        exchange = math.exp(log_a + log_b - log_sum) * -math.expm1(-cycle_affinity)
    else:
        exchange = -math.exp(log_c + log_d - log_sum) * -math.expm1(cycle_affinity)
    log_inhibition = (
        _softplus(EXTRACELLULAR_MODIFIER_HILL * (
            math.log(h_o_mM) - math.log(k.modifier_k_o_mM)
        ))
        + INTRACELLULAR_MODIFIER_HILL * (
            math.log(k.modifier_k_i_mM) - math.log(h_i_mM)
        )
    )
    modifier = math.exp(-_softplus(log_inhibition))
    exchange_per_s = 1000.0 * exchange
    return ChaNHE1Rate(exchange_per_s, modifier, exchange_per_s * modifier, affinity)


def cha_nhe1_flux_fmol_s(
    *, na_i_mM: float, h_i_mM: float, na_o_mM: float, h_o_mM: float,
    kinetics: ChaNHE1Kinetics, carrier_amount_fmol: float,
    activity_scale: float = 1.0,
) -> float:
    """Eq. 1 in amount coordinates: amount * Mod2 * turnover (fmol/s).

    Carrier amount is N/A * 1e15. No division by cell volume is needed:
    the salivary model evolves amounts. One fmol is 6.02214076e8 carriers.
    """
    for value in (carrier_amount_fmol, activity_scale):
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("NHE amount and activity must be finite and nonnegative")
    rate = cha_nhe1_rate(
        na_i_mM=na_i_mM, h_i_mM=h_i_mM, na_o_mM=na_o_mM,
        h_o_mM=h_o_mM, kinetics=kinetics,
    )
    return carrier_amount_fmol * activity_scale * rate.regulated_per_s
