"""Exact structural checks for the Task 12 AE4 candidate family.

The functions in this module do not simulate the salivary chassis.  They
encode consequences of stoichiometry, thermodynamic sign, and the historical
steady water/lumen closure that can be checked without fitting a kinetic law.

Positive AE4 cycle flux imports ``chloride`` chloride ions and exports
``cation`` monovalent cations and ``bicarbonate`` bicarbonate ions.  Positive
AE2 flux imports one chloride and exports one bicarbonate.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, log, sqrt


@dataclass(frozen=True)
class AE4Stoichiometry:
    """Integer AE4 cycle counts in the positive, chloride-inward direction."""

    chloride: int
    cation: int
    bicarbonate: int

    def __post_init__(self) -> None:
        if self.chloride <= 0 or self.bicarbonate <= 0 or self.cation < 0:
            raise ValueError(
                "anion counts must be positive and cation count nonnegative"
            )

    @property
    def charge_residual(self) -> int:
        """Net charge moved into the cell per cycle, in elementary charges."""

        return self.bicarbonate - self.chloride - self.cation

    @property
    def is_electroneutral(self) -> bool:
        return self.charge_residual == 0

    @property
    def bicarbonate_per_chloride(self) -> float:
        return self.bicarbonate / self.chloride

    @property
    def cation_per_chloride(self) -> float:
        return self.cation / self.chloride


SOURCE_STOICHIOMETRIES = {
    "1Cl_1Cat_2B": AE4Stoichiometry(1, 1, 2),
    "1Cl_2Cat_3B": AE4Stoichiometry(1, 2, 3),
    "2Cl_1Cat_3B": AE4Stoichiometry(2, 1, 3),
}


def balance_signature(
    stoichiometry: AE4Stoichiometry,
    *,
    sodium_fraction: float,
) -> tuple[float, float, float, float]:
    """AE4 source vector in ``(Cl, Na, K, HCO3)`` amount balances."""

    if not 0.0 <= sodium_fraction <= 1.0:
        raise ValueError("sodium_fraction must lie in [0, 1]")
    sodium = stoichiometry.cation * sodium_fraction
    potassium = stoichiometry.cation - sodium
    return (
        float(stoichiometry.chloride),
        -sodium,
        -potassium,
        -float(stoichiometry.bicarbonate),
    )


def fixed_state_candidate_residual(
    *,
    historical_flux: float,
    candidate_flux: float,
    historical_signature: tuple[float, float, float, float],
    candidate_signature: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    """Residual created by replacing AE4 at a historically closed state.

    If the fixed non-AE4 chassis closes with ``historical_signature`` times
    ``historical_flux``, its non-AE4 residual is the negative of that vector.
    A replacement law can change only the scalar cycle flux at the same state.
    It closes that state exactly only when its signature is collinear with the
    historical one and the scalar magnitudes agree.
    """

    return tuple(
        candidate_flux * candidate - historical_flux * historical
        for historical, candidate in zip(
            historical_signature, candidate_signature, strict=True
        )
    )


def baseline_net_k_fraction_bound(
    *, required_cycle_flux: float, amount_tolerance: float
) -> float:
    """Bound a hidden baseline K source under fixed-row closure.

    If the frozen chloride row requires ``required_cycle_flux`` events and the
    historical potassium row requires zero AE4 source, an absolute K-balance
    tolerance ``amount_tolerance`` admits at most their ratio as the net K
    source fraction.  This constrains net transport at that row, not latent K
    capacity at a state where a K branch happens to be at equilibrium.
    """

    if required_cycle_flux == 0.0:
        raise ValueError("required_cycle_flux must be nonzero")
    if amount_tolerance < 0.0:
        raise ValueError("amount_tolerance must be nonnegative")
    return amount_tolerance / abs(required_cycle_flux)


def mixed_cation_fixed_row_correction(
    *, required_cycle_flux: float, potassium_fraction: float
) -> tuple[float, float, float, float]:
    """Non-AE4 correction for replacing a Na-only event by Na/K partition.

    The result uses ``(Na, K, Cl, HCO3)`` order.  A positive K fraction in the
    replacement AE4 event removes too little Na and removes K that the frozen
    historical row did not require.  The correction restores those balances.
    """

    if not 0.0 <= potassium_fraction <= 1.0:
        raise ValueError("potassium_fraction must lie in [0, 1]")
    shifted = required_cycle_flux * potassium_fraction
    return (-shifted, shifted, 0.0, 0.0)


def nak_nhe1_repair_coefficients(
    correction: tuple[float, float, float, float],
) -> tuple[float, float]:
    """Resolve a cation-only correction into Na/K-pump and NHE1 fluxes.

    Pump and NHE1 signatures are ``(-3,+2,0,0)`` and ``(+1,0,0,0)``.
    The returned pair is ``(pump_flux, nhe1_flux)``.  Corrections with chloride
    or bicarbonate components are outside this two-path span.
    """

    na, potassium, chloride, bicarbonate = correction
    if chloride != 0.0 or bicarbonate != 0.0:
        raise ValueError("correction is outside the pump/NHE1 cation subspace")
    pump = potassium / 2.0
    nhe1 = na + 3.0 * pump
    return pump, nhe1


def required_wt_observable_factor(
    *, current_knockout_to_wt_ratio: float, target_knockout_to_wt_ratio: float
) -> float:
    """WT gain required to hit a lower KO/WT ratio while KO stays fixed."""

    if current_knockout_to_wt_ratio <= 0.0:
        raise ValueError("current knockout/WT ratio must be positive")
    if target_knockout_to_wt_ratio <= 0.0:
        raise ValueError("target knockout/WT ratio must be positive")
    return current_knockout_to_wt_ratio / target_knockout_to_wt_ratio


def required_mean_gate_elasticity(
    *,
    activation_fold: float,
    current_knockout_to_wt_ratio: float,
    target_knockout_to_wt_ratio: float,
) -> float:
    """Mean log-elasticity a WT stimulus gate must achieve.

    Knockout output is assumed independent of the WT-only gate and all outputs
    are positive.  The result is a necessary average elasticity over a gate
    change from one to ``activation_fold``; it does not assert monotonicity of
    the nonlinear chassis response.
    """

    if activation_fold <= 1.0:
        raise ValueError("activation_fold must exceed one")
    required_factor = required_wt_observable_factor(
        current_knockout_to_wt_ratio=current_knockout_to_wt_ratio,
        target_knockout_to_wt_ratio=target_knockout_to_wt_ratio,
    )
    return log(required_factor) / log(activation_fold)


def c1_stimulus_source_comparison(
    *, basal_cycle_flux: float, activation_fold: float
) -> tuple[
    tuple[float, float, float, float],
    tuple[float, float, float, float],
]:
    """Compare common C1 gating with stimulus-recruited K transport.

    Returns ``(common_gate, k_recruitment)`` intracellular source vectors in
    ``(Na, K, Cl, HCO3)`` order.  Both diagnostics have the same total anion
    source at the same state.  K recruitment instead substitutes K export for
    the recruited fraction of Na export, so it rotates the cation signature.
    """

    if activation_fold < 1.0:
        raise ValueError("activation_fold must be at least one")
    total = activation_fold * basal_cycle_flux
    recruited = (activation_fold - 1.0) * basal_cycle_flux
    common = (-total, 0.0, total, -2.0 * total)
    k_recruitment = (
        -basal_cycle_flux,
        -recruited,
        total,
        -2.0 * total,
    )
    return common, k_recruitment


def parallel_112_supported_chloride_flux(
    *, nkcc_cycle_flux: float, pump_cycle_flux: float, k_ae4_cycle_flux: float
) -> float:
    """Steady chloride efflux implied by the separate K balance.

    For parallel Na/K 1:1:2 AE4 branches, exact steady cation and anion
    balances reduce to ``C = N + 3 P - J_K``.  The negative K-branch term is a
    balance consequence, not an assumed channel interaction.
    """

    return nkcc_cycle_flux + 3.0 * pump_cycle_flux - k_ae4_cycle_flux


def recover_ae2_ae4_fluxes(
    chloride_aggregate: float,
    bicarbonate_aggregate: float,
    stoichiometry: AE4Stoichiometry,
) -> tuple[float, float]:
    """Recover signed AE2 and AE4 cycle fluxes from two balance aggregates.

    ``chloride_aggregate = C_apical - 2 J_NKCC`` and
    ``bicarbonate_aggregate = J_buffer``.  At steady state they equal

    ``J_AE2 + chloride * J_AE4`` and
    ``J_AE2 + bicarbonate * J_AE4``, respectively.
    """

    denominator = stoichiometry.bicarbonate - stoichiometry.chloride
    if denominator == 0:
        raise ValueError("AE2 and AE4 have identical anion signatures")
    j4 = (bicarbonate_aggregate - chloride_aggregate) / denominator
    j2 = (
        stoichiometry.bicarbonate * chloride_aggregate
        - stoichiometry.chloride * bicarbonate_aggregate
    ) / denominator
    return j2, j4


def forward_flux_wedge(
    chloride_aggregate: float,
    bicarbonate_aggregate: float,
    stoichiometry: AE4Stoichiometry,
    *,
    tolerance: float = 1e-12,
) -> bool:
    """Return whether both reconstructed exchanger cycle rates are nonnegative."""

    j2, j4 = recover_ae2_ae4_fluxes(
        chloride_aggregate, bicarbonate_aggregate, stoichiometry
    )
    return j2 >= -tolerance and j4 >= -tolerance


def delta_g_over_rt(
    stoichiometry: AE4Stoichiometry,
    *,
    chloride_i: float,
    chloride_o: float,
    cation_i: float,
    cation_o: float,
    bicarbonate_i: float,
    bicarbonate_o: float,
) -> float:
    """Dimensionless affinity for chloride influx/cation-base efflux.

    The expression is the ideal-solution reaction Gibbs energy used by
    Pena-Munzenmayer et al. (2016).  Negative values favor chloride influx.
    Electrical-potential terms cancel exactly for an electroneutral cycle.
    """

    values = (
        chloride_i,
        chloride_o,
        cation_i,
        cation_o,
        bicarbonate_i,
        bicarbonate_o,
    )
    if any(value <= 0.0 for value in values):
        raise ValueError("thermodynamic concentrations must be positive")
    if not stoichiometry.is_electroneutral:
        raise ValueError("voltage-free affinity requires electroneutral stoichiometry")
    return (
        stoichiometry.chloride * log(chloride_i / chloride_o)
        + stoichiometry.cation * log(cation_o / cation_i)
        + stoichiometry.bicarbonate * log(bicarbonate_o / bicarbonate_i)
    )


def equilibrium_chloride_i(
    stoichiometry: AE4Stoichiometry,
    *,
    chloride_o: float,
    cation_i: float,
    cation_o: float,
    bicarbonate_i: float,
    bicarbonate_o: float,
) -> float:
    """Intracellular chloride at zero affinity, holding other species fixed."""

    return chloride_o * (
        (cation_i / cation_o) ** stoichiometry.cation
        * (bicarbonate_i / bicarbonate_o) ** stoichiometry.bicarbonate
    ) ** (1.0 / stoichiometry.chloride)


def equilibrium_cation_i(
    stoichiometry: AE4Stoichiometry,
    *,
    chloride_i: float,
    chloride_o: float,
    cation_o: float,
    bicarbonate_i: float,
    bicarbonate_o: float,
) -> float:
    """Intracellular cation concentration at zero affinity."""

    return cation_o * (
        (chloride_i / chloride_o) ** stoichiometry.chloride
        * (bicarbonate_o / bicarbonate_i) ** stoichiometry.bicarbonate
    ) ** (1.0 / stoichiometry.cation)


def equilibrium_bicarbonate_i(
    stoichiometry: AE4Stoichiometry,
    *,
    chloride_i: float,
    chloride_o: float,
    cation_i: float,
    cation_o: float,
    bicarbonate_o: float,
) -> float:
    """Intracellular bicarbonate concentration at zero affinity."""

    return bicarbonate_o * exp(
        (
            stoichiometry.chloride * log(chloride_i / chloride_o)
            + stoichiometry.cation * log(cation_o / cation_i)
        )
        / stoichiometry.bicarbonate
    )


def reversible_rate(positive_mobility: float, delta_g_rt: float) -> float:
    """A minimal detailed-balance rate with sign controlled by affinity."""

    if positive_mobility < 0.0:
        raise ValueError("mobility must be nonnegative")
    return positive_mobility * (1.0 - exp(delta_g_rt))


def pooled_partition_error(
    *,
    forward_factor: float,
    reverse_factor: float,
    na_i: float,
    k_i: float,
    na_o: float,
    k_o: float,
) -> float:
    """Na-branch flux minus intracellular-fraction allocation of pooled flux.

    This applies to a one-cation cycle with identical Na and K branch kinetic
    factors.  A nonzero result proves that allocating the *net* pooled cycle by
    intracellular mole fractions is not the species-resolved reversible law.
    """

    cation_i = na_i + k_i
    cation_o = na_o + k_o
    if min(cation_i, cation_o) <= 0.0:
        raise ValueError("cation sums must be positive")
    na_branch = forward_factor * na_i - reverse_factor * na_o
    pooled = forward_factor * cation_i - reverse_factor * cation_o
    return na_branch - (na_i / cation_i) * pooled


def historical_flow_from_chloride_flux(
    chloride_flux: float,
    *,
    water_gain: float,
    fixed_osmotic_offset: float,
) -> float:
    """Positive steady flow implied by the seven-state water/lumen closure.

    At steady cell volume, ``Q = water_gain * (2 Cl_l + offset)`` and the
    lumen chloride balance is ``chloride_flux = Q * Cl_l``.  Eliminating
    ``Cl_l`` gives the positive quadratic root below.
    """

    if water_gain <= 0.0:
        raise ValueError("water_gain must be positive")
    discriminant = (
        (water_gain * fixed_osmotic_offset) ** 2
        + 8.0 * water_gain * chloride_flux
    )
    if discriminant < 0.0:
        raise ValueError("no real steady flow for these inputs")
    return (
        water_gain * fixed_osmotic_offset + sqrt(discriminant)
    ) / 2.0


def historical_chloride_flux_from_flow(
    flow: float,
    *,
    water_gain: float,
    fixed_osmotic_offset: float,
) -> float:
    """Inverse of :func:`historical_flow_from_chloride_flux`."""

    if water_gain <= 0.0:
        raise ValueError("water_gain must be positive")
    return flow * (flow / water_gain - fixed_osmotic_offset) / 2.0


def historical_flow_slope(
    flow: float,
    *,
    water_gain: float,
    fixed_osmotic_offset: float,
) -> float:
    """Exact ``dQ/dC`` along the positive historical steady branch."""

    denominator = 2.0 * flow - water_gain * fixed_osmotic_offset
    if denominator <= 0.0:
        raise ValueError("inputs are not on the monotone positive branch")
    return 2.0 * water_gain / denominator
