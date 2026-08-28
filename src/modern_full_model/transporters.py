"""Source-disciplined AE4 transport for the Task 13B whole-cell model.

The accepted core is the two-conformation quasi-steady reduction of the
Task-13 ``SR2_SHARED_112`` graph.  A common chloride path and separate Na and
K return paths share one finite carrier pool.  Every path is reversible, and
the rate construction enforces local detailed balance exactly.  Positive
branch current denotes basolateral chloride loading of the cell::

    Cl_o + C_i + 2 HCO3_i <-> Cl_i + C_o + 2 HCO3_o,

where ``C`` is Na or K.  The module reports amount sources in fmol/s so it can
be inserted directly into the conserved amount coordinates in ``states.py``.

The optional fast cooperative gate is the QSS limit of an SR5-like catalytic
occupancy.  It changes both directions of a branch by the same positive
factor and therefore cannot move reversal.  Mapping the 2016 apparent Hill
summary onto that gate is an explicit observation-model assumption, not a
measured microscopic site or transported-ion count.

Dynamic beta/cAMP/PKA states deliberately live outside this file.  Their
output enters as ``regulation_gain`` and, in the minimal coupling, scales
active AE4 capacity without changing ion affinity.  State-specific regulation
is not implemented because the 2021 experiments do not identify a regulated
transport transition.

Provenance tags used below:

* ``PRIMARY_MEASUREMENT``: transported Na and K, electroneutrality, the 2016
  apparent dose summaries, and the qualitative 2025 mutant hierarchy.
* ``DERIVED_CONSTRAINT``: source matrix, charge closure, affinities, reversal,
  stationary carrier occupancy, and local detailed balance.
* ``NEW_MODELING_DECISION``: two-state barrier gauge, QSS use, common-capacity
  regulatory coupling, and any fast empirical cooperative gate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Literal, Mapping


FARADAY_C_PER_MOL = 96485.33212
FMOL_TO_MOL = 1.0e-15

AE4_BRANCH_SOURCE_STOICHIOMETRY: Mapping[str, Mapping[str, float]] = (
    MappingProxyType(
        {
            "na": MappingProxyType(
                {
                    "na_i": -1.0,
                    "k_i": 0.0,
                    "cl_i": +1.0,
                    "hco3_i": -2.0,
                    "tic_i": -2.0,
                    "alkalinity_i": -2.0,
                }
            ),
            "k": MappingProxyType(
                {
                    "na_i": 0.0,
                    "k_i": -1.0,
                    "cl_i": +1.0,
                    "hco3_i": -2.0,
                    "tic_i": -2.0,
                    "alkalinity_i": -2.0,
                }
            ),
        }
    )
)


@dataclass(frozen=True)
class AE4Environment:
    """Aqueous concentrations seen by basolateral AE4 (all in mM)."""

    na_i_mM: float
    k_i_mM: float
    cl_i_mM: float
    hco3_i_mM: float
    na_o_mM: float
    k_o_mM: float
    cl_o_mM: float
    hco3_o_mM: float

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")


@dataclass(frozen=True)
class AE4FastCooperativeGate:
    """Positive SR5-QSS branch gate; it is not a transport stoichiometry.

    ``factor`` is applied symmetrically to the forward and reverse effective
    rates.  The 2016 dose curves constrain an apparent response, whereas this
    class places that response on a kinetic barrier.  Consequently the
    :meth:`from_2016_shape_sensitivity` constructor is suitable for nested
    sensitivity/model-comparison work, not for claiming microscopic parameter
    identification.
    """

    ec50_na_mM: float
    ec50_k_mM: float
    hill_na: float
    hill_k: float
    floor_na: float = 0.0
    floor_k: float = 0.0
    context: str = "user-specified fast catalytic gate"

    def __post_init__(self) -> None:
        for name in ("ec50_na_mM", "ec50_k_mM", "hill_na", "hill_k"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("floor_na", "floor_k"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")

    @classmethod
    def from_2016_shape_sensitivity(cls) -> "AE4FastCooperativeGate":
        """Return a declared empirical SR5 shape sensitivity.

        Peña-Münzenmayer et al. (2016), Fig. 9, reported apparent Na and K
        ``EC50`` values of 49 and 62 mM and Hill summaries of 2.0 and 1.8.
        Those apparent shapes are placed on a barrier only for a nested
        sensitivity test.  The reported ``Rmin`` values are assay-observation
        backgrounds and are deliberately *not* converted into AE4 gate floors.
        This mapping does not assert two transported ions, a measured
        allosteric site, or native-cell parameter identity.
        """

        return cls(
            ec50_na_mM=49.0,
            ec50_k_mM=62.0,
            hill_na=2.0,
            hill_k=1.8,
            floor_na=0.0,
            floor_k=0.0,
            context="2016 Fig. 9 apparent EC50/Hill shape sensitivity on a fast barrier",
        )

    def factor(self, branch: Literal["na", "k"], external_cation_mM: float) -> float:
        if not math.isfinite(external_cation_mM) or external_cation_mM <= 0.0:
            raise ValueError("external cation concentration must be finite and positive")
        if branch == "na":
            ec50, hill, floor = self.ec50_na_mM, self.hill_na, self.floor_na
        elif branch == "k":
            ec50, hill, floor = self.ec50_k_mM, self.hill_k, self.floor_k
        else:
            raise ValueError(f"unknown AE4 branch {branch!r}")
        log_ratio = hill * math.log(ec50 / external_cation_mM)
        if log_ratio > 700.0:
            occupancy = 0.0
        elif log_ratio < -700.0:
            occupancy = 1.0
        else:
            occupancy = 1.0 / (1.0 + math.exp(log_ratio))
        return floor + (1.0 - floor) * occupancy


@dataclass(frozen=True)
class AE4Parameters:
    """Unmeasured kinetic lumps for the QSS shared-carrier module.

    Attempt rates are in s^-1 per carrier and carrier amount is in fmol.  No
    source establishes their absolute values, so the four required values have
    no scientific defaults.  The Na/K attempt-rate fraction is the explicit
    mixed-bath routing degree of freedom that pure-cation assays do not fix.
    """

    carrier_amount_fmol: float
    common_cl_attempt_rate_s: float
    na_loaded_attempt_rate_s: float
    k_loaded_attempt_rate_s: float
    cooperative_gate: AE4FastCooperativeGate | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.carrier_amount_fmol) or self.carrier_amount_fmol < 0.0:
            raise ValueError("carrier_amount_fmol must be finite and nonnegative")
        for name in (
            "common_cl_attempt_rate_s",
            "na_loaded_attempt_rate_s",
            "k_loaded_attempt_rate_s",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")

    @classmethod
    def equal_branch_attempts(
        cls,
        *,
        carrier_amount_fmol: float,
        common_cl_attempt_rate_s: float,
        loaded_attempt_rate_s: float,
        cooperative_gate: AE4FastCooperativeGate | None = None,
    ) -> "AE4Parameters":
        """Named equal-attempt mixed-bath model, not a measured equality."""

        return cls(
            carrier_amount_fmol=carrier_amount_fmol,
            common_cl_attempt_rate_s=common_cl_attempt_rate_s,
            na_loaded_attempt_rate_s=loaded_attempt_rate_s,
            k_loaded_attempt_rate_s=loaded_attempt_rate_s,
            cooperative_gate=cooperative_gate,
        )

    @property
    def mixed_bath_na_attempt_fraction(self) -> float:
        total = self.na_loaded_attempt_rate_s + self.k_loaded_attempt_rate_s
        return self.na_loaded_attempt_rate_s / total


@dataclass(frozen=True)
class AE4BranchModulation:
    """Symmetric branch-conductance changes that preserve reversal.

    Values for mutants are qualitative same-assay activity encodings.  They
    must not be interpreted as a particular binding or flip rate.
    """

    na_scale: float = 1.0
    k_scale: float = 1.0
    label: str = "WT"

    def __post_init__(self) -> None:
        for name in ("na_scale", "k_scale"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")

    @classmethod
    def assay_variant(cls, name: str) -> "AE4BranchModulation":
        """Return Task-13 qualitative encodings of the 2025 hierarchy."""

        variants = {
            "WT": cls(1.0, 1.0, "WT"),
            "T448I": cls(0.5, 0.5, "T448I approximate assay encoding"),
            "T756A": cls(0.7, 0.7, "T756A approximate assay encoding"),
            # Na-supported activity was comparable with nontransfected cells;
            # K-supported exchange remained substantial (2025 Figs. 5-6).
            "T448I_T756A": cls(
                0.0,
                0.5,
                "T448I-T756A qualitative Na-collapse/K-retention encoding",
            ),
        }
        try:
            return variants[name]
        except KeyError as exc:
            raise KeyError(f"unknown AE4 assay variant {name!r}") from exc


@dataclass(frozen=True)
class AE4EffectiveRates:
    """Two-node multigraph rates in s^-1 per carrier.

    ``cl_o_to_i`` and branch ``*_i_to_o`` rates point in the positive
    chloride-loading cycle direction.  Reverse rates use the opposite suffix.
    """

    cl_o_to_i_s: float
    cl_i_to_o_s: float
    na_i_to_o_s: float
    na_o_to_i_s: float
    k_i_to_o_s: float
    k_o_to_i_s: float
    gate_na: float
    gate_k: float

    @property
    def outward_to_inward_total_s(self) -> float:
        return self.cl_o_to_i_s + self.na_o_to_i_s + self.k_o_to_i_s

    @property
    def inward_to_outward_total_s(self) -> float:
        return self.cl_i_to_o_s + self.na_i_to_o_s + self.k_i_to_o_s

    @property
    def relaxation_rate_s(self) -> float:
        return self.outward_to_inward_total_s + self.inward_to_outward_total_s

    @property
    def relaxation_time_s(self) -> float:
        return 1.0 / self.relaxation_rate_s


@dataclass(frozen=True)
class AE4Flux:
    """AE4 QSS result with explicit species and conserved-coordinate sources."""

    model_id: str
    units: str
    j_na_fmol_s: float
    j_k_fmol_s: float
    source_na_i_fmol_s: float
    source_k_i_fmol_s: float
    source_cl_i_fmol_s: float
    source_hco3_i_fmol_s: float
    source_tic_i_fmol_s: float
    source_alkalinity_i_fmol_s: float
    affinity_na: float
    affinity_k: float
    entropy_production_over_r_fmol_s: float
    transported_charge_equivalents_fmol_s: float
    transported_current_A: float
    occupancy_outward: float
    occupancy_inward: float
    carrier_relaxation_time_s: float
    qss_residual_s: float
    local_detailed_balance_residual: float
    regulation_gain: float
    mixed_bath_na_attempt_fraction: float
    branch_modulation_label: str
    effective_rates: AE4EffectiveRates

    @property
    def intracellular_species_sources_fmol_s(self) -> Mapping[str, float]:
        """Explicit mobile-species sources; bicarbonate is not double counted."""

        return MappingProxyType(
            {
                "na_i": self.source_na_i_fmol_s,
                "k_i": self.source_k_i_fmol_s,
                "cl_i": self.source_cl_i_fmol_s,
                "hco3_i": self.source_hco3_i_fmol_s,
            }
        )

    @property
    def intracellular_conserved_sources_fmol_s(self) -> Mapping[str, float]:
        """Sources for ``(Na,K,Cl,TIC,alkalinity)`` amount coordinates."""

        return MappingProxyType(
            {
                "na_i": self.source_na_i_fmol_s,
                "k_i": self.source_k_i_fmol_s,
                "cl_i": self.source_cl_i_fmol_s,
                "tic_i": self.source_tic_i_fmol_s,
                "alkalinity_i": self.source_alkalinity_i_fmol_s,
            }
        )


def ae4_branch_affinity(environment: AE4Environment, branch: Literal["na", "k"]) -> float:
    """Dimensionless affinity for one positive chloride-loading cycle."""

    if branch == "na":
        cation_i, cation_o = environment.na_i_mM, environment.na_o_mM
    elif branch == "k":
        cation_i, cation_o = environment.k_i_mM, environment.k_o_mM
    else:
        raise ValueError(f"unknown AE4 branch {branch!r}")
    return (
        math.log(environment.cl_o_mM / environment.cl_i_mM)
        + math.log(cation_i / cation_o)
        + 2.0 * math.log(environment.hco3_i_mM / environment.hco3_o_mM)
    )


def ae4_reversal_cl_i_mM(
    environment: AE4Environment, branch: Literal["na", "k"]
) -> float:
    """Intracellular chloride at zero affinity with other activities fixed."""

    if branch == "na":
        cation_i, cation_o = environment.na_i_mM, environment.na_o_mM
    elif branch == "k":
        cation_i, cation_o = environment.k_i_mM, environment.k_o_mM
    else:
        raise ValueError(f"unknown AE4 branch {branch!r}")
    return (
        environment.cl_o_mM
        * (cation_i / cation_o)
        * (environment.hco3_i_mM / environment.hco3_o_mM) ** 2
    )


def _symmetric_rates(attempt_rate_s: float, log_forward_reverse_ratio: float) -> tuple[float, float]:
    half = 0.5 * log_forward_reverse_ratio
    if abs(half) > 700.0:
        raise OverflowError("AE4 activity ratio lies outside the finite rate domain")
    return attempt_rate_s * math.exp(half), attempt_rate_s * math.exp(-half)


def ae4_effective_rates(
    environment: AE4Environment,
    parameters: AE4Parameters,
    *,
    branch_modulation: AE4BranchModulation | None = None,
) -> AE4EffectiveRates:
    """Construct rates whose product ratios equal the aqueous affinities."""

    modulation = branch_modulation or AE4BranchModulation()
    gate_na = 1.0
    gate_k = 1.0
    if parameters.cooperative_gate is not None:
        gate_na = parameters.cooperative_gate.factor("na", environment.na_o_mM)
        gate_k = parameters.cooperative_gate.factor("k", environment.k_o_mM)

    # DERIVED_CONSTRAINT: the internal state-energy gauge is zero.  The two
    # edge-pair log ratios sum exactly to the complete chemical affinity.
    log_cl_ratio = math.log(environment.cl_o_mM / environment.cl_i_mM)
    log_na_return_ratio = (
        math.log(environment.na_i_mM / environment.na_o_mM)
        + 2.0 * math.log(environment.hco3_i_mM / environment.hco3_o_mM)
    )
    log_k_return_ratio = (
        math.log(environment.k_i_mM / environment.k_o_mM)
        + 2.0 * math.log(environment.hco3_i_mM / environment.hco3_o_mM)
    )
    cl_forward, cl_reverse = _symmetric_rates(
        parameters.common_cl_attempt_rate_s, log_cl_ratio
    )
    na_forward, na_reverse = _symmetric_rates(
        parameters.na_loaded_attempt_rate_s * modulation.na_scale * gate_na,
        log_na_return_ratio,
    )
    k_forward, k_reverse = _symmetric_rates(
        parameters.k_loaded_attempt_rate_s * modulation.k_scale * gate_k,
        log_k_return_ratio,
    )
    return AE4EffectiveRates(
        cl_o_to_i_s=cl_forward,
        cl_i_to_o_s=cl_reverse,
        na_i_to_o_s=na_forward,
        na_o_to_i_s=na_reverse,
        k_i_to_o_s=k_forward,
        k_o_to_i_s=k_reverse,
        gate_na=gate_na,
        gate_k=gate_k,
    )


def ae4_carrier_rhs_s(occupancy_inward: float, rates: AE4EffectiveRates) -> float:
    """Reduced carrier-occupancy RHS, for QSS time-scale diagnostics only.

    This scalar RHS does not supply the transient bound-substrate bookkeeping
    needed for an exact explicit carrier-state whole-cell model.  If its
    relaxation is not independently fast, use an edge-resolved carrier model
    with bound-ion inventories instead of integrating this diagnostic alone.
    """

    if not math.isfinite(occupancy_inward) or not 0.0 <= occupancy_inward <= 1.0:
        raise ValueError("occupancy_inward must lie in [0, 1]")
    occupancy_outward = 1.0 - occupancy_inward
    return (
        rates.outward_to_inward_total_s * occupancy_outward
        - rates.inward_to_outward_total_s * occupancy_inward
    )


def evaluate_ae4_qss(
    environment: AE4Environment,
    parameters: AE4Parameters,
    *,
    regulation_gain: float = 1.0,
    branch_modulation: AE4BranchModulation | None = None,
) -> AE4Flux:
    """Evaluate the smallest shared-pool AE4 module at carrier QSS.

    ``regulation_gain`` is a dimensionless active-capacity output from the
    external beta/cAMP/PKA subsystem.  It must be frozen from regulatory/WT
    evidence before any knockout secretion prediction.
    """

    if not math.isfinite(regulation_gain) or regulation_gain < 0.0:
        raise ValueError("regulation_gain must be finite and nonnegative")
    modulation = branch_modulation or AE4BranchModulation()
    rates = ae4_effective_rates(
        environment, parameters, branch_modulation=modulation
    )
    sigma = rates.relaxation_rate_s
    occupancy_outward = rates.inward_to_outward_total_s / sigma
    occupancy_inward = rates.outward_to_inward_total_s / sigma

    j_na_per_carrier_s = (
        rates.na_i_to_o_s * occupancy_inward
        - rates.na_o_to_i_s * occupancy_outward
    )
    j_k_per_carrier_s = (
        rates.k_i_to_o_s * occupancy_inward
        - rates.k_o_to_i_s * occupancy_outward
    )
    active_carrier_fmol = parameters.carrier_amount_fmol * regulation_gain
    j_na = active_carrier_fmol * j_na_per_carrier_s
    j_k = active_carrier_fmol * j_k_per_carrier_s
    j_total = j_na + j_k

    source_na = -j_na
    source_k = -j_k
    source_cl = j_total
    source_hco3 = -2.0 * j_total
    source_tic = source_hco3
    source_alkalinity = source_hco3
    charge_rate = math.fsum((source_na, source_k, -source_cl, -source_alkalinity))
    current_A = FARADAY_C_PER_MOL * FMOL_TO_MOL * charge_rate

    affinity_na = ae4_branch_affinity(environment, "na")
    affinity_k = ae4_branch_affinity(environment, "k")
    entropy_over_r = j_na * affinity_na + j_k * affinity_k
    qss_residual = abs(ae4_carrier_rhs_s(occupancy_inward, rates))

    ldb_residuals = []
    if rates.na_i_to_o_s > 0.0 and rates.na_o_to_i_s > 0.0:
        ldb_residuals.append(
            abs(
                math.log(rates.cl_o_to_i_s / rates.cl_i_to_o_s)
                + math.log(rates.na_i_to_o_s / rates.na_o_to_i_s)
                - affinity_na
            )
        )
    if rates.k_i_to_o_s > 0.0 and rates.k_o_to_i_s > 0.0:
        ldb_residuals.append(
            abs(
                math.log(rates.cl_o_to_i_s / rates.cl_i_to_o_s)
                + math.log(rates.k_i_to_o_s / rates.k_o_to_i_s)
                - affinity_k
            )
        )
    ldb_residual = max(ldb_residuals, default=0.0)

    model_id = (
        "SR2_SHARED_112_QSS"
        if parameters.cooperative_gate is None
        else "SR5_FAST_GATE_112_QSS"
    )
    return AE4Flux(
        model_id=model_id,
        units="fmol/s",
        j_na_fmol_s=j_na,
        j_k_fmol_s=j_k,
        source_na_i_fmol_s=source_na,
        source_k_i_fmol_s=source_k,
        source_cl_i_fmol_s=source_cl,
        source_hco3_i_fmol_s=source_hco3,
        source_tic_i_fmol_s=source_tic,
        source_alkalinity_i_fmol_s=source_alkalinity,
        affinity_na=affinity_na,
        affinity_k=affinity_k,
        entropy_production_over_r_fmol_s=entropy_over_r,
        transported_charge_equivalents_fmol_s=charge_rate,
        transported_current_A=current_A,
        occupancy_outward=occupancy_outward,
        occupancy_inward=occupancy_inward,
        carrier_relaxation_time_s=rates.relaxation_time_s,
        qss_residual_s=qss_residual,
        local_detailed_balance_residual=ldb_residual,
        regulation_gain=regulation_gain,
        mixed_bath_na_attempt_fraction=parameters.mixed_bath_na_attempt_fraction,
        branch_modulation_label=modulation.label,
        effective_rates=rates,
    )


__all__ = (
    "AE4_BRANCH_SOURCE_STOICHIOMETRY",
    "AE4BranchModulation",
    "AE4EffectiveRates",
    "AE4Environment",
    "AE4FastCooperativeGate",
    "AE4Flux",
    "AE4Parameters",
    "ae4_branch_affinity",
    "ae4_carrier_rhs_s",
    "ae4_effective_rates",
    "ae4_reversal_cl_i_mM",
    "evaluate_ae4_qss",
)
