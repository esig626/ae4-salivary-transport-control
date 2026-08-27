"""Source-auditable AE4 mechanism candidates.

The common sign convention is an AE4 event rate ``J > 0`` for basolateral
chloride entry into the cell.  For a ``Cl:C:HCO3 = a:b:c`` cycle this adds

``(+a J, -b_Na J, -b_K J, -c J)``

to the intracellular ``(Cl, Na, K, HCO3)`` amount balances.  All
concentrations supplied to this module are in mM.  ``activity_scale`` is in
the chassis' amount/time unit, except for the explicitly historical and
published empirical laws where the polynomial rate constants retain their
source conventions and the fitted activity absorbs their scale.

This module is an independent implementation of documented equations.  It is
not a port of the archived MATLAB source.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from types import MappingProxyType
from typing import Any, Mapping, Optional, Sequence


@dataclass(frozen=True)
class AE4Context:
    """Concentrations plus an independent beta/cAMP/PKA activation input."""

    na_i: float
    k_i: float
    cl_i: float
    hco3_i: float
    na_e: float
    k_e: float
    cl_e: float
    hco3_e: float
    pka_activation: float = 0.0


@dataclass(frozen=True)
class AE4Parameters:
    """AE4-only parameters shared by candidate evaluators.

    ``k_forward`` and ``k_reverse`` are optional because C1 and C2 have
    different source defaults.  For C1, ``k_reverse`` denotes the historical
    pre-fold coefficient and is multiplied by ``historical_reverse_fold``.
    For thermodynamic laws, ``equilibrium_constant`` enters the affinity as
    ``log(K_eq * forward_product / reverse_product)``.
    """

    activity_scale: float = 1.0
    k_forward: Optional[float] = None
    k_reverse: Optional[float] = None
    reference_concentration: float = 1.0
    equilibrium_constant: float = 1.0
    na_weight: float = 1.0
    k_weight: float = 1.0
    ec50_na: float = 49.0
    ec50_k: float = 62.0
    hill_na: float = 2.0
    hill_k: float = 1.8
    gate_floor_na: float = 0.20
    gate_floor_k: float = 0.25
    historical_external_hco3: float = 21.0
    historical_reverse_fold: float = 118.7219
    saturation_strength: float = 1.0
    pka_stimulated_fold: float = 3.0


@dataclass(frozen=True)
class AE4Contribution:
    """AE4 contribution to the four intracellular amount balances."""

    cycle_flux: float
    cl_i: float
    na_i: float
    k_i: float
    hco3_i: float
    transported_charge: float
    affinity: Optional[float]
    valid: bool
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    @property
    def cl_term(self) -> float:
        return self.cl_i

    @property
    def na_term(self) -> float:
        return self.na_i

    @property
    def k_term(self) -> float:
        return self.k_i

    @property
    def hco3_term(self) -> float:
        return self.hco3_i


@dataclass(frozen=True)
class CandidateSpec:
    """Machine-readable provenance and complexity record for one candidate."""

    candidate_id: str
    family: str
    label: str
    stoichiometry: Optional[tuple[int, int, int]]
    cation_mode: str
    rate_law: str
    free_parameters: tuple[str, ...]
    independent_constraints: tuple[str, ...]
    evidence: tuple[str, ...]
    new_assumptions: tuple[str, ...]
    thermodynamic_sign_controlled: bool
    active: bool = True
    rejection_reason: Optional[str] = None
    equivalent_to: Optional[str] = None
    core_candidate_id: Optional[str] = None
    stimulus_mode: str = "none"
    investigation_round: int = 1
    assay_status: str = "evaluate against declared evidence"


_PRIMARY_2016 = "Peña-Münzenmayer et al. 2016, DOI 10.1085/jgp.201611571"
_PUBLISHED_2018 = "Vera-Sigüenza et al. 2018, DOI 10.1007/s11538-017-0370-6"
_HISTORICAL = "immutable historical implementation evidence (Task 11)"
_PKA_2021 = "Peña-Münzenmayer et al. 2021, DOI 10.1152/ajpgi.00145.2021"

# Predeclared source-assay sensitivity values, not fits to secretion output.
# The paper states the recombinant forskolin effect explicitly (~25%); the WT
# total-exchanger and native Ae2-knockout/Ae4-only values are approximate
# readings of its Fig. 1B and Fig. 1C, respectively.
PKA_FOLD_ASSAY_VALUES: tuple[float, ...] = (1.25, 1.6, 3.0)


def _specs() -> dict[str, CandidateSpec]:
    common_primary = (
        f"{_PRIMARY_2016}: direct Na and K transport, electroneutrality, and reversibility",
    )
    stoich_evidence = (
        f"{_PRIMARY_2016}: explicitly considered electroneutral 1:1:2, 1:2:3, and 2:1:3",
    )
    specs = [
        CandidateSpec(
            "C1",
            "C1 historical Na-only",
            "Historical Na-only empirical polynomial",
            (1, 1, 2),
            "Na only",
            "source empirical forward-minus-reverse polynomial",
            ("activity_scale",),
            ("historical WT closure only; not primary transport evidence",),
            (f"{_HISTORICAL}: Na-only law, 21 mM reverse bicarbonate, folded reverse coefficient",),
            (),
            False,
        ),
        CandidateSpec(
            "C2",
            "C2 published pooled Na+K",
            "Exact printed pooled-cation polynomial and balance partition",
            (1, 1, 2),
            "pooled Na+K; intracellular mole-fraction balance insertion",
            "2018 printed forward-minus-reverse polynomial",
            ("activity_scale",),
            ("WT closure; printed k+ and k-",),
            (f"{_PUBLISHED_2018}: Eq. 31 and Eqs. 2-4",),
            (
                "pooled term treated as a constitutive law",
                "same intracellular partition retained even if the net cycle reverses",
            ),
            False,
        ),
        CandidateSpec(
            "C3",
            "C3 explicit parallel Na/K",
            "Explicit Na-coupled and K-coupled 1:1:2 branches",
            (1, 1, 2),
            "parallel Na and K branches",
            "empirical branchwise mass action",
            ("activity_scale",),
            ("Na:K capacity weights fixed at the independent Rmax ratio 1.5:1.6",),
            common_primary,
            (
                "independent parallel cycles; common forward/reverse constants",
                "Rmax ratio is treated as a branch-capacity ratio",
            ),
            False,
        ),
        CandidateSpec(
            "C1_PKA_COMMON",
            "C1 stimulus-dependent",
            "Historical C1 with a common beta/cAMP/PKA activity gate",
            (1, 1, 2),
            "Na only",
            "historical C1 turnover multiplied by 1 + u_PKA(F_PKA - 1)",
            ("activity_scale",),
            (
                "resting WT closure",
                "2021 assay fold envelope 1.25, approximately 1.6, approximately 3.0",
            ),
            (f"{_PKA_2021}: beta-adrenergic/cAMP/PKA stimulation increases AE4 activity",),
            (
                "macroscopic activation is represented as common turnover gating",
                "full assay activation is mapped to the model stimulus step",
            ),
            False,
            core_candidate_id="C1",
            stimulus_mode="pka_common",
            investigation_round=2,
            assay_status="inadmissible as a complete mechanism: Na-only conflicts with direct K transport",
        ),
        CandidateSpec(
            "C1_PKA_K_RECRUIT",
            "C1 stimulus-dependent",
            "Historical basal C1 plus a stimulus-recruited K-coupled branch",
            (1, 1, 2),
            "basal Na branch plus beta/cAMP-recruited K branch",
            "C1 basal flux plus u_PKA(F_PKA - 1) times that flux assigned to K",
            ("activity_scale",),
            (
                "resting WT closure",
                "2021 activity-fold envelope",
                "2016 direct K transport",
            ),
            (
                f"{_PKA_2021}: PKA increases macroscopic AE4 activity",
                f"{_PRIMARY_2016}: AE4 directly transports K",
            ),
            (
                "PKA selectively recruits the K branch",
                "recruited K flux is direction-locked to the empirical C1 flux",
                "resting AE4 is taken to be Na-only",
                "full assay activation is mapped to the model stimulus step",
            ),
            False,
            core_candidate_id="C1",
            stimulus_mode="pka_k_recruit",
            investigation_round=2,
            assay_status=(
                "novel Round-2 diagnostic: PKA activation and K transport are separately "
                "supported, but selective K recruitment is unmeasured"
            ),
        ),
        CandidateSpec(
            "C3_PKA_COMMON",
            "C3 stimulus-dependent",
            "Explicit parallel Na/K C3 with a common beta/cAMP/PKA activity gate",
            (1, 1, 2),
            "parallel transported Na and K branches",
            "C3 macroscopic source multiplied by 1 + u_PKA(F_PKA - 1)",
            ("activity_scale",),
            (
                "Na:K capacity weights fixed at the independent Rmax ratio 1.5:1.6",
                "2021 assay fold envelope 1.25, approximately 1.6, approximately 3.0",
            ),
            (
                f"{_PRIMARY_2016}: direct Na and K transport",
                f"{_PKA_2021}: beta-adrenergic/cAMP/PKA stimulation increases AE4 activity",
            ),
            (
                "the same macroscopic activation factor applies to both cation branches",
                "the activity change is not assigned uniquely to turnover, active number, or stoichiometry",
                "full assay activation is mapped to the model stimulus step",
            ),
            False,
            core_candidate_id="C3",
            stimulus_mode="pka_common",
            investigation_round=2,
            assay_status=(
                "direct-source combination of measured Na/K transport and PKA activation; "
                "macroscopic gate interpretation remains unresolved"
            ),
        ),
        CandidateSpec(
            "C4",
            "C4 cooperative occupancy",
            "Measured Hill occupancy placed in the transported-cation step",
            (1, 1, 2),
            "parallel Na and K branches",
            "directional Hill occupancy difference",
            ("activity_scale",),
            ("Na: EC50 49 mM, nH 2.0; K: EC50 62 mM, nH 1.8",),
            (f"{_PRIMARY_2016}: Fig. 9 Hill fits",),
            (
                "Hill response is used as directional occupancy",
                "one cation is transported despite a cooperative response",
            ),
            False,
        ),
        CandidateSpec(
            "C5",
            "C5 cation-dependent allosteric gating",
            "Transported 1:1:2 cycle plus an additional Hill gate",
            (1, 1, 2),
            "parallel transported Na/K plus branch gate",
            "thermodynamic mass action multiplied by a positive donor-side Hill gate",
            ("activity_scale",),
            ("measured EC50, Hill slopes, and Rmin/Rmax ratios",),
            (f"{_PRIMARY_2016}: direct cation transport plus allosteric regulation explicitly left viable",),
            ("phenomenological donor-side gate; no microscopic gating scheme",),
            True,
        ),
        CandidateSpec(
            "C5_GATE_ONLY_REJECTED",
            "C5 cation gating without transport",
            "Pure allosteric 1Cl:1HCO3 exchange",
            (1, 0, 1),
            "regulatory cation only",
            "gated anion exchange",
            ("activity_scale",),
            ("cation dose-response only",),
            (f"{_PRIMARY_2016}: direct Na and K movement",),
            (),
            True,
            False,
            "direct Na and K flux measurements rule out a pure gating-only AE4 cycle",
        ),
        CandidateSpec(
            "C6_CARBONATE_UNTESTABLE",
            "C6 alternative electroneutral stoichiometry",
            "Na-dependent Cl/carbonate exchange left viable by the primary study",
            None,
            "direct cation transport; carbonate stoichiometry unspecified",
            "not constructible on the seven-state bicarbonate-only chassis",
            (),
            ("requires a measured carbonate species and explicit acid-base partition",),
            (f"{_PRIMARY_2016}: Na-dependent Cl/CO3 exchange could not be ruled out",),
            (),
            False,
            False,
            "source-viable chemistry, but the primary study gives no stoichiometry and the fixed chassis has no carbonate state",
        ),
        CandidateSpec(
            "C6_HCO3_SELF_EXCHANGE_UNTESTABLE",
            "C6 alternative electroneutral stoichiometry",
            "Low-Cl bicarbonate/Na-bicarbonate exchange mode",
            None,
            "direct Na transport; unknown x:y bicarbonate coefficients",
            "not constructible without the unknown coefficients",
            (),
            ("requires the experimentally unresolved x:y bicarbonate stoichiometry",),
            (f"{_PRIMARY_2016}: Figs. 5-6, nonphysiological Cl-free substrate substitution",),
            (),
            False,
            False,
            "special Cl-free mode with unresolved stoichiometry; not a defined physiological Cl-cycle candidate",
        ),
    ]

    for suffix, stoich in (("112", (1, 1, 2)), ("123", (1, 2, 3)), ("213", (2, 1, 3))):
        specs.append(
            CandidateSpec(
                f"C6_{suffix}",
                "C6 alternative electroneutral stoichiometry",
                f"Parallel Na/K {stoich[0]}:{stoich[1]}:{stoich[2]} cycle",
                stoich,
                "parallel Na and K branches",
                "dimensionless reversible mass action",
                ("activity_scale",),
                ("WT closure; primary transport direction and stoichiometry set",),
                common_primary + stoich_evidence,
                (
                    "pure-Na and pure-K branches in mixed-cation conditions",
                    "symmetric passive equilibrium constant K_eq=1 unless independently changed",
                ),
                True,
            )
        )
        specs.append(
            CandidateSpec(
                f"C6_{suffix}_PKA_COMMON",
                "C6 stimulus-dependent",
                (
                    "Reversible mixed-cation "
                    f"{stoich[0]}:{stoich[1]}:{stoich[2]} core with a common "
                    "beta/cAMP/PKA gate"
                ),
                stoich,
                "parallel transported Na and K branches",
                "dimensionless reversible mass action multiplied by 1 + u_PKA(F_PKA - 1)",
                ("activity_scale",),
                (
                    "Na:K capacity weights fixed at the independent Rmax ratio 1.5:1.6",
                    "passive equilibrium constant fixed at one",
                    "2021 assay fold envelope 1.25, approximately 1.6, approximately 3.0",
                ),
                (
                    f"{_PRIMARY_2016}: direct Na/K transport, "
                    f"{stoich[0]}:{stoich[1]}:{stoich[2]} is source-viable, "
                    "and exchange is reversible",
                    f"{_PKA_2021}: beta-adrenergic/cAMP/PKA stimulation increases AE4 activity",
                ),
                (
                    "pure-Na and pure-K reversible branches represent mixed-cation transport",
                    "the same positive macroscopic activation factor applies to both branches",
                    "full assay activation is mapped to the model stimulus step",
                ),
                True,
                core_candidate_id=f"C6_{suffix}",
                stimulus_mode="pka_common",
                investigation_round=2,
                assay_status=(
                    (
                        "simplest assay- and thermodynamics-admissible source combination; "
                        if suffix == "112"
                        else "assay- and thermodynamics-admissible stoichiometric alternative; "
                    )
                    + "microscopic activation mechanism remains unresolved"
                ),
            )
        )
        specs.append(
            CandidateSpec(
                f"C7_{suffix}",
                "C7 thermodynamically reversible",
                f"Affinity-controlled reversible {stoich[0]}:{stoich[1]}:{stoich[2]} cycle",
                stoich,
                "parallel Na and K branches",
                "2 sinh(affinity/2)",
                ("activity_scale",),
                ("WT closure; independently specified equilibrium constant",),
                common_primary + stoich_evidence,
                (
                    "symmetric affinity law",
                    "pure-Na and pure-K branches in mixed-cation conditions",
                ),
                True,
            )
        )
        for law, label in (("MA", "mass action"), ("SAT", "substrate-saturating")):
            specs.append(
                CandidateSpec(
                    f"C8_{suffix}_{law}",
                    "C8 saturation versus mass action",
                    f"Reversible {label} {stoich[0]}:{stoich[1]}:{stoich[2]} cycle",
                    stoich,
                    "parallel Na and K branches",
                    "dimensionless reversible mass action" if law == "MA" else "symmetric occupancy denominator",
                    ("activity_scale",) if law == "MA" else ("activity_scale", "saturation_strength"),
                    ("WT closure; cation dose-response for saturation scale",),
                    common_primary + stoich_evidence,
                    (
                        "pure-Na and pure-K branches in mixed-cation conditions",
                        (
                            "common substrate saturation denominator; sigma=1 is an "
                            "unconstrained Round-1 default"
                            if law == "SAT"
                            else "unsaturated substrate products"
                        ),
                    ),
                    True,
                    equivalent_to=f"C6_{suffix}" if law == "MA" else None,
                )
            )
    return {spec.candidate_id: spec for spec in specs}


CANDIDATE_SPECS: Mapping[str, CandidateSpec] = MappingProxyType(_specs())


def candidate_specs(*, include_rejected: bool = True) -> tuple[CandidateSpec, ...]:
    """Return the stable candidate registry in construction order."""

    return tuple(
        spec for spec in CANDIDATE_SPECS.values() if include_rejected or spec.active
    )


def _get(context: Any, name: str) -> float:
    try:
        value = float(getattr(context, name))
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError(f"context must supply numeric {name}") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return value


def _validate_parameters(parameters: AE4Parameters) -> None:
    positive = (
        "reference_concentration",
        "equilibrium_constant",
        "ec50_na",
        "ec50_k",
        "hill_na",
        "hill_k",
        "historical_external_hco3",
        "historical_reverse_fold",
        "saturation_strength",
    )
    for name in positive:
        value = float(getattr(parameters, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and strictly positive")
    if not math.isfinite(float(parameters.activity_scale)) or parameters.activity_scale < 0.0:
        raise ValueError("activity_scale must be finite and nonnegative")
    if (
        not math.isfinite(float(parameters.pka_stimulated_fold))
        or parameters.pka_stimulated_fold < 1.0
    ):
        raise ValueError("pka_stimulated_fold must be finite and at least one")
    for name in ("na_weight", "k_weight"):
        value = float(getattr(parameters, name))
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if parameters.na_weight + parameters.k_weight <= 0.0:
        raise ValueError("at least one cation branch weight must be positive")
    for name in ("gate_floor_na", "gate_floor_k"):
        value = float(getattr(parameters, name))
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must lie in [0, 1]")
    for name in ("k_forward", "k_reverse"):
        value = getattr(parameters, name)
        if value is not None and (not math.isfinite(float(value)) or float(value) <= 0.0):
            raise ValueError(f"{name} must be finite and strictly positive when supplied")


def _context_values(context: Any) -> dict[str, float]:
    return {
        name: _get(context, name)
        for name in ("na_i", "k_i", "cl_i", "hco3_i", "na_e", "k_e", "cl_e", "hco3_e")
    }


def _resolve_stimulus_activation(
    context: Any, stimulus_activation: Optional[float]
) -> float:
    """Resolve the independent beta/cAMP/PKA input on the closed interval.

    The default is deliberately neutral.  A chassis that explicitly models
    beta-adrenergic stimulation may expose ``pka_activation`` on its context;
    an explicit evaluator keyword takes precedence.  Calcium is never used as
    a proxy because the source study distinguishes muscarinic/Ca signalling
    from beta/cAMP/PKA activation.
    """

    raw = (
        getattr(context, "pka_activation", 0.0)
        if stimulus_activation is None
        else stimulus_activation
    )
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError("stimulus_activation must be numeric") from exc
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("stimulus_activation must lie in [0, 1]")
    return value


def _hill(concentration: float, ec50: float, exponent: float) -> float:
    ratio = (ec50 / concentration) ** exponent
    return 1.0 / (1.0 + ratio)


def _gate(concentration: float, ec50: float, exponent: float, floor: float) -> float:
    return floor + (1.0 - floor) * _hill(concentration, ec50, exponent)


def _safe_affinity(forward: float, reverse: float, equilibrium_constant: float = 1.0) -> float:
    if forward <= 0.0 or reverse <= 0.0:
        raise ValueError("forward and reverse products must be positive")
    return math.log(equilibrium_constant) + math.log(forward) - math.log(reverse)


def _branch_products(
    values: Mapping[str, float],
    cation: str,
    stoichiometry: tuple[int, int, int],
    reference_concentration: float,
) -> tuple[float, float]:
    a, b, c = stoichiometry
    ref = reference_concentration
    ci = values[f"{cation}_i"] / ref
    ce = values[f"{cation}_e"] / ref
    forward = (values["cl_e"] / ref) ** a * ci**b * (values["hco3_i"] / ref) ** c
    reverse = (values["cl_i"] / ref) ** a * ce**b * (values["hco3_e"] / ref) ** c
    return forward, reverse


def _branch_weights(parameters: AE4Parameters) -> tuple[float, float]:
    total = parameters.na_weight + parameters.k_weight
    return parameters.na_weight / total, parameters.k_weight / total


def _make_contribution(
    *,
    stoichiometry: tuple[int, int, int],
    na_flux: float,
    k_flux: float,
    affinity: Optional[float],
    diagnostics: Mapping[str, Any],
    published_partition: Optional[tuple[float, float]] = None,
) -> AE4Contribution:
    a, b, c = stoichiometry
    total = na_flux + k_flux
    if published_partition is None:
        na_term = -b * na_flux
        k_term = -b * k_flux
    else:
        na_term = -b * published_partition[0] * total
        k_term = -b * published_partition[1] * total
    cl_term = a * total
    hco3_term = -c * total
    charge = na_term + k_term - cl_term - hco3_term
    finite = all(
        math.isfinite(value)
        for value in (total, cl_term, na_term, k_term, hco3_term, charge)
    )
    return AE4Contribution(
        cycle_flux=total,
        cl_i=cl_term,
        na_i=na_term,
        k_i=k_term,
        hco3_i=hco3_term,
        transported_charge=charge,
        affinity=affinity,
        valid=finite and abs(charge) <= 1e-10 * max(1.0, abs(total)),
        diagnostics=MappingProxyType(dict(diagnostics)),
    )


def _multiply_contribution(
    contribution: AE4Contribution,
    multiplier: float,
    *,
    stimulus_activation: float,
    stimulated_fold: float,
) -> AE4Contribution:
    """Apply a positive macroscopic activity gate without changing a signature."""

    diagnostics = dict(contribution.diagnostics)
    for branch in ("na", "k"):
        key = f"{branch}_flux"
        if key in diagnostics:
            diagnostics[key] = multiplier * float(diagnostics[key])
    diagnostics.update(
        {
            "pka_activation": stimulus_activation,
            "pka_stimulated_fold": stimulated_fold,
            "pka_activity_multiplier": multiplier,
            "pka_interpretation": "common macroscopic activity multiplier",
        }
    )
    return AE4Contribution(
        cycle_flux=multiplier * contribution.cycle_flux,
        cl_i=multiplier * contribution.cl_i,
        na_i=multiplier * contribution.na_i,
        k_i=multiplier * contribution.k_i,
        hco3_i=multiplier * contribution.hco3_i,
        transported_charge=multiplier * contribution.transported_charge,
        affinity=contribution.affinity,
        valid=contribution.valid and math.isfinite(multiplier) and multiplier > 0.0,
        diagnostics=MappingProxyType(diagnostics),
    )


def _pka_common_gate(
    contribution: AE4Contribution,
    p: AE4Parameters,
    stimulus_activation: float,
) -> AE4Contribution:
    if stimulus_activation == 0.0:
        # Preserve the resting core contribution bit-for-bit, including its
        # provenance diagnostics and equilibrium metadata.
        return contribution
    multiplier = 1.0 + stimulus_activation * (p.pka_stimulated_fold - 1.0)
    return _multiply_contribution(
        contribution,
        multiplier,
        stimulus_activation=stimulus_activation,
        stimulated_fold=p.pka_stimulated_fold,
    )


def _pka_recruited_k_branch(
    basal: AE4Contribution,
    p: AE4Parameters,
    stimulus_activation: float,
) -> AE4Contribution:
    """Novel localization: preserve basal C1 and assign recruited cycles to K.

    This diagnostic is electroneutral but not thermodynamically certified: the
    recruited K branch is direction-locked to the historical Na-only empirical
    flux rather than controlled by an independent K chemical affinity.
    """

    if stimulus_activation == 0.0:
        return basal

    recruited_k_flux = (
        stimulus_activation
        * (p.pka_stimulated_fold - 1.0)
        * basal.cycle_flux
    )
    return _make_contribution(
        stoichiometry=(1, 1, 2),
        na_flux=basal.cycle_flux,
        k_flux=recruited_k_flux,
        affinity=None,
        diagnostics={
            **dict(basal.diagnostics),
            "na_flux": basal.cycle_flux,
            "k_flux": recruited_k_flux,
            "pka_activation": stimulus_activation,
            "pka_stimulated_fold": p.pka_stimulated_fold,
            "pka_activity_multiplier": 1.0
            + stimulus_activation * (p.pka_stimulated_fold - 1.0),
            "pka_interpretation": "novel stimulus-recruited K branch",
            "thermodynamic_valid": False,
            "warning": "recruited K direction follows C1 rather than K affinity",
        },
    )


def _empirical_c1(values: Mapping[str, float], p: AE4Parameters) -> AE4Contribution:
    kf = 1.92e-2 if p.k_forward is None else p.k_forward
    kr_raw = 1.341840733667157e-5 if p.k_reverse is None else p.k_reverse
    kr = kr_raw * p.historical_reverse_fold
    forward = kf * values["cl_e"] * values["hco3_i"] ** 2 * values["na_i"]
    reverse = (
        kr
        * values["cl_i"]
        * p.historical_external_hco3**2
        * values["na_e"]
    )
    flux = p.activity_scale * (forward - reverse)
    affinity = _safe_affinity(forward, reverse)
    return _make_contribution(
        stoichiometry=(1, 1, 2),
        na_flux=flux,
        k_flux=0.0,
        affinity=affinity,
        diagnostics={
            "forward_propensity": forward,
            "reverse_propensity": reverse,
            "thermodynamic_valid": False,
            "source_convention": "historical Na-only; reverse bicarbonate fixed at 21 mM",
        },
    )


def _empirical_c2(values: Mapping[str, float], p: AE4Parameters) -> AE4Contribution:
    kf = 1.92e-2 if p.k_forward is None else p.k_forward
    kr = 1.3e-5 if p.k_reverse is None else p.k_reverse
    ci = values["na_i"] + values["k_i"]
    ce = values["na_e"] + values["k_e"]
    forward = kf * values["cl_e"] * values["hco3_i"] ** 2 * ci
    reverse = kr * values["cl_i"] * values["hco3_e"] ** 2 * ce
    flux = p.activity_scale * (forward - reverse)
    fractions = (values["na_i"] / ci, values["k_i"] / ci)
    return _make_contribution(
        stoichiometry=(1, 1, 2),
        na_flux=flux,
        k_flux=0.0,
        affinity=_safe_affinity(forward, reverse),
        diagnostics={
            "forward_propensity": forward,
            "reverse_propensity": reverse,
            "na_partition": fractions[0],
            "k_partition": fractions[1],
            "thermodynamic_valid": False,
            "source_convention": "exact printed pooled law and intracellular balance fractions",
        },
        published_partition=fractions,
    )


def _empirical_c3(values: Mapping[str, float], p: AE4Parameters) -> AE4Contribution:
    kf = 1.92e-2 if p.k_forward is None else p.k_forward
    kr = 1.3e-5 if p.k_reverse is None else p.k_reverse
    w_na, w_k = _branch_weights(p)
    branch_fluxes: dict[str, float] = {}
    affinities: dict[str, float] = {}
    for cation, weight in (("na", w_na), ("k", w_k)):
        forward = kf * values["cl_e"] * values["hco3_i"] ** 2 * values[f"{cation}_i"]
        reverse = kr * values["cl_i"] * values["hco3_e"] ** 2 * values[f"{cation}_e"]
        branch_fluxes[cation] = p.activity_scale * weight * (forward - reverse)
        affinities[cation] = _safe_affinity(forward, reverse)
    return _make_contribution(
        stoichiometry=(1, 1, 2),
        na_flux=branch_fluxes["na"],
        k_flux=branch_fluxes["k"],
        affinity=None,
        diagnostics={
            "na_flux": branch_fluxes["na"],
            "k_flux": branch_fluxes["k"],
            "affinity_na": affinities["na"],
            "affinity_k": affinities["k"],
            "thermodynamic_valid": False,
        },
    )


def _cooperative_occupancy(values: Mapping[str, float], p: AE4Parameters) -> AE4Contribution:
    w_na, w_k = _branch_weights(p)
    branch_fluxes: dict[str, float] = {}
    affinities: dict[str, float] = {}
    for cation, weight, ec50, exponent in (
        ("na", w_na, p.ec50_na, p.hill_na),
        ("k", w_k, p.ec50_k, p.hill_k),
    ):
        occupancy_i = _hill(values[f"{cation}_i"], ec50, exponent)
        occupancy_e = _hill(values[f"{cation}_e"], ec50, exponent)
        forward = (
            values["cl_e"]
            * values["hco3_i"] ** 2
            * occupancy_i
            / p.reference_concentration**3
        )
        reverse = (
            values["cl_i"]
            * values["hco3_e"] ** 2
            * occupancy_e
            / p.reference_concentration**3
        )
        branch_fluxes[cation] = p.activity_scale * weight * (forward - reverse)
        chemical_forward, chemical_reverse = _branch_products(
            values, cation, (1, 1, 2), p.reference_concentration
        )
        affinities[cation] = _safe_affinity(
            chemical_forward, chemical_reverse, p.equilibrium_constant
        )
    return _make_contribution(
        stoichiometry=(1, 1, 2),
        na_flux=branch_fluxes["na"],
        k_flux=branch_fluxes["k"],
        affinity=None,
        diagnostics={
            "na_flux": branch_fluxes["na"],
            "k_flux": branch_fluxes["k"],
            "affinity_na": affinities["na"],
            "affinity_k": affinities["k"],
            "thermodynamic_valid": False,
            "warning": "Hill occupancy can reverse at a different point from chemical affinity",
        },
    )


def _thermodynamic_branches(
    values: Mapping[str, float],
    p: AE4Parameters,
    stoichiometry: tuple[int, int, int],
    law: str,
    *,
    allosteric_gate: bool = False,
) -> AE4Contribution:
    w_na, w_k = _branch_weights(p)
    branch_fluxes: dict[str, float] = {}
    affinities: dict[str, float] = {}
    for cation, weight, ec50, exponent, floor in (
        ("na", w_na, p.ec50_na, p.hill_na, p.gate_floor_na),
        ("k", w_k, p.ec50_k, p.hill_k, p.gate_floor_k),
    ):
        forward, reverse = _branch_products(
            values, cation, stoichiometry, p.reference_concentration
        )
        sqrt_keq = math.sqrt(p.equilibrium_constant)
        forward_adjusted = sqrt_keq * forward
        reverse_adjusted = reverse / sqrt_keq
        affinity = _safe_affinity(forward, reverse, p.equilibrium_constant)
        if abs(affinity) <= 1e-12:
            # Products can be O(1e10) in mM-normalized high-order cycles.
            # Treat round-off at the analytically exact equilibrium as zero.
            numerator = 0.0
        else:
            numerator = forward_adjusted - reverse_adjusted
        if law == "mass_action":
            turnover = numerator
        elif law == "affinity":
            if abs(affinity) > 700.0:
                turnover = math.copysign(math.inf, affinity)
            else:
                turnover = 2.0 * math.sinh(0.5 * affinity)
        elif law == "saturating":
            denominator = 1.0 + p.saturation_strength * (
                forward_adjusted + reverse_adjusted
            )
            turnover = numerator / denominator
        else:
            raise ValueError(f"unknown thermodynamic law {law!r}")
        if allosteric_gate:
            donor = values[f"{cation}_i"] if affinity >= 0.0 else values[f"{cation}_e"]
            turnover *= _gate(donor, ec50, exponent, floor)
        branch_fluxes[cation] = p.activity_scale * weight * turnover
        affinities[cation] = affinity
    return _make_contribution(
        stoichiometry=stoichiometry,
        na_flux=branch_fluxes["na"],
        k_flux=branch_fluxes["k"],
        affinity=None,
        diagnostics={
            "na_flux": branch_fluxes["na"],
            "k_flux": branch_fluxes["k"],
            "affinity_na": affinities["na"],
            "affinity_k": affinities["k"],
            "thermodynamic_valid": all(
                branch_fluxes[cation] == 0.0
                or math.copysign(1.0, branch_fluxes[cation])
                == math.copysign(1.0, affinities[cation])
                for cation in ("na", "k")
            ),
            "rate_law": law,
            "allosteric_gate": allosteric_gate,
        },
    )


def evaluate_ae4(
    candidate_id: str,
    context: Any,
    parameters: Optional[AE4Parameters] = None,
    *,
    activity_scale: Optional[float] = None,
    stimulus_activation: Optional[float] = None,
) -> AE4Contribution:
    """Evaluate a registered AE4 candidate on a duck-typed concentration context.

    The context need only expose the eight attributes named by :class:`AE4Context`.
    The optional keyword ``activity_scale`` is the genotype/scenario multiplier
    used by the chassis (one for WT, zero for knockout).  It multiplies, rather
    than replaces, the calibrated capacity stored in ``parameters``.  The
    independent ``stimulus_activation`` input lies in ``[0, 1]`` and is used
    only by explicitly registered beta/cAMP/PKA candidates; if omitted, an
    optional context attribute ``pka_activation`` is read, otherwise rest
    (zero) is assumed.
    """

    try:
        spec = CANDIDATE_SPECS[candidate_id]
    except KeyError as exc:
        raise KeyError(f"unknown AE4 candidate {candidate_id!r}") from exc
    if not spec.active:
        raise ValueError(f"{candidate_id} is analytically rejected: {spec.rejection_reason}")
    p = parameters or AE4Parameters()
    if activity_scale is not None:
        p = replace(p, activity_scale=p.activity_scale * float(activity_scale))
    _validate_parameters(p)
    values = _context_values(context)
    pka_activation = _resolve_stimulus_activation(context, stimulus_activation)

    if candidate_id == "C1_PKA_COMMON":
        return _pka_common_gate(_empirical_c1(values, p), p, pka_activation)
    if candidate_id == "C1_PKA_K_RECRUIT":
        return _pka_recruited_k_branch(
            _empirical_c1(values, p), p, pka_activation
        )
    if candidate_id == "C3_PKA_COMMON":
        return _pka_common_gate(_empirical_c3(values, p), p, pka_activation)
    if candidate_id.startswith("C6_") and candidate_id.endswith("_PKA_COMMON"):
        suffix = candidate_id.split("_")[1]
        stoichiometry = {
            "112": (1, 1, 2),
            "123": (1, 2, 3),
            "213": (2, 1, 3),
        }[suffix]
        return _pka_common_gate(
            _thermodynamic_branches(values, p, stoichiometry, "mass_action"),
            p,
            pka_activation,
        )

    if candidate_id == "C1":
        return _empirical_c1(values, p)
    if candidate_id == "C2":
        return _empirical_c2(values, p)
    if candidate_id == "C3":
        return _empirical_c3(values, p)
    if candidate_id == "C4":
        return _cooperative_occupancy(values, p)
    if candidate_id == "C5":
        return _thermodynamic_branches(values, p, (1, 1, 2), "mass_action", allosteric_gate=True)

    if candidate_id.startswith(("C6_", "C7_", "C8_")):
        suffix = candidate_id.split("_")[1]
        stoichiometry = {
            "112": (1, 1, 2),
            "123": (1, 2, 3),
            "213": (2, 1, 3),
        }[suffix]
        if candidate_id.startswith("C7_"):
            law = "affinity"
        elif candidate_id.endswith("_SAT"):
            law = "saturating"
        else:
            law = "mass_action"
        return _thermodynamic_branches(values, p, stoichiometry, law)
    raise AssertionError(f"registered candidate has no evaluator: {candidate_id}")


@dataclass(frozen=True)
class AE4Mechanism:
    """Small callable wrapper suitable for injection into the fixed chassis."""

    candidate_id: str
    parameters: AE4Parameters = field(default_factory=AE4Parameters)

    @property
    def spec(self) -> CandidateSpec:
        return CANDIDATE_SPECS[self.candidate_id]

    def evaluate(
        self,
        context: Any,
        activity_scale: Optional[float] = None,
        stimulus_activation: Optional[float] = None,
    ) -> AE4Contribution:
        return evaluate_ae4(
            self.candidate_id,
            context,
            self.parameters,
            activity_scale=activity_scale,
            stimulus_activation=stimulus_activation,
        )

    def __call__(
        self,
        context: Any,
        activity_scale: Optional[float] = None,
        stimulus_activation: Optional[float] = None,
    ) -> AE4Contribution:
        return self.evaluate(
            context,
            activity_scale,
            stimulus_activation=stimulus_activation,
        )


def get_mechanism(
    candidate_id: str, parameters: Optional[AE4Parameters] = None
) -> AE4Mechanism:
    """Construct a checked mechanism wrapper from the stable registry."""

    if candidate_id not in CANDIDATE_SPECS:
        raise KeyError(f"unknown AE4 candidate {candidate_id!r}")
    if not CANDIDATE_SPECS[candidate_id].active:
        raise ValueError(
            f"{candidate_id} is analytically rejected: "
            f"{CANDIDATE_SPECS[candidate_id].rejection_reason}"
        )
    return AE4Mechanism(candidate_id, parameters or AE4Parameters())


def equilibrium_context(
    context: Any,
    cation: str,
    stoichiometry: tuple[int, int, int],
    *,
    equilibrium_constant: float = 1.0,
) -> AE4Context:
    """Return a context with ``hco3_e`` adjusted to zero one branch affinity.

    This exact helper is intended for tests and analytical diagnostics, not for
    physiological calibration.
    """

    if cation not in ("na", "k"):
        raise ValueError("cation must be 'na' or 'k'")
    values = _context_values(context)
    a, b, c = stoichiometry
    numerator = (
        equilibrium_constant
        * values["cl_e"] ** a
        * values[f"{cation}_i"] ** b
        * values["hco3_i"] ** c
    )
    denominator = values["cl_i"] ** a * values[f"{cation}_e"] ** b
    values["hco3_e"] = (numerator / denominator) ** (1.0 / c)
    return AE4Context(**values)


__all__: Sequence[str] = (
    "AE4Context",
    "AE4Parameters",
    "AE4Contribution",
    "CandidateSpec",
    "CANDIDATE_SPECS",
    "PKA_FOLD_ASSAY_VALUES",
    "candidate_specs",
    "evaluate_ae4",
    "AE4Mechanism",
    "get_mechanism",
    "equilibrium_context",
)
