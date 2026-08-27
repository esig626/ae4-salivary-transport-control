"""Thermodynamically constrained, explicit-state AE4 cycle models.

The implementation is independent of the immutable historical MATLAB code.
It represents each candidate as a finite continuous-time Markov chain (CTMC)
whose reversible edge rates obey local detailed balance.  Reservoir
stoichiometry is attached to edges rather than inferred from a scalar flux
law, so cycle affinities, transported ions, and mutant limits can be audited
directly.

Concentrations are converted to dimensionless activities using one common
reference concentration.  A positive cycle current denotes basolateral
chloride loading of the cell for all candidates.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray


ION_CHARGES: Mapping[str, int] = MappingProxyType(
    {
        "na_i": +1,
        "na_e": +1,
        "k_i": +1,
        "k_e": +1,
        "cl_i": -1,
        "cl_e": -1,
        "hco3_i": -1,
        "hco3_e": -1,
        "co3_i": -2,
        "co3_e": -2,
    }
)


@dataclass(frozen=True)
class ReservoirActivities:
    """Dimensionless aqueous activities supplied to a cycle network."""

    na_i: float
    k_i: float
    cl_i: float
    hco3_i: float
    na_e: float
    k_e: float
    cl_e: float
    hco3_e: float
    co3_i: float = math.nan
    co3_e: float = math.nan

    def as_mapping(self) -> Mapping[str, float]:
        return MappingProxyType(
            {
                "na_i": self.na_i,
                "k_i": self.k_i,
                "cl_i": self.cl_i,
                "hco3_i": self.hco3_i,
                "na_e": self.na_e,
                "k_e": self.k_e,
                "cl_e": self.cl_e,
                "hco3_e": self.hco3_e,
                "co3_i": self.co3_i,
                "co3_e": self.co3_e,
            }
        )

    def validate(self, required: Iterable[str]) -> None:
        values = self.as_mapping()
        for name in required:
            value = values[name]
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"activity {name} must be finite and positive")


@dataclass(frozen=True)
class ReversibleEdge:
    """One reversible CTMC edge with an auditable reservoir exchange.

    ``reservoir_stoichiometry`` is positive for an aqueous species consumed
    by the forward state transition and negative for a species released by
    it.  With state standard energies ``g``, the rate construction gives

    ``log(k_forward/k_reverse) = -(g_target-g_source) + sum(nu*log(a))``.

    The dimensionless state energies telescope around every closed cycle, so
    the product-of-rates identity is fixed exactly by the aqueous affinity.
    ``barrier_multiplier`` changes both directions equally and cannot create
    a fictitious equilibrium shift.
    """

    edge_id: str
    source: str
    target: str
    attempt_rate: float = 1.0
    reservoir_stoichiometry: Mapping[str, float] = field(default_factory=dict)
    barrier_multiplier: float = 1.0

    def __post_init__(self) -> None:
        if self.source == self.target:
            raise ValueError("a reversible edge must connect two distinct states")
        if not math.isfinite(self.attempt_rate) or self.attempt_rate <= 0.0:
            raise ValueError("attempt_rate must be finite and positive")
        if not math.isfinite(self.barrier_multiplier) or self.barrier_multiplier <= 0.0:
            raise ValueError("barrier_multiplier must be finite and positive")
        unknown = set(self.reservoir_stoichiometry) - set(ION_CHARGES)
        if unknown:
            raise ValueError(f"unknown reservoir species: {sorted(unknown)}")

    def rates(
        self,
        activities: Mapping[str, float],
        state_energies: Mapping[str, float],
    ) -> tuple[float, float]:
        delta_g0 = float(state_energies[self.target] - state_energies[self.source])
        log_forward = (
            math.log(self.attempt_rate * self.barrier_multiplier) - 0.5 * delta_g0
        )
        log_reverse = (
            math.log(self.attempt_rate * self.barrier_multiplier) + 0.5 * delta_g0
        )
        for name, coefficient in self.reservoir_stoichiometry.items():
            activity = float(activities[name])
            if not math.isfinite(activity) or activity <= 0.0:
                raise ValueError(f"activity {name} must be finite and positive")
            if coefficient > 0.0:
                log_forward += coefficient * math.log(activity)
            elif coefficient < 0.0:
                log_reverse += -coefficient * math.log(activity)
        return math.exp(log_forward), math.exp(log_reverse)


@dataclass(frozen=True)
class CycleDefinition:
    """State graph plus branch currents and intracellular source signatures."""

    model_id: str
    family: str
    states: tuple[str, ...]
    edges: tuple[ReversibleEdge, ...]
    state_energies: Mapping[str, float]
    branch_current_edges: Mapping[str, tuple[str, ...]]
    branch_cycle_edge_ids: Mapping[str, tuple[str, ...]]
    intracellular_source_vectors: Mapping[str, Mapping[str, float]]
    evidence_status: str
    free_parameter_names: tuple[str, ...]
    source_constraints: tuple[str, ...]
    new_assumptions: tuple[str, ...]
    pka_edge_ids: tuple[str, ...] = ()
    pka_switch_edge_ids: tuple[str, ...] = ()
    transported_cation_bound_states: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    transported_cation_counts: Mapping[str, int] = field(default_factory=dict)
    catalytic_bound_states: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    mutant_sensitive_edges: Mapping[str, tuple[str, ...]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(set(self.states)) != len(self.states):
            raise ValueError("cycle states must be unique")
        if set(self.state_energies) != set(self.states):
            raise ValueError("state_energies must cover every state exactly")
        edge_ids = {edge.edge_id for edge in self.edges}
        if len(edge_ids) != len(self.edges):
            raise ValueError("edge IDs must be unique")
        for edge in self.edges:
            if edge.source not in self.states or edge.target not in self.states:
                raise ValueError(f"edge {edge.edge_id} references an unknown state")
        declared_current_edges = {
            edge_id
            for branch_edges in self.branch_current_edges.values()
            for edge_id in branch_edges
        }
        if not declared_current_edges <= edge_ids:
            raise ValueError("branch current edge is absent from graph")
        if set(self.branch_current_edges) != set(self.intracellular_source_vectors):
            raise ValueError("every branch needs one intracellular source vector")
        if set(self.branch_cycle_edge_ids) != set(self.branch_current_edges):
            raise ValueError("every branch needs one declared affinity cycle")
        if not {
            edge_id
            for cycle_edges in self.branch_cycle_edge_ids.values()
            for edge_id in cycle_edges
        } <= edge_ids:
            raise ValueError("branch affinity cycle contains an absent edge")
        if not set(self.pka_edge_ids) <= edge_ids:
            raise ValueError("PKA edge is absent from graph")
        if not set(self.pka_switch_edge_ids) <= edge_ids:
            raise ValueError("PKA switch edge is absent from graph")
        if not {
            state
            for states in self.catalytic_bound_states.values()
            for state in states
        } <= set(self.states):
            raise ValueError("catalytic bound state is absent from graph")
        if set(self.transported_cation_bound_states) != set(self.transported_cation_counts):
            raise ValueError("transported bound-state counts are incomplete")
        if not {
            state
            for states in self.transported_cation_bound_states.values()
            for state in states
        } <= set(self.states):
            raise ValueError("transported cation-bound state is absent from graph")
        if not {
            edge_id
            for edge_ids in self.mutant_sensitive_edges.values()
            for edge_id in edge_ids
        } <= edge_ids:
            raise ValueError("mutant-sensitive edge is absent from graph")

    @property
    def required_reservoir_species(self) -> tuple[str, ...]:
        return tuple(
            sorted(
                {
                    name
                    for edge in self.edges
                    for name in edge.reservoir_stoichiometry
                }
            )
        )


@dataclass(frozen=True)
class CycleParameters:
    """Lumped kinetic parameters retained after microscopic nonidentifiability."""

    capacity: float = 1.0
    reference_concentration_mM: float = 1.0
    na_attempt_scale: float = 1.0
    k_attempt_scale: float = 1.0
    common_attempt_scale: float = 1.0
    ec50_na_mM: float = 49.0
    ec50_k_mM: float = 62.0
    hill_na: float = 2.0
    hill_k: float = 1.8
    gate_floor_na: float = 0.20
    gate_floor_k: float = 0.25
    pka_barrier_fold: float = 1.25
    pka_basal_activity: float = 1e-6
    pka_phosphorylation_scale: float = 1.0
    pka_dephosphorylation_scale: float = 1.0
    mutant_na_scale: float = 1.0
    mutant_k_scale: float = 1.0
    edge_barrier_overrides: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        positive = (
            "reference_concentration_mM",
            "na_attempt_scale",
            "k_attempt_scale",
            "common_attempt_scale",
            "ec50_na_mM",
            "ec50_k_mM",
            "hill_na",
            "hill_k",
            "pka_barrier_fold",
            "pka_phosphorylation_scale",
            "pka_dephosphorylation_scale",
        )
        for name in positive:
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for name in (
            "capacity",
            "gate_floor_na",
            "gate_floor_k",
            "mutant_na_scale",
            "mutant_k_scale",
            "pka_basal_activity",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if self.gate_floor_na > 1.0 or self.gate_floor_k > 1.0:
            raise ValueError("gate floors must lie in [0, 1]")
        if self.mutant_na_scale > 1.0 or self.mutant_k_scale > 1.0:
            raise ValueError("mutant branch scales must lie in [0, 1]")
        for edge_id, multiplier in self.edge_barrier_overrides.items():
            if not edge_id or not math.isfinite(multiplier) or multiplier <= 0.0:
                raise ValueError("edge barrier overrides must be named and positive")


@dataclass(frozen=True)
class CycleEvaluation:
    """Stationary CTMC solution and source terms at one chemical condition."""

    model_id: str
    occupancy: Mapping[str, float]
    edge_currents: Mapping[str, float]
    branch_currents: Mapping[str, float]
    intracellular_sources: Mapping[str, float]
    branch_affinities: Mapping[str, float]
    entropy_production: float
    transported_charge_rate: float
    generator_residual: float
    local_detailed_balance_residual: float
    pka_activation: float

    @property
    def cycle_flux(self) -> float:
        return float(sum(self.branch_currents.values()))


def hill_occupancy(concentration_mM: float, ec50_mM: float, hill: float) -> float:
    """Stable Hill occupancy used only as an edge-barrier modifier."""

    if concentration_mM <= 0.0 or ec50_mM <= 0.0 or hill <= 0.0:
        raise ValueError("Hill arguments must be positive")
    log_ratio = hill * math.log(ec50_mM / concentration_mM)
    if log_ratio > 700.0:
        return 0.0
    if log_ratio < -700.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(log_ratio))


def stationary_distribution(generator: NDArray[np.float64]) -> NDArray[np.float64]:
    """Solve ``pi Q=0, sum(pi)=1`` without selecting an eigenvector phase."""

    q = np.asarray(generator, dtype=float)
    if q.ndim != 2 or q.shape[0] != q.shape[1]:
        raise ValueError("generator must be square")
    n = q.shape[0]
    system = q.T.copy()
    rhs = np.zeros(n)
    system[-1, :] = 1.0
    rhs[-1] = 1.0
    try:
        occupancy = np.linalg.solve(system, rhs)
    except np.linalg.LinAlgError:
        occupancy, *_ = np.linalg.lstsq(system, rhs, rcond=None)
    if np.min(occupancy) < -1e-10 or not np.all(np.isfinite(occupancy)):
        raise ValueError("CTMC stationary solution is not a probability vector")
    occupancy = np.maximum(occupancy, 0.0)
    occupancy /= np.sum(occupancy)
    return occupancy


def _edge_with_scales(
    edge: ReversibleEdge,
    definition: CycleDefinition,
    parameters: CycleParameters,
    activities_mM: Mapping[str, float],
    pka_activation: float,
) -> ReversibleEdge:
    multiplier = parameters.common_attempt_scale
    edge_text = f"{edge.edge_id} {edge.source} {edge.target}".lower()
    if "na" in edge_text:
        multiplier *= parameters.na_attempt_scale
    if "k" in edge_text and "common" not in edge.edge_id.lower():
        multiplier *= parameters.k_attempt_scale
    if edge.edge_id in definition.mutant_sensitive_edges.get("na", ()):
        multiplier *= parameters.mutant_na_scale
    if edge.edge_id in definition.mutant_sensitive_edges.get("k", ()):
        multiplier *= parameters.mutant_k_scale
    if edge.edge_id in definition.pka_edge_ids:
        multiplier *= 1.0 + pka_activation * (parameters.pka_barrier_fold - 1.0)
    multiplier *= float(parameters.edge_barrier_overrides.get(edge.edge_id, 1.0))
    return replace(edge, barrier_multiplier=edge.barrier_multiplier * multiplier)


def evaluate_cycle(
    definition: CycleDefinition,
    concentrations_mM: ReservoirActivities,
    parameters: CycleParameters | None = None,
    *,
    pka_activation: float = 0.0,
) -> CycleEvaluation:
    """Evaluate one state graph at CTMC quasi-steady state.

    The CTMC occupancies are algebraic relative to whole-cell dynamics.  That
    approximation is explicit: the current is state-resolved, but transporter
    relaxation is assumed faster than the seven chassis states.
    """

    p = parameters or CycleParameters()
    if not 0.0 <= pka_activation <= 1.0 or not math.isfinite(pka_activation):
        raise ValueError("pka_activation must lie in [0, 1]")
    concentrations_mM.validate(definition.required_reservoir_species)
    raw = concentrations_mM.as_mapping()
    activities = {
        name: float(value) / p.reference_concentration_mM
        for name, value in raw.items()
        if math.isfinite(float(value))
    }
    index = {state: i for i, state in enumerate(definition.states)}
    effective_energies = dict(definition.state_energies)
    for species, states in definition.transported_cation_bound_states.items():
        ec50 = p.ec50_na_mM if species == "na" else p.ec50_k_mM
        count = definition.transported_cation_counts[species]
        delta = count * math.log(ec50 / p.reference_concentration_mM)
        for state in states:
            effective_energies[state] += delta
    # Explicit SR5 regulatory-bound states use the source-reported external
    # EC50 values as effective dissociation-energy constraints.  This maps an
    # assay-level lump to a state energy without claiming a microscopic site
    # count or fitting an intracellular concentration.
    for species, states in definition.catalytic_bound_states.items():
        ec50 = p.ec50_na_mM if species == "na" else p.ec50_k_mM
        delta = math.log(ec50 / p.reference_concentration_mM)
        for state in states:
            effective_energies[state] += delta
    n = len(definition.states)
    generator = np.zeros((n, n), dtype=float)
    rate_rows: dict[str, tuple[float, float]] = {}
    for edge in definition.edges:
        scaled = _edge_with_scales(
            edge, definition, p, raw, pka_activation
        )
        forward, reverse = scaled.rates(activities, effective_energies)
        if edge.edge_id in definition.pka_switch_edge_ids:
            forward *= p.pka_phosphorylation_scale * (
                p.pka_basal_activity + pka_activation
            )
            reverse *= p.pka_dephosphorylation_scale
        source = index[edge.source]
        target = index[edge.target]
        generator[source, target] += forward
        generator[target, source] += reverse
        rate_rows[edge.edge_id] = (forward, reverse)
    generator[np.diag_indices(n)] = -np.sum(generator, axis=1)
    occupancy = stationary_distribution(generator)
    edge_currents: dict[str, float] = {}
    ldb_residual = 0.0
    for edge in definition.edges:
        forward, reverse = rate_rows[edge.edge_id]
        current = occupancy[index[edge.source]] * forward - occupancy[index[edge.target]] * reverse
        edge_currents[edge.edge_id] = float(p.capacity * current)
        expected = -(
            effective_energies[edge.target]
            - effective_energies[edge.source]
        ) + sum(
            coefficient * math.log(activities[name])
            for name, coefficient in edge.reservoir_stoichiometry.items()
        )
        if edge.edge_id in definition.pka_switch_edge_ids:
            expected += math.log(
                p.pka_phosphorylation_scale
                * (p.pka_basal_activity + pka_activation)
                / p.pka_dephosphorylation_scale
            )
        ldb_residual = max(ldb_residual, abs(math.log(forward / reverse) - expected))
    branch_currents = {
        branch: float(sum(edge_currents[edge_id] for edge_id in edge_ids))
        for branch, edge_ids in definition.branch_current_edges.items()
    }
    sources = {"na_i": 0.0, "k_i": 0.0, "cl_i": 0.0, "hco3_i": 0.0, "co3_i": 0.0}
    affinities: dict[str, float] = {}
    entropy = 0.0
    for branch, current in branch_currents.items():
        source_vector = definition.intracellular_source_vectors[branch]
        for name, coefficient in source_vector.items():
            sources[name] += coefficient * current
        # The net reservoir consumption for a fundamental cycle is the
        # negative of its intracellular source for intracellular species and
        # the corresponding extracellular product.  Candidate builders store
        # a complete cycle affinity as a synthetic ``affinity_*`` mapping.
        cycle_edges = _branch_cycle_edges(definition, branch)
        affinity = sum(
            sum(
                coefficient * math.log(activities[name])
                for name, coefficient in edge.reservoir_stoichiometry.items()
            )
            for edge in cycle_edges
        )
        affinities[branch] = float(affinity)
        entropy += current * affinity
    transported_charge = sum(
        ION_CHARGES[name] * value
        for name, value in sources.items()
        if name in ION_CHARGES
    )
    generator_residual = float(np.max(np.abs(occupancy @ generator)))
    if entropy < -1e-10 * max(1.0, abs(entropy)):
        raise ValueError("passive CTMC produced negative entropy production")
    return CycleEvaluation(
        model_id=definition.model_id,
        occupancy=MappingProxyType(
            {state: float(occupancy[i]) for state, i in index.items()}
        ),
        edge_currents=MappingProxyType(edge_currents),
        branch_currents=MappingProxyType(branch_currents),
        intracellular_sources=MappingProxyType(sources),
        branch_affinities=MappingProxyType(affinities),
        entropy_production=float(max(entropy, 0.0)),
        transported_charge_rate=float(transported_charge),
        generator_residual=generator_residual,
        local_detailed_balance_residual=float(ldb_residual),
        pka_activation=float(pka_activation),
    )


def ctmc_snapshot(
    definition: CycleDefinition,
    concentrations_mM: ReservoirActivities,
    parameters: CycleParameters | None = None,
    *,
    pka_activation: float = 0.0,
) -> Mapping[str, object]:
    """Export ordered rates/generator for an implementation-independent check."""

    p = parameters or CycleParameters()
    concentrations_mM.validate(definition.required_reservoir_species)
    raw = concentrations_mM.as_mapping()
    activities = {
        name: float(value) / p.reference_concentration_mM
        for name, value in raw.items()
        if math.isfinite(float(value))
    }
    energies = dict(definition.state_energies)
    for species, states in definition.transported_cation_bound_states.items():
        ec50 = p.ec50_na_mM if species == "na" else p.ec50_k_mM
        count = definition.transported_cation_counts[species]
        for state in states:
            energies[state] += count * math.log(ec50 / p.reference_concentration_mM)
    for species, states in definition.catalytic_bound_states.items():
        ec50 = p.ec50_na_mM if species == "na" else p.ec50_k_mM
        for state in states:
            energies[state] += math.log(ec50 / p.reference_concentration_mM)
    index = {state: i for i, state in enumerate(definition.states)}
    generator = np.zeros((len(index), len(index)), dtype=float)
    edge_rows = []
    for edge in definition.edges:
        scaled = _edge_with_scales(edge, definition, p, raw, pka_activation)
        forward, reverse = scaled.rates(activities, energies)
        if edge.edge_id in definition.pka_switch_edge_ids:
            forward *= p.pka_phosphorylation_scale * (p.pka_basal_activity + pka_activation)
            reverse *= p.pka_dephosphorylation_scale
        source = index[edge.source]
        target = index[edge.target]
        generator[source, target] += forward
        generator[target, source] += reverse
        edge_rows.append(
            {
                "edge_id": edge.edge_id,
                "source": edge.source,
                "target": edge.target,
                "forward_rate": float(forward),
                "reverse_rate": float(reverse),
                "reservoir_stoichiometry": dict(edge.reservoir_stoichiometry),
            }
        )
    generator[np.diag_indices(len(index))] = -np.sum(generator, axis=1)
    occupancy = stationary_distribution(generator)
    for row in edge_rows:
        row["stationary_current_per_transporter"] = float(
            occupancy[index[str(row["source"])]] * float(row["forward_rate"])
            - occupancy[index[str(row["target"])]] * float(row["reverse_rate"])
        )
    return MappingProxyType(
        {
            "model_id": definition.model_id,
            "state_order": definition.states,
            "state_energies": MappingProxyType(energies),
            "generator": generator.copy(),
            "occupancy": occupancy.copy(),
            "edge_rows": tuple(MappingProxyType(row) for row in edge_rows),
            "concentrations_mM": concentrations_mM.as_mapping(),
            "pka_activation": pka_activation,
        }
    )


def _branch_cycle_edges(
    definition: CycleDefinition, branch: str
) -> tuple[ReversibleEdge, ...]:
    """Return the declared positive fundamental cycle for one branch."""

    edge_by_id = {edge.edge_id: edge for edge in definition.edges}
    try:
        return tuple(edge_by_id[edge_id] for edge_id in definition.branch_cycle_edge_ids[branch])
    except KeyError as exc:
        raise KeyError(f"no affinity cycle declared for branch {branch}") from exc


def _edge(
    edge_id: str,
    source: str,
    target: str,
    **reservoir_stoichiometry: float,
) -> ReversibleEdge:
    return ReversibleEdge(
        edge_id=edge_id,
        source=source,
        target=target,
        reservoir_stoichiometry=MappingProxyType(dict(reservoir_stoichiometry)),
    )


def _shared_bicarbonate_definition(
    model_id: str,
    family: str,
    *,
    branches: Sequence[str] = ("na", "k"),
    stoichiometry: tuple[int, int, int] = (1, 1, 2),
) -> CycleDefinition:
    cl_count, cation_count, hco3_count = stoichiometry
    if hco3_count != cl_count + cation_count:
        raise ValueError("bicarbonate candidate is not electroneutral")
    states = ["Eo", f"EoCl{cl_count}", f"EiCl{cl_count}", "Ei"]
    edges = [
        _edge("common_cl_bind", "Eo", f"EoCl{cl_count}", cl_e=+cl_count),
        _edge("common_cl_flip", f"EoCl{cl_count}", f"EiCl{cl_count}"),
        _edge("common_cl_release", f"EiCl{cl_count}", "Ei", cl_i=-cl_count),
    ]
    sources: dict[str, Mapping[str, float]] = {}
    current_edges: dict[str, tuple[str, ...]] = {}
    cycle_edges: dict[str, tuple[str, ...]] = {}
    transported_states: dict[str, tuple[str, ...]] = {}
    mutant_edges: dict[str, tuple[str, ...]] = {}
    for branch in branches:
        state_i = f"Ei{branch.title()}{cation_count}B{hco3_count}"
        state_o = f"Eo{branch.title()}{cation_count}B{hco3_count}"
        states.extend((state_i, state_o))
        cation_i = f"{branch}_i"
        cation_e = f"{branch}_e"
        edges.extend(
            (
                _edge(
                    f"{branch}_loaded_bind",
                    "Ei",
                    state_i,
                    **{cation_i: +cation_count, "hco3_i": +hco3_count},
                ),
                _edge(f"{branch}_loaded_flip", state_i, state_o),
                _edge(
                    f"{branch}_loaded_release",
                    state_o,
                    "Eo",
                    **{cation_e: -cation_count, "hco3_e": -hco3_count},
                ),
            )
        )
        sources[branch] = MappingProxyType(
            {
                "na_i": -float(cation_count) if branch == "na" else 0.0,
                "k_i": -float(cation_count) if branch == "k" else 0.0,
                "cl_i": float(cl_count),
                "hco3_i": -float(hco3_count),
            }
        )
        current_edges[branch] = (f"{branch}_loaded_flip",)
        transported_states[branch] = (state_i, state_o)
        mutant_edges[branch] = (f"{branch}_loaded_flip",)
        cycle_edges[branch] = (
            "common_cl_bind",
            "common_cl_flip",
            "common_cl_release",
            f"{branch}_loaded_bind",
            f"{branch}_loaded_flip",
            f"{branch}_loaded_release",
        )
    state_energies = {state: 0.0 for state in states}
    # The composite loaded state includes an unmeasured bicarbonate binding
    # standard state.  Fixing it to the 25 mM assay bath prevents the B^b
    # mass-action factor from being silently interpreted as a unit-mM Kd.
    # This is a canonical lumped gauge, varied later through edge barriers.
    for branch_states in transported_states.values():
        for state in branch_states:
            state_energies[state] += hco3_count * math.log(25.0)
    return CycleDefinition(
        model_id=model_id,
        family=family,
        states=tuple(states),
        edges=tuple(edges),
        state_energies=MappingProxyType(state_energies),
        branch_current_edges=MappingProxyType(current_edges),
        branch_cycle_edge_ids=MappingProxyType(cycle_edges),
        intracellular_source_vectors=MappingProxyType(sources),
        evidence_status="source-compatible candidate; microscopic rates not identified",
        free_parameter_names=tuple(
            ["capacity", "Na/K lumped attempt-rate ratio"]
        ),
        source_constraints=(
            "direct Na transport",
            "direct K transport" if len(branches) > 1 else f"pure-{branches[0]} assay submodel",
            "electroneutrality",
            "reversibility",
            "local detailed balance",
        ),
        new_assumptions=(
            f"composite {cation_count}-cation-plus-{hco3_count}-bicarbonate loaded state",
            "equal unmeasured state standard energies",
            "unit attempt rates before lumped cation scaling",
        ),
        transported_cation_bound_states=MappingProxyType(transported_states),
        transported_cation_counts=MappingProxyType(
            {branch: cation_count for branch in branches}
        ),
        mutant_sensitive_edges=MappingProxyType(mutant_edges),
    )


def _shared_carbonate_definition(
    model_id: str,
    family: str,
    *,
    branches: Sequence[str] = ("na", "k"),
) -> CycleDefinition:
    """SR4a: one transported cation and one explicit carbonate."""

    states = ["Eo", "EoCl", "EiCl", "Ei"]
    edges = [
        _edge("common_cl_bind", "Eo", "EoCl", cl_e=+1),
        _edge("common_cl_flip", "EoCl", "EiCl"),
        _edge("common_cl_release", "EiCl", "Ei", cl_i=-1),
    ]
    sources: dict[str, Mapping[str, float]] = {}
    current_edges: dict[str, tuple[str, ...]] = {}
    cycle_edges: dict[str, tuple[str, ...]] = {}
    transported_states: dict[str, tuple[str, ...]] = {}
    mutant_edges: dict[str, tuple[str, ...]] = {}
    for branch in branches:
        state_i = f"Ei{branch.title()}CO3"
        state_o = f"Eo{branch.title()}CO3"
        states.extend((state_i, state_o))
        edges.extend(
            (
                _edge(
                    f"{branch}_loaded_bind",
                    "Ei",
                    state_i,
                    **{f"{branch}_i": +1, "co3_i": +1},
                ),
                _edge(f"{branch}_loaded_flip", state_i, state_o),
                _edge(
                    f"{branch}_loaded_release",
                    state_o,
                    "Eo",
                    **{f"{branch}_e": -1, "co3_e": -1},
                ),
            )
        )
        sources[branch] = MappingProxyType(
            {
                "na_i": -1.0 if branch == "na" else 0.0,
                "k_i": -1.0 if branch == "k" else 0.0,
                "cl_i": +1.0,
                "co3_i": -1.0,
            }
        )
        current_edges[branch] = (f"{branch}_loaded_flip",)
        transported_states[branch] = (state_i, state_o)
        mutant_edges[branch] = (f"{branch}_loaded_flip",)
        cycle_edges[branch] = (
            "common_cl_bind",
            "common_cl_flip",
            "common_cl_release",
            f"{branch}_loaded_bind",
            f"{branch}_loaded_flip",
            f"{branch}_loaded_release",
        )
    state_energies = {state: 0.0 for state in states}
    # Reference carbonate at 25 mM bicarbonate and pH 7.4 with pKa2 10.33.
    co3_reference_mM = 25.0 * 10.0 ** (7.4 - 10.33)
    for branch_states in transported_states.values():
        for state in branch_states:
            state_energies[state] += math.log(co3_reference_mM)
    return CycleDefinition(
        model_id=model_id,
        family=family,
        states=tuple(states),
        edges=tuple(edges),
        state_energies=MappingProxyType(state_energies),
        branch_current_edges=MappingProxyType(current_edges),
        branch_cycle_edge_ids=MappingProxyType(cycle_edges),
        intracellular_source_vectors=MappingProxyType(sources),
        evidence_status="minimal explicit carbonate alternative; carbon species unmeasured",
        free_parameter_names=("capacity", "Na/K lumped attempt-rate ratio"),
        source_constraints=("electroneutrality", "reversibility", "local detailed balance"),
        new_assumptions=(
            "explicit CO3 rather than a hidden bicarbonate power",
            "fast aqueous carbonate equilibrium requires a declared pKa2 and pH",
            "equal unmeasured state standard energies",
        ),
        transported_cation_bound_states=MappingProxyType(transported_states),
        transported_cation_counts=MappingProxyType({branch: 1 for branch in branches}),
        mutant_sensitive_edges=MappingProxyType(mutant_edges),
    )


def _sr5_catalytic_definition(model_id: str, family: str) -> CycleDefinition:
    """SR5 with explicit same-reservoir catalytic-bound states."""

    base = _shared_bicarbonate_definition(
        model_id, family, stoichiometry=(1, 1, 2)
    )
    states = list(base.states)
    edges = list(base.edges)
    current_edges = dict(base.branch_current_edges)
    catalytic_states: dict[str, tuple[str, ...]] = {}
    mutant_edges = dict(base.mutant_sensitive_edges)
    for branch in ("na", "k"):
        loaded_i = f"Ei{branch.title()}1B2"
        loaded_o = f"Eo{branch.title()}1B2"
        bound_i = f"{loaded_i}R{branch.title()}e"
        bound_o = f"{loaded_o}R{branch.title()}e"
        states.extend((bound_i, bound_o))
        # E16-07 varied the external cation reservoir.  It is used on both
        # sides of this catalytic detour and therefore cancels exactly from
        # net transport.  The reservoir placement is discrete and NEW.
        edges.extend(
            (
                _edge(
                    f"{branch}_catalyst_bind",
                    loaded_i,
                    bound_i,
                    **{f"{branch}_e": +1},
                ),
                _edge(f"{branch}_catalytic_flip", bound_i, bound_o),
                _edge(
                    f"{branch}_catalyst_release",
                    bound_o,
                    loaded_o,
                    **{f"{branch}_e": -1},
                ),
            )
        )
        current_edges[branch] = (
            f"{branch}_loaded_flip",
            f"{branch}_catalytic_flip",
        )
        catalytic_states[branch] = (bound_i, bound_o)
        mutant_edges[branch] = (
            f"{branch}_loaded_flip",
            f"{branch}_catalytic_flip",
        )
    state_energies = {state: float(base.state_energies.get(state, 0.0)) for state in states}
    # A catalytic-bound state contains its underlying transported loaded-state
    # standard energy; the additional cation Kd is added at evaluation time.
    for branch, (bound_i, bound_o) in catalytic_states.items():
        state_energies[bound_i] = base.state_energies[f"Ei{branch.title()}1B2"]
        state_energies[bound_o] = base.state_energies[f"Eo{branch.title()}1B2"]
    return CycleDefinition(
        model_id=model_id,
        family=family,
        states=tuple(states),
        edges=tuple(edges),
        state_energies=MappingProxyType(state_energies),
        branch_current_edges=MappingProxyType(current_edges),
        branch_cycle_edge_ids=base.branch_cycle_edge_ids,
        intracellular_source_vectors=base.intracellular_source_vectors,
        evidence_status="source-compatible cooperative extension; regulatory site is unmeasured",
        free_parameter_names=(
            "capacity",
            "Na/K lumped attempt-rate ratio",
            "Na/K catalytic-state effective dissociation energies",
        ),
        source_constraints=base.source_constraints + (
            "external Na/K dose-response effective EC50 and slope",
        ),
        new_assumptions=base.new_assumptions + (
            "one nontransported external-reservoir cation catalytic occupancy per branch",
        ),
        transported_cation_bound_states=base.transported_cation_bound_states,
        transported_cation_counts=base.transported_cation_counts,
        catalytic_bound_states=MappingProxyType(catalytic_states),
        mutant_sensitive_edges=MappingProxyType(mutant_edges),
    )


def _sr6_layered_definition(
    model_id: str, family: str, *, regulated_placement: str
) -> CycleDefinition:
    """SR6 with explicit unregulated/regulated CTMC layers."""

    core = _sr5_catalytic_definition("SR5_CORE", "SR5 core")
    states = tuple(f"{state}^{layer}" for layer in ("U", "P") for state in core.states)
    edges: list[ReversibleEdge] = []
    for layer in ("U", "P"):
        for edge in core.edges:
            edges.append(
                ReversibleEdge(
                    edge_id=f"{layer}::{edge.edge_id}",
                    source=f"{edge.source}^{layer}",
                    target=f"{edge.target}^{layer}",
                    attempt_rate=edge.attempt_rate,
                    reservoir_stoichiometry=edge.reservoir_stoichiometry,
                )
            )
    edges.append(_edge("pka_switch", "Ei^U", "Ei^P"))
    current_edges = {
        branch: tuple(
            f"{layer}::{edge_id}"
            for layer in ("U", "P")
            for edge_id in core.branch_current_edges[branch]
        )
        for branch in core.branch_current_edges
    }
    cycle_edges = {
        branch: tuple(f"U::{edge_id}" for edge_id in core.branch_cycle_edge_ids[branch])
        for branch in core.branch_cycle_edge_ids
    }
    catalytic_states = {
        branch: tuple(
            f"{state}^{layer}"
            for layer in ("U", "P")
            for state in branch_states
        )
        for branch, branch_states in core.catalytic_bound_states.items()
    }
    transported_states = {
        branch: tuple(
            f"{state}^{layer}"
            for layer in ("U", "P")
            for state in branch_states
        )
        for branch, branch_states in core.transported_cation_bound_states.items()
    }
    mutant_edges = {
        branch: tuple(
            f"{layer}::{edge_id}"
            for layer in ("U", "P")
            for edge_id in edge_ids
        )
        for branch, edge_ids in core.mutant_sensitive_edges.items()
    }
    if regulated_placement == "loaded_flip":
        regulated_edges = tuple(
            f"P::{branch}_{transition}"
            for branch in ("na", "k")
            for transition in ("loaded_flip", "catalytic_flip")
        )
    elif regulated_placement == "common_cl_flip":
        regulated_edges = ("P::common_cl_flip",)
    else:
        raise ValueError("unknown SR6 regulated transition placement")
    layered_energies = {
        f"{state}^{layer}": float(core.state_energies[state])
        for layer in ("U", "P")
        for state in core.states
    }
    return CycleDefinition(
        model_id=model_id,
        family=family,
        states=states,
        edges=tuple(edges),
        state_energies=MappingProxyType(layered_energies),
        branch_current_edges=MappingProxyType(current_edges),
        branch_cycle_edge_ids=MappingProxyType(cycle_edges),
        intracellular_source_vectors=core.intracellular_source_vectors,
        evidence_status="PKA-associated regulation supported; phosphorylation identity and kinetics unmeasured",
        free_parameter_names=core.free_parameter_names
        + ("regulated-layer barrier fold", "steady phosphorylation/dephosphorylation ratio"),
        source_constraints=core.source_constraints
        + ("approximately 25% recombinant forskolin response", "S173-associated activation"),
        new_assumptions=core.new_assumptions
        + (
            "U/P labels are an effective regulatory layer, not proof of direct phosphorylation",
            "one cytosol-accessible Ei switching connection",
            f"{regulated_placement} placement of the regulated barrier",
        ),
        pka_edge_ids=regulated_edges,
        pka_switch_edge_ids=("pka_switch",),
        transported_cation_bound_states=MappingProxyType(transported_states),
        transported_cation_counts=core.transported_cation_counts,
        catalytic_bound_states=MappingProxyType(catalytic_states),
        mutant_sensitive_edges=MappingProxyType(mutant_edges),
    )


def _sequential_definition(
    model_id: str,
    family: str,
    *,
    carbonate: bool,
) -> CycleDefinition:
    branch = "co3" if carbonate else "sequential"
    loaded_in = "EiK2CO3" if carbonate else "EiKHCO3"
    loaded_out = "EoK2CO3" if carbonate else "EoKHCO3"
    states = ("Eo", "EoNaCl", "EiNaCl", "Ei", loaded_in, loaded_out)
    if carbonate:
        return_ligands = {"k_i": +2, "co3_i": +1}
        release_ligands = {"k_e": -2, "co3_e": -1}
        source = {"na_i": +1.0, "k_i": -2.0, "cl_i": +1.0, "co3_i": -1.0}
    else:
        return_ligands = {"k_i": +1, "hco3_i": +1}
        release_ligands = {"k_e": -1, "hco3_e": -1}
        source = {"na_i": +1.0, "k_i": -1.0, "cl_i": +1.0, "hco3_i": -1.0}
    edges = (
        _edge(f"{branch}_nacl_bind", "Eo", "EoNaCl", na_e=+1, cl_e=+1),
        _edge(f"{branch}_nacl_flip", "EoNaCl", "EiNaCl"),
        _edge(f"{branch}_nacl_release", "EiNaCl", "Ei", na_i=-1, cl_i=-1),
        _edge(f"{branch}_return_bind", "Ei", loaded_in, **return_ligands),
        _edge(f"{branch}_return_flip", loaded_in, loaded_out),
        _edge(f"{branch}_return_release", loaded_out, "Eo", **release_ligands),
    )
    return CycleDefinition(
        model_id=model_id,
        family=family,
        states=states,
        edges=edges,
        state_energies=MappingProxyType({state: 0.0 for state in states}),
        branch_current_edges=MappingProxyType({branch: (f"{branch}_return_flip",)}),
        branch_cycle_edge_ids=MappingProxyType(
            {branch: tuple(edge.edge_id for edge in edges)}
        ),
        intracellular_source_vectors=MappingProxyType(
            {branch: MappingProxyType(source)}
        ),
        evidence_status=(
            "2025 discussion hypothesis with explicit carbonate chemistry"
            if carbonate
            else "2025 discussion hypothesis; not an established sole AE4 cycle"
        ),
        free_parameter_names=("capacity", "lumped sequential attempt-rate scale"),
        source_constraints=("electroneutrality", "reversibility", "local detailed balance"),
        new_assumptions=(
            "fixed six-state ordering proposed from a discussion-level mechanism",
            "equal unmeasured state standard energies",
        ),
    )


def candidate_definitions() -> Mapping[str, CycleDefinition]:
    """Return the frozen SR1--SR6 graph registry."""

    definitions = {
        "SR1_NA_112": _shared_bicarbonate_definition(
            "SR1_NA_112", "SR1 minimal pure-Na 1:1:2", branches=("na",)
        ),
        "SR1_K_112": _shared_bicarbonate_definition(
            "SR1_K_112", "SR1 minimal pure-K 1:1:2", branches=("k",)
        ),
        "SR2_SHARED_112": _shared_bicarbonate_definition(
            "SR2_SHARED_112", "SR2 shared-state Na/K 1:1:2"
        ),
        "SR2_SHARED_123": _shared_bicarbonate_definition(
            "SR2_SHARED_123",
            "SR2 shared-state Na/K 1:2:3",
            stoichiometry=(1, 2, 3),
        ),
        "SR2_SHARED_213": _shared_bicarbonate_definition(
            "SR2_SHARED_213",
            "SR2 shared-state Na/K 2:1:3",
            stoichiometry=(2, 1, 3),
        ),
        "SR3_SEQUENTIAL_1111": _sequential_definition(
            "SR3_SEQUENTIAL_1111", "SR3 asymmetric Cl:Na:HCO3:K 1:1:1:1", carbonate=False
        ),
        "SR4A_CARBONATE_111": _shared_carbonate_definition(
            "SR4A_CARBONATE_111", "SR4a shared-state Cl:cation:CO3 1:1:1"
        ),
        "SR4B_CARBONATE_1112": _sequential_definition(
            "SR4B_CARBONATE_1112", "SR4b asymmetric Cl:Na:CO3:K 1:1:1:2", carbonate=True
        ),
        "SR5_COORDINATION_112": _sr5_catalytic_definition(
            "SR5_COORDINATION_112",
            "SR5 shared-state Na/K with explicit catalytic cation occupancy",
        ),
        "SR6_PKA_EDGE_112": _sr6_layered_definition(
            "SR6_PKA_EDGE_112",
            "SR6 explicit U/P layers on SR5 with regulated loaded transitions",
            regulated_placement="loaded_flip",
        ),
        "SR6_PKA_COMMON_FLIP_112": _sr6_layered_definition(
            "SR6_PKA_COMMON_FLIP_112",
            "SR6 explicit U/P layers on SR5 with regulated common Cl transition",
            regulated_placement="common_cl_flip",
        ),
    }
    return MappingProxyType(definitions)


CANDIDATE_DEFINITIONS = candidate_definitions()


MUTANT_BRANCH_SCALES: Mapping[str, tuple[float, float]] = MappingProxyType(
    {
        # These are qualitative/lumped activity encodings of the primary
        # hierarchy, not microscopic parameter estimates.  They are retained
        # as explicit assumptions for sensitivity checks.
        "WT": (1.0, 1.0),
        "T448I": (0.5, 0.5),
        "T756A": (0.7, 0.7),
        # A tiny rather than exact-zero Na scale keeps the finite CTMC
        # irreducible while representing activity below assay detection.
        "T448I_T756A": (1e-12, 0.5),
    }
)


def parameters_for_mutant(
    parameters: CycleParameters, mutant: str
) -> CycleParameters:
    """Apply the source-constrained Na-versus-K mutant hierarchy."""

    try:
        na_scale, k_scale = MUTANT_BRANCH_SCALES[mutant]
    except KeyError as exc:
        raise KeyError(f"unknown mutant {mutant}") from exc
    return replace(
        parameters,
        mutant_na_scale=na_scale,
        mutant_k_scale=k_scale,
    )


def branch_affinity_from_sources(
    definition: CycleDefinition,
    branch: str,
    concentrations_mM: ReservoirActivities,
    reference_concentration_mM: float = 1.0,
) -> float:
    """Compute a fundamental-cycle affinity independently of the CTMC solve."""

    concentrations_mM.validate(definition.required_reservoir_species)
    activities = {
        name: value / reference_concentration_mM
        for name, value in concentrations_mM.as_mapping().items()
        if math.isfinite(value)
    }
    return float(
        sum(
            coefficient * math.log(activities[name])
            for edge in _branch_cycle_edges(definition, branch)
            for name, coefficient in edge.reservoir_stoichiometry.items()
        )
    )


def exact_charge_per_branch(definition: CycleDefinition) -> Mapping[str, float]:
    """Return charge equivalents moved into the cell per completed branch."""

    values = {}
    for branch, source in definition.intracellular_source_vectors.items():
        values[branch] = float(
            sum(ION_CHARGES[name] * coefficient for name, coefficient in source.items())
        )
    return MappingProxyType(values)
