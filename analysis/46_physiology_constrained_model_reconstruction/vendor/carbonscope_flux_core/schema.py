"""Minimal immutable flux-model schema for Task 46 CBM screening.

Adapted from esig626/CarbonScope commit
`a11e25f176a2cc70de22c3f36a1c0e623248393a`.
"""

from __future__ import annotations

from dataclasses import dataclass

Number = int | float


@dataclass(frozen=True, slots=True)
class FluxMetabolite:
    metabolite_id: str
    steady_state_balanced: bool = True


@dataclass(frozen=True, slots=True)
class StoichiometricTerm:
    metabolite_id: str
    coefficient: Number


@dataclass(frozen=True, slots=True)
class FluxReaction:
    reaction_id: str
    stoichiometric_terms: tuple[StoichiometricTerm, ...]
    lower_bound: Number
    upper_bound: Number


@dataclass(frozen=True, slots=True)
class ObjectiveTerm:
    reaction_id: str
    coefficient: Number


@dataclass(frozen=True, slots=True)
class LinearObjective:
    direction: str
    terms: tuple[ObjectiveTerm, ...]


@dataclass(frozen=True, slots=True)
class FluxModel:
    metabolites: tuple[FluxMetabolite, ...]
    reactions: tuple[FluxReaction, ...]
    objective: LinearObjective
