"""Minimal validation for the local Task 46 flux model."""

from __future__ import annotations

import math

from .schema import FluxModel


class CanonicalModelError(ValueError):
    pass


def validate_flux_model(model: FluxModel) -> None:
    if not isinstance(model, FluxModel):
        raise CanonicalModelError("flux model must be a FluxModel")
    if not isinstance(model.metabolites, tuple) or not isinstance(model.reactions, tuple):
        raise CanonicalModelError("metabolites and reactions must be tuples")
    metabolite_ids = [m.metabolite_id for m in model.metabolites]
    reaction_ids = [r.reaction_id for r in model.reactions]
    if any(not isinstance(x, str) or not x for x in metabolite_ids + reaction_ids):
        raise CanonicalModelError("all IDs must be nonempty strings")
    if len(set(metabolite_ids)) != len(metabolite_ids):
        raise CanonicalModelError("duplicate metabolite ID")
    if len(set(reaction_ids)) != len(reaction_ids):
        raise CanonicalModelError("duplicate reaction ID")
    metabolite_set = set(metabolite_ids)
    for metabolite in model.metabolites:
        if type(metabolite.steady_state_balanced) is not bool:
            raise CanonicalModelError("steady_state_balanced must be bool")
    for reaction in model.reactions:
        if not isinstance(reaction.stoichiometric_terms, tuple):
            raise CanonicalModelError("stoichiometric_terms must be a tuple")
        if isinstance(reaction.lower_bound, bool) or isinstance(reaction.upper_bound, bool):
            raise CanonicalModelError("reaction bounds must be numeric")
        lower = float(reaction.lower_bound)
        upper = float(reaction.upper_bound)
        if not math.isfinite(lower) or not math.isfinite(upper) or lower > upper:
            raise CanonicalModelError(f"invalid bounds for {reaction.reaction_id}")
        for term in reaction.stoichiometric_terms:
            if term.metabolite_id not in metabolite_set:
                raise CanonicalModelError(
                    f"unknown metabolite {term.metabolite_id!r} in {reaction.reaction_id!r}"
                )
            if isinstance(term.coefficient, bool) or not math.isfinite(float(term.coefficient)):
                raise CanonicalModelError("stoichiometric coefficient must be finite")
    if model.objective.direction not in {"maximise", "minimise"}:
        raise CanonicalModelError("objective direction must be 'maximise' or 'minimise'")
    if not isinstance(model.objective.terms, tuple) or not model.objective.terms:
        raise CanonicalModelError("objective must have at least one term")
    reaction_set = set(reaction_ids)
    for term in model.objective.terms:
        if term.reaction_id not in reaction_set:
            raise CanonicalModelError(f"objective references unknown reaction {term.reaction_id!r}")
        if isinstance(term.coefficient, bool) or not math.isfinite(float(term.coefficient)):
            raise CanonicalModelError("objective coefficient must be finite")
