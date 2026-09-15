"""Minimal FBA/reference-FVA/VFFVA equivalence smoke test."""

from __future__ import annotations

from . import (
    FluxMetabolite,
    FluxModel,
    FluxReaction,
    LinearObjective,
    ObjectiveTerm,
    StoichiometricTerm,
    assert_fva_equivalent,
    run_highs_fba,
    run_highs_fva_reference,
    run_highs_vffva,
)


def build_model() -> FluxModel:
    return FluxModel(
        metabolites=(FluxMetabolite("A", True),),
        reactions=(
            FluxReaction("IN", (StoichiometricTerm("A", 1.0),), 0.0, 10.0),
            FluxReaction("OUT", (StoichiometricTerm("A", -1.0),), 0.0, 10.0),
        ),
        objective=LinearObjective("maximise", (ObjectiveTerm("OUT", 1.0),)),
    )


def main() -> None:
    model = build_model()
    fba = run_highs_fba(model)
    if abs(fba.objective_value - 10.0) > 1e-8:
        raise SystemExit(f"unexpected FBA optimum: {fba.objective_value}")
    reference = run_highs_fva_reference(model, 1.0)
    fast = run_highs_vffva(model, 1.0, workers=2, chunk_size=1)
    assert_fva_equivalent(reference, fast)
    expected = {"IN": (10.0, 10.0), "OUT": (10.0, 10.0)}
    for reaction, (low, high) in expected.items():
        row = fast.ranges.loc[reaction]
        if abs(float(row.minimum) - low) > 1e-8 or abs(float(row.maximum) - high) > 1e-8:
            raise SystemExit(f"unexpected FVA range for {reaction}: {tuple(row)}")
    print("AE4 local CBM core smoke test passed")


if __name__ == "__main__":
    main()
