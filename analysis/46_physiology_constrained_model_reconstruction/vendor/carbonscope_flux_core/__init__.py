"""Local Task 46 CBM solver core.

No runtime access to CarbonScope is required or permitted by the Task 46 workflow.
"""

from .schema import (
    FluxMetabolite,
    FluxReaction,
    FluxModel,
    LinearObjective,
    ObjectiveTerm,
    StoichiometricTerm,
)
from .solver import (
    CompiledFluxLP,
    PreparedFluxRegion,
    assert_fva_equivalent,
    compile_flux_lp,
    prepare_highs_flux_region,
    run_highs_fba,
    run_highs_fva_reference,
    run_highs_vffva,
    run_prepared_highs_vffva,
)

__all__ = [
    "FluxMetabolite",
    "FluxReaction",
    "FluxModel",
    "LinearObjective",
    "ObjectiveTerm",
    "StoichiometricTerm",
    "CompiledFluxLP",
    "PreparedFluxRegion",
    "assert_fva_equivalent",
    "compile_flux_lp",
    "prepare_highs_flux_region",
    "run_highs_fba",
    "run_highs_fva_reference",
    "run_highs_vffva",
    "run_prepared_highs_vffva",
]
