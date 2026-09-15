"""Compact result records for Task 46 FBA/FVA."""

from __future__ import annotations

from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True, slots=True)
class PrimalDiagnostics:
    max_lower_bound_violation: float
    max_upper_bound_violation: float
    max_mass_balance_residual: float
    objective_recalculation_error: float


@dataclass(frozen=True)
class FBAResult:
    objective_value: float
    status: str
    objective_direction: str
    fluxes: pd.Series
    diagnostics: PrimalDiagnostics | None = None


@dataclass(frozen=True)
class FVAResult:
    ranges: pd.DataFrame
    fraction_of_optimum: float | None
    objective_value: float
    objective_direction: str
