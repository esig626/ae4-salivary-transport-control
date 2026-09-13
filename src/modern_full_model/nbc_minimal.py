"""Minimal Na+/HCO3- cotransporter for the missing-alkalinity test.

This is a deliberately small structural intervention.  It is NOT an isoform
assignment and it is NOT a claim that mouse SMG acinar cells express a
particular NBC at the required abundance.  It asks one narrower question:
can a basolateral pathway that imports Na+ and HCO3- without exporting Cl-
remove the alkalinity bottleneck that otherwise forces reverse AE2 to cancel
AE4-mediated chloride loading?

Transport law
-------------
The surrogate is an electroneutral 1 Na : 1 HCO3 cotransporter.  Positive
flux is bath -> cell.  The ideal-solution affinity is

    A = log((Na_o * HCO3_o) / (Na_i * HCO3_i))

and the bounded reversible flux is

    J_NBC = G_NBC * tanh(A / w).

Cell source terms are therefore

    dNa/dt  += J_NBC
    dTIC/dt += J_NBC
    dTA/dt  += J_NBC
    dCl/dt  += 0.

The source is exactly electroneutral because +J_NBC Na equivalents are
balanced by +J_NBC alkalinity equivalents.

The default capacity is a DERIVED FEASIBILITY SCALE, not a measured NBC
capacity.  It comes from the Task 31 R09 flux scale plus the algebraic 70/30
chloride-loading requirement, source-fixed 2.3x stimulated NHE1, and the
existing AE2 positive capacity.  The resulting required NBC flux is
0.1964897261082554 fmol/s.  At the Task 31 R09 chemical state the simple
reversible law has tanh(A/2)=0.9674029262, so the corresponding capacity is
0.20311053520552247 fmol/s.

No optimiser is implied by this module.  The reference capacity is one fixed
structural test value and should not be swept automatically.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from .model import Genotype, ModelEvaluation, WT


REFERENCE_REQUIRED_NBC_FLUX_FMOL_S = 0.1964897261082554
REFERENCE_MINIMUM_NBC_CAPACITY_FMOL_S = 0.20311053520552247
REFERENCE_AFFINITY_LOG = 4.100247118100945


@dataclass(frozen=True)
class MinimalNbcParameters:
    """One-capacity electroneutral Na/HCO3 structural surrogate."""

    capacity_fmol_s: float = REFERENCE_MINIMUM_NBC_CAPACITY_FMOL_S
    log_width: float = 2.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.capacity_fmol_s) or self.capacity_fmol_s < 0.0:
            raise ValueError("NBC capacity must be finite and nonnegative")
        if not math.isfinite(self.log_width) or self.log_width <= 0.0:
            raise ValueError("NBC log width must be finite and positive")


@dataclass(frozen=True)
class MinimalNbcFlux:
    """Basolateral NBC flux and exact conserved cell sources."""

    inward_fmol_s: float
    affinity_log: float
    na_cell_fmol_s: float
    cl_cell_fmol_s: float
    tic_cell_fmol_s: float
    alkalinity_cell_fmol_s: float
    transported_charge_equivalents_fmol_s: float
    charge_source_residual_fmol_s: float



def evaluate_minimal_nbc(
    *,
    na_i_mM: float,
    hco3_i_mM: float,
    na_o_mM: float,
    hco3_o_mM: float,
    parameters: MinimalNbcParameters | None = None,
    activity_scale: float = 1.0,
) -> MinimalNbcFlux:
    """Evaluate the reversible 1 Na : 1 HCO3 cotransporter.

    Positive flux is bath -> cell.  Because the transported pair is neutral,
    there is no membrane-current term and no direct chloride source.
    """

    p = parameters or MinimalNbcParameters()
    for value in (na_i_mM, hco3_i_mM, na_o_mM, hco3_o_mM):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("NBC concentrations must be finite and positive")
    if not math.isfinite(activity_scale) or activity_scale < 0.0:
        raise ValueError("NBC activity scale must be finite and nonnegative")

    affinity = math.log(
        (na_o_mM * hco3_o_mM) / (na_i_mM * hco3_i_mM)
    )
    flux = (
        p.capacity_fmol_s
        * activity_scale
        * math.tanh(affinity / p.log_width)
    )
    na_source = flux
    tic_source = flux
    alkalinity_source = flux
    charge_residual = na_source - alkalinity_source

    return MinimalNbcFlux(
        inward_fmol_s=flux,
        affinity_log=affinity,
        na_cell_fmol_s=na_source,
        cl_cell_fmol_s=0.0,
        tic_cell_fmol_s=tic_source,
        alkalinity_cell_fmol_s=alkalinity_source,
        transported_charge_equivalents_fmol_s=0.0,
        charge_source_residual_fmol_s=charge_residual,
    )


class MinimalNbcModel:
    """Conservation-preserving wrapper that adds only the NBC source terms.

    The wrapped model may be the core ModernFullModel or another compatible
    conservation-preserving wrapper.  No existing transporter parameter,
    electrical closure, water law, AE4 law, or regulatory state is replaced.
    """

    def __init__(
        self,
        base_model: Any,
        parameters: MinimalNbcParameters | None = None,
        *,
        activity_scale: float = 1.0,
    ) -> None:
        self.base_model = base_model
        self.nbc_parameters = parameters or MinimalNbcParameters()
        if not math.isfinite(activity_scale) or activity_scale < 0.0:
            raise ValueError("NBC activity scale must be finite and nonnegative")
        self.nbc_activity_scale = float(activity_scale)

    @property
    def parameters(self) -> Any:
        return self.base_model.parameters

    @property
    def stimulus(self) -> Any:
        return self.base_model.stimulus

    @property
    def regulatory_model(self) -> Any:
        return getattr(self.base_model, "regulatory_model", None)

    @property
    def ae4_parameters(self) -> Any:
        return getattr(self.base_model, "ae4_parameters", None)

    @property
    def ae4_evaluator(self) -> Any:
        return getattr(self.base_model, "ae4_evaluator", None)

    @property
    def layout(self) -> Any:
        return self.base_model.layout

    @property
    def state_names(self) -> tuple[str, ...]:
        return self.base_model.state_names

    def initial_state(self) -> NDArray[np.float64]:
        return self.base_model.initial_state()

    def _evaluate_nbc(self, baseline: ModelEvaluation) -> MinimalNbcFlux:
        obs = baseline.diagnostics.observables
        ci = obs.cell_concentrations_mM
        bath = self.parameters.bath
        return evaluate_minimal_nbc(
            na_i_mM=float(ci["na"]),
            hco3_i_mM=float(ci["hco3"]),
            na_o_mM=float(bath.na_mM),
            hco3_o_mM=float(obs.bath_acid_base.hco3_mM),
            parameters=self.nbc_parameters,
            activity_scale=self.nbc_activity_scale,
        )

    def evaluate(
        self,
        time_s: float,
        vector: ArrayLike,
        *,
        genotype: Genotype = WT,
    ) -> ModelEvaluation:
        baseline = self.base_model.evaluate(time_s, vector, genotype=genotype)
        nbc = self._evaluate_nbc(baseline)
        raw = np.asarray(baseline.rhs, dtype=float).copy()

        # Core state order: Na_i, K_i, Cl_i, TIC_i, TA_i, V_i, then lumen.
        raw[0] += nbc.na_cell_fmol_s
        raw[3] += nbc.tic_cell_fmol_s
        raw[4] += nbc.alkalinity_cell_fmol_s

        conservation = dict(baseline.diagnostics.conservation_residuals)
        conservation["cell_bulk_charge_rate_fmol_s"] = float(
            raw[0] + raw[1] - raw[2] - raw[4]
        )
        conservation["minimal_nbc_charge_fmol_s"] = float(
            nbc.charge_source_residual_fmol_s
        )
        conservation["minimal_nbc_carbon_source_fmol_s"] = float(
            nbc.tic_cell_fmol_s
        )

        regulatory = dict(baseline.diagnostics.regulatory)
        regulatory.update(
            {
                "minimal_nbc_model": "ELECTRONEUTRAL_1NA_1HCO3",
                "minimal_nbc_inward_fmol_s": float(nbc.inward_fmol_s),
                "minimal_nbc_affinity_log": float(nbc.affinity_log),
                "minimal_nbc_capacity_fmol_s": float(
                    self.nbc_parameters.capacity_fmol_s
                ),
                "minimal_nbc_activity_scale": self.nbc_activity_scale,
                "minimal_nbc_evidence_status": (
                    "STRUCTURAL_SURROGATE; CAPACITY_DERIVED_FROM_FEASIBILITY_NOT_MEASURED"
                ),
            }
        )

        diagnostics = replace(
            baseline.diagnostics,
            regulatory=regulatory,
            conservation_residuals=conservation,
        )
        return replace(baseline, rhs=raw, diagnostics=diagnostics)

    def rhs(
        self,
        time_s: float,
        vector: ArrayLike,
        *,
        genotype: Genotype = WT,
    ) -> NDArray[np.float64]:
        return self.evaluate(time_s, vector, genotype=genotype).rhs

    def solve_dynamics(
        self,
        time_span_s: tuple[float, float],
        *,
        initial_state: ArrayLike | None = None,
        genotype: Genotype = WT,
        method: str = "Radau",
        rtol: float = 1.0e-7,
        atol: float | ArrayLike = 1.0e-10,
        t_eval: ArrayLike | None = None,
        max_step_s: float = math.inf,
    ) -> Any:
        if method not in {"Radau", "BDF"}:
            raise ValueError("production dynamics must use method='Radau' or 'BDF'")
        y0 = (
            self.initial_state()
            if initial_state is None
            else self.layout.validate(initial_state)
        )
        return solve_ivp(
            lambda t, y: self.rhs(t, y, genotype=genotype),
            time_span_s,
            y0,
            method=method,
            rtol=rtol,
            atol=atol,
            t_eval=None if t_eval is None else np.asarray(t_eval, dtype=float),
            max_step=max_step_s,
        )


__all__ = (
    "REFERENCE_AFFINITY_LOG",
    "REFERENCE_MINIMUM_NBC_CAPACITY_FMOL_S",
    "REFERENCE_REQUIRED_NBC_FLUX_FMOL_S",
    "MinimalNbcFlux",
    "MinimalNbcModel",
    "MinimalNbcParameters",
    "evaluate_minimal_nbc",
)
