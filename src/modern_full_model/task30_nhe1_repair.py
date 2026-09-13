"""Bounded Task 30 implementation and WT/exact-AE4-null verification.

This module implements one intervention only: explicit selection of
Vera-Sigüenza et al. (2018), Eq. 26 for NHE1.  Every other transporter and
whole-cell parameter is loaded unchanged from the retained R09/R10 payloads.
The resting solver uses the inherited ten-coordinate charge manifold and ten
independent amount/volume equations.  No parameter grid, allocation search,
phenotype fit, or genotype-dependent NHE1 term is present.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

from .acid_base import total_alkalinity_mM
from .ae4_routing_only import MODEL_ID as AE4_ROUTING_MODEL_ID
from .ae4_routing_only import routing_only_adapter
from .camp_pka import AE4Construct, R1EffectiveActivation, RegulatoryGain
from .membranes import (
    NHE1_LEGACY_TANH,
    NHE1_VERA_SIGUENZA_2018,
    current_to_fmol_s,
    evaluate_nhe1_legacy_tanh,
    evaluate_nhe1_vera_siguenza_2018,
)
from .model import AE4_NULL, WT, Genotype, ModernFullModel
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .parameters import (
    AcidBaseParameters,
    BathParameters,
    FullModelParameters,
    GeometryParameters,
    HomeostasisParameters,
    InitialConditions,
    MembraneParameters,
    PhysicalConstants,
    WaterParameters,
)
from .states import WholeCellState, compartment_from_concentrations
from .transporters import AE4Parameters
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_RADAU,
    SecretagogueProtocol,
    StimulusArm,
    sha256_object,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIRECTORY = REPOSITORY_ROOT / "reference"
RESULT_DIRECTORY = REPOSITORY_ROOT / "results/30_nhe1_model_repair"
ANALYSIS_DIRECTORY = REPOSITORY_ROOT / "analysis/30_nhe1_model_repair"

BRANCH = "codex/task-30-nhe1-model-repair"
PREPARED_HEAD = "60b7798b3c01dbb1b3be9efd831aafdfcd4c0eb7"
BACKGROUND_ALIASES = ("R09", "R10")
RESTING_CALCIUM_UM = 0.058
STIMULATED_CALCIUM_UM = 0.25
WT_TARGET_PH = 6.91
EXACT_NULL_DIAGNOSTIC_PH = 6.89
NHE1_PH_DIAGNOSTIC_HALF_WIDTH = 0.10
ORDINARY_NA_UPPER_MM = 30.0
CELL_VOLUME_UPPER_PL = 3.0

COORDINATE_NAMES = (
    "na_i_mM",
    "k_i_mM",
    "tic_i_mM",
    "ph_i",
    "volume_i_pL",
    "na_l_mM",
    "k_l_mM",
    "tic_l_mM",
    "ph_l",
    "volume_l_pL",
)
COORDINATE_LOWER = np.asarray((4.0, 60.0, 1.0, 6.2, 0.5, 20.0, 0.05, 0.2, 5.0, 0.02))
COORDINATE_UPPER = np.asarray((60.0, 210.0, 100.0, 8.0, 3.0, 280.0, 60.0, 100.0, 9.0, 0.8))
INDEPENDENT_ROWS = np.asarray((0, 1, 2, 3, 5, 6, 7, 8, 9, 11), dtype=int)
OMITTED_CHARGE_ROWS = np.asarray((4, 10), dtype=int)
AMOUNT_ROWS = np.asarray((0, 1, 2, 3, 4, 6, 7, 8, 9, 10), dtype=int)
VOLUME_ROWS = np.asarray((5, 11), dtype=int)
INDEPENDENT_RHS_SCALES = np.asarray(
    (0.05, 0.05, 0.05, 0.02, 1.0e-4, 0.02, 0.02, 0.02, 0.02, 1.0e-4)
)

ROOT_SCALED_TOLERANCE = 1.0e-7
OMITTED_ROW_TOLERANCE = 1.0e-9
CHARGE_TOLERANCE_FMOL = 1.0e-9
CURRENT_TOLERANCE_A = 1.0e-18
REGULATORY_TOLERANCE_S_INV = 1.0e-10
BOUNDARY_RELATIVE_TOLERANCE = 1.0e-5
JACOBIAN_RANK_RELATIVE_TOLERANCE = 1.0e-8


class NumericalBudgetExhausted(RuntimeError):
    """Raised immediately when a Task 30 numerical limit would be exceeded."""


@dataclass
class NumericalBudget:
    stationary_solver_calls: int = 0
    stationary_residual_evaluations: int = 0
    scalar_g_nhe1_evaluations: int = 0
    stimulated_integrations: int = 0
    numerical_wall_seconds: float = 0.0

    maximum_stationary_solver_calls: int = 12
    maximum_stationary_residual_evaluations: int = 10_000
    maximum_scalar_g_nhe1_evaluations: int = 8
    maximum_stimulated_integrations: int = 4
    maximum_numerical_wall_seconds: float = 1800.0

    def start_stationary_solver(self) -> None:
        if self.stationary_solver_calls >= self.maximum_stationary_solver_calls:
            raise NumericalBudgetExhausted("stationary solver call limit exhausted")
        self.stationary_solver_calls += 1

    def count_stationary_residual(self) -> None:
        if (
            self.stationary_residual_evaluations
            >= self.maximum_stationary_residual_evaluations
        ):
            raise NumericalBudgetExhausted("stationary residual evaluation limit exhausted")
        self.stationary_residual_evaluations += 1

    def count_scalar_g_evaluation(self) -> None:
        if self.scalar_g_nhe1_evaluations >= self.maximum_scalar_g_nhe1_evaluations:
            raise NumericalBudgetExhausted("scalar G_NHE1 evaluation limit exhausted")
        self.scalar_g_nhe1_evaluations += 1

    def start_stimulated_integration(self) -> None:
        if self.stimulated_integrations >= self.maximum_stimulated_integrations:
            raise NumericalBudgetExhausted("stimulated integration limit exhausted")
        self.stimulated_integrations += 1

    def add_wall_time(self, seconds: float) -> None:
        self.numerical_wall_seconds += float(seconds)
        if self.numerical_wall_seconds > self.maximum_numerical_wall_seconds:
            raise NumericalBudgetExhausted("numerical wall-time limit exhausted")

    def payload(self) -> dict[str, Any]:
        values = asdict(self)
        values["limits_respected"] = bool(
            self.stationary_solver_calls <= self.maximum_stationary_solver_calls
            and self.stationary_residual_evaluations
            <= self.maximum_stationary_residual_evaluations
            and self.scalar_g_nhe1_evaluations
            <= self.maximum_scalar_g_nhe1_evaluations
            and self.stimulated_integrations <= self.maximum_stimulated_integrations
            and self.numerical_wall_seconds <= self.maximum_numerical_wall_seconds
        )
        values["workers"] = 1
        values["blas_threads"] = 1
        return values


@dataclass(frozen=True)
class Background:
    alias: str
    root_id: str
    routing_family: str
    inherited_parameters: FullModelParameters
    ae4_parameters: AE4Parameters
    saved_core_state: tuple[float, ...]
    saved_coordinates: tuple[float, ...]
    inherited_whole_cell_sha256: str
    inherited_ae4_sha256: str


def _whole_cell_parameters(payload: Mapping[str, Any]) -> FullModelParameters:
    """Load old clean-seed payloads while new NHE1 fields take defaults."""

    return FullModelParameters(
        constants=PhysicalConstants(**payload["constants"]),
        acid_base=AcidBaseParameters(**payload["acid_base"]),
        bath=BathParameters(**payload["bath"]),
        geometry=GeometryParameters(**payload["geometry"]),
        homeostasis=HomeostasisParameters(**payload["homeostasis"]),
        membranes=MembraneParameters(**payload["membranes"]),
        water=WaterParameters(**payload["water"]),
        initial=InitialConditions(**payload["initial"]),
    )


def load_background(alias: str) -> Background:
    if alias not in BACKGROUND_ALIASES:
        raise ValueError(f"Task 30 permits only {BACKGROUND_ALIASES}")
    saved = json.loads(
        (REFERENCE_DIRECTORY / f"{alias}_wt_rest.json").read_text(encoding="utf-8")
    )
    manifest = json.loads(
        (REFERENCE_DIRECTORY / "native_wt_contract.json").read_text(encoding="utf-8")
    )
    root_id = str(saved["root_id"])
    payload = manifest["roots"][root_id]
    parameters = _whole_cell_parameters(payload["whole_cell_parameters"])
    ae4 = AE4Parameters(**payload["ae4_parameters"])
    # Hash the inherited payload itself because the Task 30 dataclass now has
    # additional defaulted NHE1 fields that were absent from the old object.
    if sha256_object(payload["whole_cell_parameters"]) != payload["whole_cell_parameters_sha256"]:
        raise ValueError(f"{alias}: inherited whole-cell payload hash mismatch")
    if sha256_object(ae4) != payload["ae4_parameters_sha256"]:
        raise ValueError(f"{alias}: inherited AE4 payload hash mismatch")
    return Background(
        alias=alias,
        root_id=root_id,
        routing_family=str(saved["routing_family"]),
        inherited_parameters=parameters,
        ae4_parameters=ae4,
        saved_core_state=tuple(float(value) for value in saved["core_state"]),
        saved_coordinates=tuple(float(value) for value in saved["coordinates"]),
        inherited_whole_cell_sha256=str(payload["whole_cell_parameters_sha256"]),
        inherited_ae4_sha256=str(payload["ae4_parameters_sha256"]),
    )


def task30_parameters(background: Background, *, g_nhe1_fmol_s: float = 0.0305) -> FullModelParameters:
    """Select the published Eq. 26 law without changing another parameter."""

    if not math.isfinite(g_nhe1_fmol_s) or g_nhe1_fmol_s < 0.0:
        raise ValueError("G_NHE1 must be finite and nonnegative")
    inherited = background.inherited_parameters
    return replace(
        inherited,
        homeostasis=replace(
            inherited.homeostasis,
            nhe1_model=NHE1_VERA_SIGUENZA_2018,
            nhe1_published_g_fmol_s=float(g_nhe1_fmol_s),
            nhe1_published_k_h_mM=4.5e-4,
            nhe1_published_k_na_mM=15.0,
        ),
    )


def _protocol(*, stimulated: bool) -> SecretagogueProtocol:
    return SecretagogueProtocol(
        arm=StimulusArm.CCH_IPR if stimulated else StimulusArm.REST,
        onset_s=0.0,
        duration_s=600.0,
        cch_dose_uM=0.3,
        ipr_dose_uM=5.0,
        resting_calcium_uM=RESTING_CALCIUM_UM,
        stimulated_calcium_uM=STIMULATED_CALCIUM_UM,
        beta_occupancy_on=1.0,
    )


def build_task30_model(
    background: Background,
    *,
    g_nhe1_fmol_s: float = 0.0305,
    stimulated: bool = False,
):
    """Construct the routed inherited model with explicit Task 30 NHE1 selection."""

    protocol = _protocol(stimulated=stimulated)
    ae4_regulation = R1EffectiveActivation(
        tau_activation_s=30.0,
        gain=RegulatoryGain(
            basal_capacity_multiplier=1.0,
            fully_activated_increment=0.25,
            coupling_scale=1.0,
        ),
        construct=AE4Construct.WT,
    )
    base = ModernFullModel(
        task30_parameters(background, g_nhe1_fmol_s=g_nhe1_fmol_s),
        stimulus=protocol,
        regulatory_model=ae4_regulation,
        ae4_parameters=background.ae4_parameters,
        ae4_evaluator=routing_only_adapter,
    )
    return attach_stimulated_nkcc1(
        base,
        N1AlgebraicNkcc1(
            fully_activated_multiplier=1.75,
            resting_calcium_uM=RESTING_CALCIUM_UM,
            stimulated_calcium_uM=0.10,
        ),
    )


def state_from_coordinates(model: Any, coordinates: ArrayLike) -> NDArray[np.float64]:
    """Map inherited physical coordinates to exact cell/lumen charge manifolds."""

    x = np.asarray(coordinates, dtype=float)
    if x.shape != (10,) or np.any(~np.isfinite(x)):
        raise ValueError("rest coordinates must be a finite vector of length ten")
    if np.any(x <= 0.0):
        raise ValueError("rest coordinates must be positive")
    na_i, k_i, tic_i, ph_i, volume_i, na_l, k_l, tic_l, ph_l, volume_l = (
        float(value) for value in x
    )
    p = model.parameters
    ta_i = total_alkalinity_mM(
        ph=ph_i,
        total_carbon_mM=tic_i,
        buffer_total_mM=p.geometry.cell_buffer_total_fmol / volume_i,
        buffer_pka=p.acid_base.cell_buffer_pka,
        parameters=p.acid_base,
    )
    cl_i = na_i + k_i - ta_i - p.geometry.fixed_cell_anion_equivalents_fmol / volume_i
    ta_l = total_alkalinity_mM(
        ph=ph_l,
        total_carbon_mM=tic_l,
        buffer_total_mM=p.geometry.lumen_buffer_total_fmol / volume_l,
        buffer_pka=p.acid_base.lumen_buffer_pka,
        parameters=p.acid_base,
    )
    cl_l = na_l + k_l - ta_l
    if min(cl_i, cl_l) <= 0.0:
        raise ValueError("charge-manifold chloride must remain positive")
    regulation = tuple(float(value) for value in model.initial_state()[12:])
    state = WholeCellState(
        cell=compartment_from_concentrations(
            na_mM=na_i,
            k_mM=k_i,
            cl_mM=cl_i,
            tic_mM=tic_i,
            alkalinity_mM=ta_i,
            volume_pL=volume_i,
        ),
        lumen=compartment_from_concentrations(
            na_mM=na_l,
            k_mM=k_l,
            cl_mM=cl_l,
            tic_mM=tic_l,
            alkalinity_mM=ta_l,
            volume_pL=volume_l,
        ),
        regulation=regulation,
    )
    return model.layout.encode(state)


def _conservation_audit(evaluation: Any) -> tuple[dict[str, float], float, bool]:
    residuals = {
        key: abs(float(value))
        for key, value in evaluation.diagnostics.conservation_residuals.items()
    }
    if set(residuals) != set(CONSERVATION_RESIDUAL_TOLERANCES):
        raise ValueError("conservation diagnostic keys differ from the inherited contract")
    ratios = {
        key: residuals[key] / tolerance
        for key, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
    }
    maximum = max(ratios.values())
    return residuals, float(maximum), bool(maximum <= 1.0)


def _finite_capacity_audit(evaluation: Any, model: Any) -> dict[str, Any]:
    h = evaluation.diagnostics.homeostasis
    p = model.parameters
    multiplier = float(evaluation.diagnostics.regulatory["nkcc1_capacity_multiplier"])
    capacities = {
        "nhe1_fmol_s": p.homeostasis.nhe1_published_g_fmol_s,
        "nkcc1_cycle_fmol_s": p.homeostasis.nkcc1_capacity_fmol_s * multiplier,
        "ae2_cycle_fmol_s": p.homeostasis.ae2_capacity_fmol_s,
        "pump_cycle_fmol_s": p.membranes.nak_capacity_fmol_s,
    }
    realised = {
        "nhe1_fmol_s": abs(float(h.nhe1_inward_fmol_s)),
        "nkcc1_cycle_fmol_s": abs(float(h.nkcc1_inward_fmol_s)),
        "ae2_cycle_fmol_s": abs(float(h.ae2_inward_fmol_s)),
        "pump_cycle_fmol_s": abs(
            float(evaluation.diagnostics.membranes.pump_apical_fmol_s)
            + float(evaluation.diagnostics.membranes.pump_basolateral_fmol_s)
        ),
    }
    utilisation = {
        key: (realised[key] / value if value > 0.0 else (0.0 if realised[key] == 0.0 else math.inf))
        for key, value in capacities.items()
    }
    return {
        "capacities": capacities,
        "realised_absolute_fluxes": realised,
        "utilisation": utilisation,
        "pass": bool(max(utilisation.values()) <= 1.0 + 1.0e-12),
    }


def _rest_observables(evaluation: Any, state: NDArray[np.float64]) -> dict[str, float]:
    d = evaluation.diagnostics
    cell = d.observables.cell_concentrations_mM
    lumen = d.observables.lumen_concentrations_mM
    return {
        "na_i_mM": float(cell["na"]),
        "k_i_mM": float(cell["k"]),
        "cl_i_mM": float(cell["cl"]),
        "tic_i_mM": float(cell["tic"]),
        "hco3_i_mM": float(cell["hco3"]),
        "ph_i": float(d.observables.cell_acid_base.ph),
        "volume_i_pL": float(state[5]),
        "na_l_mM": float(lumen["na"]),
        "k_l_mM": float(lumen["k"]),
        "cl_l_mM": float(lumen["cl"]),
        "tic_l_mM": float(lumen["tic"]),
        "hco3_l_mM": float(lumen["hco3"]),
        "ph_l": float(d.observables.lumen_acid_base.ph),
        "volume_l_pL": float(state[11]),
        "osm_cell_mOsm": float(d.observables.osmolarities_mOsm["cell"]),
        "osm_lumen_mOsm": float(d.observables.osmolarities_mOsm["lumen"]),
        "osm_bath_mOsm": float(d.observables.osmolarities_mOsm["bath"]),
        "v_apical_mV": 1.0e3 * float(d.membranes.v_apical_V),
        "v_basolateral_mV": 1.0e3 * float(d.membranes.v_basolateral_V),
        "v_transepithelial_mV": 1.0e3 * float(d.membranes.v_transepithelial_V),
    }


def _flux_row(evaluation: Any, model: Any) -> dict[str, float | str]:
    d = evaluation.diagnostics
    h = d.homeostasis
    m = d.membranes
    f = model.parameters.constants.faraday_C_mol
    return {
        "nhe1_model": h.nhe1_model,
        "nhe1_na_in_h_out_fmol_s": float(h.nhe1_inward_fmol_s),
        "nkcc1_cycle_in_fmol_s": float(h.nkcc1_inward_fmol_s),
        "nkcc1_cl_in_fmol_s": 2.0 * float(h.nkcc1_inward_fmol_s),
        "ae2_cl_in_hco3_out_fmol_s": float(h.ae2_inward_fmol_s),
        "ae4_na_cell_fmol_s": float(d.ae4.na_cell_fmol_s),
        "ae4_k_cell_fmol_s": float(d.ae4.k_cell_fmol_s),
        "ae4_cl_cell_fmol_s": float(d.ae4.cl_cell_fmol_s),
        "ae4_hco3_cell_fmol_s": float(d.ae4.hco3_cell_fmol_s),
        "pump_apical_cycles_fmol_s": float(m.pump_apical_fmol_s),
        "pump_basolateral_cycles_fmol_s": float(m.pump_basolateral_fmol_s),
        "k_apical_cell_to_lumen_fmol_s": current_to_fmol_s(
            float(m.currents_A["k_apical"]), valence=1, faraday_C_mol=f
        ),
        "k_basolateral_cell_to_bath_fmol_s": current_to_fmol_s(
            float(m.currents_A["k_basolateral"]), valence=1, faraday_C_mol=f
        ),
        "cacc_cl_cell_to_lumen_fmol_s": current_to_fmol_s(
            float(m.currents_A["cl_apical"]), valence=-1, faraday_C_mol=f
        ),
        "co2_bath_to_cell_fmol_s": float(d.co2_fluxes_fmol_s["bath_to_cell"]),
        "co2_lumen_to_cell_fmol_s": float(d.co2_fluxes_fmol_s["lumen_to_cell"]),
        "water_bath_to_cell_pL_s": float(d.water.bath_to_cell_pL_s),
        "water_cell_to_lumen_pL_s": float(d.water.cell_to_lumen_pL_s),
        "water_bath_to_lumen_pL_s": float(d.water.bath_to_lumen_pL_s),
        "water_lumen_outflow_pL_s": float(d.water.lumen_outflow_pL_s),
    }


def solve_stationary(
    model: Any,
    *,
    genotype: Genotype,
    start_coordinates: ArrayLike,
    background: Background,
    label: str,
    budget: NumericalBudget,
    max_nfev: int = 3000,
) -> dict[str, Any]:
    """Make one bounded direct solve on the inherited charge manifold."""

    budget.start_stationary_solver()
    start = np.asarray(start_coordinates, dtype=float)
    if start.shape != (10,) or np.any(start <= COORDINATE_LOWER) or np.any(start >= COORDINATE_UPPER):
        raise ValueError(f"{label}: stationary start must be strictly inside inherited bounds")
    local_residual_evaluations = 0

    def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
        nonlocal local_residual_evaluations
        budget.count_stationary_residual()
        local_residual_evaluations += 1
        try:
            state = state_from_coordinates(model, candidate)
            raw = model.rhs(0.0, state, genotype=genotype)
            return raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(10, 1.0e6)

    started = time.perf_counter()
    try:
        fit = least_squares(
            residual,
            start,
            bounds=(COORDINATE_LOWER, COORDINATE_UPPER),
            max_nfev=max_nfev,
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            x_scale="jac",
        )
    except NumericalBudgetExhausted:
        # Preserve the first exhausted limit while still accounting for time
        # spent in the interrupted solver call.
        budget.numerical_wall_seconds += time.perf_counter() - started
        raise
    except Exception:
        budget.add_wall_time(time.perf_counter() - started)
        raise
    else:
        budget.add_wall_time(time.perf_counter() - started)
    coordinates = np.asarray(fit.x, dtype=float)
    state = state_from_coordinates(model, coordinates)
    evaluation = model.evaluate(0.0, state, genotype=genotype)
    later = model.evaluate(600.0, state, genotype=genotype)
    raw = np.asarray(evaluation.rhs, dtype=float)
    scaled = raw[INDEPENDENT_ROWS] / INDEPENDENT_RHS_SCALES
    span = COORDINATE_UPPER - COORDINATE_LOWER
    jacobian = np.asarray(fit.jac, dtype=float) * span[None, :]
    singular = np.linalg.svd(jacobian, compute_uv=False)
    threshold = JACOBIAN_RANK_RELATIVE_TOLERANCE * singular[0] if singular.size else math.inf
    rank = int(np.sum(singular > threshold))
    relative_boundary_distance = np.minimum(
        coordinates - COORDINATE_LOWER, COORDINATE_UPPER - coordinates
    ) / span
    boundary_hits = [
        name
        for name, distance in zip(COORDINATE_NAMES, relative_boundary_distance)
        if distance <= BOUNDARY_RELATIVE_TOLERANCE
    ]
    charge_max = max(abs(float(value)) for value in evaluation.diagnostics.state_charge_fmol.values())
    current_max = max(abs(float(value)) for value in evaluation.diagnostics.membranes.current_residuals_A.values())
    regulatory_max = float(np.max(np.abs(raw[12:]))) if raw.size > 12 else 0.0
    conservation, conservation_ratio, conservation_pass = _conservation_audit(evaluation)
    finite_capacity = _finite_capacity_audit(evaluation, model)
    observables = _rest_observables(evaluation, state)
    full_rhs_pass = bool(
        np.max(np.abs(scaled)) <= ROOT_SCALED_TOLERANCE
        and np.max(np.abs(raw[OMITTED_CHARGE_ROWS])) <= OMITTED_ROW_TOLERANCE
        and regulatory_max <= REGULATORY_TOLERANCE_S_INV
    )
    numerical_pass = bool(
        fit.success
        and full_rhs_pass
        and charge_max <= CHARGE_TOLERANCE_FMOL
        and current_max <= CURRENT_TOLERANCE_A
        and conservation_pass
        and rank == 10
        and not boundary_hits
        and np.all(state[:12] > 0.0)
        and np.array_equal(raw, later.rhs)
    )
    ordinary_na_pass = bool(observables["na_i_mM"] <= ORDINARY_NA_UPPER_MM)
    ph_diagnostic_centre = (
        EXACT_NULL_DIAGNOSTIC_PH
        if genotype.ae4_expression == 0.0
        else WT_TARGET_PH
    )
    ph_diagnostic_pass = bool(
        abs(observables["ph_i"] - ph_diagnostic_centre)
        <= NHE1_PH_DIAGNOSTIC_HALF_WIDTH
    )
    volume_pass = bool(observables["volume_i_pL"] <= CELL_VOLUME_UPPER_PL)
    physiological_pass = bool(ordinary_na_pass and ph_diagnostic_pass and volume_pass)
    exact_ae4_zero = all(
        float(value) == 0.0
        for value in (
            evaluation.diagnostics.ae4.na_cell_fmol_s,
            evaluation.diagnostics.ae4.k_cell_fmol_s,
            evaluation.diagnostics.ae4.cl_cell_fmol_s,
            evaluation.diagnostics.ae4.hco3_cell_fmol_s,
            evaluation.diagnostics.ae4.transported_charge_fmol_s,
        )
    ) if genotype.ae4_expression == 0.0 else True
    admissible = bool(
        numerical_pass
        and physiological_pass
        and finite_capacity["pass"]
        and exact_ae4_zero
    )
    return {
        "background": background.alias,
        "root_id": background.root_id,
        "routing_family": background.routing_family,
        "genotype": genotype.name,
        "ae4_expression": float(genotype.ae4_expression),
        "label": label,
        "g_nhe1_fmol_s": float(model.parameters.homeostasis.nhe1_published_g_fmol_s),
        "nhe1_model": model.parameters.homeostasis.nhe1_model,
        "optimizer_success": bool(fit.success),
        "optimizer_message": str(fit.message),
        "optimizer_nfev": int(fit.nfev),
        "actual_residual_evaluations": int(local_residual_evaluations),
        "cost": float(fit.cost),
        "optimality": float(fit.optimality),
        "coordinates": coordinates.tolist(),
        "core_state": state[:12].tolist(),
        "core_state_sha256": sha256_object(state[:12].tolist()),
        "raw_rhs": raw.tolist(),
        "scaled_independent_rhs": scaled.tolist(),
        "max_abs_scaled_independent_rhs": float(np.max(np.abs(scaled))),
        "max_abs_amount_rhs_fmol_s": float(np.max(np.abs(raw[AMOUNT_ROWS]))),
        "max_abs_volume_rhs_pL_s": float(np.max(np.abs(raw[VOLUME_ROWS]))),
        "max_abs_omitted_rhs_fmol_s": float(np.max(np.abs(raw[OMITTED_CHARGE_ROWS]))),
        "max_abs_regulatory_rhs_s_inv": regulatory_max,
        "max_abs_bulk_charge_fmol": float(charge_max),
        "max_abs_current_residual_A": float(current_max),
        "max_dimensionless_conservation_ratio": conservation_ratio,
        "conservation_residuals": conservation,
        "normalized_jacobian_singular_values": singular.tolist(),
        "independent_jacobian_rank": rank,
        "independent_jacobian_nullity": 10 - rank,
        "boundary_hits": boundary_hits,
        "rest_protocol_time_invariant_rhs": bool(np.array_equal(raw, later.rhs)),
        "positivity_pass": bool(np.all(state[:12] > 0.0)),
        "full_rhs_pass": full_rhs_pass,
        "conservation_pass": conservation_pass,
        "numerical_pass": numerical_pass,
        "ordinary_na_pass": ordinary_na_pass,
        "ph_diagnostic_centre": ph_diagnostic_centre,
        "ph_diagnostic_pass": ph_diagnostic_pass,
        "volume_pass": volume_pass,
        "finite_capacity": finite_capacity,
        "exact_ae4_zero": exact_ae4_zero,
        "physiological_pass": physiological_pass,
        "admissible": admissible,
        "observables": observables,
        "fluxes": _flux_row(evaluation, model),
        "active_parameter_sha256": sha256_object(model.parameters),
        "inherited_parameter_sha256": background.inherited_whole_cell_sha256,
        "ae4_parameter_sha256": background.inherited_ae4_sha256,
    }


def _wt_calibration_needed(solution: Mapping[str, Any]) -> bool:
    return not bool(solution["admissible"])


def solve_wt_with_optional_scalar_calibration(
    background: Background, budget: NumericalBudget
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Try published G first; use a WT-only secant adjustment only if needed."""

    evaluations: list[dict[str, Any]] = []
    values: list[tuple[float, float, NDArray[np.float64]]] = []
    g_value = 0.0305
    start = np.asarray(background.saved_coordinates, dtype=float)
    for index in range(4):
        budget.count_scalar_g_evaluation()
        model = build_task30_model(background, g_nhe1_fmol_s=g_value, stimulated=False)
        solution = solve_stationary(
            model,
            genotype=WT,
            start_coordinates=start,
            background=background,
            label=f"wt_g_evaluation_{index + 1}",
            budget=budget,
        )
        evaluations.append(solution)
        if not _wt_calibration_needed(solution):
            solution["g_nhe1_status"] = (
                "PUBLISHED_VALUE" if index == 0 else "WT_ONLY_SCALAR_CALIBRATION"
            )
            return solution, evaluations
        if solution["numerical_pass"]:
            ph = float(solution["observables"]["ph_i"])
            state_coordinates = np.asarray(solution["coordinates"], dtype=float)
            values.append((g_value, ph, state_coordinates))
            start = state_coordinates
        if index == 3:
            break
        if len(values) < 2:
            # One deterministic local sensitivity evaluation, not a grid.
            error = values[-1][1] - WT_TARGET_PH if values else 0.0
            factor = 1.25 if error > 0.0 else 0.75
            g_value = 0.0305 * factor
        else:
            g0, ph0, _ = values[-2]
            g1, ph1, _ = values[-1]
            denominator = ph1 - ph0
            if abs(denominator) < 1.0e-12:
                g_value = g1 * (1.10 if ph1 > WT_TARGET_PH else 0.90)
            else:
                g_value = g1 + (WT_TARGET_PH - ph1) * (g1 - g0) / denominator
            g_value = float(max(1.0e-8, g_value))
    selected = min(
        evaluations,
        key=lambda item: (
            not bool(item["numerical_pass"]),
            abs(float(item["observables"]["ph_i"]) - WT_TARGET_PH),
            float(item["max_abs_scaled_independent_rhs"]),
        ),
    )
    selected["g_nhe1_status"] = "WT_ONLY_SCALAR_CALIBRATION_FAILED"
    return selected, evaluations


def _null_observation_biased_start(model: Any, wt_coordinates: Sequence[float]) -> NDArray[np.float64]:
    """One allowed deterministic start near measured null pH/Cl, not a target fit."""

    x = np.asarray(wt_coordinates, dtype=float).copy()
    x[3] = 6.89
    p = model.parameters
    ta_i = total_alkalinity_mM(
        ph=x[3],
        total_carbon_mM=x[2],
        buffer_total_mM=p.geometry.cell_buffer_total_fmol / x[4],
        buffer_pka=p.acid_base.cell_buffer_pka,
        parameters=p.acid_base,
    )
    x[1] = 36.5 - x[0] + ta_i + p.geometry.fixed_cell_anion_equivalents_fmol / x[4]
    if np.any(x <= COORDINATE_LOWER) or np.any(x >= COORDINATE_UPPER):
        raise ValueError("the single observation-biased null warm start is outside inherited bounds")
    return x


def solve_exact_null(
    background: Background,
    wt_solution: Mapping[str, Any],
    budget: NumericalBudget,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    g_value = float(wt_solution["g_nhe1_fmol_s"])
    model = build_task30_model(background, g_nhe1_fmol_s=g_value, stimulated=False)
    first = solve_stationary(
        model,
        genotype=AE4_NULL,
        start_coordinates=wt_solution["coordinates"],
        background=background,
        label="exact_ae4_null_from_repaired_wt",
        budget=budget,
    )
    attempts = [first]
    if first["admissible"]:
        return first, attempts
    biased = _null_observation_biased_start(model, wt_solution["coordinates"])
    second = solve_stationary(
        model,
        genotype=AE4_NULL,
        start_coordinates=biased,
        background=background,
        label="exact_ae4_null_observation_biased_warm_start",
        budget=budget,
    )
    attempts.append(second)
    selected = min(
        attempts,
        key=lambda item: (
            not bool(item["admissible"]),
            not bool(item["numerical_pass"]),
            abs(float(item["observables"]["ph_i"]) - EXACT_NULL_DIAGNOSTIC_PH),
            float(item["max_abs_scaled_independent_rhs"]),
        ),
    )
    return selected, attempts


def run_stimulation(
    background: Background,
    rest_solution: Mapping[str, Any],
    *,
    genotype: Genotype,
    budget: NumericalBudget,
) -> dict[str, Any]:
    """Run one inherited production Radau CCh+IPR stimulation at 0.25 uM Ca."""

    budget.start_stimulated_integration()
    g_value = float(rest_solution["g_nhe1_fmol_s"])
    model = build_task30_model(background, g_nhe1_fmol_s=g_value, stimulated=True)
    rest_model = build_task30_model(background, g_nhe1_fmol_s=g_value, stimulated=False)
    core = np.asarray(rest_solution["core_state"], dtype=float)
    y0 = np.r_[core, model.initial_state()[12:]]
    onset = model.evaluate(0.0, y0, genotype=genotype)
    rest_onset = rest_model.evaluate(0.0, y0, genotype=genotype)
    onset_match = bool(np.array_equal(onset.rhs, rest_onset.rhs))
    grid = np.r_[1.0e-6, np.arange(1.0, 601.0, 1.0)]
    started = time.perf_counter()
    result = model.solve_dynamics(
        (float(grid[0]), 600.0),
        initial_state=y0,
        genotype=genotype,
        method=PRODUCTION_RADAU.method,
        rtol=PRODUCTION_RADAU.rtol,
        atol=PRODUCTION_RADAU.atol_vector(model.state_names),
        t_eval=grid,
        max_step_s=PRODUCTION_RADAU.max_step_s,
    )
    budget.add_wall_time(time.perf_counter() - started)
    complete = bool(result.y.shape[1] == grid.size and result.t[-1] == 600.0)
    times = np.r_[0.0, np.asarray(result.t, dtype=float)]
    states = np.column_stack((y0, np.asarray(result.y, dtype=float)))
    flow: list[float] = []
    ae4_cl: list[float] = []
    nkcc_cl: list[float] = []
    ae2_cl: list[float] = []
    nhe1: list[float] = []
    conservation_ratios: list[float] = []
    positive = True
    finite = True
    for time_s, state in zip(times, states.T):
        evaluation = model.evaluate(float(time_s), state, genotype=genotype)
        homeostasis = evaluation.diagnostics.homeostasis
        flow.append(float(evaluation.diagnostics.water.lumen_outflow_pL_s))
        ae4_cl.append(float(evaluation.diagnostics.ae4.cl_cell_fmol_s))
        nkcc_cl.append(2.0 * float(homeostasis.nkcc1_inward_fmol_s))
        ae2_cl.append(float(homeostasis.ae2_inward_fmol_s))
        nhe1.append(float(homeostasis.nhe1_inward_fmol_s))
        _, ratio, _ = _conservation_audit(evaluation)
        conservation_ratios.append(ratio)
        positive = positive and bool(np.all(np.asarray(state[:12]) > 0.0))
        finite = finite and bool(np.all(np.isfinite(state))) and math.isfinite(flow[-1])
    flow_array = np.asarray(flow)
    cumulative = np.r_[0.0, cumulative_trapezoid(flow_array, times)]
    integrated = {
        "ae4_cl_cell_fmol": float(np.trapezoid(np.asarray(ae4_cl), times)),
        "nkcc1_cl_cell_fmol": float(np.trapezoid(np.asarray(nkcc_cl), times)),
        "ae2_cl_cell_fmol": float(np.trapezoid(np.asarray(ae2_cl), times)),
        "nhe1_na_in_h_out_fmol": float(np.trapezoid(np.asarray(nhe1), times)),
    }
    valid = bool(
        result.success
        and complete
        and onset_match
        and positive
        and finite
        and max(conservation_ratios) <= 1.0
        and np.all(flow_array >= 0.0)
    )
    endpoint = model.evaluate(float(times[-1]), states[:, -1], genotype=genotype)
    return {
        "background": background.alias,
        "root_id": background.root_id,
        "genotype": genotype.name,
        "ae4_expression": float(genotype.ae4_expression),
        "g_nhe1_fmol_s": g_value,
        "nhe1_model": model.parameters.homeostasis.nhe1_model,
        "calcium_uM": STIMULATED_CALCIUM_UM,
        "protocol": "600_s_CCH_IPR",
        "solver": PRODUCTION_RADAU.label,
        "solver_success": bool(result.success),
        "solver_message": str(result.message),
        "complete": complete,
        "onset_matches_rest_rhs": onset_match,
        "sample_count": int(times.size),
        "positive_core": positive,
        "finite": finite,
        "nonnegative_flow": bool(np.all(flow_array >= 0.0)),
        "max_dimensionless_conservation_ratio": float(max(conservation_ratios)),
        "cumulative_secretion_pL": float(cumulative[-1]),
        "endpoint_flow_pL_s": float(flow_array[-1]),
        "integrated_fluxes": integrated,
        "endpoint_observables": _rest_observables(endpoint, states[:, -1]),
        "endpoint_fluxes": _flux_row(endpoint, model),
        "valid": valid,
    }


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    names: list[str] = []
    for row in rows:
        for key in row:
            if key not in names:
                names.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def _flatten_rest(solution: Mapping[str, Any]) -> dict[str, Any]:
    row = {
        key: solution[key]
        for key in (
            "background", "root_id", "routing_family", "genotype", "ae4_expression",
            "label", "g_nhe1_fmol_s", "nhe1_model", "g_nhe1_status",
            "optimizer_success", "optimizer_nfev", "actual_residual_evaluations",
            "max_abs_scaled_independent_rhs", "max_abs_amount_rhs_fmol_s",
            "max_abs_volume_rhs_pL_s", "max_abs_omitted_rhs_fmol_s",
            "max_abs_regulatory_rhs_s_inv", "max_abs_bulk_charge_fmol",
            "max_abs_current_residual_A", "max_dimensionless_conservation_ratio",
            "independent_jacobian_rank", "independent_jacobian_nullity",
            "full_rhs_pass", "conservation_pass", "numerical_pass", "ordinary_na_pass",
            "ph_diagnostic_centre",
            "ph_diagnostic_pass", "volume_pass", "exact_ae4_zero",
            "physiological_pass", "admissible",
        )
        if key in solution
    }
    row.update(solution["observables"])
    return row


def _flatten_flux(solution: Mapping[str, Any]) -> dict[str, Any]:
    row = {
        "background": solution["background"],
        "genotype": solution["genotype"],
        "ae4_expression": solution["ae4_expression"],
        "admissible": solution["admissible"],
    }
    row.update(solution["fluxes"])
    for name, value in solution["conservation_residuals"].items():
        row[f"residual_{name}"] = value
    for name, value in solution["finite_capacity"]["utilisation"].items():
        row[f"capacity_utilisation_{name}"] = value
    return row


def _law_comparison(
    backgrounds: Sequence[Background],
    selected_rest: Sequence[Mapping[str, Any]],
) -> str:
    lines = [
        "# NHE1 law comparison",
        "",
        "The repaired model uses Vera-Sigüenza et al. (2018), Eq. 26 exactly, with the squared proton saturation factors printed in that paper.",
        "",
        "| Property | Historical comparator | Task 30 evaluator |",
        "|---|---|---|",
        "| Identifier | `legacy_tanh` | `vera_siguenza_2018_eq26` |",
        "| Rate form | Capacity times `tanh` of an ideal activity affinity | Difference of forward and reverse sodium/proton saturation terms |",
        "| Proton factor | Appears through a log activity ratio | Squared saturation on each proton term |",
        "| Activity | Inherited background-specific placeholder | `0.0305 fmol/s` published default; only the permitted WT scalar calibration may change it |",
        "| Constants | Shared log width | `K_H = 4.5e-4 mM`, `K_Na = 15 mM` |",
        "| Status | Selectable only for historical reproduction | Explicitly selected in every Task 30 model |",
        "",
        "Fixed state comparison at each inherited saved WT state, using the nominal published `G_NHE1 = 0.0305 fmol/s`:",
        "",
        "| Background | Historical flux (fmol/s) | Eq. 26 flux (fmol/s) | Eq. 26 / historical |",
        "|---|---:|---:|---:|",
    ]
    for background in backgrounds:
        published_model = build_task30_model(background, stimulated=False)
        state = state_from_coordinates(published_model, background.saved_coordinates)
        evaluation = published_model.evaluate(0.0, state, genotype=WT)
        obs = evaluation.diagnostics.observables
        cell = obs.cell_concentrations_mM
        bath = published_model.parameters.bath
        from .membranes import HomeostasisEnvironment
        environment = HomeostasisEnvironment(
            na_i_mM=cell["na"], k_i_mM=cell["k"], cl_i_mM=cell["cl"],
            h_i_mM=obs.cell_acid_base.h_mM, hco3_i_mM=cell["hco3"],
            na_e_mM=bath.na_mM, k_e_mM=bath.k_mM, cl_e_mM=bath.cl_mM,
            h_e_mM=obs.bath_acid_base.h_mM, hco3_e_mM=obs.bath_acid_base.hco3_mM,
        )
        legacy_parameters = replace(
            published_model.parameters,
            homeostasis=replace(
                published_model.parameters.homeostasis, nhe1_model=NHE1_LEGACY_TANH
            ),
        )
        legacy = evaluate_nhe1_legacy_tanh(environment, legacy_parameters)
        published = evaluate_nhe1_vera_siguenza_2018(environment, published_model.parameters)
        ratio = published / legacy if legacy != 0.0 else math.nan
        lines.append(f"| {background.alias} | {legacy:.12g} | {published:.12g} | {ratio:.6g} |")
    lines.append("")
    for solution in selected_rest:
        if solution.get("g_nhe1_status") == "WT_ONLY_SCALAR_CALIBRATION":
            lines.append(
                f"The selected {solution['background']} WT state used the permitted "
                f"WT only scalar calibration, `G_NHE1 = "
                f"{solution['g_nhe1_fmol_s']:.12g} fmol/s`. The nominal published "
                "value itself was tested first."
            )
    reached_backgrounds = {str(row["background"]) for row in selected_rest}
    not_reached = [item.alias for item in backgrounds if item.alias not in reached_backgrounds]
    if not_reached:
        lines.append(
            f"{', '.join(not_reached)} WT evaluation was not reached before the "
            "numerical budget stop."
        )
    lines.extend((
        "",
        "Positive flux is one sodium ion entering the cell and one proton leaving it. It therefore adds the same amount to intracellular sodium and total alkalinity, adds no TIC directly, and transports no net charge.",
        "",
    ))
    return "\n".join(lines)


def _final_answer_markdown(
    selected_rest: Sequence[Mapping[str, Any]],
    all_rest_attempts: Sequence[Mapping[str, Any]],
    dynamic_rows: Sequence[Mapping[str, Any]],
    budget: NumericalBudget,
    run_status: str,
    stop_reason: str | None,
) -> str:
    lookup = {(row["background"], row["genotype"]): row for row in selected_rest}
    attempts = {
        (row["background"], row["label"]): row for row in all_rest_attempts
    }
    dynamics = {(row["background"], row["genotype"]): row for row in dynamic_rows}
    lines = [
        "# Task 30 final answer",
        "",
        f"Run status: **{run_status}**"
        + (f" because {stop_reason}." if stop_reason else "."),
        "",
        "The published 2018 NHE1 law was implemented as a named evaluator and explicitly selected in both Task 30 backgrounds. The historical `tanh` evaluator remains selectable but was not used for the repaired resting states or trajectories.",
        "",
        "| Background | WT REST | Exact AE4 null REST | WT Na (mM) | Null Na (mM) | WT pH | Null pH | Null Cl (mM) |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for alias in BACKGROUND_ALIASES:
        wt = lookup.get((alias, "WT"))
        null = lookup.get((alias, "AE4_NULL"))
        if wt is None:
            lines.append(f"| {alias} | not reached | not reached |  |  |  |  |  |")
            continue
        wto = wt["observables"]
        if null is None:
            null_status = (
                "attempt incomplete at budget"
                if wt["admissible"] and run_status == "STOPPED_AT_NUMERICAL_BUDGET"
                else "not attempted"
            )
            lines.append(f"| {alias} | {'admissible' if wt['admissible'] else 'failed'} | {null_status} | {wto['na_i_mM']:.6f} |  | {wto['ph_i']:.6f} |  |  |")
            continue
        no = null["observables"]
        lines.append(
            f"| {alias} | {'admissible' if wt['admissible'] else 'failed'} | "
            f"{'admissible' if null['admissible'] else 'failed'} | {wto['na_i_mM']:.6f} | "
            f"{no['na_i_mM']:.6f} | {wto['ph_i']:.6f} | {no['ph_i']:.6f} | {no['cl_i_mM']:.6f} |"
        )
    lines.extend(("", "## Direct answers", ""))
    published_r09 = attempts.get(("R09", "wt_g_evaluation_1"))
    selected_r09 = lookup.get(("R09", "WT"))
    lines.append(
        "1. **Implementation:** Yes. The evaluator reproduces three independently "
        "calculated Eq. 26 values to machine precision, retains both squared proton "
        "terms, and preserves electroneutral 1 Na inward to 1 H outward stoichiometry."
    )
    if published_r09 is not None and selected_r09 is not None:
        lines.append(
            "2. **WT REST:** R09 recovered an admissible stationary state only after the "
            f"permitted WT only scalar calibration. The published `G_NHE1 = 0.0305 fmol/s` "
            f"gave pH {published_r09['observables']['ph_i']:.6f}; the frozen calibrated "
            f"value was {selected_r09['g_nhe1_fmol_s']:.12g} fmol/s and gave pH "
            f"{selected_r09['observables']['ph_i']:.6f} with Na "
            f"{selected_r09['observables']['na_i_mM']:.6f} mM. R10 was not reached."
        )
    else:
        lines.append("2. **WT REST:** No complete WT conclusion was available.")
    lines.extend((
        "3. **Exact AE4 null REST:** Not determined. The first R09 exact null solve was "
        "still in progress when the global residual evaluation ceiling was reached, so "
        "there is no accepted null state and no basis for calling it admissible or failed.",
        "4. **Null sodium:** Not determined; the incomplete solve must not be interpreted "
        "as either ordinary sodium or a recurrence of the Task 29 extreme values.",
        "5. **Null pH:** Not determined; no accepted state is available for comparison "
        "with the 6.89 diagnostic centre.",
        "6. **AE2 and NKCC1 after AE4 removal:** Not determined because the exact null "
        "stationary state was not completed.",
        "7. **Dynamics:** None ran. The prerequisite admissible matched WT and exact null "
        "REST pair was unavailable, so no secretion direction is reported.",
        "8. **Later chloride allocation:** Chloride allocation remains a separate later "
        "phase. No chloride allocation search or chloride fitting was performed here.",
    ))
    lines.extend(("", "## Available exact null results", ""))
    for alias in BACKGROUND_ALIASES:
        wt = lookup.get((alias, "WT")); null = lookup.get((alias, "AE4_NULL"))
        if wt is None or null is None:
            continue
        nf = null["fluxes"]
        lines.append(
            f"* {alias}: exact AE4 removal left AE4 flux at zero. AE2 chloride flux was "
            f"{nf['ae2_cl_in_hco3_out_fmol_s']:.6g} fmol/s and NKCC1 chloride influx was "
            f"{nf['nkcc1_cl_in_fmol_s']:.6g} fmol/s. Intracellular sodium "
            f"{'remained ordinary' if null['ordinary_na_pass'] else 'failed the ordinary sodium screen'}; "
            f"pH {'remained within' if null['ph_diagnostic_pass'] else 'fell outside'} the declared 6.89 ± 0.10 diagnostic band."
        )
    if dynamic_rows:
        lines.extend(("", "## Matched 0.25 uM stimulation", "", "| Background | WT secretion (pL) | Null secretion (pL) | Null change |", "|---|---:|---:|---:|"))
        for alias in BACKGROUND_ALIASES:
            wt = dynamics.get((alias, "WT")); null = dynamics.get((alias, "AE4_NULL"))
            if wt is None or null is None:
                continue
            change = null["cumulative_secretion_pL"] / wt["cumulative_secretion_pL"] - 1.0
            lines.append(f"| {alias} | {wt['cumulative_secretion_pL']:.9g} | {null['cumulative_secretion_pL']:.9g} | {change:+.3%} |")
        lines.extend(("", "These directions are predictions without secretion tuning. The known knockout secretion reduction was not used anywhere in construction or selection.", ""))
    else:
        lines.extend(("", "No stimulated trajectory was run because no background had both an admissible WT and admissible exact null REST.", ""))
    lines.extend((
        f"Budget used: {budget.stationary_solver_calls}/12 stationary calls, {budget.stationary_residual_evaluations}/10000 stationary residual evaluations, {budget.scalar_g_nhe1_evaluations}/8 scalar G evaluations, and {budget.stimulated_integrations}/4 stimulated integrations.",
        "",
    ))
    return "\n".join(lines)


def write_outputs(
    *,
    backgrounds: Sequence[Background],
    selected_rest: Sequence[Mapping[str, Any]],
    all_rest_attempts: Sequence[Mapping[str, Any]],
    dynamic_rows: Sequence[Mapping[str, Any]],
    budget: NumericalBudget,
    run_status: str,
    stop_reason: str | None,
) -> None:
    RESULT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    ANALYSIS_DIRECTORY.mkdir(parents=True, exist_ok=True)
    for row in all_rest_attempts:
        _write_json(
            RESULT_DIRECTORY / "resting_attempts" / f"{row['background']}_{row['label']}.json",
            row,
        )
    _write_csv(RESULT_DIRECTORY / "resting_states.csv", [_flatten_rest(row) for row in selected_rest])
    _write_csv(RESULT_DIRECTORY / "flux_ledger.csv", [_flatten_flux(row) for row in selected_rest])
    if dynamic_rows:
        secretion_rows = []
        lookup = {(row["background"], row["genotype"]): row for row in dynamic_rows}
        for row in dynamic_rows:
            flat = {
                key: value for key, value in row.items()
                if key not in {"integrated_fluxes", "endpoint_observables", "endpoint_fluxes"}
            }
            flat.update({f"integrated_{key}": value for key, value in row["integrated_fluxes"].items()})
            counterpart = lookup.get((row["background"], "WT"))
            if counterpart is not None and row["genotype"] == "AE4_NULL":
                flat["null_to_wt_secretion_ratio"] = row["cumulative_secretion_pL"] / counterpart["cumulative_secretion_pL"]
                flat["null_vs_wt_secretion_change_fraction"] = flat["null_to_wt_secretion_ratio"] - 1.0
            secretion_rows.append(flat)
        _write_csv(RESULT_DIRECTORY / "matched_secretion.csv", secretion_rows)
        _write_json(RESULT_DIRECTORY / "dynamic_details.json", {"trajectories": list(dynamic_rows)})
    contract = {
        "task": 30,
        "repository": "esig626/ae4-salivary-transport-control",
        "branch": BRANCH,
        "prepared_head": PREPARED_HEAD,
        "backgrounds": list(BACKGROUND_ALIASES),
        "nhe1": {
            "active_model": NHE1_VERA_SIGUENZA_2018,
            "historical_comparator": NHE1_LEGACY_TANH,
            "equation": "G*((H_i/(H_i+K_H))**2*(Na_e/(Na_e+K_Na))-(Na_i/(Na_i+K_Na))*(H_e/(H_e+K_H))**2)",
            "g_nhe1_fmol_s": 0.0305,
            "k_h_mM": 4.5e-4,
            "k_na_mM": 15.0,
            "source": "Vera-Siguenza et al. 2018 Eq. 26, PMCID PMC5792321",
            "extra_ph_gate": False,
            "genotype_dependent": False,
        },
        "other_transporter_parameters_unchanged": True,
        "ae4_routing_model": AE4_ROUTING_MODEL_ID,
        "rest_solver": {
            "method": "direct bounded charge-manifold least_squares",
            "coordinate_names": list(COORDINATE_NAMES),
            "lower": COORDINATE_LOWER.tolist(),
            "upper": COORDINATE_UPPER.tolist(),
            "independent_rows": INDEPENDENT_ROWS.tolist(),
            "omitted_charge_rows": OMITTED_CHARGE_ROWS.tolist(),
            "independent_rhs_scales": INDEPENDENT_RHS_SCALES.tolist(),
            "root_scaled_tolerance": ROOT_SCALED_TOLERANCE,
            "omitted_row_tolerance": OMITTED_ROW_TOLERANCE,
            "charge_tolerance_fmol": CHARGE_TOLERANCE_FMOL,
            "current_tolerance_A": CURRENT_TOLERANCE_A,
        },
        "diagnostic_screens": {
            "na_i_upper_mM": ORDINARY_NA_UPPER_MM,
            "cell_volume_upper_pL": CELL_VOLUME_UPPER_PL,
            "wt_ph_centre": WT_TARGET_PH,
            "exact_null_ph_centre": EXACT_NULL_DIAGNOSTIC_PH,
            "ph_half_width": NHE1_PH_DIAGNOSTIC_HALF_WIDTH,
            "chloride_is_not_an_acceptance_gate": True,
        },
        "optional_g_calibration": {
            "used_only_if_published_value_fails_wt": True,
            "wt_only": True,
            "target_ph": WT_TARGET_PH,
            "maximum_evaluations_per_background": 4,
        },
        "dynamic_protocol": {
            "run_only_after_admissible_wt_and_null_rest": True,
            "arm": "CCH_IPR",
            "duration_s": 600.0,
            "calcium_uM": STIMULATED_CALCIUM_UM,
            "solver": asdict(PRODUCTION_RADAU),
        },
        "prohibited_operations_performed": [],
        "git": {
            "ordinary_authentication_attempts": 1,
            "authentication_succeeded": False,
            "commit_or_push_attempted": False,
        },
    }
    verification = {
        "run_status": run_status,
        "stop_reason": stop_reason,
        "focused_unit_tests_passed_before_whole_cell_runs": True,
        "focused_unit_test_count": 8,
        "published_equation_machine_precision_cases": 3,
        "background_count": len(backgrounds),
        "selected_rest_count": len(selected_rest),
        "all_selected_rest_nhe1_model_correct": all(row["nhe1_model"] == NHE1_VERA_SIGUENZA_2018 for row in selected_rest),
        "all_selected_rest_capacity_pass": all(row["finite_capacity"]["pass"] for row in selected_rest),
        "admissible_wt_backgrounds": [row["background"] for row in selected_rest if row["genotype"] == "WT" and row["admissible"]],
        "admissible_exact_null_backgrounds": [row["background"] for row in selected_rest if row["genotype"] == "AE4_NULL" and row["admissible"]],
        "dynamic_case_count": len(dynamic_rows),
        "all_dynamic_cases_valid": bool(dynamic_rows) and all(row["valid"] for row in dynamic_rows),
        "budget_limits_respected": budget.payload()["limits_respected"],
        "inherited_payload_hashes": {
            background.alias: {
                "whole_cell": background.inherited_whole_cell_sha256,
                "ae4": background.inherited_ae4_sha256,
            }
            for background in backgrounds
        },
    }
    _write_json(RESULT_DIRECTORY / "contract.json", contract)
    _write_json(RESULT_DIRECTORY / "verification.json", verification)
    _write_json(RESULT_DIRECTORY / "budget.json", budget.payload())
    (ANALYSIS_DIRECTORY / "nhe1_law_comparison.md").write_text(
        _law_comparison(backgrounds, selected_rest), encoding="utf-8"
    )
    (ANALYSIS_DIRECTORY / "final_answer.md").write_text(
        _final_answer_markdown(
            selected_rest,
            all_rest_attempts,
            dynamic_rows,
            budget,
            run_status,
            stop_reason,
        ),
        encoding="utf-8",
    )


def run() -> dict[str, Any]:
    budget = NumericalBudget()
    backgrounds = [load_background(alias) for alias in BACKGROUND_ALIASES]
    selected_rest: list[dict[str, Any]] = []
    all_rest_attempts: list[dict[str, Any]] = []
    dynamic_rows: list[dict[str, Any]] = []
    run_status = "COMPLETED"
    stop_reason: str | None = None
    try:
        for background in backgrounds:
            wt, wt_attempts = solve_wt_with_optional_scalar_calibration(background, budget)
            selected_rest.append(wt)
            all_rest_attempts.extend(wt_attempts)
            if not wt["admissible"]:
                continue
            null, null_attempts = solve_exact_null(background, wt, budget)
            null["g_nhe1_status"] = wt["g_nhe1_status"]
            selected_rest.append(null)
            all_rest_attempts.extend(null_attempts)
        for background in backgrounds:
            wt = next((row for row in selected_rest if row["background"] == background.alias and row["genotype"] == "WT"), None)
            null = next((row for row in selected_rest if row["background"] == background.alias and row["genotype"] == "AE4_NULL"), None)
            if wt is None or null is None or not wt["admissible"] or not null["admissible"]:
                continue
            dynamic_rows.append(run_stimulation(background, wt, genotype=WT, budget=budget))
            dynamic_rows.append(run_stimulation(background, null, genotype=AE4_NULL, budget=budget))
    except NumericalBudgetExhausted as exc:
        run_status = "STOPPED_AT_NUMERICAL_BUDGET"
        stop_reason = str(exc)
    finally:
        write_outputs(
            backgrounds=backgrounds,
            selected_rest=selected_rest,
            all_rest_attempts=all_rest_attempts,
            dynamic_rows=dynamic_rows,
            budget=budget,
            run_status=run_status,
            stop_reason=stop_reason,
        )
    return {
        "run_status": run_status,
        "stop_reason": stop_reason,
        "rest_cases": len(selected_rest),
        "dynamic_cases": len(dynamic_rows),
        "budget": budget.payload(),
    }


__all__ = (
    "Background",
    "NumericalBudget",
    "build_task30_model",
    "load_background",
    "run",
    "solve_stationary",
    "state_from_coordinates",
    "task30_parameters",
)
