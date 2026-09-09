"""WT-only calibration and resting-root geometry for the modern full model.

This module is intentionally sealed from every AE4-null secretion datum.  It
uses only the integrated model API, the matched WT resting chloride and pH
measurements, and explicitly labelled historical-lineage physiology gates.

The twelve core steady-state equations contain two exact charge dependencies.
Consequently this module does *not* call :meth:`ModernFullModel.solve_resting_root`.
It parameterizes cell and lumen electroneutrality, solves ten independent core
rows, verifies the omitted alkalinity rows afterwards, and reports the
singular values of the independent Jacobian.  This prevents a redundant-row
least-squares minimum from being mistaken for a root.

Provenance
----------
* Electroneutral coordinates, omitted-row identities, and the kinetic-scale
  gauge are ``DERIVED_CONSTRAINT``.
* WT Cl and pH targets are ``PRIMARY_MEASUREMENT`` (2015 JBC, matched AE4-line
  controls).
* Na/K, volume, resting calcium, and voltage gates are
  ``PUBLISHED_MODEL``/historical-lineage constraints, not new primary data.
* Bounds, deterministic multistart design, the common conductance profile, and
  each generation-specific minimal calibration choice are
  ``NEW_MODELING_DECISION``.

No knockout genotype is accepted by any public calibration function here.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares
from scipy.stats import qmc

from .acid_base import total_alkalinity_mM
from .model import ConstantStimulus, ModernFullModel, WT
from .parameters import FullModelParameters
from .states import WholeCellState, compartment_from_concentrations
from .transporters import AE4Parameters


CORE_INDEPENDENT_ROWS: tuple[int, ...] = (0, 1, 2, 3, 5, 6, 7, 8, 9, 11)
CORE_OMITTED_CHARGE_ROWS: tuple[int, ...] = (4, 10)
REST_COORDINATE_NAMES: tuple[str, ...] = (
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
CALIBRATED_PARAMETER_NAMES: tuple[str, ...] = (
    "ae4_carrier_amount_fmol_at_unit_rate_gauge",
    "nhe1_capacity_fmol_s",
    "common_membrane_conductance_scale",
    "cell_other_impermeant_osmoles_fmol",
)


@dataclass(frozen=True)
class BoundedValue:
    name: str
    lower: float
    preferred: float
    upper: float
    units: str
    tier: int
    provenance: str
    note: str

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (self.lower, self.preferred, self.upper)):
            raise ValueError(f"{self.name}: bounds and preferred value must be finite")
        if not self.lower < self.preferred < self.upper:
            raise ValueError(f"{self.name}: require lower < preferred < upper")
        if self.tier not in range(5):
            raise ValueError(f"{self.name}: tier must be 0, 1, 2, 3, or 4")


def _default_coordinate_bounds() -> tuple[BoundedValue, ...]:
    # Broad search domains are not measurement intervals.  The much tighter
    # post-solve gates below decide whether a root is physiologically usable.
    return (
        BoundedValue("na_i_mM", 4.0, 20.0, 60.0, "mM", 3, "PUBLISHED_MODEL", "historical-lineage rest target; broad search domain"),
        BoundedValue("k_i_mM", 60.0, 120.0, 210.0, "mM", 3, "PUBLISHED_MODEL", "historical-lineage rest target; broad search domain"),
        BoundedValue("tic_i_mM", 1.0, 20.0, 100.0, "mM", 3, "NEW_MODELING_DECISION", "unmeasured intracellular TIC search domain"),
        BoundedValue("ph_i", 6.2, 6.91, 8.0, "dimensionless (p-scale)", 2, "PRIMARY_MEASUREMENT", "search domain wider than matched-WT gate"),
        BoundedValue("volume_i_pL", 0.5, 1.3, 3.0, "pL", 3, "PUBLISHED_MODEL", "historical 1.3-pL convention; not primary uncertainty"),
        BoundedValue("na_l_mM", 20.0, 125.0, 280.0, "mM", 3, "NEW_MODELING_DECISION", "unmeasured finite-lumen rest domain"),
        BoundedValue("k_l_mM", 0.05, 5.0, 60.0, "mM", 3, "NEW_MODELING_DECISION", "unmeasured finite-lumen rest domain"),
        BoundedValue("tic_l_mM", 0.2, 25.0, 100.0, "mM", 3, "NEW_MODELING_DECISION", "unmeasured finite-lumen carbon domain"),
        BoundedValue("ph_l", 5.0, 7.4, 9.0, "dimensionless (p-scale)", 3, "NEW_MODELING_DECISION", "unmeasured local-lumen pH domain"),
        BoundedValue("volume_l_pL", 0.02, 0.10, 0.80, "pL", 3, "NEW_MODELING_DECISION", "finite local-lumen volume domain"),
    )


def _default_parameter_bounds() -> tuple[BoundedValue, ...]:
    return (
        BoundedValue(
            CALIBRATED_PARAMETER_NAMES[0], 1.0e-4, 5.0e-2, 5.0e-1, "fmol", 3,
            "NEW_MODELING_DECISION",
            "only carrier-times-attempt-rate is observable; all attempt rates use a declared 1/s gauge",
        ),
        BoundedValue(
            CALIBRATED_PARAMETER_NAMES[1], 1.0e-6, 5.0e-3, 2.0e-1, "fmol s^-1", 3,
            "NEW_MODELING_DECISION", "WT acid-base closure nuisance; not fitted to knockout pH",
        ),
        BoundedValue(
            CALIBRATED_PARAMETER_NAMES[2], 1.0, 15.0, 200.0, "dimensionless", 3,
            "NEW_MODELING_DECISION", "common scale on existing apical-Cl and total-K conductances",
        ),
        BoundedValue(
            CALIBRATED_PARAMETER_NAMES[3], 20.0, 110.0, 300.0, "fmol osmoles", 3,
            "NEW_MODELING_DECISION",
            "OTHER finite impermeant osmoles used only for WT volume closure; excludes the explicit 30-fmol buffer-site pool",
        ),
    )


@dataclass(frozen=True)
class WTCalibrationSpec:
    coordinate_bounds: tuple[BoundedValue, ...] = _default_coordinate_bounds()
    parameter_bounds: tuple[BoundedValue, ...] = _default_parameter_bounds()
    independent_rhs_scales: tuple[float, ...] = (
        5.0e-2, 5.0e-2, 5.0e-2, 2.0e-2, 1.0e-4,
        2.0e-2, 2.0e-2, 2.0e-2, 2.0e-2, 1.0e-4,
    )
    regulatory_rhs_scales_s_inv: tuple[float, ...] = ()
    wt_cl_target_mM: float = 50.10
    wt_cl_sem_mM: float = 1.50
    wt_ph_target: float = 6.91
    wt_ph_sem: float = 0.07
    volume_lineage_target_pL: float = 1.30
    mean_voltage_lineage_target_mV: float = -56.0
    resting_calcium_uM: float = 0.058
    root_scaled_tolerance: float = 1.0e-7
    omitted_row_raw_tolerance: float = 1.0e-9
    charge_tolerance_fmol: float = 1.0e-9
    current_tolerance_A: float = 1.0e-18
    boundary_relative_tolerance: float = 1.0e-5
    cluster_relative_tolerance: float = 2.0e-5
    jacobian_rank_relative_tolerance: float = 1.0e-8

    def __post_init__(self) -> None:
        if tuple(item.name for item in self.coordinate_bounds) != REST_COORDINATE_NAMES:
            raise ValueError("coordinate-bound names do not match the coordinate layout")
        if tuple(item.name for item in self.parameter_bounds) != CALIBRATED_PARAMETER_NAMES:
            raise ValueError("parameter-bound names do not match the calibration layout")
        if len(self.independent_rhs_scales) != len(CORE_INDEPENDENT_ROWS):
            raise ValueError("one positive scale is required for each independent core row")
        if min(self.independent_rhs_scales) <= 0.0:
            raise ValueError("RHS scales must be positive")
        if any(value <= 0.0 or not math.isfinite(value) for value in self.regulatory_rhs_scales_s_inv):
            raise ValueError("regulatory RHS scales must be finite and positive")

    @property
    def coordinate_lower(self) -> NDArray[np.float64]:
        return np.asarray([item.lower for item in self.coordinate_bounds], dtype=float)

    @property
    def coordinate_preferred(self) -> NDArray[np.float64]:
        return np.asarray([item.preferred for item in self.coordinate_bounds], dtype=float)

    @property
    def coordinate_upper(self) -> NDArray[np.float64]:
        return np.asarray([item.upper for item in self.coordinate_bounds], dtype=float)

    @property
    def log_parameter_lower(self) -> NDArray[np.float64]:
        return np.log([item.lower for item in self.parameter_bounds])

    @property
    def log_parameter_preferred(self) -> NDArray[np.float64]:
        return np.log([item.preferred for item in self.parameter_bounds])

    @property
    def log_parameter_upper(self) -> NDArray[np.float64]:
        return np.log([item.upper for item in self.parameter_bounds])


@dataclass(frozen=True)
class ModelVariant:
    variant_id: str
    mixed_bath_na_attempt_fraction: float
    apical_pump_fraction: float
    apical_k_fraction: float
    cooperative_gate: bool = False
    regulatory_family: str = "R0_STATIC_NEGATIVE_CONTROL"

    def __post_init__(self) -> None:
        for name in (
            "mixed_bath_na_attempt_fraction",
            "apical_pump_fraction",
            "apical_k_fraction",
        ):
            value = float(getattr(self, name))
            if not 0.0 < value < 1.0:
                raise ValueError(f"{name} must be strictly between zero and one")
        if self.cooperative_gate:
            raise ValueError("the WT calibration generation currently admits only the minimal SR2 QSS gate")


@dataclass(frozen=True)
class WTCalibrationParameters:
    ae4_carrier_amount_fmol_at_unit_rate_gauge: float
    nhe1_capacity_fmol_s: float
    common_membrane_conductance_scale: float
    cell_other_impermeant_osmoles_fmol: float

    @classmethod
    def from_log_vector(cls, vector: ArrayLike) -> "WTCalibrationParameters":
        values = np.exp(np.asarray(vector, dtype=float))
        if values.shape != (4,):
            raise ValueError("calibration parameter vector must have length four")
        return cls(*(float(value) for value in values))

    def as_log_vector(self) -> NDArray[np.float64]:
        return np.log(
            np.asarray(
                (
                    self.ae4_carrier_amount_fmol_at_unit_rate_gauge,
                    self.nhe1_capacity_fmol_s,
                    self.common_membrane_conductance_scale,
                    self.cell_other_impermeant_osmoles_fmol,
                ),
                dtype=float,
            )
        )


@dataclass(frozen=True)
class RootAttempt:
    start_id: str
    optimizer_success: bool
    admissible_state: bool
    converged_root: bool
    max_abs_scaled_independent_rhs: float
    max_abs_amount_rhs_fmol_s: float
    max_abs_volume_rhs_pL_s: float
    max_abs_regulatory_rhs_s_inv: float
    max_abs_omitted_charge_rhs_fmol_equivalent_s: float
    cost: float
    optimality: float
    nfev: int
    message: str
    boundary_hits: tuple[str, ...]
    coordinates: tuple[float, ...]
    log_parameters: tuple[float, ...] = ()


@dataclass(frozen=True)
class RootBranch:
    branch_id: str
    representative_start_id: str
    member_start_ids: tuple[str, ...]
    state: tuple[float, ...]
    coordinates: tuple[float, ...]
    parameters: Mapping[str, float]
    observables: Mapping[str, float]
    raw_rhs: tuple[float, ...]
    scaled_independent_rhs: tuple[float, ...]
    omitted_raw_rhs: tuple[float, ...]
    charge_residuals_fmol: Mapping[str, float]
    current_residuals_A: Mapping[str, float]
    conservation_residuals: Mapping[str, float]
    boundary_hits: tuple[str, ...]
    source_supported_boundary_hits: tuple[str, ...]
    disallowed_boundary_hits: tuple[str, ...]
    independent_jacobian_singular_values: tuple[float, ...]
    independent_jacobian_rank: int
    independent_jacobian_nullity: int
    passes_numerical_gate: bool
    passes_wt_gate: bool
    gate_failures: tuple[str, ...]


@dataclass(frozen=True)
class RootSearchReport:
    variant: ModelVariant
    mode: str
    attempts: tuple[RootAttempt, ...]
    branches: tuple[RootBranch, ...]
    start_design: Mapping[str, Any]
    structural_rank_note: str


def _raw_rhs_unit_maxima(raw_rhs: ArrayLike) -> tuple[float, float, float, float]:
    """Return unit-homogeneous raw-RHS diagnostics.

    The core amount rows contain Na, K, Cl, and TIC rates in fmol/s.  Total
    alkalinity rows are fmol-equivalent/s and are the two charge-dependent
    rows omitted from the numerical solve.  Volume and optional regulatory
    rows have their own units.  No maximum across these four groups is a
    physically meaningful norm.
    """

    raw = np.asarray(raw_rhs, dtype=float)
    if raw.ndim != 1 or len(raw) < 12:
        raise ValueError("raw RHS must contain the twelve-state conserved core")
    amount_rows = np.asarray((0, 1, 2, 3, 6, 7, 8, 9), dtype=int)
    volume_rows = np.asarray((5, 11), dtype=int)
    omitted_rows = np.asarray(CORE_OMITTED_CHARGE_ROWS, dtype=int)
    regulatory_max = (
        float(np.max(np.abs(raw[12:]))) if len(raw) > 12 else 0.0
    )
    return (
        float(np.max(np.abs(raw[amount_rows]))),
        float(np.max(np.abs(raw[volume_rows]))),
        regulatory_max,
        float(np.max(np.abs(raw[omitted_rows]))),
    )


def build_wt_model(
    variant: ModelVariant,
    calibration: WTCalibrationParameters,
    *,
    base_parameters: FullModelParameters | None = None,
    regulatory_model: Any | None = None,
) -> ModernFullModel:
    """Construct a WT-only resting model without changing bath chemistry."""

    base = base_parameters or FullModelParameters()
    scale = calibration.common_membrane_conductance_scale
    parameters = replace(
        base,
        homeostasis=replace(
            base.homeostasis,
            nhe1_capacity_fmol_s=calibration.nhe1_capacity_fmol_s,
        ),
        membranes=replace(
            base.membranes,
            apical_pump_fraction=variant.apical_pump_fraction,
            apical_k_fraction=variant.apical_k_fraction,
            g_k_total_S=base.membranes.g_k_total_S * scale,
            g_cl_apical_S=base.membranes.g_cl_apical_S * scale,
            g_basolateral_background_S=(
                base.membranes.g_basolateral_background_S * scale
            ),
        ),
        geometry=replace(
            base.geometry,
            cell_impermeant_osmoles_fmol=(
                calibration.cell_other_impermeant_osmoles_fmol
            ),
        ),
    )
    na_fraction = variant.mixed_bath_na_attempt_fraction
    ae4 = AE4Parameters(
        carrier_amount_fmol=(
            calibration.ae4_carrier_amount_fmol_at_unit_rate_gauge
        ),
        common_cl_attempt_rate_s=1.0,
        na_loaded_attempt_rate_s=2.0 * na_fraction,
        k_loaded_attempt_rate_s=2.0 * (1.0 - na_fraction),
    )
    return ModernFullModel(
        parameters,
        stimulus=ConstantStimulus(
            calcium_uM=0.058,
            beta_input=0.0,
        ),
        regulatory_model=regulatory_model,
        ae4_parameters=ae4,
    )


def encode_electroneutral_coordinates(
    model: ModernFullModel,
    coordinates: ArrayLike,
) -> NDArray[np.float64]:
    """Map ten physical coordinates to the exact bulk-charge manifold.

    Intracellular and luminal chloride are eliminated using the total-
    alkalinity charge identity.  The fixed intracellular anion amount is not
    varied here: it remains the registered parameter of the selected model.
    """

    x = np.asarray(coordinates, dtype=float)
    regulatory_size = len(model.layout.regulatory_names)
    if x.shape != (10 + regulatory_size,):
        raise ValueError("rest coordinate vector has the wrong size")
    (
        na_i,
        k_i,
        tic_i,
        ph_i,
        volume_i,
        na_l,
        k_l,
        tic_l,
        ph_l,
        volume_l,
    ) = (float(value) for value in x[:10])
    p = model.parameters
    ta_i = total_alkalinity_mM(
        ph=ph_i,
        total_carbon_mM=tic_i,
        buffer_total_mM=p.geometry.cell_buffer_total_fmol / volume_i,
        buffer_pka=p.acid_base.cell_buffer_pka,
        parameters=p.acid_base,
    )
    cl_i = (
        na_i
        + k_i
        - ta_i
        - p.geometry.fixed_cell_anion_equivalents_fmol / volume_i
    )
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
        regulation=tuple(float(value) for value in x[10:]),
    )
    return model.layout.encode(state)


def decode_rest_coordinates(
    model: ModernFullModel, state_vector: ArrayLike
) -> NDArray[np.float64]:
    state = model.layout.decode(state_vector)
    ci = state.cell.concentrations()
    lumen = state.lumen.concentrations()
    observables = model.evaluate(0.0, state_vector, genotype=WT).diagnostics.observables
    return np.asarray(
        (
            ci.na_mM,
            ci.k_mM,
            ci.tic_mM,
            observables.cell_acid_base.ph,
            state.cell.volume_pL,
            lumen.na_mM,
            lumen.k_mM,
            lumen.tic_mM,
            observables.lumen_acid_base.ph,
            state.lumen.volume_pL,
            *state.regulation,
        ),
        dtype=float,
    )


def _bounds_with_regulation(
    model: ModernFullModel, spec: WTCalibrationSpec
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    lower = spec.coordinate_lower.copy()
    preferred = spec.coordinate_preferred.copy()
    upper = spec.coordinate_upper.copy()
    if model.layout.regulatory_names:
        initial = np.asarray(model._initial_regulatory_state(), dtype=float)
        lower = np.r_[lower, np.zeros(len(initial))]
        preferred = np.r_[preferred, np.clip(initial, 1.0e-6, 1.0 - 1.0e-6)]
        upper = np.r_[upper, np.ones(len(initial))]
    return lower, preferred, upper


def _independent_scaled_rhs(
    model: ModernFullModel,
    coordinates: NDArray[np.float64],
    spec: WTCalibrationSpec,
) -> NDArray[np.float64]:
    state = encode_electroneutral_coordinates(model, coordinates)
    raw = model.rhs(0.0, state, genotype=WT)
    core = raw[np.asarray(CORE_INDEPENDENT_ROWS)] / np.asarray(
        spec.independent_rhs_scales
    )
    if model.layout.regulatory_names:
        if len(spec.regulatory_rhs_scales_s_inv) != len(model.layout.regulatory_names):
            raise ValueError(
                "one explicit numerical s^-1 scale is required per regulatory RHS; "
                "these scales condition the root residual and are not measured kinetics"
            )
        core = np.r_[
            core,
            raw[12:] / np.asarray(spec.regulatory_rhs_scales_s_inv, dtype=float),
        ]
    return core


def _latin_starts(
    lower: NDArray[np.float64],
    preferred: NDArray[np.float64],
    upper: NDArray[np.float64],
    *,
    count: int,
    seed: int,
) -> tuple[tuple[str, NDArray[np.float64]], ...]:
    if count < 3:
        raise ValueError("at least three non-collinear starts are required")
    starts: list[tuple[str, NDArray[np.float64]]] = [
        ("preferred", preferred.copy()),
        ("lower_quartile", lower + 0.25 * (upper - lower)),
        ("upper_quartile", lower + 0.75 * (upper - lower)),
    ]
    remaining = count - len(starts)
    if remaining:
        sampler = qmc.LatinHypercube(d=len(lower), scramble=True, seed=seed)
        unit = sampler.random(remaining)
        for index, point in enumerate(lower + unit * (upper - lower)):
            starts.append((f"latin_{index:03d}", point))
    return tuple(starts)


def _boundary_hits(
    vector: NDArray[np.float64],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    names: Sequence[str],
    tolerance: float,
) -> tuple[str, ...]:
    relative = np.minimum(vector - lower, upper - vector) / (upper - lower)
    return tuple(name for name, distance in zip(names, relative) if distance <= tolerance)


def _observable_map(model: ModernFullModel, state: NDArray[np.float64]) -> dict[str, float]:
    evaluation = model.evaluate(0.0, state, genotype=WT)
    obs = evaluation.diagnostics.observables
    mem = evaluation.diagnostics.membranes
    water = evaluation.diagnostics.water
    return {
        "na_i_mM": float(obs.cell_concentrations_mM["na"]),
        "k_i_mM": float(obs.cell_concentrations_mM["k"]),
        "cl_i_mM": float(obs.cell_concentrations_mM["cl"]),
        "tic_i_mM": float(obs.cell_concentrations_mM["tic"]),
        "hco3_i_mM": float(obs.cell_concentrations_mM["hco3"]),
        "ph_i": float(obs.cell_acid_base.ph),
        "volume_i_pL": float(state[5]),
        "na_l_mM": float(obs.lumen_concentrations_mM["na"]),
        "k_l_mM": float(obs.lumen_concentrations_mM["k"]),
        "cl_l_mM": float(obs.lumen_concentrations_mM["cl"]),
        "tic_l_mM": float(obs.lumen_concentrations_mM["tic"]),
        "ph_l": float(obs.lumen_acid_base.ph),
        "volume_l_pL": float(state[11]),
        "v_apical_mV": float(1.0e3 * mem.v_apical_V),
        "v_basolateral_mV": float(1.0e3 * mem.v_basolateral_V),
        "v_transepithelial_mV": float(1.0e3 * mem.v_transepithelial_V),
        "rest_outflow_pL_s": float(water.lumen_outflow_pL_s),
        "osm_cell_mOsm": float(obs.osmolarities_mOsm["cell"]),
        "osm_lumen_mOsm": float(obs.osmolarities_mOsm["lumen"]),
        "osm_bath_mOsm": float(obs.osmolarities_mOsm["bath"]),
    }


def _wt_gate_failures(
    observables: Mapping[str, float], spec: WTCalibrationSpec
) -> tuple[str, ...]:
    failures: list[str] = []

    def gate(name: str, lower: float, upper: float) -> None:
        value = observables[name]
        if not lower <= value <= upper:
            failures.append(f"{name}={value:.9g} outside [{lower:.9g},{upper:.9g}]")

    # Only Cl and pH are matched primary WT measurements.  The other bounds
    # are deliberately broad historical-lineage physiology screens.
    gate("cl_i_mM", spec.wt_cl_target_mM - 2 * spec.wt_cl_sem_mM, spec.wt_cl_target_mM + 2 * spec.wt_cl_sem_mM)
    gate("ph_i", spec.wt_ph_target - 2 * spec.wt_ph_sem, spec.wt_ph_target + 2 * spec.wt_ph_sem)
    gate("na_i_mM", 10.0, 35.0)
    gate("k_i_mM", 100.0, 160.0)
    gate("volume_i_pL", 1.17, 1.43)
    gate("v_apical_mV", -80.0, -30.0)
    gate("v_basolateral_mV", -90.0, -35.0)
    if min(observables[name] for name in ("na_l_mM", "k_l_mM", "cl_l_mM", "tic_l_mM", "volume_l_pL")) <= 0.0:
        failures.append("one or more luminal physical observables are nonpositive")
    return tuple(failures)


def _branch_from_fit(
    *,
    branch_id: str,
    representative_start_id: str,
    member_start_ids: Sequence[str],
    model: ModernFullModel,
    coordinates: NDArray[np.float64],
    parameters: Mapping[str, float],
    jacobian: NDArray[np.float64],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    bound_names: Sequence[str],
    spec: WTCalibrationSpec,
    boundary_vector: NDArray[np.float64] | None = None,
) -> RootBranch:
    state = encode_electroneutral_coordinates(model, coordinates)
    evaluation = model.evaluate(0.0, state, genotype=WT)
    raw = evaluation.rhs
    scaled = _independent_scaled_rhs(model, coordinates, spec)
    omitted = raw[np.asarray(CORE_OMITTED_CHARGE_ROWS)]
    # Residuals are dimensionless.  Scale every physical-coordinate column by
    # its declared bound span before the SVD so numerical rank cannot change
    # merely because mM is relabelled as M or pL as L.
    dimensionless_jacobian = np.asarray(jacobian, dtype=float) * (
        upper - lower
    )[None, :]
    singular = np.linalg.svd(dimensionless_jacobian, compute_uv=False)
    if len(singular):
        threshold = spec.jacobian_rank_relative_tolerance * singular[0]
        rank = int(np.sum(singular > threshold))
    else:
        rank = 0
    expected_rank = jacobian.shape[1]
    boundary_hits = _boundary_hits(
        (
            np.asarray(boundary_vector, dtype=float)
            if boundary_vector is not None
            else (
                np.r_[coordinates, [parameters[name] for name in CALIBRATED_PARAMETER_NAMES]]
                if parameters
                else coordinates
            )
        ),
        lower,
        upper,
        bound_names,
        spec.boundary_relative_tolerance,
    )
    charge = evaluation.diagnostics.state_charge_fmol
    currents = evaluation.diagnostics.membranes.current_residuals_A
    # A beta-free dynamic regulatory fraction has an exact, source-supported
    # equilibrium at zero.  That physical boundary is reported, but it is not
    # treated like an optimizer-driven capacity/chemistry bound hit.  No core
    # coordinate or calibrated parameter receives this exemption.
    regulatory_names = tuple(model.layout.regulatory_names)
    initial_regulatory = tuple(model._initial_regulatory_state())
    allowed_boundary_hits: list[str] = []
    for name in boundary_hits:
        if name not in regulatory_names:
            continue
        regulator_index = regulatory_names.index(name)
        bound_index = tuple(bound_names).index(name)
        initial_value = float(initial_regulatory[regulator_index])
        span = float(upper[bound_index] - lower[bound_index])
        initial_is_physical_boundary = (
            abs(initial_value - lower[bound_index]) <= spec.boundary_relative_tolerance * span
            or abs(initial_value - upper[bound_index]) <= spec.boundary_relative_tolerance * span
        )
        fitted_is_basal_equilibrium = (
            abs(coordinates[10 + regulator_index] - initial_value)
            <= spec.boundary_relative_tolerance * span
        )
        if initial_is_physical_boundary and fitted_is_basal_equilibrium:
            allowed_boundary_hits.append(name)
    disallowed_boundary_hits = tuple(
        name for name in boundary_hits if name not in allowed_boundary_hits
    )
    numerical = (
        float(np.max(np.abs(scaled))) <= spec.root_scaled_tolerance
        and float(np.max(np.abs(omitted))) <= spec.omitted_row_raw_tolerance
        and max(abs(float(value)) for value in charge.values()) <= spec.charge_tolerance_fmol
        and max(abs(float(value)) for value in currents.values()) <= spec.current_tolerance_A
        and rank == expected_rank
        and not disallowed_boundary_hits
    )
    observables = _observable_map(model, state)
    failures = list(_wt_gate_failures(observables, spec))
    if not numerical:
        failures.append("numerical/root/rank/boundary gate failed")
    return RootBranch(
        branch_id=branch_id,
        representative_start_id=representative_start_id,
        member_start_ids=tuple(member_start_ids),
        state=tuple(float(value) for value in state),
        coordinates=tuple(float(value) for value in coordinates),
        parameters=dict(parameters),
        observables=observables,
        raw_rhs=tuple(float(value) for value in raw),
        scaled_independent_rhs=tuple(float(value) for value in scaled),
        omitted_raw_rhs=tuple(float(value) for value in omitted),
        charge_residuals_fmol={key: float(value) for key, value in charge.items()},
        current_residuals_A={key: float(value) for key, value in currents.items()},
        conservation_residuals={
            key: float(value)
            for key, value in evaluation.diagnostics.conservation_residuals.items()
        },
        boundary_hits=boundary_hits,
        source_supported_boundary_hits=tuple(allowed_boundary_hits),
        disallowed_boundary_hits=disallowed_boundary_hits,
        independent_jacobian_singular_values=tuple(float(value) for value in singular),
        independent_jacobian_rank=rank,
        independent_jacobian_nullity=expected_rank - rank,
        passes_numerical_gate=numerical,
        passes_wt_gate=numerical and not failures,
        gate_failures=tuple(failures),
    )


def solve_wt_resting_branches(
    model: ModernFullModel,
    *,
    variant: ModelVariant,
    calibration: WTCalibrationParameters,
    spec: WTCalibrationSpec | None = None,
    start_count: int = 32,
    seed: int = 13021,
    max_nfev: int = 5000,
    reference_coordinates: ArrayLike | None = None,
) -> RootSearchReport:
    """Return every distinct fixed-parameter WT resting root found."""

    spec = spec or WTCalibrationSpec()
    lower, preferred, upper = _bounds_with_regulation(model, spec)
    starts = list(_latin_starts(lower, preferred, upper, count=start_count, seed=seed))
    if reference_coordinates is not None:
        reference = np.asarray(reference_coordinates, dtype=float)
        if reference.shape != lower.shape:
            raise ValueError("reference coordinate vector has the wrong size")
        if np.any(reference <= lower) or np.any(reference >= upper):
            raise ValueError("reference coordinates must be strictly interior")
        starts.insert(0, ("calibrated_reference", reference.copy()))
    scale_span = upper - lower
    attempts: list[RootAttempt] = []
    successful: list[tuple[str, Any]] = []

    def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            return _independent_scaled_rhs(model, candidate, spec)
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(len(candidate), 1.0e6)

    for start_id, start in starts:
        fit = least_squares(
            residual,
            start,
            bounds=(lower, upper),
            max_nfev=max_nfev,
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            x_scale="jac",
        )
        admissible = True
        try:
            state = encode_electroneutral_coordinates(model, fit.x)
            raw = model.rhs(0.0, state, genotype=WT)
            scaled_max = float(np.max(np.abs(residual(fit.x))))
            amount_max, volume_max, regulatory_max, omitted_max = (
                _raw_rhs_unit_maxima(raw)
            )
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            admissible = False
            scaled_max = amount_max = volume_max = regulatory_max = omitted_max = math.inf
        hits = _boundary_hits(
            fit.x,
            lower,
            upper,
            (*REST_COORDINATE_NAMES, *model.layout.regulatory_names),
            spec.boundary_relative_tolerance,
        )
        converged = bool(
            admissible
            and fit.success
            and scaled_max <= spec.root_scaled_tolerance
            and omitted_max <= spec.omitted_row_raw_tolerance
        )
        attempts.append(
            RootAttempt(
                start_id=start_id,
                optimizer_success=bool(fit.success),
                admissible_state=admissible,
                converged_root=converged,
                max_abs_scaled_independent_rhs=scaled_max,
                max_abs_amount_rhs_fmol_s=amount_max,
                max_abs_volume_rhs_pL_s=volume_max,
                max_abs_regulatory_rhs_s_inv=regulatory_max,
                max_abs_omitted_charge_rhs_fmol_equivalent_s=omitted_max,
                cost=float(fit.cost),
                optimality=float(fit.optimality),
                nfev=int(fit.nfev),
                message=str(fit.message),
                boundary_hits=hits,
                coordinates=tuple(float(value) for value in fit.x),
            )
        )
        if converged:
            successful.append((start_id, fit))

    clusters: list[list[tuple[str, Any]]] = []
    for item in successful:
        normalized = item[1].x / scale_span
        for cluster in clusters:
            representative = cluster[0][1].x / scale_span
            if float(np.max(np.abs(normalized - representative))) <= spec.cluster_relative_tolerance:
                cluster.append(item)
                break
        else:
            clusters.append([item])

    branches: list[RootBranch] = []
    parameter_map = asdict(calibration)
    for index, cluster in enumerate(clusters):
        representative = min(cluster, key=lambda item: np.max(np.abs(residual(item[1].x))))
        start_id, fit = representative
        branches.append(
            _branch_from_fit(
                branch_id=f"{variant.variant_id}:B{index:02d}",
                representative_start_id=start_id,
                member_start_ids=[item[0] for item in cluster],
                model=model,
                coordinates=np.asarray(fit.x, dtype=float),
                parameters={},
                jacobian=np.asarray(fit.jac, dtype=float),
                lower=lower,
                upper=upper,
                bound_names=(*REST_COORDINATE_NAMES, *model.layout.regulatory_names),
                spec=spec,
            )
        )
    return RootSearchReport(
        variant=variant,
        mode="fixed_parameter_root_geometry",
        attempts=tuple(attempts),
        branches=tuple(branches),
        start_design={
            "kind": "named_anchors_plus_scipy_LatinHypercube",
            "count": len(starts),
            "requested_design_count": start_count,
            "seed": seed,
            "calibrated_reference_added": reference_coordinates is not None,
            "fixed_calibration_parameters": asdict(calibration),
        },
        structural_rank_note=(
            "Ten independent core rows are solved on exact cell/lumen bulk-charge manifolds; "
            "alkalinity rows 4 and 10 are omitted as charge-dependent and then verified raw."
        ),
    )


def fit_wt_variant(
    variant: ModelVariant,
    *,
    spec: WTCalibrationSpec | None = None,
    base_parameters: FullModelParameters | None = None,
    regulatory_model: Any | None = None,
    start_count: int = 20,
    seed: int = 13022,
    max_nfev: int = 10000,
) -> RootSearchReport:
    """Fit exactly four WT-only closure parameters and retain all branches.

    The four added equations are matched WT chloride, matched WT pH, the
    historical-lineage 1.3-pL volume convention, and the mean of the two
    historical-lineage membrane-potential values.  Intracellular Na and K and
    each individual potential remain validation gates, not extra fit rows.
    """

    if variant.regulatory_family != "R0_STATIC_NEGATIVE_CONTROL" and regulatory_model is None:
        raise ValueError("a dynamic regulatory family requires its explicit model")
    spec = spec or WTCalibrationSpec()
    initial_calibration = WTCalibrationParameters.from_log_vector(
        spec.log_parameter_preferred
    )
    initial_model = build_wt_model(
        variant,
        initial_calibration,
        base_parameters=base_parameters,
        regulatory_model=regulatory_model,
    )
    coordinate_lower, coordinate_preferred, coordinate_upper = _bounds_with_regulation(
        initial_model, spec
    )
    lower = np.r_[coordinate_lower, spec.log_parameter_lower]
    preferred = np.r_[coordinate_preferred, spec.log_parameter_preferred]
    upper = np.r_[coordinate_upper, spec.log_parameter_upper]
    starts = _latin_starts(lower, preferred, upper, count=start_count, seed=seed)
    attempts: list[RootAttempt] = []
    successful: list[tuple[str, Any, ModernFullModel, WTCalibrationParameters]] = []
    coordinate_size = len(coordinate_lower)

    def unpack(candidate: NDArray[np.float64]) -> tuple[ModernFullModel, WTCalibrationParameters, NDArray[np.float64]]:
        calibration = WTCalibrationParameters.from_log_vector(candidate[coordinate_size:])
        model = build_wt_model(
            variant,
            calibration,
            base_parameters=base_parameters,
            regulatory_model=regulatory_model,
        )
        return model, calibration, candidate[:coordinate_size]

    def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            model, _, coordinates = unpack(candidate)
            state = encode_electroneutral_coordinates(model, coordinates)
            evaluation = model.evaluate(0.0, state, genotype=WT)
            obs = evaluation.diagnostics.observables
            mem = evaluation.diagnostics.membranes
            rest = _independent_scaled_rhs(model, coordinates, spec)
            mean_voltage_mV = 500.0 * (mem.v_apical_V + mem.v_basolateral_V)
            targets = np.asarray(
                (
                    (obs.cell_concentrations_mM["cl"] - spec.wt_cl_target_mM) / spec.wt_cl_sem_mM,
                    (obs.cell_acid_base.ph - spec.wt_ph_target) / spec.wt_ph_sem,
                    (state[5] - spec.volume_lineage_target_pL) / 0.13,
                    (mean_voltage_mV - spec.mean_voltage_lineage_target_mV) / 5.0,
                ),
                dtype=float,
            )
            return np.r_[rest, targets]
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(len(candidate), 1.0e6)

    all_names = (
        *REST_COORDINATE_NAMES,
        *initial_model.layout.regulatory_names,
        *CALIBRATED_PARAMETER_NAMES,
    )
    for start_id, start in starts:
        fit = least_squares(
            residual,
            start,
            bounds=(lower, upper),
            max_nfev=max_nfev,
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            x_scale="jac",
        )
        admissible = True
        try:
            model, calibration, coordinates = unpack(fit.x)
            state = encode_electroneutral_coordinates(model, coordinates)
            raw = model.rhs(0.0, state, genotype=WT)
            full_residual = residual(fit.x)
            scaled_max = float(np.max(np.abs(full_residual)))
            amount_max, volume_max, regulatory_max, omitted_max = (
                _raw_rhs_unit_maxima(raw)
            )
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            admissible = False
            model = initial_model
            calibration = initial_calibration
            coordinates = fit.x[:coordinate_size]
            scaled_max = amount_max = volume_max = regulatory_max = omitted_max = math.inf
        # Log-parameter bounds need a separate physical report, but active
        # detection is correctly performed in optimization coordinates here.
        hits = _boundary_hits(fit.x, lower, upper, all_names, spec.boundary_relative_tolerance)
        converged = bool(
            admissible
            and fit.success
            and scaled_max <= spec.root_scaled_tolerance
            and omitted_max <= spec.omitted_row_raw_tolerance
        )
        attempts.append(
            RootAttempt(
                start_id=start_id,
                optimizer_success=bool(fit.success),
                admissible_state=admissible,
                converged_root=converged,
                max_abs_scaled_independent_rhs=scaled_max,
                max_abs_amount_rhs_fmol_s=amount_max,
                max_abs_volume_rhs_pL_s=volume_max,
                max_abs_regulatory_rhs_s_inv=regulatory_max,
                max_abs_omitted_charge_rhs_fmol_equivalent_s=omitted_max,
                cost=float(fit.cost),
                optimality=float(fit.optimality),
                nfev=int(fit.nfev),
                message=str(fit.message),
                boundary_hits=hits,
                coordinates=tuple(float(value) for value in coordinates),
                log_parameters=tuple(float(value) for value in fit.x[coordinate_size:]),
            )
        )
        if converged:
            successful.append((start_id, fit, model, calibration))

    span = upper - lower
    clusters: list[list[tuple[str, Any, ModernFullModel, WTCalibrationParameters]]] = []
    for item in successful:
        normalized = item[1].x / span
        for cluster in clusters:
            representative = cluster[0][1].x / span
            if float(np.max(np.abs(normalized - representative))) <= spec.cluster_relative_tolerance:
                cluster.append(item)
                break
        else:
            clusters.append([item])

    branches: list[RootBranch] = []
    for index, cluster in enumerate(clusters):
        representative = min(cluster, key=lambda item: np.max(np.abs(residual(item[1].x))))
        start_id, fit, model, calibration = representative
        physical_parameters = asdict(calibration)
        # The branch helper expects physical bounds when checking bound hits.
        physical_lower = np.r_[coordinate_lower, np.exp(spec.log_parameter_lower)]
        physical_upper = np.r_[coordinate_upper, np.exp(spec.log_parameter_upper)]
        physical_vector = np.r_[fit.x[:coordinate_size], np.exp(fit.x[coordinate_size:])]
        physical_jacobian = np.asarray(fit.jac, dtype=float).copy()
        # Chain rule: dR/dp = dR/dlog(p) / p.
        physical_jacobian[:, coordinate_size:] /= np.exp(fit.x[coordinate_size:])[None, :]
        branches.append(
            _branch_from_fit(
                branch_id=f"{variant.variant_id}:C{index:02d}",
                representative_start_id=start_id,
                member_start_ids=[item[0] for item in cluster],
                model=model,
                coordinates=np.asarray(fit.x[:coordinate_size], dtype=float),
                parameters=physical_parameters,
                jacobian=physical_jacobian,
                lower=physical_lower,
                upper=physical_upper,
                bound_names=all_names,
                spec=spec,
            )
        )
        # Assert helper and directly computed physical vector stay aligned.
        if _boundary_hits(
            physical_vector,
            physical_lower,
            physical_upper,
            all_names,
            spec.boundary_relative_tolerance,
        ) != branches[-1].boundary_hits:
            raise AssertionError("physical boundary report is inconsistent")

    return RootSearchReport(
        variant=variant,
        mode="four_parameter_WT_only_calibration",
        attempts=tuple(attempts),
        branches=tuple(branches),
        start_design={"kind": "named_anchors_plus_scipy_LatinHypercube", "count": start_count, "seed": seed},
        structural_rank_note=(
            "Ten independent charge-manifold balances plus four WT-only calibration rows form a square system; "
            "the two omitted charge-dependent alkalinity rows are verified after fitting."
        ),
    )


def fit_wt_source_compatible_variant(
    variant: ModelVariant,
    *,
    spec: WTCalibrationSpec | None = None,
    base_parameters: FullModelParameters | None = None,
    regulatory_model: Any | None = None,
    fixed_conductance_scale: float = 10.0,
    fixed_other_impermeant_start_fmol: float = 110.0,
    start_count: int = 20,
    seed: int = 13023,
    max_nfev: int = 10000,
) -> RootSearchReport:
    """Minimal exact WT calibration used after the four-row attempt failed.

    Three positive quantities are fitted: the effective AE4 capacity gauge,
    NHE1 capacity, and the OTHER impermeant osmole pool.  The equations added
    to the ten independent steady balances are matched-WT chloride, matched-WT
    pH, and the historical-lineage 1.3-pL volume convention.  The total-K and
    apical-Cl conductance scale is fixed before this fit and must independently
    pass both resting-voltage gates.  It is not identified by this solve.

    This hierarchy is preferable to the rejected four-parameter system because
    it obtains an actual root rather than compromising steady balance against
    a voltage residual.  Equal and Na-biased routing variants are retained as
    failures if they require the nonnegative NHE1 capacity to hit its bound.
    """

    if not math.isfinite(fixed_conductance_scale) or fixed_conductance_scale <= 0.0:
        raise ValueError("fixed conductance scale must be finite and positive")
    if variant.regulatory_family != "R0_STATIC_NEGATIVE_CONTROL" and regulatory_model is None:
        raise ValueError("a dynamic regulatory family requires its explicit model")
    spec = spec or WTCalibrationSpec()
    parameter_indices = np.asarray((0, 1, 3), dtype=int)
    log_parameter_lower = spec.log_parameter_lower[parameter_indices]
    log_parameter_upper = spec.log_parameter_upper[parameter_indices]
    log_parameter_preferred = np.log(
        (
            spec.parameter_bounds[0].preferred,
            spec.parameter_bounds[1].preferred,
            fixed_other_impermeant_start_fmol,
        )
    )
    initial_calibration = WTCalibrationParameters(
        ae4_carrier_amount_fmol_at_unit_rate_gauge=spec.parameter_bounds[0].preferred,
        nhe1_capacity_fmol_s=spec.parameter_bounds[1].preferred,
        common_membrane_conductance_scale=fixed_conductance_scale,
        cell_other_impermeant_osmoles_fmol=fixed_other_impermeant_start_fmol,
    )
    initial_model = build_wt_model(
        variant,
        initial_calibration,
        base_parameters=base_parameters,
        regulatory_model=regulatory_model,
    )
    coordinate_lower, coordinate_preferred, coordinate_upper = _bounds_with_regulation(
        initial_model, spec
    )
    lower = np.r_[coordinate_lower, log_parameter_lower]
    preferred = np.r_[coordinate_preferred, log_parameter_preferred]
    upper = np.r_[coordinate_upper, log_parameter_upper]
    starts = _latin_starts(lower, preferred, upper, count=start_count, seed=seed)
    coordinate_size = len(coordinate_lower)
    attempts: list[RootAttempt] = []
    successful: list[tuple[str, Any, ModernFullModel, WTCalibrationParameters]] = []

    def unpack(candidate: NDArray[np.float64]) -> tuple[ModernFullModel, WTCalibrationParameters, NDArray[np.float64]]:
        carrier, nhe1, other_osmoles = np.exp(candidate[coordinate_size:])
        calibration = WTCalibrationParameters(
            ae4_carrier_amount_fmol_at_unit_rate_gauge=float(carrier),
            nhe1_capacity_fmol_s=float(nhe1),
            common_membrane_conductance_scale=fixed_conductance_scale,
            cell_other_impermeant_osmoles_fmol=float(other_osmoles),
        )
        model = build_wt_model(
            variant,
            calibration,
            base_parameters=base_parameters,
            regulatory_model=regulatory_model,
        )
        return model, calibration, candidate[:coordinate_size]

    def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            model, _, coordinates = unpack(candidate)
            state = encode_electroneutral_coordinates(model, coordinates)
            evaluation = model.evaluate(0.0, state, genotype=WT)
            obs = evaluation.diagnostics.observables
            return np.r_[
                _independent_scaled_rhs(model, coordinates, spec),
                (obs.cell_concentrations_mM["cl"] - spec.wt_cl_target_mM)
                / spec.wt_cl_sem_mM,
                (obs.cell_acid_base.ph - spec.wt_ph_target) / spec.wt_ph_sem,
                (state[5] - spec.volume_lineage_target_pL) / 0.13,
            ]
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(len(candidate), 1.0e6)

    fitted_names = (
        CALIBRATED_PARAMETER_NAMES[0],
        CALIBRATED_PARAMETER_NAMES[1],
        CALIBRATED_PARAMETER_NAMES[3],
    )
    all_bound_names = (
        *REST_COORDINATE_NAMES,
        *initial_model.layout.regulatory_names,
        *fitted_names,
    )
    for start_id, start in starts:
        fit = least_squares(
            residual,
            start,
            bounds=(lower, upper),
            max_nfev=max_nfev,
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            x_scale="jac",
        )
        admissible = True
        try:
            model, calibration, coordinates = unpack(fit.x)
            state = encode_electroneutral_coordinates(model, coordinates)
            raw = model.rhs(0.0, state, genotype=WT)
            full_residual = residual(fit.x)
            scaled_max = float(np.max(np.abs(full_residual)))
            amount_max, volume_max, regulatory_max, omitted_max = (
                _raw_rhs_unit_maxima(raw)
            )
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            admissible = False
            model = initial_model
            calibration = initial_calibration
            coordinates = fit.x[:coordinate_size]
            scaled_max = amount_max = volume_max = regulatory_max = omitted_max = math.inf
        hits = _boundary_hits(
            fit.x,
            lower,
            upper,
            all_bound_names,
            spec.boundary_relative_tolerance,
        )
        converged = bool(
            admissible
            and fit.success
            and scaled_max <= spec.root_scaled_tolerance
            and omitted_max <= spec.omitted_row_raw_tolerance
        )
        attempts.append(
            RootAttempt(
                start_id=start_id,
                optimizer_success=bool(fit.success),
                admissible_state=admissible,
                converged_root=converged,
                max_abs_scaled_independent_rhs=scaled_max,
                max_abs_amount_rhs_fmol_s=amount_max,
                max_abs_volume_rhs_pL_s=volume_max,
                max_abs_regulatory_rhs_s_inv=regulatory_max,
                max_abs_omitted_charge_rhs_fmol_equivalent_s=omitted_max,
                cost=float(fit.cost),
                optimality=float(fit.optimality),
                nfev=int(fit.nfev),
                message=str(fit.message),
                boundary_hits=hits,
                coordinates=tuple(float(value) for value in coordinates),
                log_parameters=tuple(float(value) for value in fit.x[coordinate_size:]),
            )
        )
        if converged:
            successful.append((start_id, fit, model, calibration))

    span = upper - lower
    clusters: list[list[tuple[str, Any, ModernFullModel, WTCalibrationParameters]]] = []
    for item in successful:
        normalized = item[1].x / span
        for cluster in clusters:
            representative = cluster[0][1].x / span
            if float(np.max(np.abs(normalized - representative))) <= spec.cluster_relative_tolerance:
                cluster.append(item)
                break
        else:
            clusters.append([item])

    branches: list[RootBranch] = []
    physical_lower = np.r_[coordinate_lower, np.exp(log_parameter_lower)]
    physical_upper = np.r_[coordinate_upper, np.exp(log_parameter_upper)]
    for index, cluster in enumerate(clusters):
        representative = min(cluster, key=lambda item: np.max(np.abs(residual(item[1].x))))
        start_id, fit, model, calibration = representative
        fitted_physical = np.exp(fit.x[coordinate_size:])
        physical_vector = np.r_[fit.x[:coordinate_size], fitted_physical]
        physical_jacobian = np.asarray(fit.jac, dtype=float).copy()
        physical_jacobian[:, coordinate_size:] /= fitted_physical[None, :]
        branches.append(
            _branch_from_fit(
                branch_id=f"{variant.variant_id}:M{index:02d}",
                representative_start_id=start_id,
                member_start_ids=[item[0] for item in cluster],
                model=model,
                coordinates=np.asarray(fit.x[:coordinate_size], dtype=float),
                parameters=asdict(calibration),
                jacobian=physical_jacobian,
                lower=physical_lower,
                upper=physical_upper,
                bound_names=all_bound_names,
                spec=spec,
                boundary_vector=physical_vector,
            )
        )
    return RootSearchReport(
        variant=variant,
        mode="minimal_three_parameter_WT_only_calibration",
        attempts=tuple(attempts),
        branches=tuple(branches),
        start_design={
            "kind": "named_anchors_plus_scipy_LatinHypercube",
            "count": start_count,
            "seed": seed,
            "fixed_conductance_scale": fixed_conductance_scale,
            "fixed_conductance_status": "Tier-3 gate value, not identified by the fit",
        },
        structural_rank_note=(
            "Ten independent charge-manifold balances plus WT Cl, WT pH, and historical-lineage volume form a square system. "
            "Na/K and both membrane potentials are out-of-objective gates."
        ),
    )


def default_variant_panel() -> tuple[ModelVariant, ...]:
    """Small deterministic panel spanning unresolved routing and topology."""

    variants: list[ModelVariant] = []
    for topology_name, pump_fraction, k_fraction in (
        ("LOW_APICAL", 0.03, 0.03),
        ("BALANCED_APICAL", 0.10, 0.10),
        ("ASYMMETRIC_APICAL", 0.10, 0.30),
    ):
        for routing_name, na_fraction in (
            ("K_HEAVY", 0.05),
            ("K_BIASED", 0.10),
            ("MODERATE_NA", 0.20),
            ("EQUAL", 0.50),
            ("NA_BIASED", 0.90),
        ):
            variants.append(
                ModelVariant(
                    variant_id=f"G2_{topology_name}_{routing_name}",
                    mixed_bath_na_attempt_fraction=na_fraction,
                    apical_pump_fraction=pump_fraction,
                    apical_k_fraction=k_fraction,
                )
            )
    return tuple(variants)


def profile_wt_parameter_ablations(
    variant: ModelVariant,
    reference_calibration: WTCalibrationParameters,
    reference_coordinates: ArrayLike,
    *,
    spec: WTCalibrationSpec | None = None,
    start_count: int = 3,
    seed: int = 13201,
    max_nfev: int = 5000,
) -> tuple[Mapping[str, Any], ...]:
    """Test whether each member of the three-parameter WT fit can be removed.

    One fitted quantity at a time is fixed to its declared *pre-fit preferred*
    value while the other two and all ten charge-manifold coordinates are
    reoptimized.  This is a generation-local minimality test: failure means
    that the selected 13-equation closure cannot absorb that removal using the
    remaining two quantities.  It is not a global identifiability claim and
    does not exclude compensation by an additional unsupported parameter.
    """

    if variant.regulatory_family != "R0_STATIC_NEGATIVE_CONTROL":
        raise ValueError("parameter ablations are defined on the R0 WT-rest scaffold")
    spec = spec or WTCalibrationSpec()
    coordinates_reference = np.asarray(reference_coordinates, dtype=float)
    if coordinates_reference.shape != (10,):
        raise ValueError("R0 reference coordinates must have length ten")
    fitted_names = (
        CALIBRATED_PARAMETER_NAMES[0],
        CALIBRATED_PARAMETER_NAMES[1],
        CALIBRATED_PARAMETER_NAMES[3],
    )
    parameter_index = {
        name: CALIBRATED_PARAMETER_NAMES.index(name) for name in fitted_names
    }
    reference_values = asdict(reference_calibration)
    results: list[Mapping[str, Any]] = []

    for ablation_index, removed_name in enumerate(fitted_names):
        active_names = tuple(name for name in fitted_names if name != removed_name)
        active_indices = np.asarray(
            [parameter_index[name] for name in active_names], dtype=int
        )
        lower = np.r_[
            spec.coordinate_lower,
            spec.log_parameter_lower[active_indices],
        ]
        preferred = np.r_[
            coordinates_reference,
            np.log([reference_values[name] for name in active_names]),
        ]
        upper = np.r_[
            spec.coordinate_upper,
            spec.log_parameter_upper[active_indices],
        ]
        starts = _latin_starts(
            lower,
            preferred,
            upper,
            count=start_count,
            seed=seed + ablation_index,
        )
        removed_pre_fit_value = spec.parameter_bounds[
            parameter_index[removed_name]
        ].preferred

        def unpack(
            candidate: NDArray[np.float64],
        ) -> tuple[ModernFullModel, WTCalibrationParameters, NDArray[np.float64]]:
            values = dict(reference_values)
            values[removed_name] = float(removed_pre_fit_value)
            for name, value in zip(active_names, np.exp(candidate[10:])):
                values[name] = float(value)
            calibration = WTCalibrationParameters(**values)
            return (
                build_wt_model(variant, calibration),
                calibration,
                candidate[:10],
            )

        def residual(candidate: NDArray[np.float64]) -> NDArray[np.float64]:
            try:
                model, _, coordinates = unpack(candidate)
                state = encode_electroneutral_coordinates(model, coordinates)
                evaluation = model.evaluate(0.0, state, genotype=WT)
                observables = evaluation.diagnostics.observables
                return np.r_[
                    _independent_scaled_rhs(model, coordinates, spec),
                    (
                        observables.cell_concentrations_mM["cl"]
                        - spec.wt_cl_target_mM
                    )
                    / spec.wt_cl_sem_mM,
                    (observables.cell_acid_base.ph - spec.wt_ph_target)
                    / spec.wt_ph_sem,
                    (state[5] - spec.volume_lineage_target_pL) / 0.13,
                ]
            except (
                ValueError,
                FloatingPointError,
                OverflowError,
                np.linalg.LinAlgError,
            ):
                return np.full(13, 1.0e6)

        fits: list[tuple[str, Any]] = []
        for start_id, start in starts:
            fit = least_squares(
                residual,
                start,
                bounds=(lower, upper),
                max_nfev=max_nfev,
                xtol=1.0e-12,
                ftol=1.0e-12,
                gtol=1.0e-12,
                x_scale="jac",
            )
            fits.append((start_id, fit))
        start_id, best = min(
            fits, key=lambda item: float(np.max(np.abs(residual(item[1].x))))
        )
        model, calibration, coordinates = unpack(best.x)
        state = encode_electroneutral_coordinates(model, coordinates)
        evaluation = model.evaluate(0.0, state, genotype=WT)
        dimensionless_residual = residual(best.x)
        normalized_jacobian = np.asarray(best.jac, dtype=float) * (
            upper - lower
        )[None, :]
        singular = np.linalg.svd(normalized_jacobian, compute_uv=False)
        threshold = (
            spec.jacobian_rank_relative_tolerance * singular[0]
            if len(singular) else math.inf
        )
        rank = int(np.sum(singular > threshold))
        hits = _boundary_hits(
            best.x,
            lower,
            upper,
            (*REST_COORDINATE_NAMES, *active_names),
            spec.boundary_relative_tolerance,
        )
        observables = evaluation.diagnostics.observables
        results.append(
            {
                "ablation_id": f"FIX_{removed_name}_AT_PREFIT_PREFERRED",
                "removed_fit_parameter": removed_name,
                "held_fixed_value": float(removed_pre_fit_value),
                "held_fixed_unit": spec.parameter_bounds[
                    parameter_index[removed_name]
                ].units,
                "active_fitted_parameters": list(active_names),
                "start_design": {
                    "kind": "reference_anchor_plus_two_noncollinear_quartile_starts",
                    "count": start_count,
                    "seed": seed + ablation_index,
                },
                "representative_start_id": start_id,
                "optimizer_success": bool(best.success),
                "nfev": int(best.nfev),
                "cost": float(best.cost),
                "max_abs_dimensionless_residual": float(
                    np.max(np.abs(dimensionless_residual))
                ),
                "dimensionless_residual_vector": [
                    float(value) for value in dimensionless_residual
                ],
                "exact_root": bool(
                    best.success
                    and np.max(np.abs(dimensionless_residual))
                    <= spec.root_scaled_tolerance
                    and not hits
                ),
                "boundary_hits": list(hits),
                "normalized_jacobian_rank": rank,
                "normalized_jacobian_column_count": normalized_jacobian.shape[1],
                "normalized_jacobian_singular_values": [
                    float(value) for value in singular
                ],
                "fitted_parameter_values": {
                    name: float(getattr(calibration, name)) for name in active_names
                },
                "coordinates": [float(value) for value in coordinates],
                "target_residuals_native": {
                    "cl_i_mM": float(
                        observables.cell_concentrations_mM["cl"]
                        - spec.wt_cl_target_mM
                    ),
                    "ph_i_dimensionless_p_scale": float(
                        observables.cell_acid_base.ph - spec.wt_ph_target
                    ),
                    "volume_i_pL": float(
                        state[5] - spec.volume_lineage_target_pL
                    ),
                },
                "interpretation_scope": (
                    "generation-local necessity relative to the declared pre-fit preferred value; "
                    "not a structural or global identifiability result"
                ),
            }
        )
    return tuple(results)


def report_to_dict(report: RootSearchReport) -> dict[str, Any]:
    return {
        "variant": asdict(report.variant),
        "mode": report.mode,
        "attempts": [asdict(item) for item in report.attempts],
        "branches": [asdict(item) for item in report.branches],
        "start_design": dict(report.start_design),
        "structural_rank_note": report.structural_rank_note,
    }


__all__ = (
    "CALIBRATED_PARAMETER_NAMES",
    "CORE_INDEPENDENT_ROWS",
    "CORE_OMITTED_CHARGE_ROWS",
    "REST_COORDINATE_NAMES",
    "BoundedValue",
    "ModelVariant",
    "RootAttempt",
    "RootBranch",
    "RootSearchReport",
    "WTCalibrationParameters",
    "WTCalibrationSpec",
    "build_wt_model",
    "decode_rest_coordinates",
    "default_variant_panel",
    "encode_electroneutral_coordinates",
    "fit_wt_variant",
    "fit_wt_source_compatible_variant",
    "profile_wt_parameter_ablations",
    "report_to_dict",
    "solve_wt_resting_branches",
)
