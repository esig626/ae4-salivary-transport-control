"""Independent transcription of the Task 13B native-source panel.

This module deliberately does not import :mod:`native_source_panel` or call a
production model builder.  It translates a compact, frozen panel definition
into the immutable parameter records consumed by the independently
transcribed equations in :mod:`modern_full_model.independent`.

No held-out phenotype target or genotype-dependent calibration appears here.
The three source classes are exact nested alternatives at source scale one:

``STRICT_NATIVE``
    Source-literal membrane conductances and the frozen WT-only transporter
    construction, with no absolute-source multiplier.

``N_ABS_NKCC``
    One scalar multiplies the complete reversible NKCC1 capacity and nothing
    else.  Stoichiometry, reversal, and any subsequently composed WT stimulus
    multiplier are unchanged.

``AN_ABS_AE4_NHE``
    One scalar multiplies the frozen AE4 carrier amount and NHE1 capacity
    together.  Their ratio, the AE4 branch routing, and thermodynamics are
    unchanged.

Amounts are fmol, time is seconds, conductances are S, and hydraulic
permeabilities are pL s^-1 (mOsm/L)^-1.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .independent import (
    IndependentAE4Kinetics,
    IndependentConstantStimulus,
    IndependentWholeCell,
    free_charge_coordinates,
    independent_total_alkalinity_mM,
    state_from_free_charge_coordinates,
)
from .parameters import FullModelParameters


STRICT_NATIVE = "STRICT_NATIVE"
N_ABS_NKCC = "N_ABS_NKCC"
AN_ABS_AE4_NHE = "AN_ABS_AE4_NHE"
SOURCE_CLASSES = (STRICT_NATIVE, N_ABS_NKCC, AN_ABS_AE4_NHE)
NATIVE_ROOT_ROWS = (0, 1, 2, 3, 5, 6, 7, 8, 9, 11)

# Literal source map frozen before the native reroot.  These are maximum
# whole-cell conductances; calcium gating is applied by the membrane equation.
NATIVE_G_CL_APICAL_S = 31.4e-9
NATIVE_G_K_TOTAL_S = 14.0e-9
NATIVE_CA_HALF_UM = 0.26
NATIVE_CA_HILL = 1.46

# Frozen WT-only construction inherited from the pre-native generation.  The
# source panel changes absolute capacities only through the declared classes.
FROZEN_AE4_CARRIER_FMOL = 0.1501264265070657
FROZEN_NHE1_CAPACITY_FMOL_S = 0.007440980590414482
FROZEN_NKCC1_CAPACITY_FMOL_S = 0.08

# Local transcription of the declared WT resting-physiology gate.  Keeping
# these bounds here makes the independent reproduction genuinely independent
# of both the production calibration helper and the production native panel.
INDEPENDENT_WT_BOUNDS = {
    "cl_i_mM": (47.10, 53.10),
    "ph_i": (6.77, 7.05),
    "na_i_mM": (10.0, 35.0),
    "k_i_mM": (100.0, 160.0),
    "volume_i_pL": (1.17, 1.43),
    "v_apical_mV": (-80.0, -30.0),
    "v_basolateral_mV": (-90.0, -35.0),
}
INDEPENDENT_WT_POSITIVE_OBSERVABLES = (
    "na_l_mM",
    "k_l_mM",
    "cl_l_mM",
    "tic_l_mM",
    "volume_l_pL",
)


def independent_native_wt_observables(
    evaluation: Any,
    state: ArrayLike,
) -> Mapping[str, float]:
    """Derive every WT-gate observable from the independent evaluation.

    Concentrations are also calculated directly from the conserved amounts,
    rather than copied from the production root-table ``observables_json``.
    The latter can therefore be compared with this mapping as an additional
    transcription/integrity check by the independent runner.
    """

    values = np.asarray(state, dtype=float)
    if values.shape != (12,) or np.any(~np.isfinite(values)):
        raise ValueError("independent WT observables require a finite length-12 state")
    volume_i = float(values[5])
    volume_l = float(values[11])
    if volume_i <= 0.0 or volume_l <= 0.0:
        raise ValueError("independent WT observable volumes must be positive")
    observables = {
        "na_i_mM": float(values[0] / volume_i),
        "k_i_mM": float(values[1] / volume_i),
        "cl_i_mM": float(values[2] / volume_i),
        "tic_i_mM": float(values[3] / volume_i),
        "ph_i": float(evaluation.cell_speciation.ph),
        "volume_i_pL": volume_i,
        "na_l_mM": float(values[6] / volume_l),
        "k_l_mM": float(values[7] / volume_l),
        "cl_l_mM": float(values[8] / volume_l),
        "tic_l_mM": float(values[9] / volume_l),
        "ph_l": float(evaluation.lumen_speciation.ph),
        "volume_l_pL": volume_l,
        "v_apical_mV": 1.0e3 * float(evaluation.voltages_V["apical"]),
        "v_basolateral_mV": 1.0e3 * float(evaluation.voltages_V["basolateral"]),
    }
    if any(not math.isfinite(value) for value in observables.values()):
        raise ValueError("independent WT observables must all be finite")
    return observables


def independent_native_wt_gate_failures(
    observables: Mapping[str, Any],
) -> tuple[str, ...]:
    """Apply the local WT gate, returning deterministic failure labels."""

    failures: list[str] = []
    for key, (lower, upper) in INDEPENDENT_WT_BOUNDS.items():
        try:
            value = float(observables[key])
        except (KeyError, TypeError, ValueError):
            failures.append(f"{key}:missing_or_invalid")
            continue
        if not math.isfinite(value) or not lower <= value <= upper:
            failures.append(f"{key}:outside_[{lower:g},{upper:g}]")
    for key in INDEPENDENT_WT_POSITIVE_OBSERVABLES:
        try:
            value = float(observables[key])
        except (KeyError, TypeError, ValueError):
            failures.append(f"{key}:missing_or_invalid")
            continue
        if not math.isfinite(value) or value <= 0.0:
            failures.append(f"{key}:not_positive")
    return tuple(failures)


def independent_hydraulic_scaled_water_parameters(
    base: Any, *, scale: float, mode: str
) -> Any:
    """Apply the independent HW/HWQ transform to water parameters.

    This low-level transform admits unit scale so exact nesting can be tested.
    Public panel definitions below enforce the stricter ``H=1 <=> H1``
    member invariant.
    """

    factor = float(scale)
    if not math.isfinite(factor) or factor <= 0.0:
        raise ValueError("independent hydraulic scale must be positive")
    if mode not in {"H1", "HW", "HWQ"}:
        raise ValueError("independent hydraulic mode must be H1, HW, or HWQ")
    if mode == "H1" and factor != 1.0:
        raise ValueError("the independent H1 transform requires scale one")
    return replace(
        base,
        apical_hydraulic_pL_s_mOsm=(
            base.apical_hydraulic_pL_s_mOsm * factor
        ),
        basolateral_hydraulic_pL_s_mOsm=(
            base.basolateral_hydraulic_pL_s_mOsm * factor
        ),
        paracellular_hydraulic_pL_s_mOsm=(
            base.paracellular_hydraulic_pL_s_mOsm * factor
        ),
        outflow_rate_s=(
            base.outflow_rate_s * factor if mode == "HWQ" else base.outflow_rate_s
        ),
    )


@dataclass(frozen=True)
class IndependentNativeRootResult:
    """Joint state/OTHER reroot from the independent equation transcription."""

    state: NDArray[np.float64]
    other_impermeant_osmoles_fmol: float
    success: bool
    max_abs_scaled_residual: float
    max_abs_amount_source_fmol_s: float
    max_abs_volume_source_pL_s: float
    nfev: int
    cost: float
    scaled_jacobian_singular_values: tuple[float, ...]
    jacobian_rank: int
    jacobian_nullity: int
    message: str


@dataclass(frozen=True)
class IndependentNativeSourceDefinition:
    """Plain definition of one independently reconstructed native member."""

    source_class: str
    source_scale: float
    pump_capacity_scale: float
    apical_pump_fraction: float
    apical_k_fraction: float
    ae4_cation_fraction: float
    hydraulic_scale: float
    cell_other_impermeant_osmoles_fmol: float
    hydraulic_mode: str = "H1"

    def __post_init__(self) -> None:
        if self.source_class not in SOURCE_CLASSES:
            raise ValueError(f"unknown native source class: {self.source_class!r}")
        positive = (
            self.source_scale,
            self.pump_capacity_scale,
            self.hydraulic_scale,
            self.cell_other_impermeant_osmoles_fmol,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in positive):
            raise ValueError("native source scales and OTHER osmoles must be positive")
        fractions = (
            self.apical_pump_fraction,
            self.apical_k_fraction,
            self.ae4_cation_fraction,
        )
        if not all(math.isfinite(value) and 0.0 < value < 1.0 for value in fractions):
            raise ValueError("native routing fractions must lie strictly in (0, 1)")
        if self.source_class == STRICT_NATIVE and self.source_scale != 1.0:
            raise ValueError("STRICT_NATIVE is defined only at exact source scale one")
        if self.hydraulic_mode not in {"H1", "HW", "HWQ"}:
            raise ValueError("independent hydraulic mode must be H1, HW, or HWQ")
        if self.hydraulic_scale == 1.0 and self.hydraulic_mode != "H1":
            raise ValueError("unit independent hydraulic scale must declare H1")
        if self.hydraulic_scale != 1.0 and self.hydraulic_mode not in {"HW", "HWQ"}:
            raise ValueError("nonunit independent hydraulic scale must declare HW or HWQ")

    @classmethod
    def from_root_row(cls, row: Mapping[str, Any]) -> "IndependentNativeSourceDefinition":
        """Parse the public root-table fields without a production dependency.

        ``calibration_json`` may either be a decoded mapping or a JSON string;
        importing ``json`` is kept local so the scientific formulas remain
        visible in the module-level dependency list.
        """

        calibration: Any = row.get("calibration_json", {})
        if isinstance(calibration, str):
            import json

            calibration = json.loads(calibration)
        if not isinstance(calibration, Mapping):
            raise ValueError("calibration_json must decode to a mapping")

        other_keys = (
            "cell_other_impermeant_osmoles_fmol",
            "other_impermeant_osmoles_fmol",
            "other_osmoles_fmol",
        )
        try:
            other = next(float(calibration[key]) for key in other_keys if key in calibration)
        except StopIteration as exc:
            raise ValueError("native root row omits fitted OTHER osmoles") from exc

        return cls(
            source_class=str(row["source_class"]),
            source_scale=float(row["source_scale"]),
            pump_capacity_scale=float(row["pump_capacity_scale"]),
            apical_pump_fraction=float(row["apical_pump_fraction"]),
            apical_k_fraction=float(row["apical_k_fraction"]),
            ae4_cation_fraction=float(row["ae4_cation_fraction"]),
            hydraulic_scale=float(row["hydraulic_scale"]),
            cell_other_impermeant_osmoles_fmol=other,
            hydraulic_mode=str(
                row.get("hydraulic_mode")
                or ("H1" if float(row["hydraulic_scale"]) == 1.0 else "HW")
            ),
        )

    def reconstruct(
        self,
        base: FullModelParameters | None = None,
    ) -> tuple[FullModelParameters, IndependentAE4Kinetics]:
        """Independently construct parameters and AE4 kinetics for the member."""

        source = FullModelParameters() if base is None else base
        n_abs = self.source_scale if self.source_class == N_ABS_NKCC else 1.0
        an_abs = self.source_scale if self.source_class == AN_ABS_AE4_NHE else 1.0

        parameters = replace(
            source,
            homeostasis=replace(
                source.homeostasis,
                nkcc1_capacity_fmol_s=FROZEN_NKCC1_CAPACITY_FMOL_S * n_abs,
                nhe1_capacity_fmol_s=FROZEN_NHE1_CAPACITY_FMOL_S * an_abs,
            ),
            membranes=replace(
                source.membranes,
                nak_capacity_fmol_s=(
                    source.membranes.nak_capacity_fmol_s
                    * self.pump_capacity_scale
                ),
                apical_pump_fraction=self.apical_pump_fraction,
                apical_k_fraction=self.apical_k_fraction,
                g_k_total_S=NATIVE_G_K_TOTAL_S,
                g_cl_apical_S=NATIVE_G_CL_APICAL_S,
                calcium_half_uM=NATIVE_CA_HALF_UM,
                calcium_hill=NATIVE_CA_HILL,
                g_apical_background_S=0.0,
                g_basolateral_background_S=0.0,
            ),
            water=independent_hydraulic_scaled_water_parameters(
                source.water,
                scale=self.hydraulic_scale,
                mode=self.hydraulic_mode,
            ),
            geometry=replace(
                source.geometry,
                cell_impermeant_osmoles_fmol=(
                    self.cell_other_impermeant_osmoles_fmol
                ),
            ),
        )

        fraction = self.ae4_cation_fraction
        ae4 = IndependentAE4Kinetics(
            carrier_amount_fmol=FROZEN_AE4_CARRIER_FMOL * an_abs,
            common_cl_attempt_rate_s=1.0,
            na_loaded_attempt_rate_s=2.0 * fraction,
            k_loaded_attempt_rate_s=2.0 * (1.0 - fraction),
        )
        return parameters, ae4

    def build_independent_model(
        self,
        *,
        regulation: Any | None = None,
        nkcc_regulation: Any | None = None,
        stimulus: Any | None = None,
    ) -> IndependentWholeCell:
        parameters, ae4 = self.reconstruct()
        return IndependentWholeCell(
            parameters,
            ae4_parameters=ae4,
            regulation=regulation,
            nkcc_regulation=nkcc_regulation,
            stimulus=stimulus,
        )


def state_from_independent_rest_coordinates(
    parameters: FullModelParameters,
    coordinates: ArrayLike,
) -> NDArray[np.float64]:
    """Map the ten physical root coordinates to exact charge-manifold amounts.

    This is a second transcription of the calibration coordinate map.  The
    coordinates are intracellular Na, K, TIC, pH, volume followed by luminal
    Na, K, TIC, pH, volume.  Concentration in mM times volume in pL is fmol.
    """

    x = np.asarray(coordinates, dtype=float)
    if x.shape != (10,) or np.any(~np.isfinite(x)):
        raise ValueError("independent rest coordinates must be finite length ten")
    na_i, k_i, tic_i, ph_i, volume_i, na_l, k_l, tic_l, ph_l, volume_l = (
        float(value) for value in x
    )
    if min(na_i, k_i, tic_i, volume_i, na_l, k_l, tic_l, volume_l) <= 0.0:
        raise ValueError("independent rest concentration/volume coordinates must be positive")
    ta_i = independent_total_alkalinity_mM(
        ph_i,
        tic_i,
        parameters.geometry.cell_buffer_total_fmol / volume_i,
        parameters.acid_base.cell_buffer_pka,
        parameters,
    )
    cl_i = (
        na_i
        + k_i
        - ta_i
        - parameters.geometry.fixed_cell_anion_equivalents_fmol / volume_i
    )
    ta_l = independent_total_alkalinity_mM(
        ph_l,
        tic_l,
        parameters.geometry.lumen_buffer_total_fmol / volume_l,
        parameters.acid_base.lumen_buffer_pka,
        parameters,
    )
    cl_l = na_l + k_l - ta_l
    if min(cl_i, cl_l) <= 0.0:
        raise ValueError("independent charge-manifold chloride must be positive")
    return np.asarray(
        (
            na_i * volume_i,
            k_i * volume_i,
            cl_i * volume_i,
            tic_i * volume_i,
            ta_i * volume_i,
            volume_i,
            na_l * volume_l,
            k_l * volume_l,
            cl_l * volume_l,
            tic_l * volume_l,
            ta_l * volume_l,
            volume_l,
        ),
        dtype=float,
    )


def refine_independent_native_root(
    definition: IndependentNativeSourceDefinition,
    candidate_state: ArrayLike,
    *,
    row_scales: Sequence[float] = (
        0.05,
        0.05,
        0.05,
        0.02,
        1.0e-4,
        0.02,
        0.02,
        0.02,
        0.02,
        1.0e-4,
    ),
    relative_bound: float = 0.01,
    max_nfev: int = 500,
) -> IndependentNativeRootResult:
    """Jointly refine ten charge coordinates and the sole OTHER nuisance.

    The local one-percent box is intentionally a reproduction check, not a
    discovery search.  It asks whether the supplied root is a root of the
    independently transcribed equations when OTHER is treated as the same
    eleventh unknown used by the production search.
    """

    candidate = np.asarray(candidate_state, dtype=float)
    if candidate.shape != (12,) or np.any(~np.isfinite(candidate)) or np.any(candidate <= 0.0):
        raise ValueError("native root candidate must be positive finite length 12")
    if not math.isfinite(relative_bound) or not 0.0 < relative_bound < 0.5:
        raise ValueError("native reproduction bound must lie in (0, 0.5)")
    scales = np.asarray(tuple(row_scales), dtype=float)
    if scales.shape != (10,) or np.any(~np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError("native root row scales must be ten positive values")

    initial = np.r_[
        free_charge_coordinates(candidate),
        definition.cell_other_impermeant_osmoles_fmol,
    ]
    lower = np.maximum(initial * (1.0 - relative_bound), 1.0e-12)
    upper = initial * (1.0 + relative_bound)
    fixed_charge = definition.reconstruct()[0].geometry.fixed_cell_anion_equivalents_fmol

    def evaluate(vector: NDArray[np.float64]):
        local_definition = replace(
            definition,
            cell_other_impermeant_osmoles_fmol=float(vector[-1]),
        )
        model = local_definition.build_independent_model(
            stimulus=IndependentConstantStimulus(0.058, 0.0)
        )
        state = state_from_free_charge_coordinates(
            vector[:10],
            fixed_cell_anion_equivalents_fmol=fixed_charge,
            regulation=(),
        )
        return model, state, model.evaluate(0.0, state)

    def residual(vector: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            _model, state, result = evaluate(vector)
            independent = result.rhs[list(NATIVE_ROOT_ROWS)] / scales
            return np.r_[independent, (state[5] - 1.3) / 0.13]
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(11, 1.0e6)

    fit = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        max_nfev=max_nfev,
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        x_scale="jac",
    )
    _model, state, result = evaluate(np.asarray(fit.x, dtype=float))
    scaled = residual(np.asarray(fit.x, dtype=float))
    coordinate_scales = np.maximum(np.abs(initial), 1.0e-12)
    normalized_jacobian = np.asarray(fit.jac, dtype=float) @ np.diag(coordinate_scales)
    singular = np.linalg.svd(normalized_jacobian, compute_uv=False)
    threshold = 1.0e-8 * singular[0] if len(singular) else math.inf
    rank = int(np.sum(singular > threshold))
    raw = np.asarray(result.rhs, dtype=float)
    amount_rows = (0, 1, 2, 3, 4, 6, 7, 8, 9, 10)
    return IndependentNativeRootResult(
        state=np.asarray(state, dtype=float),
        other_impermeant_osmoles_fmol=float(fit.x[-1]),
        success=bool(
            fit.success
            and np.max(np.abs(scaled)) <= 1.0e-7
            and rank == 11
            and np.all(state > 0.0)
        ),
        max_abs_scaled_residual=float(np.max(np.abs(scaled))),
        max_abs_amount_source_fmol_s=float(np.max(np.abs(raw[list(amount_rows)]))),
        max_abs_volume_source_pL_s=float(np.max(np.abs(raw[[5, 11]]))),
        nfev=int(fit.nfev),
        cost=float(fit.cost),
        scaled_jacobian_singular_values=tuple(float(value) for value in singular),
        jacobian_rank=rank,
        jacobian_nullity=11 - rank,
        message=str(fit.message),
    )


__all__ = [
    "AN_ABS_AE4_NHE",
    "FROZEN_AE4_CARRIER_FMOL",
    "FROZEN_NHE1_CAPACITY_FMOL_S",
    "FROZEN_NKCC1_CAPACITY_FMOL_S",
    "IndependentNativeRootResult",
    "IndependentNativeSourceDefinition",
    "INDEPENDENT_WT_BOUNDS",
    "INDEPENDENT_WT_POSITIVE_OBSERVABLES",
    "NATIVE_CA_HALF_UM",
    "NATIVE_CA_HILL",
    "NATIVE_G_CL_APICAL_S",
    "NATIVE_G_K_TOTAL_S",
    "NATIVE_ROOT_ROWS",
    "N_ABS_NKCC",
    "SOURCE_CLASSES",
    "STRICT_NATIVE",
    "independent_hydraulic_scaled_water_parameters",
    "independent_native_wt_gate_failures",
    "independent_native_wt_observables",
    "refine_independent_native_root",
    "state_from_independent_rest_coordinates",
]
