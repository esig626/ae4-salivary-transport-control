"""Native whole-cell source map and WT-only absolute-source root panel.

This module supersedes the earlier common-conductance-scale calibration for
the Task 13B pre-reveal reconstruction.  It is intentionally sealed from the
AE4-null secretion ledger.  The only fitted continuous nuisance in a native
root is the finite OTHER intracellular osmole amount; intracellular and
luminal chloride remain eliminated on exact bulk-charge manifolds.

The strict model uses the whole-cell source maxima and the Palk-lineage Ca
gate literally: 31.4 nS apical CaCC, 14 nS total Ca-activated K conductance,
``K_Ca=0.26 uM``, and ``eta=1.46``.  Total K conductance is partitioned
between membranes without a background conductance or common multiplier.

Two explicitly hypothetical absolute-capacity classes are nested around the
strict model. ``N_ABS_NKCC`` multiplies the complete reversible NKCC1
capacity. ``AN_ABS_AE4_NHE`` multiplies the AE4 carrier gauge and NHE1
capacity together.  Both preserve source stoichiometry, reversal, and the
frozen within-class ratios.  A scale of one is exactly identical to
``STRICT_NATIVE`` and must never be interpreted as an independently fitted
model.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .calibration import (
    CORE_OMITTED_CHARGE_ROWS,
    REST_COORDINATE_NAMES,
    WTCalibrationSpec,
    _independent_scaled_rhs,
    _observable_map,
    _raw_rhs_unit_maxima,
    _wt_gate_failures,
    encode_electroneutral_coordinates,
)
from .model import ConstantStimulus, ModernFullModel, WT
from .parameters import FullModelParameters
from .transporters import AE4Parameters


STRICT_NATIVE = "STRICT_NATIVE"
N_ABS_NKCC = "N_ABS_NKCC"
AN_ABS_AE4_NHE = "AN_ABS_AE4_NHE"
SOURCE_CLASSES = (STRICT_NATIVE, N_ABS_NKCC, AN_ABS_AE4_NHE)

G_CACC_MAX_S = 31.4e-9
G_K_TOTAL_MAX_S = 14.0e-9
PALK_CA_HALF_UM = 0.26
PALK_CA_HILL = 1.46

APICAL_PUMP_FRACTION_LOW = 0.056878
APICAL_PUMP_FRACTION_NOMINAL = 0.075075
APICAL_PUMP_FRACTION_HIGH = 0.094586
APICAL_K_FRACTIONS = (0.20, 0.30, 0.40)
AE4_NA_ROUTING_FRACTIONS = (0.05, 0.10, 0.20)

PRODUCTION_SOURCE_SCALES = (0.5, 1.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0)
DIAGNOSTIC_STRESS_SOURCE_SCALES = (12.0, 16.0, 24.0, 32.0)
HYDRAULIC_DIAGNOSTIC_SCALES = (1.0, 2.5, 3.5, 5.0, 6.75, 8.0)
CONDITIONAL_PUMP_CAPACITY_SCALES = (0.5, 1.0, 2.0)

# Frozen from the previous WT-only construction, before the native source
# reroot.  These are model construction values rather than measured molecular
# capacities.  The native protocol forbids refitting them independently.
FROZEN_AE4_CARRIER_FMOL = 0.1501264265070657
FROZEN_NHE1_CAPACITY_FMOL_S = 0.007440980590414482
FROZEN_NKCC1_CAPACITY_FMOL_S = 0.08
FROZEN_OTHER_OSMOLES_FMOL = 124.71681310732373
FROZEN_REFERENCE_COORDINATES = (
    14.756773495778985,
    111.10374637122905,
    6.194122572604099,
    6.909999999999942,
    1.3,
    146.13568410951905,
    4.087361120509596,
    12.891629639536415,
    7.276453648930862,
    0.10164106693766255,
)


@dataclass(frozen=True)
class NativeTopology:
    topology_id: str
    apical_pump_fraction: float
    apical_k_fraction: float


@dataclass(frozen=True)
class NativeSourcePanelMember:
    panel_id: str
    topology_id: str
    routing_id: str
    source_class: str
    source_scale: float
    pump_capacity_scale: float
    apical_pump_fraction: float
    apical_k_fraction: float
    ae4_cation_fraction: float
    hydraulic_scale: float = 1.0
    hydraulic_mode: str = "H1"
    eligibility: str = "PRODUCTION_PRE_REVEAL"
    parent_panel_id: str = "ROOT_GENERATION_NATIVE_SOURCE_MAP"

    def __post_init__(self) -> None:
        if self.source_class not in SOURCE_CLASSES:
            raise ValueError(f"unknown absolute-source class {self.source_class!r}")
        for name in ("source_scale", "pump_capacity_scale", "hydraulic_scale"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
        for name in (
            "apical_pump_fraction",
            "apical_k_fraction",
            "ae4_cation_fraction",
        ):
            value = float(getattr(self, name))
            if not 0.0 < value < 1.0:
                raise ValueError(f"{name} must lie strictly between zero and one")
        if self.source_class == STRICT_NATIVE and self.source_scale != 1.0:
            raise ValueError("STRICT_NATIVE requires source_scale=1 exactly")
        if self.hydraulic_mode not in ("H1", "HW", "HWQ"):
            raise ValueError("hydraulic_mode must be H1, HW, or HWQ")
        if self.hydraulic_scale == 1.0 and self.hydraulic_mode != "H1":
            raise ValueError("H=1 must use the public H1 identity mode")
        if self.hydraulic_scale != 1.0 and self.hydraulic_mode == "H1":
            raise ValueError("H>1 must declare HW or HWQ")
        if self.hydraulic_scale != 1.0 and "DIAGNOSTIC" not in self.eligibility:
            raise ValueError("H != 1 is diagnostic-only")


@dataclass(frozen=True)
class NativeRootAttempt:
    panel_id: str
    start_id: str
    optimizer_success: bool
    admissible_state: bool
    converged_root: bool
    passes_numerical_gate: bool
    passes_wt_gate: bool
    max_abs_scaled_residual: float
    max_abs_amount_rhs_fmol_s: float
    max_abs_volume_rhs_pL_s: float
    max_abs_omitted_charge_rhs_fmol_equivalent_s: float
    normalized_jacobian_rank: int
    normalized_jacobian_nullity: int
    boundary_hits: tuple[str, ...]
    wt_gate_failures: tuple[str, ...]
    coordinates_and_other: tuple[float, ...]
    cost: float
    optimality: float
    nfev: int
    message: str


@dataclass(frozen=True)
class NativeRootRecord:
    root_id: str
    panel_id: str
    topology_id: str
    routing_id: str
    source_class: str
    source_scale: float
    pump_capacity_scale: float
    apical_pump_fraction: float
    apical_k_fraction: float
    ae4_cation_fraction: float
    hydraulic_scale: float
    hydraulic_mode: str
    eligibility: str
    passes_numerical_gate: bool
    passes_wt_gate: bool
    member_start_ids: tuple[str, ...]
    core_state: tuple[float, ...]
    coordinates: tuple[float, ...]
    other_impermeant_osmoles_fmol: float
    observables: Mapping[str, float]
    max_abs_scaled_residual: float
    max_abs_amount_rhs_fmol_s: float
    max_abs_volume_rhs_pL_s: float
    max_abs_omitted_charge_rhs_fmol_equivalent_s: float
    max_abs_charge_residual_fmol: float
    max_abs_current_residual_A: float
    normalized_jacobian_singular_values: tuple[float, ...]
    normalized_jacobian_rank: int
    normalized_jacobian_nullity: int
    boundary_hits: tuple[str, ...]
    gate_failures: tuple[str, ...]


@dataclass(frozen=True)
class NativeRootSearchReport:
    member: NativeSourcePanelMember
    attempts: tuple[NativeRootAttempt, ...]
    roots: tuple[NativeRootRecord, ...]
    start_design: Mapping[str, Any]


def native_topologies() -> tuple[NativeTopology, ...]:
    """Return the nominal point and four corners of the declared rectangle."""

    return (
        NativeTopology("PLOW_KLOW", APICAL_PUMP_FRACTION_LOW, 0.20),
        NativeTopology("PLOW_KHIGH", APICAL_PUMP_FRACTION_LOW, 0.40),
        NativeTopology("PNOM_KMID", APICAL_PUMP_FRACTION_NOMINAL, 0.30),
        NativeTopology("PHIGH_KLOW", APICAL_PUMP_FRACTION_HIGH, 0.20),
        NativeTopology("PHIGH_KHIGH", APICAL_PUMP_FRACTION_HIGH, 0.40),
    )


def _routing_id(fraction: float) -> str:
    return f"AE4NA{int(round(100.0 * fraction)):02d}"


def declared_native_panel(
    *,
    source_class: str = STRICT_NATIVE,
    source_scale: float = 1.0,
    pump_capacity_scale: float = 1.0,
    hydraulic_scale: float = 1.0,
    hydraulic_mode: str | None = None,
    eligibility: str = "PRODUCTION_PRE_REVEAL",
) -> tuple[NativeSourcePanelMember, ...]:
    """Return the full five-topology by three-routing deterministic slice."""

    members: list[NativeSourcePanelMember] = []
    mode = hydraulic_mode or ("H1" if hydraulic_scale == 1.0 else "HW")
    for topology in native_topologies():
        for routing in AE4_NA_ROUTING_FRACTIONS:
            routing_id = _routing_id(routing)
            hydraulic_id = (
                f"H{hydraulic_scale:g}"
                if mode != "HWQ"
                else f"HWQ_H{hydraulic_scale:g}"
            )
            panel_id = (
                f"{source_class}_S{source_scale:g}_{topology.topology_id}_"
                f"{routing_id}_P{pump_capacity_scale:g}_{hydraulic_id}"
            )
            strict_pump_one_id = (
                f"{STRICT_NATIVE}_S1_{topology.topology_id}_{routing_id}_P1_H1"
            )
            strict_same_pump_id = (
                f"{STRICT_NATIVE}_S1_{topology.topology_id}_{routing_id}_"
                f"P{pump_capacity_scale:g}_H1"
            )
            if source_class != STRICT_NATIVE:
                parent_panel_id = strict_same_pump_id
            elif hydraulic_scale != 1.0:
                parent_panel_id = strict_same_pump_id
            elif pump_capacity_scale != 1.0:
                parent_panel_id = strict_pump_one_id
            else:
                parent_panel_id = "ROOT_GENERATION_NATIVE_SOURCE_MAP"
            members.append(
                NativeSourcePanelMember(
                    panel_id=panel_id,
                    topology_id=topology.topology_id,
                    routing_id=routing_id,
                    source_class=source_class,
                    source_scale=source_scale,
                    pump_capacity_scale=pump_capacity_scale,
                    apical_pump_fraction=topology.apical_pump_fraction,
                    apical_k_fraction=topology.apical_k_fraction,
                    ae4_cation_fraction=routing,
                    hydraulic_scale=hydraulic_scale,
                    hydraulic_mode=mode,
                    eligibility=eligibility,
                    parent_panel_id=parent_panel_id,
                )
            )
    return tuple(members)


def declared_absolute_source_panel(
    *,
    include_identity: bool = True,
    pump_capacity_scales: Sequence[float] = (1.0,),
) -> tuple[NativeSourcePanelMember, ...]:
    """Return production N/AN rows; nonunit rows are context hypotheses."""

    members: list[NativeSourcePanelMember] = []
    for pump_scale in pump_capacity_scales:
        for source_class in (N_ABS_NKCC, AN_ABS_AE4_NHE):
            for scale in PRODUCTION_SOURCE_SCALES:
                if scale == 1.0 and not include_identity:
                    continue
                members.extend(
                    declared_native_panel(
                        source_class=source_class,
                        source_scale=scale,
                        pump_capacity_scale=float(pump_scale),
                    )
                )
    return tuple(members)


def declared_absolute_source_stress_panel(
    *, pump_capacity_scales: Sequence[float] = (1.0,)
) -> tuple[NativeSourcePanelMember, ...]:
    """Return non-production stress rows beyond the declared source grid."""

    members: list[NativeSourcePanelMember] = []
    for pump_scale in pump_capacity_scales:
        for source_class in (N_ABS_NKCC, AN_ABS_AE4_NHE):
            for scale in DIAGNOSTIC_STRESS_SOURCE_SCALES:
                members.extend(
                    declared_native_panel(
                        source_class=source_class,
                        source_scale=scale,
                        pump_capacity_scale=float(pump_scale),
                        eligibility="DIAGNOSTIC_ONLY_ABSOLUTE_SOURCE_STRESS",
                    )
                )
    return tuple(members)


def declared_hydraulic_diagnostic_panel(
    *,
    pump_capacity_scales: Sequence[float] = (1.0,),
    hydraulic_modes: Sequence[str] = ("HW", "HWQ"),
) -> tuple[NativeSourcePanelMember, ...]:
    """Return H>1 reroots, which can diagnose but never rescue production."""

    members: list[NativeSourcePanelMember] = []
    for mode in hydraulic_modes:
        if mode not in ("HW", "HWQ"):
            raise ValueError("diagnostic hydraulic modes must be HW or HWQ")
        for pump_scale in pump_capacity_scales:
            for scale in HYDRAULIC_DIAGNOSTIC_SCALES:
                if scale == 1.0:
                    continue
                members.extend(
                    declared_native_panel(
                        pump_capacity_scale=float(pump_scale),
                        hydraulic_scale=scale,
                        hydraulic_mode=mode,
                        eligibility=f"DIAGNOSTIC_ONLY_{mode}_HYDRAULIC_TRANSFER",
                    )
                )
    return tuple(members)


def native_source_map() -> Mapping[str, Any]:
    """Machine-readable source/assumption registry for the native panel."""

    return {
        "g_cacc_max_S": G_CACC_MAX_S,
        "g_cacc_topology": "whole-cell maximum assigned wholly apical",
        "g_k_total_max_S": G_K_TOTAL_MAX_S,
        "g_k_topology": "whole-cell maximum partitioned only; no background or common scale",
        "calcium_gate": {
            "half_uM": PALK_CA_HALF_UM,
            "hill": PALK_CA_HILL,
            "formula": "Ca^eta/(Ca^eta+K_Ca^eta)",
        },
        "apical_pump_fraction": {
            "low": APICAL_PUMP_FRACTION_LOW,
            "nominal": APICAL_PUMP_FRACTION_NOMINAL,
            "high": APICAL_PUMP_FRACTION_HIGH,
            "status": "derived parotid-to-mouse-SMG transfer sensitivity",
        },
        "apical_k_fraction": {
            "values": list(APICAL_K_FRACTIONS),
            "status": "model/source bracket; not a measured interval",
        },
        "ae4_na_routing_fraction": {
            "values": list(AE4_NA_ROUTING_FRACTIONS),
            "status": "mixed-bath model assumptions; pure-cation assays do not identify it",
        },
        "source_classes": {
            STRICT_NATIVE: "literal source map; all absolute multipliers one",
            N_ABS_NKCC: "one scalar on complete reversible NKCC1 capacity",
            AN_ABS_AE4_NHE: "one common scalar on AE4 carrier gauge and NHE1 capacity",
        },
        "production_source_scales": list(PRODUCTION_SOURCE_SCALES),
        "diagnostic_stress_source_scales": list(DIAGNOSTIC_STRESS_SOURCE_SCALES),
        "hydraulic_diagnostic_scales": list(HYDRAULIC_DIAGNOSTIC_SCALES),
        "hydraulic_diagnostic_modes": {
            "HW": "common scalar on Pa, Pb, and Pt only",
            "HWQ": "same scalar on Pa, Pb, Pt, and local lumen outflow_rate_s",
            "eligibility": "diagnostic only for every H>1 row",
        },
        "frozen_capacity_construction": {
            "ae4_carrier_fmol": FROZEN_AE4_CARRIER_FMOL,
            "nhe1_capacity_fmol_s": FROZEN_NHE1_CAPACITY_FMOL_S,
            "nkcc1_capacity_fmol_s": FROZEN_NKCC1_CAPACITY_FMOL_S,
        },
        "root_fit": {
            "coordinates": [*REST_COORDINATE_NAMES, "cell_other_impermeant_osmoles_fmol"],
            "equations": "ten independent charge-manifold steady RHS plus Vi=1.3 pL",
            "continuous_nuisance": "OTHER intracellular impermeant osmoles only",
        },
        "firewall": "WT-only; no AE4-null secretion datum is imported or evaluated",
    }


def hydraulic_scaled_water_parameters(
    water_parameters: Any, *, scale: float, mode: str
) -> Any:
    """Pure HW/HWQ transform; all named modes are exact at scale one."""

    value = float(scale)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("hydraulic diagnostic scale must be finite and positive")
    if mode not in ("H1", "HW", "HWQ"):
        raise ValueError("hydraulic diagnostic mode must be H1, HW, or HWQ")
    if mode == "H1" and value != 1.0:
        raise ValueError("the H1 transform requires scale one")
    return replace(
        water_parameters,
        apical_hydraulic_pL_s_mOsm=(
            water_parameters.apical_hydraulic_pL_s_mOsm * value
        ),
        basolateral_hydraulic_pL_s_mOsm=(
            water_parameters.basolateral_hydraulic_pL_s_mOsm * value
        ),
        paracellular_hydraulic_pL_s_mOsm=(
            water_parameters.paracellular_hydraulic_pL_s_mOsm * value
        ),
        outflow_rate_s=(
            water_parameters.outflow_rate_s * value
            if mode == "HWQ"
            else water_parameters.outflow_rate_s
        ),
    )


def build_native_source_model(
    member: NativeSourcePanelMember | NativeRootRecord,
    *,
    other_impermeant_osmoles_fmol: float | None = None,
    regulatory_model: Any | None = None,
    stimulus: Any | None = None,
    base_parameters: FullModelParameters | None = None,
) -> ModernFullModel:
    """Build one literal-source full model for roots or downstream dynamics."""

    if isinstance(member, NativeRootRecord):
        root = member
        member = NativeSourcePanelMember(
            panel_id=root.panel_id,
            topology_id=root.topology_id,
            routing_id=root.routing_id,
            source_class=root.source_class,
            source_scale=root.source_scale,
            pump_capacity_scale=root.pump_capacity_scale,
            apical_pump_fraction=root.apical_pump_fraction,
            apical_k_fraction=root.apical_k_fraction,
            ae4_cation_fraction=root.ae4_cation_fraction,
            hydraulic_scale=root.hydraulic_scale,
            hydraulic_mode=root.hydraulic_mode,
            eligibility=root.eligibility,
            parent_panel_id="LOADED_ROOT_RECORD",
        )
        if other_impermeant_osmoles_fmol is None:
            other_impermeant_osmoles_fmol = root.other_impermeant_osmoles_fmol
    if other_impermeant_osmoles_fmol is None:
        raise ValueError("OTHER impermeant osmoles are required for a panel member")
    other = float(other_impermeant_osmoles_fmol)
    if not math.isfinite(other) or other <= 0.0:
        raise ValueError("OTHER impermeant osmoles must be finite and positive")
    base = base_parameters or FullModelParameters()
    n_scale = member.source_scale if member.source_class == N_ABS_NKCC else 1.0
    an_scale = member.source_scale if member.source_class == AN_ABS_AE4_NHE else 1.0
    parameters = replace(
        base,
        homeostasis=replace(
            base.homeostasis,
            nkcc1_capacity_fmol_s=FROZEN_NKCC1_CAPACITY_FMOL_S * n_scale,
            nhe1_capacity_fmol_s=FROZEN_NHE1_CAPACITY_FMOL_S * an_scale,
        ),
        membranes=replace(
            base.membranes,
            nak_capacity_fmol_s=(
                base.membranes.nak_capacity_fmol_s * member.pump_capacity_scale
            ),
            apical_pump_fraction=member.apical_pump_fraction,
            apical_k_fraction=member.apical_k_fraction,
            g_k_total_S=G_K_TOTAL_MAX_S,
            g_cl_apical_S=G_CACC_MAX_S,
            calcium_half_uM=PALK_CA_HALF_UM,
            calcium_hill=PALK_CA_HILL,
            g_apical_background_S=0.0,
            g_basolateral_background_S=0.0,
        ),
        water=hydraulic_scaled_water_parameters(
            base.water,
            scale=member.hydraulic_scale,
            mode=member.hydraulic_mode,
        ),
        geometry=replace(
            base.geometry,
            cell_impermeant_osmoles_fmol=other,
        ),
    )
    na_fraction = member.ae4_cation_fraction
    ae4 = AE4Parameters(
        carrier_amount_fmol=FROZEN_AE4_CARRIER_FMOL * an_scale,
        common_cl_attempt_rate_s=1.0,
        na_loaded_attempt_rate_s=2.0 * na_fraction,
        k_loaded_attempt_rate_s=2.0 * (1.0 - na_fraction),
    )
    return ModernFullModel(
        parameters,
        stimulus=stimulus or ConstantStimulus(calcium_uM=0.058, beta_input=0.0),
        regulatory_model=regulatory_model,
        ae4_parameters=ae4,
    )


def native_root_residual(
    member: NativeSourcePanelMember,
    coordinates_and_other: ArrayLike,
    *,
    spec: WTCalibrationSpec | None = None,
) -> NDArray[np.float64]:
    """Eleven-row root residual with OTHER as the only fitted nuisance."""

    declared = spec or WTCalibrationSpec()
    vector = np.asarray(coordinates_and_other, dtype=float)
    if vector.shape != (11,):
        raise ValueError("native root vector must have ten coordinates plus OTHER")
    model = build_native_source_model(
        member, other_impermeant_osmoles_fmol=float(vector[-1])
    )
    return np.r_[
        _independent_scaled_rhs(model, vector[:10], declared),
        (float(vector[4]) - declared.volume_lineage_target_pL) / 0.13,
    ]


def _root_bounds(
    spec: WTCalibrationSpec,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    lower = np.r_[spec.coordinate_lower, 20.0]
    preferred = np.r_[np.asarray(FROZEN_REFERENCE_COORDINATES), FROZEN_OTHER_OSMOLES_FMOL]
    upper = np.r_[spec.coordinate_upper, 300.0]
    return lower, preferred, upper


def _named_starts(
    lower: NDArray[np.float64],
    preferred: NDArray[np.float64],
    upper: NDArray[np.float64],
    *,
    count: int,
    seed: int,
    anchors: Sequence[ArrayLike] = (),
) -> tuple[tuple[str, NDArray[np.float64]], ...]:
    if count not in (3, 7, 33):
        raise ValueError("native multistart count must be exactly 3, 7, or 33")
    unique_anchors: list[NDArray[np.float64]] = []
    for anchor in anchors:
        candidate = np.asarray(anchor, dtype=float)
        if candidate.shape != preferred.shape:
            raise ValueError("native root anchor has the wrong shape")
        if np.any(candidate <= lower) or np.any(candidate >= upper):
            raise ValueError("native root anchors must be strictly interior")
        if not any(np.array_equal(candidate, prior) for prior in unique_anchors):
            unique_anchors.append(candidate.copy())
    lower_quartile = lower + 0.25 * (upper - lower)
    upper_quartile = lower + 0.75 * (upper - lower)
    if count == 3:
        primary_name = "parent_continuation_00" if unique_anchors else "frozen_reference"
        primary = unique_anchors[0] if unique_anchors else preferred.copy()
        if np.array_equal(primary, lower_quartile) or np.array_equal(
            primary, upper_quartile
        ):
            raise ValueError(
                "the three-start anchor must be distinct from both quartile starts"
            )
        starts: list[tuple[str, NDArray[np.float64]]] = [
            (primary_name, primary.copy()),
            ("lower_quartile", lower_quartile),
            ("upper_quartile", upper_quartile),
        ]
    else:
        starts = [
            ("frozen_reference", preferred.copy()),
            ("lower_quartile", lower_quartile),
            ("upper_quartile", upper_quartile),
        ]
        for index, candidate in enumerate(unique_anchors):
            if any(np.array_equal(candidate, item[1]) for item in starts):
                continue
            starts.append((f"parent_continuation_{index:02d}", candidate.copy()))
    if len(starts) < count:
        rng = np.random.default_rng(seed)
        while len(starts) < count:
            point = lower + rng.random(len(lower)) * (upper - lower)
            starts.append((f"uniform_{len(starts):03d}", point))
    return tuple(starts[:count])


def _root_components(
    normalized_points: Sequence[ArrayLike], *, tolerance: float
) -> tuple[tuple[int, ...], ...]:
    """Return deterministic connected components under the root tolerance.

    Connected components make clustering independent of discovery order.  They
    also guarantee that representatives chosen from distinct components are
    separated by more than ``tolerance`` in the max norm.
    """

    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("root clustering tolerance must be finite and positive")
    points = tuple(np.asarray(point, dtype=float) for point in normalized_points)
    if not points:
        return ()
    shape = points[0].shape
    if any(point.shape != shape or np.any(~np.isfinite(point)) for point in points):
        raise ValueError("root clustering points must share one finite shape")

    parents = list(range(len(points)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)

    for left in range(len(points)):
        for right in range(left + 1, len(points)):
            if float(np.max(np.abs(points[left] - points[right]))) <= tolerance:
                union(left, right)

    grouped: dict[int, list[int]] = {}
    for index in range(len(points)):
        grouped.setdefault(find(index), []).append(index)
    components = [tuple(sorted(indices)) for indices in grouped.values()]
    return tuple(sorted(components, key=lambda item: tuple(points[item[0]])))


def _boundary_hits(
    vector: NDArray[np.float64],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    tolerance: float,
) -> tuple[str, ...]:
    names = (*REST_COORDINATE_NAMES, "cell_other_impermeant_osmoles_fmol")
    distance = np.minimum(vector - lower, upper - vector) / (upper - lower)
    return tuple(name for name, value in zip(names, distance) if value <= tolerance)


def solve_native_member(
    member: NativeSourcePanelMember,
    *,
    spec: WTCalibrationSpec | None = None,
    start_count: int = 3,
    seed: int = 13601,
    max_nfev: int = 3000,
    anchors: Sequence[ArrayLike] = (),
) -> NativeRootSearchReport:
    """Solve and retain every detected root for one fixed native slice."""

    declared = spec or WTCalibrationSpec()
    lower, preferred, upper = _root_bounds(declared)
    starts = _named_starts(
        lower, preferred, upper, count=start_count, seed=seed, anchors=anchors
    )
    span = upper - lower

    def residual(vector: NDArray[np.float64]) -> NDArray[np.float64]:
        try:
            return native_root_residual(member, vector, spec=declared)
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            return np.full(11, 1.0e6)

    attempts: list[NativeRootAttempt] = []
    successful: list[tuple[str, Any, Mapping[str, Any]]] = []
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
        details: dict[str, Any] = {}
        try:
            model = build_native_source_model(
                member, other_impermeant_osmoles_fmol=float(fit.x[-1])
            )
            state = encode_electroneutral_coordinates(model, fit.x[:10])
            evaluation = model.evaluate(0.0, state, genotype=WT)
            scaled = residual(fit.x)
            scaled_max = float(np.max(np.abs(scaled)))
            amount_max, volume_max, _, omitted_max = _raw_rhs_unit_maxima(
                evaluation.rhs
            )
            observables = _observable_map(model, state)
            wt_failures = _wt_gate_failures(observables, declared)
            charge_max = max(
                abs(float(value))
                for value in evaluation.diagnostics.state_charge_fmol.values()
            )
            current_max = max(
                abs(float(value))
                for value in evaluation.diagnostics.membranes.current_residuals_A.values()
            )
            singular = np.linalg.svd(
                np.asarray(fit.jac, dtype=float) * span[None, :],
                compute_uv=False,
            )
            threshold = (
                declared.jacobian_rank_relative_tolerance * singular[0]
                if len(singular)
                else math.inf
            )
            rank = int(np.sum(singular > threshold))
            hits = _boundary_hits(
                np.asarray(fit.x), lower, upper, declared.boundary_relative_tolerance
            )
            numerical = bool(
                fit.success
                and scaled_max <= declared.root_scaled_tolerance
                and omitted_max <= declared.omitted_row_raw_tolerance
                and charge_max <= declared.charge_tolerance_fmol
                and current_max <= declared.current_tolerance_A
                and rank == 11
                and not hits
            )
            passes_wt = numerical and not wt_failures
            details = {
                "model": model,
                "state": state,
                "observables": observables,
                "scaled_max": scaled_max,
                "amount_max": amount_max,
                "volume_max": volume_max,
                "omitted_max": omitted_max,
                "charge_max": charge_max,
                "current_max": current_max,
                "singular": singular,
                "rank": rank,
                "hits": hits,
                "wt_failures": wt_failures,
                "numerical": numerical,
                "passes_wt": passes_wt,
            }
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError):
            admissible = False
            scaled_max = amount_max = volume_max = omitted_max = math.inf
            rank = 0
            hits = ()
            wt_failures = ("inadmissible physical state",)
            numerical = passes_wt = False
        attempt = NativeRootAttempt(
            panel_id=member.panel_id,
            start_id=start_id,
            optimizer_success=bool(fit.success),
            admissible_state=admissible,
            converged_root=bool(admissible and scaled_max <= declared.root_scaled_tolerance),
            passes_numerical_gate=numerical,
            passes_wt_gate=passes_wt,
            max_abs_scaled_residual=scaled_max,
            max_abs_amount_rhs_fmol_s=amount_max,
            max_abs_volume_rhs_pL_s=volume_max,
            max_abs_omitted_charge_rhs_fmol_equivalent_s=omitted_max,
            normalized_jacobian_rank=rank,
            normalized_jacobian_nullity=11 - rank,
            boundary_hits=tuple(hits),
            wt_gate_failures=tuple(wt_failures),
            coordinates_and_other=tuple(float(value) for value in fit.x),
            cost=float(fit.cost),
            optimality=float(fit.optimality),
            nfev=int(fit.nfev),
            message=str(fit.message),
        )
        attempts.append(attempt)
        if numerical:
            successful.append((start_id, fit, details))

    components = _root_components(
        [item[1].x / span for item in successful],
        tolerance=declared.cluster_relative_tolerance,
    )
    clusters: list[
        tuple[
            tuple[str, Any, Mapping[str, Any]],
            tuple[tuple[str, Any, Mapping[str, Any]], ...],
        ]
    ] = []
    for component in components:
        members = tuple(successful[index] for index in component)
        representative = min(
            members,
            key=lambda item: (
                float(item[2]["scaled_max"]),
                tuple(float(value) for value in item[1].x),
                item[0],
            ),
        )
        clusters.append((representative, members))
    clusters.sort(
        key=lambda item: (
            tuple(float(value) for value in item[0][1].x / span),
            item[0][0],
        )
    )

    roots: list[NativeRootRecord] = []
    for index, (representative, cluster) in enumerate(clusters):
        start_id, fit, details = representative
        failures = list(details["wt_failures"])
        if not details["numerical"]:
            failures.append("numerical/root/rank/boundary gate failed")
        roots.append(
            NativeRootRecord(
                root_id=f"{member.panel_id}:B{index:02d}",
                panel_id=member.panel_id,
                topology_id=member.topology_id,
                routing_id=member.routing_id,
                source_class=member.source_class,
                source_scale=member.source_scale,
                pump_capacity_scale=member.pump_capacity_scale,
                apical_pump_fraction=member.apical_pump_fraction,
                apical_k_fraction=member.apical_k_fraction,
                ae4_cation_fraction=member.ae4_cation_fraction,
                hydraulic_scale=member.hydraulic_scale,
                hydraulic_mode=member.hydraulic_mode,
                eligibility=member.eligibility,
                passes_numerical_gate=bool(details["numerical"]),
                passes_wt_gate=bool(details["passes_wt"]),
                member_start_ids=tuple(sorted(item[0] for item in cluster)),
                core_state=tuple(float(value) for value in details["state"][:12]),
                coordinates=tuple(float(value) for value in fit.x[:10]),
                other_impermeant_osmoles_fmol=float(fit.x[-1]),
                observables={
                    key: float(value) for key, value in details["observables"].items()
                },
                max_abs_scaled_residual=float(details["scaled_max"]),
                max_abs_amount_rhs_fmol_s=float(details["amount_max"]),
                max_abs_volume_rhs_pL_s=float(details["volume_max"]),
                max_abs_omitted_charge_rhs_fmol_equivalent_s=float(
                    details["omitted_max"]
                ),
                max_abs_charge_residual_fmol=float(details["charge_max"]),
                max_abs_current_residual_A=float(details["current_max"]),
                normalized_jacobian_singular_values=tuple(
                    float(value) for value in details["singular"]
                ),
                normalized_jacobian_rank=int(details["rank"]),
                normalized_jacobian_nullity=11 - int(details["rank"]),
                boundary_hits=tuple(details["hits"]),
                gate_failures=tuple(failures),
            )
        )
    return NativeRootSearchReport(
        member=member,
        attempts=tuple(attempts),
        roots=tuple(roots),
        start_design={
            "kind": (
                "one_parent_or_frozen_anchor_plus_both_quartiles"
                if start_count == 3
                else "frozen_reference_both_quartiles_unique_parent_anchors_seeded_uniform_fill"
            ),
            "count": len(starts),
            "start_ids": [start_id for start_id, _ in starts],
            "unique_continuation_anchor_count": len(
                [start_id for start_id, _ in starts if start_id.startswith("parent_continuation_")]
            ),
            "seed": seed,
            "max_nfev": max_nfev,
            "continuous_nuisance": "cell_other_impermeant_osmoles_fmol",
        },
    )


def root_to_csv_row(root: NativeRootRecord) -> Mapping[str, Any]:
    """Stable interchange schema used by the dynamic and independent agents."""

    return {
        "root_id": root.root_id,
        "panel_id": root.panel_id,
        "topology_id": root.topology_id,
        "routing_id": root.routing_id,
        "source_class": root.source_class,
        "source_scale": root.source_scale,
        "pump_capacity_scale": root.pump_capacity_scale,
        "apical_pump_fraction": root.apical_pump_fraction,
        "apical_k_fraction": root.apical_k_fraction,
        "ae4_cation_fraction": root.ae4_cation_fraction,
        "hydraulic_scale": root.hydraulic_scale,
        "hydraulic_mode": root.hydraulic_mode,
        "eligibility": root.eligibility,
        "passes_numerical_gate": root.passes_numerical_gate,
        "passes_wt_gate": root.passes_wt_gate,
        "core_state_json": json.dumps(root.core_state, separators=(",", ":")),
        "coordinates_json": json.dumps(root.coordinates, separators=(",", ":")),
        "calibration_json": json.dumps(
            {"cell_other_impermeant_osmoles_fmol": root.other_impermeant_osmoles_fmol},
            sort_keys=True,
            separators=(",", ":"),
        ),
        "observables_json": json.dumps(
            root.observables, sort_keys=True, separators=(",", ":")
        ),
        "member_start_ids_json": json.dumps(root.member_start_ids),
        "max_abs_scaled_residual": root.max_abs_scaled_residual,
        "max_abs_amount_rhs_fmol_s": root.max_abs_amount_rhs_fmol_s,
        "max_abs_volume_rhs_pL_s": root.max_abs_volume_rhs_pL_s,
        "max_abs_omitted_charge_rhs_fmol_equivalent_s": (
            root.max_abs_omitted_charge_rhs_fmol_equivalent_s
        ),
        "max_abs_charge_residual_fmol": root.max_abs_charge_residual_fmol,
        "max_abs_current_residual_A": root.max_abs_current_residual_A,
        "normalized_jacobian_rank": root.normalized_jacobian_rank,
        "normalized_jacobian_nullity": root.normalized_jacobian_nullity,
        "boundary_hits_json": json.dumps(root.boundary_hits),
        "gate_failures_json": json.dumps(root.gate_failures),
    }


def load_native_source_roots(path: str | Path) -> tuple[NativeRootRecord, ...]:
    """Load the stable CSV interchange without touching any held-out file."""

    source = Path(path)
    if source.is_dir():
        source = source / "native_source_wt_roots.csv"
    records: list[NativeRootRecord] = []
    with source.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            observables = json.loads(row["observables_json"])
            calibration = json.loads(row["calibration_json"])
            records.append(
                NativeRootRecord(
                    root_id=row["root_id"],
                    panel_id=row["panel_id"],
                    topology_id=row["topology_id"],
                    routing_id=row["routing_id"],
                    source_class=row["source_class"],
                    source_scale=float(row["source_scale"]),
                    pump_capacity_scale=float(row["pump_capacity_scale"]),
                    apical_pump_fraction=float(row["apical_pump_fraction"]),
                    apical_k_fraction=float(row["apical_k_fraction"]),
                    ae4_cation_fraction=float(row["ae4_cation_fraction"]),
                    hydraulic_scale=float(row["hydraulic_scale"]),
                    hydraulic_mode=(
                        row.get("hydraulic_mode")
                        or ("H1" if float(row["hydraulic_scale"]) == 1.0 else "HW")
                    ),
                    eligibility=row["eligibility"],
                    passes_numerical_gate=row["passes_numerical_gate"] == "True",
                    passes_wt_gate=row["passes_wt_gate"] == "True",
                    member_start_ids=tuple(json.loads(row["member_start_ids_json"])),
                    core_state=tuple(float(value) for value in json.loads(row["core_state_json"])),
                    coordinates=tuple(
                        float(value) for value in json.loads(row["coordinates_json"])
                    ),
                    other_impermeant_osmoles_fmol=float(
                        calibration["cell_other_impermeant_osmoles_fmol"]
                    ),
                    observables={key: float(value) for key, value in observables.items()},
                    max_abs_scaled_residual=float(row["max_abs_scaled_residual"]),
                    max_abs_amount_rhs_fmol_s=float(row["max_abs_amount_rhs_fmol_s"]),
                    max_abs_volume_rhs_pL_s=float(row["max_abs_volume_rhs_pL_s"]),
                    max_abs_omitted_charge_rhs_fmol_equivalent_s=float(
                        row["max_abs_omitted_charge_rhs_fmol_equivalent_s"]
                    ),
                    max_abs_charge_residual_fmol=float(
                        row["max_abs_charge_residual_fmol"]
                    ),
                    max_abs_current_residual_A=float(row["max_abs_current_residual_A"]),
                    normalized_jacobian_singular_values=(),
                    normalized_jacobian_rank=int(row["normalized_jacobian_rank"]),
                    normalized_jacobian_nullity=int(row["normalized_jacobian_nullity"]),
                    boundary_hits=tuple(json.loads(row["boundary_hits_json"])),
                    gate_failures=tuple(json.loads(row["gate_failures_json"])),
                )
            )
    return tuple(records)


__all__ = (
    "AN_ABS_AE4_NHE",
    "N_ABS_NKCC",
    "STRICT_NATIVE",
    "NativeRootAttempt",
    "NativeRootRecord",
    "NativeRootSearchReport",
    "NativeSourcePanelMember",
    "build_native_source_model",
    "declared_absolute_source_panel",
    "declared_absolute_source_stress_panel",
    "declared_hydraulic_diagnostic_panel",
    "declared_native_panel",
    "hydraulic_scaled_water_parameters",
    "load_native_source_roots",
    "native_root_residual",
    "native_source_map",
    "native_topologies",
    "root_to_csv_row",
    "solve_native_member",
    "_root_components",
)
