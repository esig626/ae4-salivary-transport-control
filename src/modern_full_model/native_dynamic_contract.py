"""Pre-reveal WT dynamic-contract helpers for the native-source panel.

This module contains no genotype perturbation and never opens a validation or
held-out target ledger.  It turns single-cell WT flow trajectories into the
scale, geometry, sustainment, and numerical checks that must be satisfied
before a frozen Task 13B ensemble can become eligible for phenotype reveal.

The absolute gland-geometry ceiling is derived from Kondo et al. (2015): a
generous one-SMG mass ceiling of 100 mg, 66.91 percent true acinar volume, and
the frozen 1.3 pL model cell volume.  The mass ceiling is deliberately
generous.  It is not a fitted cell count and it cannot be relaxed by a WT flow
normalization.  Kondo et al. (2019) collected one ipsilateral SMG, so the
one-gland bound is the production convention; a paired-gland value is exposed
only as a labelled sensitivity.

The 9--10 uL/min WT flow envelope is graphical rather than a digitized error
model.  For each positive per-cell flow sample ``q``, an admissible whole-
gland multiplier ``S`` must satisfy ``9/q <= S <= 10/q``.  A structural row
receives one shared interval formed by intersecting these bounds across every
predeclared regulatory member, protocol arm, production solver, and minute.
Only existence of a geometry-compatible interval is used; no member-specific
or phenotype-dependent scale is fitted.

Provenance tags
---------------
``PRIMARY_MEASUREMENT``: SMG acinar fraction and WT graphical flow envelope.
``DERIVED_CONSTRAINT``: cell-count ceiling, unit conversion, interval algebra.
``NEW_MODELING_DECISION``: generous 100-mg ceiling and paired-gland sensitivity.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from collections import Counter
from typing import Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


WT_FLOW_LOWER_UL_MIN = 9.0
WT_FLOW_UPPER_UL_MIN = 10.0
ONE_SMG_MASS_CEILING_MG = 100.0
TRUE_ACINAR_VOLUME_FRACTION = 0.6691
FROZEN_CELL_VOLUME_PL = 1.3
PL_S_CELL_TO_UL_MIN = 6.0e-5

# 100 mg tissue is treated as 100 uL for this deliberately generous upper
# bound; 1 uL = 1e6 pL.  Rounding upward to the nearest cell is conservative.
ONE_SMG_EFFECTIVE_CELL_COUNT_MAX = 51_469_231
PAIRED_SMG_EFFECTIVE_CELL_COUNT_MAX = 2 * ONE_SMG_EFFECTIVE_CELL_COUNT_MAX
ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S = (
    ONE_SMG_EFFECTIVE_CELL_COUNT_MAX * PL_S_CELL_TO_UL_MIN
)
PAIRED_SMG_SCALE_MAX_UL_MIN_PER_PL_S = (
    PAIRED_SMG_EFFECTIVE_CELL_COUNT_MAX * PL_S_CELL_TO_UL_MIN
)


@dataclass(frozen=True)
class SharedScaleInterval:
    """Intersection of all WT observation-scale constraints for one row."""

    lower_uL_min_per_pL_s: float
    upper_uL_min_per_pL_s: float
    sample_count: int
    minimum_flow_pL_s: float
    maximum_flow_pL_s: float

    def __post_init__(self) -> None:
        values = (
            self.lower_uL_min_per_pL_s,
            self.upper_uL_min_per_pL_s,
            self.minimum_flow_pL_s,
            self.maximum_flow_pL_s,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise ValueError("shared-scale bounds and flows must be finite and positive")
        if int(self.sample_count) != self.sample_count or self.sample_count <= 0:
            raise ValueError("shared-scale sample_count must be a positive integer")
        if self.maximum_flow_pL_s < self.minimum_flow_pL_s:
            raise ValueError("shared-scale flow extrema must be ordered")

    @property
    def compatible(self) -> bool:
        return bool(
            math.isfinite(self.lower_uL_min_per_pL_s)
            and math.isfinite(self.upper_uL_min_per_pL_s)
            and self.lower_uL_min_per_pL_s <= self.upper_uL_min_per_pL_s
        )

    @property
    def maximum_to_minimum_flow_ratio(self) -> float:
        return float(self.maximum_flow_pL_s / self.minimum_flow_pL_s)


@dataclass(frozen=True)
class GeometryScaleAssessment:
    """One-SMG or explicitly labelled paired-SMG interval intersection."""

    convention: str
    cell_count_max: int
    scale_max_uL_min_per_pL_s: float
    unclipped_lower_uL_min_per_pL_s: float
    unclipped_upper_uL_min_per_pL_s: float
    feasible_lower_uL_min_per_pL_s: float
    feasible_upper_uL_min_per_pL_s: float
    interval_nonempty: bool
    geometry_gate_pass: bool
    implied_cell_count_at_feasible_lower: float | None


def _positive_finite_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    vector = np.asarray(values, dtype=float).reshape(-1)
    if vector.size == 0:
        raise ValueError(f"{name} must contain at least one value")
    if np.any(~np.isfinite(vector)) or np.any(vector <= 0.0):
        raise ValueError(f"{name} must be finite and strictly positive")
    return vector


def shared_observation_scale_interval(
    flow_samples_pL_s: ArrayLike,
    *,
    lower_uL_min: float = WT_FLOW_LOWER_UL_MIN,
    upper_uL_min: float = WT_FLOW_UPPER_UL_MIN,
) -> SharedScaleInterval:
    """Intersect ``[lower/q, upper/q]`` over all declared WT samples."""

    if (
        not math.isfinite(lower_uL_min)
        or not math.isfinite(upper_uL_min)
        or lower_uL_min <= 0.0
        or upper_uL_min < lower_uL_min
    ):
        raise ValueError("WT flow bounds must be finite, positive, and ordered")
    flows = _positive_finite_vector(flow_samples_pL_s, "WT flow samples")
    lower = float(np.max(lower_uL_min / flows))
    upper = float(np.min(upper_uL_min / flows))
    return SharedScaleInterval(
        lower_uL_min_per_pL_s=lower,
        upper_uL_min_per_pL_s=upper,
        sample_count=int(flows.size),
        minimum_flow_pL_s=float(np.min(flows)),
        maximum_flow_pL_s=float(np.max(flows)),
    )


def assess_geometry_scale(
    interval: SharedScaleInterval,
    *,
    paired_smg_sensitivity: bool = False,
) -> GeometryScaleAssessment:
    """Intersect a shared WT scale interval with the independent gland bound."""

    if paired_smg_sensitivity:
        convention = "PAIRED_SMG_SENSITIVITY_NOT_PRODUCTION"
        cell_count_max = PAIRED_SMG_EFFECTIVE_CELL_COUNT_MAX
        scale_max = PAIRED_SMG_SCALE_MAX_UL_MIN_PER_PL_S
    else:
        convention = "ONE_IPSILATERAL_SMG_PRODUCTION"
        cell_count_max = ONE_SMG_EFFECTIVE_CELL_COUNT_MAX
        scale_max = ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
    lower = float(interval.lower_uL_min_per_pL_s)
    upper = float(min(interval.upper_uL_min_per_pL_s, scale_max))
    interval_nonempty = bool(interval.compatible)
    gate = bool(interval_nonempty and lower <= upper)
    return GeometryScaleAssessment(
        convention=convention,
        cell_count_max=cell_count_max,
        scale_max_uL_min_per_pL_s=float(scale_max),
        unclipped_lower_uL_min_per_pL_s=lower,
        unclipped_upper_uL_min_per_pL_s=float(
            interval.upper_uL_min_per_pL_s
        ),
        feasible_lower_uL_min_per_pL_s=lower,
        feasible_upper_uL_min_per_pL_s=upper,
        interval_nonempty=interval_nonempty,
        geometry_gate_pass=gate,
        implied_cell_count_at_feasible_lower=(
            float(lower / PL_S_CELL_TO_UL_MIN) if gate else None
        ),
    )


def sustainment_ratio(
    time_s: ArrayLike,
    flow_pL_s: ArrayLike,
    *,
    early_s: float = 60.0,
    late_s: float = 600.0,
) -> float:
    """Return the scale-free late/early WT flow ratio."""

    time = np.asarray(time_s, dtype=float).reshape(-1)
    flow = _positive_finite_vector(flow_pL_s, "WT trajectory flow")
    if time.shape != flow.shape or np.any(~np.isfinite(time)) or np.any(np.diff(time) <= 0.0):
        raise ValueError("time and flow must be aligned on a finite increasing grid")
    if not time[0] <= early_s < late_s <= time[-1]:
        raise ValueError("sustainment landmarks lie outside the trajectory")
    early = float(np.interp(early_s, time, flow))
    late = float(np.interp(late_s, time, flow))
    return float(late / early)


def sustainment_gate(
    ratio: float,
    *,
    lower: float = 0.8,
    upper: float = 1.2,
) -> bool:
    if not all(math.isfinite(value) for value in (ratio, lower, upper)):
        return False
    if lower <= 0.0 or upper < lower:
        raise ValueError("sustainment bounds must be positive and ordered")
    return bool(lower <= ratio <= upper)


def relative_difference(
    left: ArrayLike,
    right: ArrayLike,
    *,
    floor: float = 1.0e-12,
) -> float:
    """Maximum symmetric relative difference with a declared numeric floor."""

    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("relative-difference arrays must be finite and aligned")
    if not math.isfinite(floor) or floor <= 0.0:
        raise ValueError("relative-difference floor must be finite and positive")
    scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), floor)
    return float(np.max(np.abs(a - b) / scale))


def paired_arm_total_ratio(
    cch_only_cumulative_pL: float,
    costim_cumulative_pL: float,
) -> float:
    """Return the same-row CCh+IPR/CCh-only WT 10-min total ratio."""

    values = (float(cch_only_cumulative_pL), float(costim_cumulative_pL))
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("paired-arm cumulative flows must be finite and positive")
    return float(values[1] / values[0])


def regulatory_member_pairs(
    members: Iterable[object],
) -> tuple[tuple[str, str], ...]:
    """Return exact R2/R3 isomorphism pairs by gain and timing labels."""

    registry: dict[tuple[str, str, str], str] = {}
    for member in members:
        family = str(getattr(member, "family"))
        if family not in {"R2", "R3"}:
            continue
        key = (
            family,
            str(getattr(member, "gain_label")),
            str(getattr(member, "timing_label")),
        )
        if key in registry:
            raise ValueError(f"duplicate regulatory axis {key}")
        registry[key] = str(getattr(member, "member_id"))
    pairs: list[tuple[str, str]] = []
    axes = sorted({(gain, timing) for _, gain, timing in registry})
    for gain, timing in axes:
        r2 = registry.get(("R2", gain, timing))
        r3 = registry.get(("R3", gain, timing))
        if r2 is None or r3 is None:
            raise ValueError(f"incomplete R2/R3 pair for {gain}/{timing}")
        pairs.append((r2, r3))
    if len(pairs) != 6:
        raise ValueError(f"expected six R2/R3 isomorphism pairs, found {len(pairs)}")
    return tuple(pairs)


def validate_regulatory_ensemble(members: Iterable[object]) -> tuple[object, ...]:
    """Fail closed unless the exact frozen 2/6/6/6 panel is present."""

    panel = tuple(members)
    if len(panel) != 20:
        raise ValueError(f"regulatory panel must contain 20 members, found {len(panel)}")
    identifiers = [str(getattr(member, "member_id")) for member in panel]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("regulatory member IDs must be unique")
    counts = Counter(str(getattr(member, "family")) for member in panel)
    if counts != Counter({"R0": 2, "R1": 6, "R2": 6, "R3": 6}):
        raise ValueError(f"regulatory family composition is not 2/6/6/6: {counts}")
    gains = {"G125_P21_PROSE", "G170_P21_FIGURE_VISUAL"}
    timings = {"TFAST", "TREFERENCE", "TSLOW"}
    for family in ("R1", "R2", "R3"):
        axes = {
            (
                str(getattr(member, "gain_label")),
                str(getattr(member, "timing_label")),
            )
            for member in panel
            if str(getattr(member, "family")) == family
        }
        expected = {(gain, timing) for gain in gains for timing in timings}
        if axes != expected:
            raise ValueError(f"{family} does not span the exact gain/timing grid")
    r0_axes = {
        (
            str(getattr(member, "gain_label")),
            str(getattr(member, "timing_label")),
        )
        for member in panel
        if str(getattr(member, "family")) == "R0"
    }
    if r0_axes != {(gain, "STATIC") for gain in gains}:
        raise ValueError("R0 must contain exactly one static control per gain")
    regulatory_member_pairs(panel)
    return panel


def geometry_contract_payload() -> Mapping[str, object]:
    """Machine-readable provenance for every absolute geometry deduction."""

    derived_acinar_volume_pL = (
        ONE_SMG_MASS_CEILING_MG
        * TRUE_ACINAR_VOLUME_FRACTION
        * 1.0e6
    )
    return {
        "production_convention": "ONE_IPSILATERAL_SMG",
        "one_smg_mass_ceiling_mg": ONE_SMG_MASS_CEILING_MG,
        "tissue_density_conversion_mg_per_uL": 1.0,
        "true_acinar_volume_fraction": TRUE_ACINAR_VOLUME_FRACTION,
        "derived_true_acinar_volume_pL": derived_acinar_volume_pL,
        "frozen_cell_volume_pL": FROZEN_CELL_VOLUME_PL,
        "unrounded_cell_count_ceiling": (
            derived_acinar_volume_pL / FROZEN_CELL_VOLUME_PL
        ),
        "one_smg_effective_cell_count_max": ONE_SMG_EFFECTIVE_CELL_COUNT_MAX,
        "one_smg_scale_max_uL_min_per_pL_s": (
            ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
        ),
        "minimum_per_cell_flow_for_9_uL_min_pL_s": (
            WT_FLOW_LOWER_UL_MIN / ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
        ),
        "minimum_per_cell_flow_for_10_uL_min_pL_s": (
            WT_FLOW_UPPER_UL_MIN / ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
        ),
        "minimum_flow_per_frozen_cell_volume_for_9_uL_min_s_inv": (
            WT_FLOW_LOWER_UL_MIN
            / ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
            / FROZEN_CELL_VOLUME_PL
        ),
        "minimum_flow_per_frozen_cell_volume_for_10_uL_min_s_inv": (
            WT_FLOW_UPPER_UL_MIN
            / ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
            / FROZEN_CELL_VOLUME_PL
        ),
        "paired_smg_effective_cell_count_max_sensitivity": (
            PAIRED_SMG_EFFECTIVE_CELL_COUNT_MAX
        ),
        "paired_smg_scale_max_uL_min_per_pL_s_sensitivity": (
            PAIRED_SMG_SCALE_MAX_UL_MIN_PER_PL_S
        ),
        "unit_conversion_uL_min_per_pL_s_per_cell": PL_S_CELL_TO_UL_MIN,
        "provenance": {
            "acinar_fraction": "PRIMARY_MEASUREMENT:KONDO_2015",
            "one_ipsilateral_gland_protocol": "PRIMARY_MEASUREMENT:KONDO_2019",
            "100_mg_ceiling": "NEW_MODELING_DECISION:GENEROUS_UPPER_BOUND",
            "tissue_density_conversion": "NEW_MODELING_DECISION:1_MG_PER_UL",
            "cell_count_and_scale": "DERIVED_CONSTRAINT",
        },
    }


__all__ = (
    "FROZEN_CELL_VOLUME_PL",
    "GeometryScaleAssessment",
    "ONE_SMG_EFFECTIVE_CELL_COUNT_MAX",
    "ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S",
    "PAIRED_SMG_EFFECTIVE_CELL_COUNT_MAX",
    "PAIRED_SMG_SCALE_MAX_UL_MIN_PER_PL_S",
    "PL_S_CELL_TO_UL_MIN",
    "SharedScaleInterval",
    "TRUE_ACINAR_VOLUME_FRACTION",
    "WT_FLOW_LOWER_UL_MIN",
    "WT_FLOW_UPPER_UL_MIN",
    "assess_geometry_scale",
    "geometry_contract_payload",
    "paired_arm_total_ratio",
    "regulatory_member_pairs",
    "relative_difference",
    "shared_observation_scale_interval",
    "sustainment_gate",
    "sustainment_ratio",
    "validate_regulatory_ensemble",
)
