"""Gate-first secretion validation metrics.

No experimental target value appears in this module.  A knockout/WT ratio is
defined only after an independently calibrated WT capacity and the WT
physiology gate both pass.  Physical early/late timing additionally requires
a certified mapping from model time to minutes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike


class PredictiveGateError(RuntimeError):
    """Raised when a held-out metric is requested before its upstream gates."""


@dataclass(frozen=True)
class PredictiveEligibility:
    """Upstream facts required before calculating secretion validation metrics."""

    wt_capacity_valid: bool
    wt_physiology_gate_pass: bool
    chassis_independently_calibrated: bool
    physical_time_certified: bool = False

    @property
    def can_form_ratio(self) -> bool:
        return bool(
            self.wt_capacity_valid
            and self.wt_physiology_gate_pass
            and self.chassis_independently_calibrated
        )

    @property
    def can_score_physical_timing(self) -> bool:
        return bool(self.can_form_ratio and self.physical_time_certified)

    def blocking_reasons(self) -> tuple[str, ...]:
        reasons: list[str] = []
        if not self.wt_capacity_valid:
            reasons.append("no valid WT capacity calibration")
        if not self.wt_physiology_gate_pass:
            reasons.append("WT physiology gate failed")
        if not self.chassis_independently_calibrated:
            reasons.append("chassis extension is not independently calibrated")
        return tuple(reasons)


@dataclass(frozen=True)
class SecretionRatios:
    """Target-free ratios under one common, already frozen protocol."""

    endpoint_ko_wt: float
    integral_ko_wt: float


def secretion_ratios(
    time: ArrayLike,
    wt_flow: ArrayLike,
    knockout_flow: ArrayLike,
    eligibility: PredictiveEligibility,
) -> SecretionRatios:
    """Return endpoint and integral ratios only after upstream validation.

    The function does not compare either value with an experimental target.
    Model selection and target joining must remain separate operations.
    """

    if not eligibility.can_form_ratio:
        reasons = "; ".join(eligibility.blocking_reasons())
        raise PredictiveGateError(f"predictive secretion ratio blocked: {reasons}")
    t = np.asarray(time, dtype=float)
    wt = np.asarray(wt_flow, dtype=float)
    ko = np.asarray(knockout_flow, dtype=float)
    if t.ndim != 1 or wt.shape != t.shape or ko.shape != t.shape or t.size < 2:
        raise ValueError("time and flow arrays must be one-dimensional and aligned")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(wt)) or not np.all(np.isfinite(ko)):
        raise ValueError("time and flow arrays must be finite")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("time must be strictly increasing")
    wt_integral = float(np.trapezoid(wt, t))
    if wt[-1] == 0.0 or wt_integral == 0.0:
        raise ZeroDivisionError("WT endpoint and integral must be nonzero")
    return SecretionRatios(
        endpoint_ko_wt=float(ko[-1] / wt[-1]),
        integral_ko_wt=float(np.trapezoid(ko, t) / wt_integral),
    )


def require_physical_timing(eligibility: PredictiveEligibility) -> None:
    """Reject early/late scoring when model time has no minute conversion."""

    if not eligibility.can_score_physical_timing:
        reasons = list(eligibility.blocking_reasons())
        if not eligibility.physical_time_certified:
            reasons.append("physical time mapping is uncertified")
        raise PredictiveGateError(
            "physical early/late validation blocked: " + "; ".join(reasons)
        )


__all__ = (
    "PredictiveEligibility",
    "PredictiveGateError",
    "SecretionRatios",
    "require_physical_timing",
    "secretion_ratios",
)
