"""A transparent first-order steady-state pump--leak diagnostic.

This module deliberately does *not* claim to reproduce the 2018 numerical model.
The publication leaves several quantities needed for that task unresolved.  The
model below isolates the part that is independently recoverable: AE2 and AE4
have distinct published stoichiometric signatures, and their local effect on a
regular steady state is governed by the implicit-function formula.

States and activities are dimensionless ratios to the published resting values.
The positive-definite restoring matrix and activity lever scales are declared
diagnostic choices.  They permit executable rank/equivalence checks without
being confused with a calibrated physiological model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Iterable, Sequence

import numpy as np


STATE_NAMES = ("Na_i", "K_i", "Cl_i", "HCO3_i")
PARAMETER_NAMES = ("G2", "G4")
OBSERVATION_NAMES = ("Q_star", "Na_i", "K_i", "Cl_i", "HCO3_i", "volume_star")


def printed_activity_factor_check() -> dict[str, float]:
    """Evaluate printed AE2/AE4 laws at the published resting concentrations.

    The two outputs are deliberately not combined or compared: the publication
    gives incompatible activity/density units.  This calculation exposes that
    incompatibility and provides a regression check on transcription of Eqs.
    (25) and (31), not a flux calibration.
    """

    cl_e, cl_i = 102.6, 50.0
    hco3_e, hco3_i = 42.9, 12.1
    na_e, na_i = 140.2, 25.0
    k_e, k_i = 5.3, 120.0
    k_cl, k_b = 5.6, 1.0e4
    g2 = 0.01807
    g4 = 0.66
    k_plus, k_minus = 1.92e-2, 1.3e-5
    phi2 = (
        cl_e / (cl_e + k_cl) * hco3_i / (hco3_i + k_b)
        - cl_i / (cl_i + k_cl) * hco3_e / (hco3_e + k_b)
    )
    phi4 = (
        k_plus * cl_e * hco3_i**2 * (na_i + k_i)
        - k_minus * cl_i * hco3_e**2 * (na_e + k_e)
    )
    return {
        "AE2_dimensionless_bracket": float(phi2),
        "printed_G2_times_bracket": float(g2 * phi2),
        "AE4_printed_kinetic_bracket": float(phi4),
        "printed_G4_times_bracket": float(g4 * phi4),
    }


@dataclass(frozen=True)
class LinearizedPumpLeakModel:
    """First-order reduction ``F(x; theta)=0`` around the published baseline.

    ``F = A (x - 1) - S (theta - 1)``.  Here ``x`` contains concentration
    ratios, ``theta`` contains activity multipliers, ``A`` is a regular
    effective restoring Jacobian, and ``S`` contains normalized AE2/AE4
    stoichiometric signatures.  Thus ``dx/dtheta = A^{-1} S`` exactly.
    """

    baseline_concentrations_mM: np.ndarray = field(
        default_factory=lambda: np.array([25.0, 120.0, 50.0, 12.1], dtype=float)
    )
    restoring_matrix: np.ndarray = field(
        default_factory=lambda: np.array(
            [
                [1.00, 0.08, -0.05, 0.02],
                [0.08, 1.00, 0.02, 0.03],
                [-0.05, 0.02, 1.00, -0.10],
                [0.02, 0.03, -0.10, 1.00],
            ],
            dtype=float,
        )
    )
    activity_levers: np.ndarray = field(
        default_factory=lambda: np.array([0.08, 0.12], dtype=float)
    )
    q_weights: np.ndarray = field(
        default_factory=lambda: np.array([0.10, 0.10, 0.70, 0.10], dtype=float)
    )
    volume_weights: np.ndarray = field(
        default_factory=lambda: np.full(4, 0.25, dtype=float)
    )

    def __post_init__(self) -> None:
        if self.baseline_concentrations_mM.shape != (4,):
            raise ValueError("baseline_concentrations_mM must have shape (4,)")
        if self.restoring_matrix.shape != (4, 4):
            raise ValueError("restoring_matrix must have shape (4, 4)")
        if np.min(np.linalg.eigvalsh(self.restoring_matrix)) <= 0:
            raise ValueError("restoring_matrix must be positive definite")
        if self.activity_levers.shape != (2,) or np.any(self.activity_levers <= 0):
            raise ValueError("activity_levers must contain two positive values")

    @property
    def stoichiometric_matrix(self) -> np.ndarray:
        """Return normalized, scaled AE2/AE4 intracellular signatures.

        Row order is Na, K, Cl, HCO3.  AE2 imports one Cl and exports one
        bicarbonate.  AE4 imports one Cl and exports two bicarbonates plus one
        cation, partitioned using the published baseline Na/K mole fractions.
        Column normalization removes arbitrary cycle-count magnitude before the
        declared diagnostic lever is applied.
        """

        na, k = self.baseline_concentrations_mM[:2]
        rho_na = na / (na + k)
        rho_k = k / (na + k)
        raw = np.array(
            [
                [0.0, -rho_na],
                [0.0, -rho_k],
                [1.0, 1.0],
                [-1.0, -2.0],
            ],
            dtype=float,
        )
        return raw / np.linalg.norm(raw, axis=0) * self.activity_levers

    @property
    def state_sensitivity(self) -> np.ndarray:
        """Exact implicit sensitivity ``-F_x^{-1} F_theta``."""

        return np.linalg.solve(self.restoring_matrix, self.stoichiometric_matrix)

    def residual(self, x: Sequence[float], theta: Sequence[float]) -> np.ndarray:
        x_arr = np.asarray(x, dtype=float)
        theta_arr = np.asarray(theta, dtype=float)
        return self.restoring_matrix @ (x_arr - 1.0) - self.stoichiometric_matrix @ (
            theta_arr - 1.0
        )

    def steady_state(self, theta: Sequence[float]) -> np.ndarray:
        theta_arr = np.asarray(theta, dtype=float)
        if theta_arr.shape != (2,):
            raise ValueError("theta must have shape (2,)")
        return 1.0 + self.state_sensitivity @ (theta_arr - 1.0)

    def finite_difference_sensitivity(
        self, theta: Sequence[float] = (1.0, 1.0), step: float = 1e-6
    ) -> np.ndarray:
        theta_arr = np.asarray(theta, dtype=float)
        answer = np.empty((4, 2), dtype=float)
        for j in range(2):
            delta = np.zeros(2)
            delta[j] = step
            answer[:, j] = (
                self.steady_state(theta_arr + delta)
                - self.steady_state(theta_arr - delta)
            ) / (2.0 * step)
        return answer

    @property
    def observation_rows(self) -> dict[str, np.ndarray]:
        eye = np.eye(4)
        return {
            "Q_star": self.q_weights.copy(),
            "Na_i": eye[0],
            "K_i": eye[1],
            "Cl_i": eye[2],
            "HCO3_i": eye[3],
            "volume_star": self.volume_weights.copy(),
        }

    def observation_values(self, theta: Sequence[float]) -> dict[str, float]:
        x = self.steady_state(theta)
        rows = self.observation_rows
        values = {name: float(1.0 + row @ (x - 1.0)) for name, row in rows.items()}
        for i, state in enumerate(STATE_NAMES):
            values[f"{state}_mM"] = float(x[i] * self.baseline_concentrations_mM[i])
        return values

    def observation_jacobian(self, panel: Iterable[str]) -> np.ndarray:
        names = tuple(panel)
        if not names:
            raise ValueError("panel must contain at least one observation")
        rows = self.observation_rows
        h_x = np.vstack([rows[name] for name in names])
        return h_x @ self.state_sensitivity

    def panel_metrics(self, panel: Iterable[str], tolerance: float = 1e-12) -> dict[str, float | int | str]:
        names = tuple(panel)
        jacobian = self.observation_jacobian(names)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        rank = int(np.sum(singular > tolerance))
        sigma_min = float(singular[-1])
        sigma_max = float(singular[0])
        condition = float(sigma_max / sigma_min) if rank == 2 and sigma_min > 0 else float("inf")
        determinant = float(np.linalg.det(jacobian)) if jacobian.shape == (2, 2) else float("nan")
        return {
            "panel": "+".join(names),
            "n_observations": len(names),
            "rank": rank,
            "sigma_min": sigma_min,
            "sigma_max": sigma_max,
            "condition_number": condition,
            "determinant": determinant,
        }

    def all_panel_metrics(self) -> list[dict[str, float | int | str]]:
        answer: list[dict[str, float | int | str]] = []
        for size in (1, 2):
            for panel in combinations(OBSERVATION_NAMES, size):
                answer.append(self.panel_metrics(panel))
        return answer

    def q_equivalence_segment(
        self,
        reference_theta: Sequence[float] = (0.5, 0.5),
        lower: float = 0.0,
        upper: float = 1.0,
        n_points: int = 201,
    ) -> np.ndarray:
        """Points on an interior affine Q* equivalence set in the square."""

        reference = np.asarray(reference_theta, dtype=float)
        gradient = self.observation_jacobian(("Q_star",))[0]
        if abs(gradient[1]) < 1e-15:
            raise ValueError("Q* gradient cannot parameterize G4")
        g2 = np.linspace(lower, upper, n_points)
        g4 = reference[1] - gradient[0] * (g2 - reference[0]) / gradient[1]
        points = np.column_stack([g2, g4])
        mask = (g4 >= lower) & (g4 <= upper)
        return points[mask]

    def row_sign_invariance_error(self, signs: Sequence[float]) -> float:
        """Verify invariance when complete residual equations change sign.

        A consistent row sign change multiplies both ``F_x`` and ``F_theta`` by
        the same diagonal matrix and cannot change the implicit sensitivity.
        """

        sign_arr = np.asarray(signs, dtype=float)
        if sign_arr.shape != (4,) or np.any(np.abs(sign_arr) != 1):
            raise ValueError("signs must contain four values in {-1, +1}")
        d = np.diag(sign_arr)
        transformed = -np.linalg.solve(
            d @ self.restoring_matrix, d @ (-self.stoichiometric_matrix)
        )
        return float(np.max(np.abs(transformed - self.state_sensitivity)))
