"""Source-fixed Palk/Benjamin two-state NKCC1, Vera-Siguenza 2018 Eq. 17.

Concentrations are mM. The historical M-based fourth-order coefficients
are multiplied by (1e-3)**4 = 1e-12. The shape is dimensionless in the
contract's rate convention; alpha_eff is a whole-cell cycle rate in fmol/s,
derived once from accepted WT REST, not the historical membrane density.
No voltage term, dynamic state, clipping, or stimulation law is added.
"""
from dataclasses import dataclass
import math

GENERIC_NKCC1 = "generic_tanh"
PALK_NKCC1 = "palk_benjamin_2010_eq17"
A1 = 157.5
A2_MM = 2.0096e-5
A3 = 1.0306
A4_MM = 1.3852e-6


def intracellular_product_mM4(na_i_mM: float, k_i_mM: float, cl_i_mM: float) -> float:
    if not all(math.isfinite(v) and v > 0 for v in (na_i_mM, k_i_mM, cl_i_mM)):
        raise ValueError("NKCC1 concentrations must be finite and positive")
    return na_i_mM * k_i_mM * cl_i_mM**2


def palk_shape_factor(na_i_mM: float, k_i_mM: float, cl_i_mM: float) -> float:
    x = intracellular_product_mM4(na_i_mM, k_i_mM, cl_i_mM)
    return (A1 - A2_MM * x) / (A3 + A4_MM * x)


def palk_cycle_flux_fmol_s(na_i_mM: float, k_i_mM: float, cl_i_mM: float,
                          *, alpha_eff_fmol_s: float, activity_multiplier: float) -> float:
    """Signed inward 1 Na : 1 K : 2 Cl cycles; the multiplier scales both ways."""
    if not math.isfinite(alpha_eff_fmol_s) or alpha_eff_fmol_s <= 0:
        raise ValueError("Palk alpha_eff must be finite and positive")
    if not math.isfinite(activity_multiplier) or activity_multiplier < 0:
        raise ValueError("NKCC1 activity multiplier must be finite and nonnegative")
    return alpha_eff_fmol_s * activity_multiplier * palk_shape_factor(
        na_i_mM, k_i_mM, cl_i_mM)


@dataclass(frozen=True)
class Nkcc1Kinetics:
    """Explicit core-law selection, separate from the inherited activity law.

    The generic default preserves historical parameter payloads and hashes.
    Task 39 must explicitly provide the algebraically derived Palk scale.
    """
    law: str = GENERIC_NKCC1
    alpha_eff_fmol_s: float | None = None

    def __post_init__(self):
        if self.law == GENERIC_NKCC1:
            if self.alpha_eff_fmol_s is not None:
                raise ValueError("The generic comparator does not use Palk alpha_eff")
        elif self.law == PALK_NKCC1:
            if self.alpha_eff_fmol_s is None or not math.isfinite(self.alpha_eff_fmol_s) or self.alpha_eff_fmol_s <= 0:
                raise ValueError("Palk selection requires a positive frozen alpha_eff")
        else:
            raise ValueError(f"Unknown NKCC1 core law: {self.law!r}")
