"""Explicit bicarbonate/carbonate chemistry for the SR4 hypothesis."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class CarbonateChemistry:
    """Second dissociation of carbonic acid at a declared temperature.

    ``pka2`` is an effective activity-based value.  The default 10.33 is a
    25-degree-C reference used only to make SR4 chemically explicit; the
    primary AE4 experiments did not identify intracellular carbonate or its
    activity coefficient.  Any SR4 result must therefore retain this as a new
    modeling assumption rather than a fitted bicarbonate surrogate.
    """

    pka2: float = 10.33
    temperature_c: float = 25.0

    def carbonate_from_bicarbonate(self, hco3_mM: float, ph: float) -> float:
        if not math.isfinite(hco3_mM) or hco3_mM <= 0.0:
            raise ValueError("bicarbonate must be finite and positive")
        if not math.isfinite(ph):
            raise ValueError("pH must be finite")
        return hco3_mM * 10.0 ** (ph - self.pka2)

    def speciate_total_base(self, total_hco3_plus_co3_mM: float, ph: float) -> tuple[float, float]:
        """Return HCO3 and CO3 from their two-species analytical total."""

        if not math.isfinite(total_hco3_plus_co3_mM) or total_hco3_plus_co3_mM <= 0.0:
            raise ValueError("analytical carbonate base must be positive")
        ratio = 10.0 ** (ph - self.pka2)
        hco3 = total_hco3_plus_co3_mM / (1.0 + ratio)
        return hco3, total_hco3_plus_co3_mM - hco3

    def equilibrium_residual(self, hco3_mM: float, co3_mM: float, ph: float) -> float:
        """Log residual of ``CO3/HCO3 = 10^(pH-pKa2)``."""

        if hco3_mM <= 0.0 or co3_mM <= 0.0:
            raise ValueError("carbonate species must be positive")
        return math.log(co3_mM / hco3_mM) - math.log(10.0) * (ph - self.pka2)


DEFAULT_CARBONATE_CHEMISTRY = CarbonateChemistry()
