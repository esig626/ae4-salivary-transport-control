"""State definitions for the Task 13B modern whole-cell model.

The dynamic core is stored as *amounts* (fmol) and compartment volumes (pL).
Because ``1 mM * 1 pL == 1 fmol``, division by volume gives concentrations in
mM without a hidden scale factor.  This choice makes dilution an automatic
consequence of volume dynamics rather than a separately coded source term.

Intracellular and luminal pH are algebraic observables reconstructed from the
total-inorganic-carbon and total-alkalinity states in :mod:`acid_base`.
Membrane potentials are likewise algebraic variables obtained from current
closure.  Those reductions are documented modeling decisions; neither is a
return to the inherited seven-state chassis.

Provenance: ``DERIVED_CONSTRAINT`` (amount/concentration conversion and
conservation bookkeeping) and ``NEW_MODELING_DECISION`` (chosen coordinates).
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


CORE_STATE_NAMES: tuple[str, ...] = (
    "na_i_fmol",
    "k_i_fmol",
    "cl_i_fmol",
    "tic_i_fmol",
    "alk_i_fmol",
    "volume_i_pL",
    "na_l_fmol",
    "k_l_fmol",
    "cl_l_fmol",
    "tic_l_fmol",
    "alk_l_fmol",
    "volume_l_pL",
)

CELL_AMOUNT_NAMES: tuple[str, ...] = CORE_STATE_NAMES[0:5]
LUMEN_AMOUNT_NAMES: tuple[str, ...] = CORE_STATE_NAMES[6:11]


@dataclass(frozen=True)
class CompartmentAmounts:
    """Conserved or externally exchanged quantities in one compartment."""

    na_fmol: float
    k_fmol: float
    cl_fmol: float
    tic_fmol: float
    alkalinity_fmol: float
    volume_pL: float

    def concentrations(self) -> "CompartmentConcentrations":
        if not math.isfinite(self.volume_pL) or self.volume_pL <= 0.0:
            raise ValueError("compartment volume must be finite and positive")
        return CompartmentConcentrations(
            na_mM=self.na_fmol / self.volume_pL,
            k_mM=self.k_fmol / self.volume_pL,
            cl_mM=self.cl_fmol / self.volume_pL,
            tic_mM=self.tic_fmol / self.volume_pL,
            alkalinity_mM=self.alkalinity_fmol / self.volume_pL,
        )


@dataclass(frozen=True)
class CompartmentConcentrations:
    """Concentrations derived from amount states (all concentrations in mM)."""

    na_mM: float
    k_mM: float
    cl_mM: float
    tic_mM: float
    alkalinity_mM: float


@dataclass(frozen=True)
class WholeCellState:
    """Decoded core state plus optional regulatory coordinates."""

    cell: CompartmentAmounts
    lumen: CompartmentAmounts
    regulation: tuple[float, ...] = ()


@dataclass(frozen=True)
class StateLayout:
    """Stable vector layout for the conserved core and a regulatory suffix."""

    regulatory_names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        names = self.names
        if len(names) != len(set(names)):
            raise ValueError("state names must be unique")
        if any(name in CORE_STATE_NAMES for name in self.regulatory_names):
            raise ValueError("regulatory names collide with core state names")

    @property
    def names(self) -> tuple[str, ...]:
        return CORE_STATE_NAMES + self.regulatory_names

    @property
    def size(self) -> int:
        return len(self.names)

    @property
    def core_size(self) -> int:
        return len(CORE_STATE_NAMES)

    def index(self, name: str) -> int:
        try:
            return self.names.index(name)
        except ValueError as exc:
            raise KeyError(name) from exc

    def validate(self, vector: ArrayLike, *, require_positive_core: bool = True) -> NDArray[np.float64]:
        state = np.asarray(vector, dtype=float)
        if state.shape != (self.size,):
            raise ValueError(f"state must have shape ({self.size},), got {state.shape}")
        if not np.all(np.isfinite(state)):
            raise ValueError("state contains a non-finite value")
        if require_positive_core and np.any(state[: self.core_size] <= 0.0):
            bad = [
                name
                for name, value in zip(CORE_STATE_NAMES, state[: self.core_size])
                if value <= 0.0
            ]
            raise ValueError(f"core amount/volume states must be positive: {bad}")
        return state

    def decode(self, vector: ArrayLike, *, require_positive_core: bool = True) -> WholeCellState:
        y = self.validate(vector, require_positive_core=require_positive_core)
        cell = CompartmentAmounts(
            na_fmol=float(y[0]),
            k_fmol=float(y[1]),
            cl_fmol=float(y[2]),
            tic_fmol=float(y[3]),
            alkalinity_fmol=float(y[4]),
            volume_pL=float(y[5]),
        )
        lumen = CompartmentAmounts(
            na_fmol=float(y[6]),
            k_fmol=float(y[7]),
            cl_fmol=float(y[8]),
            tic_fmol=float(y[9]),
            alkalinity_fmol=float(y[10]),
            volume_pL=float(y[11]),
        )
        return WholeCellState(
            cell=cell,
            lumen=lumen,
            regulation=tuple(float(value) for value in y[self.core_size :]),
        )

    def encode(self, state: WholeCellState) -> NDArray[np.float64]:
        if len(state.regulation) != len(self.regulatory_names):
            raise ValueError(
                "regulatory coordinate count does not match the state layout"
            )
        vector = np.asarray(
            (
                state.cell.na_fmol,
                state.cell.k_fmol,
                state.cell.cl_fmol,
                state.cell.tic_fmol,
                state.cell.alkalinity_fmol,
                state.cell.volume_pL,
                state.lumen.na_fmol,
                state.lumen.k_fmol,
                state.lumen.cl_fmol,
                state.lumen.tic_fmol,
                state.lumen.alkalinity_fmol,
                state.lumen.volume_pL,
                *state.regulation,
            ),
            dtype=float,
        )
        return self.validate(vector)

    def mapping(self, vector: ArrayLike) -> Mapping[str, float]:
        y = self.validate(vector, require_positive_core=False)
        return {name: float(value) for name, value in zip(self.names, y)}


def compartment_from_concentrations(
    *,
    na_mM: float,
    k_mM: float,
    cl_mM: float,
    tic_mM: float,
    alkalinity_mM: float,
    volume_pL: float,
) -> CompartmentAmounts:
    """Construct amount states from concentrations without hidden conversion."""

    values = (na_mM, k_mM, cl_mM, tic_mM, alkalinity_mM, volume_pL)
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("concentrations and volume must be finite and positive")
    return CompartmentAmounts(
        na_fmol=na_mM * volume_pL,
        k_fmol=k_mM * volume_pL,
        cl_fmol=cl_mM * volume_pL,
        tic_fmol=tic_mM * volume_pL,
        alkalinity_fmol=alkalinity_mM * volume_pL,
        volume_pL=volume_pL,
    )


def vector_from_mapping(
    values: Mapping[str, float], *, regulatory_names: Sequence[str] = ()
) -> NDArray[np.float64]:
    """Build a state vector from an explicitly named mapping."""

    layout = StateLayout(tuple(regulatory_names))
    missing = [name for name in layout.names if name not in values]
    if missing:
        raise KeyError(f"missing state values: {missing}")
    return layout.validate([values[name] for name in layout.names])


def amount_charge_equivalents_fmol(
    compartment: CompartmentAmounts, *, fixed_anion_equivalents_fmol: float = 0.0
) -> float:
    """Net bulk charge in fmol equivalents.

    Total alkalinity is defined as bicarbonate + twice carbonate + deprotonated
    buffer + hydroxide - proton.  Therefore the modeled mobile charge is
    exactly ``Na + K - Cl - alkalinity``.  A fixed impermeant-anion term is
    included for the intracellular compartment.
    """

    return (
        compartment.na_fmol
        + compartment.k_fmol
        - compartment.cl_fmol
        - compartment.alkalinity_fmol
        - fixed_anion_equivalents_fmol
    )


__all__ = (
    "CELL_AMOUNT_NAMES",
    "CORE_STATE_NAMES",
    "LUMEN_AMOUNT_NAMES",
    "CompartmentAmounts",
    "CompartmentConcentrations",
    "StateLayout",
    "WholeCellState",
    "amount_charge_equivalents_fmol",
    "compartment_from_concentrations",
    "vector_from_mapping",
)
