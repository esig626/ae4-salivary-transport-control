"""Typed parameters and provenance registry for the Task 13B full model.

Defaults in this module are executable starting values for invariant and
integration tests.  They are not a calibrated parameter set.  Every field is
tagged with one of the provenance categories required by the Task 13B prompt;
call :func:`parameter_records` before accepting or freezing any calibration.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from enum import Enum
import math
from typing import Any, Iterator, Mapping


class Provenance(str, Enum):
    PRIMARY_MEASUREMENT = "PRIMARY_MEASUREMENT"
    PUBLISHED_MODEL = "PUBLISHED_MODEL"
    HISTORICAL_IMPLEMENTATION = "HISTORICAL_IMPLEMENTATION"
    DERIVED_CONSTRAINT = "DERIVED_CONSTRAINT"
    NEW_MODELING_DECISION = "NEW_MODELING_DECISION"


def _pfield(
    default: Any,
    *,
    provenance: Provenance,
    units: str,
    note: str,
) -> Any:
    return field(
        default=default,
        metadata={"provenance": provenance.value, "units": units, "note": note},
    )


@dataclass(frozen=True)
class PhysicalConstants:
    temperature_K: float = _pfield(
        310.15,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="K",
        note="37 C reference temperature; protocol-specific temperature must replace it",
    )
    gas_constant_J_mol_K: float = _pfield(
        8.31446261815324,
        provenance=Provenance.DERIVED_CONSTRAINT,
        units="J mol^-1 K^-1",
        note="CODATA exact/standard physical constant",
    )
    faraday_C_mol: float = _pfield(
        96485.33212,
        provenance=Provenance.DERIVED_CONSTRAINT,
        units="C mol^-1",
        note="Faraday conversion between charge current and molar flux",
    )

    @property
    def thermal_voltage_V(self) -> float:
        return self.gas_constant_J_mol_K * self.temperature_K / self.faraday_C_mol


@dataclass(frozen=True)
class AcidBaseParameters:
    carbon_pka1: float = _pfield(
        6.10,
        provenance=Provenance.PUBLISHED_MODEL,
        units="dimensionless (p-scale)",
        note="effective CO2/HCO3 pKa; temperature/context refinement required",
    )
    carbon_pka2: float = _pfield(
        10.30,
        provenance=Provenance.PUBLISHED_MODEL,
        units="dimensionless (p-scale)",
        note="effective HCO3/CO3 pKa; carbonate retained for charge closure",
    )
    water_pkw: float = _pfield(
        14.00,
        provenance=Provenance.PUBLISHED_MODEL,
        units="dimensionless (p-scale)",
        note="effective water ionization value; temperature refinement required",
    )
    cell_buffer_pka: float = _pfield(
        7.00,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="single effective non-carbonate buffer site; not identified",
    )
    lumen_buffer_pka: float = _pfield(
        7.00,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="irrelevant when luminal buffer amount is zero",
    )
    ph_lower: float = _pfield(
        3.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="numerical speciation bracket, not a physiological gate",
    )
    ph_upper: float = _pfield(
        11.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="numerical speciation bracket, not a physiological gate",
    )


@dataclass(frozen=True)
class BathParameters:
    na_mM: float = _pfield(
        145.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="representative open-bath value; protocol ledger must replace it",
    )
    k_mM: float = _pfield(
        5.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="representative open-bath value; protocol ledger must replace it",
    )
    cl_mM: float = _pfield(
        126.16159280366088,
        provenance=Provenance.DERIVED_CONSTRAINT,
        units="mM",
        note="chosen from bath cations and TA for exact reference electroneutrality",
    )
    tic_mM: float = _pfield(
        25.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="open-bath total inorganic carbon starting value",
    )
    ph: float = _pfield(
        7.40,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="representative bicarbonate-buffered bath value",
    )
    buffer_total_mM: float = _pfield(
        0.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM sites",
        note="bath non-carbonate buffer omitted from exchanged species",
    )
    untracked_osmolyte_mM: float = _pfield(
        0.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mOsm L^-1",
        note="explicit correction only if protocol osmolarity requires it",
    )


@dataclass(frozen=True)
class GeometryParameters:
    cell_buffer_total_fmol: float = _pfield(
        30.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol sites",
        note="effective finite buffer pool; must be WT-calibrated independently",
    )
    lumen_buffer_total_fmol: float = _pfield(
        1.0e-9,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol sites",
        note="positive numerical limit approximating no luminal fixed buffer",
    )
    fixed_cell_anion_equivalents_fmol: float = _pfield(
        78.06176220594485,
        provenance=Provenance.DERIVED_CONSTRAINT,
        units="fmol charge equivalents",
        note="derived from reference WT cations, chloride, and total alkalinity",
    )
    cell_impermeant_osmoles_fmol: float = _pfield(
        56.16159280366088,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol osmoles",
        note=(
            "other impermeant osmotic particles; explicitly excludes the finite "
            "buffer pool, which is added separately"
        ),
    )


@dataclass(frozen=True)
class HomeostasisParameters:
    """Reversible coarse laws for NKCC1, NHE1, AE2, and CO2 exchange."""

    nkcc1_capacity_fmol_s: float = _pfield(
        0.08,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol s^-1",
        note="placeholder capacity; no AE4-null secretion calibration",
    )
    nhe1_capacity_fmol_s: float = _pfield(
        0.02,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol s^-1",
        note="placeholder capacity; calibrated only against WT acid-base evidence",
    )
    ae2_capacity_fmol_s: float = _pfield(
        0.005,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol s^-1",
        note="small placeholder capacity; AE2-null is a validation gate",
    )
    thermodynamic_saturation_log_width: float = _pfield(
        2.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless log-affinity",
        note="sets tanh approach to capacity without moving reversal",
    )
    co2_basolateral_permeability_fmol_s_mM: float = _pfield(
        0.02,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol s^-1 mM^-1",
        note="linear neutral CO2 exchange coefficient",
    )
    co2_apical_permeability_fmol_s_mM: float = _pfield(
        0.01,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol s^-1 mM^-1",
        note="linear neutral CO2 exchange coefficient",
    )


@dataclass(frozen=True)
class MembraneParameters:
    nak_capacity_fmol_s: float = _pfield(
        0.08,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fmol cycles s^-1",
        note="total pump capacity shared between membranes",
    )
    nak_na_half_mM: float = _pfield(
        10.0,
        provenance=Provenance.PUBLISHED_MODEL,
        units="mM",
        note="effective intracellular Na half-saturation starting value",
    )
    nak_k_half_mM: float = _pfield(
        1.5,
        provenance=Provenance.PUBLISHED_MODEL,
        units="mM",
        note="effective extracellular K half-saturation starting value",
    )
    apical_pump_fraction: float = _pfield(
        0.10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fraction",
        note="apical localization is supported; quantitative fraction is not measured",
    )
    apical_k_fraction: float = _pfield(
        0.10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="fraction",
        note="apical K current is supported; quantitative fraction is not measured",
    )
    g_k_total_S: float = _pfield(
        4.0e-10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="S",
        note="total Ca-activated K conductance shared between membranes",
    )
    g_cl_apical_S: float = _pfield(
        6.0e-10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="S",
        note="apical Ca-activated Cl conductance",
    )
    calcium_half_uM: float = _pfield(
        0.30,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="uM",
        note="common test gating scale, not a frozen WT estimate",
    )
    calcium_hill: float = _pfield(
        2.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless",
        note="effective Ca gating exponent",
    )
    g_para_na_S: float = _pfield(
        3.0e-10,
        provenance=Provenance.PUBLISHED_MODEL,
        units="S",
        note="historical-lineage paracellular Na pathway; value not frozen",
    )
    g_para_k_S: float = _pfield(
        5.0e-11,
        provenance=Provenance.PUBLISHED_MODEL,
        units="S",
        note="historical-lineage paracellular K pathway; value not frozen",
    )
    g_para_cl_S: float = _pfield(
        2.5e-10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="S",
        note="charge-complete paracellular Cl pathway",
    )
    g_para_hco3_S: float = _pfield(
        5.0e-11,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="S",
        note="charge-complete paracellular bicarbonate pathway",
    )
    g_apical_background_S: float = _pfield(
        0.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="S",
        note="zero preserves the exact no-apical-K nesting limit",
    )
    g_basolateral_background_S: float = _pfield(
        5.0e-12,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="S",
        note="small K-selective regularizer; set zero only if closure remains nonsingular",
    )


@dataclass(frozen=True)
class WaterParameters:
    apical_hydraulic_pL_s_mOsm: float = _pfield(
        4.32e-3,
        provenance=Provenance.PUBLISHED_MODEL,
        units="pL s^-1 (mOsm L^-1)^-1",
        note=(
            "Palk 2010/2018 Pa=4.32e-12 L^2 mol^-1 s^-1 converted by "
            "1 mM=1e-3 mol/L and 1 L=1e12 pL"
        ),
    )
    basolateral_hydraulic_pL_s_mOsm: float = _pfield(
        5.15e-2,
        provenance=Provenance.PUBLISHED_MODEL,
        units="pL s^-1 (mOsm L^-1)^-1",
        note=(
            "Palk 2010/2018 Pb=5.15e-11 L^2 mol^-1 s^-1 converted by "
            "1 mM=1e-3 mol/L and 1 L=1e12 pL"
        ),
    )
    paracellular_hydraulic_pL_s_mOsm: float = _pfield(
        2.60e-4,
        provenance=Provenance.PUBLISHED_MODEL,
        units="pL s^-1 (mOsm L^-1)^-1",
        note=(
            "Palk 2010/2018 Pt=2.6e-13 L^2 mol^-1 s^-1 converted by "
            "1 mM=1e-3 mol/L and 1 L=1e12 pL"
        ),
    )
    lumen_dead_volume_pL: float = _pfield(
        0.10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="pL",
        note="non-draining local lumen volume",
    )
    outflow_rate_s: float = _pfield(
        0.20,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="s^-1",
        note="local compliance/outflow closure, not secretion calibration",
    )


@dataclass(frozen=True)
class InitialConditions:
    cell_na_mM: float = _pfield(
        15.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="positive reference value; replace by source-ledger WT range",
    )
    cell_k_mM: float = _pfield(
        140.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="positive reference value; replace by source-ledger WT range",
    )
    cell_cl_mM: float = _pfield(
        40.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="positive reference value; not a calibration claim",
    )
    cell_tic_mM: float = _pfield(
        20.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="reference TIC pending WT carbon measurement",
    )
    cell_ph: float = _pfield(
        7.20,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="reference pH pending evidence-ledger gate",
    )
    cell_volume_pL: float = _pfield(
        1.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="pL",
        note="single-cell reference volume pending geometry audit",
    )
    lumen_na_mM: float = _pfield(
        145.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="open-bath-like positive reference value",
    )
    lumen_k_mM: float = _pfield(
        5.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="open-bath-like positive reference value",
    )
    lumen_cl_mM: float = _pfield(
        126.16159279650835,
        provenance=Provenance.DERIVED_CONSTRAINT,
        units="mM",
        note="derived from lumen cations and TA for exact reference electroneutrality",
    )
    lumen_tic_mM: float = _pfield(
        25.0,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="mM",
        note="open-bath-like positive reference value",
    )
    lumen_ph: float = _pfield(
        7.40,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="dimensionless (p-scale)",
        note="open-bath-like positive reference value",
    )
    lumen_volume_pL: float = _pfield(
        0.10,
        provenance=Provenance.NEW_MODELING_DECISION,
        units="pL",
        note="finite local lumen reference volume pending geometry audit",
    )


@dataclass(frozen=True)
class FullModelParameters:
    constants: PhysicalConstants = field(default_factory=PhysicalConstants)
    acid_base: AcidBaseParameters = field(default_factory=AcidBaseParameters)
    bath: BathParameters = field(default_factory=BathParameters)
    geometry: GeometryParameters = field(default_factory=GeometryParameters)
    homeostasis: HomeostasisParameters = field(default_factory=HomeostasisParameters)
    membranes: MembraneParameters = field(default_factory=MembraneParameters)
    water: WaterParameters = field(default_factory=WaterParameters)
    initial: InitialConditions = field(default_factory=InitialConditions)

    def __post_init__(self) -> None:
        _validate_dataclass_numbers(self)

        strictly_positive = {
            "constants.temperature_K": self.constants.temperature_K,
            "constants.gas_constant_J_mol_K": self.constants.gas_constant_J_mol_K,
            "constants.faraday_C_mol": self.constants.faraday_C_mol,
            "acid_base.ph_lower": self.acid_base.ph_lower,
            "acid_base.ph_upper": self.acid_base.ph_upper,
            "bath.na_mM": self.bath.na_mM,
            "bath.k_mM": self.bath.k_mM,
            "bath.cl_mM": self.bath.cl_mM,
            "bath.tic_mM": self.bath.tic_mM,
            "homeostasis.thermodynamic_saturation_log_width": (
                self.homeostasis.thermodynamic_saturation_log_width
            ),
            "membranes.nak_na_half_mM": self.membranes.nak_na_half_mM,
            "membranes.nak_k_half_mM": self.membranes.nak_k_half_mM,
            "membranes.calcium_half_uM": self.membranes.calcium_half_uM,
            "membranes.calcium_hill": self.membranes.calcium_hill,
            "water.lumen_dead_volume_pL": self.water.lumen_dead_volume_pL,
            "initial.cell_na_mM": self.initial.cell_na_mM,
            "initial.cell_k_mM": self.initial.cell_k_mM,
            "initial.cell_cl_mM": self.initial.cell_cl_mM,
            "initial.cell_tic_mM": self.initial.cell_tic_mM,
            "initial.cell_volume_pL": self.initial.cell_volume_pL,
            "initial.lumen_na_mM": self.initial.lumen_na_mM,
            "initial.lumen_k_mM": self.initial.lumen_k_mM,
            "initial.lumen_cl_mM": self.initial.lumen_cl_mM,
            "initial.lumen_tic_mM": self.initial.lumen_tic_mM,
            "initial.lumen_volume_pL": self.initial.lumen_volume_pL,
        }
        for name, value in strictly_positive.items():
            if value <= 0.0:
                raise ValueError(f"{name} must be strictly positive")

        nonnegative = {
            "bath.buffer_total_mM": self.bath.buffer_total_mM,
            "bath.untracked_osmolyte_mM": self.bath.untracked_osmolyte_mM,
            "geometry.cell_buffer_total_fmol": self.geometry.cell_buffer_total_fmol,
            "geometry.lumen_buffer_total_fmol": self.geometry.lumen_buffer_total_fmol,
            "geometry.fixed_cell_anion_equivalents_fmol": (
                self.geometry.fixed_cell_anion_equivalents_fmol
            ),
            "geometry.cell_impermeant_osmoles_fmol": (
                self.geometry.cell_impermeant_osmoles_fmol
            ),
            "homeostasis.nkcc1_capacity_fmol_s": (
                self.homeostasis.nkcc1_capacity_fmol_s
            ),
            "homeostasis.nhe1_capacity_fmol_s": self.homeostasis.nhe1_capacity_fmol_s,
            "homeostasis.ae2_capacity_fmol_s": self.homeostasis.ae2_capacity_fmol_s,
            "homeostasis.co2_basolateral_permeability_fmol_s_mM": (
                self.homeostasis.co2_basolateral_permeability_fmol_s_mM
            ),
            "homeostasis.co2_apical_permeability_fmol_s_mM": (
                self.homeostasis.co2_apical_permeability_fmol_s_mM
            ),
            "membranes.nak_capacity_fmol_s": self.membranes.nak_capacity_fmol_s,
            "membranes.g_k_total_S": self.membranes.g_k_total_S,
            "membranes.g_cl_apical_S": self.membranes.g_cl_apical_S,
            "membranes.g_para_na_S": self.membranes.g_para_na_S,
            "membranes.g_para_k_S": self.membranes.g_para_k_S,
            "membranes.g_para_cl_S": self.membranes.g_para_cl_S,
            "membranes.g_para_hco3_S": self.membranes.g_para_hco3_S,
            "membranes.g_apical_background_S": self.membranes.g_apical_background_S,
            "membranes.g_basolateral_background_S": (
                self.membranes.g_basolateral_background_S
            ),
            "water.apical_hydraulic_pL_s_mOsm": (
                self.water.apical_hydraulic_pL_s_mOsm
            ),
            "water.basolateral_hydraulic_pL_s_mOsm": (
                self.water.basolateral_hydraulic_pL_s_mOsm
            ),
            "water.paracellular_hydraulic_pL_s_mOsm": (
                self.water.paracellular_hydraulic_pL_s_mOsm
            ),
            "water.outflow_rate_s": self.water.outflow_rate_s,
        }
        for name, value in nonnegative.items():
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")

        for name in ("apical_pump_fraction", "apical_k_fraction"):
            value = getattr(self.membranes, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")
        if self.acid_base.ph_lower >= self.acid_base.ph_upper:
            raise ValueError("acid-base pH bracket must be ordered")
        for name, value in (
            ("bath.ph", self.bath.ph),
            ("initial.cell_ph", self.initial.cell_ph),
            ("initial.lumen_ph", self.initial.lumen_ph),
        ):
            if not self.acid_base.ph_lower < value < self.acid_base.ph_upper:
                raise ValueError(f"{name} must lie strictly inside the pH bracket")


@dataclass(frozen=True)
class ParameterRecord:
    name: str
    value: Any
    units: str
    provenance: str
    note: str


def _validate_dataclass_numbers(instance: Any) -> None:
    if not is_dataclass(instance):
        return
    for item in fields(instance):
        value = getattr(instance, item.name)
        if is_dataclass(value):
            _validate_dataclass_numbers(value)
        elif isinstance(value, (float, int)) and not math.isfinite(float(value)):
            raise ValueError(f"{item.name} must be finite")


def parameter_records(parameters: FullModelParameters) -> tuple[ParameterRecord, ...]:
    """Flatten all tagged fields into a machine-readable provenance ledger."""

    records: list[ParameterRecord] = []

    def visit(value: Any, prefix: str) -> None:
        for item in fields(value):
            child = getattr(value, item.name)
            name = f"{prefix}.{item.name}" if prefix else item.name
            if is_dataclass(child):
                visit(child, name)
                continue
            metadata: Mapping[str, Any] = item.metadata
            records.append(
                ParameterRecord(
                    name=name,
                    value=child,
                    units=str(metadata.get("units", "UNSPECIFIED")),
                    provenance=str(
                        metadata.get("provenance", Provenance.NEW_MODELING_DECISION.value)
                    ),
                    note=str(metadata.get("note", "unannotated field")),
                )
            )

    visit(parameters, "")
    return tuple(records)


__all__ = (
    "AcidBaseParameters",
    "BathParameters",
    "FullModelParameters",
    "GeometryParameters",
    "HomeostasisParameters",
    "InitialConditions",
    "MembraneParameters",
    "ParameterRecord",
    "PhysicalConstants",
    "Provenance",
    "WaterParameters",
    "parameter_records",
)
