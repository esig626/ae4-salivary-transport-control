"""Transporter-only constraint and WT-only capacity utilities.

No function in this module accepts an AE4-null phenotype or secretion target.
The effective dose-response summaries are reconstructed only from CAL-T
values frozen in ``evidence_freeze.md``; missing intracellular assay
concentrations remain explicit assumptions.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from .adapter import transporter_reservoirs_with_carbonate
from .cycles import (
    CANDIDATE_DEFINITIONS,
    CycleDefinition,
    CycleParameters,
    ReservoirActivities,
    evaluate_cycle,
    exact_charge_per_branch,
    parameters_for_mutant,
)
from .evidence import assert_calibration_targets_allowed


DOSE_CONCENTRATIONS_MM = np.array((5.0, 25.0, 50.0, 100.0, 125.0, 150.0))
ASSAY_EC50_MM = {"na": 49.0, "k": 62.0}
ASSAY_HILL = {"na": 2.0, "k": 1.8}
ASSAY_RMIN = {"na": 0.3, "k": 0.4}
ASSAY_RMAX = {"na": 1.5, "k": 1.6}


@dataclass(frozen=True)
class DoseResponseFit:
    species: str
    apparent_ec50_mM: float
    apparent_hill: float
    fitted_minimum: float
    fitted_maximum: float
    normalized_rmse: float
    solver_success: bool


@dataclass(frozen=True)
class TransporterFitResult:
    model_id: str
    gauge_id: str
    na_attempt_scale: float
    k_attempt_scale: float
    common_barrier_multiplier: float
    na_loaded_energy_lump_mM: float
    k_loaded_energy_lump_mM: float
    pka_barrier_fold: float
    na_ec50_pred_mM: float
    na_hill_pred: float
    k_ec50_pred_mM: float
    k_hill_pred: float
    pure_k_to_na_amplitude_ratio: float
    assay_background_na: float
    assay_background_k: float
    assay_common_scale: float
    mapped_na_rmax: float
    mapped_k_rmax: float
    double_mutant_na_to_wt: float
    double_mutant_k_to_wt: float
    pka_effective_fold: float
    na_transported: bool
    k_transported: bool
    electroneutral: bool
    ldb_pass: bool
    mutant_hierarchy_pass: bool
    dose_shape_status: str
    regulation_status: str
    transporter_gate_pass: bool
    n_lumped_parameters: int
    complexity_penalty: int
    n_fitted_numeric_parameters: int
    n_fixed_assay_nuisance_values: int
    n_discrete_gauge_choices: int
    model_comparison_interpretation: str
    calibration_target_ids: str
    heldout_target_read: bool
    context_assumption: str
    failure_reason: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def nominal_assay_context(
    species: str,
    external_cation_mM: float,
    *,
    carbonate: bool = False,
) -> ReservoirActivities:
    """Declared low-external-Cl heterologous context for shape diagnostics.

    External bath values follow E16-07.  Intracellular concentrations were not
    measured in the dose-response experiment; the 25/120/50/10 mM tuple is a
    labeled salivary-reference sensitivity and cannot identify microscopic
    rates.  We therefore report apparent-shape dependence and repeat barrier
    gauges rather than calling these values estimates.
    """

    if species not in {"na", "k"}:
        raise ValueError("species must be na or k")
    na_e = external_cation_mM if species == "na" else 1e-6
    k_e = external_cation_mM if species == "k" else 1e-6
    if carbonate:
        return transporter_reservoirs_with_carbonate(
            na_i=25.0,
            k_i=120.0,
            cl_i=50.0,
            hco3_i=10.0,
            ph_i=6.91,
            na_e=na_e,
            k_e=k_e,
            cl_e=4.0,
            hco3_e=25.0,
            ph_e=7.4,
        )
    return ReservoirActivities(
        na_i=25.0,
        k_i=120.0,
        cl_i=50.0,
        hco3_i=10.0,
        na_e=na_e,
        k_e=k_e,
        cl_e=4.0,
        hco3_e=25.0,
    )


def _branch_for_species(definition: CycleDefinition, species: str) -> str | None:
    if species in definition.branch_current_edges:
        return species
    if len(definition.branch_current_edges) == 1:
        branch = next(iter(definition.branch_current_edges))
        source = definition.intracellular_source_vectors[branch]
        if abs(source.get(f"{species}_i", 0.0)) > 0.0:
            return branch
    return None


def assay_rates(
    definition: CycleDefinition,
    parameters: CycleParameters,
    species: str,
    concentrations_mM: Sequence[float] = DOSE_CONCENTRATIONS_MM,
    *,
    pka_activation: float = 0.0,
) -> np.ndarray:
    if _branch_for_species(definition, species) is None:
        return np.full(len(concentrations_mM), math.nan)
    carbonate = "co3_i" in definition.required_reservoir_species
    rates = []
    for concentration in concentrations_mM:
        evaluation = evaluate_cycle(
            definition,
            nominal_assay_context(species, float(concentration), carbonate=carbonate),
            parameters,
            pka_activation=pka_activation,
        )
        # E16-07 reports alkalinization, not a selected microscopic branch.
        # In the reverse assay direction, positive inorganic-carbon alkalinity
        # source HCO3 + 2*CO3 is the minimal model-matched observable.  The
        # factor two is charge-equivalent chemistry, not a hidden conversion of
        # carbonate into the bicarbonate state.  A wrong-sign net source is
        # retained as zero rather than hidden by abs(J_branch).
        sources = evaluation.intracellular_sources
        alkalinity_source = float(sources["hco3_i"] + 2.0 * sources["co3_i"])
        rates.append(max(0.0, alkalinity_source))
    return np.asarray(rates)


def _hill_model(concentrations: np.ndarray, values: np.ndarray) -> DoseResponseFit:
    maximum = max(float(np.max(values)), 1e-15)

    def residual(x: np.ndarray) -> np.ndarray:
        baseline, amplitude, log_ec50, log_hill = x
        ec50 = math.exp(log_ec50)
        hill = math.exp(log_hill)
        predicted = baseline + amplitude / (
            1.0 + np.exp(np.clip(hill * np.log(ec50 / concentrations), -700.0, 700.0))
        )
        return (predicted - values) / maximum

    x0 = np.array(
        (
            max(0.0, float(np.min(values))),
            max(1e-12, float(np.max(values) - np.min(values))),
            math.log(50.0),
            math.log(1.5),
        )
    )
    fit = least_squares(
        residual,
        x0,
        bounds=(
            np.array((0.0, 0.0, math.log(0.1), math.log(0.1))),
            np.array((2.0 * maximum, 10.0 * maximum, math.log(1000.0), math.log(6.0))),
        ),
        max_nfev=5000,
    )
    baseline, amplitude, log_ec50, log_hill = fit.x
    return DoseResponseFit(
        species="",
        apparent_ec50_mM=float(math.exp(log_ec50)),
        apparent_hill=float(math.exp(log_hill)),
        fitted_minimum=float(baseline),
        fitted_maximum=float(baseline + amplitude),
        normalized_rmse=float(math.sqrt(np.mean(residual(fit.x) ** 2))),
        solver_success=bool(fit.success),
    )


def fit_apparent_dose_response(
    definition: CycleDefinition, parameters: CycleParameters, species: str
) -> DoseResponseFit:
    values = assay_rates(definition, parameters, species)
    if not np.all(np.isfinite(values)) or np.max(values) <= 0.0:
        return DoseResponseFit(species, math.nan, math.nan, math.nan, math.nan, math.nan, False)
    return replace(_hill_model(DOSE_CONCENTRATIONS_MM, values), species=species)


def _gauge_overrides(definition: CycleDefinition, multiplier: float) -> Mapping[str, float]:
    return {
        edge.edge_id: multiplier
        for edge in definition.edges
        if edge.edge_id.endswith("common_cl_flip")
    }


def constrain_k_attempt_scale(
    definition: CycleDefinition,
    *,
    common_barrier_multiplier: float,
    base_parameters: CycleParameters | None = None,
) -> CycleParameters:
    """Match the reported pure-cation Rmax ratio within a declared assay gauge."""

    initial = base_parameters or CycleParameters()
    overrides = {
        **dict(initial.edge_barrier_overrides),
        **_gauge_overrides(definition, common_barrier_multiplier),
    }
    initial = replace(initial, edge_barrier_overrides=overrides)
    if "na" not in definition.branch_current_edges or "k" not in definition.branch_current_edges:
        return initial
    target_ratio = (
        (ASSAY_RMAX["k"] - ASSAY_RMIN["k"])
        / (ASSAY_RMAX["na"] - ASSAY_RMIN["na"])
    )

    def objective(log_scale: float) -> float:
        p = replace(
            initial,
            k_attempt_scale=math.exp(log_scale),
        )
        na_fit = fit_apparent_dose_response(definition, p, "na")
        k_fit = fit_apparent_dose_response(definition, p, "k")
        na_amplitude = na_fit.fitted_maximum - na_fit.fitted_minimum
        k_amplitude = k_fit.fitted_maximum - k_fit.fitted_minimum
        if na_amplitude <= 0.0 or k_amplitude <= 0.0:
            return 1e6
        return (math.log(k_amplitude / na_amplitude) - math.log(target_ratio)) ** 2

    fit = minimize_scalar(objective, bounds=(-8.0, 8.0), method="bounded")
    return replace(
        initial,
        k_attempt_scale=float(math.exp(fit.x)),
    )


def constrain_loaded_state_energies(
    definition: CycleDefinition, parameters: CycleParameters
) -> CycleParameters:
    """Fit only two effective standard-energy lumps to E16-07 EC50 values.

    The resulting ``ec50_*`` fields are *not* reported as microscopic Kd.
    Together with the fixed 25 mM bicarbonate standard state they are the
    loaded-state energy combinations needed to reproduce the apparent bath
    dependence in the declared assay context.
    """

    if _branch_for_species(definition, "na") is None or _branch_for_species(definition, "k") is None:
        return parameters

    def residual(log_values: np.ndarray) -> np.ndarray:
        candidate = replace(
            parameters,
            ec50_na_mM=math.exp(float(log_values[0])),
            ec50_k_mM=math.exp(float(log_values[1])),
        )
        na_fit = fit_apparent_dose_response(definition, candidate, "na")
        k_fit = fit_apparent_dose_response(definition, candidate, "k")
        if not na_fit.solver_success or not k_fit.solver_success:
            return np.array((1e3, 1e3))
        return np.array(
            (
                math.log(max(na_fit.apparent_ec50_mM, 1e-12) / ASSAY_EC50_MM["na"]),
                math.log(max(k_fit.apparent_ec50_mM, 1e-12) / ASSAY_EC50_MM["k"]),
            )
        )

    fit = least_squares(
        residual,
        np.log(np.array((500.0, 500.0))),
        bounds=(np.log(np.array((0.1, 0.1))), np.log(np.array((1e6, 1e6)))),
        max_nfev=200,
        xtol=1e-10,
        ftol=1e-10,
        gtol=1e-10,
    )
    return replace(
        parameters,
        ec50_na_mM=float(math.exp(fit.x[0])),
        ec50_k_mM=float(math.exp(fit.x[1])),
    )


def _pka_effective_fold(
    definition: CycleDefinition, parameters: CycleParameters, species: str = "na"
) -> float:
    if _branch_for_species(definition, species) is None:
        return math.nan
    context = nominal_assay_context(
        species, 145.0, carbonate="co3_i" in definition.required_reservoir_species
    )
    basal_evaluation = evaluate_cycle(
        definition, context, parameters, pka_activation=0.0
    )
    stimulated_evaluation = evaluate_cycle(
        definition, context, parameters, pka_activation=1.0
    )
    basal = max(
        0.0,
        basal_evaluation.intracellular_sources["hco3_i"]
        + 2.0 * basal_evaluation.intracellular_sources["co3_i"],
    )
    stimulated = max(
        0.0,
        stimulated_evaluation.intracellular_sources["hco3_i"]
        + 2.0 * stimulated_evaluation.intracellular_sources["co3_i"],
    )
    return stimulated / basal if basal > 0.0 else math.nan


def constrain_pka_fold(
    definition: CycleDefinition, parameters: CycleParameters
) -> CycleParameters:
    if not definition.pka_switch_edge_ids:
        return parameters
    target = 1.25

    def objective(log_gamma: float) -> float:
        candidate = replace(parameters, pka_barrier_fold=math.exp(log_gamma))
        fold = _pka_effective_fold(definition, candidate)
        return (fold - target) ** 2 if math.isfinite(fold) else 1e6

    fit = minimize_scalar(objective, bounds=(0.0, math.log(1e6)), method="bounded")
    return replace(parameters, pka_barrier_fold=float(math.exp(fit.x)))


def evaluate_transporter_candidate(
    model_id: str,
    *,
    gauge_id: str,
    common_barrier_multiplier: float,
) -> TransporterFitResult:
    assert_calibration_targets_allowed(
        (
            "E16-03",
            "E16-04",
            "E16-05",
            "E16-07",
            "E21-03",
            "E25-04",
        )
    )
    definition = CANDIDATE_DEFINITIONS[model_id]
    parameters = CycleParameters(
        edge_barrier_overrides=_gauge_overrides(definition, common_barrier_multiplier)
    )
    parameters = constrain_loaded_state_energies(definition, parameters)
    parameters = constrain_k_attempt_scale(
        definition,
        common_barrier_multiplier=common_barrier_multiplier,
        base_parameters=parameters,
    )
    parameters = constrain_pka_fold(definition, parameters)
    na_fit = fit_apparent_dose_response(definition, parameters, "na")
    k_fit = fit_apparent_dose_response(definition, parameters, "k")
    na_branch = _branch_for_species(definition, "na")
    k_branch = _branch_for_species(definition, "k")
    na_transported = na_branch is not None
    k_transported = k_branch is not None
    charge = exact_charge_per_branch(definition)
    electroneutral = all(abs(value) <= 1e-12 for value in charge.values())
    context_na = nominal_assay_context(
        "na", 145.0, carbonate="co3_i" in definition.required_reservoir_species
    )
    context_k = nominal_assay_context(
        "k", 145.0, carbonate="co3_i" in definition.required_reservoir_species
    )
    ldb_evaluations = [
        evaluate_cycle(definition, context_na, parameters),
        evaluate_cycle(definition, context_k, parameters),
    ]
    ldb_pass = all(
        result.local_detailed_balance_residual <= 1e-10
        and result.entropy_production >= 0.0
        for result in ldb_evaluations
    )
    if na_transported and k_transported:
        na_wt_evaluation = evaluate_cycle(definition, context_na, parameters)
        k_wt_evaluation = evaluate_cycle(definition, context_k, parameters)
        na_wt = max(
            0.0,
            na_wt_evaluation.intracellular_sources["hco3_i"]
            + 2.0 * na_wt_evaluation.intracellular_sources["co3_i"],
        )
        k_wt = max(
            0.0,
            k_wt_evaluation.intracellular_sources["hco3_i"]
            + 2.0 * k_wt_evaluation.intracellular_sources["co3_i"],
        )
        mutant_parameters = parameters_for_mutant(parameters, "T448I_T756A")
        na_mut_evaluation = evaluate_cycle(definition, context_na, mutant_parameters)
        k_mut_evaluation = evaluate_cycle(definition, context_k, mutant_parameters)
        na_mut = max(
            0.0,
            na_mut_evaluation.intracellular_sources["hco3_i"]
            + 2.0 * na_mut_evaluation.intracellular_sources["co3_i"],
        )
        k_mut = max(
            0.0,
            k_mut_evaluation.intracellular_sources["hco3_i"]
            + 2.0 * k_mut_evaluation.intracellular_sources["co3_i"],
        )
        na_mut_ratio = na_mut / na_wt if na_wt else math.nan
        k_mut_ratio = k_mut / k_wt if k_wt else math.nan
        mutant_pass = na_mut_ratio < 0.05 and k_mut_ratio > 0.10
        na_amplitude = na_fit.fitted_maximum - na_fit.fitted_minimum
        k_amplitude = k_fit.fitted_maximum - k_fit.fitted_minimum
        pure_ratio = k_amplitude / na_amplitude if na_amplitude > 0.0 else math.nan
    else:
        na_mut_ratio = k_mut_ratio = pure_ratio = math.nan
        mutant_pass = False
    amplitudes = np.array(
        (
            na_fit.fitted_maximum - na_fit.fitted_minimum,
            k_fit.fitted_maximum - k_fit.fitted_minimum,
        )
    )
    observed_amplitudes = np.array(
        (
            ASSAY_RMAX["na"] - ASSAY_RMIN["na"],
            ASSAY_RMAX["k"] - ASSAY_RMIN["k"],
        )
    )
    if np.all(amplitudes > 0.0) and np.all(np.isfinite(amplitudes)):
        assay_scale = float(np.dot(amplitudes, observed_amplitudes) / np.dot(amplitudes, amplitudes))
        mapped_na_rmax = ASSAY_RMIN["na"] + assay_scale * amplitudes[0]
        mapped_k_rmax = ASSAY_RMIN["k"] + assay_scale * amplitudes[1]
        amplitude_pass = (
            abs(mapped_na_rmax - ASSAY_RMAX["na"]) <= 0.10
            and abs(mapped_k_rmax - ASSAY_RMAX["k"]) <= 0.10
        )
    else:
        assay_scale = mapped_na_rmax = mapped_k_rmax = math.nan
        amplitude_pass = False
    dose_pass = bool(
        na_fit.solver_success
        and k_fit.solver_success
        and 0.7 * ASSAY_EC50_MM["na"] <= na_fit.apparent_ec50_mM <= 1.3 * ASSAY_EC50_MM["na"]
        and 0.7 * ASSAY_EC50_MM["k"] <= k_fit.apparent_ec50_mM <= 1.3 * ASSAY_EC50_MM["k"]
        and abs(na_fit.apparent_hill - ASSAY_HILL["na"]) <= 0.4
        and abs(k_fit.apparent_hill - ASSAY_HILL["k"]) <= 0.4
        and amplitude_pass
    )
    dose_status = (
        "compatible with reported effective EC50/Hill under declared assay context"
        if dose_pass
        else "not compatible with both reported effective Hill shapes under declared context"
    )
    pka_fold = _pka_effective_fold(definition, parameters) if definition.pka_switch_edge_ids else math.nan
    if definition.pka_switch_edge_ids:
        regulation_pass = math.isfinite(pka_fold) and abs(pka_fold - 1.25) <= 0.05
        regulation_status = (
            "assay fold constrained"
            if regulation_pass
            else "edge placement cannot attain approximately 1.25 fold within gamma<=1e6"
        )
    else:
        regulation_pass = True
        regulation_status = "unregulated core; requires nested SR6 before stimulated claim"
    # SR3/SR4b directly conflict with the pure-cation/direction evidence even
    # though their cycle is thermodynamically valid.
    direction_conflict = model_id in {"SR3_SEQUENTIAL_1111", "SR4B_CARBONATE_1112"}
    gate_pass = bool(
        na_transported
        and k_transported
        and electroneutral
        and ldb_pass
        and mutant_pass
        and dose_pass
        and regulation_pass
        and not direction_conflict
    )
    failure = []
    if not na_transported or not k_transported:
        failure.append("does not carry both directly observed cations")
    if not dose_pass:
        failure.append("effective dose shape mismatch")
    if not mutant_pass:
        failure.append("cannot express Na-collapse/K-retention mutant limit")
    if not regulation_pass:
        failure.append("PKA edge placement misses the recombinant fold")
    if direction_conflict:
        failure.append("stand-alone sequential cycle conflicts with direct Na direction/K-only operation")
    return TransporterFitResult(
        model_id=model_id,
        gauge_id=gauge_id,
        na_attempt_scale=parameters.na_attempt_scale,
        k_attempt_scale=parameters.k_attempt_scale,
        common_barrier_multiplier=common_barrier_multiplier,
        na_loaded_energy_lump_mM=parameters.ec50_na_mM,
        k_loaded_energy_lump_mM=parameters.ec50_k_mM,
        pka_barrier_fold=parameters.pka_barrier_fold,
        na_ec50_pred_mM=na_fit.apparent_ec50_mM,
        na_hill_pred=na_fit.apparent_hill,
        k_ec50_pred_mM=k_fit.apparent_ec50_mM,
        k_hill_pred=k_fit.apparent_hill,
        pure_k_to_na_amplitude_ratio=pure_ratio,
        assay_background_na=ASSAY_RMIN["na"],
        assay_background_k=ASSAY_RMIN["k"],
        assay_common_scale=assay_scale,
        mapped_na_rmax=mapped_na_rmax,
        mapped_k_rmax=mapped_k_rmax,
        double_mutant_na_to_wt=na_mut_ratio,
        double_mutant_k_to_wt=k_mut_ratio,
        pka_effective_fold=pka_fold,
        na_transported=na_transported,
        k_transported=k_transported,
        electroneutral=electroneutral,
        ldb_pass=ldb_pass,
        mutant_hierarchy_pass=mutant_pass,
        dose_shape_status=dose_status,
        regulation_status=regulation_status,
        transporter_gate_pass=gate_pass,
        n_lumped_parameters=(
            2  # Na/K loaded-state energy combinations
            + 1  # Na/K attempt-rate ratio
            + 1  # common assay output scale
            + (1 if definition.pka_switch_edge_ids else 0)
        ),
        complexity_penalty=(
            2 + 1 + 1 + (1 if definition.pka_switch_edge_ids else 0)
        ),
        n_fitted_numeric_parameters=(
            2 + 1 + 1 + (1 if definition.pka_switch_edge_ids else 0)
        ),
        n_fixed_assay_nuisance_values=2,
        n_discrete_gauge_choices=1,
        model_comparison_interpretation=(
            "conditional expressivity/compatibility comparison only; not AIC because raw "
            "curve errors/covariances and intracellular assay concentrations are unavailable"
        ),
        calibration_target_ids="E16-03;E16-04;E16-05;E16-07;E21-03;E25-04",
        heldout_target_read=False,
        context_assumption=(
            "E16 external dose grid and baths; unmeasured intracellular assay values fixed "
            "to the declared 25/120/50/10 mM sensitivity; alkalinization mapped to "
            "HCO3 source plus twice CO3 source; three common-barrier gauges"
        ),
        failure_reason="; ".join(failure),
    )


def transporter_fit_grid() -> tuple[TransporterFitResult, ...]:
    rows = []
    for model_id in CANDIDATE_DEFINITIONS:
        for gauge_id, multiplier in (("slow_common", 0.1), ("unit_common", 1.0), ("fast_common", 10.0)):
            rows.append(
                evaluate_transporter_candidate(
                    model_id,
                    gauge_id=gauge_id,
                    common_barrier_multiplier=multiplier,
                )
            )
    return tuple(rows)
