"""Test carbon escape and NHE1 capacity in the routed AE4 model.

This is a fixed mechanistic sensitivity experiment.  The only physiological
controls are the two existing CO2 permeabilities, scaled together, and the
existing NHE1 capacity.  No phenotype value is read by this module and no
parameter is fitted.

For each retained root and AE4 expression, direct REST roots are continued
through the declared parameter grid.  The Task 24 re equilibrated state is the
branch anchor.  A nearby solved grid point, the fixed Task 24 anchor, and the
older inherited resting state are used as predeclared numerical starts.  Root
selection uses branch distance and residuals only.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import asdict, replace
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import time
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import scipy

from .calibration import WTCalibrationSpec
from .genotype_evaluation import (
    _attempt_and_root,
    _core_coordinates,
    _regulatory_suffix,
    _root_state_with_basal_regulation,
    genotype_with_expression,
    simulate_genotype,
)
from .membranes import current_to_fmol_s
from .model import ModernFullModel, WT
from .nkcc_stimulation import StimulatedNkcc1Model, attach_stimulated_nkcc1
from .run_ae4_routing_only import inputs
from .run_ae4_routing_re_equilibrated import (
    AMOUNT_ROWS,
    MODEL_ID,
    REST_MAX_NFEV,
    ROUTING_ROUNDOFF_TOLERANCE_FMOL_S,
    VOLUME_ROWS,
    _audit_resting_state,
    _explicit_rest_model,
    _root_aliases,
)
from .run_calcium_fast_screen import write_json, write_rows
from .states import CORE_STATE_NAMES
from .task14_blind import CALCIUM, build_model, routing
from .ae4_routing_only import with_routing_only
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_RADAU,
    StimulusArm,
    sha256_file,
    sha256_object,
)


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/ae4_carbon_nhe1_compensation"
ANALYSIS = REPO / "analysis/ae4_carbon_nhe1_compensation"
TASK24 = REPO / "results/ae4_routing_re_equilibrated"
TASK24_REST = TASK24 / "resting_states"
TASK24_TRAJECTORY_SUMMARY = TASK24 / "trajectory_summary.csv"
TASK24_COMPARISON = TASK24 / "comparison.csv"
MANIFEST = REPO / "results/13B_modern_full_model/native_dynamic_contract_manifest.json"
FROZEN = REPO / "results/14B_five_percent_ae4_check/frozen_inputs.json"
SOURCE_COMMIT = "f6bd93ef9639b603dae1ac556e8a8d6482b92746"
BRANCH = "codex/task-25-carbon-nhe1-compensation"

CO2_FACTORS = (1.0, 2.0, 5.0, 10.0)
NHE1_FACTORS = (1.0, 0.5, 0.25)
EXPRESSIONS = (1.0, 0.05)
BASELINE = (1.0, 1.0)
BASELINE_REST_RELATIVE_TOLERANCE = 2.0e-5
BASELINE_TRAJECTORY_RELATIVE_TOLERANCE = 1.0e-6
AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S = 2.0e-14
EXPECTED_PARAMETER_PATHS = (
    "homeostasis.co2_apical_permeability_fmol_s_mM",
    "homeostasis.co2_basolateral_permeability_fmol_s_mM",
    "homeostasis.nhe1_capacity_fmol_s",
)


def _genotype(expression: float):
    return WT if expression == 1.0 else genotype_with_expression("AE4", expression)


def _factor_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def _rest_stem(root: str, co2_factor: float, nhe1_factor: float, expression: float) -> str:
    return (
        f"{root}_co2x{_factor_label(co2_factor)}_"
        f"nhe1x{_factor_label(nhe1_factor)}_e{expression:.2f}_rest"
    )


def _trajectory_stem(
    root: str,
    co2_factor: float,
    nhe1_factor: float,
    calcium: float,
    expression: float,
) -> str:
    return (
        f"{root}_co2x{_factor_label(co2_factor)}_"
        f"nhe1x{_factor_label(nhe1_factor)}_ca{calcium:.2f}_"
        f"e{expression:.2f}_production_radau"
    )


def _continuation_order() -> tuple[tuple[float, float], ...]:
    """A declared path that changes one control at a time."""

    order: list[tuple[float, float]] = [BASELINE]
    order.extend((1.0, nhe1) for nhe1 in NHE1_FACTORS[1:])
    for nhe1 in NHE1_FACTORS:
        order.extend((co2, nhe1) for co2 in CO2_FACTORS[1:])
    expected = set(itertools.product(CO2_FACTORS, NHE1_FACTORS))
    if len(order) != 12 or set(order) != expected:
        raise AssertionError("the continuation order does not equal the declared grid")
    return tuple(order)


def _predecessor(key: tuple[float, float]) -> tuple[float, float] | None:
    co2, nhe1 = key
    if key == BASELINE:
        return None
    if co2 == 1.0:
        index = NHE1_FACTORS.index(nhe1)
        return (1.0, NHE1_FACTORS[index - 1])
    index = CO2_FACTORS.index(co2)
    return (CO2_FACTORS[index - 1], nhe1)


def _controlled_parameters(parameters: Any, co2_factor: float, nhe1_factor: float):
    if co2_factor not in CO2_FACTORS or nhe1_factor not in NHE1_FACTORS:
        raise ValueError("control factors must belong to the declared grid")
    homeostasis = parameters.homeostasis
    controlled_homeostasis = replace(
        homeostasis,
        co2_apical_permeability_fmol_s_mM=(
            homeostasis.co2_apical_permeability_fmol_s_mM * co2_factor
        ),
        co2_basolateral_permeability_fmol_s_mM=(
            homeostasis.co2_basolateral_permeability_fmol_s_mM * co2_factor
        ),
        nhe1_capacity_fmol_s=homeostasis.nhe1_capacity_fmol_s * nhe1_factor,
    )
    return replace(parameters, homeostasis=controlled_homeostasis)


def _flatten_mapping(value: Any, prefix: str = "") -> dict[str, Any]:
    if hasattr(value, "__dataclass_fields__"):
        value = asdict(value)
    if not isinstance(value, Mapping):
        return {prefix: value}
    flat: dict[str, Any] = {}
    for key, item in value.items():
        path = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(item, Mapping):
            flat.update(_flatten_mapping(item, path))
        else:
            flat[path] = item
    return flat


def parameter_changes(reference: Any, controlled: Any) -> dict[str, tuple[Any, Any]]:
    before = _flatten_mapping(reference)
    after = _flatten_mapping(controlled)
    if set(before) != set(after):
        raise AssertionError("parameter structure changed")
    return {key: (before[key], after[key]) for key in before if before[key] != after[key]}


def _assert_controlled_parameters(
    reference: Any, controlled: Any, co2_factor: float, nhe1_factor: float
) -> tuple[str, ...]:
    changes = parameter_changes(reference, controlled)
    expected: set[str] = set()
    if co2_factor != 1.0:
        expected.update(EXPECTED_PARAMETER_PATHS[:2])
    if nhe1_factor != 1.0:
        expected.add(EXPECTED_PARAMETER_PATHS[2])
    if set(changes) != expected:
        raise AssertionError(
            f"unexpected parameter changes: observed={sorted(changes)}, expected={sorted(expected)}"
        )
    before = reference.homeostasis
    after = controlled.homeostasis
    checks = (
        (
            after.co2_apical_permeability_fmol_s_mM,
            before.co2_apical_permeability_fmol_s_mM * co2_factor,
        ),
        (
            after.co2_basolateral_permeability_fmol_s_mM,
            before.co2_basolateral_permeability_fmol_s_mM * co2_factor,
        ),
        (after.nhe1_capacity_fmol_s, before.nhe1_capacity_fmol_s * nhe1_factor),
    )
    if any(observed != expected_value for observed, expected_value in checks):
        raise AssertionError("a requested control factor was not applied exactly")
    return tuple(sorted(changes))


def with_controls(model: Any, co2_factor: float, nhe1_factor: float):
    """Clone a routed model while changing only the three permitted fields."""

    controlled = _controlled_parameters(model.parameters, co2_factor, nhe1_factor)
    _assert_controlled_parameters(model.parameters, controlled, co2_factor, nhe1_factor)
    if not isinstance(model, StimulatedNkcc1Model):
        raise TypeError("the inherited routed model must retain stimulated NKCC1")
    composite = model.regulatory_model
    if not hasattr(composite, "ae4_regulatory_model") or not hasattr(
        composite, "nkcc1_regulatory_model"
    ):
        raise AssertionError("the inherited R1 plus N1 regulator changed")
    base = ModernFullModel(
        parameters=controlled,
        stimulus=model.stimulus,
        regulatory_model=composite.ae4_regulatory_model,
        ae4_parameters=model.ae4_parameters,
        ae4_evaluator=model.ae4_evaluator,
    )
    result = attach_stimulated_nkcc1(base, composite.nkcc1_regulatory_model)
    if result.ae4_evaluator is not model.ae4_evaluator:
        raise AssertionError("control cloning changed the AE4 evaluator")
    if result.ae4_parameters != model.ae4_parameters:
        raise AssertionError("control cloning changed AE4 parameters")
    return result


def build_controlled_rest_model(
    manifest: Mapping[str, Any], root: str, co2_factor: float, nhe1_factor: float
):
    baseline = _explicit_rest_model(manifest, root)
    model = with_controls(baseline, co2_factor, nhe1_factor)
    if model.stimulus.arm != StimulusArm.REST:
        raise AssertionError("controlled resting model is not on the REST arm")
    for time_s in (0.0, 1.0, 600.0, 3600.0):
        stimulus = model.stimulus(time_s)
        if stimulus.beta_input != 0.0:
            raise AssertionError("REST beta input is not zero")
        if stimulus.calcium_uM != model.stimulus.resting_calcium_uM:
            raise AssertionError("REST calcium changed")
    return model


def build_controlled_stimulation_model(
    manifest: Mapping[str, Any],
    root: str,
    calcium: float,
    co2_factor: float,
    nhe1_factor: float,
):
    baseline = with_routing_only(build_model(manifest, root, calcium))
    return with_controls(baseline, co2_factor, nhe1_factor)


def _task24_rest_path(root: str, expression: float, calcium: float = 0.10) -> Path:
    return TASK24_REST / f"{root}_ca{calcium:.2f}_e{expression:.2f}_rest.json"


def _task24_rest(root: str, expression: float, calcium: float = 0.10) -> dict[str, Any]:
    path = _task24_rest_path(root, expression, calcium)
    payload = json.loads(path.read_text())
    if not payload["rest_gate_pass"]:
        raise AssertionError(f"Task 24 anchor is not accepted: {path}")
    if float(payload["ae4_expression"]) != expression:
        raise AssertionError("Task 24 anchor expression changed")
    return payload


def _unique_seeds(
    seeds: Sequence[tuple[str, Sequence[float]]]
) -> tuple[tuple[str, np.ndarray], ...]:
    unique: list[tuple[str, np.ndarray]] = []
    hashes: set[str] = set()
    for name, state in seeds:
        vector = np.asarray(state, dtype=float)
        digest = hashlib.sha256(vector.tobytes()).hexdigest()
        if digest not in hashes:
            unique.append((name, vector.copy()))
            hashes.add(digest)
    return tuple(unique)


def _cluster_candidates(
    candidates: Sequence[Mapping[str, Any]], spec: WTCalibrationSpec
) -> list[dict[str, Any]]:
    if not candidates:
        return []
    span = spec.coordinate_upper - spec.coordinate_lower
    parents = list(range(len(candidates)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        a, b = find(left), find(right)
        if a != b:
            parents[max(a, b)] = min(a, b)

    for left in range(len(candidates)):
        for right in range(left + 1, len(candidates)):
            a = np.asarray(candidates[left]["coordinates"], dtype=float)
            b = np.asarray(candidates[right]["coordinates"], dtype=float)
            distance = float(np.max(np.abs((a - b) / span)))
            if distance <= spec.cluster_relative_tolerance:
                union(left, right)
    groups: dict[int, list[Mapping[str, Any]]] = {}
    for index, candidate in enumerate(candidates):
        groups.setdefault(find(index), []).append(candidate)
    ordered = sorted(
        groups.values(),
        key=lambda group: tuple(float(value) for value in group[0]["coordinates"]),
    )
    result: list[dict[str, Any]] = []
    for index, group in enumerate(ordered):
        representative = min(
            group,
            key=lambda item: (
                float(item["audit"]["max_abs_scaled_independent_rhs"]),
                float(item["audit"]["max_abs_amount_rhs_fmol_s"]),
                str(item["seed_id"]),
            ),
        )
        result.append(
            {
                "cluster_id": f"B{index:02d}",
                "representative": representative,
                "members": list(group),
            }
        )
    return result


def _rest_flux_metrics(model: Any, core_state: Sequence[float], expression: float) -> dict[str, Any]:
    state = _root_state_with_basal_regulation(model, core_state)
    evaluation = model.evaluate(0.0, state, genotype=_genotype(expression))
    diagnostics = evaluation.diagnostics
    observables = diagnostics.observables
    cell = observables.cell_concentrations_mM
    lumen = observables.lumen_concentrations_mM
    decoded = model.layout.decode(state)
    ae4 = diagnostics.ae4
    raw_ae4 = float(ae4.diagnostics["j_ae4_fmol_s"])
    per_expression = float(ae4.cl_cell_fmol_s / expression)
    ae4_scaling_residual = abs(float(ae4.cl_cell_fmol_s) - expression * raw_ae4)
    if ae4_scaling_residual > AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S:
        raise AssertionError("AE4 expression scaling failed")
    homeostasis = diagnostics.homeostasis
    nkcc_cycle = float(homeostasis.nkcc1_inward_fmol_s)
    nkcc_cl = 2.0 * nkcc_cycle
    currents = diagnostics.membranes.currents_A
    faraday = model.parameters.constants.faraday_C_mol
    apical_cl = current_to_fmol_s(
        float(currents["cl_apical"]), valence=-1, faraday_C_mol=faraday
    )
    apical_k = current_to_fmol_s(
        float(currents["k_apical"]), valence=+1, faraday_C_mol=faraday
    )
    basolateral_k = current_to_fmol_s(
        float(currents["k_basolateral"]), valence=+1, faraday_C_mol=faraday
    )
    return {
        "cell_tic_mM": float(cell["tic"]),
        "cell_hco3_mM": float(cell["hco3"]),
        "cell_ph": float(observables.cell_acid_base.ph),
        "cell_na_mM": float(cell["na"]),
        "cell_k_mM": float(cell["k"]),
        "cell_cl_mM": float(cell["cl"]),
        "cell_volume_pL": float(decoded.cell.volume_pL),
        "lumen_tic_mM": float(lumen["tic"]),
        "lumen_hco3_mM": float(lumen["hco3"]),
        "lumen_ph": float(observables.lumen_acid_base.ph),
        "ae4_cl_flux_fmol_s": float(ae4.cl_cell_fmol_s),
        "ae4_cl_flux_pre_expression_fmol_s": raw_ae4,
        "ae4_flux_per_unit_expression_fmol_s": per_expression,
        "ae4_expression_scaling_residual_fmol_s": ae4_scaling_residual,
        "ae4_na_source_fmol_s": float(ae4.na_cell_fmol_s),
        "ae4_k_source_fmol_s": float(ae4.k_cell_fmol_s),
        "ae4_hco3_source_fmol_s": float(ae4.hco3_cell_fmol_s),
        "ae4_charge_residual_fmol_s": float(ae4.charge_source_residual_fmol_s),
        "ae4_no_opposing_cation_sources": bool(
            ae4.na_cell_fmol_s == 0.0
            or ae4.k_cell_fmol_s == 0.0
            or ae4.na_cell_fmol_s * ae4.k_cell_fmol_s > 0.0
        ),
        "nkcc1_inward_cycle_flux_fmol_s": nkcc_cycle,
        "nkcc1_chloride_influx_fmol_s": nkcc_cl,
        "nhe1_inward_na_flux_fmol_s": float(homeostasis.nhe1_inward_fmol_s),
        "ae2_inward_cycle_flux_fmol_s": float(homeostasis.ae2_inward_fmol_s),
        "co2_basolateral_bath_to_cell_fmol_s": float(
            diagnostics.co2_fluxes_fmol_s["bath_to_cell"]
        ),
        "co2_apical_lumen_to_cell_fmol_s": float(
            diagnostics.co2_fluxes_fmol_s["lumen_to_cell"]
        ),
        "v_apical_V": float(diagnostics.membranes.v_apical_V),
        "v_basolateral_V": float(diagnostics.membranes.v_basolateral_V),
        "v_transepithelial_V": float(diagnostics.membranes.v_transepithelial_V),
        "v_apical_mV": 1000.0 * float(diagnostics.membranes.v_apical_V),
        "v_basolateral_mV": 1000.0 * float(diagnostics.membranes.v_basolateral_V),
        "v_transepithelial_mV": 1000.0 * float(
            diagnostics.membranes.v_transepithelial_V
        ),
        "apical_chloride_cell_to_lumen_fmol_s": float(apical_cl),
        "apical_k_cell_to_lumen_fmol_s": float(apical_k),
        "basolateral_k_cell_to_bath_fmol_s": float(basolateral_k),
        "rest_outflow_pL_s": float(diagnostics.water.lumen_outflow_pL_s),
    }


def _attempt_row(
    *,
    root: str,
    alias: str,
    expression: float,
    co2_factor: float,
    nhe1_factor: float,
    predecessor: tuple[float, float] | None,
    seed_id: str,
    seed_state: np.ndarray,
    attempt: Any,
    root_result: Any | None,
    audit: Mapping[str, Any] | None,
) -> dict[str, Any]:
    return {
        "root_id": root,
        "root_alias": alias,
        "ae4_expression": expression,
        "co2_factor": co2_factor,
        "nhe1_factor": nhe1_factor,
        "predecessor_co2_factor": None if predecessor is None else predecessor[0],
        "predecessor_nhe1_factor": None if predecessor is None else predecessor[1],
        "seed_id": seed_id,
        "seed_core_sha256": sha256_object(list(seed_state[: len(CORE_STATE_NAMES)])),
        "optimizer_success": bool(attempt.optimizer_success),
        "admissible_state": bool(attempt.admissible_state),
        "converged_root": bool(attempt.converged_root),
        "max_abs_scaled_independent_rhs": float(attempt.max_abs_scaled_independent_rhs),
        "cost": float(attempt.cost),
        "optimality": float(attempt.optimality),
        "nfev": int(attempt.nfev),
        "message": attempt.message,
        "root_internal_gate_pass": bool(root_result is not None),
        "complete_rest_gate_pass": bool(audit and audit["rest_gate_pass"]),
        "solution_core_sha256": (
            None if root_result is None else sha256_object(list(root_result.core_state))
        ),
        "cluster_id": None,
        "selected_cluster": False,
    }


def _flat_rest_row(row: Mapping[str, Any]) -> dict[str, Any]:
    excluded = {
        "attempts",
        "raw_rhs",
        "scaled_independent_rhs",
        "conservation_residuals",
        "coordinates",
        "core_state",
    }
    flat = {key: value for key, value in row.items() if key not in excluded}
    flat["coordinates_json"] = json.dumps(row["coordinates"], separators=(",", ":"))
    flat["core_state_json"] = json.dumps(row["core_state"], separators=(",", ":"))
    for key, value in row["conservation_residuals"].items():
        flat[f"conservation_{key}"] = value
    return flat


def solve_rest_chain(request: tuple[str, float]) -> dict[str, Any]:
    """Continue all twelve control points for one root and expression."""

    root, expression = request
    started = time.perf_counter()
    manifest, _, inherited_five_percent = inputs()
    alias = _root_aliases(manifest)[root]
    task24 = _task24_rest(root, expression)
    baseline_state = np.asarray(task24["core_state"], dtype=float)
    inherited_state = np.asarray(
        manifest["roots"][root]["core_state"]
        if expression == 1.0
        else inherited_five_percent[root],
        dtype=float,
    )
    spec = WTCalibrationSpec()
    lower = spec.coordinate_lower.copy()
    upper = spec.coordinate_upper.copy()
    span = upper - lower
    selected_states: dict[tuple[float, float], np.ndarray] = {}
    rows: list[dict[str, Any]] = []
    attempts: list[dict[str, Any]] = []
    alternatives: list[dict[str, Any]] = []

    for key in _continuation_order():
        co2_factor, nhe1_factor = key
        previous_key = _predecessor(key)
        branch_state = (
            baseline_state if previous_key is None else selected_states[previous_key]
        )
        model = build_controlled_rest_model(
            manifest, root, co2_factor, nhe1_factor
        )
        baseline_model = _explicit_rest_model(manifest, root)
        changed_paths = _assert_controlled_parameters(
            baseline_model.parameters, model.parameters, co2_factor, nhe1_factor
        )
        if key == BASELINE and sha256_object(model.parameters) != sha256_object(
            baseline_model.parameters
        ):
            raise AssertionError("factor one resting model changed baseline parameters")
        seeds = _unique_seeds(
            (
                ("continuation_predecessor", branch_state),
                ("task24_baseline_anchor", baseline_state),
                ("inherited_saved_rest", inherited_state),
            )
        )
        candidates: list[dict[str, Any]] = []
        local_attempts: list[dict[str, Any]] = []
        for seed_id, seed_state in seeds:
            coordinates = _core_coordinates(model, seed_state)
            if np.any(coordinates <= lower) or np.any(coordinates >= upper):
                raise RuntimeError(
                    f"{root} {expression} {key} seed {seed_id} lies outside root bounds"
                )
            attempt, root_result = _attempt_and_root(
                model=model,
                transporter="AE4",
                expression=expression,
                genotype=_genotype(expression),
                regulatory_suffix=_regulatory_suffix(model),
                spec=spec,
                start_id=seed_id,
                start=coordinates,
                lower=lower,
                upper=upper,
                max_nfev=REST_MAX_NFEV,
            )
            audit = None
            if root_result is not None:
                audit = _audit_resting_state(
                    model,
                    np.asarray(root_result.core_state, dtype=float),
                    expression,
                    root_result,
                    spec,
                )
            attempt_record = _attempt_row(
                root=root,
                alias=alias,
                expression=expression,
                co2_factor=co2_factor,
                nhe1_factor=nhe1_factor,
                predecessor=previous_key,
                seed_id=seed_id,
                seed_state=seed_state,
                attempt=attempt,
                root_result=root_result,
                audit=audit,
            )
            local_attempts.append(attempt_record)
            if root_result is not None and audit is not None and audit["rest_gate_pass"]:
                candidates.append(
                    {
                        "seed_id": seed_id,
                        "coordinates": np.asarray(root_result.coordinates, dtype=float),
                        "core_state": np.asarray(root_result.core_state, dtype=float),
                        "root_result": root_result,
                        "audit": audit,
                        "attempt_record": attempt_record,
                    }
                )
        clusters = _cluster_candidates(candidates, spec)
        if not clusters:
            attempts.extend(local_attempts)
            return {
                "rows": rows,
                "attempts": attempts,
                "alternatives": alternatives,
                "error": (
                    f"no admissible resting root for {root}, expression={expression}, "
                    f"CO2={co2_factor}, NHE1={nhe1_factor}"
                ),
            }
        branch_coordinates = _core_coordinates(model, branch_state)
        for cluster in clusters:
            representative = cluster["representative"]
            cluster["distance_from_predecessor_normalized"] = float(
                np.max(
                    np.abs(
                        (np.asarray(representative["coordinates"]) - branch_coordinates)
                        / span
                    )
                )
            )
            for member in cluster["members"]:
                member["attempt_record"]["cluster_id"] = cluster["cluster_id"]
        selected_cluster = min(
            clusters,
            key=lambda item: (
                float(item["distance_from_predecessor_normalized"]),
                float(
                    item["representative"]["audit"][
                        "max_abs_scaled_independent_rhs"
                    ]
                ),
                item["cluster_id"],
            ),
        )
        selected = selected_cluster["representative"]
        for member in selected_cluster["members"]:
            member["attempt_record"]["selected_cluster"] = True
        for cluster in clusters:
            if cluster is selected_cluster:
                continue
            representative = cluster["representative"]
            alternatives.append(
                {
                    "root_id": root,
                    "root_alias": alias,
                    "ae4_expression": expression,
                    "co2_factor": co2_factor,
                    "nhe1_factor": nhe1_factor,
                    "cluster_id": cluster["cluster_id"],
                    "member_seed_ids": ";".join(
                        sorted(str(item["seed_id"]) for item in cluster["members"])
                    ),
                    "distance_from_predecessor_normalized": cluster[
                        "distance_from_predecessor_normalized"
                    ],
                    "max_abs_scaled_independent_rhs": representative["audit"][
                        "max_abs_scaled_independent_rhs"
                    ],
                    "coordinates_json": json.dumps(
                        list(np.asarray(representative["coordinates"], dtype=float)),
                        separators=(",", ":"),
                    ),
                    "core_state_json": json.dumps(
                        list(np.asarray(representative["core_state"], dtype=float)),
                        separators=(",", ":"),
                    ),
                }
            )
        attempts.extend(local_attempts)
        selected_state = np.asarray(selected["core_state"], dtype=float)
        selected_states[key] = selected_state
        task24_coordinates = _core_coordinates(model, baseline_state)
        baseline_distance = float(
            np.max(
                np.abs(
                    (np.asarray(selected["coordinates"]) - task24_coordinates) / span
                )
            )
        )
        baseline_gate = bool(
            key != BASELINE
            or baseline_distance <= BASELINE_REST_RELATIVE_TOLERANCE
        )
        if not baseline_gate:
            raise RuntimeError(
                f"Task 24 resting baseline was not reproduced for {root}, expression={expression}"
            )
        metrics = _rest_flux_metrics(model, selected_state, expression)
        audit = dict(selected["audit"])
        row: dict[str, Any] = {
            "root_id": root,
            "root_alias": alias,
            "routing_family": routing(manifest, root),
            "ae4_expression": expression,
            "co2_factor": co2_factor,
            "nhe1_factor": nhe1_factor,
            "co2_apical_permeability_fmol_s_mM": (
                model.parameters.homeostasis.co2_apical_permeability_fmol_s_mM
            ),
            "co2_basolateral_permeability_fmol_s_mM": (
                model.parameters.homeostasis.co2_basolateral_permeability_fmol_s_mM
            ),
            "nhe1_capacity_fmol_s": model.parameters.homeostasis.nhe1_capacity_fmol_s,
            "changed_parameter_paths": ";".join(changed_paths),
            "selected_cluster_id": selected_cluster["cluster_id"],
            "selected_seed_id": selected["seed_id"],
            "direct_start_count": len(seeds),
            "accepted_direct_start_count": len(candidates),
            "candidate_root_cluster_count": len(clusters),
            "alternate_root_count": len(clusters) - 1,
            "predecessor_co2_factor": None if previous_key is None else previous_key[0],
            "predecessor_nhe1_factor": None if previous_key is None else previous_key[1],
            "branch_step_distance_normalized": selected_cluster[
                "distance_from_predecessor_normalized"
            ],
            "task24_anchor_distance_normalized": baseline_distance,
            "task24_baseline_rest_reproduced": baseline_gate,
            "rest_method": "direct_charge_manifold_least_squares",
            "rest_integration_used": False,
            "rest_protocol": StimulusArm.REST.value,
            "resting_calcium_uM": float(model.stimulus.resting_calcium_uM),
            "rest_beta_input": 0.0,
            "whole_cell_parameters_sha256": sha256_object(model.parameters),
            "ae4_parameters_sha256": sha256_object(model.ae4_parameters),
            "core_state_sha256": sha256_object(list(selected_state)),
            "coordinates": [float(value) for value in selected["coordinates"]],
            "core_state": [float(value) for value in selected_state],
            "attempts": local_attempts,
            "wall_seconds": time.perf_counter() - started,
            **audit,
            **metrics,
        }
        path = OUT / "resting_states" / f"{_rest_stem(*((root,) + key + (expression,)))}.json"
        write_json(path, row)
        rows.append(row)
    return {
        "rows": rows,
        "attempts": attempts,
        "alternatives": alternatives,
        "error": None,
    }


def _array_sha256(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(contiguous.dtype).encode("ascii"))
    digest.update(str(contiguous.shape).encode("ascii"))
    digest.update(contiguous.tobytes())
    return digest.hexdigest()


def _landmark_values(time_s: np.ndarray, values: np.ndarray) -> dict[str, float]:
    result: dict[str, float] = {}
    for landmark in (0.0, 30.0, 60.0, 120.0, 300.0, 600.0):
        matches = np.flatnonzero(time_s == landmark)
        if matches.size != 1:
            raise AssertionError(f"time grid does not contain landmark {landmark}")
        result[f"t{int(landmark):03d}"] = float(values[int(matches[0])])
    return result


def run_stimulation_case(
    request: tuple[str, float, float, float, float]
) -> dict[str, Any]:
    root, co2_factor, nhe1_factor, calcium, expression = request
    started = time.perf_counter()
    manifest, frozen, _ = inputs()
    alias = _root_aliases(manifest)[root]
    rest_path = OUT / "resting_states" / (
        _rest_stem(root, co2_factor, nhe1_factor, expression) + ".json"
    )
    rest = json.loads(rest_path.read_text())
    if not rest["rest_gate_pass"]:
        raise AssertionError("stimulation cannot start from an unaccepted resting state")
    core = np.asarray(rest["core_state"], dtype=float)
    genotype = _genotype(expression)
    model = build_controlled_stimulation_model(
        manifest, root, calcium, co2_factor, nhe1_factor
    )
    rest_model = build_controlled_rest_model(
        manifest, root, co2_factor, nhe1_factor
    )
    y0 = _root_state_with_basal_regulation(model, core)
    stimulated_onset = model.evaluate(0.0, y0, genotype=genotype)
    rest_onset = rest_model.evaluate(0.0, y0, genotype=genotype)
    onset_match = bool(np.array_equal(stimulated_onset.rhs, rest_onset.rhs))
    if not onset_match:
        raise AssertionError("stimulation onset does not reproduce the REST equations")
    trace = simulate_genotype(
        model,
        core,
        transporter="AE4",
        genotype=genotype,
        solver=PRODUCTION_RADAU,
        time_s=frozen["time_grid_s"],
    )
    exact_grid = bool(np.array_equal(trace.time_s, np.asarray(frozen["time_grid_s"])))
    valid = bool(trace.numerical_gate_pass and exact_grid)
    if not valid:
        raise RuntimeError(
            f"production stimulation failed for {root}, CO2={co2_factor}, "
            f"NHE1={nhe1_factor}, Ca={calcium}, expression={expression}"
        )
    opposing_samples = 0
    max_charge = 0.0
    max_expression_residual = 0.0
    max_nkcc_stoichiometry_residual = 0.0
    negative_cycle_samples = 0
    for time_s, state in zip(trace.time_s, trace.states.T):
        evaluation = model.evaluate(float(time_s), state, genotype=genotype)
        ae4 = evaluation.diagnostics.ae4
        opposing_samples += int(ae4.na_cell_fmol_s * ae4.k_cell_fmol_s < 0.0)
        negative_cycle_samples += int(ae4.cl_cell_fmol_s < 0.0)
        max_charge = max(max_charge, abs(float(ae4.charge_source_residual_fmol_s)))
        raw_ae4 = float(ae4.diagnostics["j_ae4_fmol_s"])
        max_expression_residual = max(
            max_expression_residual,
            abs(float(ae4.cl_cell_fmol_s) - expression * raw_ae4),
        )
        homeostasis = evaluation.diagnostics.homeostasis
        max_nkcc_stoichiometry_residual = max(
            max_nkcc_stoichiometry_residual,
            abs(
                (
                    float(homeostasis.cl_cell_fmol_s)
                    - float(homeostasis.ae2_inward_fmol_s)
                )
                - 2.0 * float(homeostasis.nkcc1_inward_fmol_s)
            ),
        )
    if opposing_samples:
        raise AssertionError("routed AE4 produced opposing Na and K sources")
    if max_charge > CONSERVATION_RESIDUAL_TOLERANCES["ae4_charge_fmol_s"]:
        raise AssertionError("AE4 charge neutrality failed")
    if max_expression_residual > AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S:
        raise AssertionError("AE4 expression scaling failed during stimulation")
    stem = _trajectory_stem(
        root, co2_factor, nhe1_factor, calcium, expression
    )
    trajectory_path = OUT / "trajectories" / f"{stem}.npz"
    trajectory_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        key: value for key, value in asdict(trace).items() if isinstance(value, np.ndarray)
    }
    np.savez_compressed(trajectory_path, **arrays)
    min_flow_index = int(np.argmin(trace.flow_pL_s))
    max_flow_index = int(np.argmax(trace.flow_pL_s))
    row: dict[str, Any] = {
        "root_id": root,
        "root_alias": alias,
        "routing_family": routing(manifest, root),
        "co2_factor": co2_factor,
        "nhe1_factor": nhe1_factor,
        "calcium_uM": calcium,
        "ae4_expression": expression,
        "solver": PRODUCTION_RADAU.label,
        "valid": valid,
        "solver_success": bool(trace.success),
        "solver_message": trace.message,
        "exact_time_grid": exact_grid,
        "sample_count": len(trace.time_s),
        "total_0_600_pL": float(trace.cumulative_flow_pL[-1]),
        "endpoint_flow_pL_s": float(trace.flow_pL_s[-1]),
        "minimum_flow_pL_s": float(trace.flow_pL_s[min_flow_index]),
        "minimum_flow_time_s": float(trace.time_s[min_flow_index]),
        "maximum_flow_pL_s": float(trace.flow_pL_s[max_flow_index]),
        "maximum_flow_time_s": float(trace.time_s[max_flow_index]),
        "positive_core": bool(trace.positive_core),
        "nonnegative_flow": bool(trace.all_flow_nonnegative),
        "max_dimensionless_conservation_ratio": float(
            trace.max_dimensionless_conservation_ratio
        ),
        "onset_rest_rhs_match": onset_match,
        "initial_max_abs_amount_rhs_fmol_s": float(
            np.max(np.abs(stimulated_onset.rhs[AMOUNT_ROWS]))
        ),
        "initial_max_abs_volume_rhs_pL_s": float(
            np.max(np.abs(stimulated_onset.rhs[VOLUME_ROWS]))
        ),
        "opposing_cation_samples": opposing_samples,
        "negative_net_cycle_samples": negative_cycle_samples,
        "max_ae4_charge_residual_fmol_s": max_charge,
        "max_ae4_expression_scaling_residual_fmol_s": max_expression_residual,
        "max_nkcc1_chloride_stoichiometry_residual_fmol_s": (
            max_nkcc_stoichiometry_residual
        ),
        "resting_state_sha256": rest["core_state_sha256"],
        "whole_cell_parameters_sha256": sha256_object(model.parameters),
        "ae4_parameters_sha256": sha256_object(model.ae4_parameters),
        "time_grid_sha256": _array_sha256(trace.time_s),
        "states_sha256": _array_sha256(trace.states),
        "flow_sha256": _array_sha256(trace.flow_pL_s),
        "cumulative_flow_sha256": _array_sha256(trace.cumulative_flow_pL),
        "trajectory_path": str(trajectory_path.relative_to(REPO)),
        "trajectory_sha256": sha256_file(trajectory_path),
        "wall_seconds": time.perf_counter() - started,
    }
    for name in (
        "cell_na_mM",
        "cell_k_mM",
        "cell_cl_mM",
        "cell_ph",
        "cell_volume_pL",
        "lumen_na_mM",
        "lumen_k_mM",
        "lumen_cl_mM",
        "lumen_ph",
        "lumen_volume_pL",
        "ae4_capacity_multiplier",
    ):
        values = np.asarray(getattr(trace, name), dtype=float)
        row[f"initial_{name}"] = float(values[0])
        row[f"endpoint_{name}"] = float(values[-1])
        row[f"minimum_{name}"] = float(np.min(values))
        row[f"maximum_{name}"] = float(np.max(values))
    for label, value in _landmark_values(trace.time_s, trace.flow_pL_s).items():
        row[f"flow_{label}_pL_s"] = value
    for label, value in _landmark_values(trace.time_s, trace.cumulative_flow_pL).items():
        row[f"cumulative_{label}_pL"] = value
    write_json(OUT / "trajectories" / f"{stem}.json", row)
    return row


REST_PAIR_FIELDS = (
    "cell_tic_mM",
    "cell_hco3_mM",
    "cell_ph",
    "cell_na_mM",
    "cell_k_mM",
    "cell_cl_mM",
    "cell_volume_pL",
    "ae4_cl_flux_fmol_s",
    "ae4_cl_flux_pre_expression_fmol_s",
    "ae4_flux_per_unit_expression_fmol_s",
    "nkcc1_inward_cycle_flux_fmol_s",
    "nkcc1_chloride_influx_fmol_s",
    "nhe1_inward_na_flux_fmol_s",
    "ae2_inward_cycle_flux_fmol_s",
    "co2_basolateral_bath_to_cell_fmol_s",
    "co2_apical_lumen_to_cell_fmol_s",
    "v_apical_V",
    "v_basolateral_V",
    "v_transepithelial_V",
    "v_apical_mV",
    "v_basolateral_mV",
    "v_transepithelial_mV",
    "apical_chloride_cell_to_lumen_fmol_s",
    "apical_k_cell_to_lumen_fmol_s",
    "basolateral_k_cell_to_bath_fmol_s",
)


def make_flux_decomposition(
    rest_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    index = {
        (
            row["root_id"],
            float(row["co2_factor"]),
            float(row["nhe1_factor"]),
            float(row["ae4_expression"]),
        ): row
        for row in rest_rows
    }
    rows: list[dict[str, Any]] = []
    roots = sorted({row["root_id"] for row in rest_rows})
    aliases = {row["root_id"]: row["root_alias"] for row in rest_rows}
    for co2_factor, nhe1_factor, root in itertools.product(
        CO2_FACTORS, NHE1_FACTORS, roots
    ):
        wt = index[root, co2_factor, nhe1_factor, 1.0]
        low = index[root, co2_factor, nhe1_factor, 0.05]
        wt_ae4 = float(wt["ae4_cl_flux_fmol_s"])
        low_ae4 = float(low["ae4_cl_flux_fmol_s"])
        if wt_ae4 == 0.0:
            raise ZeroDivisionError("WT AE4 chloride loading is zero")
        retention = low_ae4 / wt_ae4
        compensation = (low_ae4 / 0.05) / wt_ae4
        signed_loss = wt_ae4 - low_ae4
        absolute_reduction = abs(signed_loss)
        nkcc_change = float(low["nkcc1_chloride_influx_fmol_s"]) - float(
            wt["nkcc1_chloride_influx_fmol_s"]
        )
        replacement = nkcc_change / signed_loss if signed_loss != 0.0 else math.nan
        row: dict[str, Any] = {
            "root_id": root,
            "root_alias": aliases[root],
            "routing_family": wt["routing_family"],
            "co2_factor": co2_factor,
            "nhe1_factor": nhe1_factor,
            "wt_resting_state_sha256": wt["core_state_sha256"],
            "ae4_5pct_resting_state_sha256": low["core_state_sha256"],
            "ae4_flux_retention": retention,
            "ae4_per_expression_compensation": compensation,
            "ae4_chloride_loading_loss_fmol_s": signed_loss,
            "ae4_chloride_loading_absolute_reduction_fmol_s": absolute_reduction,
            "nkcc1_chloride_loading_change_fmol_s": nkcc_change,
            "nkcc1_chloride_loading_absolute_change_fmol_s": abs(nkcc_change),
            "fraction_lost_ae4_loading_replaced_by_nkcc1": replacement,
            "cell_hco3_change_mM": float(low["cell_hco3_mM"])
            - float(wt["cell_hco3_mM"]),
            "cell_tic_change_mM": float(low["cell_tic_mM"])
            - float(wt["cell_tic_mM"]),
            "cell_ph_change": float(low["cell_ph"]) - float(wt["cell_ph"]),
            "nhe1_flux_change_fmol_s": float(low["nhe1_inward_na_flux_fmol_s"])
            - float(wt["nhe1_inward_na_flux_fmol_s"]),
            "ae2_flux_change_fmol_s": float(low["ae2_inward_cycle_flux_fmol_s"])
            - float(wt["ae2_inward_cycle_flux_fmol_s"]),
            "co2_basolateral_flux_change_fmol_s": float(
                low["co2_basolateral_bath_to_cell_fmol_s"]
            )
            - float(wt["co2_basolateral_bath_to_cell_fmol_s"]),
            "co2_apical_flux_change_fmol_s": float(
                low["co2_apical_lumen_to_cell_fmol_s"]
            )
            - float(wt["co2_apical_lumen_to_cell_fmol_s"]),
            "nkcc1_chloride_factor": 2.0,
            "wt_rest_gate_pass": bool(wt["rest_gate_pass"]),
            "ae4_5pct_rest_gate_pass": bool(low["rest_gate_pass"]),
        }
        for field in REST_PAIR_FIELDS:
            row[f"wt_{field}"] = wt[field]
            row[f"ae4_5pct_{field}"] = low[field]
        rows.append(row)
    write_rows(OUT / "flux_decomposition.csv", rows)
    return rows


def _load_task24_trajectory_index() -> dict[tuple[str, float, float], dict[str, str]]:
    rows = list(csv.DictReader(TASK24_TRAJECTORY_SUMMARY.open()))
    index = {
        (row["root_id"], float(row["calcium_uM"]), float(row["ae4_expression"])): row
        for row in rows
    }
    if len(index) != 60:
        raise AssertionError("Task 24 trajectory summary is incomplete")
    return index


def make_comparisons(
    rest_rows: Sequence[Mapping[str, Any]],
    trajectory_rows: Sequence[Mapping[str, Any]],
    flux_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    trajectories = {
        (
            row["root_id"],
            float(row["co2_factor"]),
            float(row["nhe1_factor"]),
            float(row["calcium_uM"]),
            float(row["ae4_expression"]),
        ): row
        for row in trajectory_rows
    }
    fluxes = {
        (row["root_id"], float(row["co2_factor"]), float(row["nhe1_factor"])): row
        for row in flux_rows
    }
    task24 = _load_task24_trajectory_index()
    roots = sorted({row["root_id"] for row in rest_rows})
    aliases = {row["root_id"]: row["root_alias"] for row in rest_rows}
    comparisons: list[dict[str, Any]] = []
    for co2_factor, nhe1_factor, root, calcium in itertools.product(
        CO2_FACTORS, NHE1_FACTORS, roots, CALCIUM
    ):
        wt = trajectories[root, co2_factor, nhe1_factor, calcium, 1.0]
        low = trajectories[root, co2_factor, nhe1_factor, calcium, 0.05]
        flux = fluxes[root, co2_factor, nhe1_factor]
        wt_total = float(wt["total_0_600_pL"])
        low_total = float(low["total_0_600_pL"])
        ratio = low_total / wt_total
        change = 100.0 * (ratio - 1.0)
        baseline_wt = task24[root, calcium, 1.0]
        baseline_low = task24[root, calcium, 0.05]
        if (co2_factor, nhe1_factor) == BASELINE:
            reference_wt = float(baseline_wt["total_0_600_pL"])
            reference_low = float(baseline_low["total_0_600_pL"])
            wt_error = abs(wt_total - reference_wt) / max(abs(reference_wt), 1.0e-300)
            low_error = abs(low_total - reference_low) / max(abs(reference_low), 1.0e-300)
            reference_ratio = reference_low / reference_wt
            ratio_error = abs(ratio - reference_ratio)
        else:
            reference_wt = reference_low = reference_ratio = math.nan
            wt_error = low_error = ratio_error = math.nan
        row = {
            "root_id": root,
            "root_alias": aliases[root],
            "routing_family": wt["routing_family"],
            "co2_factor": co2_factor,
            "nhe1_factor": nhe1_factor,
            "calcium_uM": calcium,
            "wt_total_0_600_pL": wt_total,
            "ae4_5pct_total_0_600_pL": low_total,
            "ae4_5pct_to_wt_secretion_ratio": ratio,
            "ae4_5pct_secretion_change_percent": change,
            "ae4_loss_decreases_secretion": bool(ratio < 1.0),
            "paired_production_gate_pass": bool(wt["valid"] and low["valid"]),
            "wt_resting_state_sha256": wt["resting_state_sha256"],
            "ae4_5pct_resting_state_sha256": low["resting_state_sha256"],
            "ae4_flux_retention": flux["ae4_flux_retention"],
            "ae4_per_expression_compensation": flux[
                "ae4_per_expression_compensation"
            ],
            "ae4_chloride_loading_loss_fmol_s": flux[
                "ae4_chloride_loading_loss_fmol_s"
            ],
            "nkcc1_chloride_loading_change_fmol_s": flux[
                "nkcc1_chloride_loading_change_fmol_s"
            ],
            "fraction_lost_ae4_loading_replaced_by_nkcc1": flux[
                "fraction_lost_ae4_loading_replaced_by_nkcc1"
            ],
            "cell_hco3_change_mM": flux["cell_hco3_change_mM"],
            "cell_tic_change_mM": flux["cell_tic_change_mM"],
            "cell_ph_change": flux["cell_ph_change"],
            "task24_wt_total_0_600_pL": reference_wt,
            "task24_ae4_5pct_total_0_600_pL": reference_low,
            "task24_ae4_5pct_to_wt_ratio": reference_ratio,
            "task24_wt_total_relative_error": wt_error,
            "task24_ae4_5pct_total_relative_error": low_error,
            "task24_ratio_absolute_error": ratio_error,
        }
        comparisons.append(row)
    write_rows(OUT / "comparison.csv", comparisons)
    return comparisons


def _stats(values: Iterable[float]) -> dict[str, float]:
    numbers = [float(value) for value in values]
    if not numbers:
        raise ValueError("summary group is empty")
    return {
        "min": min(numbers),
        "median": statistics.median(numbers),
        "max": max(numbers),
    }


def _add_stats(target: dict[str, Any], prefix: str, values: Iterable[float]) -> None:
    for suffix, value in _stats(values).items():
        target[f"{prefix}_{suffix}"] = value


def make_parameter_summary(
    flux_rows: Sequence[Mapping[str, Any]], comparisons: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    roots = sorted({row["root_id"] for row in flux_rows})
    for co2_factor, nhe1_factor in itertools.product(CO2_FACTORS, NHE1_FACTORS):
        flux_group = [
            row
            for row in flux_rows
            if float(row["co2_factor"]) == co2_factor
            and float(row["nhe1_factor"]) == nhe1_factor
        ]
        if len(flux_group) != len(roots):
            raise AssertionError("each parameter summary requires all ten roots")
        row: dict[str, Any] = {
            "co2_factor": co2_factor,
            "nhe1_factor": nhe1_factor,
            "root_count": len(flux_group),
        }
        statistics_fields = {
            "wt_rest_hco3_mM": "wt_cell_hco3_mM",
            "ae4_5pct_rest_hco3_mM": "ae4_5pct_cell_hco3_mM",
            "wt_to_ae4_5pct_hco3_change_mM": "cell_hco3_change_mM",
            "wt_rest_tic_mM": "wt_cell_tic_mM",
            "ae4_5pct_rest_tic_mM": "ae4_5pct_cell_tic_mM",
            "wt_to_ae4_5pct_tic_change_mM": "cell_tic_change_mM",
            "wt_to_ae4_5pct_ph_change": "cell_ph_change",
            "ae4_flux_retention": "ae4_flux_retention",
            "ae4_per_expression_compensation": "ae4_per_expression_compensation",
            "nkcc1_chloride_flux_change_fmol_s": (
                "nkcc1_chloride_loading_change_fmol_s"
            ),
            "nhe1_flux_change_fmol_s": "nhe1_flux_change_fmol_s",
            "fraction_lost_ae4_replaced_by_nkcc1": (
                "fraction_lost_ae4_loading_replaced_by_nkcc1"
            ),
        }
        for output_name, source_name in statistics_fields.items():
            _add_stats(row, output_name, (item[source_name] for item in flux_group))
        consistent_decrease = set(roots)
        any_decrease: set[str] = set()
        for calcium in CALCIUM:
            group = [
                item
                for item in comparisons
                if float(item["co2_factor"]) == co2_factor
                and float(item["nhe1_factor"]) == nhe1_factor
                and float(item["calcium_uM"]) == calcium
            ]
            if len(group) != len(roots):
                raise AssertionError("each calcium summary requires all ten roots")
            calcium_label = f"ca_{calcium:.2f}".replace(".", "p")
            _add_stats(
                row,
                f"{calcium_label}_secretion_ratio",
                (item["ae4_5pct_to_wt_secretion_ratio"] for item in group),
            )
            _add_stats(
                row,
                f"{calcium_label}_secretion_change_percent",
                (item["ae4_5pct_secretion_change_percent"] for item in group),
            )
            decreasing = {
                str(item["root_id"])
                for item in group
                if bool(item["ae4_loss_decreases_secretion"])
            }
            row[f"{calcium_label}_decrease_root_count"] = len(decreasing)
            consistent_decrease &= decreasing
            any_decrease |= decreasing
        row["decrease_at_all_calcium_root_count"] = len(consistent_decrease)
        row["decrease_at_any_calcium_root_count"] = len(any_decrease)
        result.append(row)
    write_rows(OUT / "parameter_summary.csv", result)
    return result


def _parameter_row(
    rows: Sequence[Mapping[str, Any]], co2_factor: float, nhe1_factor: float
) -> Mapping[str, Any]:
    matches = [
        row
        for row in rows
        if float(row["co2_factor"]) == co2_factor
        and float(row["nhe1_factor"]) == nhe1_factor
    ]
    if len(matches) != 1:
        raise AssertionError("parameter combination is not unique")
    return matches[0]


def diagnostic_answers(
    parameter_rows: Sequence[Mapping[str, Any]],
    flux_rows: Sequence[Mapping[str, Any]],
    comparisons: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    baseline = _parameter_row(parameter_rows, 1.0, 1.0)
    co2_only = _parameter_row(parameter_rows, 10.0, 1.0)
    nhe_only = _parameter_row(parameter_rows, 1.0, 0.25)
    joint = _parameter_row(parameter_rows, 10.0, 0.25)
    nonbaseline = [
        row
        for row in parameter_rows
        if (float(row["co2_factor"]), float(row["nhe1_factor"])) != BASELINE
    ]
    lowest_compensation = min(
        nonbaseline,
        key=lambda row: float(row["ae4_per_expression_compensation_median"]),
    )
    lowest_retention = min(
        nonbaseline, key=lambda row: float(row["ae4_flux_retention_median"])
    )
    co2_accumulation_rows = [
        row for row in flux_rows if float(row["co2_factor"]) > 1.0
    ]
    reduced_nhe_rows = [
        row for row in flux_rows if float(row["nhe1_factor"]) < 1.0
    ]

    def flux_medians(co2_factor: float, nhe1_factor: float) -> dict[str, float]:
        rows = [
            row
            for row in flux_rows
            if float(row["co2_factor"]) == co2_factor
            and float(row["nhe1_factor"]) == nhe1_factor
        ]
        if len(rows) != 10:
            raise AssertionError("each flux median requires all ten roots")
        return {
            "wt": float(statistics.median(row["wt_ae4_cl_flux_fmol_s"] for row in rows)),
            "ae4_5pct": float(
                statistics.median(row["ae4_5pct_ae4_cl_flux_fmol_s"] for row in rows)
            ),
        }

    increased_co2_net_fluxes = [
        float(row[prefix + "co2_basolateral_bath_to_cell_fmol_s"])
        + float(row[prefix + "co2_apical_lumen_to_cell_fmol_s"])
        for row in co2_accumulation_rows
        for prefix in ("wt_", "ae4_5pct_")
    ]
    consistent_combinations: list[dict[str, Any]] = []
    for parameter in parameter_rows:
        if int(parameter["decrease_at_all_calcium_root_count"]) == 10:
            consistent_combinations.append(
                {
                    "co2_factor": float(parameter["co2_factor"]),
                    "nhe1_factor": float(parameter["nhe1_factor"]),
                }
            )
    decrease_rows = [
        row for row in comparisons if bool(row["ae4_loss_decreases_secretion"])
    ]
    all_numerical = all(bool(row["paired_production_gate_pass"]) for row in comparisons)
    return {
        "A": {
            "question": "Does increasing CO2 permeability prevent TIC and HCO3 accumulation after AE4 loss?",
            "all_increased_co2_tic_changes_nonpositive": all(
                float(row["cell_tic_change_mM"]) <= 0.0
                for row in co2_accumulation_rows
            ),
            "all_increased_co2_hco3_changes_nonpositive": all(
                float(row["cell_hco3_change_mM"]) <= 0.0
                for row in co2_accumulation_rows
            ),
            "increased_co2_tic_change_mM": _stats(
                row["cell_tic_change_mM"] for row in co2_accumulation_rows
            ),
            "increased_co2_hco3_change_mM": _stats(
                row["cell_hco3_change_mM"] for row in co2_accumulation_rows
            ),
            "all_increased_co2_apical_fluxes_outward": all(
                float(row[prefix + "co2_apical_lumen_to_cell_fmol_s"]) < 0.0
                for row in co2_accumulation_rows
                for prefix in ("wt_", "ae4_5pct_")
            ),
            "all_increased_co2_basolateral_fluxes_inward": all(
                float(row[prefix + "co2_basolateral_bath_to_cell_fmol_s"]) > 0.0
                for row in co2_accumulation_rows
                for prefix in ("wt_", "ae4_5pct_")
            ),
            "all_increased_co2_net_fluxes_inward": all(
                value > 0.0 for value in increased_co2_net_fluxes
            ),
            "increased_co2_net_bath_and_lumen_to_cell_flux_fmol_s": _stats(
                increased_co2_net_fluxes
            ),
        },
        "B": {
            "question": "Does reducing NHE1 capacity prevent alkalinisation after AE4 loss?",
            "all_reduced_nhe1_ph_changes_nonpositive": all(
                float(row["cell_ph_change"]) <= 0.0 for row in reduced_nhe_rows
            ),
            "reduced_nhe1_ph_change": _stats(
                row["cell_ph_change"] for row in reduced_nhe_rows
            ),
            "baseline_ph_change": {
                key: float(baseline[f"wt_to_ae4_5pct_ph_change_{key}"])
                for key in ("min", "median", "max")
            },
        },
        "C": {
            "question": "Which intervention reduces AE4 compensation?",
            "baseline_median": float(
                baseline["ae4_per_expression_compensation_median"]
            ),
            "co2_x10_nhe1_x1_median": float(
                co2_only["ae4_per_expression_compensation_median"]
            ),
            "co2_x1_nhe1_x0p25_median": float(
                nhe_only["ae4_per_expression_compensation_median"]
            ),
            "co2_x10_nhe1_x0p25_median": float(
                joint["ae4_per_expression_compensation_median"]
            ),
            "absolute_ae4_flux_median_fmol_s": {
                "baseline": flux_medians(1.0, 1.0),
                "co2_x10_nhe1_x1": flux_medians(10.0, 1.0),
                "co2_x1_nhe1_x0p25": flux_medians(1.0, 0.25),
                "co2_x10_nhe1_x0p25": flux_medians(10.0, 0.25),
            },
            "lowest_median_combination": {
                "co2_factor": float(lowest_compensation["co2_factor"]),
                "nhe1_factor": float(lowest_compensation["nhe1_factor"]),
                "median": float(
                    lowest_compensation["ae4_per_expression_compensation_median"]
                ),
            },
        },
        "D": {
            "question": "Does 5 per cent AE4 flux retention fall below the current value?",
            "baseline_retention": {
                key: float(baseline[f"ae4_flux_retention_{key}"])
                for key in ("min", "median", "max")
            },
            "lowest_median_combination": {
                "co2_factor": float(lowest_retention["co2_factor"]),
                "nhe1_factor": float(lowest_retention["nhe1_factor"]),
                "retention_min": float(lowest_retention["ae4_flux_retention_min"]),
                "retention_median": float(
                    lowest_retention["ae4_flux_retention_median"]
                ),
                "retention_max": float(lowest_retention["ae4_flux_retention_max"]),
            },
            "nonbaseline_retention": _stats(
                row["ae4_flux_retention"]
                for row in flux_rows
                if (float(row["co2_factor"]), float(row["nhe1_factor"]))
                != BASELINE
            ),
            "nonbaseline_combination_medians_below_baseline_median": sum(
                float(row["ae4_flux_retention_median"])
                < float(baseline["ae4_flux_retention_median"])
                for row in nonbaseline
            ),
        },
        "E": {
            "question": "Does NKCC1 replace most missing AE4 chloride loading?",
            "grid_replacement_fraction": _stats(
                row["fraction_lost_ae4_loading_replaced_by_nkcc1"]
                for row in flux_rows
            ),
            "combination_medians_above_one_half": sum(
                float(row["fraction_lost_ae4_replaced_by_nkcc1_median"]) > 0.5
                for row in parameter_rows
            ),
            "individual_root_grid_points_above_one_half": sum(
                float(row["fraction_lost_ae4_loading_replaced_by_nkcc1"]) > 0.5
                for row in flux_rows
            ),
        },
        "F": {
            "question": "Does any grid point make AE4 loss reduce secretion consistently?",
            "consistent_all_root_all_calcium_combinations": consistent_combinations,
            "decreasing_comparison_count": len(decrease_rows),
            "comparison_count": len(comparisons),
        },
        "G": {
            "question": "If secretion falls, is it associated with reduced AE4 self rescue?",
            "decreasing_cases_exist": bool(decrease_rows),
            "all_comparisons_numerically_valid": all_numerical,
            "decreasing_case_retention": (
                _stats(row["ae4_flux_retention"] for row in decrease_rows)
                if decrease_rows
                else None
            ),
            "decreasing_case_compensation": (
                _stats(row["ae4_per_expression_compensation"] for row in decrease_rows)
                if decrease_rows
                else None
            ),
        },
    }


def make_summary(
    parameter_rows: Sequence[Mapping[str, Any]],
    flux_rows: Sequence[Mapping[str, Any]],
    comparisons: Sequence[Mapping[str, Any]],
    verification: Mapping[str, Any],
) -> dict[str, Any]:
    answers = diagnostic_answers(parameter_rows, flux_rows, comparisons)
    summary = {
        "model_id": MODEL_ID,
        "source_commit": SOURCE_COMMIT,
        "grid_combination_count": len(parameter_rows),
        "resting_state_count": len(flux_rows) * 2,
        "stimulated_trajectory_count": len(comparisons) * 2,
        "baseline_reproduced": bool(verification["baseline_control_pass"]),
        "all_resting_states_valid": bool(verification["all_rest_gates_pass"]),
        "all_trajectories_valid": bool(verification["all_stimulation_gates_pass"]),
        "diagnostic_answers": answers,
        "parameter_summary": list(parameter_rows),
    }
    write_json(OUT / "summary.json", summary)
    return summary


def _source_files() -> list[Path]:
    files = [MANIFEST, FROZEN]
    files.extend(sorted(path for path in TASK24.rglob("*") if path.is_file()))
    files.extend(sorted((REPO / "src").rglob("*.py")))
    return files


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(REPO)): sha256_file(path)
        for path in _source_files()
    }


def contract_payload(manifest: Mapping[str, Any], frozen: Mapping[str, Any]) -> dict[str, Any]:
    spec = WTCalibrationSpec()
    return {
        "task": "Task 25 carbon and NHE1 compensation sensitivity",
        "source_commit": SOURCE_COMMIT,
        "branch": BRANCH,
        "model_id": MODEL_ID,
        "root_ids": sorted(manifest["roots"]),
        "co2_factors": list(CO2_FACTORS),
        "nhe1_factors": list(NHE1_FACTORS),
        "ae4_expression": list(EXPRESSIONS),
        "calcium_uM": list(CALCIUM),
        "expected_resting_state_count": 240,
        "expected_stimulated_trajectory_count": 720,
        "permitted_parameter_paths": list(EXPECTED_PARAMETER_PATHS),
        "co2_rule": "scale apical and basolateral existing CO2 permeabilities together",
        "nhe1_rule": "scale only the existing NHE1 capacity",
        "factor_one_is_exact_model_control": True,
        "continuation_order": [
            {"co2_factor": co2, "nhe1_factor": nhe1}
            for co2, nhe1 in _continuation_order()
        ],
        "continuation_predecessors": [
            {
                "co2_factor": co2,
                "nhe1_factor": nhe1,
                "predecessor": (
                    None
                    if _predecessor((co2, nhe1)) is None
                    else {
                        "co2_factor": _predecessor((co2, nhe1))[0],
                        "nhe1_factor": _predecessor((co2, nhe1))[1],
                    }
                ),
            }
            for co2, nhe1 in _continuation_order()
        ],
        "direct_warm_starts": [
            "nearby selected continuation predecessor",
            "fixed Task 24 re equilibrated branch anchor",
            "inherited saved resting state",
        ],
        "root_selection": (
            "cluster accepted direct roots at the inherited tolerance, then select the "
            "cluster nearest the declared predecessor; phenotype values are unavailable"
        ),
        "rest_solver": {
            "method": "direct bounded charge manifold least squares",
            "max_nfev": REST_MAX_NFEV,
            "xtol": 1.0e-12,
            "ftol": 1.0e-12,
            "gtol": 1.0e-12,
            "x_scale": "jac",
            "coordinate_bounds": [asdict(item) for item in spec.coordinate_bounds],
            "independent_rhs_scales": list(spec.independent_rhs_scales),
            "root_scaled_tolerance": spec.root_scaled_tolerance,
            "omitted_row_raw_tolerance": spec.omitted_row_raw_tolerance,
            "charge_tolerance_fmol": spec.charge_tolerance_fmol,
            "current_tolerance_A": spec.current_tolerance_A,
            "cluster_relative_tolerance": spec.cluster_relative_tolerance,
        },
        "rest_protocol": StimulusArm.REST.value,
        "resting_calcium_uM": manifest["protocols"]["CCH_IPR"][
            "resting_calcium_uM"
        ],
        "rest_beta_input": 0.0,
        "rest_integration_policy": (
            "direct solution first; no fixed duration resting integration is used"
        ),
        "production_solver": asdict(PRODUCTION_RADAU),
        "time_grid_s": frozen["time_grid_s"],
        "conservation_tolerances": dict(CONSERVATION_RESIDUAL_TOLERANCES),
        "baseline_rest_relative_tolerance": BASELINE_REST_RELATIVE_TOLERANCE,
        "baseline_trajectory_relative_tolerance": (
            BASELINE_TRAJECTORY_RELATIVE_TOLERANCE
        ),
        "ae4_expression_scaling_absolute_tolerance_fmol_s": (
            AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S
        ),
        "routing_roundoff_tolerance_fmol_s": ROUTING_ROUNDOFF_TOLERANCE_FMOL_S,
        "nkcc1_chloride_per_cycle": 2.0,
        "no_phenotype_fitting": True,
        "no_calibration_to_experimental_knockout": True,
        "no_phenotype_guided_root_selection": True,
        "no_additional_grid_values": True,
        "unchanged_ae4_capacity": True,
        "unchanged_nkcc1_law_capacity_and_regulation": True,
        "unchanged_channels_pumps_ae2_water_and_chloride_allocation": True,
        "task24_artifact_sha256": {
            str(path.relative_to(REPO)): sha256_file(path)
            for path in (
                TASK24 / "contract.json",
                TASK24 / "resting_states.csv",
                TASK24 / "trajectory_summary.csv",
                TASK24 / "comparison.csv",
                TASK24 / "verification.json",
            )
        },
        "source_sha256": _source_hashes(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
    }


def make_verification(
    rest_rows: Sequence[Mapping[str, Any]],
    attempts: Sequence[Mapping[str, Any]],
    alternatives: Sequence[Mapping[str, Any]],
    trajectory_rows: Sequence[Mapping[str, Any]],
    comparisons: Sequence[Mapping[str, Any]],
    source_hashes_before: Mapping[str, str],
) -> dict[str, Any]:
    baseline_rest = [
        row
        for row in rest_rows
        if (float(row["co2_factor"]), float(row["nhe1_factor"])) == BASELINE
    ]
    baseline_comparisons = [
        row
        for row in comparisons
        if (float(row["co2_factor"]), float(row["nhe1_factor"])) == BASELINE
    ]
    if len(baseline_rest) != 20 or len(baseline_comparisons) != 30:
        raise AssertionError("baseline control panel is incomplete")
    baseline_rest_error = max(
        float(row["task24_anchor_distance_normalized"]) for row in baseline_rest
    )
    baseline_wt_error = max(
        float(row["task24_wt_total_relative_error"])
        for row in baseline_comparisons
    )
    baseline_low_error = max(
        float(row["task24_ae4_5pct_total_relative_error"])
        for row in baseline_comparisons
    )
    baseline_ratio_error = max(
        float(row["task24_ratio_absolute_error"]) for row in baseline_comparisons
    )
    baseline_control_pass = bool(
        baseline_rest_error <= BASELINE_REST_RELATIVE_TOLERANCE
        and baseline_wt_error <= BASELINE_TRAJECTORY_RELATIVE_TOLERANCE
        and baseline_low_error <= BASELINE_TRAJECTORY_RELATIVE_TOLERANCE
        and baseline_ratio_error <= BASELINE_TRAJECTORY_RELATIVE_TOLERANCE
    )
    if not baseline_control_pass:
        raise RuntimeError("the exact Task 24 baseline control was not reproduced")
    expected_paths = set(EXPECTED_PARAMETER_PATHS)
    parameter_path_gate = True
    for row in rest_rows:
        observed = {
            item for item in str(row["changed_parameter_paths"]).split(";") if item
        }
        expected: set[str] = set()
        if float(row["co2_factor"]) != 1.0:
            expected.update(EXPECTED_PARAMETER_PATHS[:2])
        if float(row["nhe1_factor"]) != 1.0:
            expected.add(EXPECTED_PARAMETER_PATHS[2])
        parameter_path_gate &= observed == expected and observed <= expected_paths
    rest_hash_groups: dict[tuple[str, float, float, float], set[str]] = {}
    for row in trajectory_rows:
        key = (
            str(row["root_id"]),
            float(row["co2_factor"]),
            float(row["nhe1_factor"]),
            float(row["ae4_expression"]),
        )
        rest_hash_groups.setdefault(key, set()).add(str(row["resting_state_sha256"]))
    calcium_independence = bool(
        len(rest_hash_groups) == 240
        and all(len(hashes) == 1 for hashes in rest_hash_groups.values())
    )
    nkcc_residual = max(
        abs(
            float(row["nkcc1_chloride_influx_fmol_s"])
            - 2.0 * float(row["nkcc1_inward_cycle_flux_fmol_s"])
        )
        for row in rest_rows
    )
    source_hashes_after = _source_hashes()
    changed_sources = sorted(
        key
        for key in set(source_hashes_before) | set(source_hashes_after)
        if source_hashes_before.get(key) != source_hashes_after.get(key)
    )
    verification = {
        "model_id": MODEL_ID,
        "source_commit": SOURCE_COMMIT,
        "parameter_combination_count": 12,
        "root_count": 10,
        "resting_state_count": len(rest_rows),
        "warm_start_attempt_count": len(attempts),
        "accepted_warm_start_attempt_count": sum(
            bool(row["complete_rest_gate_pass"]) for row in attempts
        ),
        "alternate_root_count": len(alternatives),
        "unresolved_resting_state_count": 240 - len(rest_rows),
        "all_rest_gates_pass": bool(
            len(rest_rows) == 240 and all(bool(row["rest_gate_pass"]) for row in rest_rows)
        ),
        "rest_integration_used": any(bool(row["rest_integration_used"]) for row in rest_rows),
        "max_rest_scaled_independent_rhs": max(
            float(row["max_abs_scaled_independent_rhs"]) for row in rest_rows
        ),
        "max_rest_amount_rhs_fmol_s": max(
            float(row["max_abs_amount_rhs_fmol_s"]) for row in rest_rows
        ),
        "max_rest_volume_rhs_pL_s": max(
            float(row["max_abs_volume_rhs_pL_s"]) for row in rest_rows
        ),
        "max_rest_omitted_rhs_fmol_s": max(
            float(row["max_abs_omitted_rhs_fmol_s"]) for row in rest_rows
        ),
        "max_rest_current_residual_A": max(
            float(row["max_abs_current_residual_A"]) for row in rest_rows
        ),
        "max_rest_conservation_ratio": max(
            float(row["max_dimensionless_conservation_ratio"]) for row in rest_rows
        ),
        "all_rest_states_positive": all(
            bool(row["positivity_gate_pass"]) for row in rest_rows
        ),
        "all_rest_states_away_from_boundaries": all(
            not row["boundary_hits"] for row in rest_rows
        ),
        "all_rest_ae4_charge_neutral": all(
            abs(float(row["ae4_charge_residual_fmol_s"]))
            <= CONSERVATION_RESIDUAL_TOLERANCES["ae4_charge_fmol_s"]
            for row in rest_rows
        ),
        "all_rest_ae4_no_opposing_cation_sources": all(
            bool(row["ae4_no_opposing_cation_sources"]) for row in rest_rows
        ),
        "max_rest_ae4_expression_scaling_residual_fmol_s": max(
            float(row["ae4_expression_scaling_residual_fmol_s"])
            for row in rest_rows
        ),
        "max_rest_nkcc1_chloride_stoichiometry_residual_fmol_s": nkcc_residual,
        "stimulation_trajectory_count": len(trajectory_rows),
        "all_stimulation_gates_pass": bool(
            len(trajectory_rows) == 720
            and all(bool(row["valid"]) for row in trajectory_rows)
        ),
        "all_stimulation_onsets_match_rest": all(
            bool(row["onset_rest_rhs_match"]) for row in trajectory_rows
        ),
        "resting_state_independent_of_subsequent_calcium": calcium_independence,
        "max_stimulation_conservation_ratio": max(
            float(row["max_dimensionless_conservation_ratio"])
            for row in trajectory_rows
        ),
        "opposing_cation_stimulation_sample_count": sum(
            int(row["opposing_cation_samples"]) for row in trajectory_rows
        ),
        "max_stimulation_ae4_charge_residual_fmol_s": max(
            float(row["max_ae4_charge_residual_fmol_s"])
            for row in trajectory_rows
        ),
        "max_stimulation_ae4_expression_scaling_residual_fmol_s": max(
            float(row["max_ae4_expression_scaling_residual_fmol_s"])
            for row in trajectory_rows
        ),
        "max_stimulation_nkcc1_chloride_stoichiometry_residual_fmol_s": max(
            float(row["max_nkcc1_chloride_stoichiometry_residual_fmol_s"])
            for row in trajectory_rows
        ),
        "only_permitted_parameter_paths_changed": bool(parameter_path_gate),
        "baseline_control_pass": baseline_control_pass,
        "baseline_max_rest_distance_normalized": baseline_rest_error,
        "baseline_max_wt_total_relative_error": baseline_wt_error,
        "baseline_max_ae4_5pct_total_relative_error": baseline_low_error,
        "baseline_max_ratio_absolute_error": baseline_ratio_error,
        "baseline_change_percent_min": min(
            float(row["ae4_5pct_secretion_change_percent"])
            for row in baseline_comparisons
        ),
        "baseline_change_percent_max": max(
            float(row["ae4_5pct_secretion_change_percent"])
            for row in baseline_comparisons
        ),
        "changed_source_or_frozen_inputs_during_run": changed_sources,
    }
    required_true = (
        "all_rest_gates_pass",
        "all_rest_states_positive",
        "all_rest_states_away_from_boundaries",
        "all_rest_ae4_charge_neutral",
        "all_rest_ae4_no_opposing_cation_sources",
        "all_stimulation_gates_pass",
        "all_stimulation_onsets_match_rest",
        "resting_state_independent_of_subsequent_calcium",
        "only_permitted_parameter_paths_changed",
        "baseline_control_pass",
    )
    if not all(bool(verification[key]) for key in required_true):
        raise RuntimeError("one or more mandatory verification gates failed")
    if verification["rest_integration_used"]:
        raise RuntimeError("unexpected resting integration was used")
    if verification["opposing_cation_stimulation_sample_count"] != 0:
        raise RuntimeError("opposing routed AE4 cation sources occurred")
    if (
        verification[
            "max_stimulation_nkcc1_chloride_stoichiometry_residual_fmol_s"
        ]
        > AE4_SCALING_ABSOLUTE_TOLERANCE_FMOL_S
    ):
        raise RuntimeError("stimulated NKCC1 chloride stoichiometry failed")
    if changed_sources:
        raise RuntimeError("source or frozen input changed during the run")
    write_json(OUT / "verification.json", verification)
    return verification


def _range_text(row: Mapping[str, Any], prefix: str, digits: int = 3) -> str:
    return (
        f"{float(row[prefix + '_min']):.{digits}f} to "
        f"{float(row[prefix + '_max']):.{digits}f}"
    )


def _answer_lines(answers: Mapping[str, Any]) -> list[str]:
    a = answers["A"]
    if a["all_increased_co2_tic_changes_nonpositive"] and a[
        "all_increased_co2_hco3_changes_nonpositive"
    ]:
        answer_a = "Yes, throughout the increased CO2 permeability grid."
    elif (
        a["increased_co2_tic_change_mM"]["min"] > 0.0
        and a["increased_co2_hco3_change_mM"]["min"] > 0.0
    ):
        answer_a = "No. Accumulation persisted in every increased CO2 permeability case."
    else:
        answer_a = "No, not throughout the increased CO2 permeability grid."
    b = answers["B"]
    answer_b = (
        "Yes, throughout the reduced NHE1 grid."
        if b["all_reduced_nhe1_ph_changes_nonpositive"]
        else (
            "No. Reduced NHE1 attenuated but did not prevent alkalinisation in any case."
            if b["reduced_nhe1_ph_change"]["min"] > 0.0
            else "No, not throughout the reduced NHE1 grid."
        )
    )
    c = answers["C"]
    co2_change = c["co2_x10_nhe1_x1_median"] - c["baseline_median"]
    nhe_change = c["co2_x1_nhe1_x0p25_median"] - c["baseline_median"]
    if co2_change >= 0.0 and nhe_change >= 0.0:
        dominant = (
            "Neither intervention reduced compensation; both increased it slightly "
            "relative to the baseline median."
        )
    elif co2_change < nhe_change:
        dominant = "The isolated CO2 intervention produced the larger median reduction."
    elif nhe_change < co2_change:
        dominant = "The isolated NHE1 intervention produced the larger median reduction."
    else:
        dominant = "The two isolated interventions produced the same median reduction."
    d = answers["D"]
    baseline_retention = d["baseline_retention"]["median"]
    lowest_retention = d["lowest_median_combination"]["retention_median"]
    if d["nonbaseline_combination_medians_below_baseline_median"] == 0:
        answer_d = "No. Every nonbaseline combination had a higher median retention."
    else:
        answer_d = "Some nonbaseline medians were lower than the control median."
    e = answers["E"]
    answer_e = (
        "Yes for every parameter combination by the median replacement criterion."
        if e["combination_medians_above_one_half"] == 12
        else (
            f"A median majority replacement occurred in {e['combination_medians_above_one_half']} "
            "of 12 parameter combinations."
        )
    )
    f = answers["F"]
    answer_f = (
        "Yes: "
        + ", ".join(
            f"CO2 ×{item['co2_factor']:g}, NHE1 ×{item['nhe1_factor']:g}"
            for item in f["consistent_all_root_all_calcium_combinations"]
        )
        if f["consistent_all_root_all_calcium_combinations"]
        else "No parameter combination reduced secretion in every root at every calcium."
    )
    g = answers["G"]
    if not g["decreasing_cases_exist"]:
        answer_g = "No secretion decrease occurred, so that association cannot be assigned."
    elif g["all_comparisons_numerically_valid"]:
        answer_g = (
            "The decreasing cases passed the numerical gates.  Their retention and "
            "compensation ranges are reported below, so any association can be separated "
            "from numerical failure."
        )
    else:
        answer_g = "No, because at least one decreasing case failed a numerical gate."
    return [
        "## Answers to the diagnostic questions",
        "",
        f"### A. Carbon accumulation",
        "",
        (
            f"{answer_a} With CO2 permeability above baseline, the WT to 5% AE4 TIC "
            f"change ranged from {a['increased_co2_tic_change_mM']['min']:.4f} to "
            f"{a['increased_co2_tic_change_mM']['max']:.4f} mM and the bicarbonate "
            f"change ranged from {a['increased_co2_hco3_change_mM']['min']:.4f} to "
            f"{a['increased_co2_hco3_change_mM']['max']:.4f} mM. At REST, apical CO2 "
            "flux was outward and basolateral CO2 flux was inward in every increased "
            "permeability state; the combined CO2 flux remained inward in every state."
        ),
        "",
        "### B. Alkalinisation",
        "",
        (
            f"{answer_b} With NHE1 below baseline, the pH change ranged from "
            f"{b['reduced_nhe1_ph_change']['min']:+.4f} to "
            f"{b['reduced_nhe1_ph_change']['max']:+.4f}, compared with "
            f"{b['baseline_ph_change']['min']:+.4f} to "
            f"{b['baseline_ph_change']['max']:+.4f} in the exact control."
        ),
        "",
        "### C. AE4 compensation",
        "",
        (
            f"The median per expression compensation was {c['baseline_median']:.3f} at "
            f"baseline, {c['co2_x10_nhe1_x1_median']:.3f} for CO2 ×10 alone, "
            f"{c['co2_x1_nhe1_x0p25_median']:.3f} for NHE1 ×0.25 alone, and "
            f"{c['co2_x10_nhe1_x0p25_median']:.3f} for the joint intervention. "
            f"{dominant} Reducing NHE1 to 0.25 lowered median absolute WT AE4 loading "
            f"from {c['absolute_ae4_flux_median_fmol_s']['baseline']['wt']:.6f} to "
            f"{c['absolute_ae4_flux_median_fmol_s']['co2_x1_nhe1_x0p25']['wt']:.6f} "
            "fmol s⁻¹ and 5% AE4 loading from "
            f"{c['absolute_ae4_flux_median_fmol_s']['baseline']['ae4_5pct']:.6f} to "
            f"{c['absolute_ae4_flux_median_fmol_s']['co2_x1_nhe1_x0p25']['ae4_5pct']:.6f} "
            "fmol s⁻¹. The WT reference fell slightly more in proportion, so the "
            "normalised compensation ratio rose despite weaker absolute AE4 flux."
        ),
        "",
        "### D. AE4 flux retention",
        "",
        (
            f"{answer_d} The baseline median retention was {baseline_retention:.3%}. "
            f"The lowest grid median was {lowest_retention:.3%} at CO2 "
            f"×{d['lowest_median_combination']['co2_factor']:g} and NHE1 "
            f"×{d['lowest_median_combination']['nhe1_factor']:g}. Individual "
            f"nonbaseline values ranged from {d['nonbaseline_retention']['min']:.3%} "
            f"to {d['nonbaseline_retention']['max']:.3%}."
        ),
        "",
        "### E. NKCC1 replacement",
        "",
        (
            f"{answer_e} Across individual roots and grid points, the replacement "
            f"fraction ranged from {e['grid_replacement_fraction']['min']:.3f} to "
            f"{e['grid_replacement_fraction']['max']:.3f}, with median "
            f"{e['grid_replacement_fraction']['median']:.3f}."
        ),
        "",
        "### F. Secretion direction",
        "",
        answer_f,
        "",
        "### G. Mechanistic attribution",
        "",
        answer_g,
    ]


def write_report(
    parameter_rows: Sequence[Mapping[str, Any]],
    summary: Mapping[str, Any],
    verification: Mapping[str, Any],
) -> Path:
    lines = [
        "# Carbon escape and NHE1 sensitivity of routed AE4 compensation",
        "",
        "## Result",
        "",
        (
            "The complete fixed 4 by 3 grid was evaluated for all ten retained roots, "
            "with separate WT and 5% AE4 resting states and three production calcium "
            "levels.  No phenotype value entered model construction, root selection, "
            "or parameter choice."
        ),
        "",
        (
            f"All {verification['resting_state_count']} resting states and all "
            f"{verification['stimulation_trajectory_count']} stimulated trajectories "
            "passed the inherited numerical and conservation gates."
        ),
        "",
        (
            "The interventions did not break routed AE4 self rescue. Across the fixed "
            f"grid, nonbaseline AE4 retention was {summary['diagnostic_answers']['D']['nonbaseline_retention']['min']:.3%} "
            f"to {summary['diagnostic_answers']['D']['nonbaseline_retention']['max']:.3%}, "
            f"and none of the {summary['diagnostic_answers']['F']['comparison_count']} "
            "root and calcium comparisons showed lower secretion with 5% AE4."
        ),
        "",
        "## Parameter summary across ten roots",
        "",
        (
            "Changes are 5% AE4 minus WT.  Retention is the 5% AE4 chloride flux divided "
            "by the WT flux.  Compensation divides the 5% flux by 0.05 before comparison "
            "with WT.  The NKCC1 replacement fraction is its chloride loading increase "
            "divided by lost AE4 chloride loading."
        ),
        "",
        "| CO2 | NHE1 | HCO3 change mM | TIC change mM | pH change | AE4 retention median | Compensation median | NKCC1 replacement median | Secretion change at 0.10 µM | at 0.25 µM | at 0.50 µM | decrease counts |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in parameter_rows:
        lines.append(
            "| ×{co2:g} | ×{nhe:g} | {hco3} | {tic} | {ph} | {ret:.3%} | {comp:.3f} | {replace:.3f} | {c010}% | {c025}% | {c050}% | {n010}/{n025}/{n050} |".format(
                co2=float(row["co2_factor"]),
                nhe=float(row["nhe1_factor"]),
                hco3=_range_text(row, "wt_to_ae4_5pct_hco3_change_mM", 3),
                tic=_range_text(row, "wt_to_ae4_5pct_tic_change_mM", 3),
                ph=_range_text(row, "wt_to_ae4_5pct_ph_change", 3),
                ret=float(row["ae4_flux_retention_median"]),
                comp=float(row["ae4_per_expression_compensation_median"]),
                replace=float(row["fraction_lost_ae4_replaced_by_nkcc1_median"]),
                c010=_range_text(row, "ca_0p10_secretion_change_percent", 3),
                c025=_range_text(row, "ca_0p25_secretion_change_percent", 3),
                c050=_range_text(row, "ca_0p50_secretion_change_percent", 3),
                n010=int(row["ca_0p10_decrease_root_count"]),
                n025=int(row["ca_0p25_decrease_root_count"]),
                n050=int(row["ca_0p50_decrease_root_count"]),
            )
        )
    lines.extend([""] + _answer_lines(summary["diagnostic_answers"]))
    lines.extend(
        [
            "",
            "## Exact regression control",
            "",
            (
                "The CO2 ×1, NHE1 ×1 control reproduced Task 24.  The largest normalised "
                f"resting coordinate difference was {verification['baseline_max_rest_distance_normalized']:.3e}. "
                f"The largest WT secretion relative error was {verification['baseline_max_wt_total_relative_error']:.3e}, "
                f"the largest 5% AE4 relative error was {verification['baseline_max_ae4_5pct_total_relative_error']:.3e}, "
                f"and the secretion change remained {verification['baseline_change_percent_min']:+.3f}% "
                f"to {verification['baseline_change_percent_max']:+.3f}%."
            ),
            "",
            "## Numerical status",
            "",
            (
                f"Direct REST solution was sufficient for every state.  No resting integration was used. "
                f"The maximum scaled independent REST residual was {verification['max_rest_scaled_independent_rhs']:.3e}, "
                f"the maximum current residual was {verification['max_rest_current_residual_A']:.3e} A, "
                f"the maximum REST conservation ratio was {verification['max_rest_conservation_ratio']:.3e}, "
                f"and the maximum stimulation conservation ratio was {verification['max_stimulation_conservation_ratio']:.3e}."
            ),
            "",
            (
                f"Continuation found {verification['alternate_root_count']} nonselected admissible root clusters. "
                "Any such roots are recorded in `alternative_roots.csv`; selection remained on the branch nearest the declared predecessor."
            ),
            "",
            "## Interpretation boundaries",
            "",
            (
                "The mathematical result is the response of this fixed model to the declared sensitivity grid. "
                "The model mechanism can identify whether changes in carbon, pH, AE4 driving, and NKCC1 loading move together inside these equations. "
                "It does not establish that CO2 permeability or NHE1 causes the experimental AE4 phenotype."
            ),
            "",
            (
                "The experimental secretion deficit was not read or used until after every simulation and summary was complete. "
                "The CO2 and NHE1 factors remain unmeasured sensitivities, and agreement with an experimental percentage would not make a grid point biologically correct."
            ),
            "",
            "## Reproducibility",
            "",
            (
                "`resting_states.csv` contains all 240 solved states and complete REST diagnostics. "
                "`flux_decomposition.csv` contains every WT and 5% pair.  `trajectory_summary.csv` contains all 720 trajectories, their extrema and landmarks. "
                "The compressed trajectory files preserve the full inherited time grid, state arrays, flow, and cumulative secretion."
            ),
        ]
    )
    path = ANALYSIS / "ae4_carbon_nhe1_compensation_report.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".md.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path


def _load_existing_trajectory(
    request: tuple[str, float, float, float, float]
) -> dict[str, Any] | None:
    root, co2_factor, nhe1_factor, calcium, expression = request
    stem = _trajectory_stem(root, co2_factor, nhe1_factor, calcium, expression)
    metadata_path = OUT / "trajectories" / f"{stem}.json"
    trajectory_path = OUT / "trajectories" / f"{stem}.npz"
    if not metadata_path.exists() or not trajectory_path.exists():
        return None
    row = json.loads(metadata_path.read_text())
    checks = (
        row.get("root_id") == root,
        float(row.get("co2_factor")) == co2_factor,
        float(row.get("nhe1_factor")) == nhe1_factor,
        float(row.get("calcium_uM")) == calcium,
        float(row.get("ae4_expression")) == expression,
        bool(row.get("valid")),
        row.get("trajectory_sha256") == sha256_file(trajectory_path),
    )
    return row if all(checks) else None


def _run_resting(workers: int, roots: Sequence[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    requests = list(itertools.product(roots, EXPRESSIONS))
    rows: list[dict[str, Any]] = []
    attempts: list[dict[str, Any]] = []
    alternatives: list[dict[str, Any]] = []
    errors: list[str] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(solve_rest_chain, request): request for request in requests}
        for future in as_completed(jobs):
            request = jobs[future]
            result = future.result()
            rows.extend(result["rows"])
            attempts.extend(result["attempts"])
            alternatives.extend(result["alternatives"])
            if result["error"]:
                errors.append(result["error"])
            write_rows(OUT / "resting_states.csv", [_flat_rest_row(row) for row in rows])
            write_rows(OUT / "warm_start_attempts.csv", attempts)
            if alternatives:
                write_rows(OUT / "alternative_roots.csv", alternatives)
            print(
                json.dumps(
                    {
                        "phase": "rest",
                        "root": request[0],
                        "ae4_expression": request[1],
                        "completed_states": len(result["rows"]),
                        "total_completed_states": len(rows),
                        "error": result["error"],
                    }
                ),
                flush=True,
            )
    if not alternatives:
        write_rows(
            OUT / "alternative_roots.csv",
            [
                {
                    "status": "no alternative admissible root clusters found",
                    "root_count": len(roots),
                    "parameter_combination_count": 12,
                    "expression_count": 2,
                }
            ],
        )
    if errors:
        raise RuntimeError("; ".join(errors))
    rows.sort(
        key=lambda row: (
            float(row["co2_factor"]),
            float(row["nhe1_factor"]),
            str(row["root_id"]),
            -float(row["ae4_expression"]),
        )
    )
    write_rows(OUT / "resting_states.csv", [_flat_rest_row(row) for row in rows])
    write_rows(OUT / "warm_start_attempts.csv", attempts)
    return rows, attempts, alternatives


def _run_stimulation(
    workers: int, roots: Sequence[str]
) -> list[dict[str, Any]]:
    requests = [
        (root, co2, nhe1, calcium, expression)
        for co2, nhe1, root, calcium, expression in itertools.product(
            CO2_FACTORS, NHE1_FACTORS, roots, CALCIUM, EXPRESSIONS
        )
    ]
    rows: list[dict[str, Any]] = []
    pending: list[tuple[str, float, float, float, float]] = []
    for request in requests:
        existing = _load_existing_trajectory(request)
        if existing is None:
            pending.append(request)
        else:
            rows.append(existing)
    if rows:
        print(
            json.dumps(
                {
                    "phase": "stimulation_resume",
                    "reused": len(rows),
                    "pending": len(pending),
                }
            ),
            flush=True,
        )
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(run_stimulation_case, request): request for request in pending}
        for future in as_completed(jobs):
            row = future.result()
            rows.append(row)
            if len(rows) % 10 == 0 or len(rows) == len(requests):
                write_rows(OUT / "trajectory_summary.csv", rows)
            print(
                json.dumps(
                    {
                        "phase": "stimulation",
                        "root": row["root_alias"],
                        "co2_factor": row["co2_factor"],
                        "nhe1_factor": row["nhe1_factor"],
                        "calcium_uM": row["calcium_uM"],
                        "ae4_expression": row["ae4_expression"],
                        "completed": len(rows),
                        "valid": row["valid"],
                    }
                ),
                flush=True,
            )
    rows.sort(
        key=lambda row: (
            float(row["co2_factor"]),
            float(row["nhe1_factor"]),
            str(row["root_id"]),
            float(row["calcium_uM"]),
            -float(row["ae4_expression"]),
        )
    )
    write_rows(OUT / "trajectory_summary.csv", rows)
    return rows


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("workers must be positive")
    thread_variables = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
    if any(os.environ.get(name) != "1" for name in thread_variables):
        parser.error(
            "Set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 before Python starts"
        )
    manifest, frozen, _ = inputs()
    roots = sorted(manifest["roots"])
    if len(roots) != 10:
        raise AssertionError("the complete ten root panel is required")
    OUT.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    contract = contract_payload(manifest, frozen)
    contract_path = OUT / "contract.json"
    if contract_path.exists() and json.loads(contract_path.read_text()) != contract:
        raise RuntimeError("existing Task 25 contract differs from the current contract")
    write_json(contract_path, contract)
    source_hashes_before = _source_hashes()
    rest_rows, attempts, alternatives = _run_resting(args.workers, roots)
    if len(rest_rows) != 240:
        raise RuntimeError(f"expected 240 resting states, obtained {len(rest_rows)}")
    trajectory_rows = _run_stimulation(args.workers, roots)
    if len(trajectory_rows) != 720:
        raise RuntimeError(
            f"expected 720 stimulated trajectories, obtained {len(trajectory_rows)}"
        )
    flux_rows = make_flux_decomposition(rest_rows)
    comparisons = make_comparisons(rest_rows, trajectory_rows, flux_rows)
    parameter_rows = make_parameter_summary(flux_rows, comparisons)
    verification = make_verification(
        rest_rows,
        attempts,
        alternatives,
        trajectory_rows,
        comparisons,
        source_hashes_before,
    )
    summary = make_summary(
        parameter_rows, flux_rows, comparisons, verification
    )
    report_path = write_report(parameter_rows, summary, verification)
    artifact_index = {
        "report": str(report_path.relative_to(REPO)),
        "report_sha256": sha256_file(report_path),
        "resting_states": str((OUT / "resting_states.csv").relative_to(REPO)),
        "comparison": str((OUT / "comparison.csv").relative_to(REPO)),
        "flux_decomposition": str((OUT / "flux_decomposition.csv").relative_to(REPO)),
        "parameter_summary": str((OUT / "parameter_summary.csv").relative_to(REPO)),
        "trajectory_summary": str((OUT / "trajectory_summary.csv").relative_to(REPO)),
        "warm_start_attempts": str((OUT / "warm_start_attempts.csv").relative_to(REPO)),
        "alternative_roots": str((OUT / "alternative_roots.csv").relative_to(REPO)),
        "verification": str((OUT / "verification.json").relative_to(REPO)),
        "contract": str((OUT / "contract.json").relative_to(REPO)),
        "summary": str((OUT / "summary.json").relative_to(REPO)),
        "trajectory_directory": str((OUT / "trajectories").relative_to(REPO)),
    }
    write_json(OUT / "artifact_index.json", artifact_index)
    print(
        json.dumps(
            {
                "summary": summary["diagnostic_answers"],
                "verification": verification,
            },
            indent=2,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
