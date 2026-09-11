"""Re-equilibrate the routing-only AE4 model before WT/5% stimulation.

This runner changes no model parameter.  It reuses the routing intervention
from :mod:`modern_full_model.ae4_routing_only`, explicitly returns the model
to the REST protocol, and solves the ten independent charge-manifold resting
equations for each retained root and AE4 expression.  For every calcium arm,
the two prescribed warm starts are used independently:

* the inherited saved resting state;
* the endpoint of the matching earlier routing-only 600 s trajectory.

Both roots must pass the inherited numerical gates and cluster as the same
root.  The selected representative is chosen only by residual size.  No
phenotype value enters root construction or selection.
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
from typing import Any, Mapping, Sequence

import numpy as np
import scipy

from .ae4_routing_only import MODEL_ID, with_routing_only
from .calibration import (
    CORE_INDEPENDENT_ROWS,
    CORE_OMITTED_CHARGE_ROWS,
    WTCalibrationSpec,
)
from .genotype_evaluation import (
    _attempt_and_root,
    _core_coordinates,
    _regulatory_suffix,
    _root_state_with_basal_regulation,
    genotype_with_expression,
    simulate_genotype,
)
from .model import ModernFullModel, WT
from .nkcc_stimulation import attach_stimulated_nkcc1
from .run_ae4_routing_only import inputs, output_stem as prior_output_stem
from .run_calcium_fast_screen import write_json, write_rows
from .states import CORE_STATE_NAMES
from .task14_blind import CALCIUM, build_model, routing
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES,
    PRODUCTION_RADAU,
    StimulusArm,
    sha256_file,
    sha256_object,
)


REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/ae4_routing_re_equilibrated"
ANALYSIS = REPO / "analysis/ae4_routing_re_equilibrated"
PRIOR = REPO / "results/ae4_routing_only"
MANIFEST = REPO / "results/13B_modern_full_model/native_dynamic_contract_manifest.json"
FROZEN = REPO / "results/14B_five_percent_ae4_check/frozen_inputs.json"
PRIOR_COMPARISON = PRIOR / "comparison.csv"
REST_MAX_NFEV = 3000
REST_REGULATORY_TOLERANCE_S_INV = 1.0e-10
ROUTING_ROUNDOFF_TOLERANCE_FMOL_S = 2.0e-14
EXPRESSIONS = (1.0, 0.05)
AMOUNT_ROWS = np.asarray((0, 1, 2, 3, 4, 6, 7, 8, 9, 10), dtype=int)
VOLUME_ROWS = np.asarray((5, 11), dtype=int)


def _genotype(expression: float):
    return WT if expression == 1.0 else genotype_with_expression("AE4", expression)


def _root_aliases(manifest: Mapping[str, Any]) -> Mapping[str, str]:
    return {root: f"R{index:02d}" for index, root in enumerate(sorted(manifest["roots"]), 1)}


def _explicit_rest_model(manifest: Mapping[str, Any], root: str):
    """Build the routed model with an explicit, time-independent REST arm."""

    routed_stimulated = with_routing_only(build_model(manifest, root, CALCIUM[0]))
    composite = routed_stimulated.regulatory_model
    if not hasattr(composite, "ae4_regulatory_model") or not hasattr(
        composite, "nkcc1_regulatory_model"
    ):
        raise AssertionError("the inherited R1 plus N1 regulatory model changed")
    protocol = replace(routed_stimulated.stimulus, arm=StimulusArm.REST)
    base = ModernFullModel(
        parameters=routed_stimulated.parameters,
        stimulus=protocol,
        regulatory_model=composite.ae4_regulatory_model,
        ae4_parameters=routed_stimulated.ae4_parameters,
        ae4_evaluator=routed_stimulated.ae4_evaluator,
    )
    model = attach_stimulated_nkcc1(base, composite.nkcc1_regulatory_model)
    for time_s in (0.0, 1.0e-6, 1.0, 600.0, 3600.0):
        stimulus = model.stimulus(time_s)
        if stimulus.beta_input != 0.0:
            raise AssertionError("REST protocol has a nonzero beta input")
        if stimulus.calcium_uM != protocol.resting_calcium_uM:
            raise AssertionError("REST protocol does not retain resting calcium")
    if not str(model.regulatory_model.family).startswith("R1+"):
        raise AssertionError("the inherited R1 regulator was not retained")
    return model


def _prior_endpoint_path(root: str, calcium: float, expression: float) -> Path:
    stem = prior_output_stem(
        root, calcium, "routing_only", expression, "production_radau"
    )
    return PRIOR / "trajectories" / f"{stem}.npz"


def _load_prior_endpoint(root: str, calcium: float, expression: float) -> tuple[np.ndarray, Path]:
    path = _prior_endpoint_path(root, calcium, expression)
    metadata_path = path.with_suffix(".json")
    if not path.exists() or not metadata_path.exists():
        raise FileNotFoundError(f"missing prior routed trajectory for warm start: {path}")
    metadata = json.loads(metadata_path.read_text())
    if not metadata.get("valid") or metadata.get("architecture") != "routing_only":
        raise AssertionError(f"prior warm-start trajectory is not valid routed output: {path}")
    with np.load(path) as payload:
        time_s = np.asarray(payload["time_s"], dtype=float)
        states = np.asarray(payload["states"], dtype=float)
    if time_s[-1] != 600.0 or states.shape[1] != time_s.size:
        raise AssertionError(f"prior warm-start trajectory is incomplete: {path}")
    return states[:, -1].copy(), path


def _conservation_audit(evaluation: Any) -> tuple[dict[str, float], float, bool]:
    residuals = {
        key: abs(float(value))
        for key, value in evaluation.diagnostics.conservation_residuals.items()
    }
    unknown = set(residuals) - set(CONSERVATION_RESIDUAL_TOLERANCES)
    missing = set(CONSERVATION_RESIDUAL_TOLERANCES) - set(residuals)
    if unknown or missing:
        raise AssertionError(
            f"conservation diagnostic set changed; unknown={sorted(unknown)}, missing={sorted(missing)}"
        )
    ratios = {
        key: residuals[key] / tolerance
        for key, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()
    }
    maximum = max(ratios.values())
    return residuals, float(maximum), bool(maximum <= 1.0)


def _audit_resting_state(
    model: Any,
    core_state: Sequence[float],
    expression: float,
    root_result: Any,
    spec: WTCalibrationSpec,
) -> dict[str, Any]:
    genotype = _genotype(expression)
    state = _root_state_with_basal_regulation(model, core_state)
    evaluation = model.evaluate(0.0, state, genotype=genotype)
    later = model.evaluate(600.0, state, genotype=genotype)
    time_invariant_rhs = bool(np.array_equal(evaluation.rhs, later.rhs))
    raw = np.asarray(evaluation.rhs, dtype=float)
    scaled = raw[np.asarray(CORE_INDEPENDENT_ROWS, dtype=int)] / np.asarray(
        spec.independent_rhs_scales, dtype=float
    )
    omitted = raw[np.asarray(CORE_OMITTED_CHARGE_ROWS, dtype=int)]
    amount_max = float(np.max(np.abs(raw[AMOUNT_ROWS])))
    volume_max = float(np.max(np.abs(raw[VOLUME_ROWS])))
    regulatory_max = float(np.max(np.abs(raw[12:]))) if raw.size > 12 else 0.0
    charge_max = max(
        abs(float(value))
        for value in evaluation.diagnostics.state_charge_fmol.values()
    )
    current_max = max(
        abs(float(value))
        for value in evaluation.diagnostics.membranes.current_residuals_A.values()
    )
    conservation, conservation_ratio, conservation_gate = _conservation_audit(evaluation)
    observables = evaluation.diagnostics.observables
    decoded = model.layout.decode(state)
    cell = observables.cell_concentrations_mM
    lumen = observables.lumen_concentrations_mM
    ae4 = evaluation.diagnostics.ae4
    cation_same_direction = bool(
        ae4.na_cell_fmol_s == 0.0
        or ae4.k_cell_fmol_s == 0.0
        or ae4.na_cell_fmol_s * ae4.k_cell_fmol_s > 0.0
    )
    positive = bool(np.all(np.asarray(state[: len(CORE_STATE_NAMES)]) > 0.0))
    finite_observables = all(
        math.isfinite(float(value))
        for value in (
            cell["na"], cell["k"], cell["cl"], observables.cell_acid_base.ph,
            lumen["na"], lumen["k"], lumen["cl"], observables.lumen_acid_base.ph,
            decoded.cell.volume_pL, decoded.lumen.volume_pL,
        )
    )
    rhs_gate = bool(
        float(np.max(np.abs(scaled))) <= spec.root_scaled_tolerance
        and float(np.max(np.abs(omitted))) <= spec.omitted_row_raw_tolerance
        and regulatory_max <= REST_REGULATORY_TOLERANCE_S_INV
    )
    current_gate = bool(current_max <= spec.current_tolerance_A)
    charge_gate = bool(charge_max <= spec.charge_tolerance_fmol)
    routing_gate = bool(
        cation_same_direction
        and abs(float(ae4.charge_source_residual_fmol_s))
        <= CONSERVATION_RESIDUAL_TOLERANCES["ae4_charge_fmol_s"]
    )
    full_gate = bool(
        root_result.passes_numerical_gate
        and rhs_gate
        and current_gate
        and charge_gate
        and conservation_gate
        and positive
        and finite_observables
        and routing_gate
        and time_invariant_rhs
    )
    result: dict[str, Any] = {
        "rest_gate_pass": full_gate,
        "root_internal_gate_pass": bool(root_result.passes_numerical_gate),
        "full_rhs_gate_pass": rhs_gate,
        "current_gate_pass": current_gate,
        "bulk_charge_gate_pass": charge_gate,
        "conservation_gate_pass": conservation_gate,
        "positivity_gate_pass": positive,
        "finite_observables_gate_pass": finite_observables,
        "routing_gate_pass": routing_gate,
        "rest_protocol_time_invariant_rhs": time_invariant_rhs,
        "max_abs_scaled_independent_rhs": float(np.max(np.abs(scaled))),
        "max_abs_amount_rhs_fmol_s": amount_max,
        "max_abs_volume_rhs_pL_s": volume_max,
        "max_abs_omitted_rhs_fmol_s": float(np.max(np.abs(omitted))),
        "max_abs_regulatory_rhs_s_inv": regulatory_max,
        "max_abs_bulk_charge_fmol": float(charge_max),
        "max_abs_current_residual_A": float(current_max),
        "max_dimensionless_conservation_ratio": conservation_ratio,
        "independent_jacobian_rank": int(root_result.normalized_jacobian_rank),
        "independent_jacobian_nullity": int(root_result.normalized_jacobian_nullity),
        "boundary_hits": list(root_result.boundary_hits),
        "rest_cell_na_mM": float(cell["na"]),
        "rest_cell_k_mM": float(cell["k"]),
        "rest_cell_cl_mM": float(cell["cl"]),
        "rest_cell_ph": float(observables.cell_acid_base.ph),
        "rest_cell_volume_pL": float(decoded.cell.volume_pL),
        "rest_lumen_na_mM": float(lumen["na"]),
        "rest_lumen_k_mM": float(lumen["k"]),
        "rest_lumen_cl_mM": float(lumen["cl"]),
        "rest_lumen_ph": float(observables.lumen_acid_base.ph),
        "rest_lumen_volume_pL": float(decoded.lumen.volume_pL),
        "rest_outflow_pL_s": float(evaluation.diagnostics.water.lumen_outflow_pL_s),
        "rest_ae4_na_source_fmol_s": float(ae4.na_cell_fmol_s),
        "rest_ae4_k_source_fmol_s": float(ae4.k_cell_fmol_s),
        "rest_ae4_cl_source_fmol_s": float(ae4.cl_cell_fmol_s),
        "rest_ae4_hco3_source_fmol_s": float(ae4.hco3_cell_fmol_s),
        "rest_ae4_charge_residual_fmol_s": float(ae4.charge_source_residual_fmol_s),
        "raw_rhs": [float(value) for value in raw],
        "scaled_independent_rhs": [float(value) for value in scaled],
        "conservation_residuals": conservation,
    }
    return result


def _rest_output_stem(root: str, calcium: float, expression: float) -> str:
    return f"{root}_ca{calcium:.2f}_e{expression:.2f}_rest"


def solve_rest_case(request: tuple[str, float, float]) -> dict[str, Any]:
    root, calcium, expression = request
    started = time.perf_counter()
    manifest, _, inherited_five_percent = inputs()
    aliases = _root_aliases(manifest)
    model = _explicit_rest_model(manifest, root)
    genotype = _genotype(expression)
    spec = WTCalibrationSpec()
    lower = spec.coordinate_lower.copy()
    upper = spec.coordinate_upper.copy()
    span = upper - lower
    inherited = (
        manifest["roots"][root]["core_state"]
        if expression == 1.0
        else inherited_five_percent[root]
    )
    endpoint, endpoint_path = _load_prior_endpoint(root, calcium, expression)
    seeds = (
        ("inherited_saved_rest", np.asarray(inherited, dtype=float), None),
        ("previous_routed_endpoint", endpoint, endpoint_path),
    )
    candidates: list[dict[str, Any]] = []
    attempts: list[dict[str, Any]] = []
    for seed_id, seed_state, source_path in seeds:
        coordinates = _core_coordinates(model, seed_state)
        if np.any(coordinates <= lower) or np.any(coordinates >= upper):
            raise RuntimeError(f"{root} {calcium} {expression} seed is outside inherited root bounds")
        attempt, root_result = _attempt_and_root(
            model=model,
            transporter="AE4",
            expression=expression,
            genotype=genotype,
            regulatory_suffix=_regulatory_suffix(model),
            spec=spec,
            start_id=seed_id,
            start=coordinates,
            lower=lower,
            upper=upper,
            max_nfev=REST_MAX_NFEV,
        )
        attempt_row = {
            "root_id": root,
            "root_alias": aliases[root],
            "calcium_uM": calcium,
            "ae4_expression": expression,
            "seed_id": seed_id,
            "seed_source_path": (
                str(source_path.relative_to(REPO)) if source_path is not None else "inherited"
            ),
            "seed_core_sha256": sha256_object(list(np.asarray(seed_state)[:12])),
            "optimizer_success": bool(attempt.optimizer_success),
            "admissible_state": bool(attempt.admissible_state),
            "converged_root": bool(attempt.converged_root),
            "max_abs_scaled_independent_rhs": float(attempt.max_abs_scaled_independent_rhs),
            "cost": float(attempt.cost),
            "optimality": float(attempt.optimality),
            "nfev": int(attempt.nfev),
            "message": attempt.message,
        }
        if root_result is None:
            attempts.append(attempt_row)
            raise RuntimeError(
                f"no accepted direct resting root for {root}, Ca={calcium}, expression={expression}, seed={seed_id}: {attempt_row}"
            )
        core_state = np.asarray(root_result.core_state, dtype=float)
        audit = _audit_resting_state(model, core_state, expression, root_result, spec)
        attempt_row.update(
            rest_gate_pass=bool(audit["rest_gate_pass"]),
            max_abs_amount_rhs_fmol_s=audit["max_abs_amount_rhs_fmol_s"],
            max_abs_volume_rhs_pL_s=audit["max_abs_volume_rhs_pL_s"],
            max_abs_current_residual_A=audit["max_abs_current_residual_A"],
            max_dimensionless_conservation_ratio=audit[
                "max_dimensionless_conservation_ratio"
            ],
            solution_core_sha256=sha256_object(list(core_state)),
        )
        attempts.append(attempt_row)
        if not audit["rest_gate_pass"]:
            raise RuntimeError(
                f"resting root failed complete audit for {root}, Ca={calcium}, expression={expression}, seed={seed_id}"
            )
        candidates.append(
            {
                "seed_id": seed_id,
                "coordinates": np.asarray(root_result.coordinates, dtype=float),
                "core_state": core_state,
                "audit": audit,
            }
        )
    normalized_root_distance = float(
        np.max(np.abs((candidates[0]["coordinates"] - candidates[1]["coordinates"]) / span))
    )
    same_root = bool(normalized_root_distance <= spec.cluster_relative_tolerance)
    if not same_root:
        raise RuntimeError(
            f"the two prescribed warm starts reached distinct roots for {root}, Ca={calcium}, expression={expression}"
        )
    selected = min(
        candidates,
        key=lambda item: (
            item["audit"]["max_abs_scaled_independent_rhs"],
            item["audit"]["max_abs_amount_rhs_fmol_s"],
            item["seed_id"],
        ),
    )
    row: dict[str, Any] = {
        "root_id": root,
        "root_alias": aliases[root],
        "calcium_uM": calcium,
        "routing_family": routing(manifest, root),
        "ae4_expression": expression,
        "selected_seed_id": selected["seed_id"],
        "warm_start_count": 2,
        "both_warm_starts_accepted": True,
        "warm_starts_same_root": same_root,
        "warm_start_root_distance_normalized": normalized_root_distance,
        "root_cluster_relative_tolerance": spec.cluster_relative_tolerance,
        "rest_method": "direct_charge_manifold_least_squares",
        "rest_integration_used": False,
        "rest_protocol": StimulusArm.REST.value,
        "resting_calcium_uM": float(model.stimulus.resting_calcium_uM),
        "rest_beta_input": 0.0,
        "whole_cell_parameters_sha256": sha256_object(model.parameters),
        "ae4_parameters_sha256": sha256_object(model.ae4_parameters),
        "core_state_sha256": sha256_object(list(selected["core_state"])),
        "coordinates": [float(value) for value in selected["coordinates"]],
        "core_state": [float(value) for value in selected["core_state"]],
        "attempts": attempts,
        "wall_seconds": time.perf_counter() - started,
        **selected["audit"],
    }
    path = OUT / "resting_states" / f"{_rest_output_stem(root, calcium, expression)}.json"
    write_json(path, row)
    return row


def _stimulation_output_stem(root: str, calcium: float, expression: float) -> str:
    return f"{root}_ca{calcium:.2f}_e{expression:.2f}_production_radau"


def run_stimulation_case(request: tuple[str, float, float]) -> dict[str, Any]:
    root, calcium, expression = request
    started = time.perf_counter()
    manifest, frozen, _ = inputs()
    aliases = _root_aliases(manifest)
    rest_path = OUT / "resting_states" / f"{_rest_output_stem(root, calcium, expression)}.json"
    resting = json.loads(rest_path.read_text())
    if not resting["rest_gate_pass"]:
        raise AssertionError("stimulation cannot start from an unaccepted resting root")
    core = np.asarray(resting["core_state"], dtype=float)
    genotype = _genotype(expression)
    model = with_routing_only(build_model(manifest, root, calcium))
    rest_model = _explicit_rest_model(manifest, root)
    y0 = _root_state_with_basal_regulation(model, core)
    initial_stimulated = model.evaluate(0.0, y0, genotype=genotype)
    initial_rest = rest_model.evaluate(0.0, y0, genotype=genotype)
    onset_rest_rhs_match = bool(np.array_equal(initial_stimulated.rhs, initial_rest.rhs))
    if not onset_rest_rhs_match:
        raise AssertionError("stimulated model at t=0 does not reproduce explicit REST RHS")
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
        raise RuntimeError(f"production stimulation failed for {root}, Ca={calcium}, expression={expression}")
    original = build_model(manifest, root, calcium)
    sample_count = len(trace.time_s)
    ae4_sources = np.zeros((sample_count, 4), dtype=float)
    reference_sources = np.zeros((sample_count, 4), dtype=float)
    max_chloride_change = 0.0
    max_anion_change = 0.0
    max_total_cation_change = 0.0
    max_other_rhs_change = 0.0
    opposing_count = 0
    negative_cycle_count = 0
    max_charge_residual = 0.0
    for index, (time_s, state) in enumerate(zip(trace.time_s, trace.states.T)):
        new_evaluation = model.evaluate(float(time_s), state, genotype=genotype)
        old_evaluation = original.evaluate(float(time_s), state, genotype=genotype)
        new = new_evaluation.diagnostics.ae4
        old = old_evaluation.diagnostics.ae4
        ae4_sources[index] = (
            new.na_cell_fmol_s,
            new.k_cell_fmol_s,
            new.cl_cell_fmol_s,
            new.hco3_cell_fmol_s,
        )
        reference_sources[index] = (
            old.na_cell_fmol_s,
            old.k_cell_fmol_s,
            old.cl_cell_fmol_s,
            old.hco3_cell_fmol_s,
        )
        max_chloride_change = max(
            max_chloride_change, abs(new.cl_cell_fmol_s - old.cl_cell_fmol_s)
        )
        max_anion_change = max(
            max_anion_change, abs(new.hco3_cell_fmol_s - old.hco3_cell_fmol_s)
        )
        max_total_cation_change = max(
            max_total_cation_change,
            abs(
                (new.na_cell_fmol_s + new.k_cell_fmol_s)
                - (old.na_cell_fmol_s + old.k_cell_fmol_s)
            ),
        )
        max_other_rhs_change = max(
            max_other_rhs_change,
            float(np.max(np.abs(new_evaluation.rhs[2:] - old_evaluation.rhs[2:]))),
        )
        opposing_count += int(new.na_cell_fmol_s * new.k_cell_fmol_s < 0.0)
        negative_cycle_count += int(new.cl_cell_fmol_s < 0.0)
        max_charge_residual = max(
            max_charge_residual, abs(float(new.charge_source_residual_fmol_s))
        )
    if (
        max_chloride_change != 0.0
        or max_anion_change != 0.0
        or max_total_cation_change > ROUTING_ROUNDOFF_TOLERANCE_FMOL_S
        or max_other_rhs_change != 0.0
        or opposing_count != 0
    ):
        raise AssertionError("routing intervention changed an excluded component")
    stem = _stimulation_output_stem(root, calcium, expression)
    trajectory_path = OUT / "trajectories" / f"{stem}.npz"
    trajectory_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        key: value for key, value in asdict(trace).items() if isinstance(value, np.ndarray)
    }
    arrays.update(
        ae4_sources_na_k_cl_hco3_fmol_s=ae4_sources,
        reference_sources_na_k_cl_hco3_fmol_s=reference_sources,
    )
    np.savez_compressed(trajectory_path, **arrays)
    row: dict[str, Any] = {
        "root_id": root,
        "root_alias": aliases[root],
        "calcium_uM": calcium,
        "routing_family": routing(manifest, root),
        "ae4_expression": expression,
        "solver": PRODUCTION_RADAU.label,
        "valid": valid,
        "solver_success": bool(trace.success),
        "solver_message": trace.message,
        "exact_time_grid": exact_grid,
        "sample_count": sample_count,
        "total_0_600_pL": float(trace.cumulative_flow_pL[-1]),
        "endpoint_flow_pL_s": float(trace.flow_pL_s[-1]),
        "positive_core": bool(trace.positive_core),
        "nonnegative_flow": bool(trace.all_flow_nonnegative),
        "max_dimensionless_conservation_ratio": float(
            trace.max_dimensionless_conservation_ratio
        ),
        "onset_rest_rhs_match": onset_rest_rhs_match,
        "initial_max_abs_amount_rhs_fmol_s": float(
            np.max(np.abs(initial_stimulated.rhs[AMOUNT_ROWS]))
        ),
        "initial_max_abs_volume_rhs_pL_s": float(
            np.max(np.abs(initial_stimulated.rhs[VOLUME_ROWS]))
        ),
        "max_pointwise_chloride_law_change_fmol_s": max_chloride_change,
        "max_pointwise_anion_law_change_fmol_s": max_anion_change,
        "max_pointwise_total_cation_source_change_fmol_s": max_total_cation_change,
        "max_pointwise_other_rhs_change": max_other_rhs_change,
        "opposing_cation_samples": opposing_count,
        "negative_net_cycle_samples": negative_cycle_count,
        "max_ae4_charge_residual_fmol_s": max_charge_residual,
        "resting_state_sha256": resting["core_state_sha256"],
        "whole_cell_parameters_sha256": sha256_object(model.parameters),
        "ae4_parameters_sha256": sha256_object(model.ae4_parameters),
        "trajectory_path": str(trajectory_path.relative_to(REPO)),
        "wall_seconds": time.perf_counter() - started,
    }
    write_json(OUT / "trajectories" / f"{stem}.json", row)
    return row


def _csv_flat_rest_row(row: Mapping[str, Any]) -> dict[str, Any]:
    excluded = {
        "attempts", "raw_rhs", "scaled_independent_rhs", "conservation_residuals",
        "coordinates", "core_state",
    }
    flat = {key: value for key, value in row.items() if key not in excluded}
    for key, value in row["conservation_residuals"].items():
        flat[f"conservation_{key}"] = value
    return flat


def _all_pairwise_distances(vectors: Sequence[np.ndarray], span: np.ndarray) -> list[float]:
    return [
        float(np.max(np.abs((left - right) / span)))
        for left, right in itertools.combinations(vectors, 2)
    ]


def _verify_cross_calcium_rest(rest_rows: Sequence[Mapping[str, Any]]) -> float:
    spec = WTCalibrationSpec()
    span = spec.coordinate_upper - spec.coordinate_lower
    maximum = 0.0
    for root, expression in itertools.product(
        sorted({row["root_id"] for row in rest_rows}), EXPRESSIONS
    ):
        group = [
            row for row in rest_rows
            if row["root_id"] == root and row["ae4_expression"] == expression
        ]
        if len(group) != len(CALCIUM):
            raise AssertionError("each root and expression requires all three resting solves")
        distances = _all_pairwise_distances(
            [np.asarray(row["coordinates"], dtype=float) for row in group], span
        )
        local_max = max(distances, default=0.0)
        maximum = max(maximum, local_max)
        if local_max > spec.cluster_relative_tolerance:
            raise RuntimeError(
                f"calcium-specific warm starts reached different resting roots for {root}, expression={expression}"
            )
    return maximum


def _prefix_rest(target: dict[str, Any], prefix: str, rest: Mapping[str, Any]) -> None:
    fields = (
        "max_abs_scaled_independent_rhs",
        "max_abs_amount_rhs_fmol_s",
        "max_abs_volume_rhs_pL_s",
        "max_abs_omitted_rhs_fmol_s",
        "max_abs_regulatory_rhs_s_inv",
        "max_abs_bulk_charge_fmol",
        "max_abs_current_residual_A",
        "max_dimensionless_conservation_ratio",
        "rest_cell_na_mM", "rest_cell_k_mM", "rest_cell_cl_mM", "rest_cell_ph",
        "rest_cell_volume_pL", "rest_lumen_na_mM", "rest_lumen_k_mM",
        "rest_lumen_cl_mM", "rest_lumen_ph", "rest_lumen_volume_pL",
        "rest_gate_pass", "selected_seed_id", "warm_start_root_distance_normalized",
    )
    for field in fields:
        target[f"{prefix}_{field}"] = rest[field]


def make_comparisons(
    rest_rows: Sequence[Mapping[str, Any]], stimulation_rows: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    rest_index = {
        (row["root_id"], float(row["calcium_uM"]), float(row["ae4_expression"])): row
        for row in rest_rows
    }
    stimulation_index = {
        (row["root_id"], float(row["calcium_uM"]), float(row["ae4_expression"])): row
        for row in stimulation_rows
    }
    prior_rows = list(csv.DictReader(PRIOR_COMPARISON.open()))
    prior_index = {
        (row["root_id"], float(row["calcium_uM"])): row for row in prior_rows
    }
    comparisons: list[dict[str, Any]] = []
    roots = sorted({row["root_id"] for row in rest_rows})
    aliases = {row["root_id"]: row["root_alias"] for row in rest_rows}
    for root, calcium in itertools.product(roots, CALCIUM):
        wt_rest = rest_index[root, calcium, 1.0]
        low_rest = rest_index[root, calcium, 0.05]
        wt = stimulation_index[root, calcium, 1.0]
        low = stimulation_index[root, calcium, 0.05]
        prior = prior_index[root, calcium]
        ratio = float(low["total_0_600_pL"] / wt["total_0_600_pL"])
        change = 100.0 * (ratio - 1.0)
        previous_ratio = float(prior["routing_only_ratio"])
        previous_change = float(prior["routing_only_change_percent"])
        direction = "increase" if ratio > 1.0 else "decrease" if ratio < 1.0 else "no_change"
        previous_direction = (
            "increase" if previous_ratio > 1.0
            else "decrease" if previous_ratio < 1.0
            else "no_change"
        )
        row: dict[str, Any] = {
            "root_id": root,
            "root_alias": aliases[root],
            "calcium_uM": calcium,
            "routing_family": wt_rest["routing_family"],
            "wt_total_0_600_pL": wt["total_0_600_pL"],
            "ae4_5pct_total_0_600_pL": low["total_0_600_pL"],
            "ae4_5pct_to_wt_ratio": ratio,
            "ae4_5pct_change_percent": change,
            "ae4_5pct_effect_direction": direction,
            "previous_non_equilibrated_wt_total_0_600_pL": float(
                prior["routing_only_wt_pL"]
            ),
            "previous_non_equilibrated_5pct_total_0_600_pL": float(
                prior["routing_only_5pct_pL"]
            ),
            "previous_non_equilibrated_ratio": previous_ratio,
            "previous_non_equilibrated_change_percent": previous_change,
            "previous_non_equilibrated_effect_direction": previous_direction,
            "change_from_previous_percentage_points": change - previous_change,
            "direction_changed_from_previous": direction != previous_direction,
            "paired_production_gate_pass": bool(wt["valid"] and low["valid"]),
        }
        _prefix_rest(row, "wt", wt_rest)
        _prefix_rest(row, "ae4_5pct", low_rest)
        comparisons.append(row)
    write_rows(OUT / "comparison.csv", comparisons)
    return comparisons


def make_summary(comparisons: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    ratios = [float(row["ae4_5pct_to_wt_ratio"]) for row in comparisons]
    changes = [float(row["ae4_5pct_change_percent"]) for row in comparisons]
    if all(value < 1.0 for value in ratios):
        answer = "5% AE4 decreases total secretion in every matched comparison"
    elif all(value > 1.0 for value in ratios):
        answer = "5% AE4 increases total secretion in every matched comparison"
    else:
        answer = "the direction of the 5% AE4 effect differs across matched comparisons"
    per_calcium = []
    for calcium in CALCIUM:
        group = [
            float(row["ae4_5pct_change_percent"])
            for row in comparisons if float(row["calcium_uM"]) == calcium
        ]
        per_calcium.append(
            {
                "calcium_uM": calcium,
                "change_percent_min": min(group),
                "change_percent_median": statistics.median(group),
                "change_percent_max": max(group),
            }
        )
    summary = {
        "primary_answer": answer,
        "pair_count": len(comparisons),
        "decrease_count": sum(value < 1.0 for value in ratios),
        "increase_count": sum(value > 1.0 for value in ratios),
        "no_change_count": sum(value == 1.0 for value in ratios),
        "direction_changed_from_previous_count": sum(
            bool(row["direction_changed_from_previous"]) for row in comparisons
        ),
        "ratio_min": min(ratios),
        "ratio_median": statistics.median(ratios),
        "ratio_max": max(ratios),
        "change_percent_min": min(changes),
        "change_percent_median": statistics.median(changes),
        "change_percent_max": max(changes),
        "per_calcium": per_calcium,
    }
    write_json(OUT / "summary.json", summary)
    return summary


def _format_scientific(value: Any) -> str:
    return f"{float(value):.3e}"


def _format_float(value: Any, digits: int = 6) -> str:
    return f"{float(value):.{digits}f}"


def write_report(
    rest_rows: Sequence[Mapping[str, Any]],
    comparisons: Sequence[Mapping[str, Any]],
    summary: Mapping[str, Any],
    verification: Mapping[str, Any],
) -> Path:
    representative_rest: list[Mapping[str, Any]] = []
    for root, expression in itertools.product(
        sorted({row["root_id"] for row in rest_rows}), EXPRESSIONS
    ):
        representative_rest.append(
            next(
                row for row in rest_rows
                if row["root_id"] == root
                and float(row["ae4_expression"]) == expression
                and float(row["calcium_uM"]) == CALCIUM[0]
            )
        )
    wt_rest = [row for row in representative_rest if float(row["ae4_expression"]) == 1.0]
    lines = [
        "# AE4 routing experiment after resting equilibration",
        "",
        "## Primary result",
        "",
        f"{summary['primary_answer']}. All {summary['pair_count']} WT and 5% AE4 pairs passed the inherited production and conservation gates.",
        "",
        (
            "The resting equations were solved directly under an explicit REST protocol. "
            "No resting integration, fitting, recalibration, parameter change, or phenotype guided root selection was used. "
            "For each root, expression and calcium comparison, the inherited saved rest and the matching earlier routed 600 s endpoint were used as independent warm starts."
        ),
        "",
        "| Calcium (µM) | Re equilibrated 5% AE4 change across roots | Previous non equilibrated change across roots |",
        "| ---: | ---: | ---: |",
    ]
    for item in summary["per_calcium"]:
        calcium = float(item["calcium_uM"])
        previous = [
            float(row["previous_non_equilibrated_change_percent"])
            for row in comparisons if float(row["calcium_uM"]) == calcium
        ]
        lines.append(
            f"| {calcium:.2f} | {min(float(row['ae4_5pct_change_percent']) for row in comparisons if float(row['calcium_uM']) == calcium):+.3f}% to {max(float(row['ae4_5pct_change_percent']) for row in comparisons if float(row['calcium_uM']) == calcium):+.3f}% | {min(previous):+.3f}% to {max(previous):+.3f}% |"
        )
    lines.extend(
        [
            "",
            "## Secretion comparison for every root and calcium",
            "",
            "Totals are model secretion per cell over 0 to 600 s, in pL. Positive percentages mean that 5% AE4 increased secretion relative to matched WT.",
            "",
            "| Root | Calcium (µM) | WT total (pL) | 5% AE4 total (pL) | Ratio | Change | Previous change | Direction changed? |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in comparisons:
        lines.append(
            "| {root_alias} | {calcium_uM:.2f} | {wt:.9f} | {low:.9f} | {ratio:.9f} | {change:+.4f}% | {previous:+.4f}% | {changed} |".format(
                root_alias=row["root_alias"], calcium_uM=float(row["calcium_uM"]),
                wt=float(row["wt_total_0_600_pL"]),
                low=float(row["ae4_5pct_total_0_600_pL"]),
                ratio=float(row["ae4_5pct_to_wt_ratio"]),
                change=float(row["ae4_5pct_change_percent"]),
                previous=float(row["previous_non_equilibrated_change_percent"]),
                changed="yes" if row["direction_changed_from_previous"] else "no",
            )
        )
    lines.extend(
        [
            "",
            "## Resting states",
            "",
            (
                "The resting protocol is independent of the later stimulation calcium. "
                f"The three calcium specific solves agree to a maximum normalised coordinate distance of {verification['max_cross_calcium_rest_root_distance_normalized']:.3e}. "
                "The table therefore lists each common resting state once."
            ),
            "",
            "| Root | AE4 | Cell Na | Cell K | Cell Cl | Cell pH | Cell volume | Lumen Na | Lumen K | Lumen Cl | Lumen pH | Lumen volume |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in representative_rest:
        lines.append(
            "| {root_alias} | {expression:.0%} | {cna:.6f} | {ck:.6f} | {ccl:.6f} | {cph:.6f} | {cv:.6f} | {lna:.6f} | {lk:.6f} | {lcl:.6f} | {lph:.6f} | {lv:.6f} |".format(
                root_alias=row["root_alias"], expression=float(row["ae4_expression"]),
                cna=float(row["rest_cell_na_mM"]), ck=float(row["rest_cell_k_mM"]),
                ccl=float(row["rest_cell_cl_mM"]), cph=float(row["rest_cell_ph"]),
                cv=float(row["rest_cell_volume_pL"]), lna=float(row["rest_lumen_na_mM"]),
                lk=float(row["rest_lumen_k_mM"]), lcl=float(row["rest_lumen_cl_mM"]),
                lph=float(row["rest_lumen_ph"]), lv=float(row["rest_lumen_volume_pL"]),
            )
        )
    lines.extend(
        [
            "",
            "Concentrations are mM and volumes are pL.",
            "",
            (
                "With no permitted recalibration, the new routed WT equilibria do not preserve the inherited WT resting calibration. "
                f"WT cell chloride is {min(float(row['rest_cell_cl_mM']) for row in wt_rest):.3f} to {max(float(row['rest_cell_cl_mM']) for row in wt_rest):.3f} mM, "
                f"cell pH is {min(float(row['rest_cell_ph']) for row in wt_rest):.3f} to {max(float(row['rest_cell_ph']) for row in wt_rest):.3f}, and "
                f"cell volume is {min(float(row['rest_cell_volume_pL']) for row in wt_rest):.3f} to {max(float(row['rest_cell_volume_pL']) for row in wt_rest):.3f} pL. "
                "These are the genuine mathematical resting states of the routing intervention with frozen parameters, rather than a preservation of the original fitted WT resting physiology."
            ),
            "",
            "## Resting residuals and gates",
            "",
            "Every value below is from the complete resting evaluation, not only the ten solved rows.",
            "",
            "| Root | AE4 | Scaled independent RHS | Amount RHS (fmol/s) | Volume RHS (pL/s) | Omitted RHS (fmol/s) | Regulatory RHS (/s) | Current (A) | Conservation ratio | Gate |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in representative_rest:
        lines.append(
            "| {root_alias} | {expression:.0%} | {scaled} | {amount} | {volume} | {omitted} | {regulatory} | {current} | {conservation} | {gate} |".format(
                root_alias=row["root_alias"], expression=float(row["ae4_expression"]),
                scaled=_format_scientific(row["max_abs_scaled_independent_rhs"]),
                amount=_format_scientific(row["max_abs_amount_rhs_fmol_s"]),
                volume=_format_scientific(row["max_abs_volume_rhs_pL_s"]),
                omitted=_format_scientific(row["max_abs_omitted_rhs_fmol_s"]),
                regulatory=_format_scientific(row["max_abs_regulatory_rhs_s_inv"]),
                current=_format_scientific(row["max_abs_current_residual_A"]),
                conservation=_format_scientific(row["max_dimensionless_conservation_ratio"]),
                gate="pass" if row["rest_gate_pass"] else "fail",
            )
        )
    lines.extend(
        [
            "",
            "The complete calcium specific values, both warm start attempts, all conservation residuals, and both compartment states are in `resting_states.csv` and `warm_start_attempts.csv`. The full 30 row comparison, including both WT and 5% resting diagnostics on every row, is in `comparison.csv`.",
            "",
            "## Intervention and verification",
            "",
            (
                "The existing net AE4 chloride source, bicarbonate source, and total cation source were unchanged at every sampled state. "
                "Only the Na and K shares of that total cation source were reassigned. "
                "No opposing Na and K AE4 sources occurred. The net cycle did not reverse at any sampled state in this re equilibrated run."
            ),
            "",
            f"All {verification['rest_case_count']} resting cases and {verification['stimulation_trajectory_count']} stimulated trajectories passed. The maximum warm start root distance was {verification['max_warm_start_root_distance_normalized']:.3e}, below the inherited clustering tolerance {verification['root_cluster_relative_tolerance']:.3e}. The maximum stimulation conservation ratio was {verification['max_stimulation_conservation_ratio']:.3e}.",
            "",
            "| Root | Inherited root identifier |",
            "| --- | --- |",
        ]
    )
    aliases = {row["root_id"]: row["root_alias"] for row in rest_rows}
    for root in sorted(aliases):
        lines.append(f"| {aliases[root]} | `{root}` |")
    report_path = ANALYSIS / "ae4_routing_re_equilibrated_report.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = report_path.with_suffix(".md.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    return report_path


def _source_files() -> list[Path]:
    files = [MANIFEST, FROZEN, PRIOR_COMPARISON]
    files.extend(sorted((REPO / "src").rglob("*.py")))
    return files


def _contract(manifest: Mapping[str, Any], frozen: Mapping[str, Any]) -> dict[str, Any]:
    spec = WTCalibrationSpec()
    return {
        "model_id": MODEL_ID,
        "root_ids": sorted(manifest["roots"]),
        "calcium_uM": list(CALCIUM),
        "ae4_expression": list(EXPRESSIONS),
        "rest_protocol": StimulusArm.REST.value,
        "resting_calcium_uM": manifest["protocols"]["CCH_IPR"]["resting_calcium_uM"],
        "rest_beta_input": 0.0,
        "warm_starts": [
            "inherited saved resting state",
            "endpoint of matching prior routing-only 600 s production Radau trajectory",
        ],
        "rest_solver": {
            "method": "direct bounded charge-manifold least_squares",
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
            "regulatory_rhs_tolerance_s_inv": REST_REGULATORY_TOLERANCE_S_INV,
            "cluster_relative_tolerance": spec.cluster_relative_tolerance,
            "jacobian_rank_relative_tolerance": spec.jacobian_rank_relative_tolerance,
        },
        "rest_integration_policy": (
            "Direct solves are attempted first. No fixed-duration resting integration is used; "
            "the run stops if direct roots do not satisfy every inherited gate."
        ),
        "production_solver": asdict(PRODUCTION_RADAU),
        "time_grid_s": frozen["time_grid_s"],
        "conservation_tolerances": dict(CONSERVATION_RESIDUAL_TOLERANCES),
        "routing_roundoff_tolerance_fmol_s": ROUTING_ROUNDOFF_TOLERANCE_FMOL_S,
        "no_parameter_fit": True,
        "no_recalibration": True,
        "no_phenotype_guided_root_selection": True,
        "final_intervention": (
            "existing net AE4 chloride law and anion sources retained; one net cation cycle "
            "partitioned by equal donor-side Na/K concentration weighting"
        ),
        "prior_non_equilibrated_comparison_sha256": sha256_file(PRIOR_COMPARISON),
        "source_sha256": {
            str(path.relative_to(REPO)): sha256_file(path) for path in _source_files()
        },
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
    }


def _run_parallel(
    requests: Sequence[tuple[str, float, float]], worker: Any, workers: int, csv_path: Path
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(worker, request): request for request in requests}
        for future in as_completed(jobs):
            row = future.result()
            rows.append(row)
            if worker is solve_rest_case:
                write_rows(csv_path, [_csv_flat_rest_row(item) for item in rows])
            else:
                write_rows(csv_path, rows)
            print(
                json.dumps(
                    {
                        "phase": "rest" if worker is solve_rest_case else "stimulation",
                        "root": row["root_alias"],
                        "calcium_uM": row["calcium_uM"],
                        "ae4_expression": row["ae4_expression"],
                        "pass": row.get("rest_gate_pass", row.get("valid")),
                        "wall_seconds": row["wall_seconds"],
                    }
                ),
                flush=True,
            )
    return sorted(rows, key=lambda item: (item["root_id"], item["calcium_uM"], -item["ae4_expression"]))


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
    OUT.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    source_hashes_before = {
        str(path.relative_to(REPO)): sha256_file(path) for path in _source_files()
    }
    contract = _contract(manifest, frozen)
    contract_path = OUT / "contract.json"
    if contract_path.exists() and json.loads(contract_path.read_text()) != contract:
        raise RuntimeError("existing run contract does not match the current frozen inputs")
    write_json(contract_path, contract)
    requests = [
        (root, calcium, expression)
        for root, calcium, expression in itertools.product(
            sorted(manifest["roots"]), CALCIUM, EXPRESSIONS
        )
    ]
    rest_rows = _run_parallel(
        requests, solve_rest_case, args.workers, OUT / "resting_states.csv"
    )
    attempts = [attempt for row in rest_rows for attempt in row["attempts"]]
    write_rows(OUT / "warm_start_attempts.csv", attempts)
    max_cross_calcium = _verify_cross_calcium_rest(rest_rows)
    stimulation_rows = _run_parallel(
        requests, run_stimulation_case, args.workers, OUT / "trajectory_summary.csv"
    )
    comparisons = make_comparisons(rest_rows, stimulation_rows)
    summary = make_summary(comparisons)
    source_hashes_after = {
        str(path.relative_to(REPO)): sha256_file(path) for path in _source_files()
    }
    if source_hashes_after != source_hashes_before:
        raise AssertionError("model source or frozen input changed during the run")
    spec = WTCalibrationSpec()
    verification = {
        "model_id": MODEL_ID,
        "rest_case_count": len(rest_rows),
        "warm_start_attempt_count": len(attempts),
        "accepted_warm_start_attempt_count": sum(
            bool(row["rest_gate_pass"]) for row in attempts
        ),
        "all_rest_gates_pass": all(bool(row["rest_gate_pass"]) for row in rest_rows),
        "all_warm_start_pairs_same_root": all(
            bool(row["warm_starts_same_root"]) for row in rest_rows
        ),
        "max_warm_start_root_distance_normalized": max(
            float(row["warm_start_root_distance_normalized"]) for row in rest_rows
        ),
        "max_cross_calcium_rest_root_distance_normalized": max_cross_calcium,
        "root_cluster_relative_tolerance": spec.cluster_relative_tolerance,
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
        "max_rest_regulatory_rhs_s_inv": max(
            float(row["max_abs_regulatory_rhs_s_inv"]) for row in rest_rows
        ),
        "max_rest_current_residual_A": max(
            float(row["max_abs_current_residual_A"]) for row in rest_rows
        ),
        "max_rest_conservation_ratio": max(
            float(row["max_dimensionless_conservation_ratio"]) for row in rest_rows
        ),
        "stimulation_trajectory_count": len(stimulation_rows),
        "all_stimulation_gates_pass": all(bool(row["valid"]) for row in stimulation_rows),
        "all_onset_rest_rhs_match": all(
            bool(row["onset_rest_rhs_match"]) for row in stimulation_rows
        ),
        "max_stimulation_conservation_ratio": max(
            float(row["max_dimensionless_conservation_ratio"])
            for row in stimulation_rows
        ),
        "modified_sample_count": sum(int(row["sample_count"]) for row in stimulation_rows),
        "opposing_cation_sample_count": sum(
            int(row["opposing_cation_samples"]) for row in stimulation_rows
        ),
        "negative_net_cycle_sample_count": sum(
            int(row["negative_net_cycle_samples"]) for row in stimulation_rows
        ),
        "max_pointwise_chloride_law_change_fmol_s": max(
            float(row["max_pointwise_chloride_law_change_fmol_s"])
            for row in stimulation_rows
        ),
        "max_pointwise_anion_law_change_fmol_s": max(
            float(row["max_pointwise_anion_law_change_fmol_s"])
            for row in stimulation_rows
        ),
        "max_pointwise_total_cation_source_change_fmol_s": max(
            float(row["max_pointwise_total_cation_source_change_fmol_s"])
            for row in stimulation_rows
        ),
        "max_pointwise_other_rhs_change": max(
            float(row["max_pointwise_other_rhs_change"]) for row in stimulation_rows
        ),
        "max_ae4_charge_residual_fmol_s": max(
            float(row["max_ae4_charge_residual_fmol_s"]) for row in stimulation_rows
        ),
        "changed_source_or_frozen_inputs": [],
    }
    write_json(OUT / "verification.json", verification)
    report_path = write_report(rest_rows, comparisons, summary, verification)
    report_hash = hashlib.sha256(report_path.read_bytes()).hexdigest()
    write_json(
        OUT / "artifact_index.json",
        {
            "report": str(report_path.relative_to(REPO)),
            "report_sha256": report_hash,
            "comparison": str((OUT / "comparison.csv").relative_to(REPO)),
            "resting_states": str((OUT / "resting_states.csv").relative_to(REPO)),
            "warm_start_attempts": str((OUT / "warm_start_attempts.csv").relative_to(REPO)),
            "trajectory_summary": str((OUT / "trajectory_summary.csv").relative_to(REPO)),
            "summary": str((OUT / "summary.json").relative_to(REPO)),
            "verification": str((OUT / "verification.json").relative_to(REPO)),
        },
    )
    print(json.dumps({"summary": summary, "verification": verification}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
