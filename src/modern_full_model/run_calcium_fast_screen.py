"""Task 13C: frozen WT calcium screen, with selective production confirmation.

Run from the repository root with BLAS thread counts set to one::

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python -m src.modern_full_model.run_calcium_fast_screen --workers 9

The default is the available CPU count, capped at ten. The first loose case
is timed serially and reused as one of the twenty primary trajectories.
Completed requests are cached with input and scientific payload hashes.
No resting solve, genotype driver, or target ledger is used.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

import numpy as np
import scipy

from .membranes import hill_activation
from .model import WT
from .native_source_panel import NativeRootRecord, build_native_source_model
from .native_dynamic_contract import (
    ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S, PL_S_CELL_TO_UL_MIN,
    assess_geometry_scale, geometry_contract_payload,
    shared_observation_scale_interval, sustainment_gate, sustainment_ratio,
)
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES, LOOSE_RADAU, PRODUCTION_BDF,
    PRODUCTION_RADAU, SecretagogueProtocol, SolverSpecification, StimulusArm,
    Trajectory, canonical_json_bytes, physical_time_grid,
    pre_reveal_regulatory_ensemble, sha256_file, sha256_object, simulate_wt,
    solver_comparison,
)


REPO = Path(__file__).resolve().parents[2]
INHERITED = REPO / "results/13B_modern_full_model"
OUTPUT = REPO / "results/13C_calcium_fast_screen"
BASE_COMMIT = "289a313fc0c7ba103d490681c43ab3987ca1b62e"
MANIFEST_SHA256 = "6d283961421a6718a54fd7c86e6575bd5c032a93e649923777b02295cee97f29"
MEMBER_ID = "R1_G125_P21_PROSE_TREFERENCE"
NEW_CALCIUM = (0.25, 0.50)
MINUTES = tuple(range(60, 601, 60))
SOLVERS = {s.label: s for s in (LOOSE_RADAU, PRODUCTION_RADAU, PRODUCTION_BDF)}
REQUIRED_FLOW = 9.0 / ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
DECISION_RULE = {
    "version": "TASK13C_FAST_V1",
    "near_boundary_fraction": 0.90,
    "minimum_minute_flow_boundary_pL_s": REQUIRED_FLOW,
    "near_boundary_minimum_minute_flow_pL_s": 0.90 * REQUIRED_FLOW,
    "all_roots_production_trigger": (
        "Any loose absolute geometry pass, any minimum minute flow >= 90% "
        "of 9/one_SMG_scale_max, or any inconclusive numerical result. "
        "Shape or sustainment failure cannot suppress this trigger."
    ),
    "clear_025_failure": "No production trajectories at 0.25.",
    "clear_050_failure": "Production Radau for the best flowing 0.50 root only.",
    "best_root_order": "Largest minimum minute flow, then mean flow, then root ID ascending.",
    "bdf_rule": (
        "At each production amplitude, crosscheck the passing root closest "
        "to the boundary; if none passes, the best flowing root. If production "
        "roots disagree, also crosscheck the best flowing failing root."
    ),
    "pilot_rule": "Lexicographically first root at 0.25; reuse within the twenty cases.",
    "unexpected_pilot_wall_seconds": 120.0,
    "threshold_refinement": False,
}


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(canonical_json_bytes(payload) + b"\n")
    temporary.replace(path)


def write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def preferred_member() -> Any:
    return next(m for m in pre_reveal_regulatory_ensemble() if m.member_id == MEMBER_ID)


def load_freeze() -> tuple[dict[str, Any], dict[str, Any]]:
    """Verify saved evidence; use final manifest roots, never recalibrate."""
    source = INHERITED / "native_dynamic_contract_manifest.json"
    if sha256_file(source) != MANIFEST_SHA256:
        raise ValueError("Task 13B final manifest changed")
    manifest = json.loads(source.read_text())
    if len(manifest["roots"]) != 10 or not all(manifest["native_final_ready_flags"].values()):
        raise ValueError("Expected all ten final confirmed WT roots")
    verified = {str(source.relative_to(REPO)): MANIFEST_SHA256}
    for name in ("native_dynamic_contract_profile.csv",
                 "native_dynamic_contract_group_gate.csv",
                 "native_dynamic_contract_solver_crosschecks.csv"):
        digest = sha256_file(INHERITED / name)
        if digest != manifest["wt_dynamic_report_artifact_sha256"][name]:
            raise ValueError(f"Frozen dynamics artifact changed: {name}")
        verified[str((INHERITED / name).relative_to(REPO))] = digest
    discrepancies = []
    for name, expected in manifest["equation_source_file_sha256"].items():
        actual = sha256_file(REPO / name)
        if actual != expected:
            # This orchestration file was already different at branch intake.
            # It is neither imported nor executed by this new runner.
            if name != "src/modern_full_model/run_native_dynamic_contract.py":
                raise ValueError(f"Frozen implementation changed: {name}")
            discrepancies.append({"path": name, "expected_sha256": expected,
                                  "actual_sha256": actual, "used_for_simulation": False})
        verified[name] = actual
    for name in ("native_source_wt_roots.csv", "native_source_wt_summary.json"):
        actual = sha256_file(INHERITED / name)
        expected = manifest["native_input_artifact_sha256"][name]
        if actual != expected:
            discrepancies.append({"path": str((INHERITED / name).relative_to(REPO)),
                                  "expected_sha256": expected, "actual_sha256": actual,
                                  "used_for_simulation": False})
    for root_id, payload in manifest["roots"].items():
        for key in ("native_root_object", "core_state", "whole_cell_parameters", "ae4_parameters"):
            if sha256_object(payload[key]) != payload[key + "_sha256"]:
                raise ValueError(f"Invalid frozen {key}: {root_id}")
        root = NativeRootRecord(**payload["native_root_object"])
        if root.root_id != root_id or not root.passes_wt_gate or not root.passes_numerical_gate:
            raise ValueError("Invalid WT root intake")
        if sha256_object(root.core_state) != payload["core_state_sha256"]:
            raise ValueError("Root object and frozen core state disagree")
    member = preferred_member()
    if sha256_object(asdict(member.regulatory_model)) != sha256_object(
        manifest["regulatory_members"][MEMBER_ID]["model"]
    ):
        raise ValueError("Preferred regulatory parameters changed")
    if geometry_contract_payload() != manifest["geometry"]:
        raise ValueError("One SMG geometry changed")
    expected_solvers = manifest["production_solvers"] + manifest["tolerance_solvers_reference_r3_only"]
    for spec in SOLVERS.values():
        if asdict(spec) not in expected_solvers:
            raise ValueError("Inherited solver settings changed")
    return manifest, {"verified_input_sha256": verified,
                      "inherited_provenance_discrepancies": discrepancies}


@dataclass(frozen=True)
class Request:
    root_id: str
    calcium_uM: float
    solver_label: str = LOOSE_RADAU.label

    def __post_init__(self) -> None:
        if self.calcium_uM not in NEW_CALCIUM:
            raise ValueError("Only 0.25 and 0.50 uM may be newly simulated")
        if self.solver_label not in SOLVERS:
            raise ValueError("Undeclared solver")


def primary_requests(manifest: Mapping[str, Any]) -> list[Request]:
    return [Request(root_id, ca) for ca in NEW_CALCIUM for root_id in sorted(manifest["roots"])]


def build_model(manifest: Mapping[str, Any], request: Request) -> Any:
    payload = manifest["roots"][request.root_id]
    root = NativeRootRecord(**payload["native_root_object"])
    values = dict(manifest["protocols"]["CCH_IPR"])
    values.update(arm=StimulusArm.CCH_IPR, stimulated_calcium_uM=request.calcium_uM)
    protocol = SecretagogueProtocol(**values)
    template = build_native_source_model(root, regulatory_model=preferred_member().regulatory_model,
                                         stimulus=protocol)
    for value, name in ((template.parameters, "whole_cell_parameters"),
                        (template.ae4_parameters, "ae4_parameters")):
        if sha256_object(value) != payload[name + "_sha256"]:
            raise ValueError(f"Reconstructed {name} differs from Task 13B")
    # This is a frozen NKCC normalisation parameter, NOT the new Ca input.
    return attach_stimulated_nkcc1(template, N1AlgebraicNkcc1(
        fully_activated_multiplier=manifest["stimulated_nkcc1"]["fully_activated_multiplier"],
        resting_calcium_uM=manifest["protocols"]["CCH_IPR"]["resting_calcium_uM"],
        stimulated_calcium_uM=manifest["stimulated_nkcc1"]["effective_calcium_uM"],
    ))


def summarise(trajectory: Trajectory) -> dict[str, Any]:
    """Existing physical grid, accounting, sustainment and minute scale gates."""
    t = trajectory
    row: dict[str, Any] = {
        "success": t.success, "message": t.message, "output_samples": len(t.time_s),
        "positive_core": t.positive_core,
        "all_flow_nonnegative": bool(np.all(t.flow_pL_s >= 0)),
        "physiology_domain_gate_pass": bool(t.success and t.positive_core),
        "physiology_scope": "Inherited physical domains; no new stimulated ion target ranges",
    }
    ratios = {name: t.max_abs_conservation_residuals[name] / tolerance
              for name, tolerance in CONSERVATION_RESIDUAL_TOLERANCES.items()}
    row["max_dimensionless_conservation_ratio"] = max(ratios.values())
    row["conservation_gate_pass"] = all(math.isfinite(x) and x <= 1 for x in ratios.values())
    row["numerical_gate_pass"] = bool(t.success and t.positive_core
        and row["all_flow_nonnegative"] and row["conservation_gate_pass"])
    for name, value in t.max_abs_conservation_residuals.items():
        row["max_abs_residual__" + name] = value
    row.update(one_smg_geometry_gate_pass=False, sustainment_gate_pass=False,
               scientific_gate_pass=False)
    if not t.success:
        return row
    minute_flow = np.interp(MINUTES, t.time_s, t.flow_pL_s)
    row.update({f"flow_{minute}s_pL_s": float(q) for minute, q in zip(MINUTES, minute_flow)})
    row.update(
        endpoint_flow_pL_s=float(t.flow_pL_s[-1]),
        mean_flow_0_600_pL_s=float(t.cumulative_flow_pL[-1] / 600),
        cumulative_flow_600_pL=float(t.cumulative_flow_pL[-1]),
        minimum_post_stimulation_flow_pL_s=float(min(t.flow_pL_s[t.time_s > 0])),
        maximum_post_stimulation_flow_pL_s=float(max(t.flow_pL_s[t.time_s > 0])),
        minimum_minute_flow_pL_s=float(min(minute_flow)),
        maximum_minute_flow_pL_s=float(max(minute_flow)),
        capacity_multiplier_endpoint=float(t.capacity_multiplier[-1]),
    )
    for name, values in (("na_mM", t.cell_na_mM), ("k_mM", t.cell_k_mM),
                         ("cl_mM", t.cell_cl_mM), ("ph", t.cell_ph),
                         ("volume_pL", t.cell_volume_pL)):
        row["endpoint_" + name] = float(values[-1])
        row["minimum_" + name] = float(min(values))
        row["maximum_" + name] = float(max(values))
    if np.all(t.flow_pL_s > 0):
        ratio = sustainment_ratio(t.time_s, t.flow_pL_s)
        row.update(sustainment_ratio_q600_q60=ratio, sustainment_gate_pass=sustainment_gate(ratio))
    if np.all(minute_flow > 0):
        interval = shared_observation_scale_interval(minute_flow)
        geometry = assess_geometry_scale(interval)
        row.update(
            shared_scale_lower_uL_min_per_pL_s=interval.lower_uL_min_per_pL_s,
            shared_scale_upper_uL_min_per_pL_s=interval.upper_uL_min_per_pL_s,
            shared_scale_interval_nonempty=interval.compatible,
            one_smg_geometry_gate_pass=geometry.geometry_gate_pass,
            one_smg_feasible_scale_lower=geometry.feasible_lower_uL_min_per_pL_s,
            one_smg_feasible_scale_upper=geometry.feasible_upper_uL_min_per_pL_s,
            required_cell_count_at_minimum_minute=interval.lower_uL_min_per_pL_s / PL_S_CELL_TO_UL_MIN,
            implied_cell_count_at_mean_9_5=9.5 / row["mean_flow_0_600_pL_s"] / PL_S_CELL_TO_UL_MIN,
            maximum_gland_flow_at_minimum_minute_uL_min=min(minute_flow) * ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S,
            one_smg_scale_ceiling=ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S,
        )
    row["scientific_gate_pass"] = bool(row["numerical_gate_pass"]
        and row["sustainment_gate_pass"] and row["one_smg_geometry_gate_pass"])
    return row


def decode_trajectory(payload: Mapping[str, Any]) -> Trajectory:
    values = dict(payload)
    values["solver"] = SolverSpecification(**values["solver"])
    values["state_names"] = tuple(values["state_names"])
    for key in ("time_s", "states", "flow_pL_s", "cumulative_flow_pL", "cell_na_mM",
                "cell_k_mM", "cell_cl_mM", "cell_ph", "cell_volume_pL", "capacity_multiplier"):
        values[key] = np.asarray(values[key], dtype=float)
    return Trajectory(**values)


def scientific_payload(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: scientific_payload(v) for k, v in value.items()
                if not k.endswith("wall_seconds") and k != "scientific_sha256"}
    if isinstance(value, list):
        return [scientific_payload(v) for v in value]
    return value


def frozen_parameter_hash(model: Any) -> str:
    """Hash all non calcium inputs, excluding the sole experimental variable."""
    protocol = asdict(model.stimulus)
    protocol.pop("stimulated_calcium_uM")
    return sha256_object({
        "whole_cell_parameters": model.parameters,
        "ae4_parameters": model.ae4_parameters,
        "ae4_regulation": model.regulatory_model.ae4_regulatory_model,
        "nkcc1_regulation": model.regulatory_model.nkcc1_regulatory_model,
        "non_calcium_protocol": protocol,
    })


def simulate_request(job: tuple[Mapping[str, Any], Request]) -> dict[str, Any]:
    manifest, request = job
    started = time.perf_counter()
    model = build_model(manifest, request)
    counters: dict[str, Any] = {}
    solve = model.solve_dynamics

    def measured_solve(*args: Any, **kwargs: Any) -> Any:
        tick = time.perf_counter()
        result = solve(*args, **kwargs)
        counters.update(integrator_wall_seconds=time.perf_counter() - tick,
                        nfev=int(result.nfev), njev=int(result.njev), nlu=int(result.nlu))
        return result

    # Observe the existing solver result without changing equations or options.
    model.solve_dynamics = measured_solve
    frozen = manifest["roots"][request.root_id]
    p = model.parameters.membranes
    row = {**asdict(request), "genotype": "WT", "arm": "CCH_IPR",
           "regulatory_member_id": MEMBER_ID,
           "solver_method": SOLVERS[request.solver_label].method,
           "core_state_sha256": frozen["core_state_sha256"],
           "frozen_parameter_sha256": frozen_parameter_hash(model),
           "whole_cell_parameters_sha256": frozen["whole_cell_parameters_sha256"],
           "ae4_parameters_sha256": frozen["ae4_parameters_sha256"],
           "calcium_gate_activation": hill_activation(request.calcium_uM, p.calcium_half_uM, p.calcium_hill)}
    trajectory = None
    try:
        trajectory = simulate_wt(model, frozen["core_state"], family=MEMBER_ID,
            initial_condition="native_source_rest_root", solver=SOLVERS[request.solver_label],
            time_s=physical_time_grid(duration_s=600.0, sample_step_s=5.0))
        row.update(summarise(trajectory))
        if trajectory.success:
            diagnostics = model.evaluate(600.0, trajectory.states[:, -1], genotype=WT).diagnostics
            for name in ("apical", "basolateral", "transepithelial"):
                row[f"endpoint_v_{name}_mV"] = 1000 * getattr(diagnostics.membranes, f"v_{name}_V")
    except Exception as exc:
        row.update(success=False, numerical_gate_pass=False, scientific_gate_pass=False,
                   one_smg_geometry_gate_pass=False, sustainment_gate_pass=False,
                   message=f"{type(exc).__name__}: {exc}")
    row.update(counters, wall_seconds=time.perf_counter() - started)
    result = json.loads(canonical_json_bytes({"row": row, "trajectory": trajectory,
        "request": asdict(request), "input_manifest_sha256": MANIFEST_SHA256}))
    result["scientific_sha256"] = sha256_object(scientific_payload(result))
    return result


def run_requests(manifest: Mapping[str, Any], requests: Sequence[Request], workers: int,
                 worker: Any = simulate_request) -> list[dict[str, Any]]:
    jobs = [(manifest, request) for request in requests]
    if workers == 1:
        return [worker(job) for job in jobs]
    with ProcessPoolExecutor(max_workers=min(workers, 10)) as pool:
        return list(pool.map(worker, jobs, chunksize=1))


def best_row(rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    return min(rows, key=lambda row: (-float(row.get("minimum_minute_flow_pL_s", -math.inf)),
                                      -float(row.get("mean_flow_0_600_pL_s", -math.inf)),
                                      row["root_id"]))


def production_requests(rows: Sequence[Mapping[str, Any]]) -> tuple[list[Request], list[dict[str, Any]]]:
    requests, decisions = [], []
    for ca in NEW_CALCIUM:
        group = [r for r in rows if r["calcium_uM"] == ca]
        if len(group) != 10 or len({r["root_id"] for r in group}) != 10:
            raise ValueError("Confirmation requires the complete twenty case screen")
        trigger = any(not r["numerical_gate_pass"] or r["one_smg_geometry_gate_pass"]
                      or r.get("minimum_minute_flow_pL_s", 0) >= 0.9 * REQUIRED_FLOW for r in group)
        selected = sorted(group, key=lambda r: r["root_id"]) if trigger else (
            [best_row(group)] if ca == 0.50 else [])
        requests.extend(Request(r["root_id"], ca, PRODUCTION_RADAU.label) for r in selected)
        decisions.append({"calcium_uM": ca, "all_roots_triggered": trigger,
                          "production_root_ids": [r["root_id"] for r in selected],
                          "rule_sha256": sha256_object(DECISION_RULE)})
    return requests, decisions


def bdf_requests(rows: Sequence[Mapping[str, Any]]) -> list[Request]:
    requests = []
    for ca in NEW_CALCIUM:
        group = [r for r in rows if r["calcium_uM"] == ca]
        if not group:
            continue
        passed = [r for r in group if r["scientific_gate_pass"]]
        failed = [r for r in group if not r["scientific_gate_pass"]]
        if passed:
            chosen = min(passed, key=lambda r: (r["minimum_minute_flow_pL_s"], r["root_id"]))
            selected = [chosen] + ([best_row(failed)] if failed else [])
        else:
            selected = [best_row(group)]
        requests.extend(Request(r["root_id"], ca, PRODUCTION_BDF.label) for r in selected)
    return requests


def load_reference(manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Load the frozen 0.10 rows; this function cannot call an integrator."""
    rows = []
    with (INHERITED / "native_dynamic_contract_profile.csv").open(newline="") as handle:
        for raw in csv.DictReader(handle):
            if (raw["root_id"] in manifest["roots"] and raw["regulatory_member_id"] == MEMBER_ID
                and raw["arm"] == "CCH_IPR" and raw["solver_label"] == PRODUCTION_RADAU.label
                and raw["initial_condition"] == "native_source_rest_root"):
                row = dict(raw)
                row.update(calcium_uM=0.10, data_origin="FROZEN_TASK13B_NOT_RECOMPUTED",
                           mean_flow_0_600_pL_s=float(raw["cumulative_flow_600_pL"]) / 600)
                rows.append(row)
    if len(rows) != 10 or len({r["root_id"] for r in rows}) != 10:
        raise ValueError("Expected exactly ten frozen preferred reference rows")
    return sorted(rows, key=lambda r: r["root_id"])


def available_workers() -> int:
    cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    return max(1, min(10, cores))


def cache_path(request: Request, output: Path = OUTPUT) -> Path:
    return output / "trajectories" / (request.root_id.replace(":", "_")
        + f"_ca{request.calcium_uM:.2f}_{request.solver_label}.json")


def run_cached(manifest: Mapping[str, Any], requests: Sequence[Request], workers: int) -> list[dict[str, Any]]:
    cached, pending = {}, []
    for request in requests:
        path = cache_path(request)
        if path.exists():
            result = json.loads(path.read_text())
            if (result["request"] != asdict(request) or result["input_manifest_sha256"] != MANIFEST_SHA256
                or sha256_object(scientific_payload(result)) != result["scientific_sha256"]):
                raise ValueError(f"Invalid cached trajectory: {path.name}")
            cached[request] = result
        else:
            pending.append(request)
    if pending:
        # Both production modes use the same tested stateless dispatcher.
        for request, result in zip(pending, run_requests(manifest, pending, workers)):
            write_json(cache_path(request), result)
            cached[request] = result
    return [cached[request] for request in requests]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=available_workers())
    parser.add_argument("--phase", choices=("prepare", "pilot", "screen", "confirm", "all"), default="all")
    args = parser.parse_args(argv)
    if not 1 <= args.workers <= available_workers():
        parser.error(f"workers must be between 1 and {available_workers()} available cores")
    if any(os.environ.get(k) != "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")):
        parser.error("Set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 before Python starts")
    manifest, verification = load_freeze()
    plan = {"task": "13C_FAST_CALCIUM_SCREEN", "base_commit": BASE_COMMIT,
            "task13b_manifest_sha256": MANIFEST_SHA256, **verification,
            "root_intake": "Self hashed complete root objects in final Task 13B manifest",
            "root_ids": sorted(manifest["roots"]),
            "root_hashes": {rid: {k: v for k, v in p.items() if k.endswith("_sha256")}
                            for rid, p in manifest["roots"].items()},
            "reference_calcium_uM": 0.10, "reference_recomputed": False,
            "new_calcium_uM": NEW_CALCIUM, "primary_trajectory_count": 20,
            "protocol": manifest["protocols"]["CCH_IPR"], "regulatory_member_id": MEMBER_ID,
            "screening_solver": asdict(LOOSE_RADAU),
            "geometry": manifest["geometry"], "decision_rule": DECISION_RULE,
            "decision_rule_sha256": sha256_object(DECISION_RULE),
            "holdout_accessed": False, "genotype": "WT",
            "implementation_change": "New orchestration and counter capture only",
            "sample_grid_s": physical_time_grid(600.0, 5.0).tolist(),
            "runtime_environment": {"numpy": np.__version__, "scipy": scipy.__version__,
                                    "available_workers": available_workers()}}
    target = OUTPUT / "frozen_manifest.json"
    if target.exists() and json.loads(target.read_text()) != json.loads(canonical_json_bytes(plan)):
        raise ValueError("Predeclared freeze or rule changed; do not reuse stale results")
    write_json(target, plan)
    write_rows(OUTPUT / "calcium_reference_frozen.csv", load_reference(manifest))
    print(json.dumps({"phase": args.phase, "workers": args.workers,
                      "rule_sha256": sha256_object(DECISION_RULE)}), flush=True)
    if args.phase == "prepare":
        return 0
    requests = primary_requests(manifest)
    pilot = run_cached(manifest, requests[:1], 1)[0]
    print(json.dumps({"pilot": pilot["row"]}), flush=True)
    if args.phase == "pilot":
        return 0
    if pilot["row"]["wall_seconds"] > DECISION_RULE["unexpected_pilot_wall_seconds"]:
        raise RuntimeError("Pilot exceeded declared runtime threshold; profile before continuing")
    if args.phase == "confirm" and any(not cache_path(r).exists() for r in requests):
        raise ValueError("Complete the primary screen before confirmation")
    screen_pending = sum(not cache_path(r).exists() for r in requests)
    tick = time.perf_counter()
    screen = run_cached(manifest, requests, args.workers)
    write_rows(OUTPUT / "calcium_screen.csv", [r["row"] for r in screen])
    selection, decisions = production_requests([r["row"] for r in screen])
    decision_path = OUTPUT / "confirmation_decisions.json"
    write_json(decision_path, {"decision_rule": DECISION_RULE, "decisions": decisions,
                              "radau_requests": [asdict(r) for r in selection]})
    print(json.dumps({"screen_complete": len(screen), "production_decisions": decisions}), flush=True)
    if screen_pending or not (OUTPUT / "screen_runtime.json").exists():
        write_json(OUTPUT / "screen_runtime.json", {"workers": args.workers,
            "phase_wall_seconds": time.perf_counter() - tick,
            "summed_trajectory_wall_seconds": sum(r["row"]["wall_seconds"] for r in screen)})
    if args.phase == "screen":
        return 0
    tick = time.perf_counter()
    confirmation_pending = sum(not cache_path(r).exists() for r in selection)
    production = run_cached(manifest, selection, args.workers)
    cross_requests = bdf_requests([r["row"] for r in production])
    confirmation_pending += sum(not cache_path(r).exists() for r in cross_requests)
    crosschecks = run_cached(manifest, cross_requests, args.workers)
    by_key = {(r["row"]["root_id"], r["row"]["calcium_uM"]): r for r in production}
    for result in crosschecks:
        row = result["row"]
        reference = by_key[(row["root_id"], row["calcium_uM"])]
        if reference["trajectory"] and result["trajectory"] and reference["row"]["success"] and row["success"]:
            comparison = solver_comparison(decode_trajectory(reference["trajectory"]),
                                           decode_trajectory(result["trajectory"]))
            row.update({"crosscheck_" + k: v for k, v in comparison.items()})
            row["solver_crosscheck_gate_pass"] = bool(comparison["both_success"]
                and reference["row"]["numerical_gate_pass"] and row["numerical_gate_pass"] and all(
                comparison[k] <= manifest["solver_relative_tolerance"] for k in
                ("max_relative_state_difference", "max_relative_flow_difference",
                 "endpoint_cumulative_flow_relative_difference")))
        else:
            row["solver_crosscheck_gate_pass"] = False
    write_rows(OUTPUT / "calcium_confirmatory.csv", [r["row"] for r in production + crosschecks])
    write_json(decision_path, {"decision_rule": DECISION_RULE, "decisions": decisions,
        "radau_requests": [asdict(r) for r in selection], "bdf_requests": [asdict(r) for r in cross_requests]})
    if confirmation_pending or not (OUTPUT / "confirmation_runtime.json").exists():
        write_json(OUTPUT / "confirmation_runtime.json", {"workers": args.workers,
            "phase_wall_seconds": time.perf_counter() - tick,
            "production_radau_count": len(production), "production_bdf_count": len(crosschecks),
            "total_new_scientific_trajectory_count": len(screen) + len(production) + len(crosschecks)})
    print(json.dumps({"production_radau": len(production), "production_bdf": len(crosschecks),
                      "confirmation_complete": True}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
