"""Task 16 WT-only allocation, admission and freeze. No genotype outcome reads.

The scientific modules are imported unchanged. Only whole AE4 carrier amount
and whole NKCC1 capacity are replaced in newly constructed parameter objects.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from datetime import datetime, timezone
import csv
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from .calibration import WTCalibrationSpec, _observable_map, _wt_gate_failures
from .genotype_evaluation import (
    _attempt_and_root, _cluster_roots, _core_coordinates, _regulatory_suffix,
    _scaled_residual, _relative_difference, simulate_genotype,
)
from .model import ModernFullModel, WT
from .native_dynamic_contract import sustainment_gate, sustainment_ratio
from .nkcc_stimulation import StimulatedNkcc1Model
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import build_model as inherited_model, routing, slug, trajectory_path
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU, PRODUCTION_BDF,
    attach_basal_regulation, physical_time_grid, sha256_file, sha256_object,
)

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/16_wt_chloride_allocation"
BASE = "4d4403f279fc36aec940c9e4759486d02f07e3f6"
INTAKE = "fa59ef7670218125feccc5db72134995afd8fece"
BRANCH = "codex/task-16-wt-chloride-allocation"
CONDITIONS = ("inherited_baseline", "share_0.10", "share_0.30")
TARGETS = {"share_0.10": 0.10, "share_0.30": 0.30}
CALCIUM = (0.10, 0.25, 0.50)


def now():
    return datetime.now(timezone.utc).isoformat()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def read(path):
    return json.loads(Path(path).read_text())


def relative(path):
    return str(Path(path).relative_to(REPO))


def contract():
    path = OUT / "allocation_contract.json"
    return read(path), sha256_file(path)


def check_inputs():
    ledger = read(OUT / "baseline_input_hashes.json")
    assert ledger["baseline_commit"] == BASE
    for path, expected in ledger["input_files"].items():
        source = REPO / path
        assert source.stat().st_size == expected["bytes"], path
        assert sha256_file(source) == expected["sha256"], path
        assert git("hash-object", path) == expected["git_blob_sha1"], path
        assert git("rev-parse", BASE + ":" + path) == expected["git_blob_sha1"], path
    assert len(ledger["input_files"]) == 187
    return ledger


def scales(a, n, target):
    if not (math.isfinite(a) and math.isfinite(n) and a > 0 and n > 0):
        raise ValueError("Reference AE4 and NKCC1 must both load chloride")
    if target not in TARGETS.values():
        raise ValueError("Only allocation targets 0.10 and 0.30 are licensed")
    result = (target * (a + n) / a, (1 - target) * (a + n) / n)
    if not all(math.isfinite(x) and x >= 0 for x in result):
        raise ValueError("Nonfinite or negative capacity scale")
    return result


def allocated_model(manifest, root_id, calcium, ae4_scale, nkcc1_scale):
    if calcium not in CALCIUM:
        raise ValueError("Calcium outside frozen panel")
    if not all(math.isfinite(x) and x >= 0 for x in (ae4_scale, nkcc1_scale)):
        raise ValueError("Invalid complete-capacity scale")
    template = inherited_model(manifest, root_id, calcium)
    hp = replace(template.parameters.homeostasis,
        nkcc1_capacity_fmol_s=template.parameters.homeostasis.nkcc1_capacity_fmol_s * nkcc1_scale)
    parameters = replace(template.parameters, homeostasis=hp)
    ae4 = replace(template.ae4_parameters,
        carrier_amount_fmol=template.ae4_parameters.carrier_amount_fmol * ae4_scale)
    return StimulatedNkcc1Model(ModernFullModel(parameters=parameters,
        stimulus=template.stimulus, regulatory_model=template.regulatory_model,
        ae4_parameters=ae4, ae4_evaluator=template.ae4_evaluator))


def partition(model, core):
    state = attach_basal_regulation(model, core)
    ev = model.evaluate(0.0, state, genotype=WT)
    a, h = ev.diagnostics.ae4, ev.diagnostics.homeostasis
    A, N, E = a.cl_cell_fmol_s, 2 * h.nkcc1_inward_fmol_s, h.ae2_inward_fmol_s
    jna, jk = a.diagnostics["j_na_fmol_s"], a.diagnostics["j_k_fmol_s"]
    gross = abs(jna) + abs(jk)
    residuals = dict(ev.diagnostics.conservation_residuals)
    return {**_observable_map(model, state),
        "ae4_cl_fmol_s": A, "nkcc1_cl_fmol_s": N, "ae2_signed_cl_fmol_s": E,
        "total_basolateral_cl_fmol_s": A + N + E, "positive_pool_fmol_s": A + N,
        "positive_loading_eligible": bool(A > 0 and N > 0),
        "realized_ae4_share": A / (A + N) if A > 0 and N > 0 else None,
        "ae4_j_na_fmol_s": jna, "ae4_j_k_fmol_s": jk,
        "ae4_gross_branch_fmol_s": gross, "ae4_net_cl_fmol_s": A,
        "ae4_productive_fraction": abs(A) / gross if gross else None,
        "ae4_cancellation_fraction": 1 - abs(A) / gross if gross else None,
        "ae4_branches_oppose": bool(jna * jk < 0),
        "max_dimensionless_conservation_ratio": max(abs(v) / CONSERVATION_RESIDUAL_TOLERANCES[k]
            for k, v in residuals.items()),
        "max_abs_charge_fmol": max(abs(v) for v in ev.diagnostics.state_charge_fmol.values()),
        "max_abs_current_A": max(abs(v) for v in ev.diagnostics.membranes.current_residuals_A.values()),
    }


def prepare():
    if (OUT / "allocation_contract.json").exists():
        raise FileExistsError("Allocation contract is immutable")
    assert git("branch", "--show-current") == BRANCH
    assert git("diff", "--name-only", BASE, INTAKE, "src", "model", "archive", "data", "results") == ""
    inputs = check_inputs()
    manifest, inherited = load_freeze()
    rows, definitions = [], []
    for root_id in sorted(manifest["roots"]):
        model = inherited_model(manifest, root_id, 0.10)
        values = partition(model, manifest["roots"][root_id]["core_state"])
        rows.append({"root_id": root_id, "routing_family": routing(manifest, root_id), **values})
        for condition in CONDITIONS:
            target = TARGETS.get(condition)
            eligible = values["positive_loading_eligible"]
            lam, mu = (1.0, 1.0) if target is None else (
                scales(values["ae4_cl_fmol_s"], values["nkcc1_cl_fmol_s"], target)
                if eligible else (None, None))
            definitions.append({"candidate_id": slug(root_id) + "__" + condition,
                "root_id": root_id, "routing_family": routing(manifest, root_id),
                "condition": condition, "target_share": target,
                "inherited_share": values["realized_ae4_share"],
                "ae4_capacity_scale": lam, "nkcc1_capacity_scale": mu,
                "allocation_eligible": eligible})
    spec = WTCalibrationSpec()
    value = {"contract_id": "TASK16_WT_ONLY_THREE_CONDITIONS_V1", "created_utc": now(),
        "baseline_commit": BASE, "baseline_tree": git("rev-parse", BASE + "^{tree}"),
        "branch_intake_commit": INTAKE, "branch": BRANCH,
        "task15_override": "AGENTS.md accepts diagnostic conclusions; remote Task 15 SHA/artifacts not prerequisite; no rerun",
        "task15_accepted_classification": "No AE4 sign, source-mapping or expression-scaling defect; opposing mixed-cation branches",
        "instruction_sha256": {p: sha256_file(REPO / p) for p in
            ("AGENTS.md", "prompts/16_wt_constrained_ae4_chloride_allocation.md")},
        "conditions": CONDITIONS, "modified_target_shares": list(TARGETS.values()),
        "calcium_uM": CALCIUM, "candidate_definitions": definitions,
        "allocation_rule": "L0=A0+N0 with A0,N0>0; lambda=f*L0/A0; mu=(1-f)*L0/N0; signed AE2 excluded",
        "licensed_parameters": ["ae4_parameters.carrier_amount_fmol", "parameters.homeostasis.nkcc1_capacity_fmol_s"],
        "all_other_parameters_fixed": True, "no_genotype_specific_compensation": True,
        "reference_share_role": "Predeclared sensitivity coordinate; not a measured fraction or fitted target; realized share is reported separately",
        "wt_rest_gates_source": "WTCalibrationSpec and calibration._wt_gate_failures; genotype_evaluation._attempt_and_root numerical kernel",
        "wt_rest_spec": asdict(spec),
        "continuation": {"fixed_parameter_count_fitted": 0, "state_coordinate_count": 10,
            "local_starts": 3, "local_start_fraction": 0.001, "maximum_connected_step_fraction": 0.25,
            "max_nfev": 3000, "capacity_path": "lambda(t)=lambda_target**t; mu(t)=mu_target**t",
            "initial_internal_step": 0.25, "minimum_internal_step": 0.25 / 256,
            "retry_rule": "Halve homotopy step only on failed numerical solve or disconnected root; restore step up to 0.25 after success",
            "selection": "Nearest previous root by maximum coordinate difference/span, tie by residual; retain alternate roots",
            "internal_points": "Numerical solver states only. No dynamics or genotype simulation at internal t values."},
        "wt_dynamic_gates": ["complete 600s production Radau", "positive core", "finite nonnegative flow",
            "inherited conservation/current tolerances", "bulk charge <= inherited 1e-9 fmol",
            "WT sustainment Q600/Q60 in inherited [0.8,1.2]", "required BDF agreement"],
        "nonblocking_diagnostics": ["absolute one-SMG mapping", "shared minute observation scale", "WT secretion change from inherited baseline"],
        "no_new_stimulated_ion_bands": True,
        "production_solver": asdict(PRODUCTION_RADAU), "confirmation_solver": asdict(PRODUCTION_BDF),
        "conservation_tolerances": CONSERVATION_RESIDUAL_TOLERANCES,
        "time_grid_s": physical_time_grid().tolist(),
        "bdf_rule": "Lexicographically first rest-admissible root per modified share/routing at 0.50 uM; any disagreement expands to all rest-admissible roots/all calcium in that share/routing. Inherited baseline confirmations reused.",
        "solver_relative_tolerance": 1e-4,
        "wt_retention_rule": "Keep every rest-admissible root/share/calcium passing WT dynamic gates; no ranking or phenotype-based deletion",
        "freeze_rule": "All WT decisions, parameters, roots and trajectories hashed, committed and pushed before any modified genotype evaluation",
        "post_freeze_protocol": {"ae4_expression_grid": [round(1 - i / 20, 2) for i in range(20)],
            "ae2_expression_grid": [round(1 - i / 20, 2) for i in range(21)],
            "near_loss_rule": "New connected 0.05 AE4 resting state under allocated capacities; never exact AE4 zero",
            "ae2_rule": "Inherited target-free continuation to exact AE2 zero; complete transporter deletion only",
            "all_frozen_combinations_retained": True, "baseline": "Reuse byte/hash-valid Task14/14B results",
            "bdf_rule": "Same lexicographic frozen representative per share/routing at largest retained calcium for each genotype; disagreement expands share/routing"},
        "baseline_input_hashes_sha256": sha256_file(OUT / "baseline_input_hashes.json"),
        "baseline_input_count": len(inputs["input_files"]), "inherited_verification": inherited,
        "archive_tree": git("rev-parse", BASE + ":archive"),
        "regulatory_member": "R1_G125_P21_PROSE_TREFERENCE",
        "phenotype_firewall": "WT stage does not read genotype results or experimental phenotype data; byte hashing is not outcome access",
    }
    write_rows(OUT / "inherited_partition.csv", rows)
    write_json(OUT / "allocation_contract.json", value)
    print(json.dumps({"contract_sha256": contract()[1], "candidate_count": len(definitions)}), flush=True)


def require_contract_commit():
    _, digest = contract()
    assert git("show", "HEAD:" + relative(OUT / "allocation_contract.json"))
    assert git("hash-object", relative(OUT / "allocation_contract.json")) == git("rev-parse", "HEAD:" + relative(OUT / "allocation_contract.json"))
    return digest


def solve_local(model, previous, spec):
    lower, upper = spec.coordinate_lower, spec.coordinate_upper
    span = upper - lower
    starts = [("previous_connected_root", previous)]
    starts += [(name, np.clip(previous + direction * 0.001 * span,
        lower + 1e-10 * span, upper - 1e-10 * span))
        for name, direction in (("local_minus", -1), ("local_plus", 1))]
    attempts, roots = [], []
    for name, start in starts:
        attempt, root = _attempt_and_root(model=model, transporter="AE4", expression=1.0,
            genotype=WT, regulatory_suffix=_regulatory_suffix(model), spec=spec,
            start_id=name, start=start, lower=lower, upper=upper, max_nfev=3000)
        attempts.append(asdict(attempt))
        if root is not None:
            roots.append(root)
    roots = _cluster_roots(roots, expression=1.0, transporter="AE4", span=span,
        relative_tolerance=spec.cluster_relative_tolerance)
    selected = min(roots, key=lambda r: (np.max(np.abs(np.asarray(r.coordinates)-previous)/span),
        r.max_abs_scaled_independent_rhs)) if roots else None
    distance = float(np.max(np.abs(np.asarray(selected.coordinates)-previous)/span)) if selected else None
    return {"attempts": attempts, "roots": [asdict(r) for r in roots],
        "selected_root": asdict(selected) if selected else None,
        "distance_from_previous": distance, "connected": selected is not None and distance <= .25}


def root_path(candidate):
    return OUT / "wt_roots" / (candidate["candidate_id"] + ".json")


def rest_job(candidate):
    digest = require_contract_commit()
    path = root_path(candidate)
    if path.exists():
        saved = read(path); assert saved["contract_sha256"] == digest
        return saved["row"]
    manifest, _ = load_freeze()
    root_id = candidate["root_id"]
    core = manifest["roots"][root_id]["core_state"]
    spec = WTCalibrationSpec()
    steps = []
    numerical, selected, status = True, None, "INHERITED_ROOT_REUSED"
    start = time.monotonic()
    if candidate["condition"] != "inherited_baseline":
        if not candidate["allocation_eligible"]:
            numerical, status = False, "INELIGIBLE_REFERENCE_PARTITION"
        else:
            previous = _core_coordinates(inherited_model(manifest, root_id, .10), core)
            t, step = 0., .25
            while t < 1:
                target_t = min(1., t + step)
                model = allocated_model(manifest, root_id, .10,
                    candidate["ae4_capacity_scale"] ** target_t, candidate["nkcc1_capacity_scale"] ** target_t)
                report = solve_local(model, previous, spec)
                steps.append({"homotopy_t": target_t, "scientific_candidate": target_t == 1., **report})
                if report["connected"]:
                    selected = report["selected_root"]
                    core = selected["core_state"]
                    previous = np.asarray(selected["coordinates"])
                    t = target_t
                    step = min(.25, 2 * step)
                elif step > .25 / 256:
                    step /= 2
                else:
                    numerical, status = False, "NUMERICAL_CONTINUATION_FAILURE"
                    break
            if numerical:
                status = "ENDPOINT_ROOT_SOLVED"
    model = allocated_model(manifest, root_id, .10,
        candidate["ae4_capacity_scale"] or 1., candidate["nkcc1_capacity_scale"] or 1.)
    values = partition(model, core) if numerical else {}
    failures = list(_wt_gate_failures(values, spec)) if numerical else [status]
    if numerical and values["max_dimensionless_conservation_ratio"] > 1:
        failures.append("rest conservation/current gate")
    if numerical and values["max_abs_charge_fmol"] > spec.charge_tolerance_fmol:
        failures.append("rest bulk charge gate")
    if numerical and candidate["condition"] == "inherited_baseline":
        # Verify the fixed-parameter ten-state rank without re-solving the saved root.
        coords = _core_coordinates(model, core); span = spec.coordinate_upper-spec.coordinate_lower
        suffix = _regulatory_suffix(model)
        fun = lambda x: _scaled_residual(model, x, suffix, WT, spec)
        jac = np.column_stack([(fun(coords+np.eye(10)[i]*span[i]*1e-6)
            - fun(coords-np.eye(10)[i]*span[i]*1e-6))/(2e-6) for i in range(10)])
        singular = np.linalg.svd(jac, compute_uv=False)
        rank = int(sum(singular > spec.jacobian_rank_relative_tolerance * singular[0]))
        residual = float(np.max(np.abs(fun(coords))))
        if rank != 10 or residual > spec.root_scaled_tolerance:
            failures.append("inherited fixed-parameter root rank/residual")
    else:
        rank = selected["normalized_jacobian_rank"] if numerical else None
        residual = selected["max_abs_scaled_independent_rhs"] if numerical else None
    admitted = numerical and not failures
    row = {**candidate, **values, "numerical_rest_pass": numerical,
        "wt_rest_admissible": admitted, "status": status if admitted else (
            "WT_PHYSIOLOGY_REJECTED" if numerical else status), "gate_failures": "; ".join(failures),
        "fixed_parameter_rank": rank, "max_abs_scaled_rest_residual": residual,
        "continuation_step_attempts": len(steps), "root_payload_path": relative(path)}
    payload = {"contract_sha256": digest, "row": row, "created_utc": now(),
        "core_state": core if numerical else None, "core_state_sha256": sha256_object(core) if numerical else None,
        "whole_cell_parameters": asdict(model.parameters), "whole_cell_parameters_sha256": sha256_object(model.parameters),
        "ae4_parameters": asdict(model.ae4_parameters), "ae4_parameters_sha256": sha256_object(model.ae4_parameters),
        "regulatory_model_sha256": sha256_object(model.regulatory_model),
        "continuation_steps": steps, "wall_seconds": time.monotonic()-start}
    write_json(path, payload)
    print(json.dumps({k: row[k] for k in ("candidate_id", "status", "gate_failures")}), flush=True)
    return row


def trace_diagnostics(model, trace):
    states = np.asarray(trace["states"])
    cell_charge = states[0] + states[1] - states[2] - states[4] - model.parameters.geometry.fixed_cell_anion_equivalents_fmol
    lumen_charge = states[6] + states[7] - states[8] - states[10]
    charge = max(float(np.max(abs(cell_charge))), float(np.max(abs(lumen_charge))))
    ratio = sustainment_ratio(trace["time_s"], trace["flow_pL_s"])
    return {"wt_total_0_600_pL": trace["cumulative_flow_pL"][-1],
        "numerical_gate_pass": trace["numerical_gate_pass"],
        "max_dimensionless_conservation_ratio": trace["max_dimensionless_conservation_ratio"],
        "max_abs_charge_fmol": charge, "sustainment_Q600_Q60": ratio,
        "sustainment_gate_pass": sustainment_gate(ratio),
        "wt_dynamic_gate_without_bdf": bool(trace["numerical_gate_pass"] and
            charge <= WTCalibrationSpec().charge_tolerance_fmol and sustainment_gate(ratio))}


def dynamic_job(job):
    candidate, calcium, solver_label = job
    require_contract_commit()
    inherited = candidate["condition"] == "inherited_baseline"
    path = trajectory_path(candidate["root_id"], calcium, "WT", solver_label) if inherited else (
        OUT / "wt_trajectories" / f"{candidate['candidate_id']}__ca{calcium:.2f}__{solver_label}.json")
    manifest, _ = load_freeze()
    model = allocated_model(manifest, candidate["root_id"], calcium,
        candidate["ae4_capacity_scale"], candidate["nkcc1_capacity_scale"])
    root = read(root_path(candidate))
    if path.exists():
        value = read(path)
        if not inherited:
            assert value["contract_sha256"] == contract()[1]
    else:
        assert not inherited, "Missing inherited WT trajectory; do not silently recompute"
        solver = {s.label: s for s in (PRODUCTION_RADAU, PRODUCTION_BDF)}[solver_label]
        started = time.monotonic()
        trace = asdict(simulate_genotype(model, root["core_state"], transporter="AE4", genotype=WT,
            solver=solver, time_s=physical_time_grid()))
        value = {"contract_sha256": contract()[1], "candidate_id": candidate["candidate_id"],
            "genotype": "WT", "calcium_uM": calcium, "trajectory": trace,
            "whole_cell_parameters_sha256": sha256_object(model.parameters),
            "ae4_parameters_sha256": sha256_object(model.ae4_parameters), "wall_seconds": time.monotonic()-started}
        write_json(path, value)
    assert value["whole_cell_parameters_sha256"] == sha256_object(model.parameters)
    assert value["ae4_parameters_sha256"] == sha256_object(model.ae4_parameters)
    trace = value["trajectory"]
    assert np.array_equal(np.asarray(trace["states"])[:12, 0], root["core_state"])
    row = {**candidate, "calcium_uM": calcium, "solver_label": solver_label,
        "trajectory_path": relative(path), "trajectory_sha256": sha256_file(path),
        "reused_inherited": inherited, **trace_diagnostics(model, trace)}
    print(json.dumps({"candidate": candidate["candidate_id"], "ca": calcium, "solver": solver_label,
        "wt_dynamic_gate_without_bdf": row["wt_dynamic_gate_without_bdf"]}), flush=True)
    return row


def pool(function, jobs, workers):
    if workers == 1:
        return [function(job) for job in jobs]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(function, job) for job in jobs]
        return [f.result() for f in as_completed(futures)]


def wt_stage(workers=4):
    assert not (OUT / "frozen_candidate_manifest.json").exists(), "WT checkpoint is immutable"
    check_inputs()
    c, _ = contract()
    candidates = c["candidate_definitions"]
    rows = pool(rest_job, candidates, workers)
    rows.sort(key=lambda x: x["candidate_id"])
    write_rows(OUT / "wt_rest_candidates.csv", rows)
    valid = [r for r in rows if r["wt_rest_admissible"]]
    dynamic = pool(dynamic_job, [(r, ca, "production_radau") for r in valid for ca in CALCIUM], workers)
    comparisons = []
    groups = sorted({(r["condition"], r["routing_family"]) for r in valid if r["condition"] != "inherited_baseline"})
    for condition, family in groups:
        eligible = [r for r in valid if (r["condition"], r["routing_family"]) == (condition, family)]
        representative = min(eligible, key=lambda r: r["root_id"])
        checks = [dynamic_job((representative, .50, "production_bdf"))]
        expanded = False
        while checks:
            bdf = checks.pop(0)
            radau = next(r for r in dynamic if r["candidate_id"] == bdf["candidate_id"] and r["calcium_uM"] == bdf["calcium_uM"])
            b, a = read(REPO / bdf["trajectory_path"])["trajectory"], read(REPO / radau["trajectory_path"])["trajectory"]
            errors = {key: _relative_difference(a[key], b[key], floor) for key, floor in
                (("states", 1e-10), ("flow_pL_s", 1e-12), ("cumulative_flow_pL", 1e-12))}
            agrees = a["numerical_gate_pass"] and b["numerical_gate_pass"] and max(errors.values()) <= 1e-4
            comparisons.append({"candidate_id": bdf["candidate_id"], "calcium_uM": bdf["calcium_uM"],
                "condition": condition, "routing_family": family, "agrees": agrees, "errors": errors,
                "bdf_path": bdf["trajectory_path"], "bdf_sha256": bdf["trajectory_sha256"]})
            if not agrees and not expanded:
                expanded = True
                jobs = [(r, ca, "production_bdf") for r in eligible for ca in CALCIUM
                    if (r["candidate_id"], ca) != (representative["candidate_id"], .50)]
                checks += pool(dynamic_job, jobs, workers)
    for row in dynamic:
        checks = [x for x in comparisons if x["candidate_id"] == row["candidate_id"] and x["calcium_uM"] == row["calcium_uM"]]
        row["bdf_gate_pass"] = all(x["agrees"] for x in checks)
        row["wt_admissible"] = row["wt_dynamic_gate_without_bdf"] and row["bdf_gate_pass"]
        row["status"] = "WT_ADMISSIBLE" if row["wt_admissible"] else "WT_DYNAMIC_REJECTED"
    for row in rows:
        if not row["wt_rest_admissible"]:
            dynamic += [{**row, "calcium_uM": ca, "solver_label": "production_radau", "wt_admissible": False,
                "status": "NOT_RUN_REST_REJECTED", "trajectory_path": ""} for ca in CALCIUM]
    dynamic.sort(key=lambda x: (x["candidate_id"], x["calcium_uM"]))
    write_rows(OUT / "wt_dynamic_candidates.csv", dynamic)
    write_json(OUT / "wt_solver_confirmations.json", comparisons)
    write_json(OUT / "wt_selection_decisions.json", {"contract_sha256": contract()[1],
        "wt_only": True, "decisions": dynamic, "created_utc": now()})


def freeze():
    check_inputs()
    c, digest = contract()
    path = OUT / "frozen_candidate_manifest.json"
    assert not path.exists(), "WT checkpoint is immutable"
    decisions = read(OUT / "wt_selection_decisions.json")["decisions"]
    expected = {(r["candidate_id"], ca) for r in c["candidate_definitions"] for ca in CALCIUM}
    assert {(r["candidate_id"], r["calcium_uM"]) for r in decisions} == expected
    assert len(decisions) == len(expected) == 90
    hashes = {relative(p): sha256_file(p) for p in OUT.rglob("*") if p.is_file()}
    for d in decisions:
        if d.get("trajectory_path"):
            hashes[d["trajectory_path"]] = sha256_file(REPO / d["trajectory_path"])
    for source in (Path(__file__), Path(__file__).with_name("task16_wt_diagnostics.py"),
            REPO / "tests/test_task16_wt_allocation.py"):
        hashes[relative(source)] = sha256_file(source)
    for source in (REPO / "analysis/16_wt_chloride_allocation").glob("*.md"):
        hashes[relative(source)] = sha256_file(source)
    value = {"manifest_id": "TASK16_WT_ONLY_FROZEN_CANDIDATES_V1", "created_utc": now(),
        "baseline_commit": BASE, "contract_sha256": digest, "wt_only": True,
        "pre_freeze_head": git("rev-parse", "HEAD"), "expected_root_share_calcium_count": 90,
        "all_decisions": decisions, "admissible_candidates": [d for d in decisions if d["wt_admissible"]],
        "hashes": hashes, "frozen_before_new_genotype_evaluation": True,
        "requires_pushed_checkpoint_before_genotype": True, "no_post_genotype_candidate_pruning": True}
    write_json(path, value)
    print(json.dumps({"frozen_combinations": len(value["admissible_candidates"]), "manifest_sha256": sha256_file(path)}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "wt", "freeze"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.stage == "prepare": prepare()
    elif args.stage == "wt": wt_stage(args.workers)
    else: freeze()


if __name__ == "__main__":
    main()
