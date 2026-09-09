"""Frozen Task 14 prediction. No experimental data are read by this module."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import scipy

from .calibration import WTCalibrationSpec
from .genotype_evaluation import (
    DEFAULT_EXPRESSION_GRID, GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
    continue_genotype_expression, genotype_with_expression, simulate_genotype,
    _relative_difference,
)
from .independent import (
    IndependentN1NkccRegulation, IndependentSecretagogueProtocol,
    independent_r1_regulation,
)
from .independent_genotype import continue_independent_genotype_expression
from .independent_native import IndependentNativeSourceDefinition
from .model import WT
from .native_source_panel import NativeRootRecord, build_native_source_model
from .nkcc_stimulation import N1AlgebraicNkcc1, attach_stimulated_nkcc1
from .run_calcium_fast_screen import load_freeze, preferred_member, write_json, write_rows
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU, PRODUCTION_BDF,
    SecretagogueProtocol, StimulusArm, attach_basal_regulation,
    physical_time_grid, sha256_file, sha256_object,
)

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/14_scale_free_genotype_holdout"
ANALYSIS = REPO / "analysis/14_scale_free_genotype_holdout"
BASE = "2f54e7c4f87b6da746d3a8427c7bdca40a77c3de"
BASE_TREE = "338e2891183546734ff9db81611a459136346dce"
INTAKE = "9a1c3fb4e469fa016983b79b274981b1cef21b18"
CALCIUM = (0.10, 0.25, 0.50)
MEMBER = "R1_G125_P21_PROSE_TREFERENCE"
GENOTYPES = ("WT", "AE4", "AE2")
MINUTES = tuple(range(60, 601, 60))


def now():
    return datetime.now(timezone.utc).isoformat()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def contract():
    path = OUT / "prediction_contract.json"
    value = json.loads(path.read_text())
    return value, sha256_file(path)


def stamp(value):
    return {"prediction_contract_sha256": contract()[1], **value}


def write_output(path, value):
    if (OUT / "blind_checkpoint.json").exists():
        raise RuntimeError("The pushed blind checkpoint is immutable")
    write_json(path, stamp(value))


def slug(root_id):
    return root_id.replace(":", "_")


def routing(manifest, root_id):
    return manifest["roots"][root_id]["native_root_object"]["routing_id"]


def build_model(manifest, root_id, calcium):
    if calcium not in CALCIUM:
        raise ValueError("calcium is outside the declared panel")
    payload = manifest["roots"][root_id]
    root = NativeRootRecord(**payload["native_root_object"])
    protocol = SecretagogueProtocol(**{
        **manifest["protocols"]["CCH_IPR"], "arm": StimulusArm.CCH_IPR,
        "stimulated_calcium_uM": calcium,
    })
    member = preferred_member()
    assert member.member_id == MEMBER
    base = build_native_source_model(root, regulatory_model=member.regulatory_model,
                                     stimulus=protocol)
    for value, key in ((base.parameters, "whole_cell_parameters"),
                       (base.ae4_parameters, "ae4_parameters")):
        assert sha256_object(value) == payload[key + "_sha256"]
    return attach_stimulated_nkcc1(base, N1AlgebraicNkcc1(
        fully_activated_multiplier=manifest["stimulated_nkcc1"]["fully_activated_multiplier"],
        resting_calcium_uM=manifest["protocols"]["CCH_IPR"]["resting_calcium_uM"],
        stimulated_calcium_uM=manifest["stimulated_nkcc1"]["effective_calcium_uM"],
    ))


def independent_model(manifest, root_id):
    r = manifest["roots"][root_id]["native_root_object"]
    definition = IndependentNativeSourceDefinition(**{
        key: r[key] for key in ("source_class", "source_scale", "pump_capacity_scale",
        "apical_pump_fraction", "apical_k_fraction", "ae4_cation_fraction",
        "hydraulic_scale", "hydraulic_mode")
    }, cell_other_impermeant_osmoles_fmol=r["other_impermeant_osmoles_fmol"])
    model = definition.build_independent_model(
        regulation=independent_r1_regulation(tau_activation_s=30.0,
            basal_multiplier=1.0, fully_activated_increment=0.25),
        nkcc_regulation=IndependentN1NkccRegulation(fully_activated_multiplier=1.75,
            resting_calcium_uM=0.058, stimulated_calcium_uM=0.10),
        stimulus=IndependentSecretagogueProtocol(resting_calcium_uM=0.058,
            stimulated_calcium_uM=0.10, beta_occupancy_on=1.0),
    )
    assert sha256_object(model.parameters) == manifest["roots"][root_id]["whole_cell_parameters_sha256"]
    return model


def prepare():
    if (OUT / "prediction_contract.json").exists():
        raise FileExistsError("Prediction contract already exists; it is immutable")
    manifest, verification = load_freeze()
    assert git("rev-parse", BASE + "^{tree}") == BASE_TREE
    assert set(git("diff", "--name-only", BASE, INTAKE).splitlines()) == {
        "AGENTS.md", "prompts/14_scale_free_genotype_holdout_validation.md"}
    sources = sorted(set(manifest["equation_source_file_sha256"]) | {
        "src/modern_full_model/genotype_evaluation.py",
        "src/modern_full_model/independent_genotype.py",
        "src/modern_full_model/run_calcium_fast_screen.py",
    })
    source_hashes = {}
    for path in sources:
        assert git("hash-object", path) == git("rev-parse", BASE + ":" + path)
        source_hashes[path] = sha256_file(REPO / path)
    roots = sorted(manifest["roots"])
    families = sorted({routing(manifest, r) for r in roots})
    representatives = {f: min(r for r in roots if routing(manifest, r) == f) for f in families}
    spec = WTCalibrationSpec()
    payload = {
        "contract_id": "TASK14_SCALE_FREE_BLIND_V1", "created_utc": now(),
        "frozen_scientific_commit": BASE, "frozen_scientific_tree": BASE_TREE,
        "branch_intake_commit": INTAKE, "branch_intake_tree": git("rev-parse", INTAKE + "^{tree}"),
        "root_ids": roots, "root_count": 10,
        "null_genotypes": ["AE4", "AE2"], "genotype_expression_objects": {
            g: asdict(WT if g == "WT" else genotype_with_expression(g, 0.0)) for g in GENOTYPES},
        "calcium_uM": CALCIUM, "regulatory_member_id": MEMBER,
        "regulatory_member": asdict(preferred_member()),
        "protocol": manifest["protocols"]["CCH_IPR"],
        "nkcc_normalisation_frozen": manifest["stimulated_nkcc1"],
        "production_solver": asdict(PRODUCTION_RADAU), "confirmation_solver": asdict(PRODUCTION_BDF),
        "solver_relative_tolerance": GENOTYPE_SOLVER_RELATIVE_TOLERANCE,
        "solver_comparison": "maximum symmetric relative error on states, flow and cumulative flow; floors 1e-10, 1e-12, 1e-12",
        "time_grid_s": physical_time_grid().tolist(), "reporting_times_s": MINUTES,
        "early_interval_s": [0, 180], "sustained_interval_s": [180, 600],
        "metrics": ["R_total", "D_total", "R_flow(t)", "R_cum(t)", "R_early", "R_late",
                    "delta_emergence=R_late-R_early", "AE4_effect", "AE2_effect", "effect_contrast",
                    "resting and endpoint null minus WT Cl, pH, Na, K, volume"],
        "comparison_rule": {
            "total": "signed and absolute error in reduction fraction; reported uncertainty only",
            "time_course": "unweighted RMSE over available matching relative time points; no fitted transformations",
            "minute_collections": "compare collection bins with differences of cumulative flow, not endpoint instantaneous flow",
            "early_and_sustained": "use matching integrated intervals only if target resolution permits; otherwise explicitly unavailable",
            "uncertainty": "Do not manufacture an interval from a rounded point estimate or qualitative description",
            "ionic": "separate resting and stimulated values, gland and assay contexts",
        },
        "continuation": {
            "expression_grid": DEFAULT_EXPRESSION_GRID, "coordinates": [b.name for b in spec.coordinate_bounds],
            "coordinate_bounds": [[b.lower, b.upper] for b in spec.coordinate_bounds],
            "local_start_count": 3, "local_start_fraction": 0.001,
            "maximum_connected_step_fraction": 0.25, "max_nfev": 3000,
            "distance_rule": "nearest previous root in maximum coordinate difference divided by coordinate span; tie by residual",
            "fitted_parameter_count": 0, "fixed_OTHER": True, "alternate_roots_retained": True,
            "numerical_tolerances": {k: getattr(spec, k) for k in ("root_scaled_tolerance",
                "omitted_row_raw_tolerance", "charge_tolerance_fmol", "current_tolerance_A",
                "boundary_relative_tolerance", "cluster_relative_tolerance", "jacobian_rank_relative_tolerance")},
        },
        "routing_family_by_root": {r: routing(manifest, r) for r in roots},
        "representatives_lexicographic": representatives,
        "bdf_rule": "Both nulls and WT for each representative at all calcium values; any disagreement expands to all roots in that routing family and calcium",
        "independent_rule": "Both deletions in each representative; also every branch with failure, alternate roots, unsuccessful local solves, or dynamic numerical anomaly",
        "independent_relative_tolerance": 1e-4,
        "numerical_failure_rules": ["failed continuation to exact zero blocks that null trajectory",
            "positive physical core and nonnegative finite outflow required", "exact deleted transporter sources required",
            "all inherited conservation and current tolerances required", "full resting Jacobian rank required",
            "incomplete or failed dynamics are recorded, never replaced by WT initial states",
            "unresolved numerical disagreement blocks valid biological scoring of affected cases"],
        "conservation_tolerances": CONSERVATION_RESIDUAL_TOLERANCES,
        "no_pruning": True, "all_roots_and_calcium_retained_regardless_of_agreement": True,
        "absolute_scale_role": "secondary nonblocking diagnostic; same inherited map for all genotypes; no fit",
        "root_payload_hashes": {r: {k: v for k, v in manifest["roots"][r].items() if k.endswith("_sha256")} for r in roots},
        "equation_source_sha256": source_hashes,
        "task13B_final_manifest_sha256": sha256_file(REPO / "results/13B_modern_full_model/native_dynamic_contract_manifest.json"),
        "task13C_frozen_manifest_sha256": sha256_file(REPO / "results/13C_calcium_fast_screen/frozen_manifest.json"),
        "inherited_verification": verification,
        "environment": {"numpy": np.__version__, "scipy": scipy.__version__,
            "workers": min(9, os.cpu_count() or 1), "BLAS_threads": 1},
        "pre_reveal_state": "No Task 14 genotype result exists. Exact experimental magnitude and time course remain sealed.",
    }
    write_json(OUT / "prediction_contract.json", payload)
    print(json.dumps({"contract_sha256": contract()[1], "representatives": representatives}), flush=True)


def verify_contract():
    c, digest = contract()
    for path, expected in c["equation_source_sha256"].items():
        assert sha256_file(REPO / path) == expected, path
    return c, digest


def continuation_path(root_id, transporter):
    return OUT / "continuation" / f"{slug(root_id)}_{transporter}.json"


def continuation_job(root_id):
    c, digest = verify_contract()
    manifest, _ = load_freeze()
    model = build_model(manifest, root_id, 0.10)
    core = manifest["roots"][root_id]["core_state"]
    results = []
    for g in ("AE4", "AE2"):
        path = continuation_path(root_id, g)
        if path.exists():
            data = json.loads(path.read_text()); assert data["prediction_contract_sha256"] == digest
            results.append({"root": root_id, "genotype": g, "cached": True}); continue
        start = time.monotonic()
        report = continue_genotype_expression(model, core, transporter=g)
        data = asdict(report)
        data.update(root_id=root_id, routing_family=routing(manifest, root_id))
        data["selected_conservation"] = []
        for step in report.steps:
            if step.selected_root is None:
                continue
            state = attach_basal_regulation(model, step.selected_root.core_state)
            ev = model.evaluate(0.0, state, genotype=genotype_with_expression(g, step.expression))
            residuals = dict(ev.diagnostics.conservation_residuals)
            ratio = max(abs(v) / CONSERVATION_RESIDUAL_TOLERANCES[k] for k, v in residuals.items())
            data["selected_conservation"].append({"expression": step.expression,
                "residuals": residuals, "max_dimensionless_ratio": ratio})
        anomaly = (not report.completed_to_exact_zero or any(s.alternate_root_count or
            any(not a.optimizer_success or not a.converged_root for a in s.attempts) for s in report.steps))
        data["continuation_anomaly"] = anomaly
        data["primary_wall_seconds"] = time.monotonic() - start
        if anomaly or root_id in c["representatives_lexicographic"].values():
            independent = continue_independent_genotype_expression(independent_model(manifest, root_id), core, transporter=g)
            data["independent_audit"] = asdict(independent)
            differences = []
            for primary, audit in zip(report.steps, independent.steps):
                if primary.selected_root is not None and audit.selected_root is not None:
                    differences.append(_relative_difference(primary.selected_root.core_state, audit.selected_root.core_state, 1e-10))
            data["independent_max_state_relative_difference"] = max(differences, default=None)
            data["independent_agrees"] = bool(report.completed_to_exact_zero == independent.completed_to_exact_zero
                and differences and max(differences) <= c["independent_relative_tolerance"])
        data["total_wall_seconds"] = time.monotonic() - start
        write_output(path, data)
        results.append({"root": root_id, "genotype": g, "complete": report.completed_to_exact_zero,
            "failure": report.branch_failure, "seconds": data["total_wall_seconds"]})
        print(json.dumps(results[-1]), flush=True)
    return results


def null_core(root_id, genotype):
    data = json.loads(continuation_path(root_id, genotype).read_text())
    assert data["prediction_contract_sha256"] == contract()[1]
    if not data["completed_to_exact_zero"]:
        return None
    if any(r["max_dimensionless_ratio"] > 1.0 for r in data["selected_conservation"]):
        return None
    step = data["steps"][-1]
    return next(r["core_state"] for r in step["roots"] if r["root_id"] == step["selected_root_id"])


def trajectory_path(root_id, calcium, genotype, label):
    return OUT / "trajectories" / f"{slug(root_id)}_ca{calcium:.2f}_{genotype}_{label}.json"


def trajectory_job(job):
    root_id, calcium, genotype, label = job
    c, digest = verify_contract()
    path = trajectory_path(*job)
    if path.exists():
        value = json.loads(path.read_text()); assert value["prediction_contract_sha256"] == digest
        return {"job": job, "cached": True}
    manifest, _ = load_freeze()
    core = manifest["roots"][root_id]["core_state"] if genotype == "WT" else null_core(root_id, genotype)
    base = {"root_id": root_id, "routing_family": routing(manifest, root_id), "calcium_uM": calcium,
            "genotype": genotype, "solver_label": label}
    if core is None:
        write_output(path, {**base, "status": "CONTINUATION_FAILED", "trajectory": None})
        return {**base, "status": "CONTINUATION_FAILED"}
    model = build_model(manifest, root_id, calcium)
    deletion = WT if genotype == "WT" else genotype_with_expression(genotype, 0.0)
    solver = {s.label: s for s in (PRODUCTION_RADAU, PRODUCTION_BDF)}[label]
    start = time.monotonic()
    trajectory = simulate_genotype(model, core, transporter="AE4" if genotype == "WT" else genotype,
        genotype=deletion, solver=solver, time_s=np.asarray(c["time_grid_s"]))
    t = asdict(trajectory)
    # Charge of the physical state is separate from its rate in the inherited accounting.
    charges = trajectory.states[0] + trajectory.states[1] - trajectory.states[2] - trajectory.states[4]
    charges -= model.parameters.geometry.fixed_cell_anion_equivalents_fmol
    lumen_charge = trajectory.states[6] + trajectory.states[7] - trajectory.states[8] - trajectory.states[10]
    additional = {"max_abs_cell_charge_fmol": float(np.max(np.abs(charges))),
                  "max_abs_lumen_charge_fmol": float(np.max(np.abs(lumen_charge)))}
    for which, i in (("rest", 0), ("endpoint", -1)):
        ev = model.evaluate(float(trajectory.time_s[i]), trajectory.states[:, i], genotype=deletion)
        additional[which + "_v_apical_mV"] = 1000 * ev.diagnostics.membranes.v_apical_V
        additional[which + "_v_basolateral_mV"] = 1000 * ev.diagnostics.membranes.v_basolateral_V
        additional[which + "_fluxes"] = {"ae4": asdict(ev.diagnostics.ae4),
            "homeostasis": asdict(ev.diagnostics.homeostasis)}
    data = {**base, "status": "COMPLETE" if trajectory.numerical_gate_pass else "NUMERICAL_FAILURE",
        "trajectory": t, "additional_diagnostics": additional, "wall_seconds": time.monotonic() - start,
        "whole_cell_parameters_sha256": sha256_object(model.parameters),
        "ae4_parameters_sha256": sha256_object(model.ae4_parameters)}
    write_output(path, data)
    return {**base, "status": data["status"], "seconds": data["wall_seconds"]}


def load_trajectory(root, ca, genotype, solver="production_radau"):
    value = json.loads(trajectory_path(root, ca, genotype, solver).read_text())
    assert value["prediction_contract_sha256"] == contract()[1]
    return value


def ratios(wt, null):
    if wt["solver"] != null["solver"] or wt["time_s"] != null["time_s"]:
        raise ValueError("WT and null must use the same solver and time grid")
    times = np.asarray(wt["time_s"])
    wc, nc = np.asarray(wt["cumulative_flow_pL"]), np.asarray(null["cumulative_flow_pL"])
    wq, nq = np.asarray(wt["flow_pL_s"]), np.asarray(null["flow_pL_s"])
    at = lambda values, t: float(np.interp(t, times, values))
    if wc[-1] <= 0 or at(wc, 180) <= 0 or wc[-1] - at(wc, 180) <= 0:
        raise ValueError("nonpositive WT integrated denominator")
    row = {"R_total": float(nc[-1] / wc[-1]), "R_early": at(nc, 180) / at(wc, 180),
           "R_late": (nc[-1] - at(nc, 180)) / (wc[-1] - at(wc, 180))}
    row.update(D_total=1 - row["R_total"], delta_emergence=row["R_late"] - row["R_early"])
    for minute in MINUTES:
        row[f"R_flow_{minute}s"] = at(nq, minute) / at(wq, minute)
        row[f"R_cum_{minute}s"] = at(nc, minute) / at(wc, minute)
        row[f"R_bin_{minute}s"] = (at(nc, minute) - at(nc, minute - 60)) / (at(wc, minute) - at(wc, minute - 60))
    for field in ("cell_cl_mM", "cell_ph", "cell_na_mM", "cell_k_mM", "cell_volume_pL"):
        for label, i in (("rest", 0), ("endpoint", -1)):
            row[f"{label}_delta_{field}"] = null[field][i] - wt[field][i]
    return row


def solver_crosscheck(root, ca, genotype):
    a = load_trajectory(root, ca, genotype, "production_radau")
    b = load_trajectory(root, ca, genotype, "production_bdf")
    row = {"root_id": root, "calcium_uM": ca, "genotype": genotype}
    if a["status"] == b["status"] == "CONTINUATION_FAILED":
        return {**row, "pass": None, "reason": "No exact null root; neither solver trajectory is licensed"}
    if a["status"] != "COMPLETE" or b["status"] != "COMPLETE":
        return {**row, "pass": False, "reason": "trajectory incomplete or numerical gate failed"}
    x, y = a["trajectory"], b["trajectory"]
    diffs = {k: _relative_difference(x[k], y[k], floor) for k, floor in
             (("states", 1e-10), ("flow_pL_s", 1e-12), ("cumulative_flow_pL", 1e-12))}
    return {**row, **diffs, "pass": max(diffs.values()) <= GENOTYPE_SOLVER_RELATIVE_TOLERANCE}


def parallel(function, jobs, workers):
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(function, job): job for job in jobs}
        for future in as_completed(futures):
            print(json.dumps(future.result()), flush=True)


def run_dynamics(workers):
    c, _ = verify_contract()
    jobs = [(r, ca, g, "production_radau") for ca in CALCIUM for g in GENOTYPES for r in c["root_ids"]]
    parallel(trajectory_job, jobs, workers)
    representatives = list(c["representatives_lexicographic"].values())
    jobs = [(r, ca, g, "production_bdf") for r in representatives for ca in CALCIUM for g in GENOTYPES]
    parallel(trajectory_job, jobs, workers)
    checks = [solver_crosscheck(*job[:3]) for job in jobs]
    affected = {(c["routing_family_by_root"][row["root_id"]], row["calcium_uM"]) for row in checks if row["pass"] is False}
    extra = [(r, ca, g, "production_bdf") for r in c["root_ids"] for ca in CALCIUM for g in GENOTYPES
             if (c["routing_family_by_root"][r], ca) in affected and r not in representatives]
    if extra:
        parallel(trajectory_job, extra, workers)
        checks.extend(solver_crosscheck(*job[:3]) for job in extra)
    write_output(OUT / "solver_confirmation.json", {"checks": checks,
        "expanded_family_calcium": sorted(affected),
        "all_available_pass": all(row["pass"] is not False for row in checks),
        "unavailable": sum(row["pass"] is None for row in checks)})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("prepare", "continuation", "dynamics"))
    parser.add_argument("--workers", type=int, default=min(9, os.cpu_count() or 1))
    args = parser.parse_args()
    if args.phase == "prepare": prepare()
    elif args.phase == "continuation": parallel(continuation_job, contract()[0]["root_ids"], args.workers)
    else: run_dynamics(args.workers)


if __name__ == "__main__":
    main()
