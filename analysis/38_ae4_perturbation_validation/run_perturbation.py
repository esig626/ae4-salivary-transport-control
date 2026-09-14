"""Task 38: only the two authorized, ordered frozen AE4 predictions.

Invoke with ae4_5pct first, then ae4_null after the first case passes.
No WT integration, stationary solve, search, or scientific parameter change.
Task 37's model/state factory, diagnostic gates, solver and quadrature are reused.
Existing attempts cannot be repeated, including after an interrupted process.
"""
from pathlib import Path
from dataclasses import asdict, replace
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
import traceback

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis/37_wt_nbc_validation"))
from validation_common import (models_and_state, diagnose, write_json, WT,
    CHECKPOINT, CARRIER, Q_REST, CONSERVATION_RESIDUAL_TOLERANCES, sha256_object)
import numpy as np
import scipy
from scipy.integrate import Radau, BDF, cumulative_trapezoid
from threadpoolctl import threadpool_info
from modern_full_model.validation import PRODUCTION_RADAU, PRODUCTION_BDF

OUT = ROOT / "results/38_ae4_perturbation_validation"
HERE = Path(__file__).resolve().parent
CASES = {"ae4_5pct": 0.05, "ae4_null": 0.00}
FAILURES = {"ae4_5pct": "AE4_5PCT_PERTURBATION_FAILED",
            "ae4_null": "AE4_NULL_PERTURBATION_FAILED"}
QUADRATURE = "Trapezoidal on integer-second output plus t=0 and 1e-6; identical to Task 37."
FLUXES = (
    ("positive_nkcc1_cl_loading", "nkcc1_cl_inward_fmol_s", True, 1., "fmol"),
    ("positive_ae4_cl_loading", "ae4_cl_inward_fmol_s", True, 1., "fmol"),
    ("positive_ae2_cl_loading", "ae2_cl_inward_signed_fmol_s", True, 1., "fmol"),
    ("signed_ae2_cl_flux", "ae2_cl_inward_signed_fmol_s", False, 1., "fmol"),
    ("nbc_cycles", "nbc_cycle_inward_fmol_s", False, 1., "fmol"),
    ("nbc_bicarbonate_equivalent_influx", "nbc_cycle_inward_fmol_s", False, 2., "fmol"),
    ("nhe1_flux", "nhe1_inward_fmol_s", False, 1., "fmol"),
    ("cacc_cl_export", "cacc_cl_export_fmol_s", False, 1., "fmol"),
    ("paracellular_cl_return", "paracellular_cl_return_fmol_s", False, 1., "fmol"),
    ("water_outflow", "q_out_pL_s", False, 1., "pL"),
)


class ScientificStop(Exception):
    pass


class ComputeStop(Exception):
    pass


def read_json(path):
    return json.loads(path.read_text())


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_sources():
    manifest = read_json(OUT / "frozen_inputs.json")
    branch = subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip()
    assert branch == manifest["branch"]
    for item in manifest["prepared_files"]:
        raw = (ROOT / item["path"]).read_bytes()
        digest = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        assert digest == item["git_blob_sha"], item["path"]
    return manifest


def run_case(case):
    start = time.perf_counter()
    manifest = verify_sources()
    wt = read_json(ROOT / "results/37_wt_nbc_validation/wt_verification.json")
    assert wt["complete_600_s"] and wt["first_failure"] is None
    assert wt["status"].startswith("WT_NBC_VALIDATION_PASS")
    if case == "ae4_null":
        first = read_json(OUT / "ae4_5pct_verification.json")
        assert first["status"] == "PASS" and first["complete_600_s"]
        assert first["first_failure"] is None
    budget_path = OUT / "budget.json"
    if budget_path.exists():
        budget = read_json(budget_path)
        assert not budget.get("no_further_scientific_execution", False)
    else:
        assert case == "ae4_5pct"
        budget = {"intended_perturbation_integrations": 2, "started_cases": [],
            "completed_cases": [], "integration_attempts": 0, "numerical_retries": 0,
            "numerical_execution_seconds": 0.0, "maximum_numerical_execution_seconds": 600,
            "wt_integrations": 0, "stationary_solves": 0, "optimisation_calls": 0,
            "parameter_sweeps": 0, "additional_scientific_parameter_evaluations": 0,
            "scientific_workers": 1, "blas_threads": 1}
    assert case not in budget["started_cases"]
    assert budget["started_cases"] == ([] if case == "ae4_5pct" else ["ae4_5pct"])
    assert budget["numerical_execution_seconds"] < 600
    # Durable ledger entry is written before any model evaluation/integration.
    budget["started_cases"].append(case)
    budget["status"] = "RUNNING_" + case
    write_json(budget_path, budget)

    saved, _, _, model, y0 = models_and_state()  # State loading; no root solve/evaluation.
    assert np.array_equal(y0, wt["initial_state_vector"])
    assert asdict(model.stimulus) == wt["stimulus"]
    assert asdict(PRODUCTION_RADAU) == wt["integration_attempts"][0]["solver"]
    assert dict(CONSERVATION_RESIDUAL_TOLERANCES) == wt["conservation_tolerances"]
    assert {"numpy": np.__version__, "scipy": scipy.__version__} == wt["software_versions"]
    genotype = replace(WT, ae4_expression=CASES[case])
    genotype_diff = {key: value for key, value in asdict(genotype).items()
                     if value != asdict(WT)[key]}
    assert genotype_diff == {"ae4_expression": CASES[case]}
    hashes = {"whole_cell": sha256_object(model.parameters),
              "ae4": sha256_object(model.ae4_parameters),
              "nbc": sha256_object(model.nbc_parameters)}
    assert hashes == wt["parameter_hashes"]
    assert asdict(model.nbc_parameters) == wt["nbc_parameters"]
    pools = threadpool_info()
    assert pools and all(p["num_threads"] == 1 for p in pools)
    started_utc = datetime.now(timezone.utc).isoformat()
    rows, grid, attempts = [], [], []
    minima, maxima, residual_max, ratio_max = {}, {}, {}, {}
    first_failure, previous_pass, raw_last, error = None, None, None, None
    status, solver, complete, accepted_steps = None, None, False, 0

    def check_time(*_):
        if budget["numerical_execution_seconds"] + time.perf_counter() - start >= 600:
            raise ComputeStop("Task 38 total numerical execution limit: 600 seconds")

    def alarm(*_):
        raise ComputeStop("Task 38 total numerical execution limit: 600 seconds")

    def inspect(t, y, kind):
        nonlocal first_failure, previous_pass, raw_last
        check_time()
        evaluation = model.evaluate(float(t), y, genotype=genotype)
        # Pass the actual perturbation evaluation so no default WT evaluation occurs.
        row, failures, ratios = diagnose(model, t, y, evaluation)
        rows.append(row)
        raw_last = list(map(float, y))
        if t == 0.0 or t == 1e-6 or float(t).is_integer():
            if not grid or grid[-1]["time_s"] != t:
                grid.append(row)
        for key, value in row.items():
            if key != "time_s":
                minima[key] = min(minima.get(key, value), value)
                maxima[key] = max(maxima.get(key, value), value)
        for key in CONSERVATION_RESIDUAL_TOLERANCES:
            residual_max[key] = max(residual_max.get(key, 0.0),
                abs(float(evaluation.diagnostics.conservation_residuals[key])))
            ratio_max[key] = max(ratio_max.get(key, 0.0), ratios[key])
        if failures:
            first_failure = {"time_s": float(t), "previous_checked_pass_time_s": previous_pass,
                "sample_kind": kind, "failed_gates": failures, "observables": row,
                "state_vector": raw_last, "rhs": list(map(float, evaluation.rhs)),
                "conservation_residuals": dict(evaluation.diagnostics.conservation_residuals)}
            raise ScientificStop(FAILURES[case])
        previous_pass = float(t)

    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, max(0.001,
        600 - budget["numerical_execution_seconds"] - (time.perf_counter() - start)))
    try:
        inspect(0.0, y0, "saved_WT_rest_at_exact_onset_with_perturbed_expression")
        inspect(1e-6, y0, "inherited_post_onset_start")
        for specification in (PRODUCTION_RADAU, PRODUCTION_BDF):
            if attempts:
                assert attempts[-1]["solver_status"] == "failed" and first_failure is None
                budget["numerical_retries"] += 1
                rows.clear(); grid.clear()
                minima.clear(); maxima.clear(); residual_max.clear(); ratio_max.clear()
                previous_pass = None
                inspect(0.0, y0, "saved_WT_rest_at_exact_onset_with_perturbed_expression")
                inspect(1e-6, y0, "inherited_post_onset_start")
            budget["integration_attempts"] += 1
            write_json(budget_path, budget)
            attempt = {"attempt": len(attempts) + 1, "solver": asdict(specification),
                "initial_state_sha256": sha256_object(list(map(float, y0)))}
            attempts.append(attempt)
            attempt_start = time.perf_counter()
            solver_class = Radau if specification.method == "Radau" else BDF
            solver = solver_class(lambda t, y: model.rhs(t, y, genotype=genotype),
                1e-6, y0.copy(), 600.0, rtol=specification.rtol,
                atol=specification.atol_vector(model.state_names), max_step=specification.max_step_s)
            next_time = 1.0
            while solver.status == "running":
                check_time()
                message = solver.step()
                if solver.status == "failed":
                    attempt["solver_message"] = str(message)
                    break
                accepted_steps += 1
                dense = solver.dense_output()
                while next_time <= solver.t:
                    inspect(next_time, dense(next_time), "one_second_dense_output")
                    next_time += 1.0
                if not rows or rows[-1]["time_s"] != float(solver.t):
                    inspect(float(solver.t), solver.y, "accepted_solver_endpoint")
            attempt.update({"solver_status": solver.status, "last_solver_time_s": float(solver.t),
                "nfev": int(solver.nfev), "njev": int(solver.njev), "nlu": int(solver.nlu),
                "wall_seconds": time.perf_counter() - attempt_start})
            if solver.status == "finished":
                complete = bool(solver.t == 600.0)
                break
        status = "PASS" if complete else "TASK38_NUMERICAL_SOLVER_FAILED"
    except ScientificStop:
        status = FAILURES[case]
    except ComputeStop as exc:
        status, error = "TASK38_COMPUTE_BUDGET_STOP", str(exc)
    except Exception:
        status, error = "TASK38_NUMERICAL_EXECUTION_ABORTED", traceback.format_exc()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)

    if attempts and solver is not None:
        attempts[-1].update({"solver_status": solver.status,
            "stopped_by_gate": first_failure is not None, "last_solver_time_s": float(solver.t),
            "nfev": int(solver.nfev), "njev": int(solver.njev), "nlu": int(solver.nlu)})
    assert hashes == {"whole_cell": sha256_object(model.parameters),
                      "ae4": sha256_object(model.ae4_parameters),
                      "nbc": sha256_object(model.nbc_parameters)}
    verify_sources()
    secretion, partition, integrals = None, None, []
    if complete and status == "PASS":
        times = np.array([row["time_s"] for row in grid])
        assert np.array_equal(times, np.r_[0., 1e-6, np.arange(1., 601.)])
        q = np.array([row["q_out_pL_s"] for row in grid])
        cumulative = cumulative_trapezoid(q, times, initial=0.)
        window = times >= 60
        for row, value in zip(grid, cumulative):
            row["cumulative_outflow_0_t_pL"] = float(value)
        # Use the same aggregate operation as Task 37 for the headline ratio.
        total = float(np.trapezoid(q, times))
        mean = float(np.trapezoid(q[window], times[window]) / 540)
        secretion = {"cumulative_0_600_pL": total, "mean_60_600_pL_s": mean,
            "endpoint_pL_s": float(q[-1]), "frozen_wt_cumulative_0_600_pL": wt["secretion_gate"]["cumulative_0_600_pL"],
            "ratio_to_frozen_wt": total / wt["secretion_gate"]["cumulative_0_600_pL"],
            "quadrature": QUADRATURE,
            "wt_secretion_thresholds_for_context_only": {
                "mean_60_600_gt_1p10_rest": mean > 1.10 * Q_REST,
                "cumulative_0_600_gt_1p10_rest": total > 1.10 * 600 * Q_REST,
                "endpoint_at_least_rest": float(q[-1]) >= Q_REST}}
        values = {}
        for name, field, positive, scale, unit in FLUXES:
            v = np.array([row[field] for row in grid])[window] * scale
            if positive:
                v = np.maximum(v, 0.)
            values[name] = float(np.trapezoid(v, times[window]))
            integrals.append({"quantity": name, "value": values[name], "unit": unit,
                "window_start_s": 60, "window_end_s": 600})
        pool = sum(values["positive_" + key + "_cl_loading"] for key in ("nkcc1", "ae4", "ae2"))
        partition = {"window_s": [60, 600], "positive_loading_pool_fmol": pool,
            "positive_shares": {key: values["positive_" + key + "_cl_loading"] / pool
                                for key in ("nkcc1", "ae4", "ae2")},
            "definition": "Integrate max(J_Cl, 0) separately for NKCC1, AE4 and AE2; sum these three only."}
    readouts = [row for row in grid if row["time_s"] % 60 == 0]
    if first_failure:
        readouts.append(first_failure["observables"])
    write_csv(OUT / (case + "_timeseries.csv"), readouts)
    write_csv(OUT / (case + "_integrated_fluxes.csv"), integrals)
    verification = {"status": status, "case": case, "ae4_expression": CASES[case],
        "prepared_head": manifest["prepared_head"], "frozen_wt_commit": manifest["frozen_wt_commit"],
        "branch": manifest["branch"], "started_utc": started_utc,
        "complete_600_s": complete, "first_failure": first_failure, "solver_exception": error,
        "integration_attempts": attempts, "accepted_solver_steps": accepted_steps,
        "checked_state_count": len(rows), "last_checked_time_s": rows[-1]["time_s"] if rows else None,
        "monitoring": wt["monitoring"], "onset_convention": wt["onset_convention"],
        "state_names": list(model.state_names), "initial_state_vector": list(map(float, y0)),
        "initial_state_sha256": sha256_object(list(map(float, y0))),
        "initial_state_bitwise_identical_to_task37": True, "last_checked_state_vector": raw_last,
        "checkpoint_path": str(CHECKPOINT.relative_to(ROOT)), "core_state_sha256": saved["core_state_sha256"],
        "genotype": asdict(genotype), "only_genotype_difference": genotype_diff,
        "stimulus": asdict(model.stimulus), "parameter_hashes": hashes,
        "scientific_parameters_unchanged": True, "all_prepared_blobs_unchanged": True,
        "nbc_parameters": asdict(model.nbc_parameters), "nhe1_carrier_amount_fmol": CARRIER,
        "minima": minima, "maxima": maxima, "max_abs_conservation_residuals": residual_max,
        "conservation_tolerances": dict(CONSERVATION_RESIDUAL_TOLERANCES),
        "max_conservation_ratios": ratio_max, "secretion": secretion, "chloride_partition": partition,
        "software_versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "thread_pools": [{k: v for k, v in p.items() if k in ("internal_api", "prefix", "num_threads", "version")} for p in pools],
        "source_fields_note": wt["source_fields_note"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    write_json(OUT / (case + "_verification.json"), verification)
    if complete:
        budget["completed_cases"].append(case)
    budget["numerical_execution_seconds"] += time.perf_counter() - start
    budget["status"] = ("AE4_PERTURBATION_VALIDATION_COMPLETE" if case == "ae4_null" and status == "PASS"
                        else "AE4_5PCT_PASS_AWAITING_NULL" if status == "PASS" else status)
    budget["no_further_scientific_execution"] = case == "ae4_null" or status != "PASS"
    budget["limits_respected"] = (budget["integration_attempts"] <= 4
        and len(budget["started_cases"]) <= 2 and len(attempts) <= 2
        and budget["numerical_execution_seconds"] <= 600)
    write_json(budget_path, budget)
    (HERE / (case + "_execution.log")).write_text(json.dumps({"case": case, "status": status,
        "attempts": attempts, "first_failure": first_failure, "error": error, "budget": budget}, indent=2) + "\n")
    print(json.dumps({"case": case, "status": status, "secretion": secretion,
        "chloride_partition": partition, "first_failure": first_failure,
        "max_conservation_ratio": max(ratio_max.values(), default=None), "budget": budget}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=list(CASES))
    run_case(parser.parse_args().case)
