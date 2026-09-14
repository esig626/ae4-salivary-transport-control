"""One continuous production WT integration, stopped at the first failed gate.

Low-level scipy Radau/BDF stepping permits checks before the next accepted
step. This does not restart the ODE at observation times. The exact production
method, rtol, coordinate-specific atol, max step, and inherited 1-us post-onset
start are retained. A fresh invocation refuses to repeat an existing attempt.
"""
import time
START = time.perf_counter()
from validation_common import *
import csv
import signal
import traceback
import scipy
from scipy.integrate import Radau, BDF
from threadpoolctl import threadpool_info
from modern_full_model.validation import PRODUCTION_RADAU, PRODUCTION_BDF

PHYSIOLOGY_STOP = "WT_NBC_TRAJECTORY_PHYSIOLOGY_FAILED"
SECRETION_STOP = "WT_NBC_SECRETION_GATE_FAILED"


class GateStop(Exception):
    pass


class TimeStop(Exception):
    pass


def write_csv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    budget = json.loads((RESULTS / "budget.json").read_text())
    assert budget["focused_tests_pass"] and budget["rest_nesting_pass"]
    assert budget["wt_integration_attempts"] == 0, "Existing WT attempt must not be repeated"
    assert budget["numerical_execution_seconds"] < 600
    with (RESULTS / ".wt_execution_started").open("x") as f:
        f.write("Task 37 single intended WT integration; do not rerun.\n")
    saved, _, _, model, y0 = models_and_state()
    thread_info = threadpool_info()
    assert all(item["num_threads"] == 1 for item in thread_info)
    initial_parameter_hash = sha256_object(model.parameters)
    initial_ae4_hash = sha256_object(model.ae4_parameters)
    initial_nbc_hash = sha256_object(model.nbc_parameters)
    started_utc = __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()
    logs = []
    attempts = []
    accepted_steps = 0
    rows = []
    grid_rows = []
    minima, maxima, residual_max, ratio_max = {}, {}, {}, {}
    previous_pass_time = None
    first_failure = None
    raw_last = None
    solver = None
    status = None
    complete = False
    solver_exception = None

    def check_time(*_):
        if budget["numerical_execution_seconds"] + time.perf_counter() - START >= 600:
            raise TimeStop("Task 37 600-second numerical execution limit")

    def alarm(*_):
        raise TimeStop("Task 37 600-second numerical execution limit")

    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, max(0.001, 600 - budget["numerical_execution_seconds"] - (time.perf_counter() - START)))

    def inspect(t, y, *, sample_kind):
        nonlocal previous_pass_time, first_failure, raw_last
        check_time()
        e = model.evaluate(float(t), y, genotype=WT)
        row, failures, ratios = diagnose(model, t, y, e)
        rows.append(row)
        raw_last = list(map(float, y))
        if t == 0.0 or t == 1e-6 or float(t).is_integer():
            if not grid_rows or grid_rows[-1]["time_s"] != t:
                grid_rows.append(row)
        for key, val in row.items():
            if key != "time_s":
                minima[key] = min(minima.get(key, val), val)
                maxima[key] = max(maxima.get(key, val), val)
        for key, tol in CONSERVATION_RESIDUAL_TOLERANCES.items():
            val = abs(float(e.diagnostics.conservation_residuals[key]))
            residual_max[key] = max(residual_max.get(key, 0.0), val)
            ratio_max[key] = max(ratio_max.get(key, 0.0), ratios[key])
        if failures:
            first_failure = {"time_s": float(t), "previous_checked_pass_time_s": previous_pass_time,
                "sample_kind": sample_kind, "failed_gates": failures, "observables": row,
                "state_vector": raw_last, "rhs": list(map(float, e.rhs)),
                "conservation_residuals": dict(e.diagnostics.conservation_residuals)}
            raise GateStop(PHYSIOLOGY_STOP)
        previous_pass_time = float(t)
        return row

    try:
        inspect(0.0, y0, sample_kind="saved_rest_at_exact_onset")
        inspect(1e-6, y0, sample_kind="inherited_post_onset_start")
        for specification in (PRODUCTION_RADAU, PRODUCTION_BDF):
            # BDF is reached only after Radau's explicit ODE status 'failed'.
            budget["wt_integration_attempts"] += 1
            if budget["wt_integration_attempts"] == 2:
                budget["numerical_retries"] = 1
                assert attempts[-1]["solver_status"] == "failed" and first_failure is None
                rows.clear(); grid_rows.clear()
                minima.clear(); maxima.clear(); residual_max.clear(); ratio_max.clear()
                previous_pass_time = None
                inspect(0.0, y0, sample_kind="saved_rest_at_exact_onset")
                inspect(1e-6, y0, sample_kind="inherited_post_onset_start")
            write_json(RESULTS / "budget.json", budget)
            attempt = {"attempt": budget["wt_integration_attempts"], "solver": asdict(specification),
                       "initial_state_sha256": sha256_object(list(map(float, y0)))}
            attempts.append(attempt)
            attempt_start = time.perf_counter()
            solver_class = Radau if specification.method == "Radau" else BDF
            solver = solver_class(lambda t, y: model.rhs(t, y, genotype=WT), 1e-6, y0.copy(), 600.0,
                rtol=specification.rtol, atol=specification.atol_vector(model.state_names),
                max_step=specification.max_step_s)
            next_grid_time = 1.0
            while solver.status == "running":
                check_time()
                message = solver.step()
                if solver.status == "failed":
                    attempt["solver_message"] = str(message)
                    break
                accepted_steps += 1
                dense = solver.dense_output()
                # Check integer-second dense output in chronological order and
                # every accepted endpoint. No gate is applied to trial stages.
                while next_grid_time <= solver.t:
                    inspect(next_grid_time, dense(next_grid_time), sample_kind="one_second_dense_output")
                    next_grid_time += 1.0
                if not rows or rows[-1]["time_s"] != float(solver.t):
                    inspect(float(solver.t), solver.y, sample_kind="accepted_solver_endpoint")
            attempt.update({"solver_status": solver.status, "last_solver_time_s": float(solver.t),
                "nfev": int(solver.nfev), "njev": int(solver.njev), "nlu": int(solver.nlu),
                "wall_seconds": time.perf_counter() - attempt_start})
            if solver.status == "finished":
                complete = bool(solver.t == 600.0)
                break
            logs.append("Numerical ODE solver failure: " + str(attempt.get("solver_message")))
            # No scientific-input changes and no retry after a physical failure.
        if not complete:
            status = "TASK37_NUMERICAL_SOLVER_FAILED"
    except GateStop:
        status = PHYSIOLOGY_STOP
    except TimeStop as exc:
        status = "TASK37_COMPUTE_BUDGET_STOP"
        solver_exception = str(exc)
    except Exception as exc:
        # An exception is not presumed to justify the numerical retry. Abort.
        status = "TASK37_NUMERICAL_EXECUTION_ABORTED"
        solver_exception = type(exc).__name__ + ": " + str(exc)
        logs.append(traceback.format_exc())
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)

    if attempts and solver is not None:
        attempts[-1].update({"solver_status": solver.status,
            "stopped_by_gate": status == PHYSIOLOGY_STOP,
            "last_solver_time_s": float(solver.t), "nfev": int(solver.nfev),
            "njev": int(solver.njev), "nlu": int(solver.nlu)})
    unchanged = (initial_parameter_hash == sha256_object(model.parameters)
        and initial_ae4_hash == sha256_object(model.ae4_parameters)
        and initial_nbc_hash == sha256_object(model.nbc_parameters))
    assert unchanged
    secretion = {"status": "NOT_EVALUATED_EARLY_STOP", "q_out_rest_pL_s": Q_REST}
    partition = {"status": "NOT_EVALUATED_EARLY_STOP", "reason": status,
                 "requested_window_s": [60, 600], "validation_only": True}
    integrated_rows = []
    if complete and status is None:
        times = np.asarray([row["time_s"] for row in grid_rows])
        q = np.asarray([row["q_out_pL_s"] for row in grid_rows])
        window = times >= 60
        window_times = times[window]
        mean_q = float(np.trapezoid(q[window], window_times) / 540)
        cumulative_q = float(np.trapezoid(q, times))
        endpoint_q = float(q[-1])
        checks = {"mean_60_600_gt_1p10_rest": mean_q > 1.10 * Q_REST,
                  "cumulative_0_600_gt_1p10_rest": cumulative_q > 1.10 * 600 * Q_REST,
                  "endpoint_at_least_rest": endpoint_q >= Q_REST}
        secretion = {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks,
            "q_out_rest_pL_s": Q_REST, "mean_60_600_pL_s": mean_q,
            "mean_threshold_pL_s": 1.10 * Q_REST, "cumulative_0_600_pL": cumulative_q,
            "cumulative_threshold_pL": 1.10 * 600 * Q_REST, "endpoint_pL_s": endpoint_q,
            "quadrature": "trapezoidal, inherited one-second output grid plus t=0 and 1e-6"}
        if not all(checks.values()):
            status = SECRETION_STOP
            partition["reason"] = status
        else:
            specs = (
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
            values = {}
            for name, field, positive, scale, unit in specs:
                v = np.asarray([row[field] for row in grid_rows])[window] * scale
                if positive:
                    v = np.maximum(v, 0.)
                values[name] = float(np.trapezoid(v, window_times))
                integrated_rows.append({"quantity": name, "value": values[name], "unit": unit,
                    "window_start_s": 60, "window_end_s": 600, "status": "COMPUTED"})
            pool = sum(values["positive_" + transporter + "_cl_loading"] for transporter in ("nkcc1", "ae4", "ae2"))
            shares = {name: values["positive_" + name + "_cl_loading"] / pool for name in ("nkcc1", "ae4", "ae2")}
            in_context = 0.65 <= shares["nkcc1"] <= 0.75
            status = "WT_NBC_VALIDATION_PASS" if in_context else "WT_NBC_VALIDATION_PASS_PARTITION_OUTSIDE_CONTEXT"
            partition = {"status": "COMPUTED", "window_s": [60,600], "positive_loading_pool_fmol": pool,
                "positive_shares": shares, "nkcc1_context_band": [0.65,0.75],
                "nkcc1_within_context": in_context, "validation_only": True,
                "rest_ae4_share": grid_rows[0]["ae4_cl_inward_fmol_s"] / (
                    grid_rows[0]["nkcc1_cl_inward_fmol_s"] + grid_rows[0]["ae4_cl_inward_fmol_s"]),
                "ae4_share_fold_vs_rest": shares["ae4"] / (grid_rows[0]["ae4_cl_inward_fmol_s"] / (
                    grid_rows[0]["nkcc1_cl_inward_fmol_s"] + grid_rows[0]["ae4_cl_inward_fmol_s"])),
                "note": "No predeclared numerical gate defines substantial AE4 loading; assess reported share and fold change directly."}
    if not integrated_rows:
        integrated_rows.append({"quantity": "all_requested_integrals", "value": "", "unit": "",
            "window_start_s": 60, "window_end_s": 600, "status": "NOT_EVALUATED_AFTER_STOP"})
    readouts = [row for row in grid_rows if row["time_s"] == 0 or (row["time_s"] >= 60 and row["time_s"] % 60 == 0)]
    if first_failure and (not readouts or readouts[-1]["time_s"] != first_failure["time_s"]):
        readouts.append(first_failure["observables"])
    if rows and (not readouts or readouts[-1]["time_s"] != rows[-1]["time_s"]):
        readouts.append(rows[-1])
    write_csv(RESULTS / "wt_timeseries.csv", readouts, list(rows[0]) if rows else ["time_s"])
    write_csv(RESULTS / "wt_integrated_fluxes.csv", integrated_rows,
        ["quantity", "value", "unit", "window_start_s", "window_end_s", "status"])
    write_json(RESULTS / "wt_chloride_partition.json", partition)
    verification = {"status": status, "prepared_head": PREPARED_HEAD, "branch": BRANCH,
        "started_utc": started_utc, "complete_600_s": complete,
        "first_failure": first_failure, "solver_exception": solver_exception,
        "integration_attempts": attempts, "accepted_solver_steps": accepted_steps,
        "checked_state_count": len(rows), "last_checked_time_s": rows[-1]["time_s"] if rows else None,
        "monitoring": "Every accepted solver endpoint plus dense output at each integer second; first detected failed sample stops all further stepping.",
        "onset_convention": "Inherited protocol is REST at exactly t=0 and active for 0<t<=600; inherited production integration starts at 1e-6 s with the unchanged saved state.",
        "state_names": list(model.state_names), "initial_state_vector": list(map(float, y0)),
        "last_checked_state_vector": raw_last, "checkpoint_path": str(CHECKPOINT.relative_to(ROOT)),
        "genotype": "WT", "ae4_expression": 1.0, "stimulus": asdict(model.stimulus),
        "parameter_hashes": {"whole_cell": initial_parameter_hash, "ae4": initial_ae4_hash, "nbc": initial_nbc_hash},
        "scientific_parameters_unchanged": unchanged, "nbc_parameters": asdict(model.nbc_parameters),
        "nhe1_carrier_amount_fmol": CARRIER, "minima": minima, "maxima": maxima,
        "max_abs_conservation_residuals": residual_max, "conservation_tolerances": dict(CONSERVATION_RESIDUAL_TOLERANCES),
        "max_conservation_ratios": ratio_max, "secretion_gate": secretion,
        "software_versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "thread_pools": [{k:v for k,v in pool.items() if k in ("internal_api", "prefix", "num_threads", "version")} for pool in thread_info],
        "readout_status": "COMPLETE" if complete else "PARTIAL_EARLY_STOP",
        "source_fields_note": "minimal_nbc_charge_source_fmol_s=-J_NBC and minimal_nbc_carbon_source_fmol_s=2*J_NBC are nonzero physical sources, not conservation errors."}
    write_json(RESULTS / "wt_verification.json", verification)
    budget["numerical_execution_seconds"] += time.perf_counter() - START
    budget.update({"status": status, "maximum_numerical_execution_seconds": 600,
        "intended_wt_integrations": 1, "completed_wt_integrations": int(complete),
        "parameter_sweeps": 0, "genotype_tests": 0, "scientific_parameters_changed": 0,
        "limits_respected": budget["wt_integration_attempts"] <= 2 and budget["numerical_execution_seconds"] <= 600,
        "no_further_scientific_execution": True})
    write_json(RESULTS / "budget.json", budget)
    (ANALYSIS / "numerical_execution.log").write_text(
        "Task 37 WT-only. Production continuous stiff solver; one worker, one BLAS thread.\n"
        + "\n".join(logs) + "\n" + json.dumps({"status": status, "attempts": attempts,
        "first_failure": first_failure, "budget": budget}, indent=2) + "\n")
    print(json.dumps({"status": status, "complete": complete, "first_failure": first_failure,
                      "secretion_gate": secretion, "partition": partition, "budget": budget}, indent=2))


if __name__ == "__main__":
    main()
