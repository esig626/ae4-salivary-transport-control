"""Read-only equation diagnostics of saved WT solver attempts; never re-solve."""
from __future__ import annotations

import csv
import numpy as np

from . import task16_wt_allocation as task
from .calibration import WTCalibrationSpec, _wt_gate_failures
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .validation import sha256_file


def diagnose():
    assert not (task.OUT / "frozen_candidate_manifest.json").exists()
    manifest, _ = load_freeze()
    failures = []
    for candidate in task.contract()[0]["candidate_definitions"]:
        if candidate["condition"] == "inherited_baseline":
            continue
        saved = task.read(task.root_path(candidate))
        if saved["row"]["numerical_rest_pass"]:
            continue
        connected = [s for s in saved["continuation_steps"] if s["connected"]]
        if not connected:
            raise RuntimeError("Failure before any connected point needs separate inspection")
        last = connected[-1]
        # This is the last existing internal numerical state, never a new allocation candidate.
        model = task.allocated_model(manifest, candidate["root_id"], .10,
            candidate["ae4_capacity_scale"] ** last["homotopy_t"],
            candidate["nkcc1_capacity_scale"] ** last["homotopy_t"])
        values = task.partition(model, last["selected_root"]["core_state"])
        final = saved["continuation_steps"][-1]
        coordinates = np.asarray([a["coordinates"] for a in final["attempts"]])
        spec = WTCalibrationSpec()
        hits = [b.name for i, b in enumerate(spec.coordinate_bounds)
            if np.any(np.minimum(coordinates[:, i]-b.lower, b.upper-coordinates[:, i])
                <= spec.boundary_relative_tolerance*(b.upper-b.lower))]
        row = {**candidate, "diagnostic_scope": "SAVED_INTERNAL_NUMERICAL_STATES_ONLY_NOT_ENDPOINT_MODELS",
            "last_connected_internal_t": last["homotopy_t"],
            "last_connected_rank": last["selected_root"]["normalized_jacobian_rank"],
            "last_connected_scaled_residual": last["selected_root"]["max_abs_scaled_independent_rhs"],
            "last_connected_wt_gate_failures": "; ".join(_wt_gate_failures(values, spec)),
            "final_failed_internal_t": final["homotopy_t"], "final_failed_boundary_hits": "; ".join(hits),
            "final_failed_min_scaled_residual": min(a["max_abs_scaled_independent_rhs"] for a in final["attempts"]),
            "fixed_endpoint_attempted": any(s["homotopy_t"] == 1 for s in saved["continuation_steps"]),
            "fixed_endpoint_solved": False, "internal_attempted_steps": len(saved["continuation_steps"]),
            "root_path": task.relative(task.root_path(candidate)),
            "root_sha256": sha256_file(task.root_path(candidate)),
            **{"last_connected_"+k: v for k, v in values.items()}}
        failures.append(row)
    write_rows(task.OUT / "wt_continuation_failure_diagnostics.csv", failures)
    # Reuse only the WT rows from the frozen solver summary, with no new integrations.
    source = task.REPO / "results/14_scale_free_genotype_holdout/solver_confirmation.json"
    wt_checks = [r for r in task.read(source)["checks"] if r["genotype"] == "WT"]
    assert len(wt_checks) == 6 and all(r["pass"] for r in wt_checks)
    write_json(task.OUT / "inherited_wt_solver_confirmations.json", {
        "source_path": task.relative(source), "source_sha256": sha256_file(source),
        "only_wt_checks_used": True, "checks": wt_checks})
    write_json(task.OUT / "wt_numerical_disposition.json", {
        "wt_only": True, "new_rest_parameter_fits": 0,
        "modified_endpoint_count": 20, "modified_endpoint_roots_solved": 20-len(failures),
        "modified_continuation_failures": len(failures),
        "all_failed_at_inherited_na_search_boundary": all(r["final_failed_boundary_hits"] == "na_i_mM" for r in failures),
        "all_last_connected_states_fail_wt_physiology": all(r["last_connected_wt_gate_failures"] for r in failures),
        "new_wt_ode_runs": len(list((task.OUT / "wt_trajectories").glob("*.json"))),
        "baseline_wt_trajectories_reused": 30,
        "classification_limit": "Failure to continue inside frozen search domains does not prove nonexistence of all possible endpoint equilibria. No domain expansion, disconnected-root search or additional share is licensed.",
        "post_freeze_modified_genotype_runs_licensed": 0 if len(failures) == 20 else None,
    })


if __name__ == "__main__":
    diagnose()
