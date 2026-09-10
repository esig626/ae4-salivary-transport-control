"""Post-freeze accounting for Task 16's actual empty modified-candidate set.

Reuse the validated baseline results. This module has no root solver, ODE
runner or WT selection operation, and refuses an unexpected modified set.
"""
from __future__ import annotations

from dataclasses import asdict
import csv
import json
import math
from pathlib import Path

import numpy as np

from . import task16_wt_allocation as task
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14b_five_percent import _five_percent_core, _output_path
from .validation import PRODUCTION_RADAU, sha256_file, sha256_object


def checked_checkpoint():
    receipt = task.read(task.OUT / "pushed_wt_checkpoint.json")
    path = task.OUT / "frozen_candidate_manifest.json"
    assert receipt["remote_ref_verified"] is True
    assert receipt["remote_ref_payload"]["object"]["sha"] == receipt["checkpoint_commit"]
    assert receipt["branch"] == task.BRANCH
    assert receipt["frozen_manifest_sha256"] == sha256_file(path)
    assert task.git("rev-parse", receipt["checkpoint_commit"] + ":" + task.relative(path)) == task.git("hash-object", task.relative(path))
    task.git("merge-base", "--is-ancestor", receipt["checkpoint_commit"], "HEAD")
    manifest = task.read(path)
    assert manifest["contract_sha256"] == task.contract()[1]
    for name, digest in manifest["hashes"].items():
        assert sha256_file(task.REPO / name) == digest, name
    assert all(r["condition"] == "inherited_baseline" for r in manifest["admissible_candidates"]), (
        "This evaluator is specific to the actual zero-modified-candidate checkpoint")
    assert len(manifest["admissible_candidates"]) == 30
    return manifest, receipt


def baseline_pair(decision, model_manifest, ae2_rows, input_ledger):
    assert decision["wt_admissible"] and decision["condition"] == "inherited_baseline"
    assert decision["ae4_capacity_scale"] == decision["nkcc1_capacity_scale"] == 1.
    root, ca = decision["root_id"], decision["calcium_uM"]
    wt_path = task.REPO / decision["trajectory_path"]
    near_path = _output_path(root, ca)
    for path in (wt_path, near_path):
        expected = input_ledger["input_files"][task.relative(path)]
        assert sha256_file(path) == expected["sha256"]
    wt, near = task.read(wt_path), task.read(near_path)
    assert near["status"] == "COMPLETE"
    assert near["root_id"] == root and near["calcium_uM"] == ca
    assert near["genotype"]["ae4_expression"] == .05
    assert near["genotype"]["ae2_expression"] == 1.
    assert near["wt_denominator_sha256"] == sha256_file(wt_path)
    for key in ("whole_cell_parameters_sha256", "ae4_parameters_sha256"):
        assert near[key] == wt[key] == model_manifest["roots"][root][key]
    assert near["selected_core_state_sha256"] == sha256_object(_five_percent_core(root))
    t, w = near["trajectory"], wt["trajectory"]
    assert t["solver"] == w["solver"] == asdict(PRODUCTION_RADAU)
    assert t["time_s"] == w["time_s"] == task.contract()[0]["time_grid_s"]
    assert t["numerical_gate_pass"] and w["numerical_gate_pass"]
    assert np.array_equal(np.asarray(t["states"])[:12, 0], _five_percent_core(root))
    total, wt_total = t["cumulative_flow_pL"][-1], w["cumulative_flow_pL"][-1]
    assert math.isclose(total, float(np.trapezoid(t["flow_pL_s"], t["time_s"])), rel_tol=1e-12)
    assert math.isclose(wt_total, float(np.trapezoid(w["flow_pL_s"], w["time_s"])), rel_tol=1e-12)
    assert wt_total == decision["wt_total_0_600_pL"] and wt_total > 0
    ae2 = ae2_rows[root, ca]
    assert ae2["genotype"] == "AE2" and ae2["status"] == "COMPLETE"
    assert ae2["paired_numerical_gate_pass"] == "True"
    assert ae2["prediction_contract_sha256"] == task.read(
        task.REPO / "results/14B_five_percent_ae4_check/frozen_inputs.json")["prediction_contract_sha256"]
    r4, r2 = total / wt_total, float(ae2["R_total"])
    d4, d2 = 1-r4, 1-r2
    return {"status": "BASELINE_HASH_VALID_REUSE", "genotype_numerical_status": "COMPLETE",
        "ae4_expression": .05, "ae2_expression_for_loss": 0.,
        "wt_total_0_600_pL": wt_total, "ae4_5pct_total_0_600_pL": total,
        "R_AE4_5pct": r4, "D_AE4_5pct": d4, "R_AE2": r2, "D_AE2": d2,
        "AE4_minus_AE2_reduction": d4-d2,
        "ae4_5pct_source_path": task.relative(near_path), "ae4_5pct_source_sha256": sha256_file(near_path),
        "ae2_source_path": "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv",
        "ae2_reuse_scope": "Frozen validated matched ratio; no new AE2 trajectory",
        "paired_numerical_gate_pass": True, "new_genotype_simulation": False}


def evaluate():
    freeze, receipt = checked_checkpoint()
    assert not (task.OUT / "genotype_results_frozen.json").exists(), "Genotype results are immutable"
    inputs = task.check_inputs()
    model_manifest, _ = load_freeze()
    ae2_path = task.REPO / "results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv"
    # The historical 187-file sign/slip inventory covers WT and AE4, not
    # the AE2 summary. Anchor this additional input directly to Task 14B.
    ae2_blob = task.git("hash-object", task.relative(ae2_path))
    assert ae2_blob == task.git("rev-parse", task.BASE + ":" + task.relative(ae2_path))
    with ae2_path.open() as f:
        source_rows = list(csv.DictReader(f))
    ae2_rows = {(r["root_id"], float(r["calcium_uM"])): r for r in source_rows}
    assert len(ae2_rows) == len(source_rows) == 30
    rows = []
    for decision in freeze["all_decisions"]:
        row = {k: decision.get(k) for k in ("candidate_id", "root_id", "routing_family", "condition",
            "target_share", "inherited_share", "realized_ae4_share", "ae4_capacity_scale", "nkcc1_capacity_scale",
            "calcium_uM", "na_i_mM", "k_i_mM", "cl_i_mM", "hco3_i_mM", "ph_i", "volume_i_pL",
            "ae4_cancellation_fraction", "ae4_productive_fraction", "ae4_gross_branch_fmol_s",
            "ae4_net_cl_fmol_s", "ae2_signed_cl_fmol_s", "wt_total_0_600_pL", "wt_admissible")}
        row["wt_status"] = decision["status"]
        if decision["wt_admissible"]:
            row.update(baseline_pair(decision, model_manifest, ae2_rows, inputs))
        else:
            row.update(status="NOT_EVALUATED_WT_INADMISSIBLE", genotype_numerical_status="NOT_LICENSED",
                paired_numerical_gate_pass=False, new_genotype_simulation=False,
                exclusion_reason=decision.get("gate_failures", decision["status"]))
        rows.append(row)
    assert len(rows) == 90 and sum(r["wt_admissible"] for r in rows) == 30
    common = ["candidate_id", "root_id", "routing_family", "condition", "target_share", "calcium_uM",
        "wt_admissible", "wt_status", "status", "genotype_numerical_status", "wt_total_0_600_pL",
        "paired_numerical_gate_pass", "new_genotype_simulation", "exclusion_reason"]
    for filename, keys in (("ae4_5pct_results.csv", common+["ae4_expression", "R_AE4_5pct", "D_AE4_5pct",
            "ae4_5pct_total_0_600_pL", "ae4_5pct_source_path", "ae4_5pct_source_sha256"]),
            ("ae2_results.csv", common+["ae2_expression_for_loss", "R_AE2", "D_AE2", "ae2_source_path", "ae2_reuse_scope"])):
        write_rows(task.OUT / filename, [{k: r.get(k) for k in keys} for r in rows])
    write_rows(task.OUT / "allocation_vs_genotype.csv", rows)
    valid = [r for r in rows if r["wt_admissible"]]
    ranges = {key: {"minimum": min(r[key] for r in valid), "maximum": max(r[key] for r in valid)}
        for key in ("R_AE4_5pct", "D_AE4_5pct", "R_AE2", "D_AE2", "AE4_minus_AE2_reduction")}
    result = {"created_utc": task.now(), "pushed_wt_checkpoint": receipt["checkpoint_commit"],
        "wt_manifest_sha256": receipt["frozen_manifest_sha256"], "expected_combinations": 90,
        "baseline_paired_results_reused": 30, "modified_wt_admissible_combinations": 0,
        "modified_combinations_explicitly_unavailable": 60, "new_ae4_near_loss_runs": 0,
        "new_ae2_loss_runs": 0, "new_genotype_resting_solves": 0,
        "new_exact_zero_ae4_attempts": 0, "all_frozen_candidates_accounted": True,
        "baseline_ranges": ranges, "experimental_context_not_used_in_results": True,
        "additional_baseline_inputs": {task.relative(ae2_path): {"git_blob_sha1": ae2_blob,
            "sha256": sha256_file(ae2_path), "bytes": ae2_path.stat().st_size,
            "baseline_commit": task.BASE}},
        "hashes": {task.relative(task.OUT / f): sha256_file(task.OUT / f) for f in
            ("ae4_5pct_results.csv", "ae2_results.csv", "allocation_vs_genotype.csv")}}
    write_json(task.OUT / "genotype_results_frozen.json", result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    evaluate()
