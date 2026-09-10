from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from src.modern_full_model import task18_fixed_wt as task
from src.modern_full_model import task18_postfreeze as post
from src.modern_full_model.validation import sha256_file


def rows(name: str) -> list[dict[str, str]]:
    with (task.OUT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def postfreeze_payloads() -> list[dict]:
    return [
        json.loads(path.read_text())
        for path in sorted((task.OUT / "postfreeze_payloads").glob("*.json"))
    ]


def test_pushed_wt_checkpoint_and_reuse_hashes_are_still_valid():
    frozen, receipt = post.checked_checkpoint()
    verified = post.verify_reuse_inputs()
    assert receipt["checkpoint_sha"] == post.EXPECTED_WT_CHECKPOINT
    assert len(frozen["feasible_parameterisations"]) == 30
    assert frozen["genotype_outcomes_used"] is False
    assert frozen["genotype_result_files_read"] == 0
    assert len(verified) == 71


def test_all_four_dynamic_tables_have_the_complete_cartesian_ensemble():
    expected = {
        (root_id, condition, calcium)
        for root_id in task.read(task.OUT / "contract.json")["root_ids"]
        for condition in task.CONDITIONS
        for calcium in task.CALCIUM
    }
    for name in (
        "wt_dynamic_results.csv",
        "ae4_5pct_results.csv",
        "ae2_results.csv",
        "final_comparison.csv",
    ):
        actual_rows = rows(name)
        assert len(actual_rows) == 90
        actual = {
            (row["root_id"], row["condition"], float(row["calcium_uM"]))
            for row in actual_rows
        }
        assert actual == expected


def test_every_trajectory_and_paired_comparison_passes_numerical_gates():
    for name in (
        "wt_dynamic_results.csv",
        "ae4_5pct_results.csv",
        "ae2_results.csv",
    ):
        actual = rows(name)
        assert {row["status"] for row in actual} == {"COMPLETE"}
        assert {row["numerical_gate_pass"] for row in actual} == {"True"}
        conservation = [
            float(row["max_dimensionless_conservation_ratio"])
            for row in actual
            if row["max_dimensionless_conservation_ratio"]
        ]
        assert conservation
        assert max(conservation) < 1.0
        assert all(math.isfinite(float(row["total_0_600_pL"])) for row in actual)
    comparison = rows("final_comparison.csv")
    assert {row["wt_numerical_gate_pass"] for row in comparison} == {"True"}
    assert {row["ae4_paired_numerical_gate_pass"] for row in comparison} == {"True"}
    assert {row["ae2_paired_numerical_gate_pass"] for row in comparison} == {"True"}
    assert {row["capacities_refitted_after_freeze"] for row in comparison} == {"False"}


def test_modified_continuations_are_complete_without_old_calibration_boxes():
    payloads = postfreeze_payloads()
    assert len(payloads) == 20
    ae4_volumes = []
    for payload in payloads:
        assert payload["fitted_parameter_count_after_freeze"] == 0
        assert payload["exact_ae4_zero_attempts"] == 0
        assert sha256_file(task.REPO / payload["wt_payload_path"]) == payload["wt_payload_sha256"]
        ae4 = payload["ae4_5pct_continuation"]
        ae2 = payload["ae2_loss_continuation"]
        assert ae4["complete"] is True
        assert ae4["terminal_expression"] == 0.05
        assert ae4["exact_zero_attempted"] is False
        assert ae4["selected_root"]["normalized_jacobian_rank"] == 10
        assert ae4["selected_root"]["boundary_hits"] == []
        assert ae2["complete"] is True
        assert ae2["terminal_expression"] == 0.0
        assert ae2["exact_zero_attempted"] is True
        assert ae2["selected_root"]["normalized_jacobian_rank"] == 10
        assert ae2["selected_root"]["boundary_hits"] == []
        ae4_volumes.append(ae4["selected_root"]["coordinates"][4])
    assert max(ae4_volumes) > 3.0
    assert max(ae4_volumes) < post.TASK18_CONTINUATION_LIMITS[4][1]


def test_reported_ratios_reproduce_the_matched_trajectory_integrals():
    for row in rows("final_comparison.csv"):
        wt = float(row["wt_total_0_600_pL"])
        ae4 = float(row["ae4_5pct_total_0_600_pL"])
        ae2 = float(row["ae2_total_0_600_pL"])
        assert math.isclose(float(row["R_AE4_5pct"]), ae4 / wt, rel_tol=2.0e-14)
        assert math.isclose(float(row["R_AE2"]), ae2 / wt, rel_tol=2.0e-14)
        assert math.isclose(float(row["D_AE4_5pct"]), 1.0 - ae4 / wt, rel_tol=2.0e-14, abs_tol=2.0e-14)
        assert math.isclose(float(row["D_AE2"]), 1.0 - ae2 / wt, rel_tol=2.0e-14, abs_tol=2.0e-14)


def test_postfreeze_summary_hashes_and_accounting_are_exact():
    summary = task.read(task.OUT / "postfreeze_summary.json")
    assert summary["expected_rows"] == 90
    assert summary["modified_wt_numerically_valid"] == 60
    assert summary["modified_ae4_pairs_numerically_valid"] == 60
    assert summary["modified_ae2_pairs_numerically_valid"] == 60
    assert summary["ae4_exact_zero_attempts"] == 0
    assert summary["postfreeze_parameter_refits"] == 0
    for name, expected in summary["result_hashes"].items():
        assert sha256_file(task.REPO / name) == expected
