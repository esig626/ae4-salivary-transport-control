from __future__ import annotations

import csv
import inspect
import math

import numpy as np

from src.modern_full_model.run_calcium_fast_screen import load_freeze
from src.modern_full_model import task18_fixed_wt as task


def payloads():
    return [task.read(path) for path in sorted((task.OUT / "wt_parameter_payloads").glob("*.json"))]


def test_contract_has_only_the_three_predeclared_conditions_and_no_phenotype_inputs():
    contract = task.read(task.OUT / "contract.json")
    assert tuple(contract["conditions"]) == task.CONDITIONS
    assert contract["target_shares"] == [0.10, 0.30]
    assert contract["genotype_information_in_objective"] is False
    assert contract["genotype_information_read_before_freeze"] is False
    assert len(contract["root_ids"]) == 10
    lowered = "\n".join(contract["prefreeze_input_hashes"]).lower()
    assert "ae4_5pct" not in lowered
    assert "ae2_results" not in lowered
    assert "genotype_results" not in lowered


def test_source_rank_keeps_both_voltages_and_both_current_equations():
    source = task.read(task.OUT / "source_rank.json")
    assert len(source) == 10
    for result in source.values():
        augmented = result["explicit_voltage_augmented"]
        assert augmented["variable_names"][-2:] == ["v_apical_mV", "v_basolateral_mV"]
        assert augmented["equation_names"][-2:] == [
            "apical_current_fmol_charge_s",
            "basolateral_current_fmol_charge_s",
        ]
        assert result["rank"] == 10
        assert result["nullity"] == 4
        assert result["left_nullity"] == 4
        assert augmented["baseline_max_abs_residual"] < 1.0e-10


def test_complete_ensemble_preserves_exact_inherited_states_and_full_rhs():
    manifest, _ = load_freeze()
    items = payloads()
    assert len(items) == 30
    assert {(item["row"]["root_id"], item["row"]["condition"]) for item in items} == {
        (root_id, condition)
        for root_id in manifest["roots"]
        for condition in task.CONDITIONS
    }
    for item in items:
        root_id = item["row"]["root_id"]
        model, state, baseline_evaluation = task.baseline(manifest, root_id)
        assert item["core_state"] == manifest["roots"][root_id]["core_state"]
        assert item["complete_state"] == state.tolist()
        assert item["complete_state_bytes_sha256"] == task.state_bytes_sha256(state)
        modified = task.model_with_capacities(model, item["accepted_capacities"])
        evaluation = modified.evaluate(0.0, state)
        summary = task.residual_summary(evaluation)
        assert summary["max_abs_scaled_independent_rhs"] <= task.SPEC.root_scaled_tolerance
        assert summary["max_abs_omitted_rhs"] <= task.SPEC.omitted_row_raw_tolerance
        assert summary["max_conservation_tolerance_ratio"] <= 1.0
        assert summary["max_abs_current_A"] <= 1.0e-20


def test_ae4_and_other_positive_loader_are_scaled_as_complete_transporters():
    manifest, _ = load_freeze()
    for item in payloads():
        root_id = item["row"]["root_id"]
        model, state, baseline_evaluation = task.baseline(manifest, root_id)
        modified = task.model_with_capacities(model, item["accepted_capacities"])
        evaluation = modified.evaluate(0.0, state)
        lam = item["allocation"]["ae4_capacity_multiplier"]
        rho = item["allocation"]["nkcc1_capacity_multiplier"]
        before = task._ae4_sources(baseline_evaluation)
        after = task._ae4_sources(evaluation)
        for species in ("na", "k", "cl", "tic", "alkalinity"):
            assert math.isclose(after[species], lam * before[species], rel_tol=2.0e-12, abs_tol=2.0e-14)
        assert math.isclose(
            evaluation.diagnostics.homeostasis.nkcc1_inward_fmol_s,
            rho * baseline_evaluation.diagnostics.homeostasis.nkcc1_inward_fmol_s,
            rel_tol=2.0e-12,
            abs_tol=2.0e-14,
        )
        assert math.isclose(
            evaluation.diagnostics.homeostasis.ae2_inward_fmol_s,
            baseline_evaluation.diagnostics.homeostasis.ae2_inward_fmol_s,
            rel_tol=2.0e-12,
            abs_tol=2.0e-14,
        )


def test_target_share_pool_and_nonloading_chloride_fluxes_are_preserved():
    manifest, _ = load_freeze()
    for item in payloads():
        root_id = item["row"]["root_id"]
        model, state, baseline_evaluation = task.baseline(manifest, root_id)
        modified = task.model_with_capacities(model, item["accepted_capacities"])
        evaluation = modified.evaluate(0.0, state)
        inherited = task.chloride_ledger(baseline_evaluation)
        accepted = task.chloride_ledger(evaluation)
        pool = accepted["AE4"] + accepted["NKCC1"]
        assert math.isclose(
            pool,
            item["allocation"]["inherited_positive_loading_pool_fmol_s"],
            rel_tol=2.0e-11,
            abs_tol=2.0e-13,
        )
        expected = item["allocation"]["target_share"]
        if expected is None:
            expected = inherited["AE4"] / (inherited["AE4"] + inherited["NKCC1"])
        assert math.isclose(accepted["AE4"] / pool, expected, rel_tol=2.0e-11, abs_tol=2.0e-13)
        for pathway in ("AE2", "CaCC_cell_to_lumen", "paracellular_Cl_lumen_to_bath"):
            assert math.isclose(accepted[pathway], inherited[pathway], rel_tol=2.0e-9, abs_tol=2.0e-12)


def test_historical_capacity_values_are_references_not_feasibility_bounds():
    rows = list(csv.DictReader((task.OUT / "wt_inverse_solutions.csv").open()))
    modified = [row for row in rows if row["condition"] != "inherited_baseline"]
    assert len(modified) == 20
    assert all(row["fixed_wt_feasible"] == "True" for row in modified)
    assert any(float(row["g_k_total_nS"]) > 14.0 for row in modified)
    assert any(float(row["g_cl_apical_nS"]) < 31.4 for row in modified)
    assert max(float(row["g_k_total_nS"]) for row in modified) > 1000.0
    source = inspect.getsource(task.solve_canonical)
    assert "14e-9" not in source
    assert "31.4" not in source
    assert "native_K" not in source


def test_canonical_objective_is_recorded_and_satisfied_by_every_payload():
    for item in payloads():
        ref = item["reference_capacities"]
        accepted = item["accepted_capacities"]
        logs = np.asarray(
            [math.log(accepted[name] / ref[name]) for name in task.OBJECTIVE_CAPACITIES]
        )
        assert math.isclose(
            np.max(np.abs(logs)), item["row"]["largest_abs_log_fold"], rel_tol=2.0e-11
        )
        assert math.isclose(logs @ logs, item["row"]["sum_squared_log_fold"], rel_tol=2.0e-11)
        assert item["stage1_optimisation"]["success"] is True
        assert item["stage2_optimisation"]["success"] is True
        assert item["phenotype_information_used"] is False


def test_frozen_manifest_hashes_when_present():
    path = task.OUT / "frozen_wt_manifest.json"
    if not path.exists():
        return
    frozen = task.read(path)
    assert frozen["genotype_outcomes_used"] is False
    assert frozen["genotype_result_files_read"] == 0
    assert len(frozen["all_decisions"]) == 30
    assert len(frozen["feasible_parameterisations"]) == 30
    for name, expected in frozen["hashes"].items():
        assert task.sha256_file(task.REPO / name) == expected, name
