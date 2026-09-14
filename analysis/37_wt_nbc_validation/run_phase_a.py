"""Direct Task 37 REST nesting and bookkeeping audit; zero integrations."""
import time
START = time.perf_counter()
from validation_common import *
from modern_full_model.nbc_minimal import evaluate_minimal_nbc, nbc_affinity_log


def main():
    budget = json.loads((RESULTS / "budget.json").read_text())
    assert budget["focused_tests_pass"] and budget["wt_integration_attempts"] == 0
    saved, base, rest, stim, y0 = models_and_state()
    frozen_y0 = y0.copy()
    baseline = base.evaluate(0.0, y0, genotype=WT)
    nested = rest.evaluate(0.0, y0, genotype=WT)
    delta = np.asarray(nested.rhs) - np.asarray(baseline.rhs)
    mem0, mem1 = baseline.diagnostics.membranes, nested.diagnostics.membranes
    voltage_differences = {
        name: float(getattr(mem1, name) - getattr(mem0, name))
        for name in ("v_apical_V", "v_basolateral_V", "v_transepithelial_V")
    }
    rr, physical_failures, ratios = diagnose(rest, 0.0, y0, nested)
    # A fixed-state, stimulated RHS evaluation is an implementation check,
    # not a trajectory or stationary solve.
    onset = stim.evaluate(1e-6, y0, genotype=WT)
    local = bookkeeping(stim, onset)
    r = onset.diagnostics.regulatory
    local_residuals = onset.diagnostics.conservation_residuals
    ci = onset.diagnostics.observables.cell_concentrations_mM
    bath_hco3 = onset.diagnostics.observables.bath_acid_base.hco3_mM
    p = stim.parameters
    chemical = math.log(p.bath.na_mM / ci["na"]) + 2 * math.log(bath_hco3 / ci["hco3"])
    reversal = -p.constants.thermal_voltage_V * chemical
    args = dict(na_i_mM=ci["na"], hco3_i_mM=ci["hco3"], na_o_mM=p.bath.na_mM,
                hco3_o_mM=bath_hco3, thermal_voltage_V=p.constants.thermal_voltage_V,
                faraday_C_mol=p.constants.faraday_C_mol, activation_fraction=1.0)
    forward = evaluate_minimal_nbc(**args, v_basolateral_V=reversal + 0.02)
    reverse = evaluate_minimal_nbc(**args, v_basolateral_V=reversal - 0.02)
    at_reversal = evaluate_minimal_nbc(**args, v_basolateral_V=reversal)
    isolated = []
    for label, flux in (("forward", forward), ("reverse", reverse)):
        j = flux.cycle_inward_fmol_s
        residual = flux.transported_charge_equivalents_fmol_s + flux.current_A / p.constants.faraday_C_mol * 1e15
        isolated.append({"direction": label, "cycle_fmol_s": j,
            "source_signature_per_cycle": [flux.na_cell_fmol_s / j, 0.0, flux.cl_cell_fmol_s / j,
                                           flux.tic_cell_fmol_s / j, flux.alkalinity_cell_fmol_s / j],
            "charge_equivalents_per_cycle": flux.transported_charge_equivalents_fmol_s / j,
            "source_plus_current_residual_fmol_s": residual,
            "tic_minus_2cycle_fmol_s": flux.tic_cell_fmol_s - 2*j,
            "ta_minus_2cycle_fmol_s": flux.alkalinity_cell_fmol_s - 2*j})
    checks = {
        "saved_state_unchanged": bool(np.array_equal(y0, frozen_y0) and np.array_equal(y0[:12], saved["core_state"])),
        "all_rhs_rows_bitwise_equal": bool(np.array_equal(nested.rhs, baseline.rhs)),
        "amount_rhs_within_1e_minus_12": bool(np.max(np.abs(delta[AMOUNT_ROWS])) <= 1e-12),
        "volume_rhs_within_1e_minus_14": bool(np.max(np.abs(delta[VOLUME_ROWS])) <= 1e-14),
        "voltages_within_1e_minus_12": max(map(abs, voltage_differences.values())) <= 1e-12,
        "q_out_within_1e_minus_12": abs(rr["q_out_pL_s"] - baseline.diagnostics.water.lumen_outflow_pL_s) <= 1e-12,
        "accepted_q_out_reproduced": abs(rr["q_out_pL_s"] - Q_REST) <= 1e-12,
        "rest_nbc_exact_zero": rr["nbc_activation"] == rr["nbc_cycle_inward_fmol_s"] == rr["nbc_current_A"] == 0.0,
        "rest_nhe_exact_one": rr["nhe1_multiplier"] == 1.0,
        "stim_nbc_exact_one_and_nhe_exact_2p3": r["minimal_nbc_activation_fraction"] == 1.0 and r["nhe1_stimulation_multiplier"] == 2.3,
        "rest_physical_and_conservation": not physical_failures,
        "stimulated_local_conservation": all(abs(local_residuals[k]) <= tol for k, tol in CONSERVATION_RESIDUAL_TOLERANCES.items()),
        "whole_cell_nbc_bookkeeping_roundoff": max(abs(local[k]) for k in local if k != "isolated_nbc_source") < 1e-12,
        "isolated_nbc_forward_and_reverse_bookkeeping": all(
            x["source_signature_per_cycle"] == [1.0, 0.0, 0.0, 2.0, 2.0]
            and x["charge_equivalents_per_cycle"] == -1.0
            and abs(x["source_plus_current_residual_fmol_s"]) < 1e-14
            and x["tic_minus_2cycle_fmol_s"] == x["ta_minus_2cycle_fmol_s"] == 0.0
            for x in isolated),
        "full_affinity_reversal": forward.cycle_inward_fmol_s > 0 and reverse.cycle_inward_fmol_s < 0 and abs(at_reversal.affinity_log) < 1e-14,
        "outward_current_sign_for_inward_negative_charge": forward.current_A > 0 and reverse.current_A < 0,
    }
    report = {
        "status": "PASS" if all(checks.values()) else "TASK37_IMPLEMENTATION_OR_REST_NESTING_FAILED",
        "prepared_head": PREPARED_HEAD, "branch": BRANCH,
        "checkpoint_path": str(CHECKPOINT.relative_to(ROOT)),
        "checkpoint_label": saved["label"], "core_state_sha256": saved["core_state_sha256"],
        "active_parameter_sha256": saved["active_parameter_sha256"],
        "ae4_parameter_sha256": saved["ae4_parameter_sha256"],
        "state_names": list(rest.state_names), "state_vector": list(map(float, y0)),
        "nhe1_carrier_amount_fmol": CARRIER, "nbc_parameters": asdict(rest.nbc_parameters),
        "checks": checks, "max_abs_amount_rhs_difference_fmol_s": float(np.max(np.abs(delta[AMOUNT_ROWS]))),
        "max_abs_volume_rhs_difference_pL_s": float(np.max(np.abs(delta[VOLUME_ROWS]))),
        "voltage_differences_V": voltage_differences,
        "q_out_task31_pL_s": float(baseline.diagnostics.water.lumen_outflow_pL_s),
        "q_out_task36_pL_s": rr["q_out_pL_s"],
        "q_out_difference_pL_s": float(rr["q_out_pL_s"] - baseline.diagnostics.water.lumen_outflow_pL_s),
        "task31_rhs": list(map(float, baseline.rhs)), "task36_rhs": list(map(float, nested.rhs)),
        "rest_observables": rr, "rest_gate_failures": physical_failures,
        "rest_conservation_ratios": ratios, "stimulated_fixed_state_bookkeeping": local,
        "isolated_nbc_bookkeeping": isolated,
        "reversal_voltage_V": reversal, "affinity_at_reversal": at_reversal.affinity_log,
        "notes": ["No stationary solve or dynamic integration.",
                  "The two minimal_nbc_*_source fields are physical sources; closure is tested against their declared values, not against zero.",
                  "The inherited R09 root/routing identifier contains AE4NA05; active WT expression remains exactly 1.0."],
    }
    write_json(RESULTS / "rest_nesting.json", report)
    budget["numerical_execution_seconds"] += time.perf_counter() - START
    budget["rest_nesting_pass"] = all(checks.values())
    budget["direct_rest_nesting_runs"] = 1
    write_json(RESULTS / "budget.json", budget)
    print(json.dumps({"status": report["status"], "checks": checks,
        "amount_rhs_difference": report["max_abs_amount_rhs_difference_fmol_s"],
        "volume_rhs_difference": report["max_abs_volume_rhs_difference_pL_s"],
        "q_out": rr["q_out_pL_s"], "bookkeeping": local}, indent=2))


if __name__ == "__main__":
    main()
