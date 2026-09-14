"""Saved-data-only Task 38 report. This module does not import/evaluate the model."""
from pathlib import Path
import csv
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/38_ae4_perturbation_validation"
WT = ROOT / "results/37_wt_nbc_validation"
HERE = Path(__file__).resolve().parent
CASES = ("ae4_5pct", "ae4_null")


def read_json(path):
    return json.loads(path.read_text())


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def numeric_rows(path):
    return [{k: float(v) for k, v in row.items()} for row in read_csv(path)]


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |"] +
        ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    budget = read_json(OUT / "budget.json")
    assert budget["status"] == "AE4_PERTURBATION_VALIDATION_COMPLETE"
    assert budget["started_cases"] == budget["completed_cases"] == list(CASES)
    assert budget["integration_attempts"] == 2 and budget["numerical_retries"] == 0
    assert budget["no_further_scientific_execution"]
    wt = read_json(WT / "wt_verification.json")
    v = {case: read_json(OUT / (case + "_verification.json")) for case in CASES}
    for case in CASES:
        assert v[case]["status"] == "PASS" and v[case]["first_failure"] is None
        assert max(v[case]["max_conservation_ratios"].values()) <= 1
    ts = {"wt": numeric_rows(WT / "wt_timeseries.csv"),
          **{case: numeric_rows(OUT / (case + "_timeseries.csv")) for case in CASES}}
    flux = {"wt": {r["quantity"]: float(r["value"]) for r in read_csv(WT / "wt_integrated_fluxes.csv")},
            **{case: {r["quantity"]: float(r["value"]) for r in read_csv(OUT / (case + "_integrated_fluxes.csv"))}
               for case in CASES}}
    assert all([r["time_s"] for r in rows] == list(range(0, 601, 60)) for rows in ts.values())
    wt_total = wt["secretion_gate"]["cumulative_0_600_pL"]
    wt_cum60 = wt_total - flux["wt"]["water_outflow"]
    cumulative, flow, physiology = [], [], []
    sparse_cum = 0.0
    for i, wr in enumerate(ts["wt"]):
        t = wr["time_s"]
        if t == 60:
            sparse_cum = wt_cum60
        elif t > 60:
            previous = ts["wt"][i-1]
            sparse_cum += (t - previous["time_s"]) * (wr["q_out_pL_s"] + previous["q_out_pL_s"]) / 2
        basis = ("ZERO" if t == 0 else "FROZEN_TOTAL_MINUS_FROZEN_60_600_INTEGRAL" if t == 60 else
                 "FROZEN_REPORTED_0_600_TOTAL" if t == 600 else "ESTIMATED_60S_TRAPEZOID_ANCHORED_AT_60S")
        wt_cum = wt_total if t == 600 else sparse_cum
        cr = {"time_s": t, "wt_cumulative_0_t_pL": wt_cum, "wt_cumulative_basis": basis}
        fr = {"time_s": t, "wt_q_out_pL_s": wr["q_out_pL_s"]}
        pr = {"time_s": t}
        for name, source in (("wt", wr), *((c, ts[c][i]) for c in CASES)):
            for key in ("na_i_mM", "k_i_mM", "cl_i_mM", "ph_i", "hco3_i_mM", "tic_i_mM", "volume_i_pL"):
                pr[name + "_" + key] = source[key]
        for case in CASES:
            row = ts[case][i]
            cr[case + "_cumulative_0_t_pL"] = row["cumulative_outflow_0_t_pL"]
            cr["wt_minus_" + case + "_cumulative_pL"] = wt_cum - row["cumulative_outflow_0_t_pL"]
            fr[case + "_q_out_pL_s"] = row["q_out_pL_s"]
            fr[case + "_flow_deficit_percent_of_wt"] = 100 * (1 - row["q_out_pL_s"] / wr["q_out_pL_s"])
        cumulative.append(cr); flow.append(fr); physiology.append(pr)
    write_csv(OUT / "cumulative_secretion_comparison.csv", cumulative)
    write_csv(OUT / "flow_comparison.csv", flow)
    write_csv(OUT / "intracellular_comparison.csv", physiology)
    integrated = []
    for name in flux["wt"]:
        row = {"quantity": name, "unit": "pL" if name == "water_outflow" else "fmol",
               "window_start_s": 60, "window_end_s": 600, "wt": flux["wt"][name]}
        for case in CASES:
            row[case] = flux[case][name]
            row[case + "_minus_wt"] = flux[case][name] - flux["wt"][name]
        integrated.append(row)
    write_csv(OUT / "integrated_flux_comparison.csv", integrated)
    sparse_error = sparse_cum - wt_total
    ratio = {case: v[case]["secretion"]["ratio_to_frozen_wt"] for case in CASES}
    pair = {
        "5pct_minus_null_cumulative_pL": v["ae4_5pct"]["secretion"]["cumulative_0_600_pL"] - v["ae4_null"]["secretion"]["cumulative_0_600_pL"],
        "5pct_minus_null_percentage_points_of_wt": 100 * (ratio["ae4_5pct"] - ratio["ae4_null"]),
        "null_to_5pct_cumulative_ratio": ratio["ae4_null"] / ratio["ae4_5pct"],
        "max_abs_flow_difference_on_60s_grid_pL_s": max(abs(a["q_out_pL_s"]-b["q_out_pL_s"]) for a,b in zip(ts[CASES[0]],ts[CASES[1]])),
        "max_abs_chloride_difference_on_60s_grid_mM": max(abs(a["cl_i_mM"]-b["cl_i_mM"]) for a,b in zip(ts[CASES[0]],ts[CASES[1]])),
        "endpoint_ph_null_minus_5pct": ts["ae4_null"][-1]["ph_i"] - ts["ae4_5pct"][-1]["ph_i"],
        "equivalence_threshold": None,
    }
    summary = {"status": budget["status"], "frozen_wt_cumulative_0_600_pL": wt_total,
        "perturbations": {case: {"secretion": v[case]["secretion"],
            "chloride_partition": v[case]["chloride_partition"],
            "nkcc1_positive_loading_change_percent": 100 * (flux[case]["positive_nkcc1_cl_loading"] / flux["wt"]["positive_nkcc1_cl_loading"] - 1),
            "nbc_cycles_change_percent": 100 * (flux[case]["nbc_cycles"] / flux["wt"]["nbc_cycles"] - 1),
            "nhe1_flux_change_percent": 100 * (flux[case]["nhe1_flux"] / flux["wt"]["nhe1_flux"] - 1),
            "largest_tabulated_flow_deficit_percent": max(r[case + "_flow_deficit_percent_of_wt"] for r in flow),
            "endpoint": ts[case][-1]} for case in CASES},
        "5pct_vs_null": pair, "wt_cumulative_resolution": {
            "frozen_wt_cumulative_0_60_pL": wt_cum60,
            "exact_stored_quadrature_available_at_s": [0, 60, 600],
            "estimated_from_frozen_60s_flow_at_s": list(range(120, 541, 60)),
            "60s_quadrature_minus_frozen_total_at_600_pL": sparse_error,
            "caveat": "Endpoint discrepancy is a resolution diagnostic, not a bound on intermediate errors. No WT trajectory was rerun."},
        "held_out_context": {"approximate_experimental_null_deficit_percent": 35,
            "predicted_null_deficit_percent": 100 * (1 - ratio["ae4_null"]), "used_for_fitting": False}}
    write_json(OUT / "comparison_summary.json", summary)
    blocks = ["# Task 38: frozen AE4 perturbation predictions", "`AE4_PERTURBATION_VALIDATION_COMPLETE`",
        "Both authorized trajectories passed every inherited Task 37 physical and conservation gate. Cumulative secretion was "
        f"**{100*ratio['ae4_5pct']:.4f}% of frozen WT at AE4=0.05** and **{100*ratio['ae4_null']:.4f}% at AE4=0.00**. "
        "The null prediction is a small secretion loss, with substantial NKCC1 compensation. It does not reproduce the held-out approximately 35% deficit; nothing was fitted or changed in response.",
        "Work and publication are confined to `codex/task-38-ae4-perturbation-validation`, prepared head "
        "`7aaa53414714bae459227b4b972c63d49e13b79e`. Frozen WT reference: "
        "`2ba0e4e0c2a9e940aa96ed86c3c72456afd0e3aa`. `AGENTS.md` was read first. All 97 prepared files "
        "matched their GitHub blob hashes and remained unchanged; the prepared Git tree was "
        "`d3e17198e16ff0c3a8cd06b1472a23d621880370`. Only the Task 38 instructions differ from the frozen WT tree.",
        "## Frozen setup and execution",
        "The Task 37 factory loaded the accepted Task 31 R09 WT resting checkpoint directly. All 12 core coordinates and "
        "the resting regulatory coordinate (zero) were copied bitwise. Full initial-state SHA-256: "
        "`879a650761e9044691ab6b2a9b1bfc3be8f951ece176925117596cc2baa31abe`. No genotype-specific resting state was solved. "
        "The only changed scientific input was `genotype.ae4_expression`, first 0.05 and then 0.00; all other genotype fields were identical to WT.",
        "The NHE1 carrier amount remained `2.339370005697548e-05 fmol`. The unchanged electrogenic NBC was 1 Na : 2 HCO3, "
        "with capacity `0.11570913197464398 fmol/s`, voltage-dependent affinity/current closure and stimulus recruitment. "
        "The stimulated NHE1 multiplier remained 2.3. Routed AE4 retained the 1 Cl : 1 monovalent cation : 2 HCO3 stoichiometry "
        "and inherited beta/PKA regulation. NKCC1, AE2, pump, channels, CO2, paracellular transport, water, bath and geometry were unchanged.",
        "Both trajectories used CCh 0.3 uM, isoproterenol 5 uM, stimulated Ca 0.25 uM and beta occupancy 1.0 for 600 s. "
        "The inherited onset convention is REST at exactly zero and stimulation for 0<t<=600; integration starts at 1e-6 s "
        "with the unchanged initial state. Each used one continuous production Radau integration: rtol 1e-7; amount/regulatory "
        "atol 1e-10; volume atol 1e-12 pL; max step 2 s. NumPy 2.3.5 and SciPy 1.17.0 match Task 37. No numerical retry was used.",
        "## Secretion",
        table(["Readout", "Frozen WT", "AE4=0.05", "AE4=0.00"], [
            ["Cumulative 0-600 s (pL)", f"{wt_total:.12g}", *[f"{v[c]['secretion']['cumulative_0_600_pL']:.12g}" for c in CASES]],
            ["Ratio to frozen WT", "1", *[f"{ratio[c]:.12g}" for c in CASES]],
            ["Cumulative deficit (%)", "0", *[f"{100*(1-ratio[c]):.6f}" for c in CASES]],
            ["Mean q_out 60-600 s (pL/s)", f"{wt['secretion_gate']['mean_60_600_pL_s']:.12g}", *[f"{v[c]['secretion']['mean_60_600_pL_s']:.12g}" for c in CASES]],
            ["q_out at 600 s (pL/s)", f"{wt['secretion_gate']['endpoint_pL_s']:.12g}", *[f"{v[c]['secretion']['endpoint_pL_s']:.12g}" for c in CASES]],
        ]),
        "Headline totals and flux integrals use the same integer-second trapezoidal quadrature as Task 37. "
        "Both perturbations also exceed the original WT resting secretion thresholds; Task 38's decision gates are the specified physiology/conservation checks.",
        table(["Time (s)", "WT q_out (pL/s)", "5% q_out (pL/s)", "0% q_out (pL/s)", "5% flow deficit (%)", "0% flow deficit (%)"],
            [[int(r['time_s']), f"{r['wt_q_out_pL_s']:.9g}", *[f"{r[c+'_q_out_pL_s']:.9g}" for c in CASES],
              *[f"{r[c+'_flow_deficit_percent_of_wt']:.4f}" for c in CASES]] for r in flow]),
        "## Intracellular time courses",
        "All concentrations are mM; cell volume is pL. The CSV files additionally include lumen osmolarity, voltages, NBC activation/current, NHE1 multiplier, and every requested instantaneous transporter flux."]
    for case, label in zip(CASES, ["AE4=0.05", "AE4=0.00"]):
        blocks += ["### " + label, table(["Time (s)", "Na_i", "K_i", "Cl_i", "pH_i", "HCO3_i", "TIC_i", "Volume"],
            [[int(r['time_s']), *[f"{r[k]:.7g}" for k in ('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL')]] for r in ts[case]])]
    blocks += ["## Integrated fluxes and chloride partition",
        "Window: 60-600 s. Positive chloride loading integrates max(J_Cl,0) for each transporter separately. "
        "Signed AE2 is positive for chloride entry to the cell. NBC cycles and NHE1 are signed inward fluxes; NBC bicarbonate-equivalent influx is twice the cycle integral. "
        "CaCC export is cell-to-lumen; paracellular chloride return is lumen-to-bath.",
        table(["Quantity", "Unit", "Frozen WT", "AE4=0.05", "AE4=0.00"],
              [[r['quantity'], r['unit'], *[f"{r[k]:.10g}" for k in ('wt',*CASES)]] for r in integrated]),
        table(["Positive chloride pool", "Frozen WT", "AE4=0.05", "AE4=0.00"],
              [[name.upper()+" share (%)", f"{100*read_json(WT/'wt_chloride_partition.json')['positive_shares'][name]:.6f}",
                *[f"{100*v[c]['chloride_partition']['positive_shares'][name]:.6f}" for c in CASES]] for name in ('nkcc1','ae4','ae2')]),
        "## Cumulative secretion deficit versus time",
        "Task 37 retained its 0-600 and 60-600 s integrals and q_out at 60-second intervals, but did not save its integer-second trajectory or intermediate cumulative values. "
        f"Its exact stored-quadrature 0-60 s total is therefore {wt_cum60:.12g} pL by subtraction. Values at 120-540 s below "
        "use trapezoidal integration of the frozen 60-second q_out samples anchored at that 60 s total. They are explicitly estimates. "
        "At 600 s the original frozen total is used directly. Perturbation cumulative values use the original one-second grid throughout. "
        f"The uncorrected coarse WT quadrature would exceed the stored 600 s total by {sparse_error:.9g} pL "
        f"({100*sparse_error/wt_total:.6f}%). This is a resolution diagnostic, not a bound on intermediate errors. No WT trajectory was rerun.",
        table(["Time (s)", "WT cumulative (pL)", "WT basis", "5% cumulative (pL)", "0% cumulative (pL)", "WT minus 5% (pL)", "WT minus 0% (pL)"],
            [[int(r['time_s']), f"{r['wt_cumulative_0_t_pL']:.9f}", 'estimated' if 60<r['time_s']<600 else 'stored/derived',
              *[f"{r[c+'_cumulative_0_t_pL']:.9f}" for c in CASES],
              *[f"{r['wt_minus_'+c+'_cumulative_pL']:.9f}" for c in CASES]] for r in cumulative]),
        "## Interpretation",
        f"**5% versus null:** the 5% trajectory secretes {pair['5pct_minus_null_cumulative_pL']:.9g} pL more, "
        f"a difference of {pair['5pct_minus_null_percentage_points_of_wt']:.6f} percentage points of WT. "
        f"Null/5% secretion is {pair['null_to_5pct_cumulative_ratio']:.9f}. These are distinct numerical predictions with a small secretion difference. "
        "No materiality or equivalence threshold was declared, so this is not a formal equivalence claim. "
        f"Their largest chloride difference on the reported grid is {pair['max_abs_chloride_difference_on_60s_grid_mM']:.6f} mM, "
        f"and terminal pH differs by {pair['endpoint_ph_null_minus_5pct']:.6f}.",
        "**Timing:** the null's cumulative shortfall is 0.00129703 pL at 60 s; it grows to approximately 0.00480960 pL at 120 s and "
        "0.00878875 pL at 180 s. The largest tabulated cumulative shortfall is around 480 s, about 0.0190027 pL, then it narrows to the "
        "stored-reference deficit of 0.01815625 pL at 600 s. Thus an initially small cumulative deficit develops over the first minutes, "
        "but instantaneous flow does not progressively diverge after 2-3 minutes: its largest reported deficit is 3.9061% at 120 s, "
        "3.8528% at 180 s, then shrinks. Null flow is 0.9804% above WT at 600 s. The 5% case has the same qualitative time pattern "
        "and is 1.0590% above WT at 600 s. Cumulative loss remains positive because earlier deficits outweigh the late recovery.",
        "**Chloride and pH:** loss of AE4 lowers intracellular chloride at every reported post-onset time. At 600 s, Cl_i is "
        "54.07952 mM in WT, 53.35498 mM at 5%, and 53.16615 mM in null. pH stays within 6.6-7.3 throughout, but rises to "
        "7.19750 and 7.21760 versus WT 7.03483. It is controlled within the declared gate while showing substantial relative alkalinization; "
        "the finite trajectory does not establish a new resting state or long-time stability.",
        "**NKCC1 and AE2:** positive NKCC1 chloride loading rises by 24.1302% at 5% AE4 and 26.8944% in null despite identical capacity, "
        "expression and stimulation parameters. It supplies 97.6949% and 99.7697% of the respective positive pools. AE2 reverses from net "
        "WT chloride export (-0.464995 fmol) to small net import (+0.376029 and +0.455393 fmol). It is positive by the 120 s readout in both "
        "perturbations but contributes only 0.1907% and 0.2303% of the positive pools. The dominant compensatory loading response is NKCC1.",
        "**NBC and NHE1:** NBC recruitment remains one after onset, and the NHE1 multiplier remains 2.3; changes in realized flux arise "
        "within the frozen coupled model as concentrations and membrane voltage evolve. NBC cycle influx falls by 77.0740% and 85.5233% "
        "over 60-600 s, and NHE1 flux falls by 51.0621% and 54.5756%. Both remain inward on monitored stimulated samples. "
        "Bicarbonate accumulates to 14.23208 and 15.28780 mM versus WT 7.27856 mM at 600 s. This behavior is consistent with reduced "
        "AE4 bicarbonate export and altered acid-base/voltage driving forces; it is an interpretation of the coupled trajectory, not an additional intervention.",
        "**Held-out validation:** the predicted null cumulative deficit is 1.83955%, versus the supplied experimental context of approximately 35%. "
        "The frozen model therefore misses the magnitude of that held-out phenotype. No parameter was tuned, mechanism added, or rescue attempted.",
        "## Gate evidence",
        "Task 37's unmodified `diagnose` function was applied to each perturbation evaluation. Each trajectory had 333 accepted solver steps, "
        "2336 RHS evaluations and 934 checked states, including every accepted endpoint and every integer-second dense-output point. "
        "All tracked core amounts/volumes and derived physiological readouts passed the finite/positive and specified range gates. "
        "These are monitored numerical checks, not a formal enclosure between samples.",
        table(["Readout", "AE4=0.05 min to max", "AE4=0.00 min to max", "Gate"],
            [[key, *[f"{v[c]['minima'][key]:.9g} to {v[c]['maxima'][key]:.9g}" for c in CASES], gate]
             for key,gate in [('na_i_mM','<40'),('k_i_mM','50-200'),('cl_i_mM','30-80'),('ph_i','6.6-7.3'),
                              ('hco3_i_mM','0<HCO3<100'),('tic_i_mM','>0'),('volume_i_pL','0<V<3')]]),
        table(["Conservation diagnostic", "Tolerance", "5% max absolute", "0% max absolute"],
            [[key, f"{tol:.4g}", *[f"{v[c]['max_abs_conservation_residuals'][key]:.9g}" for c in CASES]]
             for key,tol in wt['conservation_tolerances'].items()]),
        "The largest tolerance-normalized residual was 0.00592593 at 5% and 0.00588685 in null (cell speciation alkalinity). "
        "The physical NBC source fields were checked as charge=-J_NBC and carbon=2*J_NBC, rather than treated as zero residuals. "
        "No scientific gate failed. The first case's PASS was verified before the null integration started.",
        "## Budget, provenance and files",
        f"Exactly two intended and completed perturbation integrations; zero WT reruns, numerical retries, stationary solves, optimisation calls, sweeps or additional scientific parameter evaluations. "
        f"One scientific agent/worker and one thread per BLAS pool. The numerical runner used {budget['numerical_execution_seconds']:.6f} s of the 600 s limit. "
        "The persisted budget prohibits further scientific execution.",
        "[5% readouts](../../results/38_ae4_perturbation_validation/ae4_5pct_timeseries.csv), "
        "[null readouts](../../results/38_ae4_perturbation_validation/ae4_null_timeseries.csv), "
        "[flow comparison](../../results/38_ae4_perturbation_validation/flow_comparison.csv), "
        "[intracellular comparison](../../results/38_ae4_perturbation_validation/intracellular_comparison.csv), "
        "[cumulative comparison](../../results/38_ae4_perturbation_validation/cumulative_secretion_comparison.csv), "
        "[integrated flux comparison](../../results/38_ae4_perturbation_validation/integrated_flux_comparison.csv), "
        "[summary](../../results/38_ae4_perturbation_validation/comparison_summary.json), "
        "[5% verification](../../results/38_ae4_perturbation_validation/ae4_5pct_verification.json), "
        "[null verification](../../results/38_ae4_perturbation_validation/ae4_null_verification.json), "
        "[frozen input manifest](../../results/38_ae4_perturbation_validation/frozen_inputs.json), "
        "[budget](../../results/38_ae4_perturbation_validation/budget.json).",
        "The minimal runner is [run_perturbation.py](run_perturbation.py); "
        "[summarize_results.py](summarize_results.py) produces the comparisons from saved data only. "
        "Publication uses the connected GitHub integration on the specified branch. No merge to main or further scientific task is performed."]
    (HERE / "final_answer.md").write_text("\n\n".join(blocks) + "\n")
    print(json.dumps({"status": budget['status'], "ratios": ratio, "5pct_vs_null": pair,
        "wt_cumulative_60s_quadrature_endpoint_error_pL": sparse_error}, indent=2))


if __name__ == "__main__":
    main()
