"""Generate Task 13C tables and reports exclusively from saved WT results."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from . import run_calcium_fast_screen as run


ANALYSIS = run.REPO / "analysis/13C_calcium_fast_screen"
CLASSIFICATIONS = (
    "WT ABSOLUTE FLOW GATE PASSES AT 0.25 UM",
    "WT ABSOLUTE FLOW GATE FAILS AT 0.25 UM BUT PASSES AT 0.50 UM",
    "CALCIUM 0.50 UM REMAINS INSUFFICIENT TO CLOSE WT ABSOLUTE FLOW DEFICIT",
    "CALCIUM INCREASE EXPOSES PHYSIOLOGICAL OR NUMERICAL FAILURE BEFORE FLOW CLOSURE",
)


def read_rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            row = {}
            for key, value in raw.items():
                if value in ("True", "False"):
                    row[key] = value == "True"
                elif value == "":
                    row[key] = None
                else:
                    try:
                        row[key] = float(value)
                    except ValueError:
                        row[key] = value
            rows.append(row)
    return rows


def classify(rows: list[dict], confirmations: list[dict]) -> str:
    broken_crosscheck = any(r.get("solver_crosscheck_gate_pass") is False for r in confirmations
                           if r["solver_label"] == run.PRODUCTION_BDF.label)
    if broken_crosscheck:
        return CLASSIFICATIONS[3]
    for index, ca in enumerate(run.NEW_CALCIUM):
        if any(r["calcium_uM"] == ca and r["scientific_gate_pass"] for r in confirmations
               if r["solver_label"] == run.PRODUCTION_RADAU.label):
            return CLASSIFICATIONS[index]
    if any(not r["numerical_gate_pass"] or not r["sustainment_gate_pass"] for r in rows + confirmations):
        return CLASSIFICATIONS[3]
    return CLASSIFICATIONS[2]


def validate_complete_results(manifest: dict, screens: list[dict], confirmations: list[dict]) -> None:
    def requests(rows):
        return [run.Request(r["root_id"], r["calcium_uM"], r["solver_label"]) for r in rows]
    expected_screen = run.primary_requests(manifest)
    if len(screens) != 20 or set(requests(screens)) != set(expected_screen):
        raise ValueError("Report requires the complete twenty case primary screen")
    expected_radau, _ = run.production_requests(screens)
    actual_radau = [r for r in confirmations if r["solver_label"] == run.PRODUCTION_RADAU.label]
    expected_bdf = run.bdf_requests(actual_radau)
    if (set(requests(actual_radau)) != set(expected_radau)
        or len(confirmations) != len(expected_radau) + len(expected_bdf)
        or set(requests(confirmations)) != set(expected_radau + expected_bdf)):
        raise ValueError("Report requires exactly the predeclared confirmations")


def table(headers: list[str], rows: list[list[Any]]) -> str:
    return ("| " + " | ".join(headers) + " |\n| " + " | ".join("---" for _ in headers)
            + " |\n" + "\n".join("| " + " | ".join(map(str, row)) + " |" for row in rows) + "\n")


def span(rows: list[dict], key: str, scale: float = 1.0, digits: int = 6) -> str:
    values = [float(r[key]) * scale for r in rows]
    return f"{min(values):.{digits}f} to {max(values):.{digits}f}"


def main() -> None:
    manifest, verification = run.load_freeze()
    screens = read_rows(run.OUTPUT / "calcium_screen.csv")
    confirmation = read_rows(run.OUTPUT / "calcium_confirmatory.csv")
    validate_complete_results(manifest, screens, confirmation)
    classification = classify(screens, confirmation)
    if classification != CLASSIFICATIONS[2]:
        raise ValueError("This report describes the completed negative screen; review unexpected results before reporting")
    references = read_rows(run.OUTPUT / "calcium_reference_frozen.csv")
    for row in references:
        minute = [row[f"flow_minute_{i}_pL_s"] for i in range(1, 11)]
        interval = run.shared_observation_scale_interval(minute)
        row.update({f"flow_{t}s_pL_s": q for t, q in zip(run.MINUTES, minute)})
        row.update(shared_scale_lower_uL_min_per_pL_s=interval.lower_uL_min_per_pL_s,
            shared_scale_upper_uL_min_per_pL_s=interval.upper_uL_min_per_pL_s,
            minimum_minute_flow_pL_s=min(minute), maximum_minute_flow_pL_s=max(minute),
            one_smg_geometry_gate_pass=run.assess_geometry_scale(interval).geometry_gate_pass,
            shared_scale_interval_nonempty=interval.compatible,
            required_cell_count_at_minimum_minute=interval.lower_uL_min_per_pL_s / run.PL_S_CELL_TO_UL_MIN)
    by_reference = {r["root_id"]: r for r in references}
    confirmed = {(r["root_id"], r["calcium_uM"]): r for r in confirmation
                 if r["solver_label"] == run.PRODUCTION_RADAU.label}
    effective = [confirmed.get((r["root_id"], r["calcium_uM"]), r) for r in screens]
    for row in effective:
        row["mean_fold_vs_frozen_010"] = row["mean_flow_0_600_pL_s"] / by_reference[row["root_id"]]["mean_flow_0_600_pL_s"]
    groups = {0.10: references, **{ca: [r for r in effective if r["calcium_uM"] == ca] for ca in run.NEW_CALCIUM}}
    ids = {root_id: f"R{i:02}" for i, root_id in enumerate(sorted(manifest["roots"]), 1)}
    best = run.best_row(groups[0.50])
    baseline_best = by_reference[best["root_id"]]
    bdf = [r for r in confirmation if r["solver_label"] == run.PRODUCTION_BDF.label][0]
    pilot = json.loads(run.cache_path(run.primary_requests(manifest)[0]).read_text())["row"]
    screen_time = json.loads((run.OUTPUT / "screen_runtime.json").read_text())
    confirm_time = json.loads((run.OUTPUT / "confirmation_runtime.json").read_text())
    primary_wall = pilot["wall_seconds"] + screen_time["phase_wall_seconds"]
    total_wall = primary_wall + confirm_time["phase_wall_seconds"]
    prod_times = [r["wall_seconds"] for r in confirmation if r["solver_label"] == run.PRODUCTION_RADAU.label]
    approximate_nine_point_wall = math.ceil(90 / screen_time["workers"]) * (sum(prod_times) / len(prod_times))
    assessment = {
        "classification": classification,
        "all_ten_absolute_geometry_fail_at_each_new_calcium": all(not r["one_smg_geometry_gate_pass"] for r in effective),
        "root_shape_heterogeneity": {str(ca): [r["root_id"] for r in group if not r["shared_scale_interval_nonempty"]]
                                     for ca, group in groups.items()},
        "best_050_root_id": best["root_id"],
        "best_050_fraction_of_required_minute_flow": best["minimum_minute_flow_pL_s"] / run.REQUIRED_FLOW,
        "best_050_maximum_gland_flow_at_minimum_minute_uL_min": best["maximum_gland_flow_at_minimum_minute_uL_min"],
        "primary_trajectory_count": len(screens), "confirmation_trajectory_count": len(confirmation),
        "baseline_new_trajectory_count": 0,
        "all_new_trajectories_numerical_and_sustainment_pass": all(r["numerical_gate_pass"] and r["sustainment_gate_pass"]
                                                                  for r in screens + confirmation),
        "solver_crosschecks_pass": all(r["solver_crosscheck_gate_pass"] for r in confirmation
                                      if r["solver_label"] == run.PRODUCTION_BDF.label),
        "measured_primary_compute_wall_seconds": primary_wall,
        "measured_total_compute_wall_seconds": total_wall,
        "estimated_nine_point_primary_compute_wall_seconds": approximate_nine_point_wall,
        "runtime_comparison_is_estimated_not_benchmarked": True,
        "analysis_scope": "Preferred regulatory member, WT CCH_IPR, exactly two new calcium inputs",
    }
    run.write_json(run.OUTPUT / "assessment.json", assessment)
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    lines = [
        "Increasing calcium to 0.25 or 0.50 µM does not close the absolute WT secretion deficit in any of the ten frozen roots. "
        "All twenty loose trajectories and both endpoint confirmation trajectories pass the inherited numerical and sustainment gates.",
        f"The best 0.50 µM root, {ids[best['root_id']]}, reaches only {100 * assessment['best_050_fraction_of_required_minute_flow']:.2f}% "
        f"of the required minimum minute flow. Even at the unchanged generous one SMG size ceiling, its limiting minute gives "
        f"{best['maximum_gland_flow_at_minimum_minute_uL_min']:.6f} µL/min, against the 9 to 10 µL/min requirement.",
        "The comparison below uses the ten saved Task 13B production rows at 0.10 µM, the twenty new loose screen rows, "
        "and the production Radau result for the single confirmed 0.50 µM root. No 0.10 µM trajectory was recomputed. "
        "Each range is across roots. Means and cumulative flow use the inherited 122 sample grid and trapezoidal output calculation.",
    ]
    summary = []
    for ca, group in groups.items():
        gate = run.hill_activation(ca, 0.26, 1.46)
        fold = "1.000000" if ca == 0.10 else span(group, "mean_fold_vs_frozen_010")
        summary.append([f"{ca:.2f}", f"{gate:.6f}", span(group, "mean_flow_0_600_pL_s", digits=9),
                        span(group, "cumulative_flow_600_pL"), fold,
                        sum(r["shared_scale_interval_nonempty"] for r in group),
                        sum(r["one_smg_geometry_gate_pass"] for r in group)])
    lines.append(table(["Ca (µM)", "Ca gate", "Mean per cell flow (pL/s)", "600 s total (pL/cell)",
                        "Paired mean fold increase", "Minute shape passes / 10", "Absolute passes / 10"], summary))
    lines.extend([
        "Thus 0.25 µM fails the absolute gate, and 0.50 µM also fails it. The negative absolute conclusion is consistent across "
        "all ten roots, but their flow shapes differ: the five AE4NA05_P1 roots have empty minute scale intervals at both new amplitudes; "
        "the five AE4NA20_P2 roots retain nonempty minute intervals. An empty interval is a WT flow shape failure even before applying "
        "the gland size limit. It is distinct from the broader sustainment test, which all roots pass.",
        "There is no new numerical or inherited physical domain failure. The model exposes no separate source based stimulated Na, K, Cl "
        "or pH target bands, so this statement does not establish broader physiological validity. Resting calibration ranges have not been "
        "reused as new stimulated acceptance thresholds. Endpoint concentrations, pH, volume and membrane potentials are reported below.",
        "At each calcium input the following table gives the per cell flow range across roots at the unchanged minute landmarks. "
        "It retains the decline within each root instead of hiding it in an average.",
    ])
    lines.append(table(["Time (s)", "0.10 µM (pL/s)", "0.25 µM (pL/s)", "0.50 µM (pL/s)"],
        [[t] + [span(groups[ca], f"flow_{t}s_pL_s", digits=9) for ca in groups] for t in run.MINUTES]))
    lines.append("Root identifiers used in the remaining tables:")
    lines.append(table(["Root", "Frozen root ID"], [[alias, f"`{root_id}`"] for root_id, alias in ids.items()]))
    lines.append("The observation interval is the intersection of [9/q, 10/q] over all ten minute flows within each root. "
                 "Alternative roots do not share a fitted scale. Every interval fails the fixed ceiling 3088.15386. "
                 "The required cell count below is the lower scale divided by 0.00006, even when infeasible; "
                 "it is a requirement for that model, not an estimate of the actual gland count. The ceiling is 51,469,231 cells.")
    lines.append(table(["Ca (µM)", "Root", "Scale lower", "Scale upper", "Required cells (millions)", "q600/q60"],
        [[f"{ca:.2f}", ids[r["root_id"]], f"{r['shared_scale_lower_uL_min_per_pL_s']:.3f}",
          f"{r['shared_scale_upper_uL_min_per_pL_s']:.3f}", f"{r['required_cell_count_at_minimum_minute'] / 1e6:.3f}",
          f"{r['sustainment_ratio_q600_q60']:.6f}"] for ca, group in groups.items() for r in group]))
    lines.append("Endpoint intracellular states at 600 s. Potentials are apical and basolateral, in that order. "
                 "The frozen reference profile did not save endpoint potentials; those entries remain unavailable rather than being recomputed.")
    lines.append(table(["Ca (µM)", "Root", "Na (mM)", "K (mM)", "Cl (mM)", "pH", "Volume (pL)", "Va (mV)", "Vb (mV)"],
        [[f"{ca:.2f}", ids[r["root_id"]], *[f"{r[k]:.6f}" for k in
          ("endpoint_na_mM", "endpoint_k_mM", "endpoint_cl_mM", "endpoint_ph", "endpoint_volume_pL")],
          "Unavailable" if ca == 0.10 else f"{r['endpoint_v_apical_mV']:.6f}",
          "Unavailable" if ca == 0.10 else f"{r['endpoint_v_basolateral_mV']:.6f}"]
         for ca, group in groups.items() for r in group]))
    lines.extend([
        "The machine readable screen also retains sampled post stimulus flow extrema, all endpoint and trajectory domain diagnostics, "
        "each accounting residual in its inherited units, and solver counters. The minimum over t > 0 includes the inherited right limit "
        "at 0.000001 s, whose state is still the resting state; it is not the minimum of the ten minute measurements.",
        f"The predeclared 90% proximity threshold was {0.9 * run.REQUIRED_FLOW:.12f} pL/s. No root reached it at either new amplitude. "
        f"Accordingly only {ids[best['root_id']]} at 0.50 µM was confirmed with production Radau and BDF. The maximum relative "
        f"Radau versus BDF differences were {bdf['crosscheck_max_relative_state_difference']:.3g} for states, "
        f"{bdf['crosscheck_max_relative_flow_difference']:.3g} for flow and "
        f"{bdf['crosscheck_endpoint_cumulative_flow_relative_difference']:.3g} for cumulative flow, all below the inherited 0.0001 gate.",
        f"Measured numerical phase time was {primary_wall:.3f} s for the twenty case primary screen including its pilot, and "
        f"{total_wall:.3f} s including the two confirmations, on nine available worker processes. These timings exclude preparation, "
        "interpreter startup, analysis and testing. The original ninety trajectory primary screen was not executed; "
        "its runtime can only be estimated from the observed production timing. The exact primary trajectory reduction is 77.8%. "
        "The runtime report gives the assumptions and the estimated elapsed saving.",
        "A matched SMG calcium measurement remains useful for choosing the physiological input, but it is no longer the deciding "
        "measurement for the proposed amplitude only rescue within the tested range: even 0.50 µM remains far below the required flow. "
        "These two tests do not prove failure at every untested amplitude or outside the tested range. The unresolved distinction "
        "between cellular secretion capacity and the cell to gland observation map still calls for matched absolute secretion per cell. "
        "No calcium threshold was refined.",
        "This result concerns WT calcium only. It does not identify an AE4 phenotype, establish AE4 insufficiency, or assess the sealed target.",
        classification,
    ])
    (ANALYSIS / "final_answer.md").write_text("\n\n".join(lines) + "\n")
    evidence = [
        "Task 13C uses the fast screen prompt on branch `codex/task-13c-calcium-input-sweep`, read after `AGENTS.md`. "
        "The user's explicit fast screen instruction controls the two new amplitudes, preferred regulator and selective confirmation. "
        "It supersedes the broader grid, extra protocol and ensemble provisions in the older prompt.",
        f"Branch intake was `{run.BASE_COMMIT}`. The controlling final Task 13B manifest has SHA256 `{run.MANIFEST_SHA256}`. "
        "Its own hash ledger, saved dynamics profile, group gates and solver comparisons agree. All ten embedded root objects, "
        "twelve component core states, whole cell parameter payloads and AE4 parameter payloads pass their individual saved hashes.",
        "Inherited provenance discrepancies were found before simulation. The current native roots CSV is truncated at 1,048,606 bytes "
        "and contains none of the ten final roots; its hash and the current native summary hash disagree with the final manifest. "
        "The old native source hash ledger also describes a different revision. The complete final manifest contains the exact retained roots "
        "and is therefore the intake source. No missing row was regenerated and no root was solved again.",
        "The current historical dynamic runner also differs from its saved source hash. The new screen neither imports nor invokes it. "
        "Every other equation source file named in the final manifest matches its saved hash, including model, chemistry, membranes, "
        "transporters, water, regulation, numerical validation, native model construction and the NKCC wrapper. "
        "The new runner calls these unchanged implementations and verifies its reconstructed parameter objects against the frozen payloads.",
        table(["Inherited discrepant path", "Expected SHA256", "Actual SHA256"],
              [[f"`{r['path']}`", f"`{r['expected_sha256']}`", f"`{r['actual_sha256']}`"]
               for r in verification["inherited_provenance_discrepancies"]]),
        "Only the protocol's stimulated calcium changes. Resting calcium stays at 0.058 µM. The Ca gate stays "
        "C^1.46 / (0.26^1.46 + C^1.46). The NKCC normalisation reference remains 0.10 µM and its gain remains 1.75; "
        "both new inputs saturate that same existing NKCC arm. Whole cell transporter capacities, NKCC source scale 4, "
        "AE4 routing, membrane conductances, water, geometry, bath and preferred R1 regulation remain unchanged.",
        "The physical protocol remains WT CCH_IPR for 600 s, 0.3 µM CCh and 5 µM IPR, with the frozen beta step. "
        "The grid is basal zero, the right limit at 0.000001 s, then 5 s increments through 600 s. No regulatory ensemble, "
        "nearby states, alternative agonist arm, null genotype or additional calcium value is simulated.",
        "The confirmation rule was written to frozen_manifest.json before the pilot. A pass, a minimum minute flow at least "
        "90% of 9/3088.15386, or inconclusive numerical status triggers production Radau across all ten roots at that amplitude. "
        "A clear failure at 0.25 adds no confirmation; a clear failure at 0.50 adds only its best flowing root. "
        "Best means greatest minimum minute flow, then mean flow, then ascending root ID. BDF is restricted to the decisive "
        "production representative, with a second representative only for heterogeneous production classifications. "
        "The actual decision and selected requests are saved in confirmation_decisions.json.",
        "The calcium domain is motivated by the parotid context recorded in AGENTS.md and the original prompt: "
        "[Vera Sigüenza et al. (2019)](https://doi.org/10.1007/s11538-018-0534-z) and "
        "[Vera Sigüenza et al. (2020)](https://doi.org/10.1007/s11538-020-00712-3). "
        "These are not matched SMG calcium measurements and their per cell flow values are not calibration targets for this screen.",
        "The new runner never opens a target ledger, an opaque target hash path, the genotype evaluation driver or archive content. "
        "The held out AE4 phenotype remains sealed. The existing manuscript and archive are unchanged.",
        "The single independent audit checked root and parameter hashes, raw saved initial states and time grids, minute interpolation, "
        "trapezoidal totals, interval algebra, exact trajectory inventory and solver agreement. It confirmed the absolute deficit and "
        "the five versus five flow shape split. Review corrected the non calcium hash metadata to exclude the experimental calcium input "
        "and required both numerical gates in the BDF comparison. These changes required no new ODE integrations.",
        "Focused execution tests replay genuine saved solver outputs through the complete stateless worker in serial and process modes, "
        "and recompute all twenty saved trajectory summaries in both modes. This tests execution and reduction identity, not a repeated "
        "independent ODE integration. No extra trajectory is charged to this verification. Test results and file hashes are saved separately.",
    ]
    (ANALYSIS / "evidence_and_freeze.md").write_text("\n\n".join(evidence) + "\n")
    runtime = [
        "The primary screen contains exactly twenty trajectories. Its serial pilot is counted once and reused; the remaining nineteen "
        "run through a process pool with nine workers and one BLAS thread per worker. Each request runs the existing 600 s loose Radau "
        "protocol. A single 0.50 µM root then receives production Radau and BDF confirmation, for twenty two new integrations overall.",
        table(["Run", "Wall time (s)", "nfev", "njev", "nlu", "Samples"],
              [["Serial 0.25 pilot", f"{pilot['wall_seconds']:.6f}", pilot["nfev"], pilot["njev"], pilot["nlu"], pilot["output_samples"]]]
              + [[r["solver_label"] + " decisive 0.50", f"{r['wall_seconds']:.6f}", int(r["nfev"]), int(r["njev"]), int(r["nlu"]),
                  int(r["output_samples"])] for r in confirmation]),
        f"The pilot integrator itself took {pilot['integrator_wall_seconds']:.6f} s. The full pilot took "
        f"{pilot['wall_seconds']:.6f} s, far below the predeclared 120 s profiling threshold. No numerical optimisation was needed. "
        "Counter capture wraps the existing solve call and observes its return values without changing its arguments or equations.",
        f"The measured remaining screen phase took {screen_time['phase_wall_seconds']:.6f} s, including cached pilot loading and output work. "
        f"Adding the earlier pilot gives {primary_wall:.6f} s for primary computation. The confirmation phase took "
        f"{confirm_time['phase_wall_seconds']:.6f} s, giving {total_wall:.6f} s overall. These are numerical phase timings; "
        "repository retrieval, implementation, interpreter startup, documentation and tests are excluded. Per trajectory parallel wall "
        "times include CPU contention and must not be added to claim elapsed process pool time.",
        f"The original nine point primary design required ninety production Radau trajectories. It was not run. "
        f"The exact primary trajectory reduction is 90 to 20, or 77.8%; including the two confirmations, 22 trajectories are 75.6% fewer "
        f"than that old primary stage alone. A rough same machine estimate uses the observed production Radau cost "
        f"{prod_times[0]:.6f} s and ten waves of nine workers: about {approximate_nine_point_wall:.3f} s for the original primary stage, "
        f"excluding process and output overhead. On that estimate, primary computation saves about "
        f"{approximate_nine_point_wall - primary_wall:.3f} s ({100 * (1 - primary_wall / approximate_nine_point_wall):.1f}%), "
        f"or about {approximate_nine_point_wall - total_wall:.3f} s ({100 * (1 - total_wall / approximate_nine_point_wall):.1f}%) "
        "when all fast screen confirmations are counted. This is an estimate based on one production root and calcium level, "
        "not a measured counterfactual; costs can vary with calcium, root and concurrent load. No unperformed control or ensemble "
        "runs have been counted as measured savings.",
        "Recommended invocation, allowing the worker count to follow available CPU affinity and the ten root cap:\n\n"
        "```bash\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\\n"
        "python -m src.modern_full_model.run_calcium_fast_screen\n```\n\n"
        "This run used `--workers 9` with separate `--phase pilot`, `--phase screen` and `--phase confirm` calls. "
        "The default `--phase all` follows the same order. Saved completed requests are reused with hashes checked. "
        "Generate the reports from saved files with `python -m src.modern_full_model.report_calcium_fast_screen`.",
    ]
    (ANALYSIS / "runtime.md").write_text("\n\n".join(runtime) + "\n")
    print(json.dumps(assessment, indent=2))


if __name__ == "__main__":
    main()
