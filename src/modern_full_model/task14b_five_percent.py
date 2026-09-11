"""Task 14B: evaluate the already validated 5% AE4 states under stimulation.

This is a diagnostic on the frozen Task 14 model. It changes no model parameter,
does not solve a new resting state, and does not use the exact-zero failure as a
biological acceptance criterion. The primary quantity is total 0--600 s
secretion at 5% AE4 activity divided by the matching frozen WT total.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from statistics import median
import subprocess
import time

import numpy as np

from .genotype_evaluation import genotype_with_expression, simulate_genotype
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import (
    CALCIUM, build_model, continuation_path, routing, slug, trajectory_path,
    verify_contract,
)
from .validation import PRODUCTION_RADAU, sha256_file, sha256_object

REPO = Path(__file__).resolve().parents[2]
TASK14 = REPO / "results/14_scale_free_genotype_holdout"
OUT = REPO / "results/14B_five_percent_ae4_check"
ANALYSIS = REPO / "analysis/14B_five_percent_ae4_check"
SOURCE_COMMIT = "a2fde97e30d79601e7455d9f1d77d8d67609cb49"
EXPRESSION = 0.05


def _five_percent_core(root_id: str) -> tuple[float, ...]:
    data = json.loads(continuation_path(root_id, "AE4").read_text())
    if data["root_id"] != root_id or data["transporter"] != "AE4":
        raise AssertionError(f"Mismatched saved continuation for {root_id}")
    steps = [s for s in data["steps"] if abs(float(s["expression"]) - EXPRESSION) < 1.0e-10]
    if len(steps) != 1:
        raise AssertionError(f"Expected exactly one 5% AE4 continuation step for {root_id}")
    step = steps[0]
    if not step["connected_to_previous"] or step["selected_root_id"] is None:
        raise AssertionError(f"5% AE4 state is not a connected retained root for {root_id}")
    selected = [r for r in step["roots"] if r["root_id"] == step["selected_root_id"]]
    if (len(selected) != 1 or not selected[0]["passes_numerical_gate"]
            or selected[0]["normalized_jacobian_rank"] != 10):
        raise AssertionError(f"Selected 5% AE4 state does not pass the numerical gate for {root_id}")
    conservation = [r for r in data["selected_conservation"]
                    if abs(float(r["expression"]) - EXPRESSION) < 1.0e-10]
    if (len(conservation) != 1
            or not 0.0 <= float(conservation[0]["max_dimensionless_ratio"]) <= 1.0):
        raise AssertionError(f"Selected 5% AE4 state fails conservation for {root_id}")
    core = tuple(float(x) for x in selected[0]["core_state"])
    if len(core) != 12 or not all(math.isfinite(x) and x > 0.0 for x in core):
        raise AssertionError(f"Invalid saved 5% core for {root_id}")
    return core


def _verify_task14_files() -> dict[str, str]:
    """Anchor the inherited hash ledger to the user supplied scientific commit."""
    ledger = TASK14 / "final_artifact_hashes.csv"
    relative = str(ledger.relative_to(REPO))
    expected = subprocess.check_output(
        ["git", "rev-parse", f"{SOURCE_COMMIT}:{relative}"], cwd=REPO, text=True,
    ).strip()
    data = ledger.read_bytes()
    actual = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
    if actual != expected:
        raise AssertionError("Task 14 final hash ledger differs from the frozen commit")
    verified = {relative: sha256_file(ledger)}
    for row in csv.DictReader(data.decode().splitlines()):
        path = REPO / row["path"]
        if sha256_file(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise AssertionError(f"Frozen Task 14 artifact changed: {row['path']}")
        verified[row["path"]] = row["sha256"]
    return verified


def _wt_totals(contract: dict, manifest: dict) -> dict[tuple[str, float], float]:
    path = TASK14 / "blind_trajectory_summary.csv"
    totals: dict[tuple[str, float], float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (row["genotype"] != "WT" or row["solver_label"] != "production_radau"
                    or row["status"] != "COMPLETE"):
                continue
            key = (row["root_id"], float(row["calcium_uM"]))
            value = float(row["total_0_600_pL"])
            if key in totals:
                raise AssertionError(f"Duplicate WT denominator {key}")
            if not math.isfinite(value) or value <= 0.0 or row["numerical_gate_pass"] != "True":
                raise AssertionError(f"Unusable frozen WT denominator {key}")
            saved = json.loads(trajectory_path(*key, "WT", "production_radau").read_text())
            trace = saved["trajectory"]
            if (saved["root_id"] != key[0] or saved["calcium_uM"] != key[1]
                    or saved["genotype"] != "WT" or saved["status"] != "COMPLETE"
                    or saved["solver_label"] != "production_radau"
                    or trace["genotype_name"] != "WT"
                    or trace["solver"] != asdict(PRODUCTION_RADAU)
                    or not trace["success"] or not trace["numerical_gate_pass"]
                    or trace["time_s"] != contract["time_grid_s"]
                    or trace["cumulative_flow_pL"][-1] != value):
                raise AssertionError(f"Frozen WT trajectory/summary mismatch {key}")
            for field in ("whole_cell_parameters_sha256", "ae4_parameters_sha256"):
                if saved[field] != manifest["roots"][key[0]][field]:
                    raise AssertionError(f"Frozen WT parameter mismatch {key}: {field}")
            integral = float(np.trapezoid(trace["flow_pL_s"], trace["time_s"]))
            if not math.isclose(value, integral, rel_tol=1.0e-12, abs_tol=1.0e-14):
                raise AssertionError(f"Frozen WT integral mismatch {key}")
            totals[key] = value
    return totals


def _preflight() -> tuple[dict, dict, dict, dict]:
    verified = _verify_task14_files()
    contract, digest = verify_contract()
    manifest, inherited = load_freeze()
    roots = sorted(manifest["roots"])
    if len(roots) != 10 or roots != contract["root_ids"] or list(CALCIUM) != contract["calcium_uM"]:
        raise AssertionError("Task 14B requires the complete frozen root/calcium panel")
    expected = {(r, float(c)) for r in roots for c in CALCIUM}
    wt = _wt_totals(contract, manifest)
    if set(wt) != expected:
        raise AssertionError(f"Frozen WT matrix mismatch: missing={sorted(expected - set(wt))}, "
                             f"extra={sorted(set(wt) - expected)}")
    ae2 = {}
    with (TASK14 / "ae2_null_blind_predictions.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["root_id"], float(row["calcium_uM"]))
            ratio = float(row["R_total"])
            if (key in ae2 or row["status"] != "COMPLETE"
                    or row["paired_numerical_gate_pass"] != "True"
                    or not math.isfinite(ratio) or ratio <= 0.0):
                raise AssertionError(f"Invalid frozen AE2 comparison {key}")
            ae2[key] = ratio
    if set(ae2) != expected:
        raise AssertionError("Frozen AE2 comparison matrix differs from the declared panel")
    states = [{"root_id": root, "routing_family": routing(manifest, root),
               "core_state": _five_percent_core(root),
               "core_state_sha256": sha256_object(_five_percent_core(root)),
               "continuation_path": str(continuation_path(root, "AE4").relative_to(REPO))}
              for root in roots]
    required = [s["continuation_path"] for s in states] + [
        str(trajectory_path(*key, "WT", "production_radau").relative_to(REPO)) for key in expected]
    if any(path not in verified for path in required):
        raise AssertionError("A required input is missing from the frozen artifact ledger")
    frozen = {
        "source_task14_commit": SOURCE_COMMIT,
        "prediction_contract_sha256": digest,
        "verified_task14_files_sha256": verified,
        "equation_source_sha256": contract["equation_source_sha256"],
        "inherited_verification": inherited,
        "root_payload_hashes": contract["root_payload_hashes"],
        "selected_five_percent_states": states,
        "ae4_expression": EXPRESSION, "new_trajectory_count": len(expected),
        "new_resting_solves": 0, "new_wt_trajectories": 0, "new_ae2_trajectories": 0,
        "regulatory_member": contract["regulatory_member"],
        "protocol": contract["protocol"], "calcium_uM": list(CALCIUM),
        "nkcc_normalisation_frozen": contract["nkcc_normalisation_frozen"],
        "production_solver": asdict(PRODUCTION_RADAU),
        "time_grid_s": contract["time_grid_s"],
        "conservation_tolerances": contract["conservation_tolerances"],
        "runner_sha256": sha256_file(Path(__file__)),
        "requests": [{"root_id": root, "calcium_uM": ca,
                      "wt_total_0_600_pL": wt[root, ca],
                      "wt_path": str(trajectory_path(root, ca, "WT", "production_radau").relative_to(REPO))}
                     for root, ca in sorted(expected, key=lambda k: (k[1], k[0]))],
    }
    return manifest, wt, ae2, frozen


def _output_path(root_id: str, calcium: float) -> Path:
    return OUT / "trajectories" / f"{slug(root_id)}_ca{calcium:.2f}_AE4_5pct_production_radau.json"


def _finite_json(value):
    """Preserve a failed trajectory while representing nonfinite diagnostics as null."""
    if isinstance(value, np.ndarray):
        return _finite_json(value.tolist())
    if isinstance(value, dict):
        return {key: _finite_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _run_one(request: tuple[str, float]) -> dict[str, object]:
    root_id, calcium = request
    manifest, _ = load_freeze()
    model = build_model(manifest, root_id, calcium)
    core = _five_percent_core(root_id)
    genotype = genotype_with_expression("AE4", EXPRESSION)
    wt_path = trajectory_path(root_id, calcium, "WT", "production_radau")
    grid = json.loads(wt_path.read_text())["trajectory"]["time_s"]
    started = time.perf_counter()
    trajectory = simulate_genotype(
        model,
        core,
        transporter="AE4",
        genotype=genotype,
        solver=PRODUCTION_RADAU,
        time_s=grid,
    )
    complete = bool(trajectory.success and np.array_equal(trajectory.time_s, grid))
    valid = bool(complete and trajectory.numerical_gate_pass)
    total = float(trajectory.cumulative_flow_pL[-1]) if valid else None
    row = {
        "root_id": root_id,
        "routing_family": routing(manifest, root_id),
        "calcium_uM": float(calcium),
        "ae4_expression": EXPRESSION,
        "status": "COMPLETE" if valid else "NUMERICAL_FAILURE",
        "solver_label": PRODUCTION_RADAU.label,
        "message": trajectory.message,
        "ae4_5pct_total_0_600_pL": total,
        "ae4_5pct_mean_0_600_pL_s": total / 600.0 if valid else None,
        "ae4_5pct_rest_na_mM": float(trajectory.cell_na_mM[0]),
        "ae4_5pct_rest_k_mM": float(trajectory.cell_k_mM[0]),
        "ae4_5pct_rest_cl_mM": float(trajectory.cell_cl_mM[0]),
        "ae4_5pct_rest_ph": float(trajectory.cell_ph[0]),
        "ae4_5pct_endpoint_na_mM": float(trajectory.cell_na_mM[-1]),
        "ae4_5pct_endpoint_k_mM": float(trajectory.cell_k_mM[-1]),
        "ae4_5pct_endpoint_cl_mM": float(trajectory.cell_cl_mM[-1]),
        "ae4_5pct_endpoint_ph": float(trajectory.cell_ph[-1]),
        "positive_core": bool(trajectory.positive_core),
        "all_flow_nonnegative": bool(trajectory.all_flow_nonnegative),
        "numerical_gate_pass": valid,
        "max_dimensionless_conservation_ratio": float(trajectory.max_dimensionless_conservation_ratio),
        "solver_success": bool(trajectory.success),
        "complete_matching_time_grid": complete,
        "time_sample_count": len(trajectory.time_s),
        "trajectory_path": str(_output_path(root_id, calcium).relative_to(REPO)),
        "max_abs_cell_charge_fmol": float(np.max(np.abs(
            trajectory.states[0] + trajectory.states[1] - trajectory.states[2]
            - trajectory.states[4] - model.parameters.geometry.fixed_cell_anion_equivalents_fmol))),
        "max_abs_lumen_charge_fmol": float(np.max(np.abs(
            trajectory.states[6] + trajectory.states[7] - trajectory.states[8] - trajectory.states[10]))),
    }
    for key, value in trajectory.max_abs_conservation_residuals.items():
        row[f"max_abs_residual__{key}"] = float(value)
    for compartment in ("cell", "lumen"):
        values = getattr(trajectory, f"{compartment}_volume_pL")
        row[f"ae4_5pct_rest_{compartment}_volume_pL"] = float(values[0])
        row[f"ae4_5pct_endpoint_{compartment}_volume_pL"] = float(values[-1])
    write_json(_output_path(root_id, calcium), _finite_json({
        "source_task14_commit": SOURCE_COMMIT,
        "root_id": root_id, "routing_family": row["routing_family"], "calcium_uM": calcium,
        "genotype": asdict(genotype), "status": row["status"],
        "selected_core_state_sha256": sha256_object(core),
        "whole_cell_parameters_sha256": sha256_object(model.parameters),
        "ae4_parameters_sha256": sha256_object(model.ae4_parameters),
        "wt_denominator_path": str(wt_path.relative_to(REPO)),
        "wt_denominator_sha256": sha256_file(wt_path),
        "trajectory": asdict(trajectory), "numerical_diagnostics": row,
        "wall_seconds": time.perf_counter() - started,
    }))
    return _finite_json(row)


def _pair(row: dict, wt_total: float, ae2_ratio: float) -> dict:
    valid = row["numerical_gate_pass"]
    ratio = float(row["ae4_5pct_total_0_600_pL"]) / wt_total if valid else None
    reduction = 100.0 * (1.0 - ratio) if valid else None
    ae2_reduction = 100.0 * (1.0 - ae2_ratio)
    return {**row, "wt_total_0_600_pL": wt_total,
            "ae4_5pct_to_wt_total_ratio": ratio,
            "ae4_5pct_reduction_fraction": 1.0 - ratio if valid else None,
            "ae4_5pct_reduction_percent": reduction,
            "ae2_to_wt_total_ratio": ae2_ratio,
            "ae2_reduction_percent": ae2_reduction,
            "ae4_minus_ae2_reduction_percentage_points": reduction - ae2_reduction if valid else None,
            "ae4_to_ae2_effect_magnitude_ratio": abs(reduction / ae2_reduction)
            if valid and ae2_reduction != 0.0 else None}


def _statistics(rows: list[dict]) -> dict:
    valid = [r for r in rows if r["numerical_gate_pass"]]
    result = {"pair_count": len(rows), "valid_count": len(valid),
              "all_pairs_reduce_secretion": bool(valid) and len(valid) == len(rows)
              and all(r["ae4_5pct_reduction_percent"] > 0.0 for r in valid)}
    for name, field in (("ratio", "ae4_5pct_to_wt_total_ratio"),
                        ("reduction_percent", "ae4_5pct_reduction_percent"),
                        ("ae2_reduction_percent", "ae2_reduction_percent"),
                        ("ae4_to_ae2_effect_magnitude_ratio", "ae4_to_ae2_effect_magnitude_ratio")):
        source = rows if name == "ae2_reduction_percent" else valid
        values = [r[field] for r in source if r[field] is not None]
        for suffix, operation in (("min", min), ("median", median), ("max", max)):
            result[f"{name}_{suffix}"] = operation(values) if values else None
    return result


def _summarise(rows: list[dict]) -> dict:
    overall = _statistics(rows)
    families = sorted({r["routing_family"] for r in rows})
    return {
        "source_task14_commit": SOURCE_COMMIT, "ae4_expression": EXPRESSION,
        "root_count": len({r["root_id"] for r in rows}), "calcium_uM": list(CALCIUM),
        "trajectory_count": len(rows), "valid_trajectory_count": overall["valid_count"],
        "all_numerical_gates_pass": overall["valid_count"] == len(rows),
        "all_pairs_reduce_secretion": overall["all_pairs_reduce_secretion"],
        "criterion": "Descriptive total 0 to 600 s secretion comparison. Time traces and the inherited "
                     "absolute gland scale are diagnostics only; experimental magnitude is context only.",
        "overall_ratio_range": [overall["ratio_min"], overall["ratio_max"]],
        "overall_ratio_median": overall["ratio_median"],
        "overall_reduction_percent_range": [overall["reduction_percent_min"], overall["reduction_percent_max"]],
        "overall_reduction_percent_median": overall["reduction_percent_median"],
        "overall": overall,
        "per_calcium": [{"calcium_uM": ca, **_statistics([r for r in rows if r["calcium_uM"] == ca])}
                        for ca in CALCIUM],
        "per_routing_family": [{"routing_family": family,
                                **_statistics([r for r in rows if r["routing_family"] == family])}
                               for family in families],
        "per_calcium_and_routing_family": [
            {"calcium_uM": ca, "routing_family": family,
             **_statistics([r for r in rows if r["calcium_uM"] == ca and r["routing_family"] == family])}
            for ca in CALCIUM for family in families],
        "numerical_maxima": {key: max((r[key] for r in rows if r.get(key) is not None), default=None)
                             for key in ("max_dimensionless_conservation_ratio", "max_abs_cell_charge_fmol",
                                         "max_abs_lumen_charge_fmol")},
        "failed_cases": [{"root_id": r["root_id"], "calcium_uM": r["calcium_uM"],
                          "status": r["status"], "message": r.get("message")}
                         for r in rows if not r["numerical_gate_pass"]],
    }


def _write_report(rows: list[dict], summary: dict) -> None:
    number = lambda value, digits=6: "unavailable" if value is None else f"{value:.{digits}f}"
    overall = summary["overall"]
    lines = [
        "# Five percent AE4 secretion check",
        "",
        f"Across all 30 root/calcium pairs, the 5% AE4 to WT total secretion ratio ranges from "
        f"**{number(overall['ratio_min'])} to {number(overall['ratio_max'])}** "
        f"(median **{number(overall['ratio_median'])}**). This corresponds to a "
        f"**{number(overall['reduction_percent_min'], 2)}% to {number(overall['reduction_percent_max'], 2)}% reduction** "
        f"(median **{number(overall['reduction_percent_median'], 2)}%**).",
        "",
        f"Numerically valid trajectories: **{summary['valid_trajectory_count']}/30**. "
        f"Every root at every calcium reduces secretion: **{summary['all_pairs_reduce_secretion']}**. "
        "Any failed cases remain in the tables with unavailable primary ratios.",
        "",
        "| Calcium (uM) | Ratio range | Median ratio | Reduction range | Median reduction |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in summary["per_calcium"]:
        lines.append(
            f"| {item['calcium_uM']:.2f} | {number(item['ratio_min'])} to {number(item['ratio_max'])} | "
            f"{number(item['ratio_median'])} | {number(item['reduction_percent_min'], 2)}% to "
            f"{number(item['reduction_percent_max'], 2)}% | {number(item['reduction_percent_median'], 2)}% |"
        )
    lines += [
        "",
        f"The saved Task 14 AE2 percentage reductions range from "
        f"{number(overall['ae2_reduction_percent_min'])}% to {number(overall['ae2_reduction_percent_max'])}% "
        f"(median {number(overall['ae2_reduction_percent_median'])}%). Negative reductions denote a small increase. "
        f"The matched AE4 effect magnitudes are {number(overall['ae4_to_ae2_effect_magnitude_ratio_min'], 1)} "
        f"to {number(overall['ae4_to_ae2_effect_magnitude_ratio_max'], 1)} times the AE2 effects.",
        "",
        "Both routing families are retained. Their complete calcium specific summaries follow; "
        "medians describe the five retained roots in each family and are not uncertainty intervals.",
        "",
        "| Routing family | Calcium (uM) | Ratio range | Median ratio | Reduction range (%) | Median reduction (%) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in summary["per_calcium_and_routing_family"]:
        lines.append(f"| {item['routing_family']} | {item['calcium_uM']:.2f} | "
                     f"{number(item['ratio_min'])} to {number(item['ratio_max'])} | {number(item['ratio_median'])} | "
                     f"{number(item['reduction_percent_min'], 2)} to {number(item['reduction_percent_max'], 2)} | "
                     f"{number(item['reduction_percent_median'], 2)} |")
    lines += [
        "",
        "These are model predictions at 5% residual AE4 activity, a near complete functional loss state. "
        "They are not literal genetic zero. The experimental approximately 35% AE4 null reduction is "
        "contextual evidence, not a fitted target or an exact acceptance band. Detailed experimental "
        "time signatures are not acceptance criteria. The inherited absolute gland scale discrepancy "
        "remains recorded in Task 14 and is nonblocking for these ratios.",
        "",
        f"Scientific source: `{SOURCE_COMMIT}`. The saved selected connected 5% resting state was used "
        "without another resting solve for every root. Exactly 30 new production Radau trajectories "
        "use R1, frozen NKCC1 normalisation, CCh + IPR, the three declared calcium inputs, "
        "the inherited tolerances and the original 602 point WT grid over 0 to 600 s. "
        "Each denominator is the matching saved Task 14 WT production total, checked against the "
        "raw trajectory and the commit anchored artifact hashes. Model and root hashes were "
        "verified before and after execution. No model parameter, root, regulation law or observation "
        "scale was fitted or changed. No root or calcium was selected by phenotype agreement.",
        "",
        f"Maximum dimensionless conservation residual ratio: "
        f"{number(summary['numerical_maxima']['max_dimensionless_conservation_ratio'], 9)} "
        "(inherited numerical limit 1). Full solver messages, positivity, charge, current, carbon, "
        "water, ionic states, pH and flow/cumulative traces are saved alongside the pair table.",
        "",
        "| Root | Calcium (uM) | AE4 5% / WT | Reduction (%) | AE2 reduction (%) | Numerical status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(f"| {row['root_id']} | {row['calcium_uM']:.2f} | "
                     f"{number(row['ae4_5pct_to_wt_total_ratio'])} | "
                     f"{number(row['ae4_5pct_reduction_percent'])} | "
                     f"{number(row['ae2_reduction_percent'])} | {row['status']} |")
    lines += [
        "",
        "Machine readable results: `results/14B_five_percent_ae4_check/five_percent_ae4_flow.csv`, "
        "`summary.json`, `frozen_inputs.json`, `postrun_integrity.json` and `trajectories/`.",
        "",
    ]
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    (ANALYSIS / "final_answer.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    manifest, wt, ae2, frozen = _preflight()
    requests = [(r["root_id"], r["calcium_uM"]) for r in frozen["requests"]]
    if (OUT / "frozen_inputs.json").exists() or any(_output_path(*key).exists() for key in requests):
        raise FileExistsError("Task 14B output already exists; do not silently rerun or replace trajectories")
    workers = min(4, max(1, os.cpu_count() or 1))
    frozen["workers"] = workers
    frozen["numerical_library_threads"] = {key: os.environ.get(key) for key in
        ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}
    write_json(OUT / "frozen_inputs.json", frozen)
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_run_one, req): req for req in requests}
        for future in as_completed(futures):
            key = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {"root_id": key[0], "calcium_uM": key[1],
                       "routing_family": routing(manifest, key[0]), "ae4_expression": EXPRESSION,
                       "status": "EXECUTION_ERROR", "numerical_gate_pass": False,
                       "ae4_5pct_total_0_600_pL": None,
                       "message": f"{type(error).__name__}: {error}"}
                write_json(_output_path(*key), row)
            row = _pair(row, wt[key], ae2[key])
            rows.append(row)
            rows.sort(key=lambda r: (r["calcium_uM"], r["root_id"]))
            write_rows(OUT / "five_percent_ae4_flow.csv", rows)
            print(json.dumps({name: row[name] for name in (
                "root_id", "calcium_uM", "ae4_5pct_to_wt_total_ratio",
                "ae4_5pct_reduction_percent", "numerical_gate_pass")}), flush=True)
    if len(rows) != 30 or {(r["root_id"], r["calcium_uM"]) for r in rows} != set(wt):
        raise AssertionError("Task 14B results must retain every declared root/calcium pair")
    _, post_wt, post_ae2, post_frozen = _preflight()
    unchanged = (post_wt == wt and post_ae2 == ae2
                 and all(post_frozen[key] == frozen[key] for key in post_frozen))
    write_json(OUT / "postrun_integrity.json", {
        "source_task14_commit": SOURCE_COMMIT, "all_frozen_inputs_unchanged": unchanged,
        "verified_task14_file_count": len(post_frozen["verified_task14_files_sha256"]),
        "equation_file_count": len(post_frozen["equation_source_sha256"]),
        "new_trajectory_count": len(rows), "new_resting_solves": 0,
        "new_wt_trajectories": 0, "new_ae2_trajectories": 0,
    })
    if not unchanged:
        raise AssertionError("Frozen inputs changed during Task 14B")
    summary = _summarise(rows)
    write_json(OUT / "summary.json", summary)
    _write_report(rows, summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
