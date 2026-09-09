"""Task 14B: evaluate the already validated 5% AE4 states under stimulation.

This is a diagnostic on the frozen Task 14 model. It changes no model parameter,
does not solve a new resting state, and does not use the exact-zero failure as a
biological acceptance criterion. The primary quantity is total 0--600 s
secretion at 5% AE4 activity divided by the matching frozen WT total.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
from statistics import median

import numpy as np

from .genotype_evaluation import genotype_with_expression, simulate_genotype
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import CALCIUM, build_model, continuation_path, routing
from .validation import PRODUCTION_RADAU

REPO = Path(__file__).resolve().parents[2]
TASK14 = REPO / "results/14_scale_free_genotype_holdout"
OUT = REPO / "results/14B_five_percent_ae4_check"
ANALYSIS = REPO / "analysis/14B_five_percent_ae4_check"
SOURCE_COMMIT = "a2fde97e30d79601e7455d9f1d77d8d67609cb49"
EXPRESSION = 0.05


def _five_percent_core(root_id: str) -> tuple[float, ...]:
    data = json.loads(continuation_path(root_id, "AE4").read_text())
    steps = [s for s in data["steps"] if abs(float(s["expression"]) - EXPRESSION) < 1.0e-10]
    if len(steps) != 1:
        raise AssertionError(f"Expected exactly one 5% AE4 continuation step for {root_id}")
    step = steps[0]
    if not step["connected_to_previous"] or step["selected_root_id"] is None:
        raise AssertionError(f"5% AE4 state is not a connected retained root for {root_id}")
    selected = [r for r in step["roots"] if r["root_id"] == step["selected_root_id"]]
    if len(selected) != 1 or not selected[0]["passes_numerical_gate"]:
        raise AssertionError(f"Selected 5% AE4 state does not pass the numerical gate for {root_id}")
    conservation = [r for r in data["selected_conservation"]
                    if abs(float(r["expression"]) - EXPRESSION) < 1.0e-10]
    if len(conservation) != 1 or float(conservation[0]["max_dimensionless_ratio"]) > 1.0:
        raise AssertionError(f"Selected 5% AE4 state fails conservation for {root_id}")
    return tuple(float(x) for x in selected[0]["core_state"])


def _wt_totals() -> dict[tuple[str, float], float]:
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
            totals[key] = value
    return totals


def _run_one(request: tuple[str, float]) -> dict[str, object]:
    root_id, calcium = request
    manifest, _ = load_freeze()
    model = build_model(manifest, root_id, calcium)
    core = _five_percent_core(root_id)
    genotype = genotype_with_expression("AE4", EXPRESSION)
    trajectory = simulate_genotype(
        model,
        core,
        transporter="AE4",
        genotype=genotype,
        solver=PRODUCTION_RADAU,
    )
    total = float(trajectory.cumulative_flow_pL[-1])
    return {
        "root_id": root_id,
        "routing_family": routing(manifest, root_id),
        "calcium_uM": float(calcium),
        "ae4_expression": EXPRESSION,
        "ae4_5pct_total_0_600_pL": total,
        "ae4_5pct_mean_0_600_pL_s": total / 600.0,
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
        "numerical_gate_pass": bool(trajectory.numerical_gate_pass),
        "max_dimensionless_conservation_ratio": float(trajectory.max_dimensionless_conservation_ratio),
        "solver_success": bool(trajectory.success),
    }


def main() -> None:
    manifest, _ = load_freeze()
    roots = sorted(manifest["roots"])
    if len(roots) != 10:
        raise AssertionError("Task 14B requires all ten frozen roots")
    wt = _wt_totals()
    expected = {(r, float(c)) for r in roots for c in CALCIUM}
    if set(wt) != expected:
        missing = sorted(expected - set(wt))
        extra = sorted(set(wt) - expected)
        raise AssertionError(f"Frozen WT denominator matrix mismatch: missing={missing}, extra={extra}")

    requests = sorted(expected)
    workers = min(4, max(1, os.cpu_count() or 1))
    rows: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_run_one, req): req for req in requests}
        for future in as_completed(futures):
            row = future.result()
            key = (str(row["root_id"]), float(row["calcium_uM"]))
            wt_total = wt[key]
            ratio = float(row["ae4_5pct_total_0_600_pL"]) / wt_total
            row.update({
                "wt_total_0_600_pL": wt_total,
                "ae4_5pct_to_wt_total_ratio": ratio,
                "ae4_5pct_reduction_fraction": 1.0 - ratio,
                "ae4_5pct_reduction_percent": 100.0 * (1.0 - ratio),
            })
            rows.append(row)
            print(json.dumps({
                "root_id": row["root_id"],
                "calcium_uM": row["calcium_uM"],
                "ratio": ratio,
                "reduction_percent": row["ae4_5pct_reduction_percent"],
                "numerical_gate_pass": row["numerical_gate_pass"],
            }), flush=True)

    rows.sort(key=lambda r: (float(r["calcium_uM"]), str(r["root_id"])))
    if len(rows) != 30 or not all(bool(r["numerical_gate_pass"]) for r in rows):
        raise AssertionError("All thirty 5% AE4 production trajectories must pass the numerical gate")
    write_rows(OUT / "five_percent_ae4_flow.csv", rows)

    per_calcium = []
    for calcium in CALCIUM:
        group = [r for r in rows if float(r["calcium_uM"]) == float(calcium)]
        ratios = [float(r["ae4_5pct_to_wt_total_ratio"]) for r in group]
        reductions = [float(r["ae4_5pct_reduction_percent"]) for r in group]
        per_calcium.append({
            "calcium_uM": float(calcium),
            "ratio_min": min(ratios),
            "ratio_median": median(ratios),
            "ratio_max": max(ratios),
            "reduction_percent_min": min(reductions),
            "reduction_percent_median": median(reductions),
            "reduction_percent_max": max(reductions),
        })

    all_ratios = [float(r["ae4_5pct_to_wt_total_ratio"]) for r in rows]
    all_reductions = [float(r["ae4_5pct_reduction_percent"]) for r in rows]
    summary = {
        "source_task14_commit": SOURCE_COMMIT,
        "ae4_expression": EXPRESSION,
        "root_count": 10,
        "calcium_uM": list(CALCIUM),
        "trajectory_count": 30,
        "all_numerical_gates_pass": True,
        "criterion": (
            "Descriptive general-picture check only: whether 5% residual AE4 activity produces "
            "a substantial reduction in total 0--600 s secretion relative to matching frozen WT. "
            "No absolute gland-flow gate and no minute-by-minute time-shape gate is used."
        ),
        "overall_ratio_range": [min(all_ratios), max(all_ratios)],
        "overall_reduction_percent_range": [min(all_reductions), max(all_reductions)],
        "per_calcium": per_calcium,
    }
    write_json(OUT / "summary.json", summary)

    lines = [
        "# Five percent AE4 secretion check",
        "",
        "This check starts from the already retained and numerically valid Task 14 5% AE4 resting states.",
        "No model parameter, root, calcium input, regulation law or observation scale was fitted or changed.",
        "The primary endpoint is total 0--600 s secretion relative to the matching frozen WT trajectory.",
        "The inherited absolute gland scale and detailed minute-flow shape are reported nowhere as acceptance criteria.",
        "",
        f"Across all 30 root/calcium pairs the 5% AE4 to WT total secretion ratio is **{min(all_ratios):.6f} to {max(all_ratios):.6f}**,",
        f"corresponding to a **{min(all_reductions):.2f}% to {max(all_reductions):.2f}% reduction**.",
        "",
        "| Calcium (uM) | Ratio range | Median ratio | Reduction range | Median reduction |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in per_calcium:
        lines.append(
            f"| {item['calcium_uM']:.2f} | {item['ratio_min']:.6f} to {item['ratio_max']:.6f} | "
            f"{item['ratio_median']:.6f} | {item['reduction_percent_min']:.2f}% to "
            f"{item['reduction_percent_max']:.2f}% | {item['reduction_percent_median']:.2f}% |"
        )
    lines += [
        "",
        "These are model predictions at 5% residual AE4 activity, not literal genetic zero.",
        "The experimental approximately 35% knockout reduction is contextual evidence, not a fitted target or an exact acceptance band for this diagnostic.",
        "",
    ]
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    (ANALYSIS / "final_answer.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
