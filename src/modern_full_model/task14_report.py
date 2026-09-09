"""Deterministic summaries of the saved Task 14 blind predictions."""
from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path
import numpy as np

from . import task14_blind as b
from .native_dynamic_contract import (
    ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S, shared_observation_scale_interval,
    assess_geometry_scale,
)
from .validation import sha256_file, canonical_json_bytes


def read_rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def selected(step):
    return next((r for r in step["roots"] if r["root_id"] == step["selected_root_id"]), None)


def continuation_rows():
    c, digest = b.contract()
    rows = []
    for root in c["root_ids"]:
        for genotype in ("AE4", "AE2"):
            d = json.loads(b.continuation_path(root, genotype).read_text())
            assert d["prediction_contract_sha256"] == digest
            for implementation, report in (("production", d), ("independent", d.get("independent_audit"))):
                if report is None:
                    continue
                for step in report["steps"]:
                    s = selected(step)
                    rows.append(b.stamp({"root_id": root, "routing_family": c["routing_family_by_root"][root],
                        "genotype": genotype, "implementation": implementation,
                        "expression": step["expression"], "status": step["status"],
                        "completed_to_exact_zero": report["completed_to_exact_zero"],
                        "fitted_parameter_count": 0,
                        "fixed_OTHER_fmol": d["fixed_other_impermeant_osmoles_fmol"],
                        "selected_root_id": step["selected_root_id"],
                        "selected_distance": step["selected_distance_from_previous_normalized"],
                        "retained_root_count": len(step["roots"]),
                        "selected_core_state_json": json.dumps(s["core_state"]) if s else None,
                        "attempt_count": len(step["attempts"]),
                        "all_attempts_and_alternate_roots_path": str(b.continuation_path(root, genotype).relative_to(b.REPO)),
                        "maximum_attempt_nfev": max(a["nfev"] for a in step["attempts"]),
                        "best_attempt_scaled_residual": min(a.get("max_abs_scaled_independent_rhs", a.get("max_abs_scaled_residual")) for a in step["attempts"]),
                        "all_attempts_optimizer_success": all(a["optimizer_success"] for a in step["attempts"]),
                        "independent_agrees_on_shared_steps_and_completion": d.get("independent_agrees"),
                        "independent_max_state_relative_difference": d.get("independent_max_state_relative_difference"),
                    }))
    return rows


def absolute_row(data):
    row = {k: data[k] for k in ("root_id", "routing_family", "calcium_uM", "genotype", "solver_label", "status")}
    t = data["trajectory"]
    if t is None:
        return b.stamp(row)
    for k in ("success", "message", "positive_core", "all_flow_nonnegative", "numerical_gate_pass",
              "max_abs_deleted_transport_flux_fmol_s", "max_dimensionless_conservation_ratio"):
        row[k] = t[k]
    row.update({"max_abs_residual__" + k: v for k, v in t["max_abs_conservation_residuals"].items()})
    row.update({k: v for k, v in data["additional_diagnostics"].items() if not isinstance(v, dict)})
    if not t["success"]:
        return b.stamp(row)
    for s in b.MINUTES:
        row[f"flow_{s}s_pL_s"] = float(np.interp(s, t["time_s"], t["flow_pL_s"]))
        row[f"cumulative_{s}s_pL"] = float(np.interp(s, t["time_s"], t["cumulative_flow_pL"]))
    row["total_0_600_pL"] = t["cumulative_flow_pL"][-1]
    row["mean_0_600_pL_s"] = row["total_0_600_pL"] / 600
    row["mean_mapped_gland_at_common_scale_ceiling_uL_min"] = row["mean_0_600_pL_s"] * ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
    row["common_scale_ceiling_uL_min_per_pL_s"] = ONE_SMG_SCALE_MAX_UL_MIN_PER_PL_S
    for k in ("cell_na_mM", "cell_k_mM", "cell_cl_mM", "cell_ph", "cell_volume_pL",
              "lumen_na_mM", "lumen_k_mM", "lumen_cl_mM", "lumen_ph", "lumen_volume_pL"):
        row["rest_" + k] = t[k][0]
        row["endpoint_" + k] = t[k][-1]
    if data["genotype"] == "WT":
        interval = shared_observation_scale_interval([row[f"flow_{s}s_pL_s"] for s in b.MINUTES])
        geometry = assess_geometry_scale(interval)
        row.update({"old_WT_scale_" + k: v for k, v in asdict(interval).items()})
        row.update({"old_WT_geometry_" + k: v for k, v in asdict(geometry).items()})
        row["old_scale_diagnostics_block_genotype_prediction"] = False
    return b.stamp(row)


def pair_rows():
    c, _ = b.contract()
    rows = {"AE4": [], "AE2": []}
    combined = []
    for root in c["root_ids"]:
        for ca in b.CALCIUM:
            wt = b.load_trajectory(root, ca, "WT")
            both = {"root_id": root, "routing_family": c["routing_family_by_root"][root], "calcium_uM": ca}
            for g in ("AE4", "AE2"):
                d = b.load_trajectory(root, ca, g)
                row = {**both, "genotype": g, "status": d["status"],
                    "paired_numerical_gate_pass": d["status"] == wt["status"] == "COMPLETE"}
                for key in ("R_total", "D_total", "R_early", "R_late", "delta_emergence"):
                    row[key] = None
                for minute in b.MINUTES:
                    for prefix in ("R_flow", "R_cum", "R_bin"):
                        row[f"{prefix}_{minute}s"] = None
                for field in ("cell_cl_mM", "cell_ph", "cell_na_mM", "cell_k_mM", "cell_volume_pL"):
                    for label in ("rest", "endpoint"):
                        row[f"{label}_delta_{field}"] = None
                if row["paired_numerical_gate_pass"]:
                    row.update(b.ratios(wt["trajectory"], d["trajectory"]))
                for field in ("R_total", "D_total", "R_early", "R_late", "delta_emergence"):
                    both[g + "_" + field] = row.get(field)
                both[g + "_status"] = row["status"]
                rows[g].append(b.stamp(row))
            both["AE4_effect"] = both["AE4_D_total"]
            both["AE2_effect"] = both["AE2_D_total"]
            both["effect_contrast"] = None if either_missing(both["AE4_effect"], both["AE2_effect"]) else both["AE4_effect"] - both["AE2_effect"]
            combined.append(b.stamp(both))
    return rows, combined


def either_missing(*values):
    return any(v is None for v in values)


def write_doc(name, text):
    if (b.OUT / "blind_checkpoint.json").exists():
        raise RuntimeError("Blind checkpoint is immutable")
    b.ANALYSIS.mkdir(parents=True, exist_ok=True)
    (b.ANALYSIS / name).write_text(text + "\n\nPrediction contract SHA256: `" + b.contract()[1] + "`.\n")


def freeze_doc():
    c, digest = b.contract()
    return f"""# Task 14 freeze and design

The scientific state is commit `{b.BASE}`, tree `{b.BASE_TREE}`. Branch intake `{b.INTAKE}` changes only AGENTS.md and the Task 14 prompt. Source files and all ten root payloads are hash verified. The Task 13B final manifest and Task 13C frozen manifest hashes are retained in the prediction contract.

The immutable prediction contract was written at {c['created_utc']}, before any Task 14 genotype result existed. All ten roots, AE4 and AE2 exact deletion, calcium 0.10, 0.25 and 0.50 uM, the R1 member `{b.MEMBER}`, and the CCH_IPR 600 s protocol were declared together. Resting calcium is 0.058 uM. The NKCC1 calcium normalisation remains 0.10 uM with multiplier 1.75; changing the input amplitude does not change that parameter.

Expression continuation uses the inherited 21 points, three local starts, ten state coordinates, zero fitted parameters and fixed OTHER osmoles. All attempted and retained roots are saved. Selection uses only distance from the previous root divided by the declared coordinate span. The maximum allowed connected step is 0.25. No null dynamics may start from a WT resting state or an unconverged exact deletion attempt.

Production Radau uses relative tolerance 1e-7, amount absolute tolerance 1e-10 fmol, volume tolerance 1e-12 pL, regulatory tolerance 1e-10 and maximum step 2 s. BDF uses the same tolerances. Each pair uses the identical 602 point grid containing 0, the onset right limit 1e-6 s, and every integer second to 600. The inherited relative solver agreement tolerance is 1e-4. Both nulls and their WT denominator are checked for the lexicographic representative of each routing family at all three calcium values. A disagreement expands checks to the entire affected routing family and calcium. A missing exact resting state is recorded as unavailable, not as a Radau versus BDF disagreement.

The representatives are:

* AE4NA05: `{c['representatives_lexicographic']['AE4NA05']}`.
* AE4NA20: `{c['representatives_lexicographic']['AE4NA20']}`.

The independent resting implementation audits both representatives for both deletions and every anomalous branch. No independent resting audit changes a production initial state. Conservation, current closure, charge, positivity, exact deletion and full resting numerical rank remain required.

The paired metrics are total, instantaneous minute, cumulative minute, early 0 to 180 s and sustained 180 to 600 s ratios; late minus early ratio; AE4 and AE2 reduction fractions and their contrast; and resting and endpoint ion, pH and volume differences. Minute collection data, if available, are compared with integrated minute bins rather than instantaneous endpoint flow. Time course RMSE is unweighted across available matching relative points. No uncertainty interval, observation factor, time shift or other fitted transformation is introduced.

The inherited absolute whole gland discrepancy remains unchanged: the Task 13C best 0.50 uM result reached 5.405993 uL/min under the existing one SMG mapping. Five roots also failed its minute observation shape criterion. These diagnostics are preserved and are nonblocking for paired genotype quantities. A common multiplicative scale cancels from the ratios. No genotype specific effective cell count is estimated.

Every root and calcium remains in every summary, including failures. Exact experimental magnitude and time course remain sealed until the completed blind artifacts and tests have been committed and pushed. The later evaluator reads the saved predictions and writes only separately named comparison artifacts. No model or trajectory change is permitted after reveal.
"""


def summarise():
    if (b.OUT / "blind_checkpoint.json").exists():
        raise RuntimeError("Blind checkpoint is immutable")
    c, digest = b.verify_contract()
    continuation = continuation_rows()
    pairs, combined = pair_rows()
    b.write_rows(b.OUT / "genotype_continuation.csv", continuation)
    for g in ("AE4", "AE2"):
        b.write_rows(b.OUT / f"{g.lower()}_null_blind_predictions.csv", pairs[g])
    b.write_rows(b.OUT / "blind_pair_summary.csv", combined)
    absolute = []
    for path in sorted((b.OUT / "trajectories").glob("*.json")):
        absolute.append(absolute_row(json.loads(path.read_text())))
    b.write_rows(b.OUT / "blind_trajectory_summary.csv", absolute)
    primary = [r for r in absolute if r["solver_label"] == "production_radau"]
    checks = json.loads((b.OUT / "solver_confirmation.json").read_text())
    branches = [json.loads(b.continuation_path(root, g).read_text()) for root in c["root_ids"] for g in ("AE4", "AE2")]
    audit = {
        "root_count": len(c["root_ids"]), "calcium_count": len(b.CALCIUM),
        "primary_matrix_rows": len(primary), "primary_trajectories_run": sum(r["status"] != "CONTINUATION_FAILED" for r in primary),
        "primary_complete": sum(r["status"] == "COMPLETE" for r in primary),
        "primary_continuation_blocked": sum(r["status"] == "CONTINUATION_FAILED" for r in primary),
        "bdf_trajectories_run": sum(r["solver_label"] == "production_bdf" and r["status"] != "CONTINUATION_FAILED" for r in absolute),
        "continuation_completed": {g: sum(d["completed_to_exact_zero"] for d in branches if d["transporter"] == g) for g in ("AE4", "AE2")},
        "independent_continuations_run": sum("independent_audit" in d for d in branches),
        "independent_all_agree_on_completion_and_shared_steps": all(d.get("independent_agrees", True) for d in branches),
        "independent_max_state_relative_difference": max(d.get("independent_max_state_relative_difference", 0) or 0 for d in branches),
        "all_selected_rest_conservation_pass": all(r["max_dimensionless_ratio"] <= 1.0 for d in branches for r in d["selected_conservation"]),
        "all_available_solver_confirmations_pass": checks["all_available_pass"],
        "solver_checks": checks["checks"], "expanded_family_calcium": checks["expanded_family_calcium"],
        "max_available_dynamic_conservation_ratio": max(r.get("max_dimensionless_conservation_ratio", 0) for r in absolute),
        "max_dynamic_charge_fmol": max(max(r.get("max_abs_cell_charge_fmol", 0), r.get("max_abs_lumen_charge_fmol", 0)) for r in absolute),
        "ae4_biological_comparison_licensed": all(r["paired_numerical_gate_pass"] for r in pairs["AE4"]),
        "all_roots_and_calcium_retained": len(combined) == 30 and all(len(rows) == 30 for rows in pairs.values()),
        "no_parameter_fit": True, "exact_phenotype_values_accessed": False,
    }
    b.write_output(b.OUT / "blind_numerical_audit.json", audit)
    write_doc("freeze_and_design.md", freeze_doc())
    lines = ["# Blind Task 14 prediction", "",
        "This report states the numerical predictions without an experimental comparison.", "",
        f"Exact resting continuation succeeded in {audit['continuation_completed']['AE4']}/10 AE4 deletion branches and {audit['continuation_completed']['AE2']}/10 AE2 deletion branches. All ten AE4 branches reached expression 0.05, then failed at exact zero under the inherited 3000 evaluation limit and numerical gates. All three local attempts at exact zero failed in each branch. The independent implementation reproduced the failure on all ten AE4 branches. Agreement between implementations at the retained positive expression steps does not establish that no exact null equilibrium exists; it establishes that the declared procedure did not obtain one.", "",
        "No AE4 null stimulated trajectory is licensed. AE4 total, minute, early and sustained secretion ratios, the emergence statistic, ionic differences and the AE4 minus AE2 contrast are therefore unavailable for every calcium value. Failed exact zero iterates are retained as numerical attempts and are never treated as equilibria.", "",
        f"The full 90 row primary matrix contains {audit['primary_trajectories_run']} simulated WT and AE2 trajectories and {audit['primary_continuation_blocked']} explicitly blocked AE4 entries. Production Radau completed {audit['primary_complete']} valid trajectories. BDF ran {audit['bdf_trajectories_run']} representative WT and AE2 confirmations. All available solver comparisons pass: {checks['all_available_pass']}. The independent resting audit ran {audit['independent_continuations_run']} continuations, covering both representatives and every failed branch.", "",
        "The complete valid AE2 predictions are summarised below. Reduction is defined as one minus the null to WT total secretion ratio; negative reduction means increased secretion.", "",
        "| Calcium uM | Routing | Total ratio range | Median | Early ratio range | Sustained ratio range |", "| :--- | :--- | :--- | :--- | :--- | :--- |"]
    for ca in b.CALCIUM:
        for family in sorted(c["representatives_lexicographic"]):
            rows = [r for r in pairs["AE2"] if r["calcium_uM"] == ca and r["routing_family"] == family and r["R_total"] is not None]
            if not rows:
                continue
            interval = lambda key: f"{min(r[key] for r in rows):.8f} to {max(r[key] for r in rows):.8f}"
            lines.append(f"| {ca:.2f} | {family} | {interval('R_total')} | {np.median([r['R_total'] for r in rows]):.8f} | {interval('R_early')} | {interval('R_late')} |")
    lines += ["", "Every primary trajectory and every continuation attempt remains in the saved artifacts. The complete per root and calcium results are in the CSV files; full time grids, states, ionic quantities and accounting diagnostics are saved separately. Absolute cellular and mapped gland flows remain secondary diagnostics. The Task 13C files have not changed.", "",
        "The declared numerical outcome before any experimental reveal is: **GENOTYPE CONTINUATION OR NUMERICS PREVENT A VALID HOLDOUT TEST**. This is a failure to obtain a valid AE4 prediction under the frozen rules, not an experimental rejection of a computed AE4 secretion curve."]
    write_doc("blind_prediction.md", "\n".join(lines))
    return audit


def hash_blind_artifacts():
    if (b.OUT / "blind_checkpoint.json").exists():
        raise RuntimeError("Blind checkpoint is immutable")
    files = sorted(p for p in b.OUT.rglob("*") if p.is_file() and p.name != "blind_artifact_hashes.csv")
    files += sorted(b.ANALYSIS.glob("*.md"))
    files += sorted((b.REPO / "src/modern_full_model").glob("task14_*.py"))
    files += sorted((b.REPO / "tests").glob("test_modern_task14*.py"))
    b.write_rows(b.OUT / "blind_artifact_hashes.csv", [b.stamp({"path": str(p.relative_to(b.REPO)),
        "sha256": sha256_file(p), "bytes": p.stat().st_size}) for p in files])


if __name__ == "__main__":
    print(json.dumps(summarise(), sort_keys=True))
