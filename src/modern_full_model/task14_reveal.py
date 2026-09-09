"""One source evaluation followed by deterministic replay from retained records.

No model module or solver is imported. The blind hash ledger must be intact.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

from . import task14_compare as compare

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/14_scale_free_genotype_holdout"
ANALYSIS = REPO / "analysis/14_scale_free_genotype_holdout"
CLASSIFICATION = "GENOTYPE CONTINUATION OR NUMERICS PREVENT A VALID HOLDOUT TEST"


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def source_records(snapshot, checkpoint):
    """Parse each retrieved source exactly once and retain just the used rows."""
    if (OUT / "reveal_manifest.json").exists():
        raise RuntimeError("Sources have already been evaluated; use retained records for replay")
    if snapshot["blind_checkpoint_sha"] != checkpoint["blind_checkpoint_sha"] or not checkpoint["remote_ref_verified"]:
        raise ValueError("A verified pushed blind checkpoint must precede source evaluation")
    used = {
        "heldout_targets.csv": {
            "H15_AE4_NULL_TOTAL_10MIN_DEFICIT", "H15_AE4_NULL_TOTAL_10MIN_RATIO",
            "H15_AE4_NULL_EARLY_TIMECOURSE", "H15_AE4_NULL_SUSTAINED_TIMECOURSE",
            "H15_AE4_NULL_LATE_POINTWISE"},
        "validation_targets.csv": {
            "V15_AE2_SECRETION", "V15_AE4_CL_REST_NULL", "V15_AE4_PHI_REST_NULL",
            "V15_AE2_CL_REST_CONTROL", "V15_AE2_CL_REST_NULL",
            "V15_AE2_PHI_REST_CONTROL", "V15_AE2_PHI_REST_NULL"},
    }
    sources, records = [], {}
    for source in snapshot["sources"]:
        raw = source["content"].encode()
        blob = hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()
        if blob != source["git_blob_sha"]:
            raise ValueError("Retrieved source is not the recorded Git blob")
        rows = list(csv.DictReader(io.StringIO(source["content"])))
        selected = [row for row in rows if row["target_id"] in used[Path(source["path"]).name]]
        if {row["target_id"] for row in selected} != used[Path(source["path"]).name]:
            raise ValueError("A required retained target record is missing")
        records.update({row["target_id"]: row for row in selected})
        sources.append({"path": source["path"], "git_blob_sha": blob,
            "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw),
            "source_ref": checkpoint["blind_checkpoint_sha"],
            "target_ids_used": [row["target_id"] for row in selected],
            "fields_used": list(selected[0])})
    return {
        "prediction_contract_sha256": checkpoint["prediction_contract_sha256"],
        "blind_checkpoint_sha": checkpoint["blind_checkpoint_sha"],
        "blind_tree_sha": checkpoint["blind_tree_sha"],
        "blind_artifact_hashes_sha256": checkpoint["blind_artifact_hashes_sha256"],
        "pushed_blind_commit_time_utc": checkpoint["remote_commit_time_utc"],
        "reveal_time_utc": snapshot["reveal_time_utc"],
        "source_retrieval": "One GitHub connector fetch per source at the verified blind commit; no historical search or digitisation",
        "source_evaluation_count": 1, "sources": sources, "used_target_records": records,
        "direct_measurement": "Reported 35 percent total deficit with 4.7 percentage point SEM; Figure 1A and Results text, P15",
        "derived_measurement": "The 0.65 null to WT ratio and 0.047 SEM are algebraic complements, not independent observations",
        "replicates": "Six glands per genotype for the AE4 and AE2 secretion experiments; separate genetic backgrounds",
        "time_resolution": "One minute collection samples in the experiment; retained exact point values are absent. Early, sustained and late records are qualitative author descriptions",
        "uncertainty_policy": "Preserve SEM as SEM, without inventing a confidence interval or an equivalence margin",
        "numeric_time_points_used": [], "digitised_values_used": False,
        "AE4_ionic_context": "Resting isolated SMG acini only. No exact null root is available; no model score. No matched 600 s stimulated Cl or pH endpoint is retained",
        "AE2_context": "Separate acinus specific Ae2 deletion and matched Ae2 line controls; no detected secretion difference is not exact equality",
        "model_or_blind_trajectory_change_after_reveal": False,
        "outcome_fixed_before_reveal": CLASSIFICATION,
    }


def make_comparisons(manifest):
    records = manifest["used_target_records"]
    primary = records["H15_AE4_NULL_TOTAL_10MIN_DEFICIT"]
    target = {"D_total": float(primary["estimate"]) / 100,
        # The source gives a mean and SEM, not a confidence or equivalence interval.
        "D_total_interval": None, "relative_time_points": [],
        "evidence_status": primary["numeric_status"]}
    ae4 = compare.evaluate_rows(compare.read_csv(OUT / "ae4_null_blind_predictions.csv"), target)
    ae2 = compare.evaluate_rows(compare.read_csv(OUT / "ae2_null_blind_predictions.csv"),
        {"D_total": None, "D_total_interval": None, "relative_time_points": [],
         "evidence_status": "qualitative_reported_no_equivalence_margin"})
    provenance = {"prediction_contract_sha256": manifest["prediction_contract_sha256"],
        "blind_checkpoint_sha": manifest["blind_checkpoint_sha"],
        "holdout_source_sha256": manifest["sources"][0]["sha256"]}
    for row in ae4:
        row.update(provenance)
        row.update(target_R_total=0.65, target_D_total_SEM=0.047, target_R_total_SEM=0.047,
            target_uncertainty_type="SEM", target_n_WT=6, target_n_null=6,
            target_early_description=records["H15_AE4_NULL_EARLY_TIMECOURSE"]["estimate"],
            target_sustained_description=records["H15_AE4_NULL_SUSTAINED_TIMECOURSE"]["estimate"],
            target_late_description=records["H15_AE4_NULL_LATE_POINTWISE"]["estimate"],
            early_agreement="NOT_EVALUABLE_NO_EXACT_NULL_ROOT",
            sustained_agreement="NOT_EVALUABLE_NO_EXACT_NULL_ROOT",
            emergence_direction_agreement="NOT_EVALUABLE_NO_EXACT_NULL_ROOT",
            numeric_time_RMSE_status="No licensed prediction and no retained numerical target time points",
            experimental_confidence_interval_available=False,
            reported_mean_plus_minus_one_SEM_low=0.303,
            reported_mean_plus_minus_one_SEM_high=0.397,
            SEM_band_is_confidence_or_equivalence_interval=False,
            resting_Cl_null_observed_mM=36.5, resting_Cl_null_SEM_mM=1.6,
            resting_pH_null_observed=6.89, resting_pH_null_SEM=0.02,
            resting_ionic_comparison="NOT_EVALUABLE_NO_EXACT_NULL_ROOT",
            stimulated_ionic_comparison="No valid prediction or matched retained 600 s concentration and pH endpoint")
    for row in ae2:
        row.update(provenance)
        row.update(target_description=records["V15_AE2_SECRETION"]["estimate"],
            qualitative_secretion_assessment="Small model effect is qualitatively consistent; no statistical equivalence claim",
            target_n_control=6, target_n_null=6, exact_total_ratio_available=False,
            numeric_equivalence_margin_available=False,
            observed_AE2_line_rest_Cl_control_mM=53.4, observed_AE2_line_rest_Cl_null_mM=54.5,
            observed_AE2_line_rest_Cl_control_SEM_mM=1.8, observed_AE2_line_rest_Cl_null_SEM_mM=1.8,
            observed_AE2_line_rest_pH_control=6.87, observed_AE2_line_rest_pH_null=6.95,
            observed_AE2_line_rest_pH_control_SEM=0.01, observed_AE2_line_rest_pH_null_SEM=0.05,
            ionic_context="Separate AE2 genetic background; not pooled with the frozen AE4 line WT roots and not used for a combined quantitative score",
            AE4_minus_AE2_specificity="Unavailable because AE4 exact null predictions are absent")
    return ae4, ae2


def panels_and_classification(manifest, ae4, ae2):
    c = json.loads((OUT / "prediction_contract.json").read_text())
    expected = {(r, float(ca)) for r in c["root_ids"] for ca in c["calcium_uM"]}
    for rows in (ae4, ae2):
        if len(rows) != 30 or {(r["root_id"], float(r["calcium_uM"])) for r in rows} != expected:
            raise ValueError("The final comparison must retain every declared row")
    if any(compare.number(r.get("R_total")) is not None for r in ae4):
        raise ValueError("This completed blind checkpoint contains no valid AE4 predictions")
    totals = [compare.number(r["R_total"]) for r in ae2]
    combined = compare.read_csv(OUT / "blind_pair_summary.csv")
    return {
        "prediction_contract_sha256": manifest["prediction_contract_sha256"],
        "blind_checkpoint_sha": manifest["blind_checkpoint_sha"],
        "classification": CLASSIFICATION, "outcome_number": 6,
        "decision_was_fixed_before_reveal": True, "valid_AE4_pairs": 0, "valid_AE2_pairs": 30,
        "root_ids": c["root_ids"], "calcium_uM": c["calcium_uM"], "all_30_panel_rows": combined,
        "AE4_panel_summaries": compare.comparison_summary(ae4),
        "AE2_panel_summaries": compare.comparison_summary(ae2),
        "AE4_range_and_median": {"R_total_min": None, "R_total_max": None, "R_total_median": None,
            "reason": "Every exact AE4 resting branch failed under the declared rules"},
        "AE2_R_total_range": [min(totals), max(totals)],
        "AE2_R_total_median": sum(sorted(totals)[14:16]) / 2,
        "AE2_specificity": "Small AE2 secretion effect survives; AE4 versus AE2 specificity cannot be established",
        "observed_AE4_D_total": 0.35, "observed_AE4_D_total_SEM": 0.047,
        "AE4_magnitude_and_timing_agreement": "Not evaluable from this blind checkpoint",
        "absolute_scale": "Inherited one SMG mapping remains unresolved and was nonblocking; it does not explain or remove resting continuation failure",
        "uncertainty": "The declared procedure failed at exact AE4 deletion. This is not proof that no equilibrium exists, and it is not a comparison of a valid AE4 curve with experiment",
        "model_or_trajectory_changes_after_reveal": False,
    }


def documents(manifest, classification):
    lo, hi = classification["AE2_R_total_range"]
    low_increase, high_increase = 100 * (lo - 1), 100 * (hi - 1)
    checkpoint = manifest["blind_checkpoint_sha"]
    source = f"https://github.com/esig626/ae4-salivary-transport-control/blob/{checkpoint}/results/13B_modern_full_model/heldout_targets.csv"
    validation = f"https://github.com/esig626/ae4-salivary-transport-control/blob/{checkpoint}/results/13B_modern_full_model/validation_targets.csv"
    footer = f"\n\nBlind checkpoint: `{checkpoint}`. Prediction contract SHA256: `{manifest['prediction_contract_sha256']}`.\n"
    reveal = f"""# Task 14 holdout reveal

The complete blind checkpoint `{checkpoint}` was committed and pushed at {manifest['pushed_blind_commit_time_utc']}. The remote branch reference was then verified and recorded before the first source fetch at {manifest['reveal_time_utc']}. Its tree and every blind artifact hash remain unchanged. The original blind prediction already classified this task as a continuation failure.

The retained [strict holdout ledger]({source}) records a **35 ± 4.7% reduction in total 600 s AE4 null secretion**, with 4.7 percentage points reported as SEM and six glands per genotype. The derived null to WT ratio is 0.65 ± 0.047. This derived row is not independent evidence. The stated SEM is preserved; the mean ± one SEM band is 0.303 to 0.397 in reduction fraction, but it is not labelled a confidence interval, an equivalence margin or a biological acceptance law.

The time course records comparable flow over the first 2 to 3 minutes, a sustained decrease after 3 minutes, and flow below half of WT in later minutes. These are qualitative author descriptions. The retained source supplies no exact minute values, no exact minute for the late threshold and no numerical minute uncertainty. No curve was digitised, and no numeric early ratio, sustained ratio, onset time or RMSE target was invented.

The [validation ledger]({validation}) records no detected AE2 genotype difference in secretion kinetics or total over 10 minutes, with six glands per group and plotted SEM. It supplies no numerical equivalence margin. AE2 line controls are a separate genetic background and were not pooled with the AE4 line.

Independent resting SMG acinar measurements in the same ledger give AE4 null intracellular Cl of 36.5 ± 1.6 mM (SEM, n = 6) and pH of 6.89 ± 0.02 (SEM, n = 4). They are resting measurements, not stimulated endpoints. No valid exact AE4 resting state exists in the saved predictions, so neither quantity is scored. Initial SPQ uptake slopes and isolated exchanger assays are not substituted for 600 s intracellular Cl or pH measurements.

The AE2 line resting measurements, retained separately for context, are Cl 53.4 ± 1.8 mM in controls and 54.5 ± 1.8 mM after deletion, and pH 6.87 ± 0.01 and 6.95 ± 0.05, respectively. All uncertainties are SEM. These distinct genetic controls do not supply a joint quantitative score for the frozen AE4 line root family.

All source paths, complete file SHA256 values, Git blob hashes, source anchors, exact used fields, replicate counts, uncertainty labels and direct versus derived status are preserved in `reveal_manifest.json`. The source rows were evaluated once. Deterministic replay uses those retained records without reopening the sources.
""" + footer
    final = f"""# Task 14 final answer

**The frozen model did not produce a valid AE4 holdout prediction.** All ten AE4 resting continuation branches reached 5% expression but failed at exact deletion. The independent implementation reproduced every failure. All 60 exact zero attempts across the two implementations exhausted the declared 3000 evaluation budget; their scaled residuals were 0.0184222 to 0.0232276 against the required 1e-7. The failed iterates are not equilibria and were not used as initial states.

Every AE2 resting branch succeeded. The blind stage completed 30 WT and 30 AE2 production Radau trajectories over all ten roots and calcium 0.10, 0.25 and 0.50 uM. All 12 available representative BDF confirmations passed, with maximum relative discrepancy 4.65e-6 below the inherited 1e-4 tolerance. The largest dynamic accounting ratio was 0.00657 of its allowed tolerance; bulk charge error was at most 8.11e-13 fmol. Independent resting states agreed within 3.46e-11 in relative state difference wherever both implementations retained a root. All 105 blind tests passed.

The blind model predicts an AE2 null to WT total secretion ratio from **{lo:.8f} to {hi:.8f}**, a secretion increase of **{low_increase:.4f}% to {high_increase:.4f}%**. Every root and calcium remains in the CSV files and routing family summaries. This is qualitatively consistent with the reported absence of a detected AE2 secretion difference. There is no reported equivalence margin, so it is not a statistical proof of equivalence.

After the checkpoint was pushed, the experiment revealed **35 ± 4.7% lower AE4 null total secretion (SEM)**, with comparable early flow followed by a sustained deficit. The model has no valid AE4 total ratio, minute ratio, early or sustained integrated ratio, emergence statistic, or exact null ionic prediction to compare with those data. Signed and absolute errors, interval membership and time course RMSE therefore remain explicitly unavailable in all 30 AE4 comparison rows. No root or calcium value was dropped or selected to improve agreement.

The small AE2 effect survives this task. **AE4 versus AE2 specificity is not established**, because the AE4 effect and the effect contrast are unavailable. The outcome is neither an AE4 magnitude match nor a demonstrated mismatch of a valid predicted curve. It is a numerical continuation failure under the frozen rules. This does not prove that an exact null equilibrium cannot exist.

The inherited absolute one SMG mapping discrepancy remains unresolved. The Task 13C best 0.50 uM result was 5.405993 uL/min against the inherited 9 to 10 uL/min mapping, and five roots failed the minute observation shape criterion. Those results are preserved. A common multiplicative scale would cancel from valid genotype ratios, but it cannot supply the missing exact AE4 resting states.

All model equations, parameters, original roots, regulation, calcium values and saved blind trajectories remained unchanged after reveal. Any investigation of why the exact deletion continuation fails belongs in a later task. This task makes no new calibration, mechanism or biological validation claim.

{footer}
**{CLASSIFICATION}**
"""
    return {"holdout_reveal.md": reveal, "final_answer.md": final}


def evaluate(manifest):
    compare.verify_blind_ledger(REPO, OUT / "blind_artifact_hashes.csv")
    ae4, ae2 = make_comparisons(manifest)
    classification = panels_and_classification(manifest, ae4, ae2)
    files = {"reveal_manifest.json": json_bytes(manifest),
        "ae4_holdout_comparison.csv": compare.render_csv(ae4),
        "ae2_validation_comparison.csv": compare.render_csv(ae2),
        "final_classification.json": json_bytes(classification)}
    return files, documents(manifest, classification)


def write_evaluation(manifest):
    files, docs = evaluate(manifest)
    for name, data in files.items(): compare.write_final(OUT, name, data)
    for name, text in docs.items():
        if name not in ("holdout_reveal.md", "final_answer.md"):
            raise ValueError("Only final analysis documents may be written")
        path = ANALYSIS / name
        if path.is_symlink(): raise ValueError("Comparison document must not be a symlink")
        path.write_text(text)
    compare.verify_blind_ledger(REPO, OUT / "blind_artifact_hashes.csv")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-snapshot", type=Path)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.replay:
        manifest = json.loads((OUT / "reveal_manifest.json").read_text())
    else:
        checkpoint = json.loads((OUT / "blind_checkpoint.json").read_text())
        if args.source_snapshot is None: parser.error("A single retrieved source snapshot is required")
        manifest = source_records(json.loads(args.source_snapshot.read_text()), checkpoint)
    write_evaluation(manifest)
    print(CLASSIFICATION)


if __name__ == "__main__": main()
