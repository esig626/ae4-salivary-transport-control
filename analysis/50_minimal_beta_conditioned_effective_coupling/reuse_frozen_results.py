"""Task 50C: verify and read frozen artefacts using only the standard library.

No model, NumPy, SciPy, optimiser or trajectory runner is imported.
"""
import csv
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
P40 = ROOT / "results/40_ae4_equal_cation_routing"
P41 = ROOT / "results/41_ae4_loss_algebraic_design"
FINAL = P41 / "candidate_08"
FROZEN_COMMIT = "86f135e075aa91a70a752b650d7bb894353fee04"
CHECKPOINT_B = "537ce3204d11a1a6cc51b1dd39032d12cb42df5e"
CLASSIFICATION = "TARGET-CALIBRATED CONSTRUCTION"


def read(path):
    return json.loads(path.read_text())


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def main():
    if (OUT / "inherited_phenotype.json").exists():
        raise SystemExit("Frozen reuse extraction already exists; no optional repeat")
    verified = read(OUT / "verification_attempt_02.json")
    assert verified["status"] == "PASS" and verified["tests_run"] == 10
    for path, digest in verified["source_hashes"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", CHECKPOINT_B, "HEAD"], cwd=ROOT)
    assert ancestor.returncode == 0

    paths = [*FINAL.iterdir(), P41 / "final_summary.json", P41 / "completion.json",
             P41 / "final_cumulative_deficit.csv", P41 / "final_conservation_summary.csv",
             P40 / "wt_rest.json", P40 / "frozen_inputs.json", P40 / "wt_verification.json",
             P40 / "wt_timeseries.csv", P40 / "ae4_5pct_verification.json",
             P40 / "ae4_null_verification.json",
             ROOT / "src/modern_full_model/ae4_cacc_recruitment.py",
             ROOT / "src/modern_full_model/task41_selected.py"]
    hashes = []
    for path in sorted(paths):
        relative = str(path.relative_to(ROOT))
        data = path.read_bytes()
        frozen = subprocess.check_output(["git", "show", f"{FROZEN_COMMIT}:{relative}"], cwd=ROOT)
        assert data == frozen, relative
        git_blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        hashes.append({"path": relative, "sha256": hashlib.sha256(data).hexdigest(),
                       "git_blob_sha1": git_blob, "bytes": len(data),
                       "byte_identical_to_frozen_publication": True})

    f40, f41 = read(P40 / "frozen_inputs.json"), read(FINAL / "frozen_inputs.json")
    rest = read(P40 / "wt_rest.json")
    rest41 = read(FINAL / "wt_rest.json")
    assert rest["state_vector"] == rest41["state_vector"]
    state_hash = canonical_hash(rest["state_vector"])
    assert state_hash == rest["state_sha256"] == f40["frozen_rest_state_sha256"] == f41["frozen_rest_state_sha256"]
    for mode in ("rest_parameter_hashes", "stimulus_parameter_hashes"):
        shared = {k:v for k,v in f41[mode].items() if k != "ae4_cacc_recruitment"}
        assert shared == f40[mode]
    cases = {c: read(FINAL / (c + "_verification.json")) for c in ("wt", "ae4_5pct", "ae4_null")}
    wt40 = read(P40 / "wt_verification.json")
    for case, record in cases.items():
        assert record["status"] == "PASS" and record["complete_600_s"]
        assert record["initial_state_vector"] == rest["state_vector"]
        assert record["initial_state_sha256"] == state_hash
        assert record["parameter_hashes"] == f41["stimulus_parameter_hashes"]
        assert record["stimulus"] == wt40["stimulus"]
        assert record["onset_convention"] == wt40["onset_convention"]
        assert record["state_names"] == wt40["state_names"]
        assert all(v <= 1 for v in record["max_conservation_ratios"].values())
        assert len(record["integration_attempts"]) == 1
        assert record["integration_attempts"][0]["solver"] == wt40["integration_attempts"][0]["solver"]
    assert cases["wt"]["secretion"] == wt40["secretion"]

    summary = read(P41 / "final_summary.json")
    assert summary["target_used_for_design"] is True
    assert summary["independently_validated_mechanism"] is False
    assert summary["long_run_validity_established"] is False
    with (P41 / "final_cumulative_deficit.csv").open() as handle:
        cumulative_rows = list(csv.DictReader(handle))
    window = {}
    for case in ("ae4_5pct", "ae4_null"):
        rows = [r for r in cumulative_rows if r["case"] == case and float(r["time_s"]) > 0]
        assert [float(r["time_s"]) for r in rows] == list(range(60, 601, 60))
        deficits = [float(r["deficit_percent"]) for r in rows]
        window[case] = {"times_s": [float(r["time_s"]) for r in rows],
                        "saved_cumulative_deficits_percent": deficits,
                        "minimum_percent": min(deficits), "maximum_percent": max(deficits)}
    inherited = {
        "checkpoint": "50C", "status": "QUALIFIED_FROZEN_REUSE_ACCEPTED",
        "classification": CLASSIFICATION,
        "frozen_publication_commit": FROZEN_COMMIT,
        "verified_implementation_commit": CHECKPOINT_B,
        "fixed_task50_lambda": "0.89488127156712",
        "frozen_task41_b": f41["selected_parameters"]["null_stimulated_fraction"],
        "cases": summary["cases"], "endpoint": summary["endpoint"],
        "saved_cumulative_window": window, "rest_state_sha256": state_hash,
        "protocol": wt40["stimulus"],
        "state_protocol_and_shared_parameter_records_match": True,
        "parent_nesting": {"WT": "EXACT", "REST": "EXACT", "AE4_null_CCh_only": "EXACT", "AE2_KO": "EXACT"},
        "mutant_reuse_basis": "User authorised numerical equivalence for coefficient rounding; six fixed state comparisons pass; no new integration or global trajectory error bound",
        "unresolved": ["Chronic knockout resting Cl and pH adaptation", "Positive IPR only AE4 knockout uptake from genotype rest", "Microscopic identity of the effective coupling", "Long duration and parameter robustness"],
        "new_model_evaluations": 0, "new_production_trajectories": 0,
        "new_stationary_solves": 0, "new_fits_or_searches": 0,
        "long_run_validity_established": False, "independent_validation": False,
    }
    (OUT / "frozen_reuse_manifest.json").write_text(json.dumps({"frozen_publication_commit":FROZEN_COMMIT,"files":hashes}, indent=2, sort_keys=True) + "\n")
    (OUT / "inherited_phenotype.json").write_text(json.dumps(inherited, indent=2, sort_keys=True) + "\n")
    labels = {"wt":"WT", "ae4_5pct":"AE4 5%", "ae4_null":"AE4 null"}
    table = "\n".join(f"| {labels[r['case']]} | {r['cumulative_0_600_pL']:.9f} | {r['deficit_percent']:.6f}% | {r['endpoint_deficit_percent']:.6f}% |" for r in summary["cases"])
    text = f"""# Task 50C frozen result reuse

**QUALIFIED FROZEN REUSE ACCEPTED. TARGET-CALIBRATED CONSTRUCTION.**

All {len(hashes)} selected frozen artefacts are byte identical to Task 41
publication `{FROZEN_COMMIT}`. Their Git object identities and SHA-256 values
are in `output/frozen_reuse_manifest.json`. Initial states, state names,
protocols, onset conventions, solver settings and shared parameter hashes match.
The common rest hash is `{state_hash}`. All saved production conservation gates
passed. WT's secretion summary equals Task 40's exactly.

| Frozen case | Cumulative secretion, 0 to 600 s, pL | Cumulative deficit | Flow deficit at 600 s |
| --- | ---: | ---: | ---: |
{table}

These are inherited Task 41 results, not new simulations. The full stored
coefficient differs slightly from the complement of the mandated Task 50
lambda. The user authorised numerical treatment of that difference, and 50B
passed the predeclared fixed state comparisons. Accordingly the mutant outputs
are reused with that explicit precision qualification, not claimed as bit
identical Task 50 solutions or a rigorous global trajectory error bound.

The saved cumulative observations at every 60 s from 60 through 600 s show
deficits between {window['ae4_5pct']['minimum_percent']:.4f}% and {window['ae4_5pct']['maximum_percent']:.4f}% at 5% AE4,
and between {window['ae4_null']['minimum_percent']:.4f}% and {window['ae4_null']['maximum_percent']:.4f}% at null AE4.
The endpoint flow deficits are also substantial. Thus the inherited effect is
present across the recorded window and is not supported by only one selected
instant. This is a finite window observation, not a stationary, infinite time
or parameter robustness result. No experimental timing was fitted or used as
a gate. The existing table, not a newly integrated trajectory, supplies these
observations.

WT, REST, AE4 null CCh only and AE2 knockout remain exact parent initial value
problems under the declared matching states and protocols. This does not prove
agreement with every corresponding experiment. No directly matching CCh only
or AE2 Task 40/41 production cohort was identified in this frozen set, so the
formal identity and 50B software checks suffice; no additional trajectory was
run and no different Task 48 bath or chronic resting state was substituted.

Chronic established genotype resting chloride/pH adaptation remains unresolved.
For IPR only the calcium coordinate is zero, so the added coupling is off and
cannot repair the known positive uptake discrepancy from exact null genotype
rest. The microscopic coupling and long duration physiological validity also
remain unresolved. The frozen null endpoint has retained chloride and increased
pH/volume; passing the original finite window gates is not preserved overall
homeostasis. The phenotype selected parameter is not independently identified
and supplies no evidence of a molecular AE4 to TMEM16A interaction.

New 50C model evaluations, production trajectories, stationary solves and fits:
**zero**. `reuse_frozen_results.py` imports only the standard library, reads
existing results and verifies their provenance.
"""
    (HERE / "FROZEN_RESULT_REUSE.md").write_text(text)
    print(json.dumps({"status":inherited["status"],"frozen_files_verified":len(hashes),
                      "deficits_percent":{r["case"]:r["deficit_percent"] for r in summary["cases"]},
                      "new_trajectories":0}))


if __name__ == "__main__":
    main()
