"""Task 50A provenance audit. No model imports, RHS calls or numerical solves.

Compare exact decimal literals, verify frozen input hashes, and search the
recorded remote heads. This checks provenance, not model predictions.
"""
from pathlib import Path
from decimal import Decimal
import ast
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
START = "8f5fefd9fc0b563a6bd963895211b4acf76972b0"
CLASSIFICATION = "TARGET-CALIBRATED CONSTRUCTION"
PATTERN = (
    r"beta.?condition|ae4.{0,30}(cacc|tmem16|apical|secretory)|"
    r"recruitment_factor|0\.10511872843288|0\.89488127156712"
)


def git(*args):
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True)
    return result.stdout


def blob(path, ref=START):
    return git("show", f"{ref}:{path}")


def record(path, ref=START):
    return json.loads(blob(path, ref), parse_float=str)


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    OUT.mkdir(exist_ok=True)
    remote_file = OUT / "remote_heads.json"
    if remote_file.exists():
        heads = json.loads(remote_file.read_text())
    else:
        lines = git("for-each-ref", "--format=%(objectname) %(refname)",
                    "refs/remotes/origin/").decode().splitlines()
        heads = [{"sha": row.split(" ", 1)[0], "ref": row.split(" ", 1)[1]}
                 for row in lines if not row.endswith("/HEAD")]
        write("remote_heads.json", heads)

    matches = []
    source_blobs = {}
    for sha in sorted({head["sha"] for head in heads}):
        result = subprocess.run(
            ["git", "grep", "-I", "-l", "-i", "-E", PATTERN, sha, "--",
             "src", "prompts", "analysis", "results"],
            cwd=ROOT, capture_output=True)
        if result.returncode not in (0, 1):
            raise RuntimeError(result.stderr.decode())
        paths = [line.split(":", 1)[1] for line in result.stdout.decode().splitlines()]
        matches.append({"sha": sha, "matching_path_count": len(paths),
                        "matching_path_list_sha256": hashlib.sha256(
                            ("\n".join(paths) + "\n").encode()).hexdigest()})
        for path in paths:
            if path.startswith("src/"):
                key = git("rev-parse", f"{sha}:{path}").decode().strip()
                item = source_blobs.setdefault(key, {"path": path, "refs": []})
                item["refs"].append(sha)
    write("novelty_search.json", {
        "baseline": START, "pattern": PATTERN,
        "scope": ["src", "prompts", "analysis", "results"],
        "heads_file": "remote_heads.json", "head_count": len(heads),
        "unique_head_count": len(matches), "matches": matches,
        "matching_source_blobs": source_blobs,
        "limit": "Text search and source review, not a proof against arbitrary renamed implementations",
    })

    prefix = "results/41_ae4_loss_algebraic_design/candidate_08/"
    design = record(prefix + "design.json")
    frozen = record(prefix + "frozen_inputs.json")
    selected_path = "src/modern_full_model/task41_selected.py"
    text = blob(selected_path).decode()
    tree = ast.parse(text)
    selected = next(ast.get_source_segment(text, node.value) for node in tree.body
                    if isinstance(node, ast.Assign)
                    and any(isinstance(target, ast.Name)
                            and target.id == "SELECTED_NULL_CACC_RECRUITMENT"
                            for target in node.targets))
    b = design["null_stimulated_cacc_fraction"]
    assert b == frozen["selected_parameters"]["null_stimulated_fraction"] == selected

    hashes = []
    for item in frozen["execution_files"]:
        data = blob(item["path"])
        digest = hashlib.sha256(data).hexdigest()
        assert digest == item["sha256"], item["path"]
        hashes.append({**item, "matches_frozen_sha256": True})

    extra = [selected_path, prefix + "frozen_inputs.json",
             "src/modern_full_model/nbc_minimal.py",
             "src/modern_full_model/validation.py",
             "analysis/40_ae4_equal_cation_routing/validation_common.py"]
    for path in extra:
        hashes.append({"path": path, "sha256": hashlib.sha256(blob(path)).hexdigest()})
    write("audited_inputs.json", {"baseline": START, "files": hashes})

    lam = Decimal("0.89488127156712")
    reported_b = Decimal("0.10511872843288")
    stored_b = Decimal(b)
    assert Decimal(1) - lam == reported_b
    difference = Decimal(1) - lam - stored_b
    assert difference != 0
    # Exact arithmetic on provenance literals only; no parameter inference.
    result = {
        "classification": CLASSIFICATION,
        "audit_result": "STOP_EXACT_FROZEN_COEFFICIENT_IDENTITY_FAILS",
        "mandated_lambda_literal": str(lam),
        "reported_task41_b_literal": str(reported_b),
        "frozen_task41_b_literal": str(stored_b),
        "driver_and_frozen_parameter_literals_agree": True,
        "one_minus_mandated_lambda": str(Decimal(1) - lam),
        "one_minus_frozen_task41_b_for_comparison_only": str(Decimal(1) - stored_b),
        "task50_minus_frozen_task41_factor_at_a1_beta1_e0": str(difference),
        "symbolic_factor_difference_at_beta1": "(-4.46e-15) * a * (1-e)",
        "formal_parent_nesting": {"WT": True, "REST": True,
                                  "CCh_only_AE4_null": True, "AE2_KO": True},
        "formal_equivalence_to_reported_rounded_b": True,
        "formal_equivalence_to_actual_frozen_b": False,
        "frozen_execution_hashes_verified": len(frozen["execution_files"]),
        "new_model_evaluations": 0,
        "new_production_trajectories": 0,
        "new_parameter_fits_or_searches": 0,
        "scientific_source_modified": False,
        "precision_note": "No numerical magnitude or biological failure is inferred from this exact provenance mismatch",
    }
    write("coefficient_audit.json", result)
    print(json.dumps({"status": result["audit_result"],
                      "remote_heads_searched": len(heads),
                      "unique_source_blobs_for_review": len(source_blobs),
                      "frozen_execution_hashes_verified": len(frozen["execution_files"]),
                      "new_rhs_calls": 0, "new_trajectories": 0}))


if __name__ == "__main__":
    main()
