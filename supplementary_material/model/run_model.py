"""Run the Ae4 salivary secretion model used in the manuscript.

Examples
--------
Run the principal case::

    python run_model.py

Run all three auxiliary-conductance cases::

    python run_model.py --case all --output ae4_output
"""
from pathlib import Path
import argparse
import json
import sys

MODEL_DIR = Path(__file__).resolve().parent
SUPPLEMENTARY_DIR = MODEL_DIR.parent
sys.path.insert(0, str(MODEL_DIR))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=("case_01", "case_02", "case_03", "all"),
        default="case_02",
        help="Simulation case to run. case_02 is the principal manuscript case.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SUPPLEMENTARY_DIR.parent / "ae4_output",
        help="New directory for simulation output. Default: <repository>/ae4_output",
    )
    args = parser.parse_args()

    output = args.output.resolve()
    protected = (
        SUPPLEMENTARY_DIR / "data",
        SUPPLEMENTARY_DIR / "model",
        SUPPLEMENTARY_DIR / "parameters",
        SUPPLEMENTARY_DIR / "provenance",
    )
    if any(output == path or path in output.parents for path in protected):
        parser.error("Choose an output directory outside the supplied model and data directories.")
    if output.exists():
        parser.error("The output directory already exists. Choose a new output directory.")

    # Import after argument parsing so the numerical thread settings in the
    # archived runner are applied before NumPy/SciPy are initialised.
    import run_frozen_cases as runner
    from run_reporter_recovery_52D import encode_observations, patched_attempt_source
    from threadpoolctl import threadpool_limits

    freeze = json.loads(
        (SUPPLEMENTARY_DIR / "parameters/parameter_and_case_freeze.json").read_text()
    )
    cases = [
        case
        for case in freeze["case_matrix"]["cases"]
        if case["id"] in ("case_01", "case_02", "case_03")
        and (args.case == "all" or case["id"] == args.case)
    ]

    onsets = freeze["projected_onsets"]["central"]
    y0 = runner.np.concatenate([onsets[g]["state"] for g in ("WT", "KO")])
    settings = freeze["case_matrix"]["solver"]

    # Reuse the archived integration and reporting routine while bypassing only
    # the historical one-time branch/output bookkeeping.
    _, source = patched_attempt_source()
    namespace = dict(runner.__dict__, encode_observations=encode_observations)
    exec(compile(source, "manuscript_run_attempt", "exec"), namespace)
    run_attempt = namespace["run_attempt"]

    output.mkdir(parents=True)
    receipts = []

    with threadpool_limits(limits=1):
        for case in cases:
            directory = output / case["id"]
            directory.mkdir()

            model = runner.PairedModel(
                protocol=case["protocol"],
                G_aux_S=case["G_aux_S"],
                ko_supply_multiplier=case["ko_supply_multiplier"],
            )
            if not runner.dependency_audit()["pass"]:
                raise RuntimeError(
                    "A model dependency was loaded from outside the supplied model snapshot."
                )

            result = run_attempt(case, model, y0.copy(), "Radau", directory, settings)
            method = "Radau"
            if result["status"] == "solver_status_failure":
                result = run_attempt(case, model, y0.copy(), "BDF", directory, settings)
                method = "BDF"

            headline = result["headline_0_600"]
            receipt = {
                "case": case["id"],
                "method": method,
                "status": result["status"],
                "headline_0_600": headline,
            }
            receipts.append(receipt)

            if headline is not None:
                wt = headline["fluid_pL"]["WT"]
                ko = headline["fluid_pL"]["KO"]
                deficit = 100.0 * (1.0 - ko / wt)
                print(
                    f"{case['id']}: {result['status']} | "
                    f"WT={wt:.6f} pL | KO={ko:.6f} pL | "
                    f"deficit={deficit:.3f}%"
                )
            else:
                print(f"{case['id']}: {result['status']}")

    (output / "run_summary.json").write_text(
        json.dumps(receipts, indent=2, allow_nan=False) + "\n"
    )

    if any(row["status"] != "completed" for row in receipts):
        raise SystemExit("At least one requested case did not complete its numerical checks.")


if __name__ == "__main__":
    main()
