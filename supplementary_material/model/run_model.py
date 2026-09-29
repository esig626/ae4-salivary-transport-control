"""Run the Ae4 salivary secretion simulations used in the manuscript.

The default command runs the principal combined CCh + IPR control/knockout
comparison. Other article simulations are selected with --simulation.
"""
from pathlib import Path
import argparse
import csv
import json
import sys

MODEL_DIR = Path(__file__).resolve().parent
SUPPLEMENTARY_DIR = MODEL_DIR.parent
REPOSITORY_DIR = SUPPLEMENTARY_DIR.parent
sys.path.insert(0, str(MODEL_DIR))

G_AUX_CENTRAL_S = 2.32e-9
G_AUX_LOW_S = 0.0
G_AUX_HIGH_S = 4.49e-9
SIGMA_ROOT = 1.4913598957


def simulation_cases(name, sigma):
    if name == "principal":
        return [dict(id="principal", protocol="CCH_IPR",
                     G_aux_S=G_AUX_CENTRAL_S, ko_supply_multiplier=1.0)]
    if name == "auxiliary":
        return [
            dict(id="aux_0nS", protocol="CCH_IPR",
                 G_aux_S=G_AUX_LOW_S, ko_supply_multiplier=1.0),
            dict(id="aux_2p32nS", protocol="CCH_IPR",
                 G_aux_S=G_AUX_CENTRAL_S, ko_supply_multiplier=1.0),
            dict(id="aux_4p49nS", protocol="CCH_IPR",
                 G_aux_S=G_AUX_HIGH_S, ko_supply_multiplier=1.0),
        ]
    if name == "cch-only":
        return [dict(id="cch_only", protocol="CCH_ONLY",
                     G_aux_S=G_AUX_CENTRAL_S, ko_supply_multiplier=1.0)]
    if name == "ipr-only":
        return [dict(id="ipr_only", protocol="IPR_ONLY",
                     G_aux_S=G_AUX_CENTRAL_S, ko_supply_multiplier=1.0)]
    if name == "compensation":
        return [dict(id=f"sigma_{sigma:.6f}".replace(".", "p"),
                     protocol="CCH_IPR", G_aux_S=G_AUX_CENTRAL_S,
                     ko_supply_multiplier=float(sigma))]
    if name == "compensation-root":
        return [dict(id="sigma_root", protocol="CCH_IPR",
                     G_aux_S=G_AUX_CENTRAL_S,
                     ko_supply_multiplier=SIGMA_ROOT)]
    raise ValueError(name)


def run_case(case, output, runner, run_attempt, y0, settings):
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
    row = {
        "case": case["id"],
        "protocol": case["protocol"],
        "G_aux_nS": case["G_aux_S"] * 1e9,
        "sigma": case["ko_supply_multiplier"],
        "method": method,
        "status": result["status"],
        "recorded_end_s": result["recorded_end_s"],
    }
    if headline is not None:
        wt = headline["fluid_pL"]["WT"]
        ko = headline["fluid_pL"]["KO"]
        row.update(
            wt_fluid_pL=wt,
            ko_fluid_pL=ko,
            deficit_percent=100.0 * (1.0 - ko / wt),
        )
        print(
            f"{case['id']}: {result['status']} | "
            f"WT={wt:.6f} pL | KO={ko:.6f} pL | "
            f"deficit={row['deficit_percent']:.3f}%"
        )
    else:
        print(
            f"{case['id']}: {result['status']} | "
            f"recorded to {result['recorded_end_s']:.3f} s"
        )
    return row, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--simulation",
        choices=(
            "principal",
            "auxiliary",
            "cch-only",
            "ipr-only",
            "compensation",
            "compensation-root",
            "compensation-scan",
        ),
        default="principal",
        help="Article simulation to run. Default: principal.",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=1.0,
        help="Knockout NKCC1/Ae2 multiplier for --simulation compensation.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY_DIR / "ae4_output",
        help="New output directory. Default: <repository>/ae4_output",
    )
    args = parser.parse_args()

    if args.simulation == "compensation" and not (1.0 <= args.sigma <= 2.01):
        parser.error("--sigma must lie between 1.0 and 2.01.")

    output = args.output.resolve()
    protected = (
        SUPPLEMENTARY_DIR / "data",
        SUPPLEMENTARY_DIR / "model",
        SUPPLEMENTARY_DIR / "parameters",
        SUPPLEMENTARY_DIR / "provenance",
    )
    if any(output == path or path in output.parents for path in protected):
        parser.error(
            "Choose an output directory outside the supplied model and data directories."
        )
    if output.exists():
        parser.error("The output directory already exists. Choose a new output directory.")

    # The archived runner sets numerical thread limits before importing NumPy/SciPy.
    import run_frozen_cases as runner
    from run_reporter_recovery_52D import encode_observations, patched_attempt_source
    from threadpoolctl import threadpool_limits

    freeze = json.loads(
        (SUPPLEMENTARY_DIR / "parameters/parameter_and_case_freeze.json").read_text()
    )
    onsets = freeze["projected_onsets"]["central"]
    y0 = runner.np.concatenate([onsets[g]["state"] for g in ("WT", "KO")])
    settings = freeze["case_matrix"]["solver"]

    # Reuse the numerical integration and reporting routine while bypassing only
    # the historical branch/output bookkeeping attached to its original run.
    _, source = patched_attempt_source()
    namespace = dict(runner.__dict__, encode_observations=encode_observations)
    exec(compile(source, "manuscript_run_attempt", "exec"), namespace)
    run_attempt = namespace["run_attempt"]

    output.mkdir(parents=True)
    rows = []

    with threadpool_limits(limits=1):
        if args.simulation == "compensation-scan":
            sigmas = [1.0 + 0.025 * i for i in range(41)]
            for sigma in sigmas:
                case = dict(
                    id=f"sigma_{sigma:.3f}".replace(".", "p"),
                    protocol="CCH_IPR",
                    G_aux_S=G_AUX_CENTRAL_S,
                    ko_supply_multiplier=sigma,
                )
                row, _ = run_case(case, output, runner, run_attempt, y0, settings)
                rows.append(row)
        else:
            cases = simulation_cases(args.simulation, args.sigma)
            for case in cases:
                row, _ = run_case(case, output, runner, run_attempt, y0, settings)
                rows.append(row)

    (output / "run_summary.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False) + "\n"
    )

    if args.simulation == "compensation-scan":
        columns = [
            "sigma", "status", "recorded_end_s",
            "wt_fluid_pL", "ko_fluid_pL", "deficit_percent",
        ]
        with (output / "compensation_scan.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key) for key in columns})

    # IPR-only is expected to terminate when the prescribed pH acceptance limit
    # is reached, as reported in the manuscript. Other incomplete runs are errors.
    if args.simulation != "ipr-only":
        incomplete = [row for row in rows if row["status"] != "completed"]
        if incomplete:
            raise SystemExit("At least one requested simulation did not complete its checks.")


if __name__ == "__main__":
    main()
