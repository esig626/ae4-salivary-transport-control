"""Independently apply every inherited gate to nearby sensitivity roots.

The sensitivity calculation retains a compact concentration/positivity screen.
This supplement re-solves its parameter and expression perturbations and calls
the unchanged Task 40 diagnostic, including its inherited Task 37 physiological,
conservation and NBC source checks plus both equal-routing architecture checks.
No gate tolerance is relaxed.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

from core import reference, clone, Scaling, stationary, Genotype, write_json
from sensitivity_analysis import parameter_builders, STEPS


HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
EXPRESSION_STEPS = (1e-4, 5e-5)


def main():
    _, stimulus, initial, task, _ = reference()
    model = clone(stimulus)
    scale = Scaling(model, initial)
    root_path = OUT / "equilibria.json"
    sensitivity_path = OUT / "sensitivity_extended.json"
    saved = json.loads(root_path.read_text())
    sensitivity = json.loads(sensitivity_path.read_text())
    seeds = {
        row["ae4_expression"]: np.asarray(row["state"])
        for row in saved["roots"] if row["mechanism"] == "C0"
    }
    builders = parameter_builders(stimulus, model)
    assert set(builders) == {row["parameter"] for row in sensitivity["parameters"]}
    assert list(STEPS) == sensitivity["independent_log_parameter_steps"]
    cases = []

    def check_root(candidate, seed, expression, metadata):
        genotype = Genotype("independent_sensitivity_gate_check", ae4_expression=expression)
        case = dict(**metadata, expression=expression)
        try:
            root = stationary(candidate, scale, seed, genotype)
            state = np.asarray(root["state"])
            evaluation = candidate.evaluate(1., state, genotype=genotype)
            row, failures, ratios = task.diagnose(candidate, 1., state, evaluation)
            case.update(
                numerical_status="root_obtained",
                root_solver_success=root["root_success"],
                root_solver_message=root["root_message"],
                independent_rhs_max=root["residual_max"],
                full_rhs_max=float(np.max(np.abs(evaluation.rhs))),
                condition_number=root["condition_number"],
                full_inherited_gate_failures=failures,
                conservation_tolerance_ratios=ratios,
                maximum_conservation_tolerance_ratio=max(ratios.values()),
                physiological=not failures,
                pH_i=row["ph_i"],
                Na_i_mM=row["na_i_mM"],
                status="passed" if not failures else "failed_inherited_gates",
            )
            cases.append(case)
            return state
        except Exception as exc:
            case.update(status="failed_numerical_check", numerical_status="unsuccessful",
                        error_type=type(exc).__name__, error=str(exc))
            cases.append(case)
            return None

    for parameter, selection in builders.items():
        for h in STEPS:
            for sign in (1, -1):
                log_step = sign*h
                candidate = selection["build"](log_step)
                common = dict(parameter=parameter, log_parameter_step=log_step)
                wt = check_root(candidate, seeds[1.], 1.,
                                dict(**common, group="parameter_WT_null"))
                check_root(candidate, seeds[0.], 0.,
                           dict(**common, group="parameter_WT_null"))
                # These are exactly the additional roots used to validate the
                # WT expression derivative at every perturbed parameter point.
                for expression_step in EXPRESSION_STEPS:
                    for expression_sign in (1, -1):
                        check_root(candidate, wt if wt is not None else seeds[1.],
                                   1.+expression_sign*expression_step,
                                   dict(**common, group="compensation_expression",
                                        expression_offset=expression_sign*expression_step))

    groups = {}
    for name, expected in (("parameter_WT_null", 128), ("compensation_expression", 256)):
        members = [case for case in cases if case["group"] == name]
        groups[name] = dict(expected=expected, attempted=len(members),
                            passed=sum(case["status"] == "passed" for case in members),
                            failed=sum(case["status"] != "passed" for case in members))
    failed = [case for case in cases if case["status"] != "passed"]
    complete = all(group["attempted"] == group["expected"] for group in groups.values())
    result = dict(
        status="passed" if complete and not failed else "failed",
        scope="Fresh root re-solves with all Task 40 physiological, finite-value, conservation, NBC source and equal-routing architecture gates; no tolerance changes.",
        gate_source="analysis/40_ae4_equal_cation_routing/validation_common.py::diagnose -> analysis/37_wt_nbc_validation/validation_common.py::diagnose",
        parameter_count=len(builders), log_parameter_steps=list(STEPS),
        expression_validation_steps=list(EXPRESSION_STEPS), groups=groups,
        maximum_conservation_tolerance_ratio=max(
            (case.get("maximum_conservation_tolerance_ratio", 0.) for case in cases), default=0.),
        maximum_independent_rhs=max(
            (case.get("independent_rhs_max", 0.) for case in cases), default=0.),
        source_sha256={
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__), HERE / "sensitivity_analysis.py", root_path, sensitivity_path)
        },
        limitations=[
            "These are local equilibrium checks at numerical perturbations, not uncertainty intervals or finite-time trajectory gates.",
            "Gate acceptance does not remove the inherited thermodynamic opposition of the C0 source convention.",
            "A solver status flag is retained separately; root acceptance follows the inherited stationary residual check and full diagnostic gates.",
        ],
        cases=cases,
    )
    write_json(OUT / "sensitivity_full_gate_verification.json", result)
    print("SENSITIVITY_FULL_GATES", result["status"],
          "parameter_roots", groups["parameter_WT_null"]["passed"],
          "compensation_roots", groups["compensation_expression"]["passed"],
          "failed", len(failed), flush=True)
    if result["status"] != "passed":
        raise AssertionError("Full inherited sensitivity gate verification failed; cases retained in JSON.")


if __name__ == "__main__":
    main()
