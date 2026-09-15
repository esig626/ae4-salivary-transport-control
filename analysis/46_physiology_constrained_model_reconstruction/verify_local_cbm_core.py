"""Task 46G: exact analytical checks of the local CBM core only.

No AE4 phenotype, kinetic model or external CarbonScope package is used.
The derivation and independently specified endpoint witnesses are documented in
LOCAL_CBM_VERIFICATION.md. The output must be a new file.
"""

from __future__ import annotations

import argparse
import contextlib
from dataclasses import asdict, replace
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import platform
import sys

AREA = Path(__file__).resolve().parent
CORE = AREA / "vendor" / "carbonscope_flux_core"
sys.path.insert(0, str(CORE.parent))
sys.dont_write_bytecode = True

import numpy as np
import carbonscope_flux_core as core
from carbonscope_flux_core import smoke_test, solver

TOLERANCE = 1e-8
REACTIONS = ("U_A", "U_B", "T", "D_A", "D_B")
MATRIX = np.array([[1, 0, -1, -1, 0], [0, 1, 1, 0, -1]], dtype=float)
LOWER = np.array([0, 0, -2, 0, 0], dtype=float)
UPPER = np.array([6, 4, 3, 5, 7], dtype=float)
OBJECTIVE = np.array([0, 0, 0, 1, 1], dtype=float)

# These values and witnesses are derived algebraically, not with an LP solver.
CASES = (
    {
        "name": "whole_feasible_region", "fraction": None,
        "ranges": [[0, 6], [0, 4], [-2, 3], [0, 5], [0, 7]],
        "witnesses": [
            [[0, 0, 0, 0, 0], [6, 0, 1, 5, 1]],
            [[0, 0, 0, 0, 0], [0, 4, 0, 0, 4]],
            [[0, 2, -2, 2, 0], [3, 0, 3, 0, 3]],
            [[0, 0, 0, 0, 0], [5, 0, 0, 5, 0]],
            [[0, 0, 0, 0, 0], [3, 4, 3, 0, 7]],
        ],
    },
    {
        "name": "retain_80_percent", "fraction": 0.8,
        "ranges": [[4, 6], [2, 4], [-1, 3], [1, 5], [3, 7]],
        "witnesses": [
            [[4, 4, -1, 5, 3], [6, 2, 1, 5, 3]],
            [[6, 2, 1, 5, 3], [4, 4, -1, 5, 3]],
            [[4, 4, -1, 5, 3], [4, 4, 3, 1, 7]],
            [[4, 4, 3, 1, 7], [4, 4, -1, 5, 3]],
            [[4, 4, -1, 5, 3], [4, 4, 3, 1, 7]],
        ],
    },
    {
        "name": "retain_optimum", "fraction": 1.0,
        "ranges": [[6, 6], [4, 4], [1, 3], [3, 5], [5, 7]],
        "witnesses": [
            [[6, 4, 1, 5, 5], [6, 4, 1, 5, 5]],
            [[6, 4, 1, 5, 5], [6, 4, 1, 5, 5]],
            [[6, 4, 1, 5, 5], [6, 4, 3, 3, 7]],
            [[6, 4, 3, 3, 7], [6, 4, 1, 5, 5]],
            [[6, 4, 1, 5, 5], [6, 4, 3, 3, 7]],
        ],
    },
)


def check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def same(actual, expected) -> None:
    np.testing.assert_allclose(actual, expected, rtol=0, atol=TOLERANCE)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def core_hashes() -> dict[str, str]:
    return {
        str(p.relative_to(CORE)): sha256(p)
        for p in sorted(CORE.rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts
    }


def build_control() -> core.FluxModel:
    bounds = ((0, 6), (0, 4), (-2, 3), (0, 5), (0, 7))
    columns = (
        (("A", 1),), (("B", 1),), (("A", -1), ("B", 1)),
        (("A", -1),), (("B", -1),),
    )
    return core.FluxModel(
        metabolites=(core.FluxMetabolite("A"), core.FluxMetabolite("B")),
        reactions=tuple(
            core.FluxReaction(
                rid, tuple(core.StoichiometricTerm(mid, c) for mid, c in terms),
                low, high,
            )
            for rid, terms, (low, high) in zip(REACTIONS, columns, bounds)
        ),
        objective=core.LinearObjective(
            "maximise", (core.ObjectiveTerm("D_A", 1), core.ObjectiveTerm("D_B", 1))
        ),
    )


def validate_state(vector, retained=None) -> None:
    vector = np.asarray(vector, dtype=float)
    check(vector.shape == (5,) and np.isfinite(vector).all(), "invalid complete state")
    same(MATRIX @ vector, [0, 0])
    check(bool(np.all(vector >= LOWER - TOLERANCE)), "lower bound violation")
    check(bool(np.all(vector <= UPPER + TOLERANCE)), "upper bound violation")
    if retained is not None:
        check(float(OBJECTIVE @ vector) >= retained - TOLERANCE, "retention violation")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(not args.output.exists(), "refusing to overwrite a verification record")
    before = core_hashes()
    local_module_paths = {}
    for name, module in tuple(sys.modules.items()):
        if name == "carbonscope_flux_core" or name.startswith("carbonscope_flux_core."):
            path = Path(module.__file__).resolve()
            check(path.is_relative_to(CORE.resolve()), f"nonlocal core import: {name}")
            local_module_paths[name] = str(path.relative_to(AREA))
    check(not any(n == "fluxemu" or n.startswith("fluxemu.") for n in sys.modules),
          "external CarbonScope package was imported")

    smoke_stdout = io.StringIO()
    with contextlib.redirect_stdout(smoke_stdout):
        smoke_test.main()

    model = build_control()
    compiled = core.compile_flux_lp(model)
    check(compiled.reaction_ids == REACTIONS, "reaction ordering mismatch")
    check(compiled.balanced_metabolite_ids == ("A", "B"), "balance ordering mismatch")
    same(compiled.balance_matrix, MATRIX)
    same(compiled.lower_bounds, LOWER)
    same(compiled.upper_bounds, UPPER)
    same(compiled.objective, OBJECTIVE)

    fba = core.run_highs_fba(model)
    same(fba.objective_value, 10)
    check(tuple(fba.fluxes.index) == REACTIONS, "FBA ordering mismatch")
    validate_state(fba.fluxes.to_numpy(), 10)
    same(fba.fluxes.iloc[:2], [6, 4])
    same(OBJECTIVE @ fba.fluxes.to_numpy(), fba.objective_value)

    minimisation_model = replace(
        model, objective=core.LinearObjective("minimise", (core.ObjectiveTerm("T", 1),))
    )
    minimum = core.run_highs_fba(minimisation_model)
    same(minimum.objective_value, -2)
    validate_state(minimum.fluxes.to_numpy())
    same(minimum.fluxes.loc["T"], -2)

    regions = []
    for case in CASES:
        fraction = case["fraction"]
        retained = None if fraction is None else 10 * fraction
        for j, endpoint_pair in enumerate(case["witnesses"]):
            for side, witness in enumerate(endpoint_pair):
                validate_state(witness, retained)
                same(witness[j], case["ranges"][j][side])
        reference = core.run_highs_fva_reference(model, fraction)
        fast = core.run_highs_vffva(model, fraction, workers=2, chunk_size=1)
        prepared = core.prepare_highs_flux_region(model, fraction)
        reused = core.run_prepared_highs_vffva(prepared, workers=1, chunk_size=2)
        results = {"reference": reference, "vffva": fast, "prepared_vffva": reused}
        for name, result in results.items():
            check(tuple(result.ranges.index) == REACTIONS, f"{name} ordering mismatch")
            check(tuple(result.ranges.columns) == ("minimum", "maximum"),
                  f"{name} column mismatch")
            same(result.ranges.to_numpy(), case["ranges"])
            same(result.objective_value, 10)
            check(result.fraction_of_optimum == fraction, "retention metadata mismatch")
            check(result.objective_direction == "max", "objective direction mismatch")
        core.assert_fva_equivalent(reference, fast, atol=TOLERANCE, rtol=0)
        core.assert_fva_equivalent(reference, reused, atol=TOLERANCE, rtol=0)
        regions.append({
            **case,
            "observed_ranges": {name: r.ranges.to_numpy().tolist() for name, r in results.items()},
            "maximum_absolute_error": float(max(
                np.max(np.abs(r.ranges.to_numpy() - np.array(case["ranges"])))
                for r in results.values()
            )),
        })

    # Marginal FVA maxima do not in general form a jointly feasible state.
    simultaneous_maxima = np.array(CASES[0]["ranges"], dtype=float)[:, 1]
    residual = MATRIX @ simultaneous_maxima
    check(np.max(np.abs(residual)) > TOLERANCE, "expected incompatible FVA extrema")
    check(before == core_hashes(), "local core changed during verification")
    check(not any(n == "fluxemu" or n.startswith("fluxemu.") for n in sys.modules),
          "external CarbonScope package was imported")

    record = {
        "status": "passed",
        "scope": "Task 46G local solver verification; no AE4 scientific calculations",
        "core_provenance_commit": "a11e25f176a2cc70de22c3f36a1c0e623248393a",
        "core_input_ae4_commit": "3598bd866c33d41d4a6914bda56d17c73900e82a",
        "versions": {"python": platform.python_version(), **{
            p: importlib.metadata.version(p) for p in ("numpy", "pandas", "highspy")
        }},
        "core_sha256": before,
        "verification_script_sha256": sha256(Path(__file__)),
        "imports": local_module_paths,
        "external_carbonscope_imported": False,
        "core_files_unchanged": True,
        "smoke_test_stdout": smoke_stdout.getvalue().strip(),
        "tolerances": {
            "analytical_absolute": TOLERANCE,
            "analytical_relative": 0,
            "highs_primal_and_dual": solver.SOLVER_FEASIBILITY_TOLERANCE,
            "core_primal_validation": solver.FEASIBILITY_TOLERANCE,
            "core_objective_validation": solver.OBJECTIVE_TOLERANCE,
        },
        "control": {
            "metabolite_order": ["A", "B"], "reaction_order": list(REACTIONS),
            "S": MATRIX.tolist(), "lower": LOWER.tolist(), "upper": UPPER.tolist(),
            "maximisation_objective": OBJECTIVE.tolist(),
            "fba_maximum": {"expected": 10, "observed": fba.objective_value,
                            "fluxes": fba.fluxes.to_dict(), "diagnostics": asdict(fba.diagnostics)},
            "fba_minimum_transfer": {"expected": -2, "observed": minimum.objective_value,
                                     "fluxes": minimum.fluxes.to_dict(),
                                     "diagnostics": asdict(minimum.diagnostics)},
            "fva_regions": regions,
            "simultaneous_marginal_maxima": simultaneous_maxima.tolist(),
            "simultaneous_marginal_maxima_residual": residual.tolist(),
        },
        "counts": {
            "local_smoke_tests": 1, "analytical_fba_objectives": 2,
            "fva_regions": 3, "analytical_endpoint_witnesses": 30,
            "implementations_compared_with_analytical_endpoints": 3,
            "endpoint_comparisons_with_analytical_values": 90,
            "endpoint_comparisons_between_reference_and_fast_paths": 60,
            "failures": 0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"status": "passed", **record["counts"],
                      "maximum_absolute_endpoint_error": max(
                          r["maximum_absolute_error"] for r in regions
                      )}, sort_keys=True))


if __name__ == "__main__":
    main()
