"""Record the prescribed Task 50B fixed state tests; no production integration."""
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"


def main():
    attempts = sorted(OUT.glob("verification_attempt_*.json"))
    if any(json.loads(p.read_text())["status"] == "PASS" for p in attempts):
        raise SystemExit("Successful verification is frozen; no optional repeat")
    attempt = len(attempts) + 1
    path = ROOT / "tests/test_task50_effective_coupling.py"
    spec = importlib.util.spec_from_file_location("task50_verification_tests", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    count = {"core_evaluations_including_constructor_checks": 0}
    original = module.ModernFullModel.evaluate

    def counted(model, *args, **kwargs):
        count["core_evaluations_including_constructor_checks"] += 1
        return original(model, *args, **kwargs)

    stream = io.StringIO()
    with patch.object(module.ModernFullModel, "evaluate", counted):
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
            unittest.defaultTestLoader.loadTestsFromModule(module)
        )
    stem = f"verification_attempt_{attempt:02d}"
    (OUT / (stem + ".log")).write_text(stream.getvalue())
    files = [path, ROOT / "src/modern_full_model/task50_effective_coupling.py",
             Path(__file__)]
    record = {
        "checkpoint": "50B", "attempt": attempt,
        "status": "PASS" if result.wasSuccessful() else "FAIL",
        "classification": "TARGET-CALIBRATED CONSTRUCTION",
        "tests_run": result.testsRun, "test_failures": len(result.failures),
        "test_errors": len(result.errors), **count, **module.EVIDENCE,
        "new_production_trajectories": 0, "new_stationary_solves": 0,
        "new_parameter_fits_or_searches": 0,
        "tolerance_source": "DECISION_LOG.md D50-07, published before implementation",
        "source_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "error_bound_scope": "Fixed state comparisons only; no global trajectory error bound",
    }
    (OUT / (stem + ".json")).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: record[k] for k in ("status", "tests_run", "test_failures", "test_errors", "core_evaluations_including_constructor_checks", "max_absolute_differences")}))
    if not result.wasSuccessful():
        print(stream.getvalue())
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
