"""Execute only the Task 37 permitted existing tests, without a pytest dependency."""
from pathlib import Path
import importlib.util
import os
import sys
import unittest

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]


def load_test_file(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tests" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if __name__ == "__main__":
    nbc = load_test_file("test_task36_minimal_nbc")
    nhe = load_test_file("test_task31_nhe1")
    nbc_names = (
        "test_stoichiometry_is_one_na_two_bicarbonate_and_electrogenic",
        "test_voltage_enters_affinity_with_correct_reversal_sign",
        "test_rest_is_exact_zero_nbc_and_unit_nhe1_multiplier",
        "test_central_stimulus_is_full_nbc_and_source_fixed_2p3_nhe1",
        "test_reference_resting_current_closure_is_exactly_nested",
        "test_stimulated_reference_current_closure_hits_derived_cycle_scale",
    )
    nhe_names = (
        "test_independent_rational_eight_state_equation",
        "test_amount_conversion_and_shared_activity_scale",
        "test_production_uses_exact_table_s1_and_requires_salivary_amount",
    )
    suite = unittest.TestSuite(
        [unittest.FunctionTestCase(getattr(nbc, name)) for name in nbc_names]
        + [nhe.ChaEquationTests(name) for name in nhe_names]
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(0 if result.wasSuccessful() else 1)
