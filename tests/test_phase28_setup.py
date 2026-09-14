"""Seed structure checks only; no model fitting or time integration."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.check_phase28_inputs import check_inputs, git_blob_sha


class Phase28SetupTests(unittest.TestCase):
    def test_git_blob_hash(self):
        self.assertEqual(git_blob_sha(b""), "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391")

    def test_reference_inputs(self):
        result = check_inputs(check_imports=False)
        self.assertEqual(result["reference_parameter_backgrounds"], 10)
        self.assertEqual(set(result["selected_backgrounds"]), {"R09", "R10"})
        self.assertFalse(result["scientific_replay_performed"])

    def test_essential_imports(self):
        result = check_inputs(check_imports=True)
        self.assertEqual(len(result["modules_imported"]), 11)
        self.assertEqual(result["stationary_solves"], 0)
        self.assertEqual(result["stimulated_trajectories"], 0)


if __name__ == "__main__":
    unittest.main()
