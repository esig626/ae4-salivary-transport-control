"""Task 16 conservation, hard freeze, sparse design and selection checks.

No ODE integration or root solve occurs in this test module.
"""
from dataclasses import asdict
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import task16_wt_allocation as task
from modern_full_model.calibration import WTCalibrationSpec, _wt_gate_failures
from modern_full_model.model import WT
from modern_full_model.run_calcium_fast_screen import load_freeze
from modern_full_model.validation import attach_basal_regulation, sha256_file, sha256_object


def flattened(value, prefix=""):
    if isinstance(value, dict):
        return {p: v for k, item in value.items() for p, v in flattened(item, prefix+"."+k).items()}
    return {prefix: value}


class AllocationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest, _ = load_freeze()
        cls.contract, cls.digest = task.contract()
        cls.candidates = cls.contract["candidate_definitions"]

    def test_complete_exact_three_condition_design(self):
        expected = {(r, c) for r in self.manifest["roots"] for c in task.CONDITIONS}
        self.assertEqual(len(self.candidates), 30)
        self.assertEqual({(r["root_id"], r["condition"]) for r in self.candidates}, expected)
        self.assertEqual(self.contract["modified_target_shares"], [.10, .30])
        self.assertEqual(self.contract["calcium_uM"], [.10, .25, .50])
        self.assertEqual({r["routing_family"] for r in self.candidates}, {"AE4NA05", "AE4NA20"})

    def test_illegal_shares_and_nonpositive_loaders_rejected(self):
        for share in [.025, .05, .20, .40, .15]:
            with self.assertRaises(ValueError): task.scales(1, 2, share)
        for a, n in [(-1, 2), (1, 0), (float("nan"), 2), (1, float("inf"))]:
            with self.assertRaises(ValueError): task.scales(a, n, .1)

    def test_reference_loading_pool_and_signed_ae2(self):
        for row in self.candidates:
            if row["target_share"] is None: continue
            base = task.inherited_model(self.manifest, row["root_id"], .10)
            core = self.manifest["roots"][row["root_id"]]["core_state"]
            v = task.partition(base, core)
            self.assertLess(v["ae2_signed_cl_fmol_s"], 0)
            a = row["ae4_capacity_scale"] * v["ae4_cl_fmol_s"]
            n = row["nkcc1_capacity_scale"] * v["nkcc1_cl_fmol_s"]
            self.assertAlmostEqual(a+n, v["positive_pool_fmol_s"], places=14)
            self.assertAlmostEqual(a/(a+n), row["target_share"], places=14)

    def test_exactly_two_parameter_leaves_change(self):
        for row in self.candidates:
            for ca in task.CALCIUM:
                base = task.inherited_model(self.manifest, row["root_id"], ca)
                new = task.allocated_model(self.manifest, row["root_id"], ca,
                    row["ae4_capacity_scale"], row["nkcc1_capacity_scale"])
                old_flat = flattened({"parameters": asdict(base.parameters), "ae4_parameters": asdict(base.ae4_parameters)})
                new_flat = flattened({"parameters": asdict(new.parameters), "ae4_parameters": asdict(new.ae4_parameters)})
                changed = {k for k in old_flat if old_flat[k] != new_flat[k]}
                expected = set() if row["target_share"] is None else {
                    ".parameters.homeostasis.nkcc1_capacity_fmol_s", ".ae4_parameters.carrier_amount_fmol"}
                self.assertEqual(changed, expected)
                self.assertEqual(sha256_object(base.regulatory_model), sha256_object(new.regulatory_model))
                self.assertEqual(asdict(base.stimulus), asdict(new.stimulus))

    def test_complete_ae4_vector_and_nkcc_cycle_at_rest_and_onset(self):
        for row in self.candidates:
            if row["target_share"] is None: continue
            base = task.inherited_model(self.manifest, row["root_id"], .50)
            new = task.allocated_model(self.manifest, row["root_id"], .50,
                row["ae4_capacity_scale"], row["nkcc1_capacity_scale"])
            state = attach_basal_regulation(base, self.manifest["roots"][row["root_id"]]["core_state"])
            for t in [0., 1e-6, 300.]:
                old, actual = base.evaluate(t, state), new.evaluate(t, state)
                ao, an = old.diagnostics.ae4, actual.diagnostics.ae4
                names = ["na_cell_fmol_s", "k_cell_fmol_s", "cl_cell_fmol_s", "hco3_cell_fmol_s"]
                np.testing.assert_allclose([getattr(an, k) for k in names],
                    np.asarray([getattr(ao, k) for k in names])*row["ae4_capacity_scale"], rtol=1e-12, atol=1e-14)
                for k in ["source_tic_i_fmol_s", "source_alkalinity_i_fmol_s"]:
                    self.assertAlmostEqual(an.diagnostics[k], row["ae4_capacity_scale"]*ao.diagnostics[k], places=13)
                for k in ["affinity_na", "affinity_k", "effective_rates", "occupancy_inward", "occupancy_outward"]:
                    self.assertEqual(an.diagnostics[k], ao.diagnostics[k])
                self.assertLess(an.diagnostics["local_detailed_balance_residual"], 1e-12)
                self.assertGreaterEqual(an.diagnostics["entropy_production_over_r_fmol_s"], 0)
                delta = np.zeros_like(old.rhs)
                delta[:5] = (row["ae4_capacity_scale"]-1)*np.asarray([
                    ao.na_cell_fmol_s, ao.k_cell_fmol_s, ao.cl_cell_fmol_s,
                    ao.hco3_cell_fmol_s, ao.hco3_cell_fmol_s])
                dn = (row["nkcc1_capacity_scale"]-1)*old.diagnostics.homeostasis.nkcc1_inward_fmol_s
                delta[:3] += dn*np.array([1,1,2])
                np.testing.assert_allclose(actual.rhs-old.rhs, delta, rtol=1e-10, atol=2e-14)
                self.assertTrue(all(abs(v) <= task.CONSERVATION_RESIDUAL_TOLERANCES[k]
                    for k, v in actual.diagnostics.conservation_residuals.items()))

    def test_scale_one_reproduces_baseline_exactly(self):
        for root in self.manifest["roots"]:
            old = task.inherited_model(self.manifest, root, .10)
            new = task.allocated_model(self.manifest, root, .10, 1., 1.)
            state = attach_basal_regulation(old, self.manifest["roots"][root]["core_state"])
            np.testing.assert_array_equal(old.evaluate(0., state).rhs, new.evaluate(0., state).rhs)

    def test_all_187_frozen_inputs_unchanged(self):
        self.assertEqual(len(task.check_inputs()["input_files"]), 187)
        self.assertEqual(task.git("diff", "--name-only", task.BASE, "HEAD", "archive", "model"), "")

    def test_wt_contract_generation_never_reads_genotype_outcomes(self):
        original_read = Path.read_text
        def guarded_read(path, *args, **kwargs):
            s = str(path)
            forbidden = ("ae2_null_blind_predictions.csv", "ae4_5pct_results.csv", "allocation_vs_genotype.csv")
            if any(x in s for x in forbidden) or ("14B_five_percent_ae4_check/trajectories" in s):
                raise AssertionError("WT contract attempted to read genotype outcome")
            return original_read(path, *args, **kwargs)
        ledger = (task.OUT / "baseline_input_hashes.json").read_bytes()
        with TemporaryDirectory() as tmp, patch.object(task, "OUT", Path(tmp)), patch.object(Path, "read_text", guarded_read):
            (Path(tmp)/"baseline_input_hashes.json").write_bytes(ledger)
            task.prepare()
            self.assertEqual(task.contract()[0]["candidate_definitions"], self.candidates)

    def test_no_extra_scientific_candidates_or_unlicensed_dynamics(self):
        with (task.OUT/"wt_rest_candidates.csv").open() as f: rest=list(csv.DictReader(f))
        with (task.OUT/"wt_dynamic_candidates.csv").open() as f: dynamic=list(csv.DictReader(f))
        self.assertEqual(len(rest), 30)
        self.assertEqual(len(dynamic), 90)
        self.assertEqual({r["condition"] for r in dynamic}, set(task.CONDITIONS))
        for r in dynamic:
            if r["condition"] != "inherited_baseline":
                self.assertEqual(r["status"], "NOT_RUN_REST_REJECTED")
                self.assertEqual(r["trajectory_path"], "")
        self.assertFalse(list((task.OUT/"wt_trajectories").glob("*.json")))

    def test_saved_continuation_failures_have_physical_boundary_evidence(self):
        for c in self.candidates:
            if c["target_share"] is None: continue
            p = task.read(task.root_path(c))
            self.assertEqual(p["row"]["status"], "NUMERICAL_CONTINUATION_FAILURE")
            self.assertIsNone(p["core_state"])
            last = next(s for s in reversed(p["continuation_steps"]) if s["connected"])
            self.assertEqual(last["selected_root"]["normalized_jacobian_rank"], 10)
            self.assertGreater(last["selected_root"]["coordinates"][0], 35)
            self.assertLess(last["selected_root"]["coordinates"][1], 100)
            self.assertTrue(all(abs(a["coordinates"][0]-60) < 1e-7 for a in p["continuation_steps"][-1]["attempts"]))

    def test_genotype_fields_cannot_change_wt_physiology_gate(self):
        for root, r in self.manifest["roots"].items():
            obs = r["native_root_object"]["observables"]
            first = _wt_gate_failures(obs, WTCalibrationSpec())
            for effect in [-100., 0., .35, 100.]:
                altered = {**obs, "D_AE4_5pct": effect, "D_AE2": -effect}
                self.assertEqual(_wt_gate_failures(altered, WTCalibrationSpec()), first)


if __name__ == "__main__":
    unittest.main()
