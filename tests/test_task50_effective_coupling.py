"""Task 50 fixed state verification only; never integrate or refit a model."""
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from modern_full_model.ae4_cacc_recruitment import recruitment_factor
from modern_full_model.model import AE2_NULL, AE4_NULL, WT, ConstantStimulus, ModernFullModel
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model
from modern_full_model.task41_selected import Task41SelectedModel, SELECTED_NULL_CACC_RECRUITMENT
from modern_full_model.task50_effective_coupling import (
    FIXED_LAMBDA, Task50EffectiveCouplingModel, effective_coupling_factor,
)
from modern_full_model.validation import CONSERVATION_RESIDUAL_TOLERANCES, StimulusArm, sha256_object

OUT41 = ROOT / "results/41_ae4_loss_algebraic_design/candidate_08"
EVIDENCE = {"max_absolute_differences": {}, "comparisons": [], "exact_parent_comparisons": 0}


def with_protocol(template, protocol):
    core = template.base_model.base_model
    model = ModernFullModel(
        core.parameters, stimulus=protocol, regulatory_model=core.regulatory_model,
        ae4_parameters=core.ae4_parameters, ae4_evaluator=core.ae4_evaluator,
        nkcc1_kinetics=core.nkcc1_kinetics,
    )
    return MinimalNbcModel(StimulatedNkcc1Model(model), template.nbc_parameters)


class Task50Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / "analysis/40_ae4_equal_cation_routing/validation_common.py"
        spec = importlib.util.spec_from_file_location("task50_parent_fixtures", path)
        cls.parent = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.parent)
        # The only difference in the 150 file historical intake is the current
        # binding AGENTS.md. Pin it explicitly rather than restore old rules.
        intake = json.loads((path.parent / "parent_source_manifest.json").read_text())
        for item in intake["files"]:
            raw = (ROOT / item["path"]).read_bytes()
            actual = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
            expected = ("94465907d246d530e78f18cf1237fad5088cfa51"
                        if item["path"] == "AGENTS.md" else item["git_blob_sha1"])
            assert actual == expected, item["path"]
        out40 = ROOT / "results/40_ae4_equal_cation_routing"
        frozen = json.loads((out40 / "frozen_inputs.json").read_text())
        for item in frozen["task40_execution_files"]:
            assert hashlib.sha256((ROOT / item["path"]).read_bytes()).hexdigest() == item["sha256"], item["path"]
        _, _, _, cls.rest, cls.combined, _ = cls.parent.models_and_reference()
        saved_rest = json.loads((out40 / "wt_rest.json").read_text())
        assert saved_rest["status"] == "PASS" and saved_rest["admissible"]
        assert saved_rest["state_sha256"] == sha256_object(saved_rest["state_vector"]) == frozen["frozen_rest_state_sha256"]
        assert cls.parent.parameter_hashes(cls.rest) == frozen["rest_parameter_hashes"]
        assert cls.parent.parameter_hashes(cls.combined) == frozen["stimulus_parameter_hashes"]
        cls.y0 = np.asarray(saved_rest["state_vector"])
        cls.saved = {case: json.loads((OUT41 / (case + "_verification.json")).read_text())
                     for case in ("wt", "ae4_5pct", "ae4_null")}
        cls.cch = with_protocol(cls.combined, replace(cls.combined.stimulus, arm=StimulusArm.CCH_ONLY))
        cls.ipr = with_protocol(cls.combined, replace(cls.combined.stimulus, arm=StimulusArm.IPR_ONLY))
        cls.fractional = with_protocol(cls.combined, ConstantStimulus(0.13, 0.4))

    def close(self, label, first, second, atol, rtol=2e-12):
        a, b = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
        diff = float(np.max(np.abs(a - b)))
        maximum = EVIDENCE["max_absolute_differences"]
        maximum[label] = max(maximum.get(label, 0.0), diff)
        np.testing.assert_allclose(a, b, atol=atol, rtol=rtol, err_msg=label)

    def exact_parent(self, base, time, state, genotype):
        new = Task50EffectiveCouplingModel(base)
        old, got = base.evaluate(time, state, genotype=genotype), new.evaluate(time, state, genotype=genotype)
        self.assertIs(new.effective_model(1.0), base)
        np.testing.assert_array_equal(got.rhs, old.rhs)
        actual, expected = asdict(got.diagnostics), asdict(old.diagnostics)
        actual["regulatory"] = {k: v for k, v in actual["regulatory"].items() if not k.startswith("task50_")}
        self.assertEqual(actual, expected)
        self.assertEqual(got.diagnostics.regulatory["task50_effective_coupling_factor"], 1.0)
        EVIDENCE["exact_parent_comparisons"] += 1

    def test_01_fixed_constant_factor_boundaries_and_domains(self):
        self.assertEqual(FIXED_LAMBDA, 0.89488127156712)
        for a in (0.0, 0.37, 1.0):
            for beta in (0.0, 0.4, 1.0):
                self.assertEqual(effective_coupling_factor(a, beta, 1.0), 1.0)
        for e in (0.0, 0.05, 1.0):
            self.assertEqual(effective_coupling_factor(0.0, 1.0, e), 1.0)
            self.assertEqual(effective_coupling_factor(1.0, 0.0, e), 1.0)
        for inputs in ((-0.1, 1, 0), (1.1, 1, 0), (1, -0.1, 0),
                       (1, 1.1, 0), (1, 1, -0.1), (float("nan"), 1, 0),
                       (1, float("inf"), 0), (1, 1, float("nan"))):
            with self.assertRaises(ValueError): effective_coupling_factor(*inputs)

    def test_02_factor_numerical_equivalence_to_frozen_coefficient(self):
        for a in (0.0, 0.37, 1.0):
            for e in (0.0, 0.05, 1.0):
                self.close("factor", effective_coupling_factor(a, 1.0, e),
                           recruitment_factor(a, e, SELECTED_NULL_CACC_RECRUITMENT), 8e-15, 0)
        self.assertNotEqual(effective_coupling_factor(1.0, 1.0, 0.0),
                            SELECTED_NULL_CACC_RECRUITMENT)

    def test_03_wt_and_ae2_exact_parent_nesting(self):
        for base in (self.rest, self.cch, self.ipr, self.combined, self.fractional):
            for genotype in (WT, AE2_NULL):
                self.exact_parent(base, 1e-6, self.y0, genotype)

    def test_04_rest_and_ipr_exact_nesting_in_mutants(self):
        for base in (self.rest, self.ipr):
            for genotype in (AE4_NULL, replace(WT, ae4_expression=0.05)):
                self.exact_parent(base, 1e-6, self.y0, genotype)
        self.exact_parent(self.combined, 0.0, self.y0, AE4_NULL)

    def test_05_cch_only_full_rhs_and_diagnostics_nest_parent(self):
        states = (self.y0, np.asarray(self.saved["ae4_null"]["last_checked_state_vector"]))
        for state in states:
            for genotype in (WT, replace(WT, ae4_expression=0.05), AE4_NULL):
                self.exact_parent(self.cch, 1e-6, state, genotype)

    def test_06_beta_one_full_rhs_current_and_source_equivalence(self):
        new, old = Task50EffectiveCouplingModel(self.combined), Task41SelectedModel(self.combined)
        for case, e in (("wt", 1.0), ("ae4_5pct", 0.05), ("ae4_null", 0.0)):
            genotype = replace(WT, ae4_expression=e)
            states = ((1e-6, self.y0), (600.0, np.asarray(self.saved[case]["last_checked_state_vector"])))
            for time, state in states:
                a, b = new.evaluate(time, state, genotype=genotype), old.evaluate(time, state, genotype=genotype)
                self.close("rhs_native_units", a.rhs, b.rhs, 2e-13)
                ma, mb = asdict(a.diagnostics.membranes), asdict(b.diagnostics.membranes)
                for key in ma:
                    va, vb = ma[key], mb[key]
                    if isinstance(va, dict):
                        self.assertEqual(va.keys(), vb.keys())
                        va, vb = list(va.values()), list(vb.values())
                    unit = "A" if key.endswith("_A") else "V" if key.endswith("_V") else "fmol_s"
                    self.close("membrane_" + unit, va, vb, {"A":2e-25,"V":2e-14,"fmol_s":2e-13}[unit])
                for key, value in a.diagnostics.conservation_residuals.items():
                    unit = "A" if key.endswith("_A") else "native_units"
                    self.close("conservation_" + unit, value, b.diagnostics.conservation_residuals[key],
                               2e-25 if unit == "A" else 2e-13)
                for key, limit in CONSERVATION_RESIDUAL_TOLERANCES.items():
                    self.assertLessEqual(abs(a.diagnostics.conservation_residuals[key]), limit, key)
                gnew = a.diagnostics.regulatory["task50_effective_cl_conductance_S"]
                gold = b.diagnostics.regulatory["effective_cacc_conductance_S"]
                self.close("conductance_S", gnew, gold, 1e-24)
                # Both chloride amount rows, including unchanged paracellular
                # and outflow bookkeeping, are covered by the full RHS above.
                flux = -a.diagnostics.membranes.cell_sources_fmol_s["cl"]
                pump = a.diagnostics.membranes.pump_apical_fmol_s + a.diagnostics.membranes.pump_basolateral_fmol_s
                nhe = a.diagnostics.homeostasis.nhe1_inward_fmol_s
                self.assertAlmostEqual(flux, 6*pump-nhe+2*a.rhs[0]-a.rhs[4]-a.rhs[2], delta=1e-12)
                EVIDENCE["comparisons"].append({"case":case,"time_s":time,
                    "max_abs_rhs_difference":float(np.max(np.abs(a.rhs-b.rhs)))})

    def test_07_only_declared_conductance_changes_at_fixed_state(self):
        new = Task50EffectiveCouplingModel(self.combined)
        for case, expression in (("ae4_5pct", 0.05), ("ae4_null", 0.0)):
            state = np.asarray(self.saved[case]["last_checked_state_vector"])
            genotype = replace(WT, ae4_expression=expression)
            a, b = new.evaluate(600, state, genotype=genotype), self.combined.evaluate(600, state, genotype=genotype)
            for field in ("ae4", "homeostasis", "water", "observables", "co2_fluxes_fmol_s", "outflow_sources_fmol_s", "state_charge_fmol"):
                self.assertEqual(getattr(a.diagnostics, field), getattr(b.diagnostics, field), field)
            factor = a.diagnostics.regulatory["task50_effective_coupling_factor"]
            effective = new.effective_model(factor)
            modified, original = asdict(effective.parameters), asdict(self.combined.parameters)
            modified["membranes"]["g_cl_apical_S"] = original["membranes"]["g_cl_apical_S"]
            self.assertEqual(modified, original)
            for field in ("ae4_parameters", "ae4_evaluator", "nkcc1_kinetics", "nbc_parameters", "regulatory_model", "stimulus", "state_names"):
                self.assertEqual(getattr(effective,field), getattr(self.combined,field),field)
            self.assertEqual(new.state_names, self.combined.state_names)

    def test_08_existing_beta_input_controls_the_new_term(self):
        new = Task50EffectiveCouplingModel(self.fractional)
        a = new.evaluate(1.0, self.y0, genotype=AE4_NULL)
        r = a.diagnostics.regulatory
        self.assertEqual(r["task50_beta_input"], 0.4)
        # The inherited AE4 regulatory state is initially zero, but the new
        # factor uses the positive protocol beta input without a new ODE.
        self.assertLess(r["task50_effective_coupling_factor"], 1.0)
        full = effective_coupling_factor(r["task50_calcium_activation"], 1.0, 0.0)
        self.assertAlmostEqual(r["task50_effective_coupling_factor"],
                               1.0 + 0.4 * (full - 1.0), delta=2e-16)
        self.assertEqual(r["task50_classification"], "TARGET-CALIBRATED CONSTRUCTION")

    def test_09_public_integrator_dispatch_uses_amended_rhs_without_integration(self):
        new = Task50EffectiveCouplingModel(self.combined)
        with patch("modern_full_model.task41_selected.solve_ivp") as solve:
            new.solve_dynamics((1e-6, 600), initial_state=self.y0, genotype=AE4_NULL)
            callback = solve.call_args.args[0]
            np.testing.assert_array_equal(callback(1e-6, self.y0), new.rhs(1e-6, self.y0, genotype=AE4_NULL))
            self.assertFalse(np.array_equal(callback(1e-6, self.y0), self.combined.rhs(1e-6, self.y0, genotype=AE4_NULL)))
        self.assertIs(Task50EffectiveCouplingModel.solve_dynamics, Task41SelectedModel.solve_dynamics)

    def test_10_frozen_task41_execution_inputs_unchanged(self):
        manifest = json.loads((OUT41 / "frozen_inputs.json").read_text())
        for item in manifest["execution_files"]:
            self.assertEqual(hashlib.sha256((ROOT / item["path"]).read_bytes()).hexdigest(), item["sha256"])
        self.assertEqual(manifest["selected_parameters"]["null_stimulated_fraction"], SELECTED_NULL_CACC_RECRUITMENT)


if __name__ == "__main__":
    unittest.main(verbosity=2)
