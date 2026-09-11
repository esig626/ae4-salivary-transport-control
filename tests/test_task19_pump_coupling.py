"""Focused factorial/scaling tests; no new roots or ODEs are solved here."""
from dataclasses import asdict
import csv
import math
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model import task19_pump_coupling as task
from modern_full_model.validation import attach_basal_regulation, sha256_object


class Task19ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest,_=task.load_freeze()
        cls.contract,_=task.contract()

    def test_exact_panel_and_no_phenotype_selection(self):
        c=self.contract
        self.assertEqual(c["pump_scales"],[1.0,.50,.10])
        self.assertEqual(c["ae4_expressions"],[1.,.05])
        self.assertEqual(c["calcium_uM"],[.10,.25,.50])
        self.assertEqual(c["root_ids"],sorted(self.manifest["roots"]))
        self.assertEqual(len(c["root_ids"]),10)
        self.assertEqual(set(c["routing_family_by_root"].values()),{"AE4NA05","AE4NA20"})
        self.assertEqual(c["chloride_allocation"],"inherited_baseline_only")
        self.assertTrue(c["interpretation"]["no_pump_scales_added_after_results"])
        self.assertFalse(c["interpretation"]["experimental_phenotype_used_for_selection"])
        self.assertFalse(c["interpretation"]["continuous_coupling_fit"])
        self.assertFalse(c["interpretation"]["coupled_reduction_alone_is_success"])

    def test_exact_total_and_both_membrane_capacity_scaling_with_all_else_frozen(self):
        for root in self.manifest["roots"]:
            for ca in task.CALCIUM:
                base=task.build_model(self.manifest,root,ca)
                original=asdict(base.parameters)
                for s in task.PUMPS:
                    model=task.pump_model(self.manifest,root,ca,s)
                    actual=asdict(model.parameters)
                    expected=asdict(base.parameters)
                    expected["membranes"]["nak_capacity_fmol_s"] *= s
                    self.assertEqual(actual,expected)
                    p,q=base.parameters.membranes,model.parameters.membranes
                    for f in (p.apical_pump_fraction,1-p.apical_pump_fraction):
                        self.assertAlmostEqual(q.nak_capacity_fmol_s*f,p.nak_capacity_fmol_s*f*s,places=14)
                    self.assertEqual(sha256_object(model.ae4_parameters),sha256_object(base.ae4_parameters))
                    self.assertEqual(sha256_object(model.regulatory_model),sha256_object(base.regulatory_model))
                    self.assertEqual(sha256_object(model.stimulus),sha256_object(base.stimulus))
                    self.assertIs(model.ae4_evaluator,base.ae4_evaluator)
                self.assertEqual(asdict(base.parameters),original)

    def test_pump_law_flux_scaling_and_independent_3na_2k_source_accounting(self):
        for root in self.manifest["roots"]:
            base=task.pump_model(self.manifest,root,.10,1.)
            state=attach_basal_regulation(base,self.manifest["roots"][root]["core_state"])
            old=base.evaluate(0.,state,genotype=task.WT).diagnostics.membranes
            for s in task.PUMPS:
                model=task.pump_model(self.manifest,root,.10,s)
                m=model.evaluate(0.,state,genotype=task.WT).diagnostics.membranes
                self.assertAlmostEqual(m.pump_apical_fmol_s,old.pump_apical_fmol_s*s,places=14)
                self.assertAlmostEqual(m.pump_basolateral_fmol_s,old.pump_basolateral_fmol_s*s,places=14)
                conversion=1e15/model.parameters.constants.faraday_C_mol
                jka=m.currents_A["k_apical"]*conversion
                jkb=m.currents_A["k_basolateral"]*conversion
                jp=m.pump_apical_fmol_s+m.pump_basolateral_fmol_s
                self.assertAlmostEqual(m.cell_sources_fmol_s["na"],-3*jp,places=13)
                self.assertAlmostEqual(m.cell_sources_fmol_s["k"]+jka+jkb,2*jp,places=13)
                self.assertAlmostEqual(m.lumen_sources_fmol_s["na"]+m.currents_A["para_na"]*conversion,3*m.pump_apical_fmol_s,places=13)
                self.assertAlmostEqual(m.lumen_sources_fmol_s["k"]-jka+m.currents_A["para_k"]*conversion,-2*m.pump_apical_fmol_s,places=13)

    def test_all_sixty_scale_one_outputs_are_hash_valid_and_reused_without_solve(self):
        with patch.object(task,"simulate_genotype",side_effect=AssertionError("Forbidden baseline recomputation")):
            for root in self.manifest["roots"]:
                for e in task.EXPRESSIONS:
                    for ca in task.CALCIUM:
                        row=task.dynamic_job((root,e,1.,ca))
                        self.assertEqual(row["status"],"BASELINE_HASH_VALID_REUSE")
                        self.assertFalse(row["new_trajectory"])
                        self.assertTrue(row["dynamic_gate_pass"])

    def test_corrupt_reuse_hash_is_rejected(self):
        root=next(iter(self.manifest["roots"]))
        with patch.object(task,"sha256_file",return_value="bad"):
            with self.assertRaises(AssertionError):
                task.verify_reuse(self.manifest,root,.05,.10)

    def test_factorial_arithmetic_distinguishes_pump_only_and_interaction(self):
        # A severe pump-only reduction is not an interaction: AE4 still raises flow 20%.
        self.assertEqual(task.ratios(100.,120.,20.,24.),
            {"R_AE4":1.2,"R_pump_WT":.2,"R_pump_AE4":.2,"I":1.,"R_coupled":.24})
        values=task.ratios(100.,120.,20.,12.)
        self.assertEqual(values["R_AE4"],.6)
        self.assertEqual(values["I"],.5)
        self.assertEqual(values["R_coupled"],.12)
        missing=task.ratios(100.,120.,20.,None)
        self.assertEqual(missing["R_pump_WT"],.2)
        self.assertIsNone(missing["R_AE4"])
        self.assertIsNone(missing["I"])

    def test_complete_factorial_accounting_requires_intact_pump_controls(self):
        rows=[{"root_id":r,"ae4_expression":e,"pump_scale":s,"calcium_uM":ca,
            "dynamic_gate_pass":True,"Q_0_600_pL":1.,"status":"MOCK"}
            for r in self.manifest["roots"] for e in task.EXPRESSIONS for s in task.PUMPS for ca in task.CALCIUM]
        self.assertEqual(len(task.factorial(rows,self.manifest)),90)
        with self.assertRaises(AssertionError):
            task.factorial(rows[1:],self.manifest)
        with self.assertRaises(AssertionError):
            task.factorial(rows+[rows[0]],self.manifest)


if __name__=="__main__":
    unittest.main()
