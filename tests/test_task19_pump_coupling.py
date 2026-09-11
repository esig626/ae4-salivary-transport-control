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


@unittest.skipUnless((task.OUT/"results_frozen.json").exists(),"Numerical results not frozen yet")
class Task19ArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from modern_full_model.task19_reporting import decoded_results, table
        cls.manifest,_=task.load_freeze()
        cls.contract,cls.digest=task.contract()
        cls.rows=decoded_results()
        cls.rests=table(task.OUT/"resting_states.csv")
        cls.freeze=task.read(task.OUT/"results_frozen.json")

    def test_contract_precedes_all_runs_and_all_frozen_hashes_match(self):
        receipt=self.freeze
        self.assertEqual(receipt["contract_sha256"],self.digest)
        self.assertEqual(task.git("rev-parse",receipt["contract_commit"]+":"+task.relative(task.OUT/"contract.json")),
            task.git("hash-object",task.relative(task.OUT/"contract.json")))
        for name,expected in receipt["hashes"].items():
            self.assertEqual(task.sha256_file(task.REPO/name),expected,name)
        self.assertTrue(receipt["all_187_frozen_inputs_unchanged"])
        self.assertEqual(len(task.check_inputs()["input_files"]),187)

    def test_complete_panel_includes_every_failed_and_valid_control(self):
        expected={(r,e,s,c) for r in self.manifest["roots"] for e in task.EXPRESSIONS for s in task.PUMPS for c in task.CALCIUM}
        actual={(r["root_id"],r["ae4_expression"],r["pump_scale"],r["calcium_uM"]) for r in self.rows}
        self.assertEqual(len(self.rows),180)
        self.assertEqual(actual,expected)
        self.assertEqual(len(self.rests),60)
        for row in self.rows:
            if not row["dynamic_gate_pass"]:
                self.assertIsNone(row["Q_0_600_pL"])
                self.assertNotEqual(row["status"],"COMPLETE")
            if row["status"]=="NOT_RUN_REST_FAILURE":
                self.assertFalse(row["new_trajectory"])
                self.assertFalse(row.get("trajectory_path"))
        self.assertEqual(sum(r["status"]=="BASELINE_HASH_VALID_REUSE" for r in self.rows),60)

    def test_every_accepted_rest_passes_full_inherited_numerical_gates(self):
        from modern_full_model.calibration import WTCalibrationSpec
        spec=WTCalibrationSpec()
        for row in self.rests:
            if row["resting_gate_pass"]!="True":
                continue
            self.assertEqual(int(row["fixed_parameter_rank"]),10)
            self.assertLessEqual(float(row["max_abs_scaled_independent_rhs"]),spec.root_scaled_tolerance)
            self.assertLessEqual(float(row["max_abs_charge_fmol"]),spec.charge_tolerance_fmol)
            self.assertLessEqual(float(row["max_abs_current_A"]),spec.current_tolerance_A)
            self.assertLessEqual(float(row["max_dimensionless_conservation_ratio"]),1.)
        for root in self.manifest["roots"]:
            for s in task.PUMPS[1:]:
                data=task.read(task.rest_path(root,s))
                self.assertEqual(data["contract_sha256"],self.digest)
                for leg_name in ("pump_leg","ae4_leg"):
                    leg=data[leg_name]
                    for step in leg["steps"]:
                        self.assertEqual(len(step["attempts"]),3)
                        self.assertTrue(step["numerical_point_only"])
                        if step["connected"]:
                            r=step["selected_root"]
                            self.assertTrue(r["passes_numerical_gate"])
                            self.assertEqual(r["normalized_jacobian_rank"],10)
                            self.assertLessEqual(step["distance_from_previous"],.25)
                            self.assertFalse(r["gate_failures"])
                    if not data["pump_leg"]["reached_target"]:
                        self.assertEqual(data["ae4_leg"]["status"],"NOT_RUN_PARENT_REST_FAILURE")
                        self.assertEqual(data["ae4_leg"]["steps"],[])

    def test_every_available_trajectory_has_exact_states_solver_integral_and_conservation(self):
        from modern_full_model.validation import CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU
        for row in self.rows:
            if not row["dynamic_gate_pass"]:
                continue
            path=task.REPO/row["trajectory_path"]
            self.assertEqual(task.sha256_file(path),row["trajectory_sha256"])
            data=task.read(path)
            trace=data["trajectory"]
            model=task.pump_model(self.manifest,row["root_id"],row["calcium_uM"],row["pump_scale"])
            self.assertEqual(data["whole_cell_parameters_sha256"],sha256_object(model.parameters))
            self.assertEqual(data["ae4_parameters_sha256"],sha256_object(model.ae4_parameters))
            self.assertEqual(trace["solver"],asdict(PRODUCTION_RADAU))
            self.assertEqual(trace["time_s"],self.contract["time_grid_s"])
            self.assertTrue(trace["success"] and trace["numerical_gate_pass"])
            self.assertTrue(np.all(np.isfinite(trace["states"])))
            self.assertTrue(np.all(np.asarray(trace["states"])[:12]>0))
            self.assertTrue(np.all(np.asarray(trace["flow_pL_s"])>=0))
            for name,limit in CONSERVATION_RESIDUAL_TOLERANCES.items():
                self.assertLessEqual(trace["max_abs_conservation_residuals"][name],limit)
            self.assertLessEqual(float(row["max_abs_charge_fmol"]),1e-9)
            q=float(np.trapezoid(trace["flow_pL_s"],trace["time_s"]))
            self.assertTrue(math.isclose(q,row["Q_0_600_pL"],rel_tol=1e-12))
            self.assertTrue(math.isclose(q,trace["cumulative_flow_pL"][-1],rel_tol=1e-12))
            if row["pump_scale"]<1:
                rest=task.read(task.rest_path(row["root_id"],row["pump_scale"]))["pump_leg" if row["ae4_expression"]==1 else "ae4_leg"]
                self.assertTrue(np.array_equal(np.asarray(trace["states"])[:12,0],rest["selected_root"]["core_state"]))
                self.assertEqual(data["genotype"],asdict(task.genotype(row["ae4_expression"])))

    def test_all_five_comparisons_match_trajectory_integrals_without_imputation(self):
        from modern_full_model.task19_reporting import table
        recorded=table(task.OUT/"factorial_comparisons.csv")
        rebuilt=task.factorial(self.rows,self.manifest)
        self.assertEqual(len(recorded),90)
        for a,b in zip(recorded,rebuilt):
            for k in task.METRICS:
                self.assertEqual(float(a[k]) if a[k] else None,b[k])
            if b["I"] is not None:
                self.assertTrue(math.isclose(b["I"],b["R_pump_AE4"]/b["R_pump_WT"],rel_tol=1e-12))
                self.assertTrue(math.isclose(b["R_coupled"],b["R_AE4"]*b["R_pump_WT"],rel_tol=1e-12))

    def test_state_flux_diagnostics_include_mandatory_membrane_and_transporter_fields(self):
        from modern_full_model.task19_reporting import table
        rows=table(task.OUT/"state_flux_diagnostics.csv")
        for r in rows:
            if r["phase"] in ("REST","DYNAMIC"):
                for k in ("na_i_mM","k_i_mM","cl_i_mM","ph_i","volume_i_pL","v_apical_mV",
                    "v_basolateral_mV","pump_total_fmol_s","pump_apical_fmol_s","pump_basolateral_fmol_s",
                    "ae4_na_branch_fmol_s","ae4_k_branch_fmol_s","ae4_net_cl_source_fmol_s",
                    "nkcc1_inward_fmol_s","nhe1_inward_fmol_s","ae2_inward_fmol_s"):
                    self.assertTrue(math.isfinite(float(r[k])))
                self.assertTrue(math.isclose(float(r["pump_total_fmol_s"]),float(r["pump_apical_fmol_s"])+float(r["pump_basolateral_fmol_s"]),rel_tol=1e-12))
                e=float(r["ae4_expression"])
                for cation in ("na","k"):
                    actual=float(r[f"ae4_{cation}_branch_fmol_s"])
                    raw=float(r[f"ae4_{cation}_branch_before_expression_fmol_s"])
                    self.assertTrue(math.isclose(actual,e*raw,rel_tol=1e-12))
                    self.assertEqual(actual,-float(r[f"ae4_{cation}_source_fmol_s"]))
                self.assertTrue(math.isclose(float(r["ae4_na_branch_fmol_s"])+float(r["ae4_k_branch_fmol_s"]),float(r["ae4_net_cl_source_fmol_s"]),rel_tol=1e-11,abs_tol=1e-14))
        self.assertEqual(sum(r["phase"]=="REST" for r in rows),sum(r["dynamic_gate_pass"] for r in self.rows))


if __name__=="__main__":
    unittest.main()
