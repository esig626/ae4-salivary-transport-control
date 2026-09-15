"""Regression checks for real intervention and observation failure modes."""
import json
from pathlib import Path
import unittest
import numpy as np
from protocol_layer import *


class ProtocolTests(unittest.TestCase):
    def test_inherited_nesting(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location("task47_readonly", ROOT / "analysis/47_dynamic_potassium_recycling_reconstruction/common.py")
        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
        inherited, seed = mod.load_model()
        new = build_model("CCH_IPR", source_bath=False)
        np.testing.assert_array_equal(inherited.rhs(1., seed), new.rhs(1., seed))

    def test_source_bath(self):
        model = build_model()
        obs = evaluate(model, 0, original_seed(), GENOTYPES["WT"]).diagnostics.observables
        self.assertAlmostEqual(obs.bath_acid_base.hco3_mM, 25., places=8)
        self.assertEqual(model.parameters.bath.cl_mM, 128.3)
        self.assertEqual(model.parameters.bath.k_mM, 4.3)
        self.assertEqual(model.parameters.bath.untracked_osmolyte_mM, 17.)
        self.assertGreater(model.parameters.bath.tic_mM, 25.)

    def test_drug_target_masks(self):
        y = original_seed(); g = GENOTYPES["WT"]
        base = build_model("CCH_IPR")
        for drug, read in [
            ("bumetanide", lambda d: d.homeostasis.nkcc1_inward_fmol_s),
            ("eipa", lambda d: d.homeostasis.nhe1_inward_fmol_s),
            ("t16ainh_a01", lambda d: d.membranes.currents_A["cl_apical"])]:
            m = build_model("CCH_IPR", interventions=Interventions(**{drug: True}))
            out = evaluate(m, 1, y, g)
            self.assertEqual(read(out.diagnostics), 0.)
            self.assertEqual(read(m.evaluate(1, y, genotype=g).diagnostics), 0.)
            np.testing.assert_array_equal(m.rhs(1, y, genotype=g), out.rhs)
            self.assertEqual(m.task48_payload_hash, base.task48_payload_hash)
            self.assertLess(abs(out.rhs[0]+out.rhs[1]-out.rhs[2]-out.rhs[4]), 1e-10)
            if drug != "t16ainh_a01":
                self.assertEqual(m.parameters, base.parameters)
            else:
                self.assertEqual(m.parameters.homeostasis, base.parameters.homeostasis)
        self.assertEqual(g.nkcc1_expression, 1.)
        self.assertEqual(g.nhe1_expression, 1.)

    def test_actual_deletions(self):
        m = build_model("CCH_IPR"); y = original_seed()
        for g in ["AE4_KO", "AE2_KO", "DOUBLE_KO"]:
            d = evaluate(m, 1, y, GENOTYPES[g]).diagnostics
            if g in ["AE4_KO", "DOUBLE_KO"]: self.assertEqual(d.ae4.cl_cell_fmol_s, 0.)
            if g in ["AE2_KO", "DOUBLE_KO"]: self.assertEqual(d.homeostasis.ae2_inward_fmol_s, 0.)

    def test_separate_arms_and_shared_regulation(self):
        for arm, calcium, beta in [("REST", .058, 0), ("CCH_ONLY", .25, 0), ("IPR_ONLY", .058, 1), ("CCH_IPR", .25, 1)]:
            m = build_model(arm)
            self.assertEqual(m.stimulus(1).calcium_uM, calcium)
            self.assertEqual(m.stimulus(1).beta_input, beta)
            self.assertIs(m.regulatory_model.stimulus, m.stimulus)

    def test_unsupported_assays_fail_explicitly(self):
        for name in BLOCKED_PROTOCOLS:
            with self.assertRaises(UnrepresentableProtocol): require_protocol(name)
        with self.assertRaises(UnrepresentableProtocol): build_model(interventions=Interventions(ethoxyzolamide=True))

    def test_cache_separates_genotype_and_protocol(self):
        rest, stim = build_model(), build_model("CCH_IPR")
        keys = {rest_cache_key(m, g) for m in [rest, stim] for g in GENOTYPES.values()}
        self.assertEqual(len(keys), 2*len(GENOTYPES))
        drugkeys={rest_cache_key(build_model("CCH_IPR", interventions=x),GENOTYPES["WT"]) for x in
                  [Interventions(),Interventions(eipa=True),Interventions(bumetanide=True),Interventions(t16ainh_a01=True)]}
        self.assertEqual(len(drugkeys),4)
        self.assertTrue(build_model("NHE_EIPA").task48_interventions.eipa)

    def test_volume_and_optical_observation(self):
        y=np.ones(13); y[2]=60; y[5]=2
        dy=np.zeros(13); dy[2]=3; dy[5]=.2
        self.assertAlmostEqual(concentration_derivative(y,dy), -1.5)
        self.assertLess(spq_inverse_ratio(20,30,quenching_per_mM=.1), 1.)
        self.assertGreater(spq_inverse_ratio(40,30,quenching_per_mM=.1), 1.)
        with self.assertRaises(ValueError): spq_inverse_ratio(20,30,quenching_per_mM=None)
        with self.assertRaises(ValueError): uptake_slope([0,1,2],[1,2,3],window=None)
        self.assertAlmostEqual(uptake_slope([0,1,2],[1,2,3],window=(0,2)), 1.)

    def test_gland_bicarbonate_and_uncertainty(self):
        result=gland_integrals([0,1,2],[2,2,2],[3,3,3])
        self.assertEqual(result, {"water_pL":4.,"bicarbonate_fmol":12.})
        self.assertAlmostEqual(gaussian_ratio_contrast(1.,5.3,1.4,7.1,.8),1.8/np.hypot(.8,1.4))


if __name__ == "__main__":
    result = unittest.TextTestRunner(verbosity=1).run(unittest.defaultTestLoader.loadTestsFromTestCase(ProtocolTests))
    payload={"milestone":"48B","tests_run":result.testsRun,"failures":len(result.failures),"errors":len(result.errors),"passed":result.wasSuccessful(),"inference_runs":0}
    (HERE/"output/protocol_tests.json").write_text(json.dumps(payload,indent=2)+"\n")
    raise SystemExit(not result.wasSuccessful())
