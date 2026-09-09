"""Regression tests for the forensic historical-equation evaluation.

The BDF trajectory tests validate the translated harness and its frozen output;
they are not evidence that MATLAB/ode15s was executed in this environment.
"""

from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / "analysis/11_forensic_reconstruction/tools/reproduce_historical_matlab.py"
SPEC = importlib.util.spec_from_file_location("reproduce_historical_matlab", HELPER)
assert SPEC and SPEC.loader
REPRO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPRO)


class HistoricalReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result_paths = (
            ROOT / "results/11_forensic_reconstruction/reproduction_summary.json",
            ROOT / "results/11_forensic_reconstruction/reproduction_trajectories.csv",
            ROOT / "results/11_forensic_reconstruction/reproduction_residuals.csv",
        )
        cls.before = {path: path.read_bytes() for path in cls.result_paths}
        cls.summary = REPRO.write_outputs()

    def test_archive_source_hashes_and_forensic_copies(self) -> None:
        expected = {
            "Original/Parameters.m": "fa19bfabd63f0798bbbf72e269a4d469cf00b2b41434bc7d32f614163db20eef",
            "Original/Saliva_Ae4.m": "7d72a055ab240f658287643520558a592a618845945c01bbdc4cb454f84fdcf3",
            "Par.m": "06ae901e1c045e234ec08d808bf909539b5d6d8cf009ebd98797791b7d313eff",
            "Par.mat": "fea3bdbe49dd794d43204bb234bd2b31e6215d784d920f46d3b4300f7068e1ca",
            "Salivary.m": "a6524879c8d948d8be41cdc081c2423b1bda8510943ddd3ce76fab41dc87ed8e",
            "Salivary2.m": "94e8ddc12e26d4cac0437f3e3d946014c91a1dcffa8fac6072aed56ac010b923",
            "Salivary_ex.m": "5f3ecaa2473a3bb409a1ff3dac980dcf356897f42923c50d092ae847f1c6f760",
        }
        records = self.summary["archive_working_copy"]["files"]
        self.assertEqual(len(records), len(expected))
        for record in records:
            self.assertEqual(record["source_sha256"], expected[record["path"]])
            self.assertEqual(record["copy_sha256"], record["source_sha256"])

    def test_original_calibration_is_an_exact_constructed_steady_state(self) -> None:
        calibration = self.summary["original_calibration"]
        self.assertLessEqual(calibration["max_abs_residual"], 1.4e-17)
        self.assertAlmostEqual(calibration["fluxes"]["q_a"], calibration["fluxes"]["q_b"], delta=2e-18)
        self.assertAlmostEqual(calibration["CO2_i"], 6.600123698596884, delta=2e-14)
        self.assertAlmostEqual(calibration["Psi_l"], 48.800135723674394, delta=2e-13)

    def test_original_calibration_generates_every_par_mat_field(self) -> None:
        concordance = self.summary["original_to_par_concordance"]
        self.assertEqual(concordance["fields_compared"], 45)
        self.assertLess(concordance["max_abs_relative_difference"], 2.6e-15)
        fields = concordance["fields"]
        self.assertAlmostEqual(fields["a2"]["par_mat"], 2.0096e-5, delta=1e-19)
        self.assertAlmostEqual(fields["a4"]["par_mat"], 1.3852e-6, delta=1e-20)
        self.assertAlmostEqual(fields["k2"]["par_mat"], 0.00159305881398359, delta=2e-18)

    def test_outer_wt_initial_state_is_a_numerical_fixed_point(self) -> None:
        seven = self.summary["runs"]["salivary_7_state"]["wt"]
        five = self.summary["runs"]["salivary2_5_state"]["wt"]
        self.assertLess(max(abs(value) for value in seven["initial_raw_rhs"].values()), 5.3e-12)
        self.assertLess(max(abs(value) for value in seven["initial_scaled_rhs"].values()), 5.3e-8)
        self.assertLess(max(abs(value) for value in five["initial_raw_rhs"].values()), 5.3e-12)
        diagnostic = seven["initial_diagnostics"]
        self.assertAlmostEqual(diagnostic["Va"], -50.24, delta=1e-12)
        self.assertAlmostEqual(diagnostic["Vb"], -62.8, delta=1e-12)
        self.assertAlmostEqual(diagnostic["pH"], 6.91, delta=3e-10)
        self.assertAlmostEqual(diagnostic["CO2"], 6.600123698596897, delta=2e-14)
        self.assertAlmostEqual(diagnostic["cell_volume_pL_from_delta_Hi"], 1.3, delta=2e-14)
        self.assertAlmostEqual(diagnostic["q_total"], 0.0005314154767606997, delta=2e-16)

    def test_calibration_conventions_behind_the_resting_row(self) -> None:
        fields = self.summary["original_to_par_concordance"]["fields"]
        b1 = fields["b1"]["par_mat"]
        b2 = fields["b2"]["par_mat"]
        fluxes = self.summary["original_calibration"]["fluxes"]
        self.assertAlmostEqual(b1 / b2, 10.62857142857142, delta=2e-14)
        self.assertAlmostEqual(fluxes["q_a"] / b1, 0.7997659713233467, delta=2e-15)
        self.assertAlmostEqual(fluxes["q_b"] / b2, 8.50036975235099, delta=2e-14)
        self.assertEqual(fields["pco2"]["par_mat"], 1970.0)
        self.assertAlmostEqual(fields["gb"]["par_mat"], 0.199652558943136, delta=1e-15)

    def test_translated_bdf_scenario_flow_comparisons(self) -> None:
        expected = {
            "salivary_7_state": {
                "ae2_knockout": (0.9990522652444069, 0.9973818146227837),
                "ae4_knockout": (0.9985561424848957, 1.101876915912234),
            },
            "salivary2_5_state": {
                "ae2_knockout": (1.000000071188472, 1.0020228733377807),
                "ae4_knockout": (1.0590191789778476, 0.9032971494353547),
            },
        }
        for version, scenarios in expected.items():
            for scenario, (endpoint, first_sample_peak) in scenarios.items():
                comparison = self.summary["runs"][version][scenario]["comparison_to_wt"]
                self.assertAlmostEqual(comparison["endpoint_flow_ratio"], endpoint, delta=2e-8)
                self.assertAlmostEqual(comparison["peak_sample_flow_ratio"], first_sample_peak, delta=2e-8)
                peak = self.summary["runs"][version][scenario]["stimulated_peak_sample"]
                self.assertEqual(peak["time"], 101.0)

    def test_archived_protocol_endpoint_rhs_is_reported_without_steady_state_claim(self) -> None:
        expected_scaled = {
            "salivary_7_state": {
                "wt": 3.4889832451401864e-5,
                "ae2_knockout": 5.578765843887103e-5,
                "ae4_knockout": 5.7138150520742155e-5,
            },
            "salivary2_5_state": {
                "wt": 3.5180549060285156e-5,
                "ae2_knockout": 5.5840089059987404e-5,
                "ae4_knockout": 0.07222091282319474,
            },
        }
        for version, scenarios in expected_scaled.items():
            for scenario, expected in scenarios.items():
                actual = self.summary["runs"][version][scenario]["stimulated_endpoint_t200"][
                    "max_abs_scaled_rhs"
                ]
                self.assertAlmostEqual(actual, expected, delta=2e-10)
        for run in self.summary["runs"]["salivary_7_state"].values():
            self.assertLess(run["stimulated_endpoint_t200"]["max_abs_raw_rhs"], 6e-9)
        five_ae4 = self.summary["runs"]["salivary2_5_state"]["ae4_knockout"]
        self.assertGreater(five_ae4["stimulated_endpoint_t200"]["max_abs_raw_rhs"], 7e-5)
        self.assertLess(five_ae4["stimulated_endpoint_t200"]["max_abs_raw_rhs"], 7.3e-5)
        self.assertAlmostEqual(
            five_ae4["stimulated_endpoint_t200"]["max_abs_scaled_rhs"],
            0.07222091282319474,
            delta=2e-10,
        )

    def test_checked_in_machine_outputs_regenerate_byte_for_byte(self) -> None:
        for path in self.result_paths:
            self.assertEqual(path.read_bytes(), self.before[path], path.name)
        loaded = json.loads(self.result_paths[0].read_text(encoding="utf-8"))
        self.assertEqual(loaded, self.summary)


if __name__ == "__main__":
    unittest.main()
