"""Independent native-source contract tests for Task 13B.

The tests contain no held-out phenotype target and do not import the
production native-source builder.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import inspect
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model.independent import IndependentConstantStimulus, IndependentWholeCell
from modern_full_model.independent_native import (
    AN_ABS_AE4_NHE,
    FROZEN_AE4_CARRIER_FMOL,
    FROZEN_NHE1_CAPACITY_FMOL_S,
    FROZEN_NKCC1_CAPACITY_FMOL_S,
    IndependentNativeSourceDefinition,
    NATIVE_CA_HALF_UM,
    NATIVE_CA_HILL,
    NATIVE_G_CL_APICAL_S,
    NATIVE_G_K_TOTAL_S,
    N_ABS_NKCC,
    STRICT_NATIVE,
    independent_hydraulic_scaled_water_parameters,
    independent_native_wt_gate_failures,
    independent_native_wt_observables,
    refine_independent_native_root,
    state_from_independent_rest_coordinates,
)
from modern_full_model.run_independent_native import (
    CLUSTER_RELATIVE_TOLERANCE,
    PRODUCTION_ELIGIBILITY,
    ROOT_BOUNDS_LOWER,
    ROOT_BOUNDS_SPAN,
    _atomic_write_text,
    _canonical_panel_definitions,
    _connected_root_components,
    _csv_snapshot,
    _write_csv,
    audit_native_hierarchy,
)


def definition(
    source_class: str,
    source_scale: float = 1.0,
    hydraulic: float = 1.0,
    hydraulic_mode: str | None = None,
):
    return IndependentNativeSourceDefinition(
        source_class=source_class,
        source_scale=source_scale,
        pump_capacity_scale=1.0,
        apical_pump_fraction=0.075075,
        apical_k_fraction=0.30,
        ae4_cation_fraction=0.10,
        hydraulic_scale=hydraulic,
        cell_other_impermeant_osmoles_fmol=124.0,
        hydraulic_mode=(
            hydraulic_mode or ("H1" if hydraulic == 1.0 else "HW")
        ),
    )


def _synthetic_panel_row(panel_id, expected):
    if expected["kind"] == "strict":
        stage = (
            "STRICT_H1_PUMP1"
            if expected["pump_capacity_scale"] == 1.0
            else "CONDITIONAL_PUMP_AFTER_ZERO_STRICT_ROOTS"
        )
    else:
        stage = {
            "hw": "DIAGNOSTIC_HYDRAULIC_REROOT",
            "hwq": "DIAGNOSTIC_HWQ_REROOT",
            "absolute": "PRODUCTION_ABSOLUTE_SOURCE_GRID",
            "stress": "DIAGNOSTIC_ABSOLUTE_SOURCE_STRESS",
        }[expected["kind"]]
    start_ids = ("frozen_reference", "lower_quartile", "upper_quartile")
    return {
        "panel_id": panel_id,
        **{key: value for key, value in expected.items() if key != "kind"},
        "stage": stage,
        "attempt_count": 3,
        "converged_attempt_count": 0,
        "numerical_root_count": 0,
        "wt_passing_root_count": 0,
        "screen_start_count": 3,
        "confirmation_start_count": 0,
        "geometry_start_count": 0,
        "start_design_json": json.dumps(
            {
                "kind": "one_parent_or_frozen_anchor_plus_both_quartiles",
                "count": 3,
                "start_ids": start_ids,
                "actual_persisted_attempt_set": True,
                "both_quartiles_present": True,
                "legacy_priority_truncation": False,
            }
        ),
    }


def _synthetic_attempt(panel, start_id, vector=None, numerical=False):
    endpoint = (
        np.asarray(ROOT_BOUNDS_LOWER + 0.5 * ROOT_BOUNDS_SPAN)
        if vector is None
        else np.asarray(vector, dtype=float)
    )
    return {
        "attempt_id": f"{panel['panel_id']}|{start_id}",
        "panel_id": panel["panel_id"],
        "source_class": panel["source_class"],
        "source_scale": panel["source_scale"],
        "pump_capacity_scale": panel["pump_capacity_scale"],
        "hydraulic_scale": panel["hydraulic_scale"],
        "hydraulic_mode": panel["hydraulic_mode"],
        "eligibility": panel["eligibility"],
        "start_id": start_id,
        "admissible_state": "True",
        "converged_root": "True" if numerical else "False",
        "passes_numerical_gate": "True" if numerical else "False",
        "max_abs_scaled_residual": "0" if numerical else "1",
        "coordinates_and_other_json": json.dumps(endpoint.tolist()),
    }


def _complete_zero_pass_hierarchy():
    definitions = _canonical_panel_definitions()
    panels = [
        _synthetic_panel_row(panel_id, expected)
        for panel_id, expected in definitions.items()
        if expected["kind"] in {"strict", "hw", "hwq", "absolute"}
    ]
    attempts = [
        _synthetic_attempt(panel, start_id)
        for panel in panels
        for start_id in ("frozen_reference", "lower_quartile", "upper_quartile")
    ]
    return panels, attempts


class TestIndependentNativeSourceContract(unittest.TestCase):
    def test_source_is_independent_of_production_native_builder_and_holdout(self) -> None:
        import modern_full_model.independent_native as module

        source = inspect.getsource(module)
        self.assertNotIn("native_source_panel", source.replace(":mod:`native_source_panel`", ""))
        self.assertNotIn("heldout_targets", source)

    def test_literal_native_conductance_and_gate_map(self) -> None:
        parameters, ae4 = definition(STRICT_NATIVE).reconstruct()
        membranes = parameters.membranes
        self.assertEqual(membranes.g_cl_apical_S, 31.4e-9)
        self.assertEqual(membranes.g_k_total_S, 14.0e-9)
        self.assertEqual(membranes.calcium_half_uM, 0.26)
        self.assertEqual(membranes.calcium_hill, 1.46)
        self.assertEqual(membranes.g_apical_background_S, 0.0)
        self.assertEqual(membranes.g_basolateral_background_S, 0.0)
        self.assertEqual(ae4.carrier_amount_fmol, FROZEN_AE4_CARRIER_FMOL)
        self.assertEqual(
            parameters.homeostasis.nhe1_capacity_fmol_s,
            FROZEN_NHE1_CAPACITY_FMOL_S,
        )
        self.assertEqual(
            parameters.homeostasis.nkcc1_capacity_fmol_s,
            FROZEN_NKCC1_CAPACITY_FMOL_S,
        )
        self.assertEqual(NATIVE_G_CL_APICAL_S, membranes.g_cl_apical_S)
        self.assertEqual(NATIVE_G_K_TOTAL_S, membranes.g_k_total_S)
        self.assertEqual(NATIVE_CA_HALF_UM, membranes.calcium_half_uM)
        self.assertEqual(NATIVE_CA_HILL, membranes.calcium_hill)

    def test_source_scale_one_is_exact_nested_identity(self) -> None:
        strict_parameters, strict_ae4 = definition(STRICT_NATIVE).reconstruct()
        for source_class in (N_ABS_NKCC, AN_ABS_AE4_NHE):
            with self.subTest(source_class=source_class):
                parameters, ae4 = definition(source_class, 1.0).reconstruct()
                self.assertEqual(parameters, strict_parameters)
                self.assertEqual(ae4, strict_ae4)

    def test_n_abs_changes_only_whole_nkcc_capacity(self) -> None:
        base_parameters, base_ae4 = definition(STRICT_NATIVE).reconstruct()
        parameters, ae4 = definition(N_ABS_NKCC, 4.0).reconstruct()
        expected = asdict(base_parameters)
        expected["homeostasis"]["nkcc1_capacity_fmol_s"] *= 4.0
        self.assertEqual(asdict(parameters), expected)
        self.assertEqual(ae4, base_ae4)

    def test_an_abs_changes_only_ae4_and_nhe_common_capacity(self) -> None:
        base_parameters, base_ae4 = definition(STRICT_NATIVE).reconstruct()
        parameters, ae4 = definition(AN_ABS_AE4_NHE, 4.0).reconstruct()
        expected = asdict(base_parameters)
        expected["homeostasis"]["nhe1_capacity_fmol_s"] *= 4.0
        self.assertEqual(asdict(parameters), expected)
        self.assertEqual(ae4.carrier_amount_fmol, base_ae4.carrier_amount_fmol * 4.0)
        self.assertEqual(ae4.common_cl_attempt_rate_s, base_ae4.common_cl_attempt_rate_s)
        self.assertEqual(ae4.na_loaded_attempt_rate_s, base_ae4.na_loaded_attempt_rate_s)
        self.assertEqual(ae4.k_loaded_attempt_rate_s, base_ae4.k_loaded_attempt_rate_s)

    def test_hydraulic_diagnostic_scales_only_three_water_permeabilities(self) -> None:
        base_parameters, base_ae4 = definition(STRICT_NATIVE).reconstruct()
        parameters, ae4 = definition(STRICT_NATIVE, hydraulic=3.5).reconstruct()
        expected = asdict(base_parameters)
        for name in (
            "apical_hydraulic_pL_s_mOsm",
            "basolateral_hydraulic_pL_s_mOsm",
            "paracellular_hydraulic_pL_s_mOsm",
        ):
            expected["water"][name] *= 3.5
        self.assertEqual(asdict(parameters), expected)
        self.assertEqual(ae4, base_ae4)

    def test_hwq_scales_outflow_in_addition_to_three_permeabilities(self) -> None:
        base_parameters, base_ae4 = definition(STRICT_NATIVE).reconstruct()
        parameters, ae4 = definition(
            STRICT_NATIVE, hydraulic=3.5, hydraulic_mode="HWQ"
        ).reconstruct()
        expected = asdict(base_parameters)
        for name in (
            "apical_hydraulic_pL_s_mOsm",
            "basolateral_hydraulic_pL_s_mOsm",
            "paracellular_hydraulic_pL_s_mOsm",
            "outflow_rate_s",
        ):
            expected["water"][name] *= 3.5
        self.assertEqual(asdict(parameters), expected)
        self.assertEqual(ae4, base_ae4)

    def test_hydraulic_modes_are_exact_identity_at_scale_one(self) -> None:
        base_parameters, _base_ae4 = definition(STRICT_NATIVE).reconstruct()
        for mode in ("HW", "HWQ"):
            with self.subTest(mode=mode):
                transformed = independent_hydraulic_scaled_water_parameters(
                    base_parameters.water, scale=1.0, mode=mode
                )
                self.assertEqual(transformed, base_parameters.water)

    def test_public_member_enforces_hydraulic_mode_scale_invariant(self) -> None:
        with self.assertRaisesRegex(ValueError, "unit independent"):
            definition(STRICT_NATIVE, hydraulic=1.0, hydraulic_mode="HW")
        with self.assertRaisesRegex(ValueError, "nonunit independent"):
            definition(STRICT_NATIVE, hydraulic=2.0, hydraulic_mode="H1")

    def test_root_row_parser_requires_explicit_fitted_other_osmoles(self) -> None:
        row = {
            "source_class": N_ABS_NKCC,
            "source_scale": "2",
            "pump_capacity_scale": "1",
            "apical_pump_fraction": "0.075075",
            "apical_k_fraction": "0.30",
            "ae4_cation_fraction": "0.10",
            "hydraulic_scale": "1",
            "hydraulic_mode": "H1",
            "calibration_json": '{"cell_other_impermeant_osmoles_fmol": 123.5}',
        }
        parsed = IndependentNativeSourceDefinition.from_root_row(row)
        self.assertEqual(parsed.source_class, N_ABS_NKCC)
        self.assertEqual(parsed.source_scale, 2.0)
        self.assertEqual(parsed.cell_other_impermeant_osmoles_fmol, 123.5)
        with self.assertRaisesRegex(ValueError, "omits fitted OTHER"):
            IndependentNativeSourceDefinition.from_root_row(
                {**row, "calibration_json": "{}"}
            )

    def test_strict_native_rejects_a_hidden_nonunit_source_scale(self) -> None:
        with self.assertRaisesRegex(ValueError, "exact source scale one"):
            definition(STRICT_NATIVE, 2.0)

    def test_local_wt_gate_is_derived_from_independent_state_and_voltage(self) -> None:
        state = np.asarray(
            (
                20.0 * 1.3,
                120.0 * 1.3,
                50.0 * 1.3,
                5.0 * 1.3,
                10.0,
                1.3,
                140.0 * 0.1,
                5.0 * 0.1,
                130.0 * 0.1,
                6.0 * 0.1,
                0.5,
                0.1,
            )
        )
        evaluation = SimpleNamespace(
            cell_speciation=SimpleNamespace(ph=6.90),
            lumen_speciation=SimpleNamespace(ph=7.10),
            voltages_V={"apical": -0.050, "basolateral": -0.060},
        )
        observables = independent_native_wt_observables(evaluation, state)
        self.assertEqual(independent_native_wt_gate_failures(observables), ())
        self.assertEqual(observables["cl_i_mM"], 50.0)
        failures = independent_native_wt_gate_failures(
            {**observables, "cl_i_mM": 60.0}
        )
        self.assertTrue(any(item.startswith("cl_i_mM:") for item in failures))

    def test_csv_snapshot_hashes_the_only_bytes_it_parses(self) -> None:
        payload = b"panel_id,value\nP0,1\n"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "panel.csv"
            with patch.object(Path, "read_bytes", return_value=payload) as reader:
                rows, digest = _csv_snapshot(path)
        self.assertEqual(reader.call_count, 1)
        self.assertEqual(rows, [{"panel_id": "P0", "value": "1"}])
        self.assertEqual(digest, hashlib.sha256(payload).hexdigest())

    def test_independent_outputs_publish_by_atomic_sibling_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            csv_path = directory / "independent.csv"
            summary_path = directory / "summary.json"
            _write_csv(csv_path, [{"root_id": "R0", "passes": True}])
            _atomic_write_text(summary_path, '{"status":"PASS"}\n')
            self.assertEqual(
                csv_path.read_text(encoding="utf-8"),
                "root_id,passes\nR0,True\n",
            )
            self.assertEqual(
                summary_path.read_text(encoding="utf-8"),
                '{"status":"PASS"}\n',
            )
            self.assertEqual(
                sorted(path.name for path in directory.iterdir()),
                ["independent.csv", "summary.json"],
            )

    def test_complete_zero_pass_hierarchy_has_exact_current_three_start_sets(self) -> None:
        panels, attempts = _complete_zero_pass_hierarchy()
        audit = audit_native_hierarchy(panels, (), attempts)
        self.assertEqual(audit["status"], "PASS", audit["issues"][:5])
        self.assertEqual(audit["licensed_pump_contexts"], [0.5, 1.0, 2.0])
        self.assertEqual(audit["required_panel_count"], 1125)

    def test_stale_deep_marker_and_duplicate_attempt_fail_hierarchy(self) -> None:
        panels, attempts = _complete_zero_pass_hierarchy()
        panels[0]["geometry_start_count"] = 33
        attempts.append(dict(attempts[0]))
        audit = audit_native_hierarchy(panels, (), attempts)
        issue_names = {item["issue"] for item in audit["issues"]}
        self.assertEqual(audit["status"], "FAIL")
        self.assertIn("stale_or_inexact_3_7_33_lineage", issue_names)
        self.assertIn("duplicate_attempt_id", issue_names)

    def test_later_pump1_pass_excludes_stale_conditional_contexts(self) -> None:
        panels, attempts = _complete_zero_pass_hierarchy()
        pump1 = next(
            panel
            for panel in panels
            if panel["source_class"] == STRICT_NATIVE
            and panel["pump_capacity_scale"] == 1.0
        )
        pump1["wt_passing_root_count"] = 1
        audit = audit_native_hierarchy(panels, (), attempts)
        self.assertEqual(audit["licensed_pump_contexts"], [1.0])
        self.assertTrue(
            any(
                item["issue"] == "unlicensed_pump_context_panel"
                for item in audit["issues"]
            )
        )
        by_id = {panel["panel_id"]: panel for panel in panels}
        self.assertTrue(
            all(
                float(by_id[panel_id]["pump_capacity_scale"]) == 1.0
                for panel_id in audit["eligible_production_panel_ids"]
            )
        )

    def test_legacy_displaced_quartile_is_diagnostic_only_and_visible(self) -> None:
        panel_id, expected = next(
            (panel_id, expected)
            for panel_id, expected in _canonical_panel_definitions().items()
            if expected["kind"] == "hw"
            and expected["pump_capacity_scale"] == 1.0
        )
        panel = _synthetic_panel_row(panel_id, expected)
        legacy_ids = (
            "frozen_reference",
            "continuation_00",
            "lower_quartile",
        )
        panel["start_design_json"] = json.dumps(
            {
                "kind": "LEGACY_PRIORITY_TRUNCATION_UPPER_QUARTILE_DISPLACED",
                "count": 3,
                "start_ids": legacy_ids,
                "actual_persisted_attempt_set": True,
                "both_quartiles_present": False,
                "legacy_priority_truncation": True,
            }
        )
        attempts = tuple(
            _synthetic_attempt(panel, start_id) for start_id in legacy_ids
        )
        audit = audit_native_hierarchy((panel,), (), attempts)
        self.assertEqual(
            audit["accepted_legacy_diagnostic_start_design_count"], 1
        )
        self.assertEqual(
            audit["accepted_legacy_diagnostic_start_design_panel_ids"],
            [panel_id],
        )
        panel_issues = {
            item["issue"]
            for item in audit["issues"]
            if item["id"] == panel_id
        }
        self.assertNotIn("legacy_start_design_not_eligible", panel_issues)
        self.assertNotIn("canonical_quartile_starts_missing", panel_issues)

        production_id, production_expected = next(
            (candidate_id, candidate)
            for candidate_id, candidate in _canonical_panel_definitions().items()
            if candidate["kind"] == "absolute"
            and candidate["pump_capacity_scale"] == 1.0
        )
        production = _synthetic_panel_row(production_id, production_expected)
        production["start_design_json"] = panel["start_design_json"]
        production_attempts = tuple(
            _synthetic_attempt(production, start_id) for start_id in legacy_ids
        )
        rejected = audit_native_hierarchy(
            (production,), (), production_attempts
        )
        rejected_issues = {
            item["issue"]
            for item in rejected["issues"]
            if item["id"] == production_id
        }
        self.assertEqual(
            rejected["accepted_legacy_diagnostic_start_design_count"], 0
        )
        self.assertIn("legacy_start_design_not_eligible", rejected_issues)
        self.assertIn("canonical_quartile_starts_missing", rejected_issues)

    def test_cluster_representatives_must_remain_mutually_separated(self) -> None:
        panel_id, expected = next(
            (panel_id, expected)
            for panel_id, expected in _canonical_panel_definitions().items()
            if expected["kind"] == "strict"
            and expected["pump_capacity_scale"] == 1.0
        )
        panel = _synthetic_panel_row(panel_id, expected)
        panel["converged_attempt_count"] = 3
        panel["numerical_root_count"] = 2
        base = ROOT_BOUNDS_LOWER + 0.5 * ROOT_BOUNDS_SPAN
        endpoint_a = base.copy()
        endpoint_b = base.copy()
        endpoint_c = base.copy()
        endpoint_b[0] += 0.9 * CLUSTER_RELATIVE_TOLERANCE * ROOT_BOUNDS_SPAN[0]
        endpoint_c[0] += 1.8 * CLUSTER_RELATIVE_TOLERANCE * ROOT_BOUNDS_SPAN[0]
        self.assertEqual(
            _connected_root_components(
                {
                    "a": endpoint_a,
                    "b": endpoint_b,
                    "c": endpoint_c,
                }
            ),
            (("a", "b", "c"),),
        )
        attempts = [
            _synthetic_attempt(panel, "frozen_reference", endpoint_a, numerical=True),
            _synthetic_attempt(panel, "lower_quartile", endpoint_b, numerical=True),
            _synthetic_attempt(panel, "upper_quartile", endpoint_c, numerical=True),
        ]

        def root(root_id, representative, members):
            return {
                "root_id": root_id,
                "panel_id": panel_id,
                "source_class": panel["source_class"],
                "source_scale": panel["source_scale"],
                "pump_capacity_scale": panel["pump_capacity_scale"],
                "hydraulic_scale": panel["hydraulic_scale"],
                "hydraulic_mode": panel["hydraulic_mode"],
                "eligibility": PRODUCTION_ELIGIBILITY,
                "passes_numerical_gate": "True",
                "passes_wt_gate": "False",
                "root_search_stage": panel["stage"],
                "root_search_start_count": 3,
                "panel_confirmation_start_count": 0,
                "panel_geometry_start_count": 0,
                "member_start_ids_json": json.dumps(members),
                "coordinates_json": json.dumps(representative[:10].tolist()),
                "calibration_json": json.dumps(
                    {
                        "cell_other_impermeant_osmoles_fmol": float(
                            representative[-1]
                        )
                    }
                ),
            }

        roots = (
            root(
                f"{panel_id}:B00",
                endpoint_b,
                ("frozen_reference", "lower_quartile"),
            ),
            root(f"{panel_id}:B01", endpoint_c, ("upper_quartile",)),
        )
        audit = audit_native_hierarchy((panel,), roots, attempts)
        self.assertTrue(
            any(
                item["issue"]
                == "stored_clusters_do_not_match_deterministic_connected_components"
                for item in audit["issues"]
            )
        )
        self.assertTrue(
            any(
                item["issue"].startswith(
                    "stored_root_representatives_not_separated:"
                )
                for item in audit["issues"]
            )
        )

    def test_n_abs_wt_candidate_jointly_reroots_with_other_and_rank_eleven(self) -> None:
        candidate = np.asarray(
            (
                24.375004480505822,
                138.68089788712203,
                66.94776908416362,
                6.830518062268592,
                18.04637107751941,
                1.3,
                13.534346011579268,
                2.243735106264099,
                15.265233459193988,
                0.589897675936328,
                0.5128476586493772,
                0.10495597693809071,
            ),
            dtype=float,
        )
        declared = IndependentNativeSourceDefinition(
            source_class=N_ABS_NKCC,
            source_scale=4.0,
            pump_capacity_scale=1.0,
            apical_pump_fraction=0.075075,
            apical_k_fraction=0.30,
            ae4_cation_fraction=0.05,
            hydraulic_scale=1.0,
            cell_other_impermeant_osmoles_fmol=124.69936937127669,
        )
        result = refine_independent_native_root(
            declared,
            candidate,
            max_nfev=100,
        )
        self.assertTrue(result.success, result.message)
        self.assertLess(result.max_abs_scaled_residual, 1.0e-7)
        self.assertEqual(result.jacobian_rank, 11)
        self.assertEqual(result.jacobian_nullity, 0)
        np.testing.assert_allclose(result.state, candidate, rtol=0.0, atol=1.0e-10)

    def test_physiological_coordinate_map_recovers_known_native_root(self) -> None:
        declared = IndependentNativeSourceDefinition(
            source_class=N_ABS_NKCC,
            source_scale=4.0,
            pump_capacity_scale=1.0,
            apical_pump_fraction=0.075075,
            apical_k_fraction=0.30,
            ae4_cation_fraction=0.05,
            hydraulic_scale=1.0,
            cell_other_impermeant_osmoles_fmol=124.69936937127669,
        )
        parameters, _ae4 = declared.reconstruct()
        coordinates = np.asarray(
            (
                18.75000344654294,
                106.67761375932463,
                5.254244663283532,
                6.839866620908583,
                1.3,
                128.95259904600414,
                21.37786881434668,
                5.620429566238851,
                6.921669926083349,
                0.10495597693809071,
            )
        )
        expected = np.asarray(
            (
                24.375004480505822,
                138.68089788712203,
                66.94776908416362,
                6.830518062268592,
                18.04637107751941,
                1.3,
                13.534346011579268,
                2.243735106264099,
                15.265233459193988,
                0.589897675936328,
                0.5128476586493772,
                0.10495597693809071,
            )
        )
        actual = state_from_independent_rest_coordinates(parameters, coordinates)
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2.0e-11)


class TestProductionNativeBuilderAgainstIndependentTranscription(unittest.TestCase):
    def test_representative_source_builds_and_rhs_match_second_transcription(self) -> None:
        # Production objects are imported only here as comparison targets.  The
        # independent builder and equation implementation above do not depend
        # on either production module.
        from modern_full_model.native_source_panel import (
            NativeSourcePanelMember,
            build_native_source_model,
        )

        cases = (
            (STRICT_NATIVE, 1.0, 1.0, 1.0, "H1", "PRODUCTION_PRE_REVEAL"),
            (N_ABS_NKCC, 4.0, 1.0, 1.0, "H1", "PRODUCTION_PRE_REVEAL"),
            (AN_ABS_AE4_NHE, 4.0, 1.0, 1.0, "H1", "PRODUCTION_PRE_REVEAL"),
            (STRICT_NATIVE, 1.0, 2.0, 3.5, "HW", "DIAGNOSTIC_ONLY_HW_HYDRAULIC_TRANSFER"),
            (STRICT_NATIVE, 1.0, 2.0, 3.5, "HWQ", "DIAGNOSTIC_ONLY_HWQ_HYDRAULIC_TRANSFER"),
        )
        for source_class, source_scale, pump_scale, hydraulic, mode, eligibility in cases:
            with self.subTest(
                source_class=source_class,
                source_scale=source_scale,
                pump_scale=pump_scale,
                hydraulic=hydraulic,
                mode=mode,
            ):
                member = NativeSourcePanelMember(
                    panel_id="INDEPENDENT_COMPARISON",
                    topology_id="PNOM_KMID",
                    routing_id="AE4NA10",
                    source_class=source_class,
                    source_scale=source_scale,
                    pump_capacity_scale=pump_scale,
                    apical_pump_fraction=0.075075,
                    apical_k_fraction=0.30,
                    ae4_cation_fraction=0.10,
                    hydraulic_scale=hydraulic,
                    hydraulic_mode=mode,
                    eligibility=eligibility,
                )
                independent_definition = IndependentNativeSourceDefinition(
                    source_class=source_class,
                    source_scale=source_scale,
                    pump_capacity_scale=pump_scale,
                    apical_pump_fraction=0.075075,
                    apical_k_fraction=0.30,
                    ae4_cation_fraction=0.10,
                    hydraulic_scale=hydraulic,
                    cell_other_impermeant_osmoles_fmol=124.0,
                    hydraulic_mode=mode,
                )
                production = build_native_source_model(
                    member,
                    other_impermeant_osmoles_fmol=124.0,
                )
                parameters, ae4 = independent_definition.reconstruct()
                self.assertEqual(production.parameters, parameters)
                production_ae4 = production.ae4_parameters
                self.assertEqual(
                    production_ae4.carrier_amount_fmol,
                    ae4.carrier_amount_fmol,
                )
                self.assertEqual(
                    production_ae4.na_loaded_attempt_rate_s,
                    ae4.na_loaded_attempt_rate_s,
                )
                self.assertEqual(
                    production_ae4.k_loaded_attempt_rate_s,
                    ae4.k_loaded_attempt_rate_s,
                )

                independent = IndependentWholeCell(
                    parameters,
                    ae4_parameters=ae4,
                    stimulus=IndependentConstantStimulus(0.058, 0.0),
                )
                state = production.initial_state()
                expected = production.evaluate(0.0, state)
                actual = independent.evaluate(0.0, state)
                np.testing.assert_allclose(
                    actual.rhs,
                    expected.rhs,
                    rtol=5.0e-13,
                    atol=5.0e-14,
                )
                self.assertLess(
                    max(abs(value) for value in actual.current_residuals_A.values()),
                    1.0e-20,
                )


if __name__ == "__main__":
    unittest.main()
