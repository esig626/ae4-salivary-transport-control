"""Focused WT-only tests for the native source-map reconstruction panel."""

from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from modern_full_model.membranes import hill_activation
from modern_full_model.native_source_panel import (
    AN_ABS_AE4_NHE,
    FROZEN_AE4_CARRIER_FMOL,
    FROZEN_NHE1_CAPACITY_FMOL_S,
    FROZEN_NKCC1_CAPACITY_FMOL_S,
    G_CACC_MAX_S,
    G_K_TOTAL_MAX_S,
    N_ABS_NKCC,
    PALK_CA_HALF_UM,
    PALK_CA_HILL,
    STRICT_NATIVE,
    build_native_source_model,
    declared_absolute_source_panel,
    declared_hydraulic_diagnostic_panel,
    declared_native_panel,
    hydraulic_scaled_water_parameters,
    native_root_residual,
    native_source_map,
    _named_starts,
    _root_components,
)
from modern_full_model.run_native_source_panel import (
    ATTEMPT_FIELDS,
    CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY,
    LEGACY_HW_DIAGNOSTIC_ELIGIBILITY,
    LEGACY_HW_ELIGIBILITY_MIGRATION,
    PANEL_FIELDS,
    ROOT_FIELDS,
    _bool,
    _canonical_hydraulic_mode,
    _canonicalize_current_generation,
    _commit_ledger_status,
    _confirmation_panel_rows,
    _current_production_lineage_status,
    _eligible_production,
    _panel_seed,
    _prune_unlicensed_pump_contexts,
    _read_rows,
    _replace_attempt_rows_atomically,
    _replace_root_rows_atomically,
    _require_complete_upstream_grids,
    _upstream_grid_status,
    _write_rows,
    execute_stage,
    persist_reports,
    run_panel,
    scale_one_identity_rows,
    write_hashes,
    write_summary,
)


def _member(source_class: str = STRICT_NATIVE, source_scale: float = 1.0):
    return next(
        item
        for item in declared_native_panel(
            source_class=source_class, source_scale=source_scale
        )
        if item.topology_id == "PNOM_KMID" and item.ae4_cation_fraction == 0.10
    )


def _screened_panel_row(member):
    return {
        **asdict(member),
        "stage": "SYNTHETIC_HIERARCHY_TEST",
        "attempt_count": 3,
        "converged_attempt_count": 0,
        "numerical_root_count": 0,
        "wt_passing_root_count": 0,
        "best_max_abs_scaled_residual": 1.0,
        "best_attempt_id": "lower_quartile",
        "failure_union_json": "[]",
        "start_design_json": "{}",
        "screen_start_count": 3,
        "confirmation_start_count": 0,
        "geometry_start_count": 0,
    }


def _complete_upstream_rows():
    pumps = (0.5, 1.0, 2.0)
    members = [
        member
        for pump in pumps
        for member in declared_native_panel(pump_capacity_scale=pump)
    ]
    members.extend(
        declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=pumps,
            hydraulic_modes=("HW",),
        )
    )
    members.extend(
        declared_hydraulic_diagnostic_panel(
            pump_capacity_scales=pumps,
            hydraulic_modes=("HWQ",),
        )
    )
    members.extend(
        declared_absolute_source_panel(
            include_identity=False,
            pump_capacity_scales=pumps,
        )
    )
    return [_screened_panel_row(member) for member in members]


def _synthetic_attempt(start_id: str):
    return SimpleNamespace(
        start_id=start_id,
        optimizer_success=False,
        admissible_state=True,
        converged_root=False,
        passes_numerical_gate=False,
        passes_wt_gate=False,
        max_abs_scaled_residual=1.0,
        max_abs_amount_rhs_fmol_s=1.0,
        max_abs_volume_rhs_pL_s=1.0,
        max_abs_omitted_charge_rhs_fmol_equivalent_s=1.0,
        normalized_jacobian_rank=0,
        normalized_jacobian_nullity=11,
        boundary_hits=(),
        wt_gate_failures=("synthetic",),
        coordinates_and_other=tuple(float(index + 1) for index in range(11)),
        cost=1.0,
        optimality=1.0,
        nfev=1,
        message="synthetic",
    )


def _synthetic_report(member, count: int):
    standard = ["frozen_reference", "lower_quartile", "upper_quartile"]
    start_ids = standard + [f"uniform_{index:03d}" for index in range(3, count)]
    attempts = tuple(_synthetic_attempt(start_id) for start_id in start_ids)
    return SimpleNamespace(
        member=member,
        attempts=attempts,
        roots=(),
        start_design={"count": count, "start_ids": start_ids},
    )


class TestNativeSourceMap(unittest.TestCase):
    def test_literal_conductances_and_palk_gate_replace_common_scale(self) -> None:
        model = build_native_source_model(
            _member(), other_impermeant_osmoles_fmol=125.0
        )
        membrane = model.parameters.membranes
        self.assertEqual(membrane.g_cl_apical_S, G_CACC_MAX_S)
        self.assertEqual(membrane.g_k_total_S, G_K_TOTAL_MAX_S)
        self.assertEqual(membrane.calcium_half_uM, PALK_CA_HALF_UM)
        self.assertEqual(membrane.calcium_hill, PALK_CA_HILL)
        self.assertEqual(membrane.g_apical_background_S, 0.0)
        self.assertEqual(membrane.g_basolateral_background_S, 0.0)
        self.assertAlmostEqual(
            hill_activation(0.058, PALK_CA_HALF_UM, PALK_CA_HILL),
            0.10062055983703325,
        )

    def test_declared_panel_is_five_topologies_by_three_ae4_routings(self) -> None:
        panel = declared_native_panel()
        self.assertEqual(len(panel), 15)
        self.assertEqual(len({item.topology_id for item in panel}), 5)
        self.assertEqual(
            {item.ae4_cation_fraction for item in panel}, {0.05, 0.10, 0.20}
        )
        self.assertTrue(all(item.hydraulic_scale == 1.0 for item in panel))

    def test_unit_absolute_source_rows_are_exact_nested_identity(self) -> None:
        models = [
            build_native_source_model(
                _member(source_class), other_impermeant_osmoles_fmol=125.0
            )
            for source_class in (STRICT_NATIVE, N_ABS_NKCC, AN_ABS_AE4_NHE)
        ]
        reference = models[0]
        for candidate in models[1:]:
            self.assertEqual(asdict(candidate.parameters), asdict(reference.parameters))
            self.assertEqual(asdict(candidate.ae4_parameters), asdict(reference.ae4_parameters))
            state = reference.initial_state()
            np.testing.assert_array_equal(candidate.rhs(0.0, state), reference.rhs(0.0, state))

    def test_absolute_source_classes_change_only_the_declared_capacities(self) -> None:
        strict = build_native_source_model(
            _member(), other_impermeant_osmoles_fmol=125.0
        )
        n_abs = build_native_source_model(
            _member(N_ABS_NKCC, 4.0), other_impermeant_osmoles_fmol=125.0
        )
        an_abs = build_native_source_model(
            _member(AN_ABS_AE4_NHE, 4.0), other_impermeant_osmoles_fmol=125.0
        )
        self.assertEqual(
            n_abs.parameters.homeostasis.nkcc1_capacity_fmol_s,
            4.0 * FROZEN_NKCC1_CAPACITY_FMOL_S,
        )
        self.assertEqual(
            n_abs.parameters.homeostasis.nhe1_capacity_fmol_s,
            FROZEN_NHE1_CAPACITY_FMOL_S,
        )
        self.assertEqual(
            n_abs.ae4_parameters.carrier_amount_fmol,
            FROZEN_AE4_CARRIER_FMOL,
        )
        self.assertEqual(
            an_abs.parameters.homeostasis.nkcc1_capacity_fmol_s,
            FROZEN_NKCC1_CAPACITY_FMOL_S,
        )
        self.assertEqual(
            an_abs.parameters.homeostasis.nhe1_capacity_fmol_s,
            4.0 * FROZEN_NHE1_CAPACITY_FMOL_S,
        )
        self.assertEqual(
            an_abs.ae4_parameters.carrier_amount_fmol,
            4.0 * FROZEN_AE4_CARRIER_FMOL,
        )
        self.assertEqual(
            an_abs.ae4_parameters.mixed_bath_na_attempt_fraction,
            strict.ae4_parameters.mixed_bath_na_attempt_fraction,
        )

    def test_hydraulic_scaling_is_diagnostic_only_and_h1_is_not_duplicated(self) -> None:
        diagnostics = declared_hydraulic_diagnostic_panel()
        self.assertEqual(len(diagnostics), 150)
        self.assertNotIn(1.0, {item.hydraulic_scale for item in diagnostics})
        self.assertEqual({item.hydraulic_mode for item in diagnostics}, {"HW", "HWQ"})
        self.assertTrue(all("DIAGNOSTIC" in item.eligibility for item in diagnostics))

    def test_hw_and_hwq_share_h1_identity_but_only_hwq_scales_outflow(self) -> None:
        strict_member = _member()
        strict = build_native_source_model(
            strict_member, other_impermeant_osmoles_fmol=125.0
        )
        for mode in ("HW", "HWQ"):
            with self.assertRaises(ValueError):
                replace(
                    strict_member,
                    panel_id=f"{mode}_H1_INVALID_PUBLIC_MEMBER",
                    hydraulic_mode=mode,
                    eligibility=f"DIAGNOSTIC_ONLY_{mode}_IDENTITY",
                )
            identity_water = hydraulic_scaled_water_parameters(
                strict.parameters.water, scale=1.0, mode=mode
            )
            self.assertEqual(asdict(identity_water), asdict(strict.parameters.water))

        hw_member = next(
            item
            for item in declared_hydraulic_diagnostic_panel(
                hydraulic_modes=("HW",)
            )
            if item.hydraulic_scale == 3.5
            and item.topology_id == "PNOM_KMID"
            and item.ae4_cation_fraction == 0.10
        )
        hwq_member = next(
            item
            for item in declared_hydraulic_diagnostic_panel(
                hydraulic_modes=("HWQ",)
            )
            if item.hydraulic_scale == 3.5
            and item.topology_id == "PNOM_KMID"
            and item.ae4_cation_fraction == 0.10
        )
        hw = build_native_source_model(hw_member, other_impermeant_osmoles_fmol=125.0)
        hwq = build_native_source_model(hwq_member, other_impermeant_osmoles_fmol=125.0)
        self.assertEqual(
            hw.parameters.water.apical_hydraulic_pL_s_mOsm,
            hwq.parameters.water.apical_hydraulic_pL_s_mOsm,
        )
        self.assertEqual(
            hw.parameters.water.outflow_rate_s,
            strict.parameters.water.outflow_rate_s,
        )
        self.assertEqual(
            hwq.parameters.water.outflow_rate_s,
            3.5 * strict.parameters.water.outflow_rate_s,
        )

    def test_native_root_has_only_other_as_continuous_nuisance(self) -> None:
        source_map = native_source_map()
        self.assertEqual(
            source_map["root_fit"]["continuous_nuisance"],
            "OTHER intracellular impermeant osmoles only",
        )
        vector = np.asarray(
            (
                14.756773495778985,
                111.10374637122905,
                6.194122572604099,
                6.91,
                1.3,
                146.13568410951905,
                4.087361120509596,
                12.891629639536415,
                7.276453648930862,
                0.10164106693766255,
                124.71681310732373,
            )
        )
        self.assertEqual(native_root_residual(_member(), vector).shape, (11,))

    def test_absolute_panel_retains_all_structural_rows(self) -> None:
        panel = declared_absolute_source_panel(include_identity=False)
        self.assertEqual(len(panel), 210)
        self.assertEqual({item.source_class for item in panel}, {N_ABS_NKCC, AN_ABS_AE4_NHE})
        self.assertNotIn(1.0, {item.source_scale for item in panel})

    def test_production_confirmation_excludes_hw_and_hwq_diagnostics(self) -> None:
        base = {
            "wt_passing_root_count": "1",
            "hydraulic_scale": "1.0",
            "hydraulic_mode": "H1",
            "eligibility": "PRODUCTION_PRE_REVEAL",
        }
        rows = (
            {**base, "panel_id": "production"},
            {
                **base,
                "panel_id": "hw",
                "hydraulic_scale": "2.5",
                "hydraulic_mode": "HW",
                "eligibility": "DIAGNOSTIC_ONLY_HW_HYDRAULIC_TRANSFER",
            },
            {
                **base,
                "panel_id": "hwq",
                "hydraulic_scale": "2.5",
                "hydraulic_mode": "HWQ",
                "eligibility": "DIAGNOSTIC_ONLY_HWQ_HYDRAULIC_TRANSFER",
            },
            {**base, "panel_id": "failed", "wt_passing_root_count": "0"},
        )
        self.assertEqual(
            [row["panel_id"] for row in _confirmation_panel_rows(rows)],
            ["production"],
        )

    def test_rerun_replaces_root_clusters_panel_atomically(self) -> None:
        existing = (
            {"root_id": "P0:B00", "panel_id": "P0"},
            {"root_id": "P0:B01", "panel_id": "P0"},
            {"root_id": "P1:B00", "panel_id": "P1"},
        )
        incoming = ({"root_id": "P0:B00", "panel_id": "P0", "new": True},)
        merged = _replace_root_rows_atomically(
            existing, incoming, incoming_panel_ids={"P0"}
        )
        self.assertEqual({row["root_id"] for row in merged}, {"P0:B00", "P1:B00"})
        self.assertTrue(next(row for row in merged if row["root_id"] == "P0:B00")["new"])

    def test_rerun_replaces_attempts_panel_atomically(self) -> None:
        existing = (
            {"attempt_id": "P0|old", "panel_id": "P0"},
            {"attempt_id": "P1|keep", "panel_id": "P1"},
        )
        incoming = ({"attempt_id": "P0|new", "panel_id": "P0"},)
        merged = _replace_attempt_rows_atomically(
            existing, incoming, incoming_panel_ids={"P0"}
        )
        self.assertEqual(
            {row["attempt_id"] for row in merged}, {"P0|new", "P1|keep"}
        )

    def test_hydraulic_mode_backfill_is_total_for_legacy_rows(self) -> None:
        rows = (
            {"hydraulic_scale": "1.0", "hydraulic_mode": ""},
            {"hydraulic_scale": "2.5"},
            {"hydraulic_scale": "3.5", "hydraulic_mode": "HWQ"},
        )
        self.assertEqual(
            [_canonical_hydraulic_mode(row) for row in rows],
            ["H1", "HW", "HWQ"],
        )

    def test_three_and_seven_start_designs_preserve_both_quartiles(self) -> None:
        lower = np.zeros(11)
        preferred = np.full(11, 0.5)
        upper = np.ones(11)
        anchor = np.full(11, 0.6)
        three = _named_starts(
            lower, preferred, upper, count=3, seed=1, anchors=(anchor, anchor)
        )
        self.assertEqual(
            [name for name, _ in three],
            ["parent_continuation_00", "lower_quartile", "upper_quartile"],
        )
        seven = _named_starts(
            lower, preferred, upper, count=7, seed=1, anchors=(anchor, anchor)
        )
        names = [name for name, _ in seven]
        self.assertEqual(len(names), 7)
        self.assertIn("frozen_reference", names)
        self.assertIn("lower_quartile", names)
        self.assertIn("upper_quartile", names)
        self.assertEqual(sum(name.startswith("parent_continuation_") for name in names), 1)

    def test_multistart_contract_rejects_nonhierarchy_counts_and_duplicate_anchor(self) -> None:
        lower = np.zeros(11)
        preferred = np.full(11, 0.5)
        upper = np.ones(11)
        with self.assertRaisesRegex(ValueError, "exactly 3, 7, or 33"):
            _named_starts(lower, preferred, upper, count=4, seed=1)
        with self.assertRaisesRegex(ValueError, "distinct from both quartile"):
            _named_starts(
                lower,
                preferred,
                upper,
                count=3,
                seed=1,
                anchors=(np.full(11, 0.25),),
            )

    def test_root_components_are_order_independent_and_transitive(self) -> None:
        tolerance = 1.0
        points = (
            np.asarray((0.0,)),
            np.asarray((0.9,)),
            np.asarray((1.8,)),
            np.asarray((4.0,)),
        )
        components = _root_components(points, tolerance=tolerance)
        self.assertEqual({frozenset(item) for item in components}, {frozenset((0, 1, 2)), frozenset((3,))})
        reordered = tuple(points[index] for index in (3, 2, 0, 1))
        reordered_components = _root_components(reordered, tolerance=tolerance)
        coordinate_components = {
            frozenset(float(reordered[index][0]) for index in component)
            for component in reordered_components
        }
        self.assertEqual(
            coordinate_components,
            {frozenset((0.0, 0.9, 1.8)), frozenset((4.0,))},
        )

    def test_panel_seed_depends_on_identity_not_batch_order(self) -> None:
        members = declared_native_panel()[:3]
        forward = {member.panel_id: _panel_seed(42, member.panel_id) for member in members}
        reverse = {
            member.panel_id: _panel_seed(42, member.panel_id)
            for member in reversed(members)
        }
        self.assertEqual(forward, reverse)
        self.assertEqual(len(set(forward.values())), len(forward))

    def test_shallow_rerun_downgrades_current_geometry_provenance(self) -> None:
        member = _member(N_ABS_NKCC, 4.0)
        initial = _screened_panel_row(member)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_rows(output / "native_source_panel.csv", [initial], PANEL_FIELDS)
            _write_rows(output / "native_source_root_attempts.csv", [], ATTEMPT_FIELDS)
            _write_rows(output / "native_source_wt_roots.csv", [], ROOT_FIELDS)
            persist_reports(
                output,
                (_synthetic_report(member, 7),),
                stage="SEVEN_START_CONFIRMATION",
            )
            persist_reports(
                output,
                (_synthetic_report(member, 33),),
                stage="THIRTY_THREE_START_ROOT_GEOMETRY",
            )
            geometry = _read_rows(output / "native_source_panel.csv")[0]
            self.assertEqual(int(geometry["confirmation_start_count"]), 7)
            self.assertEqual(int(geometry["geometry_start_count"]), 33)
            self.assertEqual(
                {row["root_search_stage"] for row in _read_rows(output / "native_source_root_attempts.csv")},
                {"THIRTY_THREE_START_ROOT_GEOMETRY"},
            )

            persist_reports(
                output,
                (_synthetic_report(member, 7),),
                stage="SEVEN_START_CONFIRMATION",
            )
            downgraded = _read_rows(output / "native_source_panel.csv")[0]
            self.assertEqual(int(downgraded["confirmation_start_count"]), 7)
            self.assertEqual(int(downgraded["geometry_start_count"]), 0)
            self.assertEqual(int(downgraded["attempt_count"]), 7)
            self.assertEqual(
                {row["root_search_stage"] for row in _read_rows(output / "native_source_root_attempts.csv")},
                {"SEVEN_START_CONFIRMATION"},
            )

    def test_legacy_blank_current_screen_provenance_is_migrated_honestly(self) -> None:
        member = _member(N_ABS_NKCC, 4.0)
        panel = {
            **_screened_panel_row(member),
            "screen_start_count": 0,
            "confirmation_start_count": 7,
            "geometry_start_count": 33,
        }
        attempts = [
            {
                "attempt_id": f"{member.panel_id}|{start_id}",
                "panel_id": member.panel_id,
                "start_id": start_id,
                "root_search_stage": "",
                "root_search_start_count": "",
            }
            for start_id in ("frozen_reference", "lower_quartile", "upper_quartile")
        ]
        root = {
            "root_id": f"{member.panel_id}:B00",
            "panel_id": member.panel_id,
            "passes_wt_gate": "False",
            "root_search_stage": "",
            "root_search_start_count": "",
            "panel_confirmation_start_count": "",
            "panel_geometry_start_count": "",
        }
        panels, migrated_attempts, roots = _canonicalize_current_generation(
            [panel], attempts, [root]
        )
        self.assertEqual(panels[0]["screen_start_count"], 3)
        self.assertEqual(panels[0]["confirmation_start_count"], 0)
        self.assertEqual(panels[0]["geometry_start_count"], 0)
        self.assertEqual(
            {row["root_search_start_count"] for row in migrated_attempts}, {3}
        )
        self.assertEqual(roots[0]["root_search_start_count"], 3)
        self.assertEqual(roots[0]["panel_confirmation_start_count"], 0)
        self.assertEqual(roots[0]["panel_geometry_start_count"], 0)

    def test_canonicalize_migrates_only_exact_legacy_hw_diagnostics(self) -> None:
        member = declared_hydraulic_diagnostic_panel(
            hydraulic_modes=("HW",)
        )[0]
        legacy_panel = {
            **_screened_panel_row(member),
            "eligibility": LEGACY_HW_DIAGNOSTIC_ELIGIBILITY,
        }
        legacy_attempt = {
            "attempt_id": f"{member.panel_id}|frozen_reference",
            "panel_id": member.panel_id,
            "start_id": "frozen_reference",
            "eligibility": LEGACY_HW_DIAGNOSTIC_ELIGIBILITY,
        }
        legacy_root = {
            "root_id": f"{member.panel_id}:B00",
            "panel_id": member.panel_id,
            "passes_wt_gate": "False",
            "eligibility": LEGACY_HW_DIAGNOSTIC_ELIGIBILITY,
        }
        panels, attempts, roots = _canonicalize_current_generation(
            [legacy_panel], [legacy_attempt], [legacy_root]
        )
        self.assertEqual(
            panels[0]["eligibility"], CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
        )
        self.assertEqual(
            panels[0]["eligibility_migration"],
            LEGACY_HW_ELIGIBILITY_MIGRATION,
        )
        self.assertEqual(
            attempts[0]["eligibility"], CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
        )
        self.assertEqual(
            roots[0]["eligibility"], CANONICAL_HW_DIAGNOSTIC_ELIGIBILITY
        )

        malformed = {
            **legacy_panel,
            "panel_id": f"{member.panel_id}_NOT_DECLARED",
        }
        panels, _, _ = _canonicalize_current_generation([malformed], [], [])
        self.assertEqual(
            panels[0]["eligibility"], LEGACY_HW_DIAGNOSTIC_ELIGIBILITY
        )
        self.assertEqual(panels[0]["eligibility_migration"], "")

    def test_duplicate_member_ids_and_fail_open_eligibility_are_rejected(self) -> None:
        member = _member()
        with self.assertRaises(ValueError):
            run_panel(
                (member, member),
                start_count=3,
                seed=1,
                max_nfev=1,
                workers=1,
            )
        self.assertTrue(_eligible_production({"eligibility": "PRODUCTION_PRE_REVEAL"}))
        self.assertFalse(_eligible_production({"eligibility": "PRODUCTION_REJECTED"}))
        self.assertTrue(_bool("True"))
        self.assertFalse(_bool("False"))
        with self.assertRaisesRegex(ValueError, "canonical case-sensitive"):
            _bool("TRUE")

    def test_scale_one_alias_artifact_covers_all_licensed_contexts(self) -> None:
        repository = Path(__file__).resolve().parents[1]
        rows = scale_one_identity_rows(
            repository / "results" / "13B_modern_full_model"
        )
        self.assertEqual(len({row["alias_panel_id"] for row in rows}), 90)
        self.assertEqual({row["source_class"] for row in rows}, {N_ABS_NKCC, AN_ABS_AE4_NHE})
        self.assertTrue(all(row["source_scale"] == 1.0 for row in rows))
        self.assertTrue(all(row["exact_identity"] is True for row in rows))

    def test_upstream_hierarchy_requires_every_declared_panel_id(self) -> None:
        rows = _complete_upstream_rows()
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            path = output / "native_source_panel.csv"
            _write_rows(path, rows, PANEL_FIELDS)
            status = _require_complete_upstream_grids(output)
            self.assertTrue(status["complete"])
            self.assertEqual(
                status["expected"],
                {"strict": 45, "hw": 225, "hwq": 225, "absolute": 630},
            )

            missing_hwq = next(
                row["panel_id"]
                for row in rows
                if row["hydraulic_mode"] == "HWQ"
            )
            _write_rows(
                path,
                [row for row in rows if row["panel_id"] != missing_hwq],
                PANEL_FIELDS,
            )
            status = _upstream_grid_status(output)
            self.assertFalse(status["complete"])
            self.assertEqual(status["missing_panel_ids"]["hwq"], [missing_hwq])
            with self.assertRaisesRegex(RuntimeError, "hierarchy is incomplete"):
                _require_complete_upstream_grids(output)

            missing_absolute = next(
                row["panel_id"]
                for row in rows
                if row["source_class"] == N_ABS_NKCC
            )
            _write_rows(
                path,
                [row for row in rows if row["panel_id"] != missing_absolute],
                PANEL_FIELDS,
            )
            status = _upstream_grid_status(output)
            self.assertFalse(status["complete"])
            self.assertEqual(
                status["missing_panel_ids"]["absolute"], [missing_absolute]
            )

    def test_upstream_hierarchy_rejects_semantic_substitution_and_stale_pumps(self) -> None:
        rows = _complete_upstream_rows()
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            path = output / "native_source_panel.csv"

            altered = [dict(row) for row in rows]
            target = next(row for row in altered if row["source_class"] == N_ABS_NKCC)
            target["ae4_cation_fraction"] = 0.15
            _write_rows(path, altered, PANEL_FIELDS)
            status = _upstream_grid_status(output)
            self.assertFalse(status["complete"])
            self.assertIn(target["panel_id"], status["semantic_mismatches"])

            stale = [dict(row) for row in rows]
            pump_one = next(
                row
                for row in stale
                if row["source_class"] == STRICT_NATIVE
                and float(row["pump_capacity_scale"]) == 1.0
                and float(row["hydraulic_scale"]) == 1.0
            )
            pump_one["wt_passing_root_count"] = 1
            _write_rows(path, stale, PANEL_FIELDS)
            status = _upstream_grid_status(output)
            self.assertEqual(status["licensed_pump_contexts"], [1.0])
            self.assertFalse(status["complete"])
            self.assertTrue(status["unexpected_panel_ids"]["strict"])

    def test_stale_conditional_pump_rows_are_pruned_from_all_tables(self) -> None:
        rows = _complete_upstream_rows()
        pump_one = next(
            row
            for row in rows
            if row["source_class"] == STRICT_NATIVE
            and float(row["pump_capacity_scale"]) == 1.0
            and float(row["hydraulic_scale"]) == 1.0
        )
        pump_one["wt_passing_root_count"] = 1
        stale_panel = next(
            row for row in rows if float(row["pump_capacity_scale"]) == 0.5
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_rows(output / "native_source_panel.csv", rows, PANEL_FIELDS)
            _write_rows(
                output / "native_source_root_attempts.csv",
                [
                    {
                        "attempt_id": f"{stale_panel['panel_id']}|lower_quartile",
                        "panel_id": stale_panel["panel_id"],
                    }
                ],
                ATTEMPT_FIELDS,
            )
            _write_rows(
                output / "native_source_wt_roots.csv",
                [{"root_id": f"{stale_panel['panel_id']}:B00", "panel_id": stale_panel["panel_id"]}],
                ROOT_FIELDS,
            )
            removed = _prune_unlicensed_pump_contexts(output)
            self.assertIn(stale_panel["panel_id"], removed)
            self.assertEqual(
                {float(row["pump_capacity_scale"]) for row in _read_rows(output / "native_source_panel.csv")},
                {1.0},
            )
            self.assertEqual(_read_rows(output / "native_source_root_attempts.csv"), [])
            self.assertEqual(_read_rows(output / "native_source_wt_roots.csv"), [])

    def test_native_hash_ledger_is_a_fail_closed_commit_marker(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_rows(output / "native_source_panel.csv", [], PANEL_FIELDS)
            _write_rows(output / "native_source_root_attempts.csv", [], ATTEMPT_FIELDS)
            _write_rows(output / "native_source_wt_roots.csv", [], ROOT_FIELDS)
            write_hashes(output)
            self.assertTrue(_commit_ledger_status(output)["complete"])
            _write_rows(
                output / "native_source_panel.csv",
                [_screened_panel_row(_member())],
                PANEL_FIELDS,
            )
            status = _commit_ledger_status(output)
            self.assertFalse(status["complete"])
            self.assertIn(
                "commit_artifact_stale:native_source_panel.csv", status["issues"]
            )

    def test_current_production_lineage_requires_current_33_stage_rows(self) -> None:
        member = _member(N_ABS_NKCC, 4.0)
        panel = {
            **_screened_panel_row(member),
            "stage": "THIRTY_THREE_START_ROOT_GEOMETRY",
            "attempt_count": 33,
            "numerical_root_count": 1,
            "wt_passing_root_count": 1,
            "confirmation_start_count": 7,
            "geometry_start_count": 33,
        }
        start_ids = ["frozen_reference", "lower_quartile", "upper_quartile"] + [
            f"uniform_{index:03d}" for index in range(3, 33)
        ]
        attempts = [
            {
                "attempt_id": f"{member.panel_id}|{start_id}",
                "panel_id": member.panel_id,
                "start_id": start_id,
                "root_search_stage": "THIRTY_THREE_START_ROOT_GEOMETRY",
                "root_search_start_count": 33,
                "hydraulic_scale": 1.0,
                "hydraulic_mode": "H1",
                "eligibility": "PRODUCTION_PRE_REVEAL",
            }
            for start_id in start_ids
        ]
        root = {
            **asdict(member),
            "root_id": f"{member.panel_id}:B00",
            "passes_numerical_gate": "True",
            "passes_wt_gate": "True",
            "root_search_stage": "THIRTY_THREE_START_ROOT_GEOMETRY",
            "root_search_start_count": 33,
            "panel_confirmation_start_count": 7,
            "panel_geometry_start_count": 33,
            "member_start_ids_json": '["frozen_reference"]',
        }
        status = _current_production_lineage_status(
            [panel], attempts, [root], [root]
        )
        self.assertTrue(status["complete"], status["issues"])
        attempts[0]["root_search_stage"] = "SEVEN_START_CONFIRMATION"
        status = _current_production_lineage_status(
            [panel], attempts, [root], [root]
        )
        self.assertFalse(status["complete"])
        self.assertTrue(
            any("production_attempt_not_current_geometry" in item for item in status["issues"])
        )

    def test_readiness_rejects_committed_but_incomplete_tables(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_rows(output / "native_source_panel.csv", [], PANEL_FIELDS)
            _write_rows(output / "native_source_root_attempts.csv", [], ATTEMPT_FIELDS)
            _write_rows(output / "native_source_wt_roots.csv", [], ROOT_FIELDS)
            write_hashes(output)
            summary = write_summary(output, completed_stage="readiness")
            self.assertFalse(summary["exact_upstream_grids_complete"])
            self.assertFalse(summary["native_pre_reveal_freeze_ready"])
            _write_rows(
                output / "native_source_wt_roots.csv",
                [
                    {
                        "root_id": "malformed:B00",
                        "panel_id": "malformed",
                        "eligibility": "PRODUCTION_PRE_REVEAL",
                        "passes_numerical_gate": "TRUE",
                        "passes_wt_gate": "TRUE",
                    }
                ],
                ROOT_FIELDS,
            )
            write_hashes(output)
            summary = write_summary(output, completed_stage="readiness")
            self.assertTrue(summary["boolean_parse_issues"])
            self.assertFalse(summary["native_pre_reveal_freeze_ready"])

    def test_readiness_distinguishes_all_numerical_from_wt_production_roots(self) -> None:
        panels = _complete_upstream_rows()
        target = next(
            row
            for row in panels
            if row["source_class"] == N_ABS_NKCC
            and float(row["source_scale"]) == 4.0
        )
        numerical_only = next(
            row
            for row in panels
            if row["source_class"] == N_ABS_NKCC
            and float(row["source_scale"]) == 5.0
        )
        target.update(
            {
                "stage": "THIRTY_THREE_START_ROOT_GEOMETRY",
                "attempt_count": 33,
                "numerical_root_count": 1,
                "wt_passing_root_count": 1,
                "confirmation_start_count": 7,
                "geometry_start_count": 33,
            }
        )
        numerical_only["numerical_root_count"] = 1
        standard_ids = ("frozen_reference", "lower_quartile", "upper_quartile")
        attempts = []
        for panel in panels:
            start_ids = (
                (*standard_ids, *(f"uniform_{index:03d}" for index in range(3, 33)))
                if panel["panel_id"] == target["panel_id"]
                else standard_ids
            )
            for start_id in start_ids:
                attempts.append(
                    {
                        "attempt_id": f"{panel['panel_id']}|{start_id}",
                        "panel_id": panel["panel_id"],
                        "source_class": panel["source_class"],
                        "source_scale": panel["source_scale"],
                        "pump_capacity_scale": panel["pump_capacity_scale"],
                        "hydraulic_scale": panel["hydraulic_scale"],
                        "hydraulic_mode": panel["hydraulic_mode"],
                        "eligibility": panel["eligibility"],
                        "start_id": start_id,
                        "root_search_stage": panel["stage"],
                        "root_search_start_count": len(start_ids),
                    }
                )

        def root_row(panel, *, wt):
            return {
                **{
                    key: panel[key]
                    for key in (
                        "panel_id",
                        "topology_id",
                        "routing_id",
                        "source_class",
                        "source_scale",
                        "pump_capacity_scale",
                        "apical_pump_fraction",
                        "apical_k_fraction",
                        "ae4_cation_fraction",
                        "hydraulic_scale",
                        "hydraulic_mode",
                        "eligibility",
                    )
                },
                "root_id": f"{panel['panel_id']}:B00",
                "passes_numerical_gate": "True",
                "passes_wt_gate": "True" if wt else "False",
                "member_start_ids_json": '["frozen_reference"]',
                "root_search_stage": panel["stage"],
                "root_search_start_count": panel["attempt_count"],
                "panel_confirmation_start_count": panel["confirmation_start_count"],
                "panel_geometry_start_count": panel["geometry_start_count"],
            }

        roots = [root_row(target, wt=True), root_row(numerical_only, wt=False)]
        production_wt_ids = [roots[0]["root_id"]]
        production_numerical_ids = sorted(row["root_id"] for row in roots)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            _write_rows(output / "native_source_panel.csv", panels, PANEL_FIELDS)
            _write_rows(
                output / "native_source_root_attempts.csv", attempts, ATTEMPT_FIELDS
            )
            _write_rows(output / "native_source_wt_roots.csv", roots, ROOT_FIELDS)
            write_hashes(output)
            hashes = {
                name: hashlib.sha256((output / name).read_bytes()).hexdigest()
                for name in (
                    "native_source_panel.csv",
                    "native_source_wt_roots.csv",
                    "native_source_root_attempts.csv",
                )
            }
            snapshot_set = hashlib.sha256(
                json.dumps(hashes, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            independent = {
                "status": "INDEPENDENT_NATIVE_ROOT_REPRODUCTION_PASS",
                "root_table_sha256": hashes["native_source_wt_roots.csv"],
                "panel_table_sha256": hashes["native_source_panel.csv"],
                "attempt_table_sha256": hashes["native_source_root_attempts.csv"],
                "artifact_snapshot_sha256": hashes,
                "artifact_snapshot_set_sha256": snapshot_set,
                "artifact_hash_basis": "SHA256_OF_THE_EXACT_BYTES_PARSED_ONCE_FOR_EACH_CSV",
                "row_count": 2,
                "production_numerical_pass_count": 2,
                "independent_pass_count": 2,
                "production_wt_pass_count": 1,
                "independent_wt_pass_count": 1,
                "eligible_production_root_ids": production_numerical_ids,
                "eligible_production_wt_root_ids": production_wt_ids,
                "gate_disagreement_count": 0,
                "wt_gate_disagreement_count": 0,
                "observable_integrity_disagreement_count": 0,
                "duplicate_root_ids": [],
                "nested_gate_failure_count": 0,
                "nested_gate_failures": [],
                "hierarchy_audit": {
                    "status": "PASS",
                    "issue_count": 0,
                    "issues": [],
                    "licensed_pump_contexts": [0.5, 1.0, 2.0],
                    "eligible_production_root_ids": production_numerical_ids,
                    "eligible_production_wt_root_ids": production_wt_ids,
                },
                "optimizer_endpoint_reproduction": {
                    "panel_table_sha256": hashes["native_source_panel.csv"],
                    "attempt_table_sha256": hashes[
                        "native_source_root_attempts.csv"
                    ],
                    "convergence_gate_disagreement_count": 0,
                    "panel_start_design_disagreement_count": 0,
                    "root_cluster_count_disagreement_panel_ids": [],
                    "duplicate_panel_ids": [],
                    "duplicate_attempt_ids": [],
                    "all_declared_confirmation_panels_have_exactly_7_current_starts": True,
                    "all_declared_geometry_panels_have_exactly_33_current_starts": True,
                },
            }
            (output / "independent_native_summary.json").write_text(
                json.dumps(independent), encoding="utf-8"
            )
            summary = write_summary(output, completed_stage="readiness")
            self.assertTrue(summary["independent_root_reproduction_complete"])
            self.assertTrue(summary["native_pre_reveal_freeze_ready"])
            self.assertEqual(len(summary["production_numerical_root_ids"]), 2)
            self.assertEqual(len(summary["production_root_ids"]), 1)

    def test_geometry_rejects_duplicate_failed_and_unconfirmed_intake(self) -> None:
        rows = _complete_upstream_rows()
        target = next(
            row for row in rows if row["source_class"] == N_ABS_NKCC
        )
        target_id = target["panel_id"]
        target["wt_passing_root_count"] = 1
        target["confirmation_start_count"] = 7
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            path = output / "native_source_panel.csv"

            def run_geometry(panel_ids):
                execute_stage(
                    "geometry",
                    output=output,
                    workers=1,
                    start_count=33,
                    max_nfev=1,
                    seed=1,
                    panel_ids=panel_ids,
                )

            _write_rows(path, rows, PANEL_FIELDS)
            with patch(
                "modern_full_model.run_native_source_panel.run_panel",
                side_effect=AssertionError("rejected geometry must not solve"),
            ):
                with self.assertRaisesRegex(RuntimeError, "duplicate --panel-id"):
                    run_geometry((target_id, target_id))

                target["wt_passing_root_count"] = 0
                _write_rows(path, rows, PANEL_FIELDS)
                with self.assertRaisesRegex(RuntimeError, "rejects diagnostic/H>1"):
                    run_geometry((target_id,))

                target["wt_passing_root_count"] = 1
                target["confirmation_start_count"] = 0
                _write_rows(path, rows, PANEL_FIELDS)
                with self.assertRaisesRegex(RuntimeError, "rejects diagnostic/H>1"):
                    run_geometry((target_id,))


if __name__ == "__main__":
    unittest.main()
