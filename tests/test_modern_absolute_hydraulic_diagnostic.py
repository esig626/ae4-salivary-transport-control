"""Focused tests for the isolated NABS-S4 by H2.5 WT diagnostic."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from modern_full_model.absolute_hydraulic_diagnostic import (
    ATTEMPT_NAME,
    COMMIT_NAME,
    DIAGNOSTIC_ELIGIBILITY,
    GEOMETRY_STAGE,
    HASH_NAME,
    HYDRAULIC_SCALE,
    PANEL_NAME,
    ROOT_NAME,
    SCREEN_STAGE,
    SOURCE_MAP_NAME,
    SUMMARY_NAME,
    AbsoluteHydraulicDiagnosticRoot,
    build_absolute_hydraulic_diagnostic_model,
    declared_absolute_hydraulic_diagnostic_panel,
    diagnostic_member_from_root,
    load_absolute_hydraulic_diagnostic_contract,
    run_absolute_hydraulic_diagnostic,
    _parent_anchor_maps,
    _validate_reports,
)
from modern_full_model.native_source_panel import (
    FROZEN_NKCC1_CAPACITY_FMOL_S,
    N_ABS_NKCC,
    NativeRootAttempt,
    NativeRootRecord,
    NativeRootSearchReport,
    build_native_source_model,
    declared_native_panel,
    root_to_csv_row,
)


def _start_ids(count: int) -> tuple[str, ...]:
    base = ("frozen_reference", "lower_quartile", "upper_quartile")
    return base + tuple(f"uniform_{index:03d}" for index in range(3, count))


def _root(
    member,
    *,
    wt_pass: bool,
    start_ids: tuple[str, ...],
) -> NativeRootRecord:
    return NativeRootRecord(
        root_id=f"{member.panel_id}:B00",
        panel_id=member.panel_id,
        topology_id=member.topology_id,
        routing_id=member.routing_id,
        source_class=member.source_class,
        source_scale=member.source_scale,
        pump_capacity_scale=member.pump_capacity_scale,
        apical_pump_fraction=member.apical_pump_fraction,
        apical_k_fraction=member.apical_k_fraction,
        ae4_cation_fraction=member.ae4_cation_fraction,
        hydraulic_scale=member.hydraulic_scale,
        hydraulic_mode=member.hydraulic_mode,
        eligibility=member.eligibility,
        passes_numerical_gate=True,
        passes_wt_gate=wt_pass,
        member_start_ids=start_ids,
        core_state=tuple(float(index + 1) for index in range(12)),
        coordinates=tuple(float(index + 1) for index in range(10)),
        other_impermeant_osmoles_fmol=125.0,
        observables={"cl_i_mM": 50.1, "ph_i": 6.91},
        max_abs_scaled_residual=1.0e-12,
        max_abs_amount_rhs_fmol_s=1.0e-13,
        max_abs_volume_rhs_pL_s=1.0e-15,
        max_abs_omitted_charge_rhs_fmol_equivalent_s=1.0e-13,
        max_abs_charge_residual_fmol=1.0e-14,
        max_abs_current_residual_A=1.0e-25,
        normalized_jacobian_singular_values=(1.0,) * 11,
        normalized_jacobian_rank=11,
        normalized_jacobian_nullity=0,
        boundary_hits=(),
        gate_failures=() if wt_pass else ("cl_i_mM outside",),
    )


def _report(member, *, count: int, wt_pass: bool) -> NativeRootSearchReport:
    ids = _start_ids(count)
    attempts = tuple(
        NativeRootAttempt(
            panel_id=member.panel_id,
            start_id=start_id,
            optimizer_success=True,
            admissible_state=True,
            converged_root=True,
            passes_numerical_gate=True,
            passes_wt_gate=wt_pass,
            max_abs_scaled_residual=1.0e-12,
            max_abs_amount_rhs_fmol_s=1.0e-13,
            max_abs_volume_rhs_pL_s=1.0e-15,
            max_abs_omitted_charge_rhs_fmol_equivalent_s=1.0e-13,
            normalized_jacobian_rank=11,
            normalized_jacobian_nullity=0,
            boundary_hits=(),
            wt_gate_failures=() if wt_pass else ("cl_i_mM outside",),
            coordinates_and_other=tuple(float(index + 1) for index in range(11)),
            cost=0.0,
            optimality=0.0,
            nfev=1,
            message="synthetic",
        )
        for start_id in ids
    )
    return NativeRootSearchReport(
        member=member,
        attempts=attempts,
        roots=(_root(member, wt_pass=wt_pass, start_ids=ids),),
        start_design={"count": count, "start_ids": list(ids)},
    )


class _SyntheticRunner:
    def __init__(self, *, no_passes: bool = False) -> None:
        self.calls: list[tuple[int, tuple[str, ...]]] = []
        self.no_passes = no_passes

    def __call__(self, members, *, start_count, **kwargs):
        members = tuple(members)
        self.calls.append((start_count, tuple(member.panel_id for member in members)))
        if start_count == 3:
            passing = set() if self.no_passes else {members[0].panel_id, members[1].panel_id}
        elif start_count == 7:
            passing = {members[0].panel_id}
        elif start_count == 33:
            passing = {members[0].panel_id}
        else:
            raise AssertionError("unexpected start count")
        return tuple(
            _report(member, count=start_count, wt_pass=member.panel_id in passing)
            for member in members
        )


class TestAbsoluteHydraulicDiagnostic(unittest.TestCase):
    def test_exact_namespaced_cartesian_contract_and_h1_parents(self) -> None:
        members = declared_absolute_hydraulic_diagnostic_panel()
        self.assertEqual(len(members), 90)
        self.assertEqual(len({member.panel_id for member in members}), 90)
        self.assertTrue(all(member.panel_id.startswith("ABS_HYD_DIAG_") for member in members))
        self.assertEqual(Counter(member.hydraulic_mode for member in members), {"HW": 45, "HWQ": 45})
        self.assertEqual(Counter(member.pump_capacity_scale for member in members), {0.5: 30, 1.0: 30, 2.0: 30})
        self.assertEqual(len({member.topology_id for member in members}), 5)
        self.assertEqual(len({member.routing_id for member in members}), 3)
        for member in members:
            self.assertEqual(member.source_class, N_ABS_NKCC)
            self.assertEqual(member.source_scale, 4.0)
            self.assertEqual(member.hydraulic_scale, 2.5)
            self.assertEqual(member.eligibility, DIAGNOSTIC_ELIGIBILITY)
            self.assertEqual(
                member.parent_panel_id,
                f"N_ABS_NKCC_S4_{member.topology_id}_{member.routing_id}_P{member.pump_capacity_scale:g}_H1",
            )

    def test_hw_and_hwq_change_only_declared_water_direction(self) -> None:
        members = declared_absolute_hydraulic_diagnostic_panel()
        hw = next(member for member in members if member.hydraulic_mode == "HW")
        hwq = next(
            member
            for member in members
            if member.hydraulic_mode == "HWQ"
            and member.topology_id == hw.topology_id
            and member.routing_id == hw.routing_id
            and member.pump_capacity_scale == hw.pump_capacity_scale
        )
        h1 = next(
            member
            for member in declared_native_panel(
                source_class=N_ABS_NKCC,
                source_scale=4.0,
                pump_capacity_scale=hw.pump_capacity_scale,
            )
            if member.topology_id == hw.topology_id and member.routing_id == hw.routing_id
        )
        models = [
            build_native_source_model(member, other_impermeant_osmoles_fmol=125.0)
            for member in (h1, hw, hwq)
        ]
        base, hw_model, hwq_model = models
        self.assertEqual(
            base.parameters.homeostasis.nkcc1_capacity_fmol_s,
            4.0 * FROZEN_NKCC1_CAPACITY_FMOL_S,
        )
        for candidate in (hw_model, hwq_model):
            self.assertEqual(asdict(candidate.parameters.membranes), asdict(base.parameters.membranes))
            self.assertEqual(asdict(candidate.parameters.homeostasis), asdict(base.parameters.homeostasis))
            self.assertEqual(asdict(candidate.ae4_parameters), asdict(base.ae4_parameters))
            self.assertEqual(
                candidate.parameters.water.apical_hydraulic_pL_s_mOsm,
                HYDRAULIC_SCALE * base.parameters.water.apical_hydraulic_pL_s_mOsm,
            )
        self.assertEqual(
            hw_model.parameters.water.outflow_rate_s,
            base.parameters.water.outflow_rate_s,
        )
        self.assertEqual(
            hwq_model.parameters.water.outflow_rate_s,
            HYDRAULIC_SCALE * base.parameters.water.outflow_rate_s,
        )

    def test_corresponding_h1_root_is_used_as_parent_anchor(self) -> None:
        member = declared_absolute_hydraulic_diagnostic_panel()[0]
        parent = declared_native_panel(
            source_class=N_ABS_NKCC,
            source_scale=4.0,
            pump_capacity_scale=member.pump_capacity_scale,
        )[0]
        parent = parent.__class__(**{**asdict(parent), "panel_id": member.parent_panel_id})
        root = _root(parent, wt_pass=False, start_ids=("frozen_reference",))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "native_source_wt_roots.csv"
            row = root_to_csv_row(root)
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=tuple(row))
                writer.writeheader()
                writer.writerow(row)
            anchors, ids = _parent_anchor_maps(Path(temporary), (member,))
        self.assertEqual(ids[member.panel_id], (root.root_id,))
        self.assertEqual(
            anchors[member.panel_id],
            ((*root.coordinates, root.other_impermeant_osmoles_fmol),),
        )

    def test_exact_start_validation_rejects_shallow_or_duplicate_sets(self) -> None:
        member = declared_absolute_hydraulic_diagnostic_panel()[0]
        valid = _report(member, count=3, wt_pass=False)
        _validate_reports((valid,), (member,), stage=SCREEN_STAGE)
        shallow = NativeRootSearchReport(
            member=member,
            attempts=valid.attempts[:2],
            roots=valid.roots,
            start_design={"count": 2, "start_ids": list(_start_ids(3)[:2])},
        )
        with self.assertRaisesRegex(RuntimeError, "requires 3 unique starts"):
            _validate_reports((shallow,), (member,), stage=SCREEN_STAGE)

    def test_full_fake_hierarchy_is_isolated_atomic_and_loadable(self) -> None:
        runner = _SyntheticRunner()
        snapshot = {
            "ledger_sha256": "a" * 64,
            "artifacts": {},
            "native_root_table_sha256": "b" * 64,
            "dynamic_failure_gate": {"sha256": "c" * 64},
        }
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            protected_names = (
                "native_source_panel.csv",
                "native_source_root_attempts.csv",
                "native_source_wt_roots.csv",
                "native_source_wt_summary.json",
                "native_source_hashes.csv",
            )
            for name in protected_names:
                (output / name).write_bytes(f"sentinel:{name}".encode())
            before = {
                name: hashlib.sha256((output / name).read_bytes()).hexdigest()
                for name in protected_names
            }
            with patch(
                "modern_full_model.absolute_hydraulic_diagnostic._native_snapshot",
                return_value=snapshot,
            ), patch(
                "modern_full_model.absolute_hydraulic_diagnostic._parent_anchor_maps",
                return_value=({}, {}),
            ):
                summary = run_absolute_hydraulic_diagnostic(
                    output=output,
                    native_output=output,
                    workers=1,
                    panel_runner=runner,
                )
            after = {
                name: hashlib.sha256((output / name).read_bytes()).hexdigest()
                for name in protected_names
            }
            self.assertEqual(before, after)
            self.assertEqual([count for count, _ in runner.calls], [3, 7, 33])
            self.assertEqual([len(ids) for _, ids in runner.calls], [90, 2, 1])
            self.assertEqual(summary["screen_attempt_count"], 270)
            self.assertEqual(summary["confirmation_attempt_count"], 14)
            self.assertEqual(summary["geometry_attempt_count"], 33)
            self.assertEqual(summary["attempt_history_row_count"], 317)
            self.assertEqual(summary["geometry_wt_passing_root_count"], 1)
            for name in (
                SOURCE_MAP_NAME,
                PANEL_NAME,
                ATTEMPT_NAME,
                ROOT_NAME,
                SUMMARY_NAME,
                HASH_NAME,
                COMMIT_NAME,
            ):
                self.assertTrue((output / name).exists())
            contract = load_absolute_hydraulic_diagnostic_contract(output)
            self.assertEqual(len(contract.roots), 1)
            self.assertIsInstance(contract.roots[0], AbsoluteHydraulicDiagnosticRoot)
            self.assertFalse(contract.production_promotion_authorized)
            self.assertFalse(contract.heldout_reveal_authorized)
            member = diagnostic_member_from_root(contract.roots[0])
            self.assertEqual(member.eligibility, DIAGNOSTIC_ELIGIBILITY)
            with self.assertRaises(TypeError):
                diagnostic_member_from_root(contract.roots[0].record)  # type: ignore[arg-type]
            with self.assertRaises((AttributeError, TypeError)):
                build_native_source_model(  # type: ignore[arg-type]
                    contract.roots[0],
                    other_impermeant_osmoles_fmol=125.0,
                )
            model = build_absolute_hydraulic_diagnostic_model(contract.roots[0])
            self.assertEqual(model.parameters.water.outflow_rate_s > 0.0, True)

    def test_empty_pass_set_is_explicit_and_never_promoted(self) -> None:
        runner = _SyntheticRunner(no_passes=True)
        snapshot = {
            "ledger_sha256": "a" * 64,
            "artifacts": {},
            "native_root_table_sha256": "b" * 64,
            "dynamic_failure_gate": {"sha256": "c" * 64},
        }
        with tempfile.TemporaryDirectory() as temporary, patch(
            "modern_full_model.absolute_hydraulic_diagnostic._native_snapshot",
            return_value=snapshot,
        ), patch(
            "modern_full_model.absolute_hydraulic_diagnostic._parent_anchor_maps",
            return_value=({}, {}),
        ):
            output = Path(temporary)
            summary = run_absolute_hydraulic_diagnostic(
                output=output,
                native_output=output,
                workers=1,
                panel_runner=runner,
            )
            contract = load_absolute_hydraulic_diagnostic_contract(output)
        self.assertEqual([count for count, _ in runner.calls], [3])
        self.assertTrue(summary["hierarchy_complete"])
        self.assertFalse(summary["diagnostic_dynamic_contract_available"])
        self.assertEqual(contract.roots, ())

    def test_loader_fails_closed_without_commit_or_after_tamper(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with self.assertRaisesRegex(RuntimeError, "commit marker is missing"):
                load_absolute_hydraulic_diagnostic_contract(output)

        runner = _SyntheticRunner()
        snapshot = {
            "ledger_sha256": "a" * 64,
            "artifacts": {},
            "native_root_table_sha256": "b" * 64,
            "dynamic_failure_gate": {"sha256": "c" * 64},
        }
        with tempfile.TemporaryDirectory() as temporary, patch(
            "modern_full_model.absolute_hydraulic_diagnostic._native_snapshot",
            return_value=snapshot,
        ), patch(
            "modern_full_model.absolute_hydraulic_diagnostic._parent_anchor_maps",
            return_value=({}, {}),
        ):
            output = Path(temporary)
            run_absolute_hydraulic_diagnostic(
                output=output,
                native_output=output,
                workers=1,
                panel_runner=runner,
            )
            with (output / ROOT_NAME).open("ab") as handle:
                handle.write(b"\n")
            with self.assertRaisesRegex(RuntimeError, "artifact is stale"):
                load_absolute_hydraulic_diagnostic_contract(output)


if __name__ == "__main__":
    unittest.main()
