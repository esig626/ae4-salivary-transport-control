"""Focused tests for the WT-only cross-root G5 repair profiler."""

from __future__ import annotations

from pathlib import Path
import unittest

from modern_full_model.run_g5_cross_root_profile import (
    declared_calcium_probes,
    declared_cross_root_nkcc_panel,
    load_retained_wt_roots,
)
from modern_full_model.run_g5_joint_repair_profile import (
    declared_joint_repair_profiles,
)
from modern_full_model.nkcc_stimulation import N1AlgebraicNkcc1


class TestDeclaredRepairPanel(unittest.TestCase):
    def test_cross_root_panel_spans_nested_family_gain_and_timing_axes(self) -> None:
        panel = declared_cross_root_nkcc_panel()
        self.assertEqual(len(panel), 16)
        self.assertEqual(
            {member.family for member in panel},
            {"N0_BASELINE", "N1_CCH_ALGEBRAIC", "N2_CCH_EFFECTIVE_DYNAMIC"},
        )
        self.assertEqual(
            {
                member.tau_activation_s
                for member in panel
                if member.family == "N2_CCH_EFFECTIVE_DYNAMIC"
            },
            {10.0, 30.0, 120.0, 300.0},
        )
        self.assertEqual(
            {member.fully_activated_multiplier for member in panel},
            {1.0, 1.75, 3.0, 8.0},
        )
        self.assertTrue(
            all(
                "OUTSIDE_P15" in member.gain_evidence_status
                for member in panel
                if member.fully_activated_multiplier > 1.75
            )
        )

    def test_only_historical_calcium_anchor_can_select_a_candidate(self) -> None:
        probes = declared_calcium_probes()
        eligible = [probe for probe in probes if probe.eligible_for_selection]
        self.assertEqual(len(eligible), 1)
        self.assertEqual(eligible[0].stimulated_calcium_uM, 0.55)
        self.assertEqual(
            {probe.stimulated_calcium_uM for probe in probes}, {0.35, 0.55, 0.75}
        )
        self.assertTrue(
            all(
                "TIER_3" in probe.evidence_status
                for probe in probes
                if not probe.eligible_for_selection
            )
        )

    def test_joint_profiles_normalize_nkcc_to_each_effective_ca_endpoint(self) -> None:
        profiles = declared_joint_repair_profiles()
        self.assertEqual(
            {
                (profile.effective_calcium_uM, profile.nkcc1_fully_activated_multiplier)
                for profile in profiles
            },
            {(0.10, 1.75), (0.15, 1.75), (0.12, 2.0)},
        )
        for profile in profiles:
            law = N1AlgebraicNkcc1(
                fully_activated_multiplier=profile.nkcc1_fully_activated_multiplier,
                resting_calcium_uM=0.058,
                stimulated_calcium_uM=profile.effective_calcium_uM,
            )
            output = law.evaluate(
                1.0e-6, (), calcium_uM=profile.effective_calcium_uM
            )
            self.assertAlmostEqual(
                output.capacity_multiplier,
                profile.nkcc1_fully_activated_multiplier,
            )


class TestFrozenRootIntake(unittest.TestCase):
    def test_exactly_nine_rest_passing_roots_are_reconstructed(self) -> None:
        repository = Path(__file__).resolve().parents[1]
        roots = load_retained_wt_roots(
            repository / "results" / "13B_modern_full_model"
        )
        self.assertEqual(len(roots), 9)
        self.assertEqual(len({root.root_id for root in roots}), 9)
        self.assertEqual(sum(root.selected_reference for root in roots), 1)
        self.assertTrue(all(len(root.core_state) == 12 for root in roots))
        self.assertEqual(
            {root.calibration.common_membrane_conductance_scale for root in roots},
            {10.0},
        )


if __name__ == "__main__":
    unittest.main()
