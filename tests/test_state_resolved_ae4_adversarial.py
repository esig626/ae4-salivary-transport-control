"""Adversarial regression checks for the Task-13 AE4 investigation.

These tests encode failures that could otherwise turn back into apparently
successful fits: the 2016 observable is total alkalinization, the explicit
carrier states are quasi-steady rather than additional dynamic states, and a
remote null root must not be mistaken for a physiological alternate branch.
"""

from __future__ import annotations

import csv
from pathlib import Path
import unittest

import numpy as np

from ae4_mechanism_reconstruction.chassis import AE4_KNOCKOUT, FixedChassis
from state_resolved_ae4.adapter import StateResolvedAE4Adapter
from state_resolved_ae4.cycles import (
    CANDIDATE_DEFINITIONS,
    CycleParameters,
    evaluate_cycle,
)
from state_resolved_ae4.fitting import assay_rates, nominal_assay_context


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
RESULT_ROOT = REPOSITORY_ROOT / "results" / "13_state_resolved_ae4"


def _rows(name: str) -> list[dict[str, str]]:
    with (RESULT_ROOT / name).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if any(None in row for row in rows):
        raise AssertionError(f"{name} is not rectangular")
    return rows


class AssayObservableAudit(unittest.TestCase):
    def test_shared_cycle_uses_total_net_hco3_source_not_one_branch(self) -> None:
        definition = CANDIDATE_DEFINITIONS["SR2_SHARED_112"]
        parameters = CycleParameters()
        concentration = 50.0
        evaluation = evaluate_cycle(
            definition,
            nominal_assay_context("na", concentration),
            parameters,
        )
        reported = float(assay_rates(definition, parameters, "na", (concentration,))[0])
        expected = max(0.0, float(evaluation.intracellular_sources["hco3_i"]))

        self.assertAlmostEqual(reported, expected, places=14)
        # This guards the original leakage-by-observable bug: in a shared
        # network the other branch can carry current even when its external
        # cation is nearly zero, and the stoichiometric HCO3 coefficient is
        # not one.  A selected branch magnitude is therefore not the BCECF
        # alkalinization observable.
        self.assertGreater(
            abs(reported - abs(float(evaluation.branch_currents["na"]))),
            0.05,
        )


class QssNestedModelAudit(unittest.TestCase):
    def test_state_resolved_adapter_adds_no_dynamic_carrier_state(self) -> None:
        definition = CANDIDATE_DEFINITIONS["SR2_SHARED_112"]
        chassis = FixedChassis(
            StateResolvedAE4Adapter(definition, CycleParameters(capacity=1e-6))
        )
        state = np.asarray((118.7, 5.6, 28.7, 25.0, 120.0, 50.0, 10.0))

        first = chassis.raw_rhs(0.0, state)
        second = chassis.raw_rhs(0.0, state)
        self.assertEqual(first.shape, (7,))
        np.testing.assert_array_equal(first, second)
        # The carrier occupancies are algebraic diagnostics.  Consequently
        # the embedded model is exactly representable by static scalar branch
        # functions J_Na(x) and J_K(x); it cannot supply transporter memory.
        _raw, diagnostics = chassis.evaluate(state, calcium=0.05)
        self.assertIn("occupancy", diagnostics["ae4_diagnostics"])
        self.assertEqual(
            len(diagnostics["ae4_diagnostics"]["occupancy"]),
            len(definition.states),
        )


class AlternateRootAudit(unittest.TestCase):
    def test_remote_exact_null_root_is_explicitly_nonphysiological(self) -> None:
        # Found by the independent log-coordinate scipy.root route from seed
        # 1307.  A second independent start drifts along the same asymptotic
        # branch to a different enormous height.  It is an exact positive
        # mathematical root, but lies far outside the declared 2--200 um
        # physiological root-search domain.
        remote = np.asarray(
            (
                118.63335900824212,
                5.579964665609247,
                1.2512524902179945e14,
                16.87458314982318,
                128.04629969010358,
                55.31008452157387,
                89.61082531822470,
            )
        )
        residual = FixedChassis().raw_rhs(
            0.0, remote, scenario=AE4_KNOCKOUT
        )

        self.assertLess(float(np.linalg.norm(residual, ord=np.inf)), 3e-15)
        self.assertGreater(float(remote[2]), 200.0)
        self.assertGreater(float(remote[2]), 1e10)


class FinalLedgerFirewallAudit(unittest.TestCase):
    def test_every_round_five_row_keeps_strict_secretion_sealed(self) -> None:
        rows = _rows("chassis_reconstruction_scan.csv")
        self.assertEqual(len(rows), 33)
        self.assertTrue(all(row["secretion_target_read"] == "False" for row in rows))

        attempts = [row for row in rows if row["record_type"] == "wt_capacity_attempt"]
        self.assertEqual(len(attempts), 5)
        self.assertTrue(all(row["valid_wt_capacity_available"] == "False" for row in attempts))
        self.assertTrue(all(float(row["max_abs_raw_rhs"]) > 1e-8 for row in attempts))
        self.assertTrue(all(float(row["height_um"]) > 199.99 for row in attempts))

        resolved = [row for row in rows if row["record_type"] == "resolved_null_root"]
        self.assertEqual(len(resolved), 8)
        self.assertTrue(all(float(row["max_abs_raw_rhs"]) <= 2.5e-12 for row in resolved))
        self.assertGreater(min(float(row["normalized_ionic_distance"]) for row in resolved), 24.0)

    def test_strict_references_are_never_joined_to_a_prediction(self) -> None:
        rows = _rows("dynamic_trajectory_status.csv")
        self.assertTrue(all(row["strict_target_used_in_fit"] == "false" for row in rows))
        self.assertTrue(
            all(row["strict_target_joined_to_prediction"] == "false" for row in rows)
        )
        candidates = [row for row in rows if row["record_type"] == "validation_candidate"]
        self.assertEqual(len(candidates), 10)
        for row in candidates:
            self.assertEqual(row["status"], "blocked_before_dynamics")
            self.assertEqual(row["endpoint_ko_wt_ratio"], "")
            self.assertEqual(row["integral_ko_wt_ratio"], "")

    def test_adversarial_machine_ledger_is_rectangular_and_target_blind(self) -> None:
        rows = _rows("adversarial_checks.csv")
        self.assertGreaterEqual(len(rows), 24)
        self.assertEqual(len({row["check_id"] for row in rows}), len(rows))
        self.assertTrue(all(row["strict_secretion_read"] == "False" for row in rows))


if __name__ == "__main__":
    unittest.main()
