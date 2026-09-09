"""Pre-registered adversarial invariants for Task 13B.

These tests do not inspect the held-out AE4-null secretion values.  The sealed
file is checked only as opaque bytes against its pre-reveal digest.
"""

from __future__ import annotations

import csv
from dataclasses import replace
import hashlib
from pathlib import Path
import unittest

import numpy as np

from modern_full_model import membranes, model, parameters, water
from modern_full_model.transporters import AE4Parameters


REPOSITORY = Path(__file__).resolve().parents[1]
RESULTS = REPOSITORY / "results" / "13B_modern_full_model"
SEALED_HELDOUT_SHA256 = (
    "939f684ab92371f35b68fb52fbac59ebcb654fb5bf9a657feb4fd543ae1a0587"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


class TestPreRevealFirewall(unittest.TestCase):
    def test_sealed_heldout_artifact_is_byte_identical_to_freeze(self) -> None:
        # Deliberately hash opaque bytes; do not parse or print target values.
        self.assertEqual(
            _sha256(RESULTS / "heldout_targets.csv"),
            SEALED_HELDOUT_SHA256,
        )

    def test_calibration_ledger_contains_no_ae4_null_secretion_row(self) -> None:
        with (RESULTS / "calibration_targets.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        offenders = []
        for row in rows:
            genotype = row.get("genotype", "").casefold()
            measurement = row.get("measurement", "").casefold()
            is_ae4_null = "ae4" in genotype and any(
                token in genotype for token in ("null", "knockout", " ko")
            )
            is_secretion = any(
                token in measurement for token in ("saliva", "secretion", "flow")
            )
            if is_ae4_null and is_secretion:
                offenders.append(row.get("target_id", "UNNAMED"))
        self.assertEqual(offenders, [])

    def test_selection_modules_do_not_depend_on_sealed_or_legacy_holdout_files(self) -> None:
        sensitive_names = {
            "calibration.py",
            "camp_pka.py",
            "membranes.py",
            "model.py",
            "parameters.py",
            "states.py",
            "transporters.py",
            "water.py",
        }
        forbidden_fragments = (
            "heldout_targets.csv",
            "HELDOUT_AE4_KO_FLOW_RATIO",
            "13_state_resolved_ae4/heldout_validation",
            "13_state_resolved_ae4/evidence_freeze",
        )
        offenders: list[str] = []
        for path in sorted((REPOSITORY / "src" / "modern_full_model").glob("*.py")):
            if path.name not in sensitive_names:
                continue
            source = path.read_text(encoding="utf-8")
            if any(fragment in source for fragment in forbidden_fragments):
                offenders.append(path.name)
        self.assertEqual(offenders, [])


class TestAdversarialConservationAndSigns(unittest.TestCase):
    def setUp(self) -> None:
        ae4 = AE4Parameters.equal_branch_attempts(
            carrier_amount_fmol=1.0e-4,
            common_cl_attempt_rate_s=1.0,
            loaded_attempt_rate_s=1.0,
        )
        self.full = model.ModernFullModel(ae4_parameters=ae4)
        self.reference = self.full.initial_state()

    def test_cell_charge_rate_is_an_identity_even_off_electroneutral_manifold(self) -> None:
        # This deliberately changes cell bulk charge.  The zero rate proves a
        # conserved charge manifold and prevents an unconstrained root search
        # from claiming branch uniqueness across different charge leaves.
        perturbed = self.reference.copy()
        perturbed[0] += 1.0
        evaluation = self.full.evaluate(0.0, perturbed)
        self.assertGreater(
            abs(evaluation.diagnostics.state_charge_fmol["cell"]), 0.5
        )
        self.assertLess(
            abs(
                evaluation.diagnostics.conservation_residuals[
                    "cell_bulk_charge_rate_fmol_s"
                ]
            ),
            1.0e-12,
        )

    def test_conservation_and_two_membrane_closure_survive_local_state_perturbations(self) -> None:
        rng = np.random.default_rng(1302)
        for _ in range(12):
            state = self.reference.copy()
            # Scale all amounts and the matching volume together within each
            # compartment, then perturb cations slightly.  This keeps pH in
            # the declared bracket while exercising non-reference currents.
            state[:6] *= float(rng.uniform(0.85, 1.15))
            state[6:12] *= float(rng.uniform(0.85, 1.15))
            state[0] *= float(rng.uniform(0.98, 1.02))
            state[7] *= float(rng.uniform(0.98, 1.02))
            diagnostics = self.full.evaluate(0.0, state).diagnostics
            self.assertLess(abs(diagnostics.membranes.current_residuals_A["apical"]), 1e-23)
            self.assertLess(
                abs(diagnostics.membranes.current_residuals_A["basolateral"]),
                1e-23,
            )
            self.assertLess(
                abs(diagnostics.conservation_residuals["carbon_accounting_fmol_s"]),
                1e-12,
            )
            self.assertLess(
                abs(diagnostics.conservation_residuals["water_volume_accounting_pL_s"]),
                1e-12,
            )

    def test_reversible_homeostasis_fluxes_follow_their_affinities(self) -> None:
        p = parameters.FullModelParameters()
        rng = np.random.default_rng(1303)
        for _ in range(40):
            values = np.exp(rng.uniform(np.log(0.1), np.log(180.0), size=10))
            flux = membranes.evaluate_homeostasis(
                membranes.HomeostasisEnvironment(*map(float, values)), p
            )
            for name, rate in (
                ("NKCC1", flux.nkcc1_inward_fmol_s),
                ("NHE1", flux.nhe1_inward_fmol_s),
                ("AE2", flux.ae2_inward_fmol_s),
            ):
                self.assertGreaterEqual(rate * flux.affinities[name], -1.0e-14)

    def test_water_flux_signs_point_toward_the_hyperosmotic_compartment(self) -> None:
        p = parameters.FullModelParameters()
        cell_high = water.evaluate_water_fluxes(
            water.OsmoticEnvironment(310.0, 300.0, 290.0, 0.1), p
        )
        self.assertGreater(cell_high.bath_to_cell_pL_s, 0.0)
        self.assertLess(cell_high.cell_to_lumen_pL_s, 0.0)
        self.assertGreater(cell_high.bath_to_lumen_pL_s, 0.0)

        lumen_high = water.evaluate_water_fluxes(
            water.OsmoticEnvironment(300.0, 320.0, 290.0, 0.1), p
        )
        self.assertGreater(lumen_high.cell_to_lumen_pL_s, 0.0)
        self.assertGreater(lumen_high.bath_to_lumen_pL_s, 0.0)

    def test_extreme_apical_partitions_preserve_current_and_charge_closure(self) -> None:
        base = parameters.FullModelParameters()
        for fraction in (0.0, 1.0):
            configured = replace(
                base,
                membranes=replace(
                    base.membranes,
                    apical_pump_fraction=fraction,
                    apical_k_fraction=fraction,
                    g_apical_background_S=0.0,
                ),
            )
            full = model.ModernFullModel(parameters=configured)
            diagnostics = full.evaluate(0.0, full.initial_state()).diagnostics.membranes
            for residual in diagnostics.current_residuals_A.values():
                self.assertLess(abs(residual), 1.0e-23)
            for residual in diagnostics.charge_residuals_fmol_s.values():
                self.assertLess(abs(residual), 1.0e-12)


class TestAssumptionLabelsAndRootDiscipline(unittest.TestCase):
    def test_apical_fractions_are_explicit_unmeasured_modeling_decisions(self) -> None:
        records = {
            record.name: record
            for record in parameters.parameter_records(parameters.FullModelParameters())
        }
        for name in (
            "membranes.apical_pump_fraction",
            "membranes.apical_k_fraction",
        ):
            self.assertEqual(
                records[name].provenance,
                parameters.Provenance.NEW_MODELING_DECISION.value,
            )
            self.assertIn("not measured", records[name].note.casefold())

    def test_bundled_root_helper_is_never_production_eligible(self) -> None:
        full = model.ModernFullModel()
        result = full.solve_resting_root(start_count=1, max_nfev=5)
        self.assertFalse(result.production_eligible)
        self.assertIn("charge manifold monitored but not parameterized", result.limitations)
        self.assertIn("branch uniqueness not established", result.limitations)


if __name__ == "__main__":
    unittest.main()
