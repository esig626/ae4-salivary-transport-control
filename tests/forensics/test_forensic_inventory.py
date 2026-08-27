"""Regression checks for Phase 11 archive and inventory claims."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "archive/legacy-2017"
ANALYSIS = ROOT / "analysis/11_forensic_reconstruction"
RESULTS = ROOT / "results/11_forensic_reconstruction"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_helper(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ForensicInventoryTests(unittest.TestCase):
    def test_archive_ledger_accounts_for_every_file(self) -> None:
        with (RESULTS / "archive_hashes.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        actual = {
            str(path.relative_to(ARCHIVE))
            for path in ARCHIVE.rglob("*")
            if path.is_file()
        }
        recorded = {row["path"] for row in rows}
        self.assertEqual(len(rows), 67)
        self.assertEqual(len(recorded), 67)
        self.assertEqual(recorded, actual)

    def test_all_materialized_archive_hash_evidence_is_unchanged(self) -> None:
        with (RESULTS / "archive_hashes.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        exact = [row for row in rows if row["integrity_basis"] == "byte-exact local"]
        remote = [row for row in rows if row["integrity_basis"] != "byte-exact local"]
        self.assertEqual(len(exact), 65)
        self.assertEqual(len(remote), 2)
        for row in exact:
            path = ARCHIVE / row["path"]
            self.assertEqual(path.stat().st_size, int(row["authoritative_size_bytes"]))
            self.assertEqual(sha256(path), row["sha256"])
            self.assertEqual(row["local_transport_sha256"], row["sha256"])
        exceptions = {row["path"]: row for row in remote}
        paper = exceptions[
            "A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion.pdf"
        ]
        self.assertEqual(paper["authoritative_size_bytes"], "1108802")
        self.assertEqual(paper["git_blob_sha"], "04ba5733c02e2cbb0f4ed642b9a8944c0e4ba267")
        log = exceptions["Ae4_Basis.log"]
        self.assertEqual(log["authoritative_size_bytes"], "51957")
        self.assertEqual(log["git_blob_sha"], "aebe7db9e214d3c3c5a3c40d840dfcd13f57cb7f")

    def test_par_mat_dump_regenerates_exactly(self) -> None:
        helper = load_helper("dump_par_mat", ANALYSIS / "tools/dump_par_mat.py")
        generated = helper.build_inventory()
        expected = json.loads((RESULTS / "par_mat_inventory.json").read_text(encoding="utf-8"))
        self.assertEqual(generated, expected)
        self.assertEqual(generated["mat_file"]["stored_variables"][0]["field_count"], 45)
        self.assertEqual(len(generated["fields"]), 45)
        self.assertEqual(generated["comparison_to_par_m"]["numeric_mismatches"], [])

    def test_figure_mapping_accounts_for_all_46_historical_plot_files(self) -> None:
        with (RESULTS / "figure_mapping.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        actual = {
            str(path.relative_to(ROOT))
            for path in ARCHIVE.iterdir()
            if path.name.startswith("Figure") and path.suffix.lower() in {".eps", ".pdf"}
        }
        recorded = {row["path"] for row in rows}
        self.assertEqual(len(rows), 46)
        self.assertEqual(len(recorded), 46)
        self.assertEqual(recorded, actual)

    def test_required_reports_and_closed_conditional_gates(self) -> None:
        required = {
            "archive_audit.md",
            "code_lineage.md",
            "paper_code_concordance.md",
            "figure_and_manuscript_provenance.md",
            "reproduction.md",
            "blocker_resolution.md",
            "final_assessment.md",
        }
        self.assertTrue(required.issubset({path.name for path in ANALYSIS.iterdir()}))
        assessment = (ANALYSIS / "final_assessment.md").read_text(encoding="utf-8")
        self.assertIn("**`UNLIKELY`**", assessment)
        self.assertIn("**`STOP`**", assessment)
        self.assertFalse((ROOT / "src/reconstructed_2018").exists())


if __name__ == "__main__":
    unittest.main()
