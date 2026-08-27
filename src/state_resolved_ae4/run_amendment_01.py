"""Post-reveal carbonate-map correction and KO-blind finite gauge grid."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from .cycles import CANDIDATE_DEFINITIONS, ReservoirActivities, evaluate_cycle
from .fitting import evaluate_transporter_candidate, transporter_fit_grid
from .solvers import parameters_from_fit


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(repository_root: str | Path = ".") -> dict[str, object]:
    root = Path(repository_root).resolve()
    result_dir = root / "results/13_state_resolved_ae4"
    marker_path = result_dir / "stage_a_reveal_marker.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    summary_path = result_dir / "transporter_fit_summary.csv"
    original_hash = _hash(summary_path)
    declared_hash = marker["file_hashes"][
        "results/13_state_resolved_ae4/transporter_fit_summary.csv"
    ]
    if original_hash != declared_hash:
        raise RuntimeError("pre-reveal transporter summary no longer matches its marker")
    preserved = result_dir / "transporter_fit_summary_pre_reveal.csv"
    shutil.copyfile(summary_path, preserved)

    # Correction: BCECF alkalinization is mapped to total inorganic-carbon
    # alkalinity source HCO3 + 2 CO3, not HCO3 alone.  No knockout value is
    # used anywhere in the recalculation.
    fits = transporter_fit_grid()
    amended_rows = [row.as_dict() for row in fits]
    amended_path = result_dir / "transporter_fit_summary_amended.csv"
    _write_csv(amended_path, amended_rows)
    _write_csv(summary_path, amended_rows)

    definition = CANDIDATE_DEFINITIONS["SR5_COORDINATION_112"]
    baseline = ReservoirActivities(
        na_i=25.0,
        k_i=120.0,
        cl_i=50.0,
        hco3_i=10.0,
        na_e=140.2,
        k_e=5.3,
        cl_e=102.6,
        hco3_e=21.0,
    )
    scan_rows: list[dict[str, object]] = []
    for multiplier in np.logspace(-4.0, 4.0, 17):
        fit = evaluate_transporter_candidate(
            "SR5_COORDINATION_112",
            gauge_id=f"finite_grid_{multiplier:.8g}",
            common_barrier_multiplier=float(multiplier),
        )
        parameters = parameters_from_fit(fit, capacity=1.0)
        evaluation = evaluate_cycle(definition, baseline, parameters)
        sources = evaluation.intracellular_sources
        cl_source = float(sources["cl_i"])
        scan_rows.append(
            {
                "model_core": "SR5_COORDINATION_112; basal core inherited by SR6",
                "common_barrier_multiplier": float(multiplier),
                "transporter_summary_pass": fit.transporter_gate_pass,
                "na_ec50_pred_mM": fit.na_ec50_pred_mM,
                "na_hill_pred": fit.na_hill_pred,
                "k_ec50_pred_mM": fit.k_ec50_pred_mM,
                "k_hill_pred": fit.k_hill_pred,
                "na_source_per_capacity": float(sources["na_i"]),
                "k_source_per_capacity": float(sources["k_i"]),
                "cl_source_per_capacity": cl_source,
                "hco3_source_per_capacity": float(sources["hco3_i"]),
                "na_over_cl": float(sources["na_i"] / cl_source) if cl_source else np.nan,
                "k_over_cl": float(sources["k_i"] / cl_source) if cl_source else np.nan,
                "hco3_over_cl": float(sources["hco3_i"] / cl_source) if cl_source else np.nan,
                "cl_loading_direction_matches_fixed_wt_demand": cl_source > 0.0,
                "entropy_production": evaluation.entropy_production,
                "heldout_ko_data_read": False,
            }
        )
    gauge_path = result_dir / "barrier_gauge_source_scan.csv"
    _write_csv(gauge_path, scan_rows)

    amendment = {
        "amendment_id": "A01_post_reveal_carbonate_observable_and_finite_gauge_grid",
        "chronology": "performed after the original Stage-A reveal marker",
        "reason": (
            "correct SR4 BCECF projection from HCO3-only to total alkalinity "
            "HCO3+2CO3 and add a KO-blind 17-point log-spaced microscopic-barrier grid"
        ),
        "knockout_ionic_values_used": False,
        "knockout_secretion_values_used": False,
        "pre_reveal_summary_preserved_path": str(preserved.relative_to(root)),
        "pre_reveal_summary_hash": _hash(preserved),
        "original_marker_declared_hash": declared_hash,
        "amended_summary_path": str(amended_path.relative_to(root)),
        "amended_summary_hash": _hash(amended_path),
        "required_summary_current_hash": _hash(summary_path),
        "gauge_scan_hash": _hash(gauge_path),
        "survivors_before": 5,
        "survivors_after": sum(row.transporter_gate_pass for row in fits),
        "stage_a_outcome_changed": False,
        "interpretation": (
            "SR4a remains a chemistry alternative whose one-cation alkalinity curve "
            "misses the cooperative Hill slopes in this declared context; it is still "
            "not embeddable in the bicarbonate-only fixed chassis. The bounded null "
            "root and WT-gate failure are unchanged."
        ),
    }
    amendment_path = result_dir / "amendment_01.json"
    amendment_path.write_text(json.dumps(amendment, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return amendment


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
