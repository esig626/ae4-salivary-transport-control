"""Export frozen CTMC rates for Agent-H-style independent reproduction."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from .cycles import CANDIDATE_DEFINITIONS, ReservoirActivities, ctmc_snapshot
from .fitting import evaluate_transporter_candidate, nominal_assay_context
from .solvers import parameters_from_fit


SELECTED = (
    ("SR2_SHARED_123", "slow_common", 0.1),
    ("SR5_COORDINATION_112", "unit_common", 1.0),
    ("SR5_COORDINATION_112", "fast_common", 10.0),
    ("SR6_PKA_EDGE_112", "unit_common", 1.0),
    ("SR6_PKA_EDGE_112", "fast_common", 10.0),
)


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(repository_root: str | Path = ".") -> dict[str, object]:
    root = Path(repository_root).resolve()
    out = root / "results/13_state_resolved_ae4"
    condition_rows: list[dict[str, object]] = []
    state_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    physiology = ReservoirActivities(
        na_i=25.0,
        k_i=120.0,
        cl_i=50.0,
        hco3_i=10.0,
        na_e=140.2,
        k_e=5.3,
        cl_e=102.6,
        hco3_e=21.0,
    )
    for model_id, gauge_id, multiplier in SELECTED:
        fit = evaluate_transporter_candidate(
            model_id,
            gauge_id=gauge_id,
            common_barrier_multiplier=multiplier,
        )
        if not fit.transporter_gate_pass:
            raise RuntimeError(f"selected row no longer passes: {model_id}/{gauge_id}")
        parameters = parameters_from_fit(fit, capacity=1.0)
        contexts = [
            ("physiology_basal", physiology, 0.0),
            ("na_assay_150", nominal_assay_context("na", 150.0), 0.0),
            ("k_assay_150", nominal_assay_context("k", 150.0), 0.0),
        ]
        if CANDIDATE_DEFINITIONS[model_id].pka_switch_edge_ids:
            contexts.append(("na_assay_150_pka", nominal_assay_context("na", 150.0), 1.0))
        for condition_id, reservoirs, pka in contexts:
            snapshot = ctmc_snapshot(
                CANDIDATE_DEFINITIONS[model_id], reservoirs, parameters, pka_activation=pka
            )
            key = f"{model_id}|{gauge_id}|{condition_id}"
            condition_rows.append(
                {
                    "snapshot_id": key,
                    "model_id": model_id,
                    "gauge_id": gauge_id,
                    "condition_id": condition_id,
                    "pka_activation": pka,
                    "concentrations_mM_json": json.dumps(dict(snapshot["concentrations_mM"]), sort_keys=True),
                    "parameters_json": json.dumps(
                        {
                            "capacity": parameters.capacity,
                            "na_attempt_scale": parameters.na_attempt_scale,
                            "k_attempt_scale": parameters.k_attempt_scale,
                            "ec50_na_energy_lump_mM": parameters.ec50_na_mM,
                            "ec50_k_energy_lump_mM": parameters.ec50_k_mM,
                            "pka_barrier_fold": parameters.pka_barrier_fold,
                            "pka_basal_activity": parameters.pka_basal_activity,
                            "pka_phosphorylation_scale": parameters.pka_phosphorylation_scale,
                            "pka_dephosphorylation_scale": parameters.pka_dephosphorylation_scale,
                            "edge_barrier_overrides": dict(parameters.edge_barrier_overrides),
                        },
                        sort_keys=True,
                    ),
                    "state_count": len(snapshot["state_order"]),
                    "edge_count": len(snapshot["edge_rows"]),
                    "heldout_ko_data_read": False,
                }
            )
            for index, state in enumerate(snapshot["state_order"]):
                state_rows.append(
                    {
                        "snapshot_id": key,
                        "state_index": index,
                        "state_id": state,
                        "effective_state_energy": snapshot["state_energies"][state],
                        "agent_e_occupancy": snapshot["occupancy"][index],
                    }
                )
            for edge in snapshot["edge_rows"]:
                edge_rows.append(
                    {
                        "snapshot_id": key,
                        "edge_id": edge["edge_id"],
                        "source_state": edge["source"],
                        "target_state": edge["target"],
                        "forward_rate": edge["forward_rate"],
                        "reverse_rate": edge["reverse_rate"],
                        "agent_e_stationary_current": edge[
                            "stationary_current_per_transporter"
                        ],
                        "reservoir_stoichiometry_json": json.dumps(
                            edge["reservoir_stoichiometry"], sort_keys=True
                        ),
                    }
                )
    conditions_path = out / "ctmc_frozen_conditions.csv"
    states_path = out / "ctmc_frozen_states.csv"
    edges_path = out / "ctmc_frozen_edges.csv"
    _write(conditions_path, condition_rows)
    _write(states_path, state_rows)
    _write(edges_path, edge_rows)
    manifest = {
        "status": "frozen_for_independent_reproduction",
        "condition_count": len(condition_rows),
        "state_row_count": len(state_rows),
        "edge_row_count": len(edge_rows),
        "conditions_sha256": _hash(conditions_path),
        "states_sha256": _hash(states_path),
        "edges_sha256": _hash(edges_path),
        "strict_holdout_read": False,
        "column_semantics": (
            "row-oriented CTMC: each reversible edge gives source->target kf and "
            "target->source kr; construct Q with row sums zero and solve pi Q=0, sum pi=1"
        ),
    }
    manifest_path = out / "ctmc_snapshot_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
