"""Build the final Task 13B pre-reveal G5 ensemble manifest and gate.

The freeze retains all 18 WT-passing root/profile groups.  It does not select
one calcium/NKCC1 pair or one R0--R3 kinetic member, because the available WT
evidence does not identify those axes uniquely.  Every dependency is hashed;
the held-out target ledger is neither opened nor hashed by this module.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, is_dataclass
from enum import Enum
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from .run_g5_cross_root_profile import load_retained_wt_roots
from .run_g5_joint_repair_profile import declared_joint_repair_profiles
from .validation import (
    LOOSE_RADAU,
    PRODUCTION_BDF,
    PRODUCTION_RADAU,
    TIGHT_RADAU,
    SecretagogueProtocol,
    physical_time_grid,
    pre_reveal_regulatory_ensemble,
    sha256_file,
    sha256_object,
)


MANIFEST_FILENAME = "g5_final_pre_reveal_manifest.json"
GATE_FILENAME = "g5_final_pre_reveal_gate.json"
HASH_FILENAME = "g5_final_pre_reveal_hashes.csv"


def _jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "__dict__"):
        return {
            str(key): _jsonable(item)
            for key, item in vars(value).items()
            if not str(key).startswith("_")
        }
    return repr(value)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path) -> Mapping[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for name in row:
            if name not in fields:
                fields.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _typed_passing_group(row: Mapping[str, str]) -> Mapping[str, Any]:
    return {
        "root_id": row["root_id"],
        "selected_reference": row["selected_reference"] == "True",
        "profile_id": row["profile_id"],
        "effective_calcium_uM": float(row["effective_calcium_uM"]),
        "nkcc1_fully_activated_multiplier": float(
            row["nkcc1_fully_activated_multiplier"]
        ),
        "source_defensibility_status": row["source_defensibility_status"],
        "regulatory_member_count": int(row["regulatory_member_count"]),
        "solver_count": int(row["solver_count"]),
        "all_integrations_success": row["all_integrations_success"] == "True",
        "all_numerical_gates_pass": row["all_numerical_gates_pass"] == "True",
        "all_individual_shape_gates_pass": (
            row["all_individual_shape_gates_pass"] == "True"
        ),
        "common_scale_feasible_intersection": (
            row["common_scale_feasible_intersection"] == "True"
        ),
        "common_scale_feasible_lower_uL_min_per_pL_s": float(
            row["common_scale_feasible_lower_uL_min_per_pL_s"]
        ),
        "common_scale_feasible_upper_uL_min_per_pL_s": float(
            row["common_scale_feasible_upper_uL_min_per_pL_s"]
        ),
        "common_scale_multiplier_uL_min_per_pL_s": float(
            row["common_scale_multiplier_uL_min_per_pL_s"]
        ),
        "common_scale_implied_effective_cell_count": float(
            row["common_scale_implied_effective_cell_count"]
        ),
        "implied_cell_count_status": row["implied_cell_count_status"],
        "ensemble_min_ratio": float(row["ensemble_min_ratio"]),
        "ensemble_max_ratio": float(row["ensemble_max_ratio"]),
        "ensemble_endpoint_cl_range_mM": json.loads(
            row["ensemble_endpoint_cl_range_mM"]
        ),
        "ensemble_endpoint_ph_range": json.loads(
            row["ensemble_endpoint_ph_range"]
        ),
        "ensemble_endpoint_volume_range_pL": json.loads(
            row["ensemble_endpoint_volume_range_pL"]
        ),
        "max_abs_radau_bdf_ratio_difference": float(
            row["max_abs_radau_bdf_ratio_difference"]
        ),
        "g5_joint_group_gate_pass": row["g5_joint_group_gate_pass"] == "True",
    }


def build_final_g5_freeze(
    repository_root: Path,
    results_directory: Path,
) -> Mapping[str, Any]:
    group_path = results_directory / "wt_g5_joint_ca_nkcc_group_gate.csv"
    all_group_rows = _read_rows(group_path)
    passing = [
        _typed_passing_group(row)
        for row in all_group_rows
        if row["g5_joint_group_gate_pass"] == "True"
    ]
    failed = [
        {
            "root_id": row["root_id"],
            "profile_id": row["profile_id"],
            "all_numerical_gates_pass": row["all_numerical_gates_pass"] == "True",
            "all_individual_shape_gates_pass": (
                row["all_individual_shape_gates_pass"] == "True"
            ),
            "common_scale_feasible_intersection": (
                row["common_scale_feasible_intersection"] == "True"
            ),
            "ensemble_max_ratio": float(row["ensemble_max_ratio"]),
        }
        for row in all_group_rows
        if row["g5_joint_group_gate_pass"] != "True"
    ]
    if len(passing) != 18 or len(failed) != 9:
        raise AssertionError("final G5 freeze requires exactly 18 passing and 9 failed groups")
    if not all(
        group["all_integrations_success"]
        and group["all_numerical_gates_pass"]
        and group["all_individual_shape_gates_pass"]
        and group["common_scale_feasible_intersection"]
        and group["g5_joint_group_gate_pass"]
        for group in passing
    ):
        raise AssertionError("a retained final G5 group does not pass every WT gate")

    ae2_summary_path = results_directory / "wt_g5_ae2_independent_summary.json"
    nearby_summary_path = results_directory / "wt_g5_final_candidate_nearby_summary.json"
    ae2_summary = _read_json(ae2_summary_path)
    nearby_summary = _read_json(nearby_summary_path)
    if not bool(ae2_summary["all_paired_numerical_gates_pass"]):
        raise AssertionError("AE2 independent numerical validation did not pass")
    if not bool(nearby_summary["all_candidate_nearby_gates_pass"]):
        raise AssertionError("selected-root cell/lumen H23 check did not pass")

    roots = load_retained_wt_roots(results_directory)
    root_payloads = {
        root.root_id: {
            "selected_reference": root.selected_reference,
            "variant": _jsonable(root.variant),
            "calibration": _jsonable(root.calibration),
            "initial_core_state": list(root.core_state),
            "initial_core_state_sha256": sha256_object(root.core_state),
        }
        for root in roots
    }
    final_profile_ids = sorted({group["profile_id"] for group in passing})
    profile_registry = {
        profile.profile_id: profile for profile in declared_joint_repair_profiles()
    }
    profile_payloads = {
        profile_id: {
            **_jsonable(profile_registry[profile_id]),
            "stimulus_protocol": _jsonable(
                SecretagogueProtocol(
                    stimulated_calcium_uM=(
                        profile_registry[profile_id].effective_calcium_uM
                    )
                )
            ),
            "stimulus_interpretation": (
                "WT_CALIBRATED_EFFECTIVE_CHANNEL_DRIVING_CA_COORDINATE; "
                "NOT_MEASURED_FREE_INTRACELLULAR_CA"
            ),
            "nkcc1_law": "N1_CCH_ALGEBRAIC_CAPACITY_ONLY",
            "nkcc1_normalization": {
                "resting_calcium_uM": 0.058,
                "stimulated_calcium_uM": (
                    profile_registry[profile_id].effective_calcium_uM
                ),
                "capacity_multiplier_at_stimulated_endpoint": (
                    profile_registry[profile_id].nkcc1_fully_activated_multiplier
                ),
            },
        }
        for profile_id in final_profile_ids
    }
    regulatory_payloads = {
        member.member_id: {
            **_jsonable(member.metadata()),
            "state_names": list(member.regulatory_model.state_names),
            "model_parameters": _jsonable(member.regulatory_model),
        }
        for member in pre_reveal_regulatory_ensemble()
    }
    grid = physical_time_grid(duration_s=600.0, sample_step_s=5.0)

    relative_source_files = (
        "src/modern_full_model/acid_base.py",
        "src/modern_full_model/calibration.py",
        "src/modern_full_model/camp_pka.py",
        "src/modern_full_model/freeze_g5_joint.py",
        "src/modern_full_model/membranes.py",
        "src/modern_full_model/model.py",
        "src/modern_full_model/nkcc_stimulation.py",
        "src/modern_full_model/parameters.py",
        "src/modern_full_model/run_analysis.py",
        "src/modern_full_model/run_g5_ae2_validation.py",
        "src/modern_full_model/run_g5_cross_root_profile.py",
        "src/modern_full_model/run_g5_joint_repair_profile.py",
        "src/modern_full_model/run_g5_nearby_robustness.py",
        "src/modern_full_model/states.py",
        "src/modern_full_model/transporters.py",
        "src/modern_full_model/validation.py",
        "src/modern_full_model/water.py",
    )
    source_hashes = {
        name: sha256_file(repository_root / name) for name in relative_source_files
    }
    source_tree_hash = sha256_object(
        [{"path": name, "sha256": source_hashes[name]} for name in sorted(source_hashes)]
    )

    dependency_artifacts = (
        "wt_selected_candidate.json",
        "wt_root_table.csv",
        "wt_root_states.csv",
        "wt_parameter_candidates.csv",
        "dynamic_frozen_manifest.json",
        "wt_dynamic_gate.json",
        "regulatory_sustainment_summary.json",
        "wt_nkcc_stimulation_summary.json",
        "wt_g5_cross_root_repair_summary.json",
        "wt_g5_joint_ca_nkcc_regulatory_profile.csv",
        "wt_g5_joint_ca_nkcc_group_gate.csv",
        "wt_g5_joint_ca_nkcc_summary.json",
        "wt_g5_ae2_independent_validation.csv",
        "wt_g5_ae2_independent_summary.json",
        "wt_g5_final_candidate_nearby_robustness.csv",
        "wt_g5_final_candidate_nearby_summary.json",
    )
    dependency_hashes = {
        name: sha256_file(results_directory / name) for name in dependency_artifacts
    }
    selected_root_id = next(root.root_id for root in roots if root.selected_reference)
    selected_root_groups = [
        group for group in passing if group["root_id"] == selected_root_id
    ]

    gate = {
        "status": "PASS_PRE_REVEAL_ENSEMBLE_NONUNIQUE",
        "classification": (
            "PHYSIOLOGICALLY_CREDIBLE_WT_G5_ENSEMBLE_WITH_UNIDENTIFIED_"
            "EFFECTIVE_CA_NKCC_GAIN_PAIR"
        ),
        "pre_reveal_holdout_eligible": True,
        "passing_group_count": len(passing),
        "passing_root_count": len({group["root_id"] for group in passing}),
        "passing_profile_ids": final_profile_ids,
        "selected_reference_root": selected_root_id,
        "selected_reference_passing_groups": selected_root_groups,
        "regulatory_member_count": len(regulatory_payloads),
        "regulatory_families": ["R0", "R1", "R2", "R3"],
        "regulatory_gains": [1.25, 1.70],
        "regulatory_timing_profiles": ["STATIC", "TFAST", "TREFERENCE", "TSLOW"],
        "regulatory_kinetics_identified": False,
        "unidentified_regulatory_kinetics_is_reveal_blocker": False,
        "all_matrix_integrations_success": True,
        "all_retained_groups_numerical_gate_pass": True,
        "all_retained_groups_fixed_flow_shape_gate_pass": True,
        "all_retained_groups_common_scale_intersection_pass": True,
        "radau_bdf_gate_pass": True,
        "ae2_independent_gate_pass": True,
        "ae2_gate_interpretation": (
            "QUALITATIVELY_CONSISTENT_WITH_NO_DETECTED_LARGE_SECRETION_ROLE; "
            "P15_REPORTED_NO_EQUIVALENCE_MARGIN"
        ),
        "selected_root_cell_lumen_nearby_gate_pass": True,
        "effective_cell_count_status": (
            "WARNING_ONLY_NO_PRIMARY_GEOMETRY_RANGE"
        ),
        "effective_calcium_status": (
            "WT_CALIBRATED_EFFECTIVE_COORDINATE_NOT_MEASURED_FREE_CA"
        ),
        "candidate_uniqueness": False,
        "retained_nonuniqueness_reason": (
            "CA010_M1P75_AND_CA012_M2P0_BOTH_PASS_ALL_AVAILABLE_WT_GATES; "
            "AVAILABLE_CA_AND_ASSAY_CONTEXT_EVIDENCE_DOES_NOT_SELECT_ONE"
        ),
        "firewall_status": "PRE_REVEAL_FREEZE_COMPLETE_WITHOUT_TARGET_LEDGER_ACCESS",
    }

    manifest = {
        "manifest_type": "TASK_13B_FINAL_PRE_REVEAL_G5_ENSEMBLE_FREEZE",
        "status": gate["status"],
        "gate": gate,
        "equation_source_file_sha256": source_hashes,
        "equation_source_tree_sha256": source_tree_hash,
        "roots": root_payloads,
        "passing_groups": passing,
        "failed_groups_preserved": failed,
        "profiles": profile_payloads,
        "regulatory_models": regulatory_payloads,
        "physical_time": {
            "unit": "s",
            "duration_s": 600.0,
            "stimulated_right_limit_s": 1.0e-6,
            "output_sample_step_s": 5.0,
            "time_grid_sha256": sha256_object(grid),
            "production_max_step_s": 2.0,
        },
        "solver_specifications": {
            item.label: _jsonable(item)
            for item in (
                PRODUCTION_RADAU,
                PRODUCTION_BDF,
                LOOSE_RADAU,
                TIGHT_RADAU,
            )
        },
        "flow_scale": {
            "tier": "TIER_2_WT_OBSERVATION_SCALE",
            "source_envelope_uL_min": [9.0, 10.0],
            "rule": (
                "INTERSECT_ALL_20_MEMBER_AND_BOTH_SOLVER_MULTIPLIER_INTERVALS; "
                "USE_PREFERRED_R1_MEAN_SCALE_IF_INSIDE_ELSE_INTERSECTION_MIDPOINT"
            ),
            "frozen_scale_by_root_profile": {
                f"{group['root_id']}::{group['profile_id']}": {
                    "multiplier_uL_min_per_pL_s": group[
                        "common_scale_multiplier_uL_min_per_pL_s"
                    ],
                    "feasible_lower_uL_min_per_pL_s": group[
                        "common_scale_feasible_lower_uL_min_per_pL_s"
                    ],
                    "feasible_upper_uL_min_per_pL_s": group[
                        "common_scale_feasible_upper_uL_min_per_pL_s"
                    ],
                    "implied_effective_cell_count": group[
                        "common_scale_implied_effective_cell_count"
                    ],
                }
                for group in passing
            },
        },
        "independent_validation": {
            "ae2_summary": ae2_summary,
            "selected_root_nearby_summary": nearby_summary,
        },
        "supersession": {
            "baseline_g5": (
                "FAILED_FIXED_0P55_CA_STIMULUS_FLOW_SHAPE; PRESERVED_IN_"
                "dynamic_frozen_manifest.json"
            ),
            "regulatory_only_sustainment": (
                "FAILED_TO_REPAIR_FLOW_SHAPE; PRESERVED_IN_"
                "regulatory_sustainment_summary.json"
            ),
            "historical_ca_nkcc_profiles": (
                "NO_NOMINAL_0P55_CA_NKCC_PROFILE_PASSED; OUT_OF_SCALE_PROBES_"
                "PRESERVED"
            ),
            "ca015_n1_m1p75": (
                "FAILED_9_OF_9_ROOT_GROUPS; PRESERVED_IN_failed_groups_preserved"
            ),
            "replacement": (
                "JOINT_WT_CALIBRATED_EFFECTIVE_CA_PLUS_SOURCE_MOTIVATED_"
                "ALGEBRAIC_NKCC1_CAPACITY"
            ),
        },
        "dependency_artifact_sha256": dependency_hashes,
        "firewall": (
            "NO_HELDOUT_TARGET_FILE_OPEN_HASH_PARSE_OR_MODEL_SELECTION"
        ),
    }
    gate_path = results_directory / GATE_FILENAME
    manifest_path = results_directory / MANIFEST_FILENAME
    _write_json(gate_path, gate)
    _write_json(manifest_path, manifest)
    hash_rows = [
        {
            "artifact": name,
            "sha256": sha256_file(results_directory / name),
            "scope": "FINAL_PRE_REVEAL_G5",
        }
        for name in (*dependency_artifacts, GATE_FILENAME, MANIFEST_FILENAME)
    ]
    hash_path = results_directory / HASH_FILENAME
    _write_rows(hash_path, hash_rows)
    return {
        "gate": gate,
        "manifest_sha256": sha256_file(manifest_path),
        "gate_sha256": sha256_file(gate_path),
        "hash_table_sha256": sha256_file(hash_path),
        "equation_source_tree_sha256": source_tree_hash,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path("."))
    parser.add_argument(
        "--results-directory",
        type=Path,
        default=Path("results/13B_modern_full_model"),
    )
    args = parser.parse_args(argv)
    result = build_final_g5_freeze(
        args.repository_root.resolve(), args.results_directory.resolve()
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ("build_final_g5_freeze",)
