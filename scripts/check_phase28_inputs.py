"""Check the small Task 28 seed without Git, a solver, or network access."""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def git_blob_sha(data: bytes) -> str:
    header = b"blob " + str(len(data)).encode("ascii") + b"\0"
    return hashlib.sha1(header + data).hexdigest()


def check_inputs(root: Path = ROOT, *, check_original_code: bool = False,
                 check_imports: bool = True) -> dict[str, Any]:
    root = root.resolve()
    manifest = json.loads((root / "reference/import_manifest.json").read_text(encoding="utf-8"))
    checked = 0
    for entry in manifest["files"]:
        path = root / entry["path"]
        if not path.is_file():
            raise FileNotFoundError(f"Missing seed file: {entry['path']}")
        if entry["kind"] == "reference" or check_original_code:
            actual = git_blob_sha(path.read_bytes())
            if actual != entry["git_blob_sha"]:
                raise ValueError(f"Seed identity mismatch: {entry['path']}")
            checked += 1
    native = json.loads((root / "reference/native_wt_contract.json").read_text(encoding="utf-8"))
    roots = native.get("roots", {})
    if len(roots) != 10:
        raise ValueError(f"Expected ten retained parameter backgrounds, found {len(roots)}")
    for root_id, payload in roots.items():
        for key in ("whole_cell_parameters", "ae4_parameters"):
            if not isinstance(payload.get(key), dict):
                raise ValueError(f"Missing {key} for {root_id}")
    selected = {}
    for alias, family in (("R09", "AE4NA05"), ("R10", "AE4NA20")):
        state = json.loads((root / f"reference/{alias}_wt_rest.json").read_text(encoding="utf-8"))
        root_id = state["root_id"]
        if root_id not in roots or state["root_alias"] != alias:
            raise ValueError(f"Incorrect reference mapping for {alias}")
        if state["routing_family"] != family or state["ae4_expression"] != 1.0:
            raise ValueError(f"Incorrect WT reference for {alias}")
        if len(state["core_state"]) != 12 or not state["rest_gate_pass"]:
            raise ValueError(f"Incomplete saved resting state for {alias}")
        for key in ("whole_cell_parameters_sha256", "ae4_parameters_sha256"):
            if state[key] != roots[root_id][key]:
                raise ValueError(f"Parameter provenance mismatch for {alias}: {key}")
        selected[alias] = root_id
    contract = json.loads((root / "reference/routing_numerical_contract.json").read_text(encoding="utf-8"))
    if contract["model_id"] != "SR2_NET_CHLORIDE_DONOR_CATION_ROUTING":
        raise ValueError("Incorrect reference routing law")
    if contract["resting_calcium_uM"] != 0.058:
        raise ValueError("Unexpected resting calcium")
    imported = []
    if check_imports:
        sys.path.insert(0, str(root / "src"))
        for entry in manifest["files"]:
            if entry["kind"] == "code":
                module_name = entry["path"][4:-3].replace("/", ".")
                module = importlib.import_module(module_name)
                if Path(module.__file__).resolve() != (root / entry["path"]).resolve():
                    raise ValueError(f"Import resolved outside this checkout: {module_name}")
                imported.append(module_name)
    return {"status": "SEED_INPUT_CHECK_PASSED", "files_hash_checked": checked,
            "reference_parameter_backgrounds": len(roots), "selected_backgrounds": selected,
            "modules_imported": imported, "scientific_replay_performed": False,
            "stationary_solves": 0, "stimulated_trajectories": 0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-original-code", action="store_true",
                        help="Verify initial scientific code identities before any edits")
    args = parser.parse_args()
    try:
        print(json.dumps(check_inputs(check_original_code=args.check_original_code), indent=2))
    except (OSError, ValueError, KeyError, ImportError) as exc:
        raise SystemExit(f"Task 28 seed check failed: {exc}") from exc


if __name__ == "__main__":
    main()
