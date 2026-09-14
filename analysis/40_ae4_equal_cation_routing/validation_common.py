"""Frozen Task 39 inputs with only the Task 40 AE4 evaluator selected."""
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
from dataclasses import asdict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("task39_validation_common",
    ROOT / "analysis/39_palk_nkcc1_full_validation/validation_common.py")
parent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)
import numpy as np
from modern_full_model.model import ModernFullModel, WT
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model
from modern_full_model.ae4_routing_only import routing_only_adapter
from modern_full_model.ae4_equal_cation_routing import equal_routing_adapter, MODEL_ID
from modern_full_model.validation import sha256_object, CONSERVATION_RESIDUAL_TOLERANCES

HERE = ROOT / "analysis/40_ae4_equal_cation_routing"
OUT = ROOT / "results/40_ae4_equal_cation_routing"
PREPARED_HEAD = "a6a58ebabcd2aedb1702b52dfdd86359c60cec32"
BRANCH = "codex/task-40-ae4-equal-cation-routing"
write_json = parent.write_json
parameter_hashes = parent.parameter_hashes


def select_equal(template):
    old = template.base_model.base_model
    assert old.ae4_evaluator is routing_only_adapter
    core = ModernFullModel(old.parameters, stimulus=old.stimulus,
        regulatory_model=old.regulatory_model, ae4_parameters=old.ae4_parameters,
        ae4_evaluator=equal_routing_adapter, nkcc1_kinetics=old.nkcc1_kinetics)
    return MinimalNbcModel(StimulatedNkcc1Model(core), template.nbc_parameters)


def models_and_reference():
    saved, _, _, _, old_rest, old_stim, y31 = parent.models_and_state()
    return saved, old_rest, old_stim, select_equal(old_rest), select_equal(old_stim), y31


def verify_parent_sources():
    manifest = json.loads((HERE / "parent_source_manifest.json").read_text())
    for entry in manifest["files"]:
        data = (ROOT / entry["path"]).read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        assert actual == entry["git_blob_sha1"], entry["path"]
    return manifest


def verify_frozen_sources():
    verify_parent_sources()
    manifest = json.loads((OUT / "frozen_inputs.json").read_text())
    for entry in manifest["task40_execution_files"]:
        assert hashlib.sha256((ROOT / entry["path"]).read_bytes()).hexdigest() == entry["sha256"], entry["path"]
    return manifest


def frozen_models_and_state():
    manifest = verify_frozen_sources()
    saved, _, _, rest, stim, _ = models_and_reference()
    frozen = json.loads((OUT / "wt_rest.json").read_text())
    assert frozen["admissible"] and frozen["status"] == "PASS"
    assert frozen["state_sha256"] == sha256_object(frozen["state_vector"])
    assert frozen["state_sha256"] == manifest["frozen_rest_state_sha256"]
    assert parameter_hashes(stim) == manifest["stimulus_parameter_hashes"]
    assert parameter_hashes(rest) == manifest["rest_parameter_hashes"]
    return saved, rest, stim, np.array(frozen["state_vector"], dtype=float), frozen


def diagnose(model, t, y, evaluation):
    row, failures, ratios = parent.diagnose(model, t, y, evaluation)
    a = evaluation.diagnostics.ae4
    j = a.cl_cell_fmol_s
    row.update({"ae4_na_cell_source_fmol_s": float(a.na_cell_fmol_s),
                "ae4_k_cell_source_fmol_s": float(a.k_cell_fmol_s),
                "ae4_na_branch_export_fmol_s": float(-a.na_cell_fmol_s),
                "ae4_k_branch_export_fmol_s": float(-a.k_cell_fmol_s),
                "ae4_tic_cell_source_fmol_s": float(a.tic_cell_fmol_s),
                "ae4_ta_cell_source_fmol_s": float(a.alkalinity_cell_fmol_s),
                "ae2_positive_cl_loading_fmol_s": float(max(0., row["ae2_cl_inward_signed_fmol_s"]))})
    if model.ae4_evaluator is not equal_routing_adapter:
        failures.append({"gate": "equal_routing_evaluator_active"})
    if not (a.na_cell_fmol_s == a.k_cell_fmol_s == -0.5*j
            and a.tic_cell_fmol_s == a.alkalinity_cell_fmol_s == -2*j):
        failures.append({"gate": "equal_routing_source_signature"})
    return row, failures, ratios
