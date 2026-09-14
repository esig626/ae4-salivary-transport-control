"""Task 39 fixed inputs, isolated model selection, and inherited diagnostics."""
from pathlib import Path
from dataclasses import asdict
import hashlib
import importlib.util
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("task37_validation_common",
    ROOT / "analysis/37_wt_nbc_validation/validation_common.py")
inherited = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inherited)
import numpy as np
from modern_full_model.model import ModernFullModel, WT
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model
from modern_full_model.nkcc1_palk2010 import Nkcc1Kinetics, PALK_NKCC1, intracellular_product_mM4
from modern_full_model.validation import sha256_object, CONSERVATION_RESIDUAL_TOLERANCES

HERE = ROOT / "analysis/39_palk_nkcc1_full_validation"
OUT = ROOT / "results/39_palk_nkcc1_full_validation"
PREPARED_HEAD = "4aa1a1454dfb042d55523d103f10f27d3de77132"
BRANCH = "codex/task-39-palk-nkcc1-full-validation"
TASK38_COMMIT = "d4bce14c5caed4e05587e237f002814f0f27e226"
REST_CONCENTRATIONS = (11.636125748680639, 116.76320932871538, 60.30496692587399)
REST_CYCLE = 0.1281222023866353
Q_REST = inherited.Q_REST
TASK37_CUMULATIVE = 0.986995055625111
write_json = inherited.write_json


def resting_algebra():
    """Literal Eq. 17, independent of the new evaluator; no numerical solve."""
    na, k, cl = REST_CONCENTRATIONS
    x = na * k * cl**2
    shape = (157.5 - 2.0096e-5 * x) / (1.0306 + 1.3852e-6 * x)
    return {"x_rest_mM4": x, "shape_rest": shape,
            "alpha_eff_fmol_s": REST_CYCLE / shape,
            "rest_cycle_fmol_s": REST_CYCLE}


def select_palk(template, alpha):
    """Keep the identical Task 38 architecture and select only the NKCC1 core."""
    old = template.base_model.base_model
    core = ModernFullModel(old.parameters, stimulus=old.stimulus,
        regulatory_model=old.regulatory_model, ae4_parameters=old.ae4_parameters,
        ae4_evaluator=old.ae4_evaluator,
        nkcc1_kinetics=Nkcc1Kinetics(PALK_NKCC1, alpha))
    return MinimalNbcModel(StimulatedNkcc1Model(core), template.nbc_parameters)


def parameter_hashes(model):
    return {"whole_cell": sha256_object(model.parameters),
            "ae4": sha256_object(model.ae4_parameters),
            "nbc": sha256_object(model.nbc_parameters),
            "regulation": sha256_object(model.regulatory_model),
            "stimulus": sha256_object(model.stimulus),
            "nkcc1_core": sha256_object(model.nkcc1_kinetics)}


def models_and_state():
    freeze = json.loads((OUT / "alpha_eff_freeze.json").read_text())
    assert freeze["sha256"] == sha256_object(freeze["payload"])
    saved, task31, old_rest, old_stim, y0 = inherited.models_and_state()
    alpha = freeze["payload"]["alpha_eff_fmol_s"]
    return saved, task31, old_rest, old_stim, select_palk(old_rest, alpha), select_palk(old_stim, alpha), y0


def diagnose(model, t, y, e):
    row, failures, ratios = inherited.diagnose(model, t, y, e)
    h, m, r = e.diagnostics.homeostasis, e.diagnostics.membranes, e.diagnostics.regulatory
    row.update({
        "nkcc1_cycle_inward_fmol_s": float(h.nkcc1_inward_fmol_s),
        "nkcc1_activity_multiplier": float(r["nkcc1_capacity_multiplier"]),
        "nkcc1_x_i_mM4": intracellular_product_mM4(row["na_i_mM"], row["k_i_mM"], row["cl_i_mM"]),
        "nbc_hco3_inward_fmol_s": float(2 * r["minimal_nbc_cycle_inward_fmol_s"]),
        "pump_total_cycles_fmol_s": float(m.pump_apical_fmol_s + m.pump_basolateral_fmol_s),
    })
    return row, failures, ratios


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_frozen_sources():
    manifest = json.loads((OUT / "frozen_inputs.json").read_text())
    for item in manifest["execution_files"]:
        assert file_sha256(ROOT / item["path"]) == item["sha256"], item["path"]
    return manifest
