"""Task 42 inputs and gates; no stationary solve occurs on import."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
from dataclasses import asdict
import hashlib
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
spec = importlib.util.spec_from_file_location('task40_parent_common',
    ROOT / 'analysis/40_ae4_equal_cation_routing/validation_common.py')
task40 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(task40)
parent = task40.parent
import numpy as np
from modern_full_model.model import WT
from modern_full_model.ae4_catalan2025 import CLASSES, with_mechanism_class
from modern_full_model.validation import sha256_object, CONSERVATION_RESIDUAL_TOLERANCES

HERE = ROOT / 'analysis/42_catalan_2025_ae4_mechanism_classes'
RESULTS = ROOT / 'results/42_catalan_2025_ae4_mechanism_classes'
CLASS_ID = os.environ.get('AE4_TASK42_CLASS', 'C0')
assert CLASS_ID in CLASSES
OUT = RESULTS / CLASS_ID
PREPARED_HEAD = 'a3c53ce07b875a3daf4e26ded0ce9f1d55a1d76f'
BRANCH = 'codex/task-42-catalan-2025-ae4-mechanism-classes'
MODEL_ID = 'CATALAN_2025_SOURCE_' + CLASS_ID
write_json = parent.write_json
parameter_hashes = parent.parameter_hashes


def verify_parent_sources():
    manifest = json.loads((HERE / 'parent_source_manifest.json').read_text())
    for entry in manifest['files']:
        data = (ROOT / entry['path']).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        assert actual == entry['git_blob_sha1'], entry['path']
    return manifest


def parent_models_and_state():
    verify_parent_sources()
    _, _, _, rest, stim, _ = task40.models_and_reference()
    frozen = json.loads((task40.OUT / 'wt_rest.json').read_text())
    inputs = json.loads((task40.OUT / 'frozen_inputs.json').read_text())
    assert frozen['admissible'] and frozen['status'] == 'PASS'
    assert frozen['state_sha256'] == inputs['frozen_rest_state_sha256'] == sha256_object(frozen['state_vector'])
    assert parameter_hashes(rest) == inputs['rest_parameter_hashes']
    assert parameter_hashes(stim) == inputs['stimulus_parameter_hashes']
    return rest, stim, np.array(frozen['state_vector'], dtype=float), frozen


def models_and_reference():
    old_rest, old_stim, y40, frozen = parent_models_and_state()
    return (frozen, old_rest, old_stim, with_mechanism_class(old_rest, CLASS_ID),
            with_mechanism_class(old_stim, CLASS_ID), y40)


def verify_frozen_sources():
    verify_parent_sources()
    manifest = json.loads((OUT / 'frozen_inputs.json').read_text())
    for entry in manifest['task42_execution_files']:
        assert hashlib.sha256((ROOT / entry['path']).read_bytes()).hexdigest() == entry['sha256'], entry['path']
    return manifest


def frozen_models_and_state():
    manifest = verify_frozen_sources()
    saved, _, _, rest, stim, _ = models_and_reference()
    frozen = json.loads((OUT / 'wt_rest.json').read_text())
    assert frozen['admissible'] and frozen['status'] == 'PASS'
    assert frozen['state_sha256'] == sha256_object(frozen['state_vector']) == manifest['frozen_rest_state_sha256']
    assert parameter_hashes(stim) == manifest['stimulus_parameter_hashes']
    return saved, rest, stim, np.array(frozen['state_vector'], dtype=float), frozen


def diagnose(model, t, y, e):
    # Task 40's generic gates come from Task 39. Its equal routing signature is
    # deliberately replaced with the exact predeclared class signature here.
    row, failures, ratios = parent.diagnose(model, t, y, e)
    a = e.diagnostics.ae4
    coefficients = CLASSES[CLASS_ID].conserved_coefficients
    fields = ('na_cell_fmol_s', 'k_cell_fmol_s', 'cl_cell_fmol_s',
              'tic_cell_fmol_s', 'alkalinity_cell_fmol_s')
    actual = np.array([getattr(a, f) for f in fields])
    expected = np.array(coefficients) * a.cl_cell_fmol_s
    if not np.array_equal(actual, expected):
        failures.append({'gate': 'predeclared_source_signature'})
    row.update(dict(zip(('ae4_na_cell_source_fmol_s', 'ae4_k_cell_source_fmol_s',
                         'ae4_cl_source_fmol_s', 'ae4_tic_cell_source_fmol_s',
                         'ae4_ta_cell_source_fmol_s'), map(float, actual))))
    row.update({'co3_i_mM': e.diagnostics.observables.cell_acid_base.co3_mM,
                'ae4_hco3_cell_source_fmol_s': float(a.hco3_cell_fmol_s),
                'ae4_co3_cell_source_fmol_s': float(a.co3_cell_fmol_s),
                'source_signature_residual_fmol_s': float(np.max(np.abs(actual-expected))),
                'ae4_charge_residual_fmol_s': float(a.charge_source_residual_fmol_s)})
    # Independent assembly of all five conserved cellular balances, using the
    # explicit TIC and TA contribution rather than a bicarbonate alias.
    bookkeeping = parent.inherited.bookkeeping(model, e)
    row['cell_balance_residual_fmol_s'] = bookkeeping['cell_rhs_assembly_max_abs_fmol_s']
    if row['cell_balance_residual_fmol_s'] > 1e-10:
        failures.append({'gate': 'independent_cell_balance'})
    if model.base_model.base_model.mechanism.name != CLASS_ID:
        failures.append({'gate': 'predeclared_class_active'})
    return row, failures, ratios
