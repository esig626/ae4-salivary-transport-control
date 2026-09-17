"""Authorised output-schema repair; every 52C scientific file remains unchanged.

Only serialization and branch/replay bookkeeping are adapted in memory.
The original failed attempt is preserved under output/; this attempt has its
own output/reporter_recovery_52D/ namespace and cannot overwrite it.
"""
from pathlib import Path
import difflib
import hashlib
import inspect
import json
import subprocess
import sys

TASK = Path(__file__).resolve().parent
ROOT = TASK.parents[1]
SOURCE_OUT = TASK/'output'
RECOVERY_OUT = SOURCE_OUT/'reporter_recovery_52D'
BRANCH = 'analysis/task-52D-recovery-from-52C'
AUTHORISATION = 'f2b56f903d33a0ac5f5c292db92cead8b6157f1f'
SCIENTIFIC_BASE = '36f7c3d0bb2353bf55ac0ed9e22582720097c829'
FAILURE_COMMIT = '9c0fde42b5eb16746413801fae555c84cb32e201'
FREEZE_SHA = '9cebb8fcf85557643f6c7838ff57f9fb9603f5a7479cbe6fb63ba8c40b417ba2'
RUNNER_SHA = 'd5988680b99638fa91ed0e6e82a70f6b1d6bca581ee102438060d6e1025909ce'

# This import sets one-thread environment variables before numerical imports.
# Importing the frozen runner constructs/evaluates no scientific model.
import run_frozen_cases as frozen
np = frozen.np


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def encode_observations(samples, times):
    """Lossless name-based union encoding; missing has no physical value."""
    rows = [row for t in times for row in samples[t]['rows']]
    columns = tuple(dict.fromkeys(key for row in rows for key in row))
    present = np.array([[[key in row for key in columns]
                         for row in samples[t]['rows']] for t in times], dtype=bool)
    values = np.array([[[row.get(key, np.nan) for key in columns]
                        for row in samples[t]['rows']] for t in times], dtype=float)
    assert np.isnan(values[~present]).all()
    # Fields consumed by frozen quadrature summaries/extrema are unconditional
    # diagnostics. Missing optional current fields must never enter them.
    required = ('q_out_pL_s', 'J_fmol_s', 'N_fmol_s', 'A_fmol_s', 'E_fmol_s',
                'ph_i', 'cl_i_mM', 'tic_i_mM', 'na_i_mM', 'k_i_mM',
                'volume_i_pL', 'E_Cl_V', 'V_a_V', 'chloride_driving_force_V',
                'J_aux_fmol_s', 'B_fmol_s', 'H_fmol_s')
    for key in required:
        if any(key in row for row in rows):
            assert key in columns and present[:, :, columns.index(key)].all(), key
    return columns, values, present


def patched_attempt_source():
    original = inspect.getsource(frozen.run_attempt)
    edits = [
        ("    columns=tuple(samples[times[0]]['rows'][0])\n"
         "    assert all(tuple(r)==columns for s in samples.values() for r in s['rows'])\n"
         "    values=np.array([[list(r.values()) for r in samples[t]['rows']] for t in times])",
         "    columns,values,observations_present=encode_observations(samples,times)"),
        ("columns=np.array(columns),observations=values,paired_states=states,sample_kind=kinds,",
         "columns=np.array(columns),observations=values,observations_present=observations_present,\n"
         "        missing_diagnostic_encoding=np.array('NaN with observations_present=False'),\n"
         "        paired_states=states,sample_kind=kinds,"),
        ("endpoint={g:{k:float(values[-1,i,index[k]]) for k in columns} for i,g in enumerate(GENOTYPES)}",
         "endpoint={g:{k:(float(values[-1,i,index[k]]) if observations_present[-1,i,index[k]] else None)\n"
         "        for k in columns} for i,g in enumerate(GENOTYPES)}"),
        ("'finite_trajectory_arrays':bool(np.isfinite(values).all() and np.isfinite(states).all()),",
         "'finite_trajectory_arrays':bool(np.isfinite(values[observations_present]).all() and\n"
         "            np.isnan(values[~observations_present]).all() and np.isfinite(states).all()),\n"
         "        'diagnostic_schema':{'encoding':'union of original names; NaN+presence mask in NPZ; null in JSON',\n"
         "            'column_count':len(columns),'missing_value_count':int((~observations_present).sum()),\n"
         "            'original_named_values_preserved':True},"),
    ]
    adapted = original
    for old, new in edits:
        assert adapted.count(old) == 1, old
        adapted = adapted.replace(old, new)
    # All model/solver/observe/physical-gate calls precede this exact boundary.
    boundary = '    columns='
    original_prefix = original[:original.index(boundary)]
    assert adapted.startswith(original_prefix)
    return original, adapted


def verify_recovery_inputs():
    assert git('branch', '--show-current').decode().strip() == BRANCH
    assert git('rev-parse', 'HEAD').decode().strip() == AUTHORISATION
    assert git('merge-base', SCIENTIFIC_BASE, AUTHORISATION).decode().strip() == SCIENTIFIC_BASE
    assert sha((TASK/'run_frozen_cases.py').read_bytes()) == RUNNER_SHA
    assert sha((SOURCE_OUT/'parameter_and_case_freeze.json').read_bytes()) == FREEZE_SHA
    authpath = str((TASK/'AUTHORIZATION_52D_REPORTER_RECOVERY.md').relative_to(ROOT))
    assert (ROOT/authpath).read_bytes() == git('show', AUTHORISATION+':'+authpath)
    # Exact original failure evidence stays untouched; status/ledger are appendable.
    failures = ('execution_52D.json', 'execution_52D.log', 'failure_52D.json',
                'cases/case_01/started.json', 'launch_52D_recovery.py')
    for name in failures:
        p = SOURCE_OUT/name
        assert p.read_bytes() == git('show', FAILURE_COMMIT+':'+str(p.relative_to(ROOT)))
    # Retain all frozen integrity gates, changing only actual branch and
    # authorisation-HEAD versus immutable scientific-parent bookkeeping.
    source = inspect.getsource(frozen.verify_inputs)
    old_branch = "'analysis/task-52-chloride-reservoir-final-test'"
    old_parent = "publication['commit_sha'] == git('rev-parse','HEAD')"
    assert source.count(old_branch) == source.count(old_parent) == 1
    source = source.replace(old_branch, repr(BRANCH)).replace(old_parent,
                           "publication['commit_sha'] == SCIENTIFIC_BASE")
    namespace = dict(frozen.__dict__, OUT=SOURCE_OUT, SCIENTIFIC_BASE=SCIENTIFIC_BASE)
    exec(compile(source, str(__file__)+':input_bookkeeping', 'exec'), namespace)
    freeze = namespace['verify_inputs']('52D')
    selected = [c for c in freeze['case_matrix']['cases'] if c['checkpoint']=='52D']
    assert [c['id'] for c in selected] == ['case_01','case_02','case_03']
    assert [c['G_aux_S'] for c in selected] == [0.0,2.32e-9,4.49e-9]
    assert all(c['projection']=='central' and c['protocol']=='CCH_IPR'
               and c['ko_supply_multiplier']==1.0 for c in selected)
    assert all(v==0 for v in frozen.COUNTS.values())
    return freeze


def main():
    assert sys.argv[1:]==[], 'Only the authorised52D stage is available'
    freeze = verify_recovery_inputs()
    original, adapted = patched_attempt_source()
    RECOVERY_OUT.mkdir(exist_ok=False)
    (RECOVERY_OUT/'reporting_patch.diff').write_text(''.join(difflib.unified_diff(
        original.splitlines(True), adapted.splitlines(True),
        fromfile='52C/run_attempt', tofile='authorised_serialisation/run_attempt')))
    frozen.dump(RECOVERY_OUT/'recovery_receipt.json', {
        'authorisation_commit':AUTHORISATION, 'scientific_base_commit':SCIENTIFIC_BASE,
        'preserved_failure_commit':FAILURE_COMMIT, 'freeze_sha256':FREEZE_SHA,
        'runner_sha256':RUNNER_SHA, 'launcher_sha256':sha(Path(__file__).read_bytes()),
        'immutable_inputs_verified':len(freeze['immutable_artifact_sha256']),
        'scientific_evaluations_before_execution':dict(frozen.COUNTS),
        'changes':'union-schema serialization and authorised branch/parent/replay bookkeeping only',
        'original_failure_evidence_unchanged':True})
    frozen.__dict__['encode_observations'] = encode_observations
    exec(compile(adapted, str(__file__)+':serialisation', 'exec'), frozen.__dict__)
    frozen.OUT = RECOVERY_OUT
    # The frozen main records the input-freeze hash using OUT. Outputs move,
    # but immutable inputs remain in SOURCE_OUT; adapt only that record path.
    original_main = inspect.getsource(frozen.main)
    old_hash_path = "sha(OUT/'parameter_and_case_freeze.json')"
    assert original_main.count(old_hash_path) == 1
    adapted_main = original_main.replace(old_hash_path,
                                        "sha(SOURCE_OUT/'parameter_and_case_freeze.json')")
    frozen.__dict__['SOURCE_OUT'] = SOURCE_OUT
    (RECOVERY_OUT/'stage_bookkeeping_patch.diff').write_text(''.join(difflib.unified_diff(
        original_main.splitlines(True), adapted_main.splitlines(True),
        fromfile='52C/main', tofile='authorised_output_namespace/main')))
    exec(compile(adapted_main, str(__file__)+':output_bookkeeping', 'exec'), frozen.__dict__)
    def verified_52D_only(stage):
        assert stage=='52D'
        return freeze
    frozen.verify_inputs = verified_52D_only
    sys.argv = [str(TASK/'run_frozen_cases.py'), '52D']
    frozen.main()


if __name__=='__main__':
    main()
