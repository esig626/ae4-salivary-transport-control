"""Task 20 post-freeze comparison for the certified empty WT feasible set.

No simulation is launched from an infeasible reference state. Legacy outputs
are reused after Git/SHA256 validation; all unavailable pooled outcomes remain
explicit. This report refuses to run if a feasible WT candidate is present.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

from . import task20_fixed_wt as task
from . import task18_fixed_wt as inverse
from .pooled_ae4 import MODEL_ID
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .validation import sha256_file


CHECKPOINT = 'a2902286c148c886ee2d8188dabbd020d0a2d32f'
MANIFEST_BLOB = '8aba7f2a2c884207951c2d88dcfae1ed26bd28ea'
UNAVAILABLE = 'NOT_RUN_STRUCTURALLY_INFEASIBLE_WT'
LEGACY = task.REPO/'results/18_fixed_wt_inverse_rebalance'
CLASSIFICATION = 'TASK 20 NUMERICALLY INCONCLUSIVE'


def blob_sha(path):
    data = Path(path).read_bytes()
    return hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()


def checked_checkpoint():
    path = task.OUT/'frozen_wt_manifest.json'
    receipt = task.read(task.OUT/'pushed_wt_checkpoint.json')
    if receipt['checkpoint_sha'] != CHECKPOINT or receipt['branch'] != task.BRANCH:
        raise AssertionError('Wrong Task 20 WT checkpoint')
    if not receipt['remote_ref_verified'] or not receipt['remote_ref_verified_utc']:
        raise AssertionError('WT checkpoint has not been remotely verified')
    if receipt['remote_manifest_git_blob_sha'] != MANIFEST_BLOB or blob_sha(path) != MANIFEST_BLOB:
        raise AssertionError('WT manifest differs from pushed checkpoint')
    frozen = task.read(path)
    if frozen['genotype_outcomes_used'] or frozen['genotype_result_files_read']:
        raise AssertionError('Phenotype firewall violated')
    if len(frozen['all_decisions']) != 20:
        raise AssertionError('Incomplete WT decisions')
    for name,digest in frozen['hashes'].items():
        if sha256_file(task.REPO/name) != digest:
            raise AssertionError(f'Frozen artifact changed: {name}')
    return frozen


def csv_rows(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def inherited_file(path, verified):
    name = task.relative(path)
    expected = task.git('rev-parse', task.INTAKE+':'+name)
    if blob_sha(path) != expected:
        raise AssertionError(f'Inherited Git blob changed: {name}')
    verified[name] = {'sha256':sha256_file(path), 'git_blob_sha1':expected}


def verify_legacy():
    verified = {}
    contract_path = task.REPO/'results/19_ae4_loss_nak_pump_coupling/contract.json'
    inherited_file(contract_path, verified)
    frozen_inputs = task.read(contract_path)['input_files']
    if len(frozen_inputs) != 187:
        raise AssertionError('Unexpected inherited frozen-input count')
    for name, item in frozen_inputs.items():
        path = task.REPO/name
        if sha256_file(path) != item['sha256'] or blob_sha(path) != item['git_blob_sha1']:
            raise AssertionError(f'Frozen inherited input changed: {name}')
        verified[name] = item
    summary_path = LEGACY/'postfreeze_summary.json'
    inherited_file(summary_path, verified)
    for name, digest in task.read(summary_path)['result_hashes'].items():
        inherited_file(task.REPO/name, verified)
        if sha256_file(task.REPO/name) != digest:
            raise AssertionError(f'Legacy result digest failed: {name}')
    tables = {}
    for name in ('wt_dynamic_results.csv','ae4_5pct_results.csv','ae2_results.csv','final_comparison.csv'):
        rows = [r for r in csv_rows(LEGACY/name) if r['condition'] in ('inherited_baseline','share_0.10')]
        if len(rows) != 60:
            raise AssertionError('Expected exactly 60 legacy comparator rows')
        tables[name] = rows
        for row in rows:
            if row.get('source_path'):
                path = task.REPO/row['source_path']
                if sha256_file(path) != row['source_sha256']:
                    raise AssertionError('Legacy trajectory source hash failed')
                inherited_file(path, verified)
    manifest, _ = load_freeze()
    for root in manifest['roots']:
        path = LEGACY/'postfreeze_payloads'/f'{inverse.slug(root)}__share_0.10.json'
        inherited_file(path, verified)
        payload = task.read(path)
        if sha256_file(task.REPO/payload['wt_payload_path']) != payload['wt_payload_sha256']:
            raise AssertionError('Legacy WT parameter payload hash failed')
        if payload['fitted_parameter_count_after_freeze'] != 0:
            raise AssertionError('Legacy comparator refitted after freeze')
    return tables, verified


def key(row):
    return row['root_id'], row['condition'], float(row['calcium_uM'])


def legacy_comparison(tables):
    lookup = {name:{key(r):r for r in rows} for name,rows in tables.items()}
    rows = []
    for r in tables['final_comparison.csv']:
        q = {name:float(lookup[name][key(r)]['total_0_600_pL']) for name in
             ('wt_dynamic_results.csv','ae4_5pct_results.csv','ae2_results.csv')}
        qwt = q['wt_dynamic_results.csv']
        if not qwt > 0 or not all(math.isfinite(v) for v in q.values()):
            raise AssertionError('Legacy comparator has invalid flow integral')
        if not all(r[name] == 'True' for name in ('wt_numerical_gate_pass',
                   'ae4_paired_numerical_gate_pass','ae2_paired_numerical_gate_pass')):
            raise AssertionError('Invalid legacy numerical gate')
        for ratio, table in (('R_AE4_5pct','ae4_5pct_results.csv'),('R_AE2','ae2_results.csv')):
            if not math.isclose(float(r[ratio]), q[table]/qwt, rel_tol=2e-14):
                raise AssertionError('Legacy ratio disagrees with saved integrals')
        rows.append({
            'architecture':'legacy_SR2', 'condition':'legacy_inherited' if r['condition']=='inherited_baseline' else 'legacy_P10',
            'legacy_condition':r['condition'], 'root_id':r['root_id'], 'routing_family':r['routing_family'],
            'calcium_uM':float(r['calcium_uM']), 'status':'HASH_VALID_LEGACY_REUSE',
            'fixed_wt_feasible':True, 'new_simulation':False,
            'wt_total_0_600_pL':qwt, 'ae4_5pct_total_0_600_pL':q['ae4_5pct_results.csv'],
            'ae2_total_0_600_pL':q['ae2_results.csv'], 'R_AE4':float(r['R_AE4_5pct']),
            'R_AE2':float(r['R_AE2']), 'ae4_direction':'increase' if float(r['R_AE4_5pct'])>1 else 'decrease',
            'paired_gate_pass':True, 'comparison_to_pooled_available':False,
            'source_path':task.relative(LEGACY/'final_comparison.csv'),
            'source_sha256':sha256_file(LEGACY/'final_comparison.csv'),
        })
    return rows


def unavailable_tables(frozen):
    rows = []
    for decision in frozen['all_decisions']:
        for calcium in task.CALCIUM:
            rows.append({
                'candidate_id':decision['candidate_id'], 'root_id':decision['root_id'],
                'routing_family':decision['routing_family'], 'condition':decision['condition'],
                'calcium_uM':calcium, 'architecture':MODEL_ID,
                'status':UNAVAILABLE, 'reason':decision['status'],
                'fixed_wt_feasible':False, 'simulation_attempted':False,
                'solver_success':None, 'numerical_gate_pass':None,
                'total_0_600_pL':None, 'secretion_ratio':None,
                'capacity_refits':0, 'wt_checkpoint_sha':CHECKPOINT,
            })
    write_rows(task.OUT/'wt_dynamic_results.csv', [{**r,'genotype':'WT','ae4_expression':1.,'ae2_expression':1.} for r in rows])
    write_rows(task.OUT/'ae4_5pct_results.csv', [{**r,'genotype':'AE4_5pct','ae4_expression':.05,'ae2_expression':1.} for r in rows])
    write_rows(task.OUT/'ae2_results.csv', [{**r,'genotype':'AE2_loss','ae4_expression':1.,'ae2_expression':0.} for r in rows])
    return [{
        'architecture':MODEL_ID, 'condition':r['condition'], 'legacy_condition':None,
        'root_id':r['root_id'], 'routing_family':r['routing_family'], 'calcium_uM':r['calcium_uM'],
        'status':UNAVAILABLE, 'fixed_wt_feasible':False, 'new_simulation':False,
        'wt_total_0_600_pL':None, 'ae4_5pct_total_0_600_pL':None, 'ae2_total_0_600_pL':None,
        'R_AE4':None, 'R_AE2':None, 'ae4_direction':'not_evaluable',
        'paired_gate_pass':None, 'comparison_to_pooled_available':False,
        'source_path':task.relative(task.OUT/'frozen_wt_manifest.json'),
        'source_sha256':sha256_file(task.OUT/'frozen_wt_manifest.json'),
    } for r in rows]


def state_flux_summary(frozen):
    manifest, _ = load_freeze()
    rows = []
    for decision in frozen['all_decisions']:
        model, state, evaluation = inverse.baseline(manifest,decision['root_id'])
        d = evaluation.diagnostics
        ci = d.observables.cell_concentrations_mM
        p = task.read(task.OUT/'wt_parameter_payloads'/f"{decision['candidate_id']}.json")
        f = p['pooled_reference_flux']
        # These fields are inherited target diagnostics, not pooled WT or
        # genotype resting states; the accepted fields remain explicitly null.
        inherited = {'inherited_cl_i_mM':ci['cl'],'inherited_ph_i':d.observables.cell_acid_base.ph,
            'inherited_na_i_mM':ci['na'],'inherited_k_i_mM':ci['k'],
            'inherited_hco3_i_mM':ci['hco3'],
            'inherited_NKCC1_cycles_fmol_s':d.homeostasis.nkcc1_inward_fmol_s,
            'inherited_NHE1_fmol_s':d.homeostasis.nhe1_inward_fmol_s,
            'inherited_AE2_cl_fmol_s':d.homeostasis.ae2_inward_fmol_s,
            'inherited_pump_apical_cycles_fmol_s':d.membranes.pump_apical_fmol_s,
            'inherited_pump_basolateral_cycles_fmol_s':d.membranes.pump_basolateral_fmol_s,
            'reference_pooled_na_source_fmol_s':f['source_na_i_fmol_s'],
            'reference_pooled_k_source_fmol_s':f['source_k_i_fmol_s'],
            'reference_pooled_cl_source_fmol_s':f['j_ae4_fmol_s'],
            'reference_pooled_is_accepted_wt':False}
        for genotype in ('WT','AE4_5pct','AE2_loss'):
            rows.append({'candidate_id':decision['candidate_id'],'root_id':decision['root_id'],
                'condition':decision['condition'],'genotype':genotype,'status':UNAVAILABLE,
                'accepted_resting_cl_i_mM':None,'accepted_resting_ph_i':None,
                'accepted_resting_na_i_mM':None,'accepted_resting_k_i_mM':None,
                'delta_cl_i_mM':None,'delta_ph_i':None,'delta_na_i_mM':None,'delta_k_i_mM':None,
                'NKCC1_cycles_fmol_s':None,'NHE1_fmol_s':None,'AE2_cl_fmol_s':None,
                'pump_apical_cycles_fmol_s':None,'pump_basolateral_cycles_fmol_s':None,
                'lower_cl_after_ae4_loss':None, **inherited})
    write_rows(task.OUT/'state_flux_summary.csv',rows)


def main():
    frozen = checked_checkpoint()
    if frozen['feasible_parameterisations']:
        raise RuntimeError('This empty-set report cannot replace simulation of feasible candidates')
    if not all(r['status']==task.INFEASIBLE for r in frozen['all_decisions']):
        raise AssertionError('Every omitted simulation needs a structural certificate')
    tables, verified = verify_legacy()
    legacy = legacy_comparison(tables)
    pooled = unavailable_tables(frozen)
    write_rows(task.OUT/'final_comparison.csv',legacy+pooled)
    state_flux_summary(frozen)
    write_json(task.OUT/'legacy_reuse_manifest.json',{
        'created_utc':inverse.now(),'wt_checkpoint_sha':CHECKPOINT,
        'conditions_reused':['inherited_baseline','share_0.10'],
        'legacy_comparison_rows':len(legacy),'legacy_trajectories_recomputed':0,
        'frozen_inherited_inputs_unchanged':187,'verified_files':verified})
    ranges = {}
    for condition in ('legacy_inherited','legacy_P10'):
        selected = [r for r in legacy if r['condition']==condition]
        ranges[condition] = {name:[min(r[name] for r in selected),max(r[name] for r in selected)]
                             for name in ('R_AE4','R_AE2')}
    write_json(task.OUT/'postfreeze_summary.json',{
        'created_utc':inverse.now(),'wt_checkpoint_sha':CHECKPOINT,
        'primary_classification':CLASSIFICATION,
        'scientific_status':task.INFEASIBLE,
        'classification_limitation': 'The mandated six-label taxonomy has no structural-infeasibility label. Label 6 denotes unavailable phenotype comparisons only. The WT sign obstruction is proved, not numerical nonconvergence.',
        'wt_decisions':20,'feasible_P0':0,'feasible_P10':0,'certified_infeasible':20,
        'new_wt_trajectories':0,'new_ae4_trajectories':0,'new_ae2_trajectories':0,
        'genotype_continuations':0,'exact_ae4_zero_attempts':0,'capacity_refits_after_freeze':0,
        'pooled_dynamics_not_run_rows':180,'legacy_comparison_rows':60,'four_condition_comparison_rows':120,
        'pooled_no_slip_by_construction':True,'pooled_simulated_state_count':0,
        'pooled_violation_in_assays_or_reference_evaluations':False,
        'legacy_ranges':ranges,'ae4_direction_P0':None,'ae4_direction_P10':None,
        'ae4_loss_resting_cl_direction':None,'pooled_ae2_neutrality':None,
        'robust_across_all_roots':'structural infeasibility only; no pooled phenotype claim',
        'experimental_deficit_used_as_target':False,
    })
    print('Reused 60 legacy comparisons; recorded 180 unavailable pooled dynamics rows. No new simulations.')


if __name__ == '__main__':
    main()
