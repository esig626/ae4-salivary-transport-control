"""Post-push baseline reuse for Task 17's certified empty modified set."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from . import task17_fixed_wt as task
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import trajectory_path
from .task16_postfreeze import baseline_pair
from .validation import sha256_file


def checked_checkpoint():
    receipt=task.read(task.OUT/'pushed_fixed_wt_checkpoint.json')
    path=task.OUT/'frozen_fixed_wt_manifest.json'
    assert receipt['branch']==task.BRANCH
    assert receipt['remote_ref_verified_utc']
    blob=task.git('hash-object',task.relative(path))
    assert blob==receipt['remote_manifest_git_blob_sha']
    assert blob==task.git('rev-parse',receipt['checkpoint_sha']+':'+task.relative(path))
    task.git('merge-base','--is-ancestor',receipt['checkpoint_sha'],'HEAD')
    freeze=task.read(path)
    for name,digest in freeze['hashes'].items():
        assert sha256_file(task.REPO/name)==digest,name
    assert len(freeze['all_decisions'])==30
    assert len(freeze['feasible_parameterizations'])==10
    assert all(r['condition']=='inherited_baseline'
        for r in freeze['all_decisions'] if r['fixed_wt_feasible']), (
        'This reuse runner cannot replace required simulations for a modified feasible endpoint')
    return freeze,receipt


def evaluate():
    freeze,receipt=checked_checkpoint()
    assert not (task.OUT/'genotype_results_frozen.json').exists()
    inputs=task.check_inputs()
    models,_=load_freeze()
    ae2_path=task.REPO/'results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv'
    ae2_blob=task.git('hash-object',task.relative(ae2_path))
    assert ae2_blob==task.git('rev-parse',task.BASE+':'+task.relative(ae2_path))
    with ae2_path.open() as handle:
        ae2_source=list(csv.DictReader(handle))
    ae2_rows={(r['root_id'],float(r['calcium_uM'])):r for r in ae2_source}
    assert len(ae2_rows)==len(ae2_source)==30
    rows=[]
    for decision in freeze['all_decisions']:
        for calcium in (.10,.25,.50):
            row={k:decision[k] for k in ('candidate_id','root_id','routing_family',
                'condition','target_share','fixed_wt_feasible','ae4_capacity_multiplier','rho')}
            row.update(calcium_uM=calcium,wt_status=decision['status'])
            if decision['fixed_wt_feasible']:
                wt_path=trajectory_path(decision['root_id'],calcium,'WT','production_radau')
                expected=inputs['input_files'][task.relative(wt_path)]['sha256']
                assert sha256_file(wt_path)==expected
                wt=task.read(wt_path)
                inherited_decision={**decision,'calcium_uM':calcium,'wt_admissible':True,
                    'ae4_capacity_scale':1.,'nkcc1_capacity_scale':1.,
                    'trajectory_path':task.relative(wt_path),
                    'wt_total_0_600_pL':wt['trajectory']['cumulative_flow_pL'][-1]}
                # This existing reader rechecks production solver, exact core,
                # parameter hashes, 600s grid, paired denominator and integral.
                row.update(baseline_pair(inherited_decision,models,ae2_rows,inputs))
                row['wt_source_path']=task.relative(wt_path)
                row['wt_source_sha256']=expected
                row['ae2_source_sha256']=sha256_file(ae2_path)
            else:
                row.update(status='NOT_EVALUATED_FIXED_WT_INFEASIBLE',
                    genotype_numerical_status='NOT_ELIGIBLE',
                    paired_numerical_gate_pass=None,new_genotype_simulation=False,
                    exclusion_reason=decision['feasibility_failure_kind'])
            rows.append(row)
    assert len(rows)==90
    common=['candidate_id','root_id','routing_family','condition','target_share','calcium_uM',
        'fixed_wt_feasible','wt_status','status','genotype_numerical_status',
        'paired_numerical_gate_pass','new_genotype_simulation','wt_total_0_600_pL',
        'wt_source_path','wt_source_sha256','exclusion_reason']
    for file,keys in (
        ('ae4_5pct_results.csv',common+['ae4_expression','R_AE4_5pct','D_AE4_5pct',
            'ae4_5pct_total_0_600_pL','ae4_5pct_source_path','ae4_5pct_source_sha256']),
        ('ae2_results.csv',common+['ae2_expression_for_loss','R_AE2','D_AE2',
            'ae2_source_path','ae2_source_sha256','ae2_reuse_scope'])):
        write_rows(task.OUT/file,[{k:r.get(k) for k in keys} for r in rows])
    write_rows(task.OUT/'paired_baseline_results.csv',rows)
    valid=[r for r in rows if r['fixed_wt_feasible']]
    summary={key:{'minimum':min(r[key] for r in valid),
        'median':float(np.median([r[key] for r in valid])),
        'maximum':max(r[key] for r in valid)}
        for key in ('R_AE4_5pct','R_AE2','AE4_minus_AE2_reduction')}
    result={'created_utc':task.now(),'pushed_wt_checkpoint':receipt['checkpoint_sha'],
        'remote_checkpoint_verified_utc':receipt['remote_ref_verified_utc'],
        'baseline_pairs_reused':len(valid),'modified_combinations_not_evaluated':60,
        'new_genotype_simulations':0,'new_WT_state_solves':0,
        'new_genotype_state_solves':0,'exact_AE4_zero_attempts':0,
        'every_frozen_decision_accounted':True,'baseline_summary':summary,
        'source_status':'Hash-valid inherited Task14/14B; modified endpoints have no loss prediction',
        'additional_input':{'path':task.relative(ae2_path),'baseline':task.BASE,
            'git_blob_sha1':ae2_blob,'sha256':sha256_file(ae2_path)},
        'hashes':{task.relative(task.OUT/name):sha256_file(task.OUT/name)
            for name in ('ae4_5pct_results.csv','ae2_results.csv','paired_baseline_results.csv')}}
    write_json(task.OUT/'genotype_results_frozen.json',result)
    print(json.dumps(result['baseline_summary']))


if __name__=='__main__':
    evaluate()
