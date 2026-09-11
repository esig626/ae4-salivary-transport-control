"""WT candidate inventory, production dynamics, and pre-genotype ensemble freeze."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import argparse
import numpy as np
from .task21_calibration import *
from .task21_contract import CALCIUM
from .task14_blind import build_model
from .run_calcium_fast_screen import summarise
from .validation import PRODUCTION_RADAU, physical_time_grid, simulate_wt


def equivalent(a,b):
    def flatten(value):
        if isinstance(value,dict):
            return [v for k in sorted(value) for v in flatten(value[k])]
        if isinstance(value,(list,tuple)):return [v for item in value for v in flatten(item)]
        return [float(value)]
    for field in ('complete_state','whole_cell_parameters','ae4_parameters'):
        u,v=np.array(flatten(a[field])),np.array(flatten(b[field]))
        if u.shape!=v.shape or not np.all(abs(u-v)<=1e-6*np.maximum(np.maximum(abs(u),abs(v)),1e-12)):
            return False
    return True


def all_outputs():
    for folder,collection in (('seed_attempts','attempts'),('refinements','attempts'),
                               ('exact_calibration','stages')):
        for path in sorted((OUT/folder).glob('*.json')):
            if '.progress.' in path.name:continue
            result=read(path)
            for i,item in enumerate(result[collection]):
                yield result['root_id'],str(path.relative_to(REPO))+f'#{i}',item.get('candidate')
    for folder in ('flux_search','mean_flux_search','joint_tie_break'):
        for path in sorted((OUT/folder).glob('*.json')):
            result=read(path)
            yield result['root_id'],str(path.relative_to(REPO)),result.get('candidate')
    for path in sorted((OUT/'profile_calibration').glob('*.json')):
        result=read(path)
        for i,item in enumerate((result.get('payload') or {}).get('versions',[])):
            yield result['root_id'],str(path.relative_to(REPO))+f'#{i}',item.get('candidate')


def inventory():
    require_stage0();candidates=[];decisions=[]
    for root,source,payload in all_outputs():
        if payload is None:
            decisions.append(dict(root_id=root,source=source,resting_gate_pass=False,
                decision='INVALID_CAPACITY_OR_STATE',gate_failures=['invalid_candidate']))
            continue
        original_gate=payload['resting_gate_pass']
        original_closure=payload['closure']
        model=candidate_model({**payload,'root_id':root},.10)
        ev=model.evaluate(0.,payload['complete_state'],genotype=WT)
        summary=inverse.residual_summary(ev)
        failures=physiology_failures(payload['observables'])
        for name,tol in (('max_abs_scaled_independent_rhs',SPEC.root_scaled_tolerance),
            ('max_abs_omitted_rhs',SPEC.omitted_row_raw_tolerance),('max_abs_current_A',SPEC.current_tolerance_A),
            ('max_abs_state_charge_fmol',SPEC.charge_tolerance_fmol),('max_conservation_tolerance_ratio',1.)):
            if not np.isfinite(summary[name]) or summary[name]>tol:failures.append(name)
        ae4=ev.diagnostics.ae4
        if ae4.diagnostics['affinity']<=0:failures.append('nonpositive_affinity')
        if ae4.cl_cell_fmol_s<=0:failures.append('nonpositive_AE4_loading')
        if ae4.na_cell_fmol_s*ae4.k_cell_fmol_s<0:failures.append('opposing_cation_sources')
        if not np.all(np.isfinite(payload['complete_state'])) or min(payload['complete_state'][:12])<=0:
            failures.append('invalid_core')
        payload={**payload,'initial_numeric_type_gate_pass':original_gate,'initial_numeric_type_closure':original_closure,
            'resting_gate_pass':not failures,'resting_gate_failures':failures,'closure':summary,
            'production_diagnostics':asdict(ev.diagnostics),'complete_rhs':ev.rhs.tolist(),
            'validation_representation':'JSON-reloaded Python float capacities and state'}
        decisions.append(dict(root_id=root,source=source,resting_gate_pass=payload['resting_gate_pass'],
            decision='RESTING_GATE_PASS' if payload['resting_gate_pass'] else 'RESTING_GATE_FAILURE',
            initial_numeric_type_gate_pass=original_gate,
            gate_failures=payload['resting_gate_failures'],objective=payload['objective'],
            closure=payload['closure']))
        # Preserve every previously tested rest witness, including replay failures.
        if not payload['resting_gate_pass'] and not original_gate:continue
        match=next((c for c in candidates if equivalent(c,payload)),None)
        if match is not None:
            match['seed_provenance'].append(dict(root_id=root,source=source));continue
        identity=sha256_object({k:payload[k] for k in ('complete_state','whole_cell_parameters','ae4_parameters')})
        candidates.append({**payload,'solution_id':'WT_'+identity[:12],'root_id':root,
                           'seed_provenance':[dict(root_id=root,source=source)]})
    candidates.sort(key=lambda c:tuple(c['objective']))
    for c in candidates:write_json(OUT/'wt_candidates'/f"{c['solution_id']}.json",c)
    write_json(OUT/'wt_candidate_inventory.json',dict(created_utc=inverse.now(),
        candidate_ids=[c['solution_id'] for c in candidates],recorded_attempts=len(decisions),
        resting_gate_pass_count=sum(d['resting_gate_pass'] for d in decisions),decisions=decisions,
        deduplication_relative_tolerance=1e-6,genotype_evaluations=0))
    print('Distinct WT witnesses:',len(candidates),'replay-valid:',sum(c['resting_gate_pass'] for c in candidates),flush=True)
    return candidates


def candidate_model(candidate,calcium):
    manifest,_=load_freeze()
    template=select_ae4_mechanism(build_model(manifest,candidate['root_id'],calcium),MODEL_ID)
    model=with_capacities(template,candidate['calibrated_capacities'])
    if sha256_object(asdict(model.parameters))!=sha256_object(candidate['whole_cell_parameters']):
        raise AssertionError('Candidate parameters changed during WT model reconstruction')
    return model


def dynamic(job):
    candidate,ca=job;model=candidate_model(candidate,ca)
    path=OUT/'wt_trajectories'/f"{candidate['solution_id']}__Ca{ca:.2f}.json"
    if path.exists():
        saved=read(path)
        if saved['whole_cell_parameters_sha256']!=sha256_object(model.parameters):
            raise AssertionError('Cached WT parameters changed')
        if saved['trajectory'] is not None:
            if not np.array_equal(np.asarray(saved['trajectory']['states'])[:,0],candidate['complete_state']):
                raise AssertionError('Cached WT initial state changed')
        saved['replay_inventory_candidate_sha256']=sha256_object(candidate)
        saved['initial_state_sha256']=sha256_object(candidate['complete_state'])
        saved['replay_inventory_metadata_only_update']=True
        write_json(path,saved)
        return saved['summary']
    row=dict(solution_id=candidate['solution_id'],root_id=candidate['root_id'],calcium_uM=ca,
             genotype='WT',ae4_expression=1.,ae2_expression=1.,started_utc=inverse.now())
    trajectory=None;grid=physical_time_grid(600.,1.)
    try:
        trajectory=simulate_wt(model,candidate['complete_state'],family='TASK21_WT',
            initial_condition=candidate['solution_id'],solver=PRODUCTION_RADAU,time_s=grid)
        summary=summarise(trajectory);states=trajectory.states
        charge=max(float(max(abs(states[0]+states[1]-states[2]-states[4]-model.parameters.geometry.fixed_cell_anion_equivalents_fmol))),
                   float(max(abs(states[6]+states[7]-states[8]-states[10]))))
        fractions=bool(np.all(np.isfinite(states[12:])) and np.all(states[12:]>=0) and np.all(states[12:]<=1))
        passed=bool(summary['numerical_gate_pass'] and np.array_equal(trajectory.time_s,grid)
                    and charge<=SPEC.charge_tolerance_fmol and fractions)
        row.update(dynamic_gate_pass=passed,solver_success=trajectory.success,
            solver_message=trajectory.message,max_abs_charge_fmol=charge,regulatory_fractions_valid=fractions,
            conservation_gate_pass=summary['conservation_gate_pass'],
            max_conservation_tolerance_ratio=summary['max_dimensionless_conservation_ratio'],
            positive_core=trajectory.positive_core,output_samples=len(trajectory.time_s),
            Q_0_600_pL=float(trajectory.cumulative_flow_pL[-1]) if passed else None,
            raw_integrated_flow_pL=float(trajectory.cumulative_flow_pL[-1]))
    except Exception as exc:
        row.update(dynamic_gate_pass=False,solver_success=False,
            solver_message=f'{type(exc).__name__}: {exc}',Q_0_600_pL=None)
    row.update(completed_utc=inverse.now(),status='VALID' if row['dynamic_gate_pass'] else 'DYNAMIC_GATE_FAILURE')
    write_json(path,dict(summary=row,trajectory=asdict(trajectory) if trajectory is not None else None,
        solver=asdict(PRODUCTION_RADAU),candidate_sha256=sha256_object(candidate),
        whole_cell_parameters_sha256=sha256_object(model.parameters)))
    print(row['solution_id'],ca,row['status'],row.get('max_conservation_tolerance_ratio'),flush=True)
    return row


def dynamics(workers):
    candidates=inventory();jobs=[(c,ca) for c in candidates for ca in CALCIUM]
    with ProcessPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(dynamic,jobs))
    write_rows(OUT/'wt_dynamic_results.csv',rows or [dict(status='NOT_RUN_NO_REST_VALID_WT',dynamic_gate_pass=False)])
    write_json(OUT/'wt_dynamic_results.json',rows)


def freeze():
    require_stage0();manifest,_=load_freeze();inventory=read(OUT/'wt_candidate_inventory.json')
    for folder in ('seed_attempts','refinements','exact_calibration','mean_flux_search',
                   'profile_calibration','joint_tie_break'):
        roots={read(p)['root_id'] for p in (OUT/folder).glob('*.json') if '.progress.' not in p.name}
        if roots!=set(manifest['roots']):raise AssertionError(f'Incomplete seed ensemble: {folder}')
    rows=read(OUT/'wt_dynamic_results.json');retained=[];table=[]
    for cid in inventory['candidate_ids']:
        c=read(OUT/'wt_candidates'/f'{cid}.json');panel=[r for r in rows if r['solution_id']==cid]
        valid=c['resting_gate_pass'] and len(panel)==3 and all(r['dynamic_gate_pass'] for r in panel)
        if valid:retained.append(cid)
        table.append(dict(solution_id=cid,root_id=c['root_id'],**c['observables'],
            measured_objective=c['objective'][0],max_abs_log_fold=c['objective'][1],
            sum_squared_log_fold=c['objective'][2],resting_gate_pass=c['resting_gate_pass'],wt_dynamic_gate_pass=valid,
            initial_numeric_type_gate_pass=c['initial_numeric_type_gate_pass'],
            admissible=valid,positive_affinity=c['production_diagnostics']['ae4']['diagnostics']['affinity'],
            positive_AE4_loading_fmol_s=c['production_diagnostics']['ae4']['cl_cell_fmol_s']))
    for root in sorted(manifest['roots']):
        if not any(r['root_id']==root for r in table):
            table.append(dict(root_id=root,solution_id='',resting_gate_pass=False,
                wt_dynamic_gate_pass=False,admissible=False,status='NO_REST_VALID_CANDIDATE'))
    write_rows(OUT/'wt_calibration_table.csv',table)
    flux_rows=[];fold_rows=[]
    for cid in inventory['candidate_ids']:
        c=read(OUT/'wt_candidates'/f'{cid}.json');d=c['production_diagnostics'];o=d['observables']
        row=dict(solution_id=cid,root_id=c['root_id'],resting_replay_gate_pass=c['resting_gate_pass'],
            volume_cell_pL=c['complete_state'][5],volume_lumen_pL=c['complete_state'][11],
            pH_cell=o['cell_acid_base']['ph'],pH_lumen=o['lumen_acid_base']['ph'])
        for side in ('cell','lumen'):
            row.update({f'{ion}_{side}_mM':v for ion,v in o[f'{side}_concentrations_mM'].items()})
        row.update({f'AE4_{k}':v for k,v in d['ae4']['diagnostics'].items() if isinstance(v,(int,float))})
        row.update({f'homeostasis_{k}':v for k,v in d['homeostasis'].items() if isinstance(v,(int,float))})
        row.update({f'current_A_{k}':v for k,v in d['membranes']['currents_A'].items()})
        row.update({f'CO2_fmol_s_{k}':v for k,v in d['co2_fluxes_fmol_s'].items()})
        row.update({f'outflow_fmol_s_{k}':v for k,v in d['outflow_sources_fmol_s'].items()})
        for name in ('apical','basolateral','transepithelial'):
            row[f'V_{name}_mV']=1000*d['membranes'][f'v_{name}_V']
        row.update(AE4_na_cell_fmol_s=d['ae4']['na_cell_fmol_s'],
                   AE4_k_cell_fmol_s=d['ae4']['k_cell_fmol_s'],
                   no_opposing_cation_cycle=d['ae4']['na_cell_fmol_s']*d['ae4']['k_cell_fmol_s']>=0)
        flux_rows.append(row)
        for name in CAPACITIES:
            fold_rows.append(dict(solution_id=cid,root_id=c['root_id'],capacity=name,
                reference=c['reference_capacities'][name],calibrated=c['calibrated_capacities'][name],
                fold=c['capacity_folds'][name],log_fold=math.log(c['capacity_folds'][name])))
    write_rows(OUT/'resting_flux_diagnostics.csv',flux_rows)
    write_rows(OUT/'capacity_folds.csv',fold_rows)
    write_rows(OUT/'wt_seed_decisions.csv',[dict(root_id=root,
        saved_solver_outputs=sum(d['root_id']==root for d in inventory['decisions']),
        replay_rest_valid_outputs=sum(d['root_id']==root and d['resting_gate_pass'] for d in inventory['decisions']),
        admissible_WT_solutions=sum(r['root_id']==root and r['admissible'] for r in table),
        classification='NO_VALIDATED_WT_NUMERICAL_FAILURE_NOT_INFEASIBILITY_PROOF' if not retained else 'SEE_CALIBRATION_TABLE')
        for root in sorted(manifest['roots'])])
    write_json(OUT/'wt_validation_summary.json',dict(created_utc=inverse.now(),
        seed_count=len(manifest['roots']),saved_solver_outputs=inventory['recorded_attempts'],
        initial_rest_witnesses=len(inventory['candidate_ids']),
        replay_rest_valid_outputs=inventory['resting_gate_pass_count'],WT_dynamic_conditions=len(rows),
        WT_integrations_complete=sum(r['solver_success'] for r in rows),
        WT_dynamic_valid=sum(r['dynamic_gate_pass'] for r in rows),admissible_WT_count=len(retained),
        task21_genotype_evaluations=0,archive_tree_sha=git('rev-parse','HEAD:archive'),
        focused_test_transcript='results/21_pooled_ae4_wt_physiology/validation_wt_tests.txt',
        broad_test_limitation_and_ordering_deviation='results/21_pooled_ae4_wt_physiology/validation_scope_deviation.json'))
    files=[p for p in OUT.rglob('*') if p.is_file() and not p.name.endswith('.progress.json')]
    code=list((REPO/'src/modern_full_model').glob('task21_*.py'))
    code += [REPO/'tests/test_task21_pooled_wt_physiology.py',
             REPO/'analysis/21_pooled_ae4_wt_physiology/methods.md',
             REPO/'analysis/21_pooled_ae4_wt_physiology/wt_physiology_contract.md']
    frozen=dict(created_utc=inverse.now(),stage0_checkpoint_sha=STAGE0_SHA,
        seed_root_ids=sorted(manifest['roots']),tested_candidate_ids=inventory['candidate_ids'],
        rest_valid_candidate_ids=[cid for cid in inventory['candidate_ids']
            if read(OUT/'wt_candidates'/f'{cid}.json')['resting_gate_pass']],
        retained_solution_ids=retained,retained_count=len(retained),wt_dynamic_conditions=len(rows),
        wt_dynamic_valid=sum(r['dynamic_gate_pass'] for r in rows),genotype_evaluations=0,
        genotype_permission='REQUIRES_REMOTE_WT_CHECKPOINT' if retained else 'STOP_NO_ADMISSIBLE_WT',
        classification='WT_ENSEMBLE_REQUIRES_CHECKPOINT' if retained else 'NO_VALIDATED_WT_NUMERICAL_FAILURE_NOT_INFEASIBILITY_PROOF',
        file_sha256={str(p.relative_to(REPO)):sha256_file(p) for p in files+code})
    write_json(OUT/'wt_manifest.json',frozen)
    print('Frozen WT retained:',len(retained),'of',len(inventory['candidate_ids']),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['inventory','dynamic','freeze'])
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if args.stage=='inventory':inventory()
    elif args.stage=='dynamic':dynamics(args.workers)
    else:freeze()
