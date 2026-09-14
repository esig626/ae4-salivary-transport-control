"""One gated Task 41 integration per invocation, ordered WT -> 5% -> null.

The common ledger permits at most three intended integrations and one BDF
retry after an explicit Radau numerical failure. No retry follows a scientific
failure. Every accepted endpoint and integer-second dense sample is checked.
"""
import time
START = time.perf_counter()
from candidate_common import *
from dataclasses import replace
from datetime import datetime, timezone
import argparse
import csv
import signal
import traceback
import scipy
from scipy.integrate import Radau, BDF, cumulative_trapezoid
from threadpoolctl import threadpool_info
from modern_full_model.validation import PRODUCTION_RADAU, PRODUCTION_BDF

CASES = {"wt":1.,"ae4_5pct":0.05,"ae4_null":0.}
FAILURES = {"wt":"TASK41_RECRUITMENT_WT_FAILED",
    "ae4_5pct":"TASK41_RECRUITMENT_AE4_5PCT_FAILED","ae4_null":"TASK41_RECRUITMENT_AE4_NULL_FAILED"}
QUADRATURE = "Trapezoidal on integer seconds plus t=0 and 1e-6; identical to Tasks 37-38."
FLUXES = (
    ("positive_nkcc1_cl_loading","nkcc1_cl_inward_fmol_s",True,"fmol"),
    ("positive_ae4_cl_loading","ae4_cl_inward_fmol_s",True,"fmol"),
    ("positive_ae2_cl_loading","ae2_cl_inward_signed_fmol_s",True,"fmol"),
    ("signed_ae2_cl_flux","ae2_cl_inward_signed_fmol_s",False,"fmol"),
    ("nkcc1_cycles","nkcc1_cycle_inward_fmol_s",False,"fmol"),
    ("ae4_signed_cl_flux","ae4_cl_inward_fmol_s",False,"fmol"),
    ("ae4_na_cell_source","ae4_na_cell_source_fmol_s",False,"fmol"),
    ("ae4_k_cell_source","ae4_k_cell_source_fmol_s",False,"fmol"),
    ("ae4_tic_cell_source","ae4_tic_cell_source_fmol_s",False,"fmol"),
    ("ae4_ta_cell_source","ae4_ta_cell_source_fmol_s",False,"fmol"),
    ("nbc_cycles","nbc_cycle_inward_fmol_s",False,"fmol"),
    ("nbc_bicarbonate_equivalent_influx","nbc_hco3_inward_fmol_s",False,"fmol"),
    ("nbc_current_charge","nbc_current_A",False,"C"),
    ("nhe1_flux","nhe1_inward_fmol_s",False,"fmol"),
    ("pump_cycles","pump_total_cycles_fmol_s",False,"fmol"),
    ("cacc_cl_export","cacc_cl_export_fmol_s",False,"fmol"),
    ("paracellular_cl_return","paracellular_cl_return_fmol_s",False,"fmol"),
    ("water_outflow","q_out_pL_s",False,"pL"),
)


class ScientificStop(Exception): pass
class ComputeStop(Exception): pass


def write_csv(path, rows):
    if not rows: return
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def run_case(case):
    manifest=verify_frozen_sources()
    budget=json.loads((OUT/'budget.json').read_text())
    assert budget['source_tests_pass'] and budget['wt_rest_pass']
    assert not budget['no_further_scientific_execution']
    preceding=list(CASES)[:list(CASES).index(case)]
    assert budget['started_cases']==preceding and budget['completed_cases']==preceding
    for prior in preceding:
        assert json.loads((OUT/(prior+'_verification.json')).read_text())['status']=='PASS'
    assert budget['numerical_execution_seconds']<900
    saved,rest_model,model,y0,rest_record=frozen_models_and_state()
    Q_REST=rest_record['observables']['q_out_pL_s']
    frozen_wt=json.loads((ROOT/'results/37_wt_nbc_validation/wt_verification.json').read_text())
    assert np.array_equal(y0,rest_record['state_vector'])
    assert asdict(model.stimulus)==frozen_wt['stimulus']
    assert asdict(PRODUCTION_RADAU)==frozen_wt['integration_attempts'][0]['solver']
    assert dict(CONSERVATION_RESIDUAL_TOLERANCES)==frozen_wt['conservation_tolerances']
    assert {'numpy':np.__version__,'scipy':scipy.__version__}==frozen_wt['software_versions']
    hashes=parameter_hashes(model)
    assert hashes==manifest['stimulus_parameter_hashes']
    genotype=replace(WT,ae4_expression=CASES[case])
    diff={k:v for k,v in asdict(genotype).items() if v!=asdict(WT)[k]}
    assert diff==({} if case=='wt' else {'ae4_expression':CASES[case]})
    pools=threadpool_info()
    assert pools and all(p['num_threads']==1 for p in pools)
    # Record intent before the first case evaluation or integration attempt.
    budget['started_cases'].append(case)
    budget['status']='RUNNING_'+case
    write_json(OUT/'budget.json',budget)
    rows,grid,attempts=[],[],[]
    minima,maxima,residual_max,ratio_max={},{},{},{}
    first_failure,previous_pass,raw_last,error=None,None,None,None
    status,solver,complete,accepted_steps=None,None,False,0
    started_utc=datetime.now(timezone.utc).isoformat()

    def check_time(*_):
        if budget['numerical_execution_seconds']+time.perf_counter()-START>=900:
            raise ComputeStop('Task 41 cumulative 900-second numerical execution limit')

    def alarm(*_):
        raise ComputeStop('Task 41 cumulative 900-second numerical execution limit')

    def inspect(t,y,kind):
        nonlocal first_failure,previous_pass,raw_last
        check_time()
        e=model.evaluate(float(t),y,genotype=genotype)
        row,failures,ratios=diagnose(model,t,y,e)
        # Verify the parent NKCC1 law and selected recruitment on every monitored state.
        from modern_full_model.nkcc1_palk2010 import palk_cycle_flux_fmol_s
        demand=palk_cycle_flux_fmol_s(row['na_i_mM'],row['k_i_mM'],row['cl_i_mM'],
            alpha_eff_fmol_s=model.nkcc1_kinetics.alpha_eff_fmol_s,
            activity_multiplier=row['nkcc1_activity_multiplier']*genotype.nkcc1_expression)
        expected=demand
        if abs(row['nkcc1_cycle_inward_fmol_s']-expected)>1e-14:
            failures.append({'gate':'Palk_law_unchanged','actual':row['nkcc1_cycle_inward_fmol_s'],'expected':expected})
        from modern_full_model.ae4_cacc_recruitment import recruitment_factor
        factor=recruitment_factor(0. if t==0 else 1.,genotype.ae4_expression,model.cacc_recruitment.null_stimulated_fraction)
        if row['cacc_recruitment_factor']!=factor:
            failures.append({'gate':'AE4_dependent_CaCC_recruitment_law'})
        if row['nkcc1_activity_multiplier']!=(1. if t==0 else 1.75):
            failures.append({'gate':'inherited_N1_multiplier_unchanged'})
        rows.append(row);raw_last=list(map(float,y))
        if t==0. or t==1e-6 or float(t).is_integer():
            if not grid or grid[-1]['time_s']!=t: grid.append(row)
        for k,v in row.items():
            if k!='time_s':
                minima[k]=min(minima.get(k,v),v);maxima[k]=max(maxima.get(k,v),v)
        for k in CONSERVATION_RESIDUAL_TOLERANCES:
            residual_max[k]=max(residual_max.get(k,0.),abs(float(e.diagnostics.conservation_residuals[k])))
            ratio_max[k]=max(ratio_max.get(k,0.),ratios[k])
        if failures:
            first_failure={'time_s':float(t),'previous_checked_pass_time_s':previous_pass,
                'sample_kind':kind,'failed_gates':failures,'observables':row,
                'state_vector':raw_last,'rhs':list(map(float,e.rhs)),
                'conservation_residuals':dict(e.diagnostics.conservation_residuals)}
            raise ScientificStop(FAILURES[case])
        previous_pass=float(t)

    signal.signal(signal.SIGALRM,alarm)
    signal.setitimer(signal.ITIMER_REAL,max(0.001,900-budget['numerical_execution_seconds']-(time.perf_counter()-START)))
    try:
        inspect(0.,y0,'frozen_task41_WT_rest_at_exact_onset')
        inspect(1e-6,y0,'inherited_post_onset_start')
        for specification in (PRODUCTION_RADAU,PRODUCTION_BDF):
            if attempts:
                assert attempts[-1]['solver_status']=='failed' and first_failure is None
                if budget['numerical_retries']>=1: break
                budget['numerical_retries']+=1
                rows.clear();grid.clear();minima.clear();maxima.clear();residual_max.clear();ratio_max.clear()
                previous_pass=None
                inspect(0.,y0,'frozen_task41_WT_rest_at_exact_onset')
                inspect(1e-6,y0,'inherited_post_onset_start')
            assert budget['integration_attempts']<4
            budget['integration_attempts']+=1
            write_json(OUT/'budget.json',budget)
            attempt={'attempt':len(attempts)+1,'global_attempt':budget['integration_attempts'],
                'solver':asdict(specification),'initial_state_sha256':sha256_object(list(map(float,y0)))}
            attempts.append(attempt);attempt_start=time.perf_counter()
            solver_class=Radau if specification.method=='Radau' else BDF
            solver=solver_class(lambda t,y:model.rhs(t,y,genotype=genotype),1e-6,y0.copy(),600.,
                rtol=specification.rtol,atol=specification.atol_vector(model.state_names),max_step=specification.max_step_s)
            next_time=1.
            while solver.status=='running':
                check_time();message=solver.step()
                if solver.status=='failed':
                    attempt['solver_message']=str(message);break
                accepted_steps+=1;dense=solver.dense_output()
                while next_time<=solver.t:
                    inspect(next_time,dense(next_time),'one_second_dense_output');next_time+=1
                if not rows or rows[-1]['time_s']!=float(solver.t):
                    inspect(float(solver.t),solver.y,'accepted_solver_endpoint')
            attempt.update(solver_status=solver.status,last_solver_time_s=float(solver.t),
                nfev=int(solver.nfev),njev=int(solver.njev),nlu=int(solver.nlu),wall_seconds=time.perf_counter()-attempt_start)
            if solver.status=='finished':
                complete=bool(solver.t==600.);break
        status='PASS' if complete else FAILURES[case]
    except ScientificStop:
        status=FAILURES[case]
    except ComputeStop as exc:
        status,error='TASK41_RECRUITMENT_COMPUTE_BUDGET_STOP',str(exc)
    except Exception:
        status,error='TASK41_RECRUITMENT_NUMERICAL_EXECUTION_ABORTED',traceback.format_exc()
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)

    if attempts and solver is not None:
        attempts[-1].update(solver_status=solver.status,stopped_by_scientific_gate=first_failure is not None,
            last_solver_time_s=float(solver.t),nfev=int(solver.nfev),njev=int(solver.njev),nlu=int(solver.nlu))
    assert parameter_hashes(model)==hashes
    verify_frozen_sources()
    secretion,partition,integrals=None,None,[]
    if complete and status=='PASS':
        times=np.asarray([r['time_s'] for r in grid]);q=np.asarray([r['q_out_pL_s'] for r in grid])
        assert np.array_equal(times,np.r_[0.,1e-6,np.arange(1.,601.)])
        cumulative=cumulative_trapezoid(q,times,initial=0.)
        for r,value in zip(grid,cumulative): r['cumulative_outflow_0_t_pL']=float(value)
        window=times>=60
        total=float(np.trapezoid(q,times));mean=float(np.trapezoid(q[window],times[window])/540)
        checks={'mean_60_600_gt_1p10_rest':mean>1.10*Q_REST,
            'cumulative_0_600_gt_1p10_rest':total>1.10*600*Q_REST,'endpoint_at_least_rest':float(q[-1])>=Q_REST}
        secretion={'cumulative_0_600_pL':total,'mean_60_600_pL_s':mean,'endpoint_pL_s':float(q[-1]),
            'q_rest_pL_s':Q_REST,'activation_checks':checks,'activation_required':case=='wt',
            'quadrature':QUADRATURE}
        if case=='wt' and not all(checks.values()): status=FAILURES[case]
        values={}
        for name,field,positive,unit in FLUXES:
            v=np.asarray([r[field] for r in grid])[window]
            if positive: v=np.maximum(v,0)
            values[name]=float(np.trapezoid(v,times[window]))
            integrals.append({'quantity':name,'value':values[name],'unit':unit,'window_start_s':60,'window_end_s':600})
        pool=sum(values['positive_'+k+'_cl_loading'] for k in ('nkcc1','ae4','ae2'))
        shares={k:values['positive_'+k+'_cl_loading']/pool for k in ('nkcc1','ae4','ae2')}
        partition={'window_s':[60,600],'positive_loading_pool_fmol':pool,'positive_shares':shares,
            'nkcc1_context_band':[0.65,0.75],'nkcc1_within_context':0.65<=shares['nkcc1']<=0.75,
            'validation_only':True,'definition':'Integrate max(J_Cl,0) for NKCC1, AE4 and AE2 separately, then normalize their sum.'}
    readouts=[r for r in grid if r['time_s']%60==0]
    if first_failure: readouts.append(first_failure['observables'])
    write_csv(OUT/(case+'_timeseries.csv'),readouts)
    write_csv(OUT/(case+'_dense_timeseries.csv'),grid)
    write_csv(OUT/(case+'_integrated_fluxes.csv'),integrals)
    verification={'status':status,'case':case,'ae4_expression':CASES[case],
        'prepared_head':PREPARED_HEAD,'branch':BRANCH,'started_utc':started_utc,
        'complete_600_s':complete,'first_failure':first_failure,'solver_exception':error,
        'integration_attempts':attempts,'accepted_solver_steps':accepted_steps,'checked_state_count':len(rows),
        'last_checked_time_s':rows[-1]['time_s'] if rows else None,
        'monitoring':frozen_wt['monitoring'],'onset_convention':frozen_wt['onset_convention'],
        'state_names':list(model.state_names),'initial_state_vector':list(map(float,y0)),
        'initial_state_sha256':sha256_object(list(map(float,y0))),
        'initial_state_bitwise_identical_to_frozen_task41_wt_rest':True,'last_checked_state_vector':raw_last,
        'core_state_sha256':rest_record['core_state_sha256'],'genotype':asdict(genotype),
        'only_genotype_difference':diff,'stimulus':asdict(model.stimulus),'parameter_hashes':hashes,
        'parent_parameters_unchanged_except_derived_CaCC_recruitment':True,
        'AE4_expression_effects':'AE4 transport and explicit stimulated CaCC recruitment dependency','ae4_routing_model':MODEL_ID,'nbc_parameters':asdict(model.nbc_parameters),
        'nhe1_carrier_amount_fmol':parent.inherited.CARRIER,'nkcc1_kinetics':asdict(model.nkcc1_kinetics),
        'minima':minima,'maxima':maxima,'max_abs_conservation_residuals':residual_max,
        'conservation_tolerances':dict(CONSERVATION_RESIDUAL_TOLERANCES),'max_conservation_ratios':ratio_max,
        'secretion':secretion,'chloride_partition':partition,
        'software_versions':{'numpy':np.__version__,'scipy':scipy.__version__},
        'thread_pools':[{k:v for k,v in p.items() if k in ('internal_api','prefix','num_threads','version')} for p in pools],
        'source_fields_note':frozen_wt['source_fields_note']}
    write_json(OUT/(case+'_verification.json'),verification)
    if status=='PASS': budget['completed_cases'].append(case)
    budget['numerical_execution_seconds']+=time.perf_counter()-START
    budget['status']=('TASK41_RECRUITMENT_VALIDATION_COMPLETE' if case=='ae4_null' and status=='PASS'
        else case.upper()+'_PASS' if status=='PASS' else status)
    budget['no_further_scientific_execution']=case=='ae4_null' or status!='PASS'
    budget['limits_respected']=(budget['integration_attempts']<=4 and budget['numerical_retries']<=1
        and len(budget['started_cases'])<=3 and budget['numerical_execution_seconds']<=900)
    write_json(OUT/'budget.json',budget)
    log={'case':case,'status':status,'attempts':attempts,'first_failure':first_failure,
         'error':error,'secretion':secretion,'chloride_partition':partition,'budget':budget}
    write_json(HERE/(case+'_execution.log'),log)
    print(json.dumps(log,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',choices=list(CASES))
    run_case(parser.parse_args().case)
