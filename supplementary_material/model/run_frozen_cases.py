"""Execute only the eight immutable Task52 cases, once, at gated checkpoints."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import time
import traceback
import numpy as np
from scipy.integrate import Radau, BDF
from threadpoolctl import threadpool_limits, threadpool_info
from task52_model import PairedModel, COUNTS, dependency_audit
from diagnostics import diagnose, flux_values

TASK = Path(__file__).resolve().parent
ROOT = TASK.parents[1]
OUT = TASK/'output'
EPS = 1e-6
GENOTYPES = ('WT','KO')
FLUXES = ('N','E','A','H','B','P','J','J_tmem16a','J_aux','K_apical','K_basolateral','J_para_Cl','C','Q')
RULES = {n:np.polynomial.legendre.leggauss(n) for n in (4,8)}

def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()

def verify_inputs(stage):
    assert git('branch','--show-current') == 'analysis/task-52-chloride-reservoir-final-test'
    prior = {'52D':'52C','52E':'52D'}[stage]
    publication = json.loads((OUT/f'publication_{prior}.json').read_text())
    assert publication['pass'] and publication['commit_sha'] == git('rev-parse','HEAD')
    freeze = json.loads((OUT/'parameter_and_case_freeze.json').read_text())
    for path, expected in freeze['immutable_artifact_sha256'].items():
        assert sha(ROOT/path) == expected, f'Frozen input changed: {path}'
    prefix=freeze['append_only_ledger_prefix']
    data=(ROOT/prefix['path']).read_bytes()
    assert hashlib.sha256(data[:prefix['length_bytes']]).hexdigest()==prefix['sha256']
    receipt = json.loads((OUT/'checkpoints/52C.json').read_text())
    fp = str((OUT/'parameter_and_case_freeze.json').relative_to(ROOT))
    assert sha(OUT/'parameter_and_case_freeze.json') == receipt['artifacts'][fp]
    assert dependency_audit()['pass']
    assert freeze['case_matrix']['duration_s'] == 600.0
    assert freeze['case_matrix']['solver']['integration_start_s']==EPS
    return freeze

def flux_array(model, pair):
    values=[flux_values(model,e) for e in pair]
    return np.array([[f[k] for k in FLUXES] for f in values])

def quadrature(model, dense, a, b, order):
    nodes, weights = RULES[order]
    total = np.zeros((2,len(FLUXES)))
    for x,w in zip(nodes,weights):
        t = (a+b)/2+(b-a)*x/2
        total += w*flux_array(model,model.evaluate(t,dense(t)))
    return total*(b-a)/2

def deficit(wt,ko):
    return float(1-ko/wt) if wt > 0 else None

def run_attempt(case, model, y0, method, directory, settings):
    started = time.monotonic()
    counts0 = dict(COUNTS)
    samples = {}
    max_ratios = {}
    max_abs = {}
    failure = None
    secondary_errors = []
    status = 'running'
    solver = None
    integrals = {n:np.zeros((2,len(FLUXES))) for n in RULES}
    integral_times = [EPS]
    integral_history = [integrals[8].copy()]
    accepted_times=[EPS]
    accepted_states=[y0.copy()]
    last_valid_t = 0.0

    def observe(t,y,kind):
        nonlocal failure,last_valid_t
        pair = model.evaluate(t,y)
        rows=[]; faults=[]
        for i,(g,e) in enumerate(zip(GENOTYPES,pair)):
            row,bad,ratios = diagnose(model,t,y[13*i:13*(i+1)],e)
            rows.append(row)
            faults.extend(dict(genotype=g,**b) for b in bad)
            for k,v in ratios.items(): max_ratios[g+'__'+k]=max(max_ratios.get(g+'__'+k,0.),v)
            for k,v in row.items():
                if 'residual' in k or k in ('cell_charge_fmol','independent_cell_assembly_max_abs_fmol_s'):
                    max_abs[g+'__'+k]=max(max_abs.get(g+'__'+k,0.),abs(v))
        fw,fk = [flux_values(model,e) for e in pair]
        for k in ('N','E'):
            if fk[k] != case['ko_supply_multiplier']*fw[k]:
                faults.append({'gate':'exact_imposed_'+k})
        if fk['A'] != 0.0: faults.append({'gate':'exact_zero_KO_AE4'})
        for k in ('G_aux_effective_S','G_tmem16a_effective_S'):
            if rows[0][k] != rows[1][k]: faults.append({'gate':'shared_'+k})
        previous=samples.get(float(t))
        samples[float(t)] = {'rows':rows,'state':y.copy(), 'kind':kind | (previous['kind'] if previous else 0)}
        if faults:
            samples[float(t)]['kind'] |= 8
            failure = {'time_s':float(t),'failures':faults}
            return False
        last_valid_t = float(t)
        return True

    try:
        if not observe(0.,y0,4|2) or not observe(EPS,y0,4|2):
            status='physical_or_conservation_failure'
        else:
            atol=np.tile([settings['amount_atol_fmol']]*5+[settings['volume_atol_pL']]+
                         [settings['amount_atol_fmol']]*5+[settings['volume_atol_pL'],settings['regulatory_atol']],2)
            solver = {'Radau':Radau,'BDF':BDF}[method](model.rhs,EPS,y0.copy(),600.,
                       rtol=settings['rtol'],atol=atol,max_step=settings['max_step_s'])
            while solver.status == 'running':
                a=float(solver.t)
                message=solver.step()
                if solver.status == 'failed':
                    status='solver_status_failure';failure={'time_s':float(solver.t),'message':message};break
                b=float(solver.t);dense=solver.dense_output()
                accepted_times.append(b);accepted_states.append(solver.y.copy())
                points={b:1}
                for k in range(int(np.floor(a))+1,int(np.floor(b))+1):
                    if a<k<=b: points[float(k)]=points.get(float(k),0)|2
                limit=b
                for t,kind in sorted(points.items()):
                    y=solver.y.copy() if t==b else dense(t)
                    if not observe(t,y,kind):
                        limit=t;status='physical_or_conservation_failure';break
                # Commit both rules atomically so a quadrature exception cannot
                # silently advance one integral or erase an earlier gate failure.
                increments={n:quadrature(model,dense,a,limit,n) for n in RULES}
                for n in RULES: integrals[n] += increments[n]
                integral_times.append(limit);integral_history.append(integrals[8].copy())
                if status != 'running': break
            if status=='running': status='completed' if solver.status=='finished' else 'solver_status_failure'
    except Exception as exc:
        error={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
        if failure is None:
            status='exception';failure=error
        else:
            secondary_errors.append(error)

    times=np.array(sorted(samples))
    if not len(times):
        raise RuntimeError('No onset state could be evaluated; original exception: '+str(failure))
    columns=tuple(samples[times[0]]['rows'][0])
    assert all(tuple(r)==columns for s in samples.values() for r in s['rows'])
    values=np.array([[list(r.values()) for r in samples[t]['rows']] for t in times])
    states=np.array([samples[t]['state'] for t in times])
    kinds=np.array([samples[t]['kind'] for t in times])
    index={k:i for i,k in enumerate(columns)}
    T=float(times[-1]);delta=states[:,2]-states[:,15]
    # The frozen epsilon convention holds the onset state fixed. Report its
    # independent right-limit flux quadrature explicitly; never hide it in delta.
    onset=EPS*flux_array(model,model.evaluate(EPS,y0))
    primary=integrals[8]+onset
    diff=primary[0]-primary[1]
    cross=integrals[8]-integrals[4]
    at=lambda k:float(diff[FLUXES.index(k)])
    onsetdiff=onset[0]-onset[1]
    onset_bound=float(sum(abs(onsetdiff[FLUXES.index(k)])*factor for k,factor in (('J',1),('N',2),('A',1),('E',1))))
    rhs=float(delta[0]+2*at('N')+at('A')+at('E')-delta[-1])
    identity_error=at('J')-rhs
    cross_terms={k:float(cross[0,FLUXES.index(k)]-cross[1,FLUXES.index(k)]) for k in ('J','N','A','E')}
    cross_cells={g:{k:float(cross[i,FLUXES.index(k)]) for k in ('J','N','A','E')} for i,g in enumerate(GENOTYPES)}
    cross_terms['two_N']=2*cross_terms['N']
    for g in GENOTYPES: cross_cells[g]['two_N']=2*cross_cells[g]['N']
    quadrature_coverage_complete=(integral_times[-1]==T)
    budget_pass=(quadrature_coverage_complete and abs(identity_error)<=1e-5 and onset_bound<=1e-5 and
                 max(abs(v) for v in cross_terms.values())<=1e-6 and
                 max(abs(v) for d in cross_cells.values() for v in d.values())<=1e-6)
    available=float(delta[0]+2*at('N')+at('A')+at('E'))
    budget={'T_s':T,'delta_0_fmol':float(delta[0]),'delta_T_fmol':float(delta[-1]),
        'quadrature_end_s':float(integral_times[-1]),'quadrature_coverage_complete':quadrature_coverage_complete,
        'D_J_fmol':at('J'),'D_N_cycles_fmol':at('N'),'two_D_N_fmol':2*at('N'),
        'D_A_fmol':at('A'),'D_E_fmol':at('E'),'RHS_fmol':rhs,
        'identity_error_fmol':identity_error,'absolute_tolerance_fmol':1e-5,
        'quadrature_8_minus_4_differential_fmol':cross_terms,
        'quadrature_8_minus_4_each_cell_fmol':cross_cells,
        'quadrature_tolerance_fmol':1e-6,'onset_held_state_interval_s':[0.,EPS],
        'onset_right_limit_integrals':{g:dict(zip(FLUXES,onset[i].tolist())) for i,g in enumerate(GENOTYPES)},
        'onset_identity_error_bound_fmol':onset_bound,
        'available_reservoir_plus_differential_supply_fmol':available,
        'released_fraction':at('J')/available if available!=0 else None,'pass':bool(budget_pass)}
    if not budget_pass and status=='completed': status='reservoir_audit_failure'
    grid=(kinds&2)!=0
    grid[-1]=True  # failed/end point included only in explicitly partial summary
    gt=times[grid];gv=values[grid]
    def window(a,b):
        mask=(gt>=a)&(gt<=b);tt=gt[mask];vv=gv[mask]
        if len(tt)<2 or tt[0]!=a or tt[-1]!=b: return None
        result={'window_s':[float(a),float(b)]}
        for name,key in (('fluid_pL','q_out_pL_s'),('chloride_fmol','J_fmol_s')):
            q=np.trapz(vv[:,:,index[key]],tt,axis=0)
            result[name]={g:float(q[i]) for i,g in enumerate(GENOTYPES)}
            result[name]['KO_deficit_fraction']=deficit(q[0],q[1])
            result['mean_'+name+'_per_s']={g:float(q[i]/(b-a)) for i,g in enumerate(GENOTYPES)}
        return result
    partial=window(0.,T) if T>0 else None
    headline=window(0.,600.) if status=='completed' else None
    broad=window(60.,600.) if status=='completed' else None
    observations={str(t):window(0.,float(t)) for t in range(60,601,60) if t<=T}
    trapdiff={k:float(np.trapz(gv[:,0,index[k+'_fmol_s']]-gv[:,1,index[k+'_fmol_s']],gt)) for k in ('J','N','A','E')}
    budget['headline_trapezoidal_identity_error_fmol']=float(trapdiff['J']-(delta[0]+2*trapdiff['N']+trapdiff['A']+trapdiff['E']-delta[-1]))
    np.savez_compressed(directory/f'{method}_trajectory.npz',time_s=times,
        columns=np.array(columns),observations=values,paired_states=states,sample_kind=kinds,
        delta_chloride_fmol=delta,quadrature_time_s=np.array(integral_times),
        quadrature_flux_names=np.array(FLUXES),quadrature8_cumulative=np.array(integral_history),
        raw_accepted_time_s=np.array(accepted_times),raw_accepted_paired_states=np.array(accepted_states),
        raw_accepted_validated=np.array(accepted_times)<=last_valid_t)
    endpoint={g:{k:float(values[-1,i,index[k]]) for k in columns} for i,g in enumerate(GENOTYPES)}
    extrema={g:{k:{'min':float(values[:,i,index[k]].min()),'max':float(values[:,i,index[k]].max())}
        for k in ('ph_i','cl_i_mM','tic_i_mM','na_i_mM','k_i_mM','volume_i_pL','E_Cl_V','V_a_V',
                  'chloride_driving_force_V','q_out_pL_s','J_fmol_s','J_aux_fmol_s','B_fmol_s','H_fmol_s')}
        for i,g in enumerate(GENOTYPES)}
    summary={'case':case,'method':method,'status':status,'completed_600s':status=='completed',
        'first_failure':failure,'last_valid_sample_s':last_valid_t,'recorded_end_s':T,
        'secondary_errors':secondary_errors,
        'sample_count':len(times),'accepted_steps_recorded':len(integral_times)-1,
        'raw_solver_accepted_steps':len(accepted_times)-1,
        'raw_endpoint_policy':'An accepted endpoint beyond the first failed dense sample is retained as unvalidated state only; no further integration or diagnostic evaluation.',
        'runtime_s':time.monotonic()-started,
        'evaluation_counts':{k:COUNTS[k]-counts0[k] for k in COUNTS},
        'solver_counters':{k:int(getattr(solver,k,0)) for k in ('nfev','njev','nlu')},
        'reservoir_budget':budget,'headline_0_600':headline,'broad_60_600':broad,
        'partial_0_T':partial,'cumulative_observations':observations,
        'endpoint_at_recorded_T':endpoint,'sampled_extrema':extrema,
        'max_inherited_conservation_tolerance_ratios':max_ratios,
        'max_absolute_residuals':max_abs,
        'gauss8_integrals_0_T':{g:dict(zip(FLUXES,primary[i].tolist())) for i,g in enumerate(GENOTYPES)} if quadrature_coverage_complete else None,
        'gauss8_partial_integrals':{'end_s':float(integral_times[-1]),
            'values':{g:dict(zip(FLUXES,primary[i].tolist())) for i,g in enumerate(GENOTYPES)}},
        'finite_trajectory_arrays':bool(np.isfinite(values).all() and np.isfinite(states).all()),
        'trajectory_sha256':sha(directory/f'{method}_trajectory.npz')}
    dump(directory/f'{method}_summary.json',summary)
    return summary

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['52D','52E']);args=parser.parse_args()
    freeze=verify_inputs(args.stage)
    stage_path=OUT/f'execution_{args.stage}.json'
    assert not stage_path.exists(),'Stage already attempted; no rerun permitted'
    record={'stage':args.stage,'parent_commit':git('rev-parse','HEAD'),
        'freeze_sha256':sha(OUT/'parameter_and_case_freeze.json'),'cases':[], 'complete':False}
    dump(stage_path,record)
    with threadpool_limits(limits=1):
        record['threadpools']=threadpool_info()
        assert all(p['num_threads']==1 for p in record['threadpools'])
        for case in freeze['case_matrix']['cases']:
            if case['checkpoint']!=args.stage: continue
            directory=OUT/'cases'/case['id'];directory.mkdir(parents=True,exist_ok=False)
            dump(directory/'started.json',{'case':case,'freeze_sha256':record['freeze_sha256']})
            states=freeze['projected_onsets'][case['projection']]
            y0=np.r_[states['WT']['state'],states['KO']['state']]
            model=PairedModel(protocol=case['protocol'],G_aux_S=case['G_aux_S'],ko_supply_multiplier=case['ko_supply_multiplier'])
            attempts=[]
            result=run_attempt(case,model,y0,'Radau',directory,freeze['case_matrix']['solver']);attempts.append('Radau')
            if result['status']=='solver_status_failure':
                result=run_attempt(case,model,y0,'BDF',directory,freeze['case_matrix']['solver']);attempts.append('BDF')
            entry={'id':case['id'],'status':result['status'],'selected_attempt':attempts[-1],
                   'attempts':attempts,'summary_sha256':sha(directory/f'{attempts[-1]}_summary.json')}
            record['cases'].append(entry);dump(stage_path,record)
            print(json.dumps({'stage':args.stage,**entry,'T_s':result['recorded_end_s'],
                              'mass_budget_pass':result['reservoir_budget']['pass']}),flush=True)
        record['complete']=True;record['evaluation_counts']=dict(COUNTS)
        record['dependency_audit']=dependency_audit();assert record['dependency_audit']['pass']
        dump(stage_path,record)

if __name__=='__main__': main()
