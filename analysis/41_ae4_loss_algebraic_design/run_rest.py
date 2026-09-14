"""Source-test, solve one amended WT rest, and freeze before trajectories."""
from validation_common import *
from reduced_inverse import water_eliminated_state,COORD_ROWS,ION_ROWS,ION_SCALES
from inverse_balance import COORDINATE_LOWER,COORDINATE_UPPER,expit,root
from threadpoolctl import threadpool_info
from datetime import datetime,timezone
import importlib.util
import io
import time
import unittest


def main():
    clock=time.perf_counter()
    assert not (OUT/'budget.json').exists()
    spec=importlib.util.spec_from_file_location('task41_tests',ROOT/'tests/test_task41_nkcc1_capacity.py')
    tests=importlib.util.module_from_spec(spec);spec.loader.exec_module(tests)
    stream=io.StringIO();tested=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(tests))
    (HERE/'source_tests.log').write_text(stream.getvalue());print(stream.getvalue())
    if not tested.wasSuccessful():raise RuntimeError('Source tests failed; no rest solve')
    _,_,m,stim,y0,reference=models_and_reference()
    pools=threadpool_info();assert all(p['num_threads']==1 for p in pools)
    lo=COORDINATE_LOWER[COORD_ROWS];hi=COORDINATE_UPPER[COORD_ROWS]
    initial=np.array(reference['coordinates'])[COORD_ROWS];z0=np.log((initial-lo)/(hi-initial))
    calls=0
    def equations(z):
        nonlocal calls
        calls+=1
        y,_=water_eliminated_state(m,lo+(hi-lo)*expit(z),activation=0.)
        return m.rhs(0.,y)[ION_ROWS]/ION_SCALES
    fit=root(equations,z0,method='hybr',options={'xtol':1e-11,'maxfev':1900})
    y,x=water_eliminated_state(m,lo+(hi-lo)*expit(fit.x),activation=0.)
    e=m.evaluate(0.,y);row,failures,ratios=diagnose(m,0.,y,e)
    if not fit.success:failures.append({'gate':'stationary_solver_success'})
    if np.max(np.abs(e.rhs[ION_ROWS]/ION_SCALES))>1e-7:failures.append({'gate':'stationary_residual'})
    if np.max(np.abs(e.rhs[[4,10]]))>1e-9:failures.append({'gate':'omitted_charge_rows'})
    if np.max(np.abs(e.rhs[[5,11]]))>1e-12:failures.append({'gate':'stationary_water'})
    record={'status':'PASS' if not failures else 'TASK41_WT_REST_FAILED','admissible':not failures,
        'state_vector':y.tolist(),'state_names':list(m.state_names),'coordinates':x.tolist(),
        'state_sha256':sha256_object(y.tolist()),'core_state_sha256':sha256_object(y[:12].tolist()),
        'sole_initial_guess_task40_state_sha256':reference['state_sha256'],
        'initial_state_vector':y0.tolist(),'solver_success':bool(fit.success),'solver_message':str(fit.message),
        'actual_residual_calls':calls,'rest_solve_calls':1,'observables':row,'rhs':e.rhs.tolist(),
        'failed_gates':failures,'conservation_ratios':ratios,'diagnostics':asdict(e.diagnostics)}
    write_json(OUT/'wt_rest.json',record)
    designs=[read_json(OUT/p) for p in ('capacity_design.json','constrained_design.json','constrained_design_retry.json',
        'constrained_design_scaled.json','reduced_capacity_design.json','two_parameter_design.json','two_parameter_design_scaled.json')]
    budget={'status':'WT_REST_FROZEN' if not failures else record['status'],
        'source_tests_pass':True,'wt_rest_pass':not failures,'stationary_solves':1,
        'algebraic_design_solver_calls':len(designs),'algebraic_least_squares_solver_calls':2,
        'algebraic_design_residual_calls':sum(d['actual_residual_calls'] for d in designs),
        'parameter_selection_uses_target':True,'selected_production_parameters':asdict(m.nkcc1_capacity),
        'integration_attempts':0,'numerical_retries':0,'started_cases':[],'completed_cases':[],
        'numerical_execution_seconds':time.perf_counter()-clock+sum(d['wall_seconds'] for d in designs),
        'no_further_scientific_execution':bool(failures),'scientific_workers':1,'blas_threads':1,
        'unselected_algebraic_candidates_simulated':0,'genotype_specific_rest_solves':0}
    write_json(OUT/'budget.json',budget)
    files=[ROOT/'src/modern_full_model/nkcc1_capacity_limited.py',ROOT/'tests/test_task41_nkcc1_capacity.py',*HERE.glob('*.py'),OUT/'capacity_design.json']
    write_json(OUT/'frozen_inputs.json',{'prepared_head':PREPARED_HEAD,'branch':BRANCH,
        'frozen_utc':datetime.now(timezone.utc).isoformat(),'frozen_rest_state_sha256':record['state_sha256'],
        'selected_parameters':asdict(m.nkcc1_capacity),'stimulus_parameter_hashes':parameter_hashes(stim),
        'rest_parameter_hashes':parameter_hashes(m),'execution_files':[{'path':str(p.relative_to(ROOT)),
            'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)],
        'long_run_design_null_physiology_failure_explicit':True})
    print(json.dumps({'status':record['status'],'failed_gates':failures,'rest':row,'budget':budget},indent=2))


if __name__=='__main__':main()
