"""One recorded source-test preflight; one literal correction at most."""
import time
START=time.perf_counter()
from validation_common import *
import io
import unittest
from threadpoolctl import threadpool_info


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    path=OUT/'budget.json'
    if path.exists():
        budget=json.loads(path.read_text())
        assert budget['status']=='SOURCE_CORRECTION_ALLOWED'
        assert budget['focused_test_runs']==1 and budget['stationary_solves']==0
        budget['implementation_correction_rounds']=1
    else:
        budget={'status':'SOURCE_PREFLIGHT_RUNNING','focused_test_runs':0,
            'implementation_correction_rounds':0,'stationary_solves':0,
            'stationary_numerical_retries':0,'stationary_residual_evaluations':0,
            'integration_attempts':0,'numerical_retries':0,
            'started_cases':[],'completed_cases':[],
            'optimisation_calls':0,'parameter_sweeps':0,'calibration_calls':0,
            'genotype_specific_rest_solves':0,'scientific_workers':1,'blas_threads':1,
            'maximum_numerical_execution_seconds':900,'numerical_execution_seconds':0.,
            'source_tests_pass':False,'wt_rest_pass':False,'no_further_scientific_execution':False}
    budget['focused_test_runs']+=1
    write_json(path,budget)
    spec=importlib.util.spec_from_file_location('task40_source_tests',ROOT/'tests/test_task40_equal_cation_routing.py')
    tests=importlib.util.module_from_spec(spec);spec.loader.exec_module(tests)
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(tests))
    log=stream.getvalue();print(log)
    (HERE/f"focused_test_{budget['focused_test_runs']}.log").write_text(log)
    passed=result.wasSuccessful()
    pools=threadpool_info();assert pools and all(p['num_threads']==1 for p in pools)
    audit={'status':'PASS' if passed else 'FAIL','tests_run':result.testsRun,
        'failed_test_names':[str(t) for t,_ in result.failures+result.errors],
        'thread_pools':[{k:p[k] for k in ('internal_api','num_threads','version')} for p in pools]}
    if passed:
        _,old,_,new,stim,y=models_and_reference()
        audit['reference_preflight']=tests.reference_preflight(old,new,y)
        audit['rest_parameter_hashes']=parameter_hashes(new)
        audit['stimulus_parameter_hashes']=parameter_hashes(stim)
        audit['nonrouting_parent_files_verified']=len(verify_parent_sources()['files'])
        audit['old_donor_weighted_evaluator_retained']=True
        audit['independent_bookkeeping']=parent.inherited.bookkeeping(new,new.evaluate(0.,y))
    write_json(OUT/'source_verification.json',audit)
    budget['source_tests_pass']=passed
    budget['status']='SOURCE_PREFLIGHT_PASS' if passed else ('SOURCE_CORRECTION_ALLOWED' if budget['focused_test_runs']==1 else 'TASK40_EQUAL_ROUTING_SOURCE_FAILED')
    budget['no_further_scientific_execution']=budget['status']=='TASK40_EQUAL_ROUTING_SOURCE_FAILED'
    budget['numerical_execution_seconds']+=time.perf_counter()-START
    write_json(path,budget)
    print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
