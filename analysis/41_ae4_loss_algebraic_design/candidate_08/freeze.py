"""Verify exact WT nesting and freeze its inherited rest; no new root solve."""
from candidate8_common import *
from threadpoolctl import threadpool_info
from datetime import datetime,timezone
import importlib.util,io,time,unittest


def main():
    clock=time.perf_counter();assert not (OUT/'budget.json').exists()
    spec=importlib.util.spec_from_file_location('recruitment_tests',ROOT/'tests/test_task41_ae4_cacc_recruitment.py')
    tests=importlib.util.module_from_spec(spec);spec.loader.exec_module(tests)
    stream=io.StringIO();result=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(tests))
    (HERE/'source_tests.log').write_text(stream.getvalue());print(stream.getvalue())
    if not result.wasSuccessful():raise RuntimeError('Source verification failed before production')
    oldrest,oldstim,rest,stim,y,reference=models_and_reference()
    from modern_full_model.ae4_cacc_recruitment import recruitment_factor
    from dataclasses import replace
    for expression in (1.,.05,0.):
        g=replace(WT,ae4_expression=expression)
        test=stim.evaluate(1e-6,y,genotype=g)
        r=stim.cacc_recruitment.null_stimulated_fraction
        assert test.diagnostics.regulatory['ae4_dependent_cacc_recruitment_factor']==recruitment_factor(1.,expression,r)
        assert all(abs(test.diagnostics.conservation_residuals[k])<=v for k,v in CONSERVATION_RESIDUAL_TOLERANCES.items())
        assert test.diagnostics.water==oldstim.evaluate(1e-6,y,genotype=g).diagnostics.water
    e=rest.evaluate(0.,y);np.testing.assert_array_equal(e.rhs,oldrest.evaluate(0.,y).rhs)
    row,failures,ratios=diagnose(rest,0.,y,e);assert not failures
    assert max(abs(e.rhs[:12]))<1e-10
    assert all(p['num_threads']==1 for p in threadpool_info())
    record={**reference,'observables':row,'diagnostics':asdict(e.diagnostics),
        'rest_provenance':'Bitwise identical Task 40 WT equations and frozen state; no new solve',
        'task41_candidate8_rest_solves':0,'status':'PASS'}
    write_json(OUT/'wt_rest.json',record)
    design=read_json(OUT/'design.json')
    write_json(OUT/'budget.json',{'status':'WT_REST_FROZEN','source_tests_pass':True,'wt_rest_pass':True,
        'stationary_solves':0,'algebraic_design_solver_calls':design['root_calls'],
        'algebraic_design_residual_calls':design['actual_residual_calls'],
        'selected_production_parameters':asdict(stim.cacc_recruitment),'parameter_selection_uses_target':True,
        'integration_attempts':0,'numerical_retries':0,'started_cases':[],'completed_cases':[],
        'numerical_execution_seconds':time.perf_counter()-clock+design['wall_seconds'],
        'no_further_scientific_execution':False,'scientific_workers':1,'blas_threads':1,
        'genotype_specific_rest_solves':0})
    files=[ROOT/'src/modern_full_model/ae4_cacc_recruitment.py',ROOT/'tests/test_task41_ae4_cacc_recruitment.py',*HERE.glob('*.py'),OUT/'design.json']
    write_json(OUT/'frozen_inputs.json',{'prepared_head':PREPARED_HEAD,'branch':BRANCH,
        'frozen_utc':datetime.now(timezone.utc).isoformat(),'frozen_rest_state_sha256':record['state_sha256'],
        'selected_parameters':asdict(stim.cacc_recruitment),'stimulus_parameter_hashes':parameter_hashes(stim),
        'rest_parameter_hashes':parameter_hashes(rest),'execution_files':[{'path':str(p.relative_to(ROOT)),
            'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)],
        'long_run_validity_not_established':True})
    print(json.dumps({'status':'PASS','rest_sha256':record['state_sha256'],'new_rest_solves':0,'selected_b':design['null_stimulated_cacc_fraction']},indent=2))

if __name__=='__main__':main()
