"""Verify saved outputs without importing the model or assuming outcomes."""
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/42_catalan_2025_ae4_mechanism_classes'
PANEL = ('C0', 'C1', 'C2', 'C3a', 'C3b', 'C4a', 'C4b')
CASES = ('wt', 'ae4_5pct', 'ae4_null')


def read(path):
    return json.loads(path.read_text())


def check_files(entries):
    for item in entries:
        assert hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest() == item['sha256'], item['path']


def main():
    complete = read(OUT/'panel_complete.json')
    assert complete['status'] == 'ALL_PREDECLARED_CLASSES_RECORDED'
    assert complete['stationary_solves'] == 7
    assert complete['genotype_rest_solves'] == complete['optimisation_calls'] == complete['parameter_sweeps'] == 0
    assert not complete['interpretation_unlocked']
    preflight = read(OUT/'preflight.json')
    check_files(preflight['execution_files'])
    rows, hashes = [], None
    attempts = 0
    for name in PANEL:
        directory = OUT/name
        rest, budget = read(directory/'wt_rest.json'), read(directory/'budget.json')
        assert rest['initial_state_vector'] == preflight['reference_state_vector']
        assert budget['stationary_solves'] == 1 and budget['stationary_numerical_retries'] <= 1
        assert budget['workers'] == budget['blas_threads'] == 1
        assert budget['numerical_retries'] <= 1
        assert budget['no_further_scientific_execution']
        if not rest['admissible']:
            assert rest['status'] != 'PASS' and rest['failed_gates']
            assert budget['integration_attempts'] == 0 and not budget['started_cases']
            rows.append({'class': name, 'case': 'rest', 'status': rest['status']})
            continue
        assert rest['status'] == 'PASS' and all(rest['gates'].values())
        assert budget['started_cases'] == list(CASES)
        freeze = read(directory/'frozen_inputs.json')
        check_files(freeze['task42_execution_files'])
        if hashes is None: hashes = freeze['stimulus_parameter_hashes']
        assert hashes == freeze['stimulus_parameter_hashes']
        for case, expression in zip(CASES, (1., .05, 0.)):
            v = read(directory/(case+'_verification.json'))
            attempts += len(v['integration_attempts'])
            assert v['initial_state_vector'] == rest['state_vector']
            assert v['initial_state_sha256'] == rest['state_sha256'] == freeze['frozen_rest_state_sha256']
            assert v['parameter_hashes'] == hashes
            assert v['only_genotype_difference'] == ({} if case=='wt' else {'ae4_expression': expression})
            assert v['ae4_routing_model'] == 'CATALAN_2025_SOURCE_'+name
            assert all(p['num_threads']==1 for p in v['thread_pools'])
            if v['integration_attempts']:
                assert v['integration_attempts'][0]['solver']['method'] == 'Radau'
                assert len(v['integration_attempts']) <= 2
                if len(v['integration_attempts']) == 2:
                    assert v['integration_attempts'][0]['solver_status'] == 'failed'
                    assert v['integration_attempts'][1]['solver']['method'] == 'BDF'
            if v['status'] == 'PASS':
                assert v['complete_600_s'] and v['first_failure'] is None
                assert max(v['max_conservation_ratios'].values()) <= 1
                assert v['maxima']['source_signature_residual_fmol_s'] == 0.
                assert v['maxima']['cell_balance_residual_fmol_s'] <= 1e-10
                with (directory/(case+'_all_window_integrals.csv')).open() as f:
                    integrals = list(csv.DictReader(f))
                for quantity in {r['quantity'] for r in integrals}:
                    windows = {(int(r['window_start_s']),int(r['window_end_s'])):float(r['value'])
                               for r in integrals if r['quantity']==quantity}
                    assert math.isclose(windows[0,60]+windows[60,600], windows[0,600], rel_tol=2e-14, abs_tol=2e-14)
                total = next(float(r['value']) for r in integrals if r['quantity']=='water_outflow' and r['window_start_s']=='0' and r['window_end_s']=='600')
                assert math.isclose(total,v['secretion']['cumulative_0_600_pL'],rel_tol=1e-15)
                coefficients = next(c['source_coefficients_na_k_cl_tic_ta'] for c in preflight['classes'] if c['class_id']==name)
                for start,end in ((0,60),(60,600),(0,600)):
                    values = {r['quantity']:float(r['value']) for r in integrals if int(r['window_start_s'])==start and int(r['window_end_s'])==end}
                    for field,nu in zip(('ae4_na_cell_source','ae4_k_cell_source','ae4_signed_cl_flux','ae4_tic_cell_source','ae4_ta_cell_source'),coefficients):
                        assert math.isclose(values[field],nu*values['ae4_signed_cl_flux'],rel_tol=2e-14,abs_tol=2e-14)
            else:
                activation_failure = (case=='wt' and v['complete_600_s'] and v['secretion'] is not None and not all(v['secretion']['activation_checks'].values()))
                solver_failure = any(a['solver_status']=='failed' for a in v['integration_attempts'])
                assert v['first_failure'] is not None or v['solver_exception'] or activation_failure or solver_failure
            rows.append({'class': name, 'case': case, 'status': v['status'],
                'checked_states': v['checked_state_count'], 'last_checked_time_s': v['last_checked_time_s']})
    assert attempts == complete['integration_attempts']
    report = {'status': 'SAVED_PREDICTIONS_VERIFIED', 'predictions': rows,
        'stationary_solves': complete['stationary_solves'], 'production_integrations': attempts,
        'source_test_count': 7, 'numerical_retries': complete['numerical_retries'],
        'all_execution_hashes_match_preflight': True, 'interpretation_unlocked': False,
        'checks': ['All seven classes recorded, including failures', 'Each class starts its rest solve from Task 40',
                   'Genotypes share their class WT rest and unchanged scientific parameters',
                   'Only AE4 expression differs between genotypes', 'Inherited gates and retry limits preserved',
                   'Successful window integrals and source identities reconcile', 'Incomplete trajectories retain failure records']}
    (OUT/'saved_prediction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    passed = sum(r['status']=='PASS' and r['case']!='rest' for r in rows)
    failed = sum(r['status']!='PASS' and r['case']!='rest' for r in rows)
    print(f'SAVED_PREDICTIONS_VERIFIED: {passed} successful trajectories, {failed} recorded failures.')


if __name__ == '__main__': main()
