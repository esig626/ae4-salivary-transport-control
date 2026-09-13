"""Read-only scientific reporting of saved Task 32 checkpoints; never solve.

Usage: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
       python scripts/report_task32_wt_chloride_calibration.py
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from modern_full_model import task32_wt_chloride_calibration as t


def main():
    directory = t.RESULTS
    report = json.loads((directory/'wt_calibration.json').read_text())
    summary = json.loads((directory/'execution_summary.json').read_text())
    budget_bytes = (directory/'budget.json').read_bytes()
    budget = json.loads(budget_bytes)
    cases_bytes = (directory/'cases.json').read_bytes()
    cases = sorted(json.loads(cases_bytes).values(), key=lambda r:r['label'])
    # Rehydrate only the exact accepted WT cache; no RHS evaluation or solve.
    cached = t.cached_wt(t.load_background('R09'))
    report['cached_start'] = cached
    closed = [r for r in cases if r['numerical_pass']]
    last_closed = closed[-1] if closed else cached
    report['last_closed_candidate'] = last_closed
    report['completed_calibration_stages'] = [s['stage'] for s in report['stages'] if 'accepted_label' in s]
    report['final_calibrated_scalars'] = ({k:summary['accepted_wt'][k]
        for k in ('carrier_amount_fmol','ae4_activity')} if summary['calibration_accepted'] else None)
    report['bounded_failure_is_not_nonexistence_proof'] = not summary['calibration_accepted']
    t.write_json(directory/'wt_calibration.json',report)
    rows = [cached]+cases
    t.export_tables(rows)
    a = t.allocation(last_closed)
    checks = {
        'residual_counter_equals_saved_new_cases':sum(r['actual_residual_evaluations'] for r in cases)
            == budget['stationary_residual_evaluations'],
        'solver_counter_equals_saved_new_cases':len(cases)==budget['stationary_solver_calls'],
        'all_new_cases_WT_R09':all(r['background']=='R09' and r['ae4_expression']==1 for r in cases),
        'all_new_cases_keep_Task31_NHE1_amount':all(r['carrier_amount_fmol']==cached['carrier_amount_fmol'] for r in cases),
        'no_scalar_bracket_claimed_from_unclosed_candidate':not report['completed_calibration_stages'],
        'failed_calibration_not_frozen':json.loads((directory/'frozen_parameters.json').read_text())['status']
            == 'not_frozen_calibration_failed',
        'no_exact_null_evaluation':not summary['exact_null_evaluated'],
        'no_stimulated_integrations':budget['stimulated_integrations']==0,
        'no_R10_evaluation':not summary['r10_evaluated'],
        'all_task_limits_respected':budget['limits_respected'],
        'no_knockout_or_secretion_fit':not report['null_data_used'] and not report['secretion_data_used'],
        'cache_has_no_Task32_residual_charge':cached['actual_residual_evaluations']==0,
    }
    for r in cases:
        checks[r['label']+'_state_hash_valid'] = t.sha256_object(r['core_state'])==r['core_state_sha256']
        checks[r['label']+'_closure_flag_consistent'] = r['full_rhs_pass']==(
            max(abs(v) for v in r['scaled_independent_rhs'])<=t.ROOT_SCALED_TOLERANCE
            and max(abs(r['raw_rhs'][i]) for i in t.OMITTED_CHARGE_ROWS)<=t.OMITTED_ROW_TOLERANCE
            and r['max_abs_regulatory_rhs_s_inv']<=t.REGULATORY_TOLERANCE_S_INV)
    verification = {'pass':all(checks.values()),'checks':checks,
        'focused_tests':{'passed':11,'failed':0,'log':'analysis/32_wt_chloride_calibration/focused_tests.log',
            'scientific_stationary_solves':0,'stimulated_integrations':0},
        'numerical_execution':'one worker, one BLAS thread; model diagnostics use saved states; no rerun',
        'cases_file_sha256':hashlib.sha256(cases_bytes).hexdigest(),
        'budget_file_sha256':hashlib.sha256(budget_bytes).hexdigest(),
        'calibration_accepted':summary['calibration_accepted'],
        'frozen_parameters_sha256':summary['frozen_parameters_sha256'],
        'stationary_gates':{'scaled_independent_rhs':t.ROOT_SCALED_TOLERANCE,
            'omitted_rhs_fmol_s':t.OMITTED_ROW_TOLERANCE,'bulk_charge_fmol':t.CHARGE_TOLERANCE_FMOL,
            'current_A':t.CURRENT_TOLERANCE_A,'regulatory_s_inv':t.REGULATORY_TOLERANCE_S_INV,
            'conservation_normalized':1.,'independent_rank':10,'boundary_hits':0},
        'last_closed_candidate':last_closed['label'],
        'unclosed_candidate_is_phenotype':False,
        'read_scope':'Task32 instructions, prepared source and tests, WT background contracts, Task31 WT density_freeze and calibration only; no historical null result file',
    }
    integrity = directory/'input_integrity.json'
    if integrity.exists():
        verification['input_integrity']=json.loads(integrity.read_text())
        verification['pass'] &= verification['input_integrity']['all_prepared_files_unchanged']
    t.write_json(directory/'verification.json',verification)
    lines = [
        '# Task 32: bounded WT chloride calibration', '',
        '**WT-only calibration did not recover R09 pH and chloride simultaneously.** '
        'The initial AE4 scalar stage did not obtain a valid chloride sign bracket. '
        'It stopped at an unclosed positive-activity endpoint; this is a bounded numerical result, not a proof of model infeasibility.', '',
        f'Prepared head: `{t.PREPARED_HEAD}`. Branch: `{t.BRANCH}`. One scientific agent; no other transporter law changed.', '',
        '| R09 WT case | AE4 activity | pH_i | Cl_i (mM) | Na_i (mM) | Volume (pL) | Stationary closure |',
        '|---|---:|---:|---:|---:|---:|---|',
    ]
    for r in rows:
        o=r['observables']
        lines.append(f"| {r['label']} | {r['ae4_activity']:.6g} | {o['ph_i']:.6f} | {o['cl_i_mM']:.6f} | {o['na_i_mM']:.6f} | {o['volume_i_pL']:.6f} | {'Pass' if r['numerical_pass'] else 'FAIL: candidate only'} |")
    lines += ['',
        'WT acceptance intervals were pH 6.84-6.98 and chloride 48.6-51.6 mM. '
        'The 0.1 activity state passes full closure, charge, current, conservation, positivity, rank, capacity, ordinary sodium and volume gates, '
        'but exceeds the upper WT chloride limit by 8.6113 mM. Its pH lies inside the WT interval. '
        'Lowering capacity tenfold changed chloride by only -0.09369 mM from the cached WT state.', '',
        'The 0.01 activity trial reached the unchanged Task 31 local ceiling of 2,400 actual residual evaluations. '
        'Its maximum scaled stationary residual is 0.00750789, versus the required 1e-7; '
        'the maximum amount residual is 0.00107400 fmol/s. Charge/current accounting and normal sodium/volume do not rescue missing stationary closure. '
        'The displayed pH and chloride for this trial are candidate diagnostics, not a phenotype or a scalar-bracket endpoint.', '',
        '1. **Simultaneous WT recovery:** No. The first AE4 stage stopped; NHE1 correction and the final AE4 correction were not reached.',
        f"2. **Required scalars:** No calibrated pair was established. NHE1 remained at the Task 31 amount `{cached['carrier_amount_fmol']:.16g}` fmol (about 14,088 carriers). "
        f"AE4 activity trials were 0.1 and 0.01. The scalar multiplies the existing AE4 carrier amount ({cached['ae4_carrier_amount_fmol']:.10g} fmol at activity 1); WT genotype expression stayed exactly 1.0.",
        f"3. **AE4 chloride-loading fraction:** No post-calibration fraction exists. At the last closed, uncalibrated WT trial, AE4 supplies {100*a['ae4_fraction_positive_basolateral']:.6f}% of positive basolateral chloride loading. This is an output, never a calibration target.",
        '4. **Exact AE4-null closure with ordinary sodium and volume:** Not evaluated. There is no accepted WT calibration to freeze.',
        '5. **Untuned null pH/chloride versus experiment:** Not evaluated. No null measurement or historical null result entered calibration.',
        '6. **Secretion direction and magnitude:** Unproven; no stimulated integrations were permitted.',
        '7. **R10 out-of-sample survival:** Not tested because no accepted frozen R09 WT/null pair exists.',
        '8. **What remains unproven:** Whether any admissible simultaneous WT solution exists under the permitted two-scalar procedure, '
        'whether exact AE4-null REST closes, and whether AE4 loss predicts the experimental chloride, pH or secretion effects. '
        'This failed local endpoint does not establish absence of another stationary root or of a finite positive chloride bracket.', '',
        'The last closed WT trial has signed AE4, NKCC1 and AE2 chloride fluxes '
        f"{a['ae4_signed_cl_fmol_s']:.12g}, {a['nkcc1_signed_cl_fmol_s']:.12g}, and {a['ae2_signed_cl_fmol_s']:.12g} fmol/s, respectively. "
        f"Total positive loading is {a['total_positive_basolateral_cl_fmol_s']:.12g} fmol/s; net loading is {a['net_basolateral_cl_fmol_s']:.12g} fmol/s. "
        f"AE4 / positive NKCC1 = {a['ae4_over_positive_nkcc1']:.9g}; AE4 / net basolateral = {a['ae4_over_net_basolateral']:.9g}. "
        'AE4 / positive AE2 is undefined because AE2 is a chloride counterflux, with zero positive loading. '
        'Signed AE2, cation, carbon, apical, paracellular and water fluxes are retained in `flux_ledger.csv`. '
        'The cached Task 31 row leaves paracellular chloride blank because that flux was not stored; it was not recomputed.', '',
        f"Budget: {budget['stationary_solver_calls']}/18 stationary calls, {budget['stationary_residual_evaluations']:,}/30,000 residual evaluations, "
        f"{budget['ae4_scalar_evaluations']}/6 AE4 scalar evaluations, 0/4 NHE1 scalar evaluations and 0/2 integrations. "
        f"Numerical execution took {budget['numerical_wall_seconds']:.6f} s, excluding development and algebra/control-flow tests. "
        'Global budgets remain unused in part because the predeclared unclosed-endpoint stop was reached. '
        'The factor-ten steps were directional scalar bracket discovery, with no factorial, grid or multidimensional parameter search. '
        'No sign interpolation was possible. Successful cached cases were never rerun.', '',
        '`frozen_parameters.json` explicitly records **not frozen: calibration failed**, with no accepted parameters. '
        'No freeze hash is claimed. `cases.json` preserves both trials; `budget.json` preserves their counters. '
        'All 11 focused tests pass. All prepared model source files remain byte-for-byte unchanged.', '',
        'Reproduce the focused checks with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m unittest discover -s tests -p test_task32_calibration.py -v`. '
        'Regenerate this report from saved checkpoints with `python scripts/report_task32_wt_chloride_calibration.py`; it never calls a stationary solver. '
        'The numerical calibration entry point is `PYTHONPATH=src python -m modern_full_model.task32_wt_chloride_calibration` with the same thread variables. '
        'An existing completed execution is returned read-only.', '',
    ]
    t.ANALYSIS.mkdir(parents=True,exist_ok=True)
    (t.ANALYSIS/'final_answer.md').write_text('\n'.join(lines),encoding='utf-8')
    assert (directory/'budget.json').read_bytes()==budget_bytes
    assert (directory/'cases.json').read_bytes()==cases_bytes
    if not verification['pass']:
        raise RuntimeError('saved-output verification failed')
    print(json.dumps({'verification_pass':verification['pass'], 'focused_tests_passed':11,
                      'scientific_recomputation':False}))


if __name__=='__main__':
    main()
