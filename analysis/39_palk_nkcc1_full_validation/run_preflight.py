"""Ordered Task 39 source, algebraic scale and exact REST gates; no integration."""
import time
START = time.perf_counter()
from validation_common import *
from decimal import Decimal, localcontext
import io
import math
import unittest
import argparse
from datetime import datetime, timezone
from modern_full_model.nkcc1_palk2010 import palk_shape_factor


def main(*, literal_test_correction=False):
    OUT.mkdir(parents=True, exist_ok=True)
    previous_seconds = 0.
    if literal_test_correction:
        previous = json.loads((OUT/'budget.json').read_text())
        assert previous['integration_attempts']==0 and previous['implementation_correction_rounds']==0
        assert previous['focused_test_runs']==1 and not previous['source_tests_pass']
        assert previous['status']=='PALK_NKCC1_SOURCE_OR_CONSERVATION_FAILED'
        previous_seconds=previous['numerical_execution_seconds']
        # Preserve the initial coding-test failure, not a new scientific condition.
        for old,new in ((OUT/'source_verification.json',OUT/'initial_test_assertion.json'),
                        (HERE/'focused_test.log',HERE/'initial_test_assertion.log')):
            assert not new.exists()
            new.write_bytes(old.read_bytes())
    else:
        assert not (OUT / "budget.json").exists(), "Do not rerun the Task 39 gate sequence"
    budget = {"stationary_solves": 0, "optimisation_calls": 0, "parameter_sweeps": 0,
        "integration_attempts": 0, "numerical_retries": 0, "started_cases": [],
        "completed_cases": [], "scientific_workers": 1, "blas_threads": 1,
        "maximum_intended_integrations": 3, "maximum_integration_attempts": 4,
        "maximum_numerical_execution_seconds": 900,
        "numerical_execution_seconds": previous_seconds, "source_tests_pass": False,
        "rest_nesting_pass": False, "implementation_correction_rounds": int(literal_test_correction),
        "no_further_scientific_execution": False,
        "root_scope": "No stationary/resting root solves. Inherited acid-base speciation and NBC current-closure brentq evaluations are unchanged."}
    write_json(OUT / "budget.json", budget)
    algebra = resting_algebra()
    with localcontext() as ctx:
        ctx.prec = 60
        na,k,cl = map(lambda v: Decimal(str(v)), REST_CONCENTRATIONS)
        x = na*k*cl**2
        shape = (Decimal('157.5')-Decimal('2.0096e-5')*x)/(Decimal('1.0306')+Decimal('1.3852e-6')*x)
        decimal_audit = {"x_rest_mM4":str(x), "shape_rest":str(shape),
                         "alpha_eff_fmol_s":str(Decimal(str(REST_CYCLE))/shape)}
    expected = {"x_rest_mM4":4941065.334966328,"shape_rest":7.391062769440939,
                "alpha_eff_fmol_s":0.017334746894096052}
    checks = {key:math.isclose(algebra[key],v,rel_tol=2e-15) for key,v in expected.items()}
    checks["independent_decimal_agreement"] = all(math.isclose(algebra[k],float(v),rel_tol=2e-15)
        for k,v in decimal_audit.items())
    checks["implemented_shape"] = math.isclose(palk_shape_factor(*REST_CONCENTRATIONS),algebra["shape_rest"],rel_tol=2e-15)
    source = {"algebra":algebra,"expected":expected,"decimal_60_digit_audit":decimal_audit,
        "relative_tolerance":2e-15,"algebra_checks":checks,
        "source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/#APP1",
        "source_equation":"Eq. 17: (157.5 - 2.0096e-5 Na_i K_i Cl_i^2)/(1.0306 + 1.3852e-6 Na_i K_i Cl_i^2)",
        "conversion":"X_M = 1e-12 X_mM; 2.0096e7*1e-12 = 2.0096e-5; 1.3852e6*1e-12 = 1.3852e-6",
        "reversal_x_mM4":157.5/2.0096e-5, "inward_cycle_source_na_k_cl_tic_ta":[1,1,2,0,0]}
    if not all(checks.values()):
        status = "PALK_NKCC1_SOURCE_REPRODUCTION_FAILED"
    else:
        suite = unittest.defaultTestLoader.discover(str(ROOT/'tests'),pattern='test_task39_palk_nkcc1.py')
        log = io.StringIO()
        tests = unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
        (HERE/'focused_test.log').write_text(log.getvalue())
        source["tests_run"] = tests.testsRun
        source["tests_pass"] = tests.wasSuccessful()
        budget["focused_test_runs"] = 1+int(literal_test_correction)
        budget["source_tests_pass"] = tests.wasSuccessful()
        status = "PASS" if tests.wasSuccessful() else "PALK_NKCC1_SOURCE_OR_CONSERVATION_FAILED"
    source["status"] = status
    write_json(OUT/'source_verification.json',source)
    if status != 'PASS':
        budget.update(status=status,no_further_scientific_execution=True,
            numerical_execution_seconds=previous_seconds+time.perf_counter()-START)
        write_json(OUT/'budget.json',budget)
        print(json.dumps(source,indent=2)); return

    payload = {**algebra, "law":PALK_NKCC1, "rest_concentrations_mM":list(REST_CONCENTRATIONS),
        "coefficients_mM":{"a1":157.5,"a2":2.0096e-5,"a3":1.0306,"a4":1.3852e-6}}
    write_json(OUT/'alpha_eff_freeze.json',{"payload":payload,"sha256":sha256_object(payload),
        "derived_from":"Accepted Task 31 R09 WT REST only; literal division, no fit",
        "frozen_utc":datetime.now(timezone.utc).isoformat()})
    saved,task31,old_rest,old_stim,rest,stim,y0 = models_and_state()
    old = old_rest.evaluate(0.,y0,genotype=WT)
    direct31 = task31.evaluate(0.,y0,genotype=WT)
    new = rest.evaluate(0.,y0,genotype=WT)
    row,failures,ratios = diagnose(rest,0.,y0,new)
    delta = np.asarray(new.rhs)-np.asarray(old.rhs)
    vols = {name:float(getattr(new.diagnostics.membranes,name)-getattr(old.diagnostics.membranes,name))
        for name in ('v_apical_V','v_basolateral_V','v_transepithelial_V')}
    water_delta = row['q_out_pL_s']-old.diagnostics.water.lumen_outflow_pL_s
    actual_ci = new.diagnostics.observables.cell_concentrations_mM
    onset = stim.evaluate(1e-6,y0,genotype=WT)
    old_onset = old_stim.evaluate(1e-6,y0,genotype=WT)
    onset_row,onset_failures,_ = diagnose(stim,1e-6,y0,onset)
    assembly = inherited.bookkeeping(stim,onset)
    unchanged_hashes = {k:parameter_hashes(stim)[k] == parameter_hashes(old_stim)[k]
        for k in ('whole_cell','ae4','nbc','regulation','stimulus')}
    checks = {
        'exact_saved_core':bool(np.array_equal(y0[:12],saved['core_state'])),
        'rest_concentrations_exact':tuple(actual_ci[k] for k in ('na','k','cl'))==REST_CONCENTRATIONS,
        'task31_and_task38_rest_rhs_bitwise_equal':bool(np.array_equal(direct31.rhs,old.rhs)),
        'amount_rhs_within_1e_minus_12':float(np.max(np.abs(delta[inherited.AMOUNT_ROWS])))<=1e-12,
        'volume_rhs_within_1e_minus_14':float(np.max(np.abs(delta[inherited.VOLUME_ROWS])))<=1e-14,
        'regulatory_rhs_exact':bool(np.array_equal(new.rhs[12:],old.rhs[12:])),
        'voltages_within_1e_minus_12':max(map(abs,vols.values()))<=1e-12,
        'q_out_within_1e_minus_12':abs(water_delta)<=1e-12,
        'accepted_q_out_reproduced':abs(row['q_out_pL_s']-Q_REST)<=1e-12,
        'rest_cycle_reproduced':abs(row['nkcc1_cycle_inward_fmol_s']-REST_CYCLE)<=1e-15,
        'no_rest_physiology_or_conservation_failure':not failures,
        'onset_conservation_and_physiology':not onset_failures,
        'rest_N1_one_and_stimulated_N1_1p75':row['nkcc1_activity_multiplier']==1. and onset_row['nkcc1_activity_multiplier']==1.75,
        'N1_diagnostics_unchanged':all(new.diagnostics.regulatory[k]==old.diagnostics.regulatory[k]
            and onset.diagnostics.regulatory[k]==old_onset.diagnostics.regulatory[k]
            for k in old.diagnostics.regulatory if k.startswith('nkcc1')),
        'fixed_onset_only_nkcc_source_changes':bool(np.allclose(
            np.asarray(onset.rhs[:5])-np.asarray(old_onset.rhs[:5]),
            (onset.diagnostics.homeostasis.nkcc1_inward_fmol_s-old_onset.diagnostics.homeostasis.nkcc1_inward_fmol_s)*np.array([1,1,2,0,0]),rtol=0,atol=1e-12)),
        'whole_cell_bookkeeping':all(abs(v)<1e-12 for k,v in assembly.items() if k!='isolated_nbc_source'),
        'inherited_parameters_unchanged':all(unchanged_hashes.values()),
        'Palk_selected_all_layers':all(m.nkcc1_kinetics.law==PALK_NKCC1
            for m in (rest,rest.base_model,rest.base_model.base_model,stim,stim.base_model,stim.base_model.base_model)),
    }
    report = {'status':'PASS' if all(checks.values()) else 'PALK_NKCC1_REST_NESTING_FAILED',
        'checks':checks,'all_rhs_bitwise_equal':bool(np.array_equal(new.rhs,old.rhs)),
        'max_abs_amount_rhs_difference_fmol_s':float(np.max(np.abs(delta[inherited.AMOUNT_ROWS]))),
        'max_abs_volume_rhs_difference_pL_s':float(np.max(np.abs(delta[inherited.VOLUME_ROWS]))),
        'voltage_differences_V':vols,'q_out_difference_pL_s':water_delta,
        'old_rhs':list(map(float,old.rhs)),'palk_rhs':list(map(float,new.rhs)),
        'state_names':list(rest.state_names),'state_vector':list(map(float,y0)),
        'core_state_sha256':saved['core_state_sha256'],'rest_observables':row,
        'rest_conservation_ratios':ratios,'rest_gate_failures':failures,
        'onset_gate_failures':onset_failures,'onset_bookkeeping':assembly,
        'parameter_hashes':parameter_hashes(stim),'unchanged_inherited_hashes':unchanged_hashes,
        'stationary_solves':0,'integrations':0}
    write_json(OUT/'rest_nesting.json',report)
    budget.update(rest_nesting_pass=all(checks.values()),direct_rest_nesting_runs=1,
        status='PREFLIGHT_PASS' if all(checks.values()) else report['status'],
        no_further_scientific_execution=not all(checks.values()),
        numerical_execution_seconds=previous_seconds+time.perf_counter()-START)
    write_json(OUT/'budget.json',budget)
    print(json.dumps({'source':source,'rest_nesting':report,'budget':budget},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--literal-test-correction',action='store_true')
    main(literal_test_correction=parser.parse_args().literal_test_correction)
