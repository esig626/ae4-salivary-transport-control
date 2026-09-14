"""One square WT stationary root from the accepted Task 40 state; no fitting.

Use the inherited ten-coordinate charge manifold, residual scales, numerical
coordinate bounds, and acceptance tolerances. MINPACK hybr solves ten fixed
equations. Logit coordinates enforce the inherited bounds without changing
the scientific problem. No least-squares or parameter optimiser is called.
"""
import time
START=time.perf_counter()
from validation_common import *
from datetime import datetime, timezone
import signal
import traceback
from scipy.optimize import root
from scipy.special import expit
from threadpoolctl import threadpool_info
from modern_full_model.task30_nhe1_repair import (
    state_from_coordinates, COORDINATE_NAMES, COORDINATE_LOWER, COORDINATE_UPPER,
    INDEPENDENT_ROWS, INDEPENDENT_RHS_SCALES, OMITTED_CHARGE_ROWS,
    ROOT_SCALED_TOLERANCE, OMITTED_ROW_TOLERANCE, REGULATORY_TOLERANCE_S_INV,
    CHARGE_TOLERANCE_FMOL, CURRENT_TOLERANCE_A, BOUNDARY_RELATIVE_TOLERANCE,
    JACOBIAN_RANK_RELATIVE_TOLERANCE, AMOUNT_ROWS, VOLUME_ROWS, _flux_row)


def main():
    verify_parent_sources()
    budget=json.loads((OUT/'budget.json').read_text())
    assert budget['source_tests_pass'] and not budget['no_further_scientific_execution']
    assert budget['stationary_solves']==0 and not budget['started_cases']
    assert not (OUT/'wt_rest.json').exists()
    saved,old_rest,_,model,stim,y40=models_and_reference()
    before=parameter_hashes(model)
    start=np.asarray(saved['coordinates'],dtype=float)
    span=COORDINATE_UPPER-COORDINATE_LOWER
    assert np.all(start>COORDINATE_LOWER) and np.all(start<COORDINATE_UPPER)
    reconstructed=state_from_coordinates(model,start)
    assert np.max(np.abs(reconstructed-y40))<1e-12
    pools=threadpool_info();assert pools and all(p['num_threads']==1 for p in pools)
    budget['stationary_solves']=1;budget['status']='WT_REST_RUNNING'
    write_json(OUT/'budget.json',budget)
    z0=np.log((start-COORDINATE_LOWER)/(COORDINATE_UPPER-start))
    calls=0

    def alarm(*_): raise TimeoutError('Task 42 cumulative numerical execution limit')

    def coordinates(z): return COORDINATE_LOWER+span*expit(z)

    def residual(z):
        nonlocal calls
        calls+=1
        if calls>2400: raise RuntimeError('Inherited 2400 stationary-residual ceiling')
        # Use the literal accepted state at the sole initial guess.
        y=y40 if np.array_equal(z,z0) else state_from_coordinates(model,coordinates(z))
        return np.asarray(model.rhs(0.,y,genotype=WT))[INDEPENDENT_ROWS]/INDEPENDENT_RHS_SCALES

    signal.signal(signal.SIGALRM,alarm)
    signal.setitimer(signal.ITIMER_REAL,max(.001,900-budget['numerical_execution_seconds']-(time.perf_counter()-START)))
    result={'prepared_head':PREPARED_HEAD,'branch':BRANCH,'started_utc':datetime.now(timezone.utc).isoformat(),
        'sole_initial_guess_task40_core_sha256':saved['core_state_sha256'],
        'initial_state_vector':y40.tolist(),'initial_coordinates':start.tolist(),
        'coordinate_reconstruction_max_abs':float(np.max(np.abs(reconstructed-y40))),
        'solver':{'method':'scipy.optimize.root/hybr','unknowns':10,'equations':10,
            'coordinate_transform':'logit of inherited bounded charge-manifold coordinates',
            'xtol':1e-12,'max_actual_residual_evaluations':2400,'maxfev':2387},
        'retry_used':False,'optimisation_calls':0,'genotype':asdict(WT),
        'scientific_parameter_hashes':before,'status':'TASK42_SOURCE_CLASS_WT_REST_FAILED',
        'admissible':False,'failed_gates':[]}
    try:
        fit=root(residual,z0,method='hybr',options={'xtol':1e-12,'maxfev':2387})
        x=coordinates(fit.x);y=state_from_coordinates(model,x)
        e=model.evaluate(0.,y,genotype=WT)
        later=model.evaluate(600.,y,genotype=WT)
        row,failures,ratios=diagnose(model,0.,y,e)
        raw=np.asarray(e.rhs);scaled=raw[INDEPENDENT_ROWS]/INDEPENDENT_RHS_SCALES
        # Independent final Jacobian audit, not a second root solve.
        jac=np.empty((10,10))
        for k in range(10):
            dx=np.sqrt(np.finfo(float).eps)*max(abs(x[k]),1.)
            if x[k]+dx>=COORDINATE_UPPER[k]:dx=-dx
            xx=x.copy();xx[k]+=dx
            yy=state_from_coordinates(model,xx)
            rr=np.asarray(model.rhs(0.,yy,genotype=WT))[INDEPENDENT_ROWS]/INDEPENDENT_RHS_SCALES
            jac[:,k]=(rr-scaled)/dx
        singular=np.linalg.svd(jac*span[None,:],compute_uv=False)
        rank=int(np.sum(singular>JACOBIAN_RANK_RELATIVE_TOLERANCE*singular[0]))
        boundary=[key for key,d in zip(COORDINATE_NAMES,np.minimum(x-COORDINATE_LOWER,COORDINATE_UPPER-x)/span)
                  if d<=BOUNDARY_RELATIVE_TOLERANCE]
        checks={'root_solver_success':bool(fit.success),
            'independent_rhs':bool(np.max(np.abs(scaled))<=ROOT_SCALED_TOLERANCE),
            'omitted_charge_rows':bool(np.max(np.abs(raw[OMITTED_CHARGE_ROWS]))<=OMITTED_ROW_TOLERANCE),
            'regulatory_rhs':bool(np.max(np.abs(raw[12:]))<=REGULATORY_TOLERANCE_S_INV),
            'bulk_charge':bool(max(abs(v) for v in e.diagnostics.state_charge_fmol.values())<=CHARGE_TOLERANCE_FMOL),
            'current':bool(max(abs(v) for v in e.diagnostics.membranes.current_residuals_A.values())<=CURRENT_TOLERANCE_A),
            'independent_jacobian_full_rank':rank==10,'no_coordinate_boundary_hits':not boundary,
            'rest_protocol_time_invariance':bool(np.array_equal(raw,later.rhs)),
            'production_physiology_and_conservation':not failures}
        failed=failures+[{'gate':key} for key,value in checks.items() if not value]
        oldrow=parent.diagnose(old_rest,0.,y40,old_rest.evaluate(0.,y40))[0]
        shifts={key:row[key]-oldrow[key] for key in oldrow if key!='time_s'}
        result.update(status='PASS' if not failed else result['status'],admissible=not failed,
            solver_success=bool(fit.success),solver_message=str(fit.message),solver_nfev=int(fit.nfev),
            coordinates=x.tolist(),state_names=list(model.state_names),state_vector=y.tolist(),
            core_state=y[:12].tolist(),state_sha256=sha256_object(y.tolist()),
            core_state_sha256=sha256_object(y[:12].tolist()),core_shift_from_task40=(y[:12]-y40[:12]).tolist(),
            raw_rhs=raw.tolist(),scaled_independent_rhs=scaled.tolist(),
            max_abs_scaled_independent_rhs=float(np.max(np.abs(scaled))),
            max_abs_amount_rhs_fmol_s=float(np.max(np.abs(raw[AMOUNT_ROWS]))),
            max_abs_volume_rhs_pL_s=float(np.max(np.abs(raw[VOLUME_ROWS]))),
            normalized_jacobian_singular_values=singular.tolist(),independent_jacobian_rank=rank,
            boundary_hits=boundary,gates=checks,failed_gates=failed,observables=row,
            observed_shift_from_task40=shifts,flux_ledger=_flux_row(e,model),
            diagnostics=asdict(e.diagnostics),conservation_tolerances=dict(CONSERVATION_RESIDUAL_TOLERANCES),
            conservation_ratios=ratios,stationary_tolerances={
                'independent_rhs_scaled':ROOT_SCALED_TOLERANCE,'omitted_row_fmol_s':OMITTED_ROW_TOLERANCE,
                'regulatory_s_inv':REGULATORY_TOLERANCE_S_INV,'bulk_charge_fmol':CHARGE_TOLERANCE_FMOL,
                'current_A':CURRENT_TOLERANCE_A,'rank_relative':JACOBIAN_RANK_RELATIVE_TOLERANCE,
                'boundary_relative':BOUNDARY_RELATIVE_TOLERANCE},
            independent_bookkeeping=parent.inherited.bookkeeping(model,e),
            final_audit_rhs_evaluations=12)
    except Exception:
        result['exception']=traceback.format_exc()
        result['failed_gates'].append({'gate':'stationary_numerical_execution'})
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
    assert parameter_hashes(model)==before
    verify_parent_sources()
    result['actual_root_residual_evaluations']=calls
    result['wall_seconds']=time.perf_counter()-START
    write_json(OUT/'wt_rest.json',result)
    budget['stationary_residual_evaluations']=calls+result.get('final_audit_rhs_evaluations',0)
    budget['wt_rest_pass']=result['admissible']
    budget['numerical_execution_seconds']+=result['wall_seconds']
    budget['status']='WT_REST_FROZEN' if result['admissible'] else result['status']
    budget['no_further_scientific_execution']=not result['admissible']
    write_json(OUT/'budget.json',budget)
    if result['admissible']:
        files=[ROOT/'src/modern_full_model/ae4_catalan2025.py',
               ROOT/'tests/test_task42_catalan2025.py',*HERE.glob('*.py')]
        write_json(OUT/'frozen_inputs.json',{'prepared_head':PREPARED_HEAD,'branch':BRANCH,
            'scientific_parent':'fe51473d7cc201e2e04007aa62fbc1b45098732e',
            'frozen_utc':datetime.now(timezone.utc).isoformat(),
            'frozen_rest_state_sha256':result['state_sha256'],'ae4_routing_model':MODEL_ID,
            'rest_parameter_hashes':before,'stimulus_parameter_hashes':parameter_hashes(stim),
            'task42_execution_files':[{'path':str(p.relative_to(ROOT)),
                'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(set(files))]})
    print(json.dumps({k:result[k] for k in ('status','admissible','actual_root_residual_evaluations','failed_gates','wall_seconds')},indent=2))
    # Detailed observables are retained in the class rest record.


if __name__=='__main__':main()
