"""Independent raw-state integrations and quadrature for inverse results.

Unlike inverse_analysis, this integrates only the 13 production states and
then integrates dense-output flow using scipy.integrate.quad.  It also checks
the inherited gates at accepted endpoints and the original observation grid.
"""
from core import *
from scipy.integrate import quad
from modern_full_model.task41_selected import Task41SelectedModel
from inverse_analysis import DETAIL,GRID,SELECTED_RHO,TARGET


def main():
    _,stim,y0,task,_=reference()
    crossings=json.loads((DETAIL/'crossings.json').read_text())
    selected=json.loads((DETAIL/'selected_reproduction.json').read_text())
    wt=selected['wt']['cumulative_adaptive_pL']
    checks=[]
    for name,rho,expression in [('selected_wt',SELECTED_RHO,1.),('selected_5pct',SELECTED_RHO,.05),
                                 ('selected_null',SELECTED_RHO,0.),('crossing_null',crossings['adaptive']['rho'],0.)]:
        m=Task41SelectedModel(stim) if name.startswith('selected') else cacc(stim,rho)
        g=replace(WT,ae4_expression=expression)
        sol=solve_ivp(lambda t,y:m.rhs(t,y,genotype=g),(0.,600.),y0,
            method='Radau',rtol=2e-9,atol=1e-11,max_step=1.,dense_output=True)
        assert sol.success
        q=lambda t:observe(m,t,sol.sol(t),g)['q_pL_s']
        volume,error=quad(q,0.,600.,epsabs=1e-10,epsrel=2e-9,
            points=[.01,.1,1.,5.,30.,60.,120.,300.],limit=200)
        sample=float(np.trapezoid([q(t) for t in GRID],GRID))
        failed=[];max_ratio=0.
        for t in np.unique(np.r_[sol.t,GRID]):
            y=sol.sol(t);ev=m.evaluate(t,y,genotype=g)
            _,ff,ratios=task.diagnose(m,t,y,ev)
            failed.extend(ff);max_ratio=max(max_ratio,max(ratios.values()))
        expected=(crossings['adaptive']['check'] if name=='crossing_null' else
            selected[{'selected_wt':'wt','selected_5pct':'ae4_5pct','selected_null':'ae4_null'}[name]])
        state_error=float(max(abs(sol.y[:,-1]-np.array(expected['final_state']))))
        rec=dict(case=name,rho=rho,expression=expression,raw_state_count=len(y0),
            adaptive_quad_pL=volume,quad_error_estimate_pL=error,
            augmented_integral_difference_pL=volume-expected['cumulative_adaptive_pL'],
            sampled_difference_pL=sample-expected['cumulative_sampled_pL'],
            endpoint_state_max_absolute_error=state_error,
            physiological=not failed,max_conservation_ratio=max_ratio)
        assert not failed
        assert abs(rec['augmented_integral_difference_pL'])<3e-8
        assert abs(rec['sampled_difference_pL'])<3e-8
        assert state_error<2e-5
        checks.append(rec)
    sensitivity=json.loads((DETAIL/'threshold_sensitivities.json').read_text())
    assert len(sensitivity)==16
    assert all(r['maximum_validation_relative_error']<2e-3 for r in sensitivity)
    assert all(r['difference_step_relative_error']<2e-3 for r in sensitivity)
    assert all(all(e['perturbed_gates_pass'] for e in r['estimates']) for r in sensitivity)
    long=json.loads((DETAIL/'long_time_summary.json').read_text())
    assert all(r['physiological_600s'] and not r['physiological'] for r in long.values())
    assert all(600<r['precise_ph_gate_crossing_s']<3600 for r in long.values())
    result=dict(status='PASS',independent_raw_state_integrations=checks,
        all_inherited_gates_checked=True,sensitivity_rows=len(sensitivity),
        threshold_tolerance_checks={k:v['solver_tolerance_rho_difference'] for k,v in crossings.items()},
        diagnostic_scope='Long extension checked separately; no retrospective rejection or alteration of Task 41',
        verification_method='Raw 13-state Radau integration plus independent adaptive dense-flow quadrature; original sampled trapezoid also reproduced')
    write_json(DETAIL/'independent_verification.json',result)
    print('INVERSE independent verification PASS:',len(checks),'raw-state trajectories;',len(sensitivity),'local sensitivity rows')


if __name__=='__main__':main()
