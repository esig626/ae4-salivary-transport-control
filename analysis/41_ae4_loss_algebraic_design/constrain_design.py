"""Resolve the inverse problem using a pH admissibility constraint, no ODE."""
from inverse_balance import *


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--numerical-retry',action='store_true');parser.add_argument('--scaled-equation-solve',action='store_true');args=parser.parse_args()
    filename='constrained_design_scaled.json' if args.scaled_equation_solve else 'constrained_design_retry.json' if args.numerical_retry else 'constrained_design.json'
    start_time=time.perf_counter()
    attempted=read_json(OUT/'capacity_design.json')
    assert not (OUT/filename).exists()
    if args.numerical_retry:
        assert not read_json(OUT/'constrained_design.json')['solver_success']
    _,_,model,_,rest=old.frozen_models_and_state()
    span=COORDINATE_UPPER-COORDINATE_LOWER
    starts=[np.array(s['coordinates']) for s in attempted['states']]
    z0=np.r_[*[np.log((x-COORDINATE_LOWER)/(COORDINATE_UPPER-x)) for x in starts],np.log(attempted['capacity_fmol_s'])]
    # This margin is a physiology requirement, not another secretion target.
    target_null_ph=7.29
    calls=0

    def residual(z,details=False):
        nonlocal calls
        calls+=not details
        cap=float(np.exp(z[-1]));equations=[];records=[]
        for j,(name,g) in enumerate((('wt',WT),('ae4_null',AE4_NULL))):
            x=COORDINATE_LOWER+span*expit(z[j*10:(j+1)*10])
            y=state_from_coordinates(model,x);y[12:]=1.
            e=model.evaluate(600.,y,genotype=g)
            n=e.diagnostics.homeostasis.nkcc1_inward_fmol_s
            limited=float(np.clip(n,-cap,cap));delta=limited-n
            rhs=np.array(e.rhs);rhs[:3]+=np.array([delta,delta,2*delta])
            equations.extend(rhs[INDEPENDENT_ROWS]/INDEPENDENT_RHS_SCALES)
            row=old.parent.inherited.row_for(model,600.,y,e)
            row['nkcc1_cl_inward_fmol_s']=2*limited
            records.append({'case':name,'coordinates':x.tolist(),'state_vector':y.tolist(),
                'rhs':rhs.tolist(),'observables':row,'unrestricted_nkcc1_fmol_s':float(n),
                'limited_nkcc1_fmol_s':limited,'q_out_pL_s':row['q_out_pL_s']})
        equations.append(records[1]['observables']['ph_i']-target_null_ph)
        return (np.array(equations),records) if details else np.array(equations)

    options={'xtol':1e-11,'maxfev':1900}
    if args.numerical_retry:options.update(factor=.1,eps=1e-8)
    if args.scaled_equation_solve:
        from scipy.optimize import least_squares
        fit=least_squares(residual,z0,xtol=1e-12,ftol=1e-12,gtol=1e-12,x_scale='jac',max_nfev=2500)
        options={'method':'least_squares for square algebraic equations','x_scale':'jac','xtol':1e-12,'ftol':1e-12,'gtol':1e-12,'max_nfev':2500}
    else:
        fit=root(residual,z0,method='hybr',options=options)
    err,states=residual(fit.x,True);cap=float(np.exp(fit.x[-1]))
    failures=[]
    for state in states:
        x=state['observables']
        for key,ok in [('na_i_mM',0<x['na_i_mM']<40),('k_i_mM',50<=x['k_i_mM']<=200),
            ('cl_i_mM',30<=x['cl_i_mM']<=80),('ph_i',6.6<=x['ph_i']<=7.3),
            ('volume_i_pL',0<x['volume_i_pL']<3),('hco3_i_mM',0<x['hco3_i_mM']<100)]:
            if not ok:failures.append({'case':state['case'],'gate':key,'value':x[key]})
    deficit=1-states[1]['q_out_pL_s']/states[0]['q_out_pL_s']
    result={'solver_success':bool(fit.success),'solver_message':str(fit.message),
        'max_abs_scaled_design_residual':float(np.max(np.abs(err))),
        'actual_residual_calls':calls,'root_calls':1,'solver_options':options,'numerical_retry':args.numerical_retry,'capacity_fmol_s':cap,
        'design_constraint':'null steady pH=7.29, within original upper limit 7.3',
        'sustained_null_deficit':deficit,'sustained_deficit_in_requested_range':.20<=deficit<=.35,
        'physiology_failures':failures,'states':states,'production_trajectories_run':0,
        'changes_task40_rest':cap<rest['observables']['nkcc1_cycle_inward_fmol_s'],
        'wall_seconds':time.perf_counter()-start_time,
        'rejected_midpoint_design':{'capacity_fmol_s':attempted['capacity_fmol_s'],
            'null_ph':attempted['states'][1]['observables']['ph_i'],
            'null_cl_mM':attempted['states'][1]['observables']['cl_i_mM']}}
    write_json(OUT/filename,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
