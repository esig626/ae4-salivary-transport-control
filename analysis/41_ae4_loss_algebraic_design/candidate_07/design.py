"""Derive residual CaCC recruitment from full stationary salt balances."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from reduced_inverse import *
from modern_full_model.model import ModernFullModel
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model

DIRECTORY=Path(__file__).resolve().parent
OUTPUT=OUT/'candidate_07'


def with_factor(template,factor):
    core=template.base_model.base_model
    parameters=replace(core.parameters,membranes=replace(core.parameters.membranes,
        g_cl_apical_S=core.parameters.membranes.g_cl_apical_S*factor))
    new=ModernFullModel(parameters,stimulus=core.stimulus,
        regulatory_model=core.regulatory_model,ae4_parameters=core.ae4_parameters,
        ae4_evaluator=core.ae4_evaluator,nkcc1_kinetics=core.nkcc1_kinetics)
    return MinimalNbcModel(StimulatedNkcc1Model(new),template.nbc_parameters)


def main():
    start=time.perf_counter();assert not (OUTPUT/'design.json').exists()
    _,_,template,_,_=old.frozen_models_and_state()
    prior=read_json(OUT/'candidate_06/design.json')
    lo=COORDINATE_LOWER[COORD_ROWS];hi=COORDINATE_UPPER[COORD_ROWS]
    results=[];qwt=None
    for j,(case,genotype) in enumerate((('wt',WT),('ae4_null',AE4_NULL))):
        base=prior['states'][j];c=np.array(base['coordinates'])[COORD_ROWS]
        initial=np.log((c-lo)/(hi-c))
        if j==1:
            b=base['cacc_conductance_factor'];initial=np.r_[initial,np.log(b/(1-b))]
        calls=0
        def residual(z,details=False):
            nonlocal calls
            calls+=not details
            factor=1. if j==0 else float(expit(z[8]));m=with_factor(template,factor)
            y,x=water_eliminated_state(m,lo+(hi-lo)*expit(z[:8]));e=m.evaluate(600.,y,genotype=genotype)
            row,failures,ratios=old.diagnose(m,600.,y,e)
            error=e.rhs[ION_ROWS]/ION_SCALES
            if j==1:error=np.r_[error,row['q_out_pL_s']/qwt-.7]
            rec={'case':case,'cacc_conductance_factor':factor,'coordinates':x.tolist(),
                'state_vector':y.tolist(),'rhs':e.rhs.tolist(),'observables':row,
                'failures':failures,'conservation_ratios':ratios}
            return (error,rec) if details else error
        options={'xtol':1e-10,'maxfev':1000,'factor':.1,'eps':1e-8}
        fit=root(residual,initial,method='hybr',options=options)
        error,record=residual(fit.x,True)
        record.update(solver_success=bool(fit.success),solver_message=str(fit.message),
            max_abs_scaled_residual=float(max(abs(error))),actual_residual_calls=calls,
            root_calls=1,solver_options=options)
        results.append(record)
        if not fit.success or max(abs(error))>1e-7:break
        qwt=record['observables']['q_out_pL_s'] if j==0 else qwt
    passed=len(results)==2 and all(r['solver_success'] and r['max_abs_scaled_residual']<1e-7 for r in results)
    result={'solver_success':passed,'states':results,'constraints':{'null_sustained_deficit':.30},
        'null_stimulated_cacc_fraction':results[-1]['cacc_conductance_factor'] if passed else None,
        'root_calls':len(results),'actual_residual_calls':sum(r['actual_residual_calls'] for r in results),
        'production_trajectories_run':0,'wall_seconds':time.perf_counter()-start}
    write_json(OUTPUT/'design.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
