"""Algebraic pH-dependent CaCC feedback; original Task 40 otherwise."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from reduced_inverse import *
from modern_full_model.model import ModernFullModel
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model

DIRECTORY=Path(__file__).resolve().parent
OUTPUT=OUT/'candidate_06'
WIDTH=.02


def gated_at_ph(template,ph,ph_half):
    core=template.base_model.base_model
    factor=float(expit((ph_half-ph)/WIDTH))
    parameters=replace(core.parameters,membranes=replace(core.parameters.membranes,
        g_cl_apical_S=core.parameters.membranes.g_cl_apical_S*factor))
    new=ModernFullModel(parameters,stimulus=core.stimulus,
        regulatory_model=core.regulatory_model,ae4_parameters=core.ae4_parameters,
        ae4_evaluator=core.ae4_evaluator,nkcc1_kinetics=core.nkcc1_kinetics)
    return MinimalNbcModel(StimulatedNkcc1Model(new),template.nbc_parameters),factor


def main():
    start=time.perf_counter();assert not (OUTPUT/'design.json').exists()
    _,_,template,_,_=old.frozen_models_and_state()
    lo=COORDINATE_LOWER[COORD_ROWS];hi=COORDINATE_UPPER[COORD_ROWS]
    starts=[]
    for case,genotype in (('wt',WT),('ae4_null',AE4_NULL)):
        y=read_json(ROOT/'results/40_ae4_equal_cation_routing'/f'{case}_verification.json')['last_checked_state_vector']
        starts.append(state_coordinates(template,y,genotype)[COORD_ROWS])
    initial=np.r_[*[np.log((c-lo)/(hi-c)) for c in starts],7.15]
    calls=0
    def residual(z,details=False):
        nonlocal calls
        calls+=not details
        rows=[];states=[]
        for j,(case,genotype) in enumerate((('wt',WT),('ae4_null',AE4_NULL))):
            c=lo+(hi-lo)*expit(z[j*8:(j+1)*8]);m,factor=gated_at_ph(template,c[3],z[16])
            y,x=water_eliminated_state(m,c);e=m.evaluate(600.,y,genotype=genotype)
            row,failures,ratios=old.diagnose(m,600.,y,e)
            rows.extend(e.rhs[ION_ROWS]/ION_SCALES)
            states.append({'case':case,'coordinates':x.tolist(),'state_vector':y.tolist(),
                'rhs':e.rhs.tolist(),'observables':row,'failures':failures,
                'cacc_conductance_factor':factor,'conservation_ratios':ratios})
        rows.append(states[1]['observables']['q_out_pL_s']/states[0]['observables']['q_out_pL_s']-.70)
        return (np.array(rows),states) if details else np.array(rows)
    options={'xtol':1e-10,'maxfev':2500,'factor':.1,'eps':1e-8}
    try:
        fit=root(residual,initial,method='hybr',options=options)
        error,states=residual(fit.x,True)
        result={'solver_success':bool(fit.success),'solver_message':str(fit.message),
            'max_abs_scaled_residual':float(max(abs(error))),'ph_half':float(fit.x[16]),
            'width_ph':WIDTH,'constraints':{'null_sustained_flow_deficit':.30},'states':states}
    except Exception as exc:
        result={'solver_success':False,'exception':repr(exc)}
    result.update(root_calls=1,actual_residual_calls=calls,solver_options=options,
        production_trajectories_run=0,wall_seconds=time.perf_counter()-start)
    write_json(OUTPUT/'design.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
