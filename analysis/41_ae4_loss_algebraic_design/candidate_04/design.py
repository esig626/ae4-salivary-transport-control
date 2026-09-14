"""Shared NKCC1 ceiling and pump Na affinity inverse stationary design."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from reduced_inverse import *
from modern_full_model.model import ModernFullModel
from modern_full_model.nbc_minimal import MinimalNbcModel
from modern_full_model.nkcc_stimulation import StimulatedNkcc1Model
from modern_full_model.nkcc1_capacity_limited import CapacityLimitedNkcc1Model,Nkcc1Capacity

DIRECTORY=Path(__file__).resolve().parent
OUTPUT=OUT/'candidate_04'


def amended(template,capacity,pump_na_half):
    core=template.base_model.base_model
    parameters=replace(core.parameters,membranes=replace(core.parameters.membranes,
        nak_na_half_mM=pump_na_half))
    new=ModernFullModel(parameters,stimulus=core.stimulus,
        regulatory_model=core.regulatory_model,ae4_parameters=core.ae4_parameters,
        ae4_evaluator=core.ae4_evaluator,nkcc1_kinetics=core.nkcc1_kinetics)
    return CapacityLimitedNkcc1Model(MinimalNbcModel(StimulatedNkcc1Model(new),template.nbc_parameters),Nkcc1Capacity(capacity))


def main():
    start=time.perf_counter();assert not (OUTPUT/'design.json').exists()
    _,_,template,_,_=old.frozen_models_and_state()
    prior=read_json(OUT/'capacity_design.json')
    lo=COORDINATE_LOWER[COORD_ROWS];hi=COORDINATE_UPPER[COORD_ROWS]
    initial=np.r_[*[np.log((np.array(s['coordinates'])[COORD_ROWS]-lo)/(hi-np.array(s['coordinates'])[COORD_ROWS])) for s in prior['states']],np.log(prior['capacity_fmol_s']),np.log(10.)]
    calls=0
    def residual(z,details=False):
        nonlocal calls
        calls+=not details
        cap,half=np.exp(z[16:]);m=amended(template,cap,half)
        rows=[];states=[]
        for j,(case,genotype) in enumerate((('wt',WT),('ae4_null',AE4_NULL))):
            c=lo+(hi-lo)*expit(z[j*8:(j+1)*8]);y,x=water_eliminated_state(m,c)
            e=m.evaluate(600.,y,genotype=genotype)
            row,failures,ratios=old.diagnose(m,600.,y,e)
            rows.extend(e.rhs[ION_ROWS]/ION_SCALES)
            states.append({'case':case,'coordinates':x.tolist(),'state_vector':y.tolist(),
                'rhs':e.rhs.tolist(),'observables':row,'failures':failures,'conservation_ratios':ratios})
        rows.extend([states[1]['observables']['ph_i']-7.28,
            states[1]['observables']['q_out_pL_s']/states[0]['observables']['q_out_pL_s']-.67])
        return (np.array(rows),states) if details else np.array(rows)
    options={'xtol':1e-10,'maxfev':2500,'factor':.1,'eps':1e-8}
    try:
        fit=root(residual,initial,method='hybr',options=options)
        error,states=residual(fit.x,True)
        result={'solver_success':bool(fit.success),'solver_message':str(fit.message),
            'max_abs_scaled_residual':float(max(abs(error))),
            'capacity_fmol_s':float(np.exp(fit.x[16])),'pump_na_half_mM':float(np.exp(fit.x[17])),
            'constraints':{'null_pH':7.28,'null_sustained_flow_deficit':.33},'states':states}
    except Exception as exc:
        result={'solver_success':False,'exception':repr(exc)}
    result.update(root_calls=1,actual_residual_calls=calls,solver_options=options,
        production_trajectories_run=0,wall_seconds=time.perf_counter()-start)
    write_json(OUTPUT/'design.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
