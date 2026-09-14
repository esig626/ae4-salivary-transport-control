"""Algebraic co-recruitment and finite NKCC1 capacity selection."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from reduced_inverse import *
from modern_full_model.nbc_minimal import MinimalNbcModel

DIRECTORY=Path(__file__).resolve().parent
OUTPUT=OUT/'candidate_03'


def main():
    start=time.perf_counter();assert not (OUTPUT/'design.json').exists()
    _,_,template,_,_=old.frozen_models_and_state()
    prior=read_json(OUT/'reduced_capacity_design.json')
    lo=COORDINATE_LOWER[COORD_ROWS];hi=COORDINATE_UPPER[COORD_ROWS]
    initial=np.r_[*[np.log((np.array(s['coordinates'])[COORD_ROWS]-lo)/(hi-np.array(s['coordinates'])[COORD_ROWS])) for s in prior['states']],np.log(prior['capacity_fmol_s']),np.log(2.)]
    calls=0
    def residual(z,details=False):
        nonlocal calls
        calls+=not details
        cap,gain=np.exp(z[16:]);m=MinimalNbcModel(template.base_model,replace(template.nbc_parameters,capacity_fmol_s=template.nbc_parameters.capacity_fmol_s*gain))
        rows=[];states=[]
        for j,(case,genotype) in enumerate((('wt',WT),('ae4_null',AE4_NULL))):
            c=lo+(hi-lo)*expit(z[j*8:(j+1)*8]);y,x=water_eliminated_state(m,c)
            e=m.evaluate(600.,y,genotype=genotype)
            n=e.diagnostics.homeostasis.nkcc1_inward_fmol_s
            dn=float(np.clip(n,-cap,cap))-n
            da=(gain-1)*e.diagnostics.ae4.cl_cell_fmol_s
            rhs=np.array(e.rhs);rhs[:5]+=[dn-.5*da,dn-.5*da,2*dn+da,-2*da,-2*da]
            row,failures,ratios=old.diagnose(m,600.,y,e)
            row['nkcc1_cycle_inward_fmol_s']=float(np.clip(n,-cap,cap))
            row['nkcc1_cl_inward_fmol_s']=2*row['nkcc1_cycle_inward_fmol_s']
            for key in ('ae4_cl_inward_fmol_s','ae4_na_cell_source_fmol_s','ae4_k_cell_source_fmol_s','ae4_tic_cell_source_fmol_s','ae4_ta_cell_source_fmol_s'):row[key]*=gain
            rows.extend(rhs[ION_ROWS]/ION_SCALES)
            states.append({'case':case,'coordinates':x.tolist(),'state_vector':y.tolist(),
                'rhs':rhs.tolist(),'observables':row,'failures':failures,'conservation_ratios':ratios})
        rows.extend([(states[1]['observables']['cl_i_mM']-40.)/40.,
            states[1]['observables']['q_out_pL_s']/states[0]['observables']['q_out_pL_s']-.67])
        return (np.array(rows),states) if details else np.array(rows)
    options={'xtol':1e-10,'maxfev':2500,'factor':.1,'eps':1e-8}
    try:
        fit=root(residual,initial,method='hybr',options=options)
        error,states=residual(fit.x,True)
        result={'solver_success':bool(fit.success),'solver_message':str(fit.message),
            'max_abs_scaled_residual':float(max(abs(error))),
            'capacity_fmol_s':float(np.exp(fit.x[16])),'stimulated_ae4_nbc_gain':float(np.exp(fit.x[17])),
            'constraints':{'null_Cl_i_mM':40.,'null_sustained_flow_deficit':.33},'states':states}
    except Exception as exc:
        result={'solver_success':False,'exception':repr(exc)}
    result.update(root_calls=1,actual_residual_calls=calls,solver_options=options,
        production_trajectories_run=0,wall_seconds=time.perf_counter()-start)
    write_json(OUTPUT/'design.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
