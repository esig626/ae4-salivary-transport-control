"""Eliminate stationary water exactly, then solve eight ionic equations/case."""
from inverse_balance import *
from modern_full_model.acid_base import total_alkalinity_mM

COORD_ROWS=np.array([0,1,2,3,5,6,7,8])
ION_ROWS=np.array([0,1,2,3,6,7,8,9])
ION_SCALES=np.array([.05,.05,.05,.02,.02,.02,.02,.02])


def water_eliminated_state(model,coordinates8,activation=1.):
    p=model.parameters;w=p.water;g=p.geometry
    ni,ki,ti,hi,nl,kl,tl,hl=coordinates8
    ta_i=total_alkalinity_mM(ph=hi,total_carbon_mM=ti,buffer_total_mM=0.,buffer_pka=p.acid_base.cell_buffer_pka,parameters=p.acid_base)
    ta_l=total_alkalinity_mM(ph=hl,total_carbon_mM=tl,buffer_total_mM=0.,buffer_pka=p.acid_base.lumen_buffer_pka,parameters=p.acid_base)
    bi=g.cell_buffer_total_fmol/(1+10**(p.acid_base.cell_buffer_pka-hi))
    bl=g.lumen_buffer_total_fmol/(1+10**(p.acid_base.lumen_buffer_pka-hl))
    ai=2*(ni+ki)+ti-ta_i
    beta_i=g.cell_impermeant_osmoles_fmol+g.cell_buffer_total_fmol-bi-g.fixed_cell_anion_equivalents_fmol
    al=2*(nl+kl)+tl-ta_l
    # O_l=al-bl/V_l. The positive conducting lumen uses the larger root.
    bath=p.bath
    ob=bath.na_mM+bath.k_mM+bath.cl_mM+bath.tic_mM+bath.untracked_osmolyte_mM
    hydraulic=w.apical_hydraulic_pL_s_mOsm*w.basolateral_hydraulic_pL_s_mOsm/(w.apical_hydraulic_pL_s_mOsm+w.basolateral_hydraulic_pL_s_mOsm)+w.paracellular_hydraulic_pL_s_mOsm
    k=hydraulic/w.outflow_rate_s
    a=w.lumen_dead_volume_pL+k*(al-ob)
    vl=(a+np.sqrt(a*a-4*k*bl))/2
    ol=al-bl/vl
    oi=(w.basolateral_hydraulic_pL_s_mOsm*ob+w.apical_hydraulic_pL_s_mOsm*ol)/(w.basolateral_hydraulic_pL_s_mOsm+w.apical_hydraulic_pL_s_mOsm)
    vi=beta_i/(oi-ai)
    x=np.array([ni,ki,ti,hi,vi,nl,kl,tl,hl,vl])
    y=state_from_coordinates(model,x);y[12:]=activation
    return y,x


def solve_design(*,allow_ae4_recruitment=False, target_deficit=.25, null_ph=7.299, scaled=False):
    clock=time.perf_counter()
    filename='two_parameter_design.json' if allow_ae4_recruitment else 'reduced_capacity_design.json'
    if scaled:filename=filename.replace('.json','_scaled.json')
    assert not (OUT/filename).exists()
    _,_,model,_,rest=old.frozen_models_and_state()
    attempt=read_json(OUT/'constrained_design_scaled.json')
    lower=COORDINATE_LOWER[COORD_ROWS];upper=COORDINATE_UPPER[COORD_ROWS]
    initial=np.r_[*[np.log((np.array(s['coordinates'])[COORD_ROWS]-lower)/(upper-np.array(s['coordinates'])[COORD_ROWS])) for s in attempt['states']],np.log(attempt['capacity_fmol_s'])]
    if allow_ae4_recruitment:initial=np.r_[initial,np.log(2.)]
    calls=0

    def residual(z,details=False):
        nonlocal calls
        calls+=not details
        cap=float(np.exp(z[16]));gain=float(np.exp(z[17])) if allow_ae4_recruitment else 1.
        rows=[];records=[]
        for j,(name,genotype) in enumerate((('wt',WT),('ae4_null',AE4_NULL))):
            c=lower+(upper-lower)*expit(z[j*8:(j+1)*8])
            y,x=water_eliminated_state(model,c)
            e=model.evaluate(600.,y,genotype=genotype)
            n=e.diagnostics.homeostasis.nkcc1_inward_fmol_s
            limited=float(np.clip(n,-cap,cap));dn=limited-n
            da=(gain-1)*e.diagnostics.ae4.cl_cell_fmol_s
            rhs=np.array(e.rhs);rhs[:5]+=[dn-.5*da,dn-.5*da,2*dn+da,-2*da,-2*da]
            rows.extend(rhs[ION_ROWS]/ION_SCALES)
            row=old.parent.inherited.row_for(model,600.,y,e)
            row['nkcc1_cl_inward_fmol_s']=2*limited
            row['ae4_cl_inward_fmol_s']*=gain
            records.append({'case':name,'coordinates':x.tolist(),'state_vector':y.tolist(),
                'rhs':rhs.tolist(),'observables':row,'unrestricted_nkcc1_fmol_s':float(n),
                'limited_nkcc1_fmol_s':limited,'q_out_pL_s':row['q_out_pL_s']})
        rows.append(records[1]['observables']['ph_i']-(null_ph if allow_ae4_recruitment else 7.29))
        if allow_ae4_recruitment:rows.append(records[1]['q_out_pL_s']/records[0]['q_out_pL_s']-(1-target_deficit))
        return (np.array(rows),records) if details else np.array(rows)

    if scaled:
        from scipy.optimize import least_squares
        fit=least_squares(residual,initial,x_scale='jac',xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=1500)
    else:
        fit=root(residual,initial,method='hybr',options={'xtol':1e-11,'maxfev':1900})
    error,states=residual(fit.x,True)
    result={'solver_success':bool(fit.success),'solver_message':str(fit.message),
        'root_calls':1,'least_squares_equation_solver':scaled,'actual_residual_calls':calls,'equations':len(initial),
        'stationary_water_exactly_eliminated':True,'max_abs_scaled_design_residual':float(np.max(np.abs(error))),
        'max_abs_water_rhs':max(abs(s['rhs'][k]) for s in states for k in (5,11)),
        'capacity_fmol_s':float(np.exp(fit.x[16])),
        'stimulated_ae4_extra_gain':float(np.exp(fit.x[17])) if allow_ae4_recruitment else 1.,
        'sustained_null_deficit':1-states[1]['q_out_pL_s']/states[0]['q_out_pL_s'],
        'constraints':{'null_ph':null_ph if allow_ae4_recruitment else 7.29,'null_sustained_deficit':target_deficit if allow_ae4_recruitment else None},
        'states':states,'production_trajectories_run':0,'wall_seconds':time.perf_counter()-clock}
    write_json(OUT/filename,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--allow-ae4-recruitment',action='store_true');parser.add_argument('--scaled',action='store_true')
    args=parser.parse_args();solve_design(allow_ae4_recruitment=args.allow_ae4_recruitment,scaled=args.scaled)
