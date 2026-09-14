"""Sequential null capacity and WT AE4-only gain algebraic solves."""
from design import *


def main():
    start=time.perf_counter();assert not (OUTPUT/'ae4_only_design.json').exists()
    _,_,template,_,_=old.frozen_models_and_state()
    prior=read_json(OUT/'capacity_design.json')
    lo=COORDINATE_LOWER[COORD_ROWS];hi=COORDINATE_UPPER[COORD_ROWS]
    results=[];cap=None;qnull=None
    for j,(case,genotype) in enumerate((('ae4_null',AE4_NULL),('wt',WT))):
        base=prior['states'][1-j]
        c=np.array(base['coordinates'])[COORD_ROWS]
        initial=np.r_[np.log((c-lo)/(hi-c)),np.log(prior['capacity_fmol_s'] if j==0 else 2.)]
        calls=0
        def residual(z,details=False):
            nonlocal calls
            calls+=not details
            C=float(np.exp(z[8])) if j==0 else cap
            G=1. if j==0 else float(np.exp(z[8]))
            y,x=water_eliminated_state(template,lo+(hi-lo)*expit(z[:8]))
            e=template.evaluate(600.,y,genotype=genotype)
            n=e.diagnostics.homeostasis.nkcc1_inward_fmol_s
            dn=float(np.clip(n,-C,C))-n;da=(G-1)*e.diagnostics.ae4.cl_cell_fmol_s
            rhs=np.array(e.rhs);rhs[:5]+=[dn-.5*da,dn-.5*da,2*dn+da,-2*da,-2*da]
            row,failures,ratios=old.diagnose(template,600.,y,e)
            eq=(row['cl_i_mM']-34.)/34. if j==0 else row['q_out_pL_s']/(qnull/.65)-1
            r=np.r_[rhs[ION_ROWS]/ION_SCALES,eq]
            record={'case':case,'capacity_fmol_s':C,'stimulated_ae4_gain':G,
                'coordinates':x.tolist(),'state_vector':y.tolist(),'rhs':rhs.tolist(),
                'observables_before_explicit_source_substitution':row,'failures':failures,
                'limited_nkcc1_cycles_fmol_s':float(np.clip(n,-C,C)),
                'amended_ae4_cl_fmol_s':G*e.diagnostics.ae4.cl_cell_fmol_s}
            return (r,record) if details else r
        options={'xtol':1e-10,'maxfev':1800,'factor':.1,'eps':1e-8}
        try:
            fit=root(residual,initial,method='hybr',options=options)
            error,record=residual(fit.x,True)
            record.update(solver_success=bool(fit.success),solver_message=str(fit.message),
                max_abs_scaled_residual=float(max(abs(error))))
        except Exception as exc:
            record={'case':case,'solver_success':False,'exception':repr(exc)}
        record.update(actual_residual_calls=calls,root_calls=1,solver_options=options)
        results.append(record)
        if not record['solver_success'] or record['max_abs_scaled_residual']>1e-7:break
        if j==0:
            cap=record['capacity_fmol_s'];qnull=record['observables_before_explicit_source_substitution']['q_out_pL_s']
    result={'states':results,'constraints':{'null_Cl_i_mM':34.,'null_sustained_deficit':.35},
        'root_calls':len(results),'actual_residual_calls':sum(r['actual_residual_calls'] for r in results),
        'production_trajectories_run':0,'wall_seconds':time.perf_counter()-start}
    write_json(OUTPUT/'ae4_only_design.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':main()
