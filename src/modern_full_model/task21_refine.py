"""WT-only equality-constrained follow-up to Task 21 feasibility restoration."""
from concurrent.futures import ProcessPoolExecutor
import math
import time
import numpy as np
from scipy.optimize import minimize
from .task21_calibration import (WTProblem, require_stage0, read, OUT, load_freeze,
                                 slug, write_json, inverse)
from .task21_flux_search import coordinates, inner


def refine_seed(root):
    require_stage0();manifest,_=load_freeze();p=WTProblem(manifest,root)
    prior=read(OUT/'seed_attempts'/f'{slug(root)}.json')
    path=OUT/'refinements'/f'{slug(root)}.json'
    results=read(path)['attempts'] if path.exists() else []
    flux=read(OUT/'flux_search'/f'{slug(root)}.json')
    fz,_=coordinates(p,flux['outer_coordinates']);_,fz=inner(p,fz)
    fz[10:25]=np.clip(fz[10:25],-16,39)
    prior['attempts'].append(dict(kind='analytic_flux_start',coordinates=fz.tolist()))
    for attempt in prior['attempts']:
        if any(r['start_kind']==attempt['kind'] for r in results):continue
        z=np.asarray(attempt['coordinates']);lower,upper=p.bounds(40.)
        x=np.r_[z,max(abs(z[10:25]))+1.]
        def eq(x):
            try:
                d=p.fast(x[:27])
                scale=1000. if attempt['kind']=='analytic_flux_start' else 1.
                return np.r_[d['equalities']/scale,(d['cell']['cl']-50.1)/1.5,(d['cell']['ph']-6.91)/.07]
            except (ValueError,OverflowError):return np.full(14,1e15)
        def ine(x):
            try:return np.r_[p.constraints(x[:27]),x[27]-x[10:25],x[27]+x[10:25]]
            except (ValueError,OverflowError):return np.full(35,-1e15)
        start=time.monotonic()
        fit=minimize(lambda x:x[27],x,method='SLSQP',bounds=list(zip(np.r_[lower,0.],np.r_[upper,40.])),
            constraints=[{'type':'eq','fun':eq},{'type':'ineq','fun':ine}],
            options={'maxiter':3000,'ftol':1e-12})
        try:payload=p.production(fit.x[:27])
        except (ValueError,OverflowError):
            payload=p.production(z)
            fit.message=str(fit.message)+'; final state invalid; initial candidate preserved'
        result=dict(start_kind=attempt['kind'],completed_utc=inverse.now(),elapsed_s=time.monotonic()-start,
            success=bool(fit.success),message=fit.message,nit=fit.nit,nfev=fit.nfev,
            max_abs_equality=float(max(abs(eq(fit.x)))),minimum_inequality=float(min(ine(fit.x))),
            coordinates=fit.x[:27].tolist(),candidate=payload)
        results.append(result)
        print(root,'minimax',attempt['kind'],fit.message,payload['objective'],
              payload['resting_gate_failures'],flush=True)
        write_json(path,dict(root_id=root,attempts=results,
            genotype_results_used=False,measured_means_constrained=True))
    return results


if __name__=='__main__':
    require_stage0();manifest,_=load_freeze()
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(refine_seed,sorted(manifest['roots'])))
