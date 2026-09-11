"""WT measured means and joint state/capacity minimax profile using linear closure."""
from concurrent.futures import ProcessPoolExecutor
import time
import numpy as np
from scipy.optimize import differential_evolution, linprog, minimize
from .task21_flux_search import coordinates, linear_system, inner
from .task21_calibration import *


def profile(p,x,return_payload=False):
    z,pen=coordinates(p,x[:10])
    if z is None:return (100+pen,None) if return_payload else 100+pen
    A,b=linear_system(p,z);rs=np.maximum(abs(b),1);A/=rs[:,None];b/=rs
    norms=np.linalg.norm(A,axis=0);As=A/norms
    lower=np.exp(-x[10])*norms;upper=np.exp(x[10])*norms
    result=linprog(np.zeros(15),A_eq=As,b_eq=b,bounds=list(zip(lower,upper)),method='highs',
        options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not result.success:
        loss=40+inner(p,z)[0]
        return (loss,None) if return_payload else loss
    if not return_payload:return x[10]
    logs=np.log(np.maximum(result.x/norms,1e-100))
    def eq(v):return A@np.exp(v)-b
    tie=minimize(lambda v:float(v@v),logs,jac=lambda v:2*v,method='SLSQP',
        bounds=[(-x[10]-1e-7,x[10]+1e-7)]*15,
        constraints=[{'type':'eq','fun':eq,'jac':lambda v:A*np.exp(v)[None,:]}],
        options={'ftol':1e-12,'maxiter':2000})
    versions=[]
    for label,v in (('minimax_profile',logs),('capacity_tie_break',tie.x)):
        zz=z.copy();zz[10:25]=v
        versions.append(dict(stage=label,coordinates=zz.tolist(),candidate=p.production(zz)))
    return x[10],dict(versions=versions,tie_break=dict(success=bool(tie.success),message=tie.message,
        nfev=tie.nfev,nit=tie.nit,max_abs_scaled_linear_residual=float(max(abs(eq(tie.x)))),
        tie_break_scope='Capacities at the joint minimax state/voltages; no global tie-break certificate'))


def search(root):
    require_stage0();manifest,_=load_freeze();p=WTProblem(manifest,root)
    old=read(OUT/'mean_flux_search'/f'{slug(root)}.json')
    bounds=old['bounds']+[(0,32)]
    seed=int(hashlib.sha256((root+'profile').encode()).hexdigest()[:8],16)
    rng=np.random.default_rng(seed);lo,hi=np.array(bounds).T
    population=rng.uniform(lo,hi,size=(100,11))
    population[:20,:10]=old['outer_coordinates'];population[:20,10]=np.linspace(20,32,20)
    started=inverse.now();begin=time.monotonic()
    fit=differential_evolution(lambda x:profile(p,x),bounds,init=population,seed=seed,
        maxiter=250,tol=1e-7,polish=False)
    score,payload=profile(p,fit.x,True)
    result=dict(root_id=root,started_utc=started,completed_utc=inverse.now(),elapsed_s=time.monotonic()-begin,
        seed=seed,nfev=fit.nfev,success=bool(fit.success),message=fit.message,profile_score=float(score),
        outer_coordinates=fit.x.tolist(),bounds=bounds,payload=payload,genotype_results_used=False,
        measured_objective=0.,global_optimality_certified=False)
    write_json(OUT/'profile_calibration'/f'{slug(root)}.json',result)
    print(root,'profile',score,[(r['stage'],r['candidate']['resting_gate_failures']) for r in payload['versions']] if payload else 'no feasible profile',flush=True)
    return result


if __name__=='__main__':
    require_stage0();manifest,_=load_freeze()
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(search,sorted(manifest['roots'])))
