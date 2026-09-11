"""Final WT log-square minimization with every state and voltage free."""
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import minimize
from .task21_calibration import *


def solve(root):
    require_stage0();manifest,_=load_freeze();p=WTProblem(manifest,root)
    profile=read(OUT/'profile_calibration'/f'{slug(root)}.json')
    if not profile['payload']:return None
    prior=profile['payload']['versions'][-1];start=np.array(prior['coordinates'])
    limit=profile['profile_score']+2e-7;lower,upper=p.bounds(limit)
    lower[6]=8.19  # Positive-carbon robustness branch, not a new physiology screen.
    start=np.clip(start,lower+1e-12,upper-1e-12)
    def eq(z):
        d=p.fast(z)
        return np.r_[d['equalities']/1000.,(d['cell']['cl']-50.1)/1.5,(d['cell']['ph']-6.91)/.07]
    def objective(z):return float(z[10:25]@z[10:25])
    def grad(z):
        result=np.zeros(27);result[10:25]=2*z[10:25];return result
    started=inverse.now()
    fit=minimize(objective,start,jac=grad,method='SLSQP',bounds=list(zip(lower,upper)),
        constraints=[{'type':'eq','fun':eq},{'type':'ineq','fun':p.constraints}],
        options={'maxiter':2200,'ftol':1e-12,'eps':1e-9})
    payload=p.production(fit.x)
    result=dict(root_id=root,started_utc=started,completed_utc=inverse.now(),
        success=bool(fit.success),message=fit.message,nit=fit.nit,nfev=fit.nfev,
        minimax_log_limit=limit,objective=payload['objective'],coordinates=fit.x.tolist(),
        candidate=payload,all_state_coordinates_free=True,genotype_results_used=False,
        max_abs_scaled_equality=float(max(abs(eq(fit.x)))))
    write_json(OUT/'joint_tie_break'/f'{slug(root)}.json',result)
    print(root,'joint tie',fit.message,payload['resting_gate_failures'],flush=True)
    return result


if __name__=='__main__':
    require_stage0();manifest,_=load_freeze()
    with ProcessPoolExecutor(max_workers=4) as pool:list(pool.map(solve,sorted(manifest['roots'])))
