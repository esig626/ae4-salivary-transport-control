"""Lexicographic WT calibration in exact material/water balance coordinates.

Five independent pathway fluxes parameterize the ten material equations.
Capacity signs and folds, including the two voltage-current identities, are
then enforced without an ill-conditioned equality-constrained inverse solve.
"""
from concurrent.futures import ProcessPoolExecutor
import time
import numpy as np
from scipy.optimize import linprog, minimize
from .task21_flux_search import coordinates, linear_system
from .task21_calibration import *


def exact(p,x):
    z,penalty=coordinates(p,x[:10])
    if z is None:return None,penalty
    d=p.fast(z);ci=d['cell'];li=d['lumen'];q=d['flow'];part=ci['na']/(ci['na']+ci['k'])
    n,e,j=x[10:13]*1000*q;fraction=x[13];ka=x[14]*1000*q
    h=e+2*j;pt=(n+h-part*j)/3;pa=fraction*pt;pb=pt-pa
    kb=n-(1-part)*j+2*pt-ka;cac=2*n+e+j
    coa=q*(li['alkalinity']-li['tic']);cob=h-coa
    flux=np.array([j,n,e,h,pa,pb,ka,kb,cac,cob,coa,3*pa-q*li['na'],
                   ka-2*pa-q*li['k'],cac-q*li['cl'],-q*li['alkalinity']])
    ratios=flux/d['pathway_fluxes']
    return (z,ratios),0.


def search(root):
    require_stage0();manifest,_=load_freeze();p=WTProblem(manifest,root)
    old=read(OUT/'flux_search'/f'{slug(root)}.json')
    state=np.array(old['outer_coordinates']);z,_=coordinates(p,state)
    A,b=linear_system(p,z);rs=np.maximum(abs(b),1);A/=rs[:,None];b/=rs
    norm=np.linalg.norm(A,axis=0);A/=norm
    fit=linprog(np.ones(15),A_eq=A,b_eq=b,bounds=list(zip(np.exp(-6)*norm,np.exp(30)*norm)),
        method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not fit.success:raise RuntimeError(fit.message)
    z[10:25]=np.log(fit.x/norm);d=p.fast(z);q=d['flow'];j,n,e,h,pa,pb,ka=d['pathway_fluxes'][:7]
    x=np.r_[state,[n/q/1000,e/q/1000,j/q/1000,pa/(pa+pb),ka/q/1000],31.]
    bounds=old['bounds']+[(-100,100),(-100,100),(1e-14,100),(1e-12,1-1e-12),(-100,100),(0,40)]
    def inequalities(x):
        out,pen=exact(p,x)
        if out is None:return np.full(30,-100-pen)
        _,ratios=out
        signed=np.sign(ratios)*np.log1p(abs(ratios))
        return np.r_[signed-np.log1p(np.exp(-x[15])),np.log1p(np.exp(x[15]))-signed]
    def measured(x):return ((x[2]-50.1)/1.5)**2+((x[3]-6.91)/.07)**2
    constraints=[{'type':'ineq','fun':inequalities}]
    records=[]
    def solve(label,objective,constraints,x):
        begin=time.monotonic()
        fit=minimize(objective,x,method='SLSQP',bounds=bounds,constraints=constraints,
                     options={'maxiter':3500,'ftol':1e-12})
        out,pen=exact(p,fit.x);payload=None;coord=None
        if out is not None:
            coord,ratios=out
            if np.all(ratios>0):
                coord[10:25]=np.log(ratios);payload=p.production(coord)
        record=dict(stage=label,success=bool(fit.success),message=fit.message,nit=fit.nit,nfev=fit.nfev,
            elapsed_s=time.monotonic()-begin,objective=float(objective(fit.x)),
            measured_objective=float(measured(fit.x)),minimum_inequality=float(min(inequalities(fit.x))),
            exact_coordinates=fit.x.tolist(),coordinates=coord.tolist() if coord is not None else None,
            candidate=payload)
        records.append(record)
        print(root,label,record['measured_objective'],fit.x[15],fit.message,
              payload['resting_gate_failures'] if payload else 'nonpositive capacity',flush=True)
        write_json(OUT/'exact_calibration'/f'{slug(root)}.json',dict(root_id=root,stages=records,
            genotype_results_used=False,stage0_checkpoint_sha=STAGE0_SHA,completed_utc=inverse.now()))
        return fit.x
    x=solve('1_measured_fit',measured,constraints,x)
    target=measured(x)
    constraints=constraints+[{'type':'ineq','fun':lambda x:target+1e-10-measured(x)}]
    x=solve('2_minimax_capacity',lambda x:x[15],constraints,x)
    maximum=x[15]
    constraints=constraints+[{'type':'ineq','fun':lambda x:maximum+1e-7-x[15]}]
    def sumlogs(x):
        out,pen=exact(p,x)
        if out is None:return 1e8+pen
        _,r=out
        return float(np.sum(np.log(np.maximum(r,1e-30))**2))
    solve('3_sum_squared_logs',sumlogs,constraints,x)
    return records


if __name__=='__main__':
    require_stage0();manifest,_=load_freeze()
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(search,sorted(manifest['roots'])))
