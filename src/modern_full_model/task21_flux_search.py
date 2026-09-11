"""Independent WT robustness search with exact water elimination and NNLS.

Capacity scales do not enter the outer search. Inner nonnegative flux fitting
is normalized by outflow, preventing a vanishing-rate near-root from looking
like a feasible steady state. This is a numerical search, not a certificate.
"""
from concurrent.futures import ProcessPoolExecutor
import math
import numpy as np
from scipy.optimize import differential_evolution, nnls
from scipy.special import logit, expit
from .acid_base import _carbon_fractions
from .task21_calibration import *


def coordinates(p,x):
    na,k,cl,ph,u,phl,v,odds,vtmargin,va=x
    b=p.p.bath; g=p.p.geometry; a=p.p.acid_base;w=p.p.water
    ci0,ci1,ci2=_carbon_fractions(ph,a)
    beta=g.cell_buffer_total_fmol/(1+10**(a.cell_buffer_pka-ph))
    fixed=g.fixed_cell_anion_equivalents_fmol+beta
    free=1e3*(10**(ph-a.water_pkw)-10**(-ph))
    ticmin=p.bath_ab.hco3_mM*math.sqrt(cl*(b.na_mM+b.k_mM)/(b.cl_mM*(na+k)))/ci1*(1+1e-6)
    ticmax=(na+k-cl-fixed/5-free)/(ci1+2*ci2)
    if ticmax<=ticmin:return None,1+(ticmin-ticmax)/100
    tic=ticmin+u*u*(ticmax-ticmin)
    vi=fixed/(na+k-cl-(ci1+2*ci2)*tic-free)
    if vi<.3:return None,1+(.3-vi)
    oi=na+k+cl+tic+(g.cell_impermeant_osmoles_fmol+g.cell_buffer_total_fmol)/vi
    ob=b.na_mM+b.k_mM+b.cl_mM+b.tic_mM+b.untracked_osmolyte_mM
    ol=ob+(1+w.basolateral_hydraulic_pL_s_mOsm/w.apical_hydraulic_pL_s_mOsm)*(oi-ob)
    q=(w.basolateral_hydraulic_pL_s_mOsm+w.paracellular_hydraulic_pL_s_mOsm*(1+w.basolateral_hydraulic_pL_s_mOsm/w.apical_hydraulic_pL_s_mOsm))*(oi-ob)
    if q<=0:return None,1+abs(q)
    vl=w.lumen_dead_volume_pL+q/w.outflow_rate_s
    l0,l1,l2=_carbon_fractions(phl,a)
    bt=g.lumen_buffer_total_fmol/vl/(1+10**(a.lumen_buffer_pka-phl))+1e3*(10**(phl-a.water_pkw)-10**(-phl))
    tlmax=(ol-bt)/(1+l1+2*l2)
    co2min=tic*ci0;co2max=tlmax*l0
    if co2max<=co2min:return None,1+(co2min-co2max)/10
    tl=(co2min+(co2max-co2min)*v)/l0
    tal=(l1+2*l2)*tl+bt
    if tal<=tl:return None,1+(tl-tal)/100
    cll=(ol-tal-tl)/2
    if cll<=0:return None,1+abs(cll)/100
    vt=p.p.constants.thermal_voltage_V*math.log(tl*l1/p.bath_ab.hco3_mM)+vtmargin
    z=np.r_[[math.log(cl),ph,math.log(tic),math.log(vi),math.log(na/k),
             math.log(cll),phl,math.log(tl),math.log(vl),odds],np.zeros(15),[va*100,(va+vt)*100]]
    return z,0.


def linear_system(p,z):
    # Exclude the two exactly eliminated water rows; keep both current rows.
    rows=[0,1,2,3,5,6,7,8,10,11]
    refz=z.copy();refz[10:25]=0.
    d=p.fast(refz);q=d['flow'];ci=d['cell'];li=d['lumen']
    part=ci['na']/(ci['na']+ci['k'])
    S=np.zeros((14,15))
    for column,entries in enumerate((
        {0:-part,1:-(1-part),2:1,3:-2,4:-2}, {0:1,1:1,2:2},
        {2:1,3:-1,4:-1},{0:1,4:1},{0:-3,1:2,6:3,7:-2,12:1},
        {0:-3,1:2,13:1},{1:-1,7:1,12:1},{1:-1,13:1},
        {2:-1,8:1,12:-1},{3:1},{3:1,9:-1},{6:-1,12:-1,13:1},
        {7:-1,12:-1,13:1},{8:-1,12:1,13:-1},{9:-1,10:-1,12:1,13:-1})):
        for row,value in entries.items():S[row,column]=value
    select=list(CORE_INDEPENDENT_ROWS)+[12,13]
    scales=np.r_[SPEC.independent_rhs_scales,.02,.02]
    A=(S[select]/scales[:,None])[rows]*d['pathway_fluxes'][None,:]/q
    out=np.zeros(14)
    out[6:11]=[li[k] for k in ('na','k','cl','tic','alkalinity')]
    b=(out[select]/scales)[rows]
    return A,b


def inner(p,z):
    A,b=linear_system(p,z)
    scale=np.linalg.norm(A,axis=0)
    if min(scale)<=0:return 1e6,None
    A=A/scale
    # Require a positive AE4 contribution in this numerical search. No WT share target.
    floor=np.zeros(15);floor[0]=1e-6
    y,res=nnls(A,b-A@floor,maxiter=150)
    y+=floor
    capmult=y/scale
    candidate=z.copy();candidate[10:25]=np.log(np.maximum(capmult,1e-100))
    return float(np.linalg.norm(A@y-b)/max(np.linalg.norm(b),1e-30)),candidate


def search(root, measured_means=False):
    require_stage0();manifest,_=load_freeze();p=WTProblem(manifest,root)
    bounds=[(2,60),(60,200),(47.1,53.1),(6.77,7.05),(0,1),(8.19,10.99),
            (1e-5,.99999),(-12,12),(1e-6,.30),(-.30,.30)]
    if measured_means:bounds[2:4]=[(50.1,50.1),(6.91,6.91)]
    def objective(x):
        z,penalty=coordinates(p,x)
        if z is None:return penalty
        try:return inner(p,z)[0]
        except (ValueError,RuntimeError):return 1e6
    started=inverse.now()
    rng=int(hashlib.sha256(root.encode()).hexdigest()[:8],16)
    fit=differential_evolution(objective,bounds,seed=rng,popsize=10,maxiter=400 if measured_means else 220,
        tol=1e-7,polish=False,updating='immediate')
    z,penalty=coordinates(p,fit.x)
    residual,candidate=inner(p,z) if z is not None else (None,None)
    payload=p.production(candidate) if candidate is not None else None
    result=dict(root_id=root,started_utc=started,completed_utc=inverse.now(),
        method='differential_evolution_of_outflow_normalized_NNLS',seed=rng,nfev=fit.nfev,
        message=fit.message,success=bool(fit.success),objective=float(fit.fun),
        bounds=bounds,outer_coordinates=fit.x.tolist(),normalized_flux_residual=residual,
        coordinates=candidate.tolist() if candidate is not None else None,candidate=payload,
        genotype_results_used=False,structural_certificate=False,measured_means=measured_means,
        linear_system_method='analytic_stoichiometric_columns')
    write_json(OUT/('mean_flux_search' if measured_means else 'flux_search')/f'{slug(root)}.json',result)
    print(root,'flux search',residual,payload['resting_gate_failures'] if payload else penalty,flush=True)
    return result


def search_mean(root):return search(root,measured_means=True)


if __name__=='__main__':
    import sys
    require_stage0();manifest,_=load_freeze()
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(search_mean if '--means' in sys.argv else search,sorted(manifest['roots'])))
